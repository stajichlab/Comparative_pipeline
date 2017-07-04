#!/usr/bin/bash
#SBATCH --time 2:00:00 --mem 1G --ntasks 1 -J prepareBLASTP
# BASH script to setup pairwise BLASTP runs

CPU=1
DBSIZE=1000000
E=1e-3
SIZE=5000
INDIR=input
TARGET=pair_compare
SPLIT=split
EXT=m9
JOBSFILE=jobs.cmds
# can override these variables by defining them in config.txt

if [ ! -f config.txt ]; then
 if [ -f ../config.txt ]; then
  ln -s ../config.txt .
 fi
fi

if [ -f config.txt ]; then
 source config.txt
fi

for l in $INDIR/*.fasta
do
 left=$(basename $l .fasta)
 if [ ! -d $SPLIT/$left ]; then
  mkdir -p $SPLIT/$left
  bp_dbsplit.pl --prefix $SPLIT/$left/$left --size ${SIZE} ${l} 
 fi
done

# remove existing jobs set
rm -f $JOBSFILE
for db in $INDIR/*.fasta
do
 d=$(basename $db .fasta)
 for query in $INDIR/*.fasta
 do
   q=$(basename $query .fasta)
   for qi in $SPLIT/$q/$q.*
   do
    qname=$(basename $qi)
    if [ ! -d $TARGET/${q}-${d} ]; then
     mkdir -p $TARGET/${q}-${d}
    fi
    DEST=$TARGET/${q}-${d}/${qname}__${d}
    if [ ! -f $DEST.$EXT.gz -a ! -s $DEST.$EXT.gz ]; then
     echo "$DEST $qi $db" >> $JOBSFILE
    fi
   done
 done
done
