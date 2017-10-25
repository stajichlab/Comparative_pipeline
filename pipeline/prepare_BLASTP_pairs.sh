#!/usr/bin/bash
#SBATCH --time 2:00:00 --mem 1G --ntasks 1 -J prepareBLASTP
# BASH script to setup pairwise BLASTP runs

module load ncbi-blast/2.6.0+
CPU=1
DBSIZE=1000000
E=1e-3
SIZE=5000
INDIR=input
TARGET=pair_compare
SPLIT=split
BLASTEXT=m9
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
 if [ ! -f $db.phr ]; then
  makeblastdb -dbtype prot -in $db
 fi
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
    if [ ! -f $DEST.$BLASTEXT.gz -a ! -s $DEST.$BLASTEXT.gz ]; then
     echo "no $DEST.$BLASTEXT.gz"
     echo "$DEST $qi $db" >> $JOBSFILE
    fi
   done
 done
done
