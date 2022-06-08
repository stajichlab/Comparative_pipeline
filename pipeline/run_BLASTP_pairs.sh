#!/usr/bin/bash -l

#SBATCH --mem 2G --nodes 1 --ntasks 8 --cpus-per-task 1 -p intel --time 24:00:00 --output blastp.%A_%a.out

# THIS SCRIPT IS FOR RUNNING BLASTP for all pairs of split files
# THIS COULD BE BETTER DONE WITH DIAMOND AND MMSEQS2 FOR SPEED
module load ncbi-blast
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

# THIS SHOULD BE UPDATED SO THAT THESE VARIABLES CAN ALSO BE RESET IN A blast_config.txt file perhaps?
DBSIZE=10000000
EVALUE=0.0001
BLAST=blastp
JOBS=jobs.cmds
OUTDIR=pair_compare
JOBSIZE=5
CPUS=$SLURM_CPUS_ON_NODE
N=0
if [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ $1 ]; then
 N=$1
fi

# allow customizing variables with local config file
if [ -f config.txt ]; then
 source config.txt
fi

MAX=$(wc -l $JOBS | awk '{print $1}')
B=$(expr $N \* $JOBSIZE)
E=$(expr $B + $JOBSIZE - 1)
if [ "$E" -gt "$MAX" ]; then
 E=$MAX
fi
if [[ $B == 0 ]]; then 
 B=1
fi

for R in $(seq $B 1 $E);
do
 sed -n ${R}p $JOBS | while read DEST QUERY TARGET
 do
  echo "$DEST.m9"
  destdir=$(dirname $DEST)
  if [ ! -d $destdir ]; then
   mkdir -p $destdir
  fi 
  if [ ! -f $DEST.m9.gz -a ! -s $DEST.m9.gz ]; then
   time $BLAST -db $TARGET -query $QUERY -dbsize $DBSIZE -out $DEST.m9 -outfmt 6 -num_threads $CPU -evalue $EVALUE
   pigz $DEST.m9
  fi
 done
done
