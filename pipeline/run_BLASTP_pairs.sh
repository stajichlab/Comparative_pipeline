#!/usr/bin/bash

#SBATCH --mem 2G --nodes 1 --ntasks 8 --cpus-per-task 1 -p intel --time 24:00:00 --output blastp.%A_%a.out

# THIS SCRIPT IS FOR RUNNING BLASTP for all pairs of split files

module load ncbi-blast/2.6.0+
module load pigz
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

DBSIZE=10000000
EVALUE=0.0001
BLAST=blastp
JOBS=jobs.cmds
OUTDIR=pair_compare
INFASTA=in
JOBSIZE=5
CPUS=$SLURM_CPUS_ON_NODE
N=0
if [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ $1 ]; then
 N=$1
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
  if [ ! -f $DEST.m9.gz -a ! -s $DEST.m9.gz ]; then
   time $BLAST -db $TARGET -query $QUERY -dbsize $DBSIZE -out $DEST.m9 -outfmt 6 -num_threads $CPU -evalue $EVALUE
   pigz $DEST.m9
  fi
# if [ ! -f $DEST.bpo ]; then
#  orthomclBlastParser $DEST.m9 $INFASTA > $DEST.bpo
# fi
 done
done
