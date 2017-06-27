#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 2 --cpus-per-task 1 -p short --time 2:00:00 --output blastp.%A_%a.out

# THIS SCRIPT IS FOR RUNNING BLASTP for all pairs of split files

module load ncbi-blast/2.6.0+
#module load orthomcl
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

DBSIZE=10000000
E=1e-3
BLAST=blastp
JOBS=jobs.cmds
OUTDIR=pair_compare
INFASTA=in
JOBSIZE=10
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
  if [ ! -f $DEST.m9 -a ! -s $DEST.m9 ]; then
   time $BLAST -db $TARGET -query $QUERY -dbsize $DBSIZE -out $DEST.m9 -outfmt 6 -num_threads $CPU -evalue $E
  fi
# if [ ! -f $DEST.bpo ]; then
#  orthomclBlastParser $DEST.m9 $INFASTA > $DEST.bpo
# fi
 done
done
