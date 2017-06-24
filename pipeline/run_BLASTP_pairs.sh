#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 2 --cpus-per-task 1 -p short --time 2:00:00 --output blastp.%A_%a.out

# THIS SCRIPT IS FOR RUNNING BLASTP for all pairs of split files

module load ncbi-blast/2.6.0+
module load orthomcl
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
CPUS=$SLURM_CPUS_ON_NODE

N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
 N=$1
fi

if [ ! $N ]; then
 N=1
fi

sed -n ${N}p $JOBS | while read DEST QUERY TARGET
do
 echo "$DEST.m9"
 if [ ! -f $DEST.m9 -a ! -s $DEST.m9 ]; then
  time $BLAST -db $TARGET -query $QUERY -dbsize $DBSIZE -out $DEST.m9 -outfmt 6 -num_threads $CPU -evalue $E
 fi
 if [ ! -f $DEST.bpo ]; then
  orthomclBlastParser $DEST.m9 $INFASTA > $DEST.bpo
 fi
done
