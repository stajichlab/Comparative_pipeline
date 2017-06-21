#!/usr/bin/bash

#SBATCH --ntasks 1 --nodes 1 --time 1:00:00 --job-name=parse.orthomcl.split
#SBATCH --output=parseblast.%A_%a.out -p batch

module load orthomcl

N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
 N=$1
 if [ ! $N ]; then
  echo "need an array or cmdline input job id default to 1"
  N=1
 fi
fi

PROTEINS=input
PREFIX=db
OUT=BLASTP

if [ -f config.txt ]; then
 source config.txt
fi

if [ ! -f $OUT/${PREFIX}.${N}.BLASTP.bpo ]; then
 orthomclBlastParser $OUT/${PREFIX}.${N}.BLASTP.tab $PROTEINS > $OUT/${PREFIX}.${N}.BLASTP.bpo
fi
