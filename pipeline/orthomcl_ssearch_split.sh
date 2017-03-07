#!/usr/bin/bash

#SBATCH --ntasks 4 --nodes 1 --time 6:00:00 --job-name=fasta.ortho.split
#SBATCH --output=fasta.ortho.%A_%a.out

module load fasta

N=${SLURM_ARRAY_TASK_ID}
if [ ! $N ]; then
 N=$1
 if [ ! $N ]; then
  echo "need an array or cmdline input job id default to 1"
  N=1
 fi
fi

FOLDER=split
DB=goodProteins.fasta
PREF=db
OUT=FASTA
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=1
fi

if [ ! -f $OUT/${PREF}.${N}.tab ]; then
 ssearch36 -T $CPU -E 0.01 -m 8c -d 0 $FOLDER/${PREF}.$N $DB > $OUT/${PREF}.${N}.tab
fi
