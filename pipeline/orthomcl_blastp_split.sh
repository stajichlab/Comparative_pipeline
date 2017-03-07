#!/usr/bin/bash

#SBATCH --ntasks 4 --nodes 1 --time 6:00:00 --job-name=blastp.ortho.split
#SBATCH --output=blastp.ortho.%A_%a.out

module load ncbi-blast/2.6.0+


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
OUT=BLASTP
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=1
fi

if [ ! -f $OUT/${PREF}.${N}.BLASTP.tab ]; then
 blastp -num_threads $CPU -db $DB -query $FOLDER/${PREF}.$N -outfmt 6 -evalue 0.01 -out $OUT/${PREF}.${N}.BLASTP.tab
fi
