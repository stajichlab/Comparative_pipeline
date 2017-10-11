#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 2 --mem-per-cpu=1G --time 12:00:00
#SBATCH --job-name=CAZY
#SBATCH --output=domains.CAZY.%A_%a.log

DOMAINS=domains
EXT=aa.fasta
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

mkdir -p $DOMAINS/CAZY

if [ ! $EXT ]; then
 EXT=aa.fasta
fi

if [ ! $PROTEINS ]; then
 PROTEINS=pep
fi

module load hmmer/3
module load db-cazy

if [ ! $CAZY_DB ]; then
 echo "Need a CAZY_DB env variable either from config.txt or 'module load db-cazy'"
 exit
fi
CPUS=$SLURM_CPUS_ON_NODE
if [ ! $CPUS ]; then
 CPUS=1
fi

IN=${SLURM_ARRAY_TASK_ID}

if [ ! $IN ]; then
 IN=$1
 if [ ! $IN ]; then
   IN=1
   echo "defaulting to IN value is 1 - specify with --array or cmdline"
 fi
fi

TOTAL=$(ls $PROTEINS/*.${EXT} | wc -l)
if [ $IN -gt $TOTAL ]; then
 echo "Only $TOTAL files in folder $PROTEINS, skipping $IN"
 exit
fi
INFILE=$(ls $PROTEINS/*.${EXT} | sed -n ${IN}p)
OUT=$DOMAINS/CAZY/$(basename ${INFILE} .${EXT})

if [ ! -f ${OUT}.hmmscan ]; then
 hmmscan --cpu $CPUS --domtbl ${OUT}.domtbl -o ${OUT}.hmmscan $CAZY_DB $INFILE
 bash $CAZY_FOLDER/hmmscan-parser.sh ${OUT}.domtbl | sort > ${OUT}.tsv
fi
