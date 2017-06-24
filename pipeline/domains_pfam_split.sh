#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 2 --mem-per-cpu=1G
#SBATCH --job-name=domains.Pfam
#SBATCH --time=1-0:00:00
#SBATCH --output=domains.Pfam.%A_%a.log


PROTEINS=split
if [ -f splitconfig.txt ]; then
 source splitconfig.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

if [ $DOMAINS ]; then
 OUTDIR=$DOMAINS
else
 OUTDIR=splitdomains
fi

mkdir -p $OUTDIR

module load db-pfam
module load hmmer/3

if [ ! $PFAM_DB ]; then
 echo "Need a PFAM_DB env variable either from config.txt or 'module load db-pfam'"
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

INFILE=$PROTEINS/$PREFIX.$IN
if [ ! -f $INFILE ]; then
 echo "No file $INFILE"
 exit
fi
OUT=$DOMAINS/$PREFIX.$IN

if [ ! -f ${OUT}.hmmscan ]; then
 hmmscan --cut_ga --cpu $CPUS --domtbl ${OUT}.domtbl -o ${OUT}.hmmscan $PFAM_DB/Pfam-A.hmm $INFILE
fi
