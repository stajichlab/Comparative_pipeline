#!/usr/bin/bash


#SBATCH --nodes 1 --ntasks 2 --mem-per-cpu=1G
#SBATCH --job-name=domains.Pfam
#SBATCH --output=domains.Pfam.%A_%a.log


if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

if [ $DOMAINS ]; then
 OUTDIR=$DOMAINS
else
 OUTDIR=domains
fi

mkdir -p $OUTDIR/Pfam

if [ ! $EXT ]; then
 EXT=aa.fasta
fi

if [ ! $PROTEINS ]; then
 PROTEINS=pep
fi
module unload python
module load antismash/4.0-branch


