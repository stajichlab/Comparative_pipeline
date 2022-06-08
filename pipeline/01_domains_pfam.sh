#!/usr/bin/bash -l
#SBATCH --nodes 1 --ntasks 48 --mem 32gb  -p short -C ryzen
#SBATCH --job-name=Pfam.domains
#SBATCH --output=logs/domains.Pfam.%a.log

EXT=aa.fasta
PROTEINS=pep
DOMAINS=domains

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi
mkdir -p $DOMAINS/Pfam

module load db-pfam
module load hmmer/3.3.2-mpi

echo "running $PFAM_DB"

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

TOTAL=$(ls $PROTEINS/*.${EXT} | wc -l)
if [ $IN -gt $TOTAL ]; then
 echo "Only $TOTAL files in folder $PROTEINS, skipping $IN"
 exit
fi
INFILE=$(ls $PROTEINS/*.${EXT} | sed -n ${IN}p)
OUT=$DOMAINS/Pfam/$(basename ${INFILE} .${EXT})

if [ ! -f ${OUT}.hmmscan ]; then
 srun hmmsearch --mpi --cut_ga --domtbl ${OUT}.domtbl -o ${OUT}.hmmscan $PFAM_DB/Pfam-A.hmm $INFILE
 #hmmscan --cut_ga --cpu $CPUS --domtbl ${OUT}.domtbl -o ${OUT}.hmmscan $PFAM_DB/Pfam-A.hmm $INFILE
fi
