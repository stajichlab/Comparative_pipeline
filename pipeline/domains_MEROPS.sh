#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 4 --mem-per-cpu=1G --time 2:00:00 -p short
#SBATCH --job-name=MEROPS.domains
#SBATCH --output=logs/domains.MEROPS.%a.log
mkdir -p logs

OUTEXT=blasttab
PROTEINS=pep
EXT=aa.fasta
DOMAINS=domains
MEROPS_CUTOFF=1e-10
MEROPS_MAX_TARGETS=10
if [ -f config.txt ]; then
    source config.txt
else
    echo "need config file to set some project-specific variables"
    exit
fi

mkdir -p $DOMAINS/MEROPS

module load db-merops/120
module load ncbi-blast/2.9.0+

if [ ! $MEROPS_DB ]; then
    echo "Need a MEROPS_DB env variable either from config.txt or 'module load db-merops'"
    exit
fi

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
    CPUS=1
fi

IN=${SLURM_ARRAY_TASK_ID}

if [ -z $IN ]; then
    IN=$1
    if [ -z $IN ]; then
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
OUT=$DOMAINS/MEROPS/$(basename ${INFILE} .${EXT}).${OUTEXT}

if [ ! -f ${OUT} ]; then
    blastp -query $INFILE -db $MEROPS_DB/merops_scan.lib -out ${OUT} -num_threads $CPUS -seg yes -soft_masking true -max_target_seqs $MEROPS_MAX_TARGETS -evalue $MEROPS_CUTOFF -outfmt 6 -use_sw_tback
fi
