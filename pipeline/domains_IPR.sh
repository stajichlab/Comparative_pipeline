#!/bin/bash -l
#SBATCH --nodes 1 --ntasks 24 --mem 96G --out logs/IPR.%A_%a.log -J IPR
mkdir -p logs
module load funannotate
module load iprscan

DOMAINS=domains
EXT=aa.fasta
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi


mkdir -p $DOMAINS/IPR

CPUS=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPUS=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
	IN=1
	echo "defaulting to IN value is 1 - specify with --array or cmdline"
    fi
fi
TOTAL=$(ls $PROTEINS | grep -c -P "${EXT}$" )
if [ $N -gt $TOTAL ]; then
    echo "Only $TOTAL files in folder $PROTEINS, skipping $IN"
    exit
fi
INFILE=$(ls $PROTEINS/*.${EXT} | sed -n ${IN}p)
XML=$DOMAINS/IPR/$(basename ${INFILE} .${EXT}).xml
TAB=$DOMAINS/IPR/$(basename ${INFILE} .${EXT}).tab
IPRPATH=$(which interproscan.sh)
if [ ! -s $XML ]; then    
     funannotate iprscan -i $INFILE -o $XML -m local -c $CPUS --iprscan_path $IPRPATH
fi
if [ ! -s $TAB ]; then
    $IPRPATH -mode convert -i $XML -o $TAB -f TSV
fi
