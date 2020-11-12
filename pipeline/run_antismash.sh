#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 24 --mem 24G -p intel --time 3-0:00:00
#SBATCH --job-name=AntiSMASH
#SBATCH --output=AntiSMASH.%A_%a.log

CPUS=2

if [ $SLURM_CPUS_ON_NODE ]; then
 CPUS=$SLURM_CPUS_ON_NODE
fi

GENBANK=gbk
EXT=gbk

N=1
if [ $1 ]; then
 N=$1
elif [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
fi

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

OUTDIR=secondary_metabolite
if [ ! -d $OUTDIR ]; then
 mkdir -p $OUTDIR
fi
TOTAL=$(ls $GENBANK/*.${EXT} | wc -l)
if [ $N -gt $TOTAL ]; then
 echo "Only $TOTAL files in folder $GENBANK, skipping $N"
 exit
elif [[ $N == 0 ]]; then
 echo "N must be between 1 and $TOTAL"
 exit
fi
INFILE=$(ls $GENBANK/*.${EXT} | sed -n ${N}p)
echo "INFILE=$INFILE"
OUT=$OUTDIR/$(basename ${INFILE} .${EXT})
module unload perl
module unload perl
module load antismash/4.1.0
module unload python/3
source activate antismash
CPU=$SLURM_CPUS_ON_NODE

antismash --taxon fungi -c $CPUS --outputfolder $OUT   --clusterblast  --subclusterblast --smcogs --knownclusterblast \
 --borderpredict --asf \
 --full-hmmer --cassis --clusterblast --smcogs --subclusterblast --knownclusterblast $INFILE 
