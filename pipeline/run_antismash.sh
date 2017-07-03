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
module unload openmpi
module unload python
module load antismash/4.0-branch
export PYTHONHOME=$VIRTUAL_ENV
echo "$VIRTUAL_ENV"
srun -J AntiSmash --out AntiSmash.%A_$N.out --ntasks 24 --nodes 1 --mem 24G run_antismash.py -c $CPUS --taxon fungi --input-type nucl --clusterblast  --subclusterblast --smcogs --knownclusterblast --cassis --borderpredict --asf --outputfolder $OUT $INFILE

