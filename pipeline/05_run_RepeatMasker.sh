#!/usr/bin/bash -l
#SBATCH --mem 4G --ntasks 16 --nodes 1 --time 1-0:00:00 -J RM --out=RM.%A_%a.out
module load RepeatMasker

OUTDIR=Repeats
CPUS=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPUS=${SLURM_ARRAY_TASK_ID}
fi

RMSPECIES=fungi
if [ -f config.txt ]; then
 source ./config.txt
fi

if [ ! $GENOME ]; then
 echo "No Genome dir defined in config.txt"
 exit
fi

mkdir -p $OUTDIR

N=1
if [ $1 ]; then
 N=$1
elif [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
fi

GENOMEFILE=$(ls $GENOME/*.fasta | sed -n ${N}p)
STEM=$(basename $GENOMEFILE .fasta)
pushd $OUTDIR
if [ ! -e ${STEM}.fasta ]; then 
 ln -s ../$GENOMEFILE .
fi
if [ ! -f ${STEM}.RM.out ]; then
 RepeatMasker -pa $CPUS -e ncbi -species $RMSPECIES -s ${STEM}.fasta > ${STEM}.RM.out
fi
