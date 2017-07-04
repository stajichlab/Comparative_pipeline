#!/usr/bin/bash
#SBATCH --mem 64G --ntasks 8 --nodes 1 -p intel --time 24:00:00

TMPDIR=/scratch
INDIR=pair_compare
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

if [ ! -d $TMPDIR ]; then
 TMPDIR=/tmp
fi

if [ ! $PROJECT ]; then
 echo "Need a project variable in config.txt"
 exit
fi
if [ ! -f $PROJECT.BLASTP ]; then

rm -f $PROJECT.BLASTP
for file in $(ls input/*.fasta)
do
 base=$(basename $file .fasta)
 find -L pair_compare/$base-* -name "*.m9.gz"  | xargs zcat | time sort -T $TMPDIR --parallel $CPU --buffer-size=64G -k1,1 -k12,12nr >> $PROJECT.BLASTP
done

fi
