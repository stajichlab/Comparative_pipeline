#!/usr/bin/bash

#SBATCH --time 2-0:00:00 --ntasks 8 --mem 32G
module load orthagogue
module load mcl

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi
ORTHOFOLDER=orthologs
CPU=2
if  [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [ ! $PROJECT ]; then
 echo "Need a project variable in config.txt"
 exit
fi
mkdir -p $ORTHOFOLDER
pushd $ORTHOFOLDER

if [ ! -f all.abc ]; then
 orthAgogue -i ../$PROJECT.BLASTP -s '|' -e 6 -c $CPU
fi

mcl all.abc -te $CPU --abc -I 1.5 -o $PROJECT.I15.mcl.out
