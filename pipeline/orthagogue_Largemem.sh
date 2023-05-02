#!/usr/bin/bash -l
#SBATCH --time 2-0:00:00 --ntasks 8 --mem 350G
module load orthagogue
module load mcl

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi
 #THIS REALLY SHOULD BE USING ORTHOFINDER INSTEAD
CPU=8
if  [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [ ! $PROJECT ]; then
 echo "Need a project variable in config.txt"
 exit
fi

if [ ! -f all.abc ]; then
 orthAgogue -i $PROJECT.BLASTP -dbs 100000000 -s '|' -e 6 -c $CPU
fi

mcl all.abc -te $CPU --abc -I 1.5 -o $PROJECT.I15.mcl.out
