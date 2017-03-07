#!/usr/bin/bash

#SBATCH --time 2-0:00:0
#SBATCH --job-name OrthomclLoad
#SBATCH --ntasks 1
#SBATCH -p batch

module load orthomcl
module load mcl
if [ -f config.txt ]; then
 source config.txt
else
 echo "need a config file to determine prefix"
 exit
fi

if [ ! $PREFIX ]; then
 echo "Need PREFIX defined in config.txt"
 exit
fi
orthomclLoadBlast dbconfig.txt $PREFIX.BLASTP.bpo 
orthomclPairs dbconfig.txt  pairs.log cleanup=yes
orthomclDumpPairsFiles dbconfig.txt

mcl mclInput --abc -I 1.5 -o $PREFIX.OrthoMCL.out

# orthomclMclToGroups
