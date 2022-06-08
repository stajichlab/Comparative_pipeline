#!/usr/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --job-name=buildTreeRX
#SBATCH --time=4-0:00:00
#SBATCH --mem 12G
#SBATCH --output=buildTree.RAxML.%A.out

module load RAxML

CPU=2

RUNFOLDER=phylo
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

if [ ! $OUTGROUP  ]; then
 echo "need an OUTGROUP in config.txt"
 exit
fi
 
mkdir -p $RUNFOLDER
count=`wc -l expected | awk '{print $1}'`
datestr=`date +%Y_%b_%d`
str=$PROJECT.${count}sp
IN=all_${count}.singlecopy
if [ ! -f phylo/$str.fasaln ]; then
 cp $IN.fasaln phylo/$str.fasaln
 cp $IN.partitions.txt phylo/$str.partitions
 cp $IN.phy phylo/$str.phy
fi
cd phylo
raxmlHPC-PTHREADS-AVX -T $CPU -f a -x 227 -p 771 -o $OUTGROUP -m PROTGAMMAAUTO -s $str.fasaln -n $datestr.$str -N autoMRE 
