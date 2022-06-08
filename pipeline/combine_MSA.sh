#!/usr/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=2:00:00
#SBATCH --job-name=combineThemAll

module load RAxML
EXPECTEDNAMES=expected
MARKERS=singlecopy
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some of the pipeline variables"
 exit
fi

if [ ! -f $EXPECTEDNAMES ]; then
 head -q -n1 $ORTHOINPUT/*.fasta | awk -F\| '{print $1}' | awk '{print $1}' > $EXPECTEDNAMES
fi

ALN=$PROJECT.I${INFLATION}.single_copy_orthologs
count=$(wc -l $EXPECTEDNAMES | awk '{print $1}')
#echo "count is $count"

perl scripts/combine_fasaln.pl -ext trim -o all_${count}.${MARKERS}.fasaln -of fasta -d $ALN -expected $EXPECTEDNAMES > all_${count}.${MARKERS}.partitions.txt
convertFasta2Phylip.sh all_${count}.${MARKERS}.fasaln > all_${count}.${MARKERS}.phy
