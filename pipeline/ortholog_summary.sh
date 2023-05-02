#!/usr/bin/bash -l
#SBATCH --time 2:00:00 -p short --mem 2G --ntasks 1 -J orthologSum --out ortholog_summary.log
COMPAREFOLDER=./Comparative_pipeline
INFLATION=15
ORTHOINPUT=input
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

# REPLACE WITH ORTHOFINDER INSTEAD OF ORTHOMCL
if [ ! -f goodProteins.fasta ]; then
 module load hmmer/2
 cat $ORTHOINPUT/*.fasta | sreformat fasta - > goodProteins.fasta
fi 
perl $COMPAREFOLDER/scripts/orthomcl_to_table.pl -r $PROJECT.I${INFLATION}.ortholog_counts.tab -i $PROJECT.I${INFLATION}.mcl.out -db goodProteins.fasta --pfam ./domains/Pfam --cazy ./domains/CAZY --merops ./domains/MEROPS
