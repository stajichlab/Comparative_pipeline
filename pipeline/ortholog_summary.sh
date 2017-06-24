#!/usr/bin/bash
#SBATCH --time 2:00:00 -p short --mem 2G --ntasks 1 -J orthologSum --out ortholog_summary.log
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi
if [ ! -f goodProteins.fasta ]; then
 module load hmmer/2
 cat in/*.fasta | sreformat fasta - > goodProteins.fasta
fi 
perl ../Comparative_pipeline/scripts/orthomcl_to_table.pl -r $PROJECT.I15.ortholog_counts.tab -i $PROJECT.I15.mcl.out -db goodProteins.fasta --pfam ../domains/Pfam
