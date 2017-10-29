#!/usr/bin/bash
#SBATCH --time 2:00:00 -p short --mem 4G --ntasks 1 -J orthSingleCopy --out ortholog_singlecopy.log
COMPAREFOLDER=./Comparative_pipeline
INFLATION=15
ORTHOINPUT=input
ALLOWED_MISSING_SINGLECOPY=1
ORTHODIR=orthologs
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

module load cdbfasta
ALLPROTEINS=goodProteins.fasta
if [ ! -f $ALLPROTEINS ]; then
 module load hmmer/2
 cat $ORTHOINPUT/*.fasta | sreformat fasta - > $ALLPROTEINS
fi 
if [ ! -f $ALLPROTEINS.cidx ]; then
 cdbfasta $ALLPROTEINS
fi
perl $COMPAREFOLDER/scripts/orthomcl_extract_goodorth.pl -r $PROJECT.I${INFLATION}.single_copy_orthologs.tab -i $ORTHODIR/$PROJECT.I${INFLATION}.mcl.out -db $ALLPROTEINS --missing $ALLOWED_MISSING_SINGLECOPY -o $PROJECT.I${INFLATION}.single_copy_orthologs
