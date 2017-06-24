#!/usr/bin/bash

if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

perl ../Comparative_pipeline/scripts/orthomcl_to_table.pl -r $PROJECT.I15.ortholog_counts.tab -i $PROJECT.I15.mcl.out -db input --pfam ../domains/Pfam
