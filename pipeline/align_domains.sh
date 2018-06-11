#!/bin/bash

module load hmmer
module load trimal
ODIR=domain_aln
Pfamdb="/srv/projects/db/PFAM/2017-06-11-Pfam31.0/Pfam-A.hmm"
for dom in AMP-binding
do
    perl ../Comparative_pipeline/scripts/extract_domain.pl -i $dom
    cat $ODIR/$dom/*.fasta > domain_aln/$dom.fas
    if [ ! -f $ODIR$dom.hmm ]; then
	hmmfetch $Pfamdb $dom > $ODIR/$dom.hmm
    fi
    ALN=$ODIR/$dom.aln
    hmmalign --trim --amino $ODIR/$dom.hmm $ODIR/$dom.fas > $ALN
    esl-reformat --replace=x:- --gapsym=- -o $dom.tmp afa $ALN
    perl -i -p -e 'if (! /^>/) { s/[ZBzbXx\*]/-/g }' $dom.tmp 
    trimal -resoverlap 0.50 -seqoverlap 60 -in $dom.tmp -out $dom.trim1
    trimal -automated1 -in $dom.trim1 -out $ODIR/$dom.trim
    unlink $dom.trim1 $dom.tmp
done
    
