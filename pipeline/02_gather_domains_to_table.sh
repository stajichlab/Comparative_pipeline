#!/usr/bin/bash
#SBATCH -p short -n 1 --mem 8gb --out logs/gather_MEROPS.log
module load hmmer/3
module load miniconda3
perl scripts/gather_MEROPS_family_counts.pl --db pep --sfetch $(which esl-sfetch)
perl scripts/gather_hmmscanMPI_domain_counts.pl
perl scripts/gather_CAZY_counts.pl
perl scripts/gather_CAZY_rundbcan_counts.pl
