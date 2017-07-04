for file in *.fasta; do pre=$(head -q -n1 $file | awk -F\| '{print $1}' | perl -p -e 's/^>//'); mv $file $pre.fasta; done
