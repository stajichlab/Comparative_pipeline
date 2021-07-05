#!/usr/bin/bash
module load db-merops

if [[ ! -z $MEROPS_DB && -f $MEROPS_DB/merops_scan.lib ]]; then
    grep ">" $MEROPS_DB/merops_scan.lib | perl -p -e 's/^>(\S+).+\[([^]]+)\]\#([^\#]+)\#.+/$1\t$3\t$2/' > lib/merops_lib.families.tab
fi
