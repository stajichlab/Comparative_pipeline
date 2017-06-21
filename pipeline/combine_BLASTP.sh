#!/usr/bin/bash
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

if [ ! $PROJECT ]; then
 echo "Need a project variable in config.txt"
 exit
fi

find . -name "*.m9" | xargs cat | sort -k1,12nr > $PROJECT.BLASTP
