#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 2G --time 2:00:00 -p short -J MSA

module load muscle
module load trimal
MSATOOL="muscle -quiet"
N=0
if [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ $1 ]; then
 N=$1
fi
EXT=fa
ALNEXT=fasaln
LIST=files.txt
JOBSIZE=10
if [ ! -f $LIST ]; then
 ls *.${EXT} > $LIST
fi

MAX=$(wc -l $LIST | awk '{print $1}')

B=$(expr $N \* $JOBSIZE)
E=$(expr $B + $JOBSIZE - 1)
if [ "$E" -gt "$MAX" ]; then
 E=$MAX
fi
if [[ $B == 0 ]]; then
 B=1
fi

for R in $(seq $B 1 $E);
do
 sed -n ${R}p $LIST | while read FILENAME
 do
 STEM=$(basename $FILENAME .${EXT})
 DEST=$STEM.$ALNEXT
 TRIM=$STEM.trim
 if [ ! -f $DEST ]; then
   $MSATOOL -in $FILENAME -out $DEST
   trimal -in $DEST -out $TRIM -automated1
   #perl -i -p -e 's/>(\S+)\|/>$1 /' $TRIM
 fi
 done
done
