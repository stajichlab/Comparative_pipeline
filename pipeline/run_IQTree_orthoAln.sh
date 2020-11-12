#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 2 --mem 2G --time 2:00:00 -p short -J IQTREE.Orthologs

CPU=2
module load IQ-TREE
N=0
if [ ${SLURM_ARRAY_TASK_ID} ]; then
 N=${SLURM_ARRAY_TASK_ID}
elif [ $1 ]; then
 N=$1
fi
EXT=trim
ALNEXT=treefile
LIST=alnfiles.txt
JOBSIZE=5
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
 DEST=$FILENAME.$ALNEXT
 if [ ! -f $DEST ]; then
  iqtree-omp -s $FILENAME -bb 1000 -nt $CPU 
 fi
 done
done
