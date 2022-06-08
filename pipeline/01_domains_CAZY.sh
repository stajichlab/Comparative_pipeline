#!/usr/bin/bash -l
#SBATCH --nodes 1 --ntasks 24 -p short -C xeon --mem=8G --time 2:00:00
#SBATCH --job-name=CAZY
#SBATCH --output=logs/domains.CAZY.%a.log
module load workspace/scratch
export OMPI_TMPDIR=$SCRATCH
hostname
DOMAINS=domains
EXT=aa.fasta
if [ -f config.txt ]; then
 source config.txt
else
 echo "need config file to set some project-specific variables"
 exit
fi

mkdir -p $DOMAINS/CAZY

if [ ! $EXT ]; then
 EXT=aa.fasta
fi

if [ ! $PROTEINS ]; then
 PROTEINS=pep
fi

module load db-cazy

if [ ! $CAZY_DB ]; then
 echo "Need a CAZY_DB env variable either from config.txt or 'module load db-cazy'"
g exit
fi
CPUS=${SLURM_CPUS_ON_NODE}
if [ ! $CPUS ]; then
 CPUS=1
fi
IN=${SLURM_ARRAY_TASK_ID}

if [ -z $IN ]; then
 IN=$1
 if [ -z $IN ]; then
   IN=1
   echo "defaulting to IN value is 1 - specify with --array or cmdline"
 fi
fi

TOTAL=$(ls $PROTEINS/*.${EXT} | wc -l)
if [ $IN -gt $TOTAL ]; then
 echo "Only $TOTAL files in folder $PROTEINS, skipping $IN"
 exit
fi
INFILE=$(ls $PROTEINS/*.${EXT} | sed -n ${IN}p)
OUT=$DOMAINS/CAZY/$(basename ${INFILE} .${EXT})

if [ ! -f ${OUT}.hmmscan ]; then
    module load hmmer/3.3.2-mpi
    # hmmscan --cpu $CPUS --domtbl ${OUT}.domtbl -o ${OUT}.hmmscan $CAZY_DB $INFILE
    srun hmmscan --mpi --domtbl ${OUT}.domtbl -o ${OUT}.hmmscan $CAZY_DB $INFILE
    bash $CAZY_FOLDER/hmmscan-parser.sh ${OUT}.domtbl | sort > ${OUT}.tsv
    module unload hmmer
fi

if [[ ! -d $OUT.run_dbcan || ! -f $OUT.run_dbcan/overview.txt ]]; then
    module load run_dbcan    
    module load hmmer/3
    run_dbcan --db_dir $CAZY_FOLDER --out_dir $OUT.run_dbcan --tools all \
	--stp_cpu $CPUS --eCAMI_jobs $CPUS --hmm_cpu $CPUS --dia_cpu $CPUS --tf_cpu $CPUS \
	$INFILE protein 
fi
