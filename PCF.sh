#!/bin/bash/
#
# Author: Rick Field
# Email: richard.field@uga.edu
#
# Purpose: to submit many PCF.py jobs

INFILE=$1
LOCAL_DB_FASTAS='/scratch/rdf85141/PPL/blast_db/Acomosus_Athaliana_Creinhardtii_Gmax_Ptrichocarpa_Osativa_Slycopersicum_Vvinifera_Zmays.fa'
SUB_SCRIPT_DIR=$PWD'/submission_scripts'
if [ ! -d $SUB_SCRIPT_DIR ]
then
  mkdir $SUB_SCRIPT_DIR
fi

while read SPECIES PREFIX HAP DIR GENE_BED PROTEINS BAM_LIST CP_GENOME
do
  OUT_FILE='PCF_'$SPECIES'_'$HAP'.out'
  DATA_DIR=$DIR'data/'$HAP'/'
  OUT_DIR=$DIR'analysis/'$HAP'/'
  SUB_SCRIPT=$PREFIX'_PCF.sh'
  {
    echo '#!/bin/bash'
    echo '#SBATCH --job-name='$PREFIX'_PCF'
    echo '#SBATCH --partition=batch'
    echo '#SBATCH --ntasks=1'
    echo '#SBATCH --nodes=1'
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH --mem=120gb'
    echo '#SBATCH --time=24:00:00'
    echo '#SBATCH --output=PCF_'$PREFIX'.out'
    echo '#SBATCH --error=PCF_'$PREFIX'.err'
    echo ''
    echo 'cd $SLURM_SUBMIT_DIR'
    echo ''
    echo 'ml pybedtools/0.8.1-foss-2019b'
    echo 'ml Biopython/1.75-intel-2019b-Python-3.7.4'
    echo 'ml pandas/0.25.3-intel-2019b-Python-3.7.4'
    echo ''
    echo 'python PCF.py -dp '$DATA_DIR' -op '$OUT_DIR' \'
    echo '-bf '$BAM_LIST' -bed '$GENE_BED' \'
    echo '--reference_pep_fastas '$PROTEINS' -p '$PREFIX' \'
    echo '-ol 0.5 --local_blast --evalue 1e-20 --chrom_ignore '$CP_GENOME' \'
    echo '--make_local_db_fastas '$LOCAL_DB_FASTAS' \'
    echo '> '$OUT_FILE
  } > $SUB_SCRIPT
  sbatch $SUB_SCRIPT
  echo $SUB_SCRIPT' Submitted!'
  mv $SUB_SCRIPT $SUB_SCRIPT_DIR  
done < $INFILE
