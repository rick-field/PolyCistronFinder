#!/bin/bash/
#
# Author: Rick Field
#
# email: richard.field@uga.edu
#
# Purpose: to submit many PCF.py jobs



ml pybedtools/0.8.1-foss-2019b
ml Biopython/1.75-intel-2019b-Python-3.7.4
ml pandas/0.25.3-intel-2019b-Python-3.7.4

INFILE=$1

LOCAL_DB_FASTAS="/scratch/rdf85141/PPL/blast_db/Acomosus_Athaliana_Creinhardtii_Gmax_Ptrichocarpa_Osativa_Slycopersicum_Vvinifera_Zmays.fa"

SUB_SCRIPT_DIR=$PWD'/submission_scripts'
if [ -! $SUB_SCRIPT_DIR ]
then
  mkdir $SUB_SCRIPT_DIR
fi

while read SPECIES PREFIX HAP DIR SAMPLE GENE_BED PROTEINS BAM_LIST
do
  OUT_FILE="PCF_$SPECIES_$HAP.out"
  OUT_DIR="$DIRdata/$HAP/"
  SUB_SCRIPT="$PREFIX_$SAMPLE_PCF.sh"
  {
    echo '#!/bin/bash'
    echo "#SBATCH --job-name=$PREFIX_PCF"
    echo '#SBATCH --partition=batch'
    echo '#SBATCH --ntasks=1'
    echo '#SBATCH --nodes=1'
    echo '#SBATCH --cpus-per-task=1'
    echo '#SBATCH --mem=50gb'
    echo '#SBATCH --time=4:00:00'
    echo "#SBATCH --output=PCF.$PREFIX_$SAMPLE.out"
    echo "#SBATCH --error=PCF.$PREFIX_$SAMPLE.err"
    echo ''
    echo 'cd $SLURM_SUBMIT_DIR'
    echo ''
    echo 'ml pybedtools/0.8.1-foss-2019b'
    echo 'ml Biopython/1.75-intel-2019b-Python-3.7.4'
    echo 'ml pandas/0.25.3-intel-2019b-Python-3.7.4'
    echo ''
    echo "python PCF.py -dp ./ -op ./ \"
    echo "-bf $BAM_LIST -bed $DIR/data/$HAP/$GENE_BED \"
    echo "--reference_pep_fastas $DIR/data/$HAP/$PROTEINS -p $PREFIX \"
    echo "-ol 0.5 --local_blast --evalue 1e-20 \"
    echo "--make_local_db_fastas $LOCAL_DB_FASTAS"
    echo '> '$OUT_FILE
  } > $SUB_SCRIPT
#  echo -e $SPECIES'\t'$PREFIX'\t'$HAP'\t'$DIR'\t'$SAMPLE_ID'\t'$READ_FILE'\t'$OUT_FILE >> $OUT_LIST
  sbatch $SUB_SCRIPT
  echo $SUB_SCRIPT' Submitted!'
  mv $SUB_SCRIPT $SUB_SCRIPT_DIR  
