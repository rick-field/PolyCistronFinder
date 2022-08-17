#!/bin/bash/


ml pybedtools/0.8.1-foss-2019b
ml Biopython/1.75-intel-2019b-Python-3.7.4
ml pandas/0.25.3-intel-2019b-Python-3.7.4


INFILE=$1

while read SPECIES PREFIX HAP DIR SAMPLE GENE_BED PROTEINS BAM_LIST
do
  OUT_FILE=${IN_FILE%.*}'_sorted.bam'
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
    echo "-bf $BAM_LIST -bed $GENE_BED \"
    echo "--reference_pep_fastas $PROTEINS -p $PREFIX \"
    echo "-f 7 -ol 0.5 --local_blast --evalue 1e-20 \"
    echo "--make_local_db_fastas $LOCAL_DB_FASTAS"
    echo '> '$OUT_FILE
  } > $SUB_SCRIPT
