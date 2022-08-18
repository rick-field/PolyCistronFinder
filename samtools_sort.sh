#!/bin/zsh

# Author: Rick Field
# Email: richard.field@uga.edu
# Purpose: to submit many samtools sort jobs

IN_LIST=$1
OUT_LIST='./samtools_sort_out.txt'
if [ -e $OUT_LIST ]
then
  rm $OUT_LIST
fi

SUB_SCRIPT_DIR=$PWD'/submission_scripts/'
if [ ! -d $SUB_SCRIPT_DIR ]
then
  mkdir	$SUB_SCRIPT_DIR
fi

while read SPECIES PREFIX HAP DIR SAMPLE_ID READ_FILE IN_FILE
do
  OUT_FILE=${IN_FILE%.*}'_sorted.bam'
  OUT_DIR=$DIR'data/'$HAP'/'
  SUB_SCRIPT=$PREFIX'_'$HAP'_'$SAMPLE_ID'_samtools_sort.sh'
  {
    echo '#!/bin/bash'
    echo '#SBATCH --job-name='$SUB_SCRIPT
    echo '#SBATCH --partition=batch'
    echo '#SBATCH --ntasks=1'
    echo '#SBATCH --nodes=1'
    echo '#SBATCH --cpus-per-task=4'
    echo '#SBATCH --mem=50gb'
    echo '#SBATCH --time=12:00:00'
    echo '#SBATCH --output='${SUB_SCRIPT%.*}'.out'
    echo '#SBATCH --error='${SUB_SCRIPT%.*}'.err'
    echo ''
    echo 'cd $SLURM_SUBMIT_DIR'
    echo ''
    echo 'ml SAMtools/1.14-GCC-8.3.0'
    echo ''
    echo 'samtools sort -@ 4 \'
    echo '-O bam \'
    echo '-o '$OUT_DIR$OUT_FILE' \'
    echo $OUT_DIR$IN_FILE
  } > $SUB_SCRIPT
  echo -e $SPECIES'\t'$PREFIX'\t'$HAP'\t'$DIR'\t'$SAMPLE_ID'\t'$READ_FILE'\t'$OUT_FILE >> $OUT_LIST
  sbatch $SUB_SCRIPT
  echo $SUB_SCRIPT' Submitted!'
  mv $SUB_SCRIPT $SUB_SCRIPT_DIR
done < $IN_LIST
