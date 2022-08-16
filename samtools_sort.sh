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
if [ -! $SUB_SCRIPT_DIR ]
then
  mkdir	$SUB_SCRIPT_DIR
fi

while read species prefix hap dir sample_id in_file
do
  OUT_FILE=${in_file%.*}'_sorted.bam'
  OUT_DIR=$dir'data/'$hap'/'
  SUB_SCRIPT=$prefix'_'$hap'_'$sample_id'_samtools_sort.sh'
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
    echo '-o '$OUT_FILE' \'
    echo $in_file
  } > $SUB_SCRIPT
  echo -e $species'\t'$prefix'\t'$hap'\t'$dir'\t'$sample_id'\t'$OUT_FILE >> $OUT_LIST
  sbatch $SUB_SCRIPT
  echo $SUB_SCRIPT' Submitted!'
  mv $SUB_SCRIPT $SUB_SCRIPT_DIR
done < $IN_LIST
