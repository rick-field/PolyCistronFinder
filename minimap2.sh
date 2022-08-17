#!/bin/bash

# Author: Rick Field
# Email: richard.field@uga.edu
# Purpose: to submit many minimap2 jobs

INFILE=$1

OUT_FILE="./minimap2_out.txt"
if [ -e $OUT_FILE ]
then
  rm $OUT_FILE
fi

SUB_SCRIPT_DIR=$PWD'/submission_scripts'
if [ ! -d $SUB_SCRIPT_DIR ]
then
  mkdir -p $SUB_SCRIPT_DIR
fi
 
while read species prefix hap dir genome read_list
do
  {
    while read read_file sample_ID
      do
        OUT_DIR=$dir'data/'$hap'/'
        SAM=$prefix'_'$hap'_'$sample_ID'_aln.sam'
        SUB_SCRIPT=$prefix'_'$hap'_'$sample_ID'_minimap2.sh'
        {
          echo '#!/bin/bash'
          echo '#SBATCH --job-name='${SUB_SCRIPT%.*}
          echo '#SBATCH --partition=batch'
          echo '#SBATCH --ntasks=1'
          echo '#SBATCH --nodes=1'
          echo '#SBATCH --cpus-per-task=8'
          echo '#SBATCH --mem=50gb'
          echo '#SBATCH --time=36:00:00'
          echo '#SBATCH --output='${SUB_SCRIPT%.*}'.out'
          echo '#SBATCH --error='${SUB_SCRIPT%.*}'.err'
          echo ''
          echo 'cd $SLURM_SUBMIT_DIR'
          echo ''
          echo 'module load minimap2/2.17-GCC-8.3.0'
          echo ''
          echo 'minimap2 -ax splice:hq -t 8 -uf \'
          echo $dir'data/'$hap'/'$genome' \'
          echo $dir'data/reads/'$read_file' \'
          echo '> '$OUT_DIR$SAM
        } > $SUB_SCRIPT
        echo -e $species'\t'$prefix'\t'$hap'\t'$dir'\t'$sample_ID'\t'$read_file'\t'$OUT_DIR$SAM >> $OUT_FILE
        sbatch $SUB_SCRIPT
        echo $SUB_SCRIPT' Submitted!'
        mv $SUB_SCRIPT $SUB_SCRIPT_DIR
      done < $dir'/data/reads/'$read_list
  }
done < $INFILE

