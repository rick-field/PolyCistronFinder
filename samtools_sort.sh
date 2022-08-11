#!/bin/zsh

# Author: Rick Field
# Email: richard.field@uga.edu
# Purpose: to submit many samtools sort jobs

INFILE=$1

if [ -e './samtools_sort.out' ]
then
  rm './samtools_sort.out'
fi

while read SPECIES PREFIX HAP DIR SAMPLE_ID SAM
SORTED_SAM=$DIR$HAP'analyses/'$SAM
do
  {
    echo '#!/bin/bash'
    echo '#SBATCH --job-name='$PREFIX'_'$SAMPLE_ID'samtools_sort'
    echo '#SBATCH --partition=batch'
    echo '#SBATCH --ntasks=1'
    echo '#SBATCH --nodes=1'
    echo '#SBATCH --cpus-per-task=8'
    echo '#SBATCH --mem=50gb'
    echo '#SBATCH --time=24:00:00'
    echo '#SBATCH --output=SAMtools_sort.'$PREFIX'_'$SAMPLE_ID'.out'
    echo '#SBATCH --error=SAMtools_sort.'$PREFIX'_'$SAMPLE_ID'.err'
    echo ''
    echo 'cd $SLURM_SUBMIT_DIR'
    echo ''
    echo 'ml samtools/1.10-GCC-8.3.0'
    echo ''
    echo 'samtools sort -@ 4 \'
    echo '-O bam \'
    echo '-o '$SORTED_SAM' \'
    echo $DIR$HAP'analyses/'$SAM
  } > $PREFIX'_'$SAMPLE_ID'_samtools_sort.sh'
  echo -e $SPECIES'\t'$PREFIX'\t'$HAP'\t'$DIR'\t'$DIR'analyses/'$HAP'/'$PREFIX'_'$SAMPLE_ID'_samtools_sort_aln.sam' >> './samtools_sort.out'
  sbatch $PREFIX'_'$SAMPLE_ID'_samtools_sort.sh'
  echo $PREFIX'_'$SAMPLE_ID'_samtools_sort.sh Submitted!'
done < $INFILE
