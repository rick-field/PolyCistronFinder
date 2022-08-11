#!/bin/zsh

# Author: Rick Field
# Email: richard.field@uga.edu
# Purpose: to submit many minimap2 jobs

INFILE=$1

if [ -e './minimap2.out' ]
then
  rm './minimap2.out'
fi

while read species prefix hap dir genome read_list
do
  {
    while read read_file sample_ID
      do
        {
          echo '#!/bin/bash'
          echo '#SBATCH --job-name='$prefix'_'$sample_ID'minimap2'
          echo '#SBATCH --partition=batch'
          echo '#SBATCH --ntasks=1'
          echo '#SBATCH --nodes=1'
          echo '#SBATCH --cpus-per-task=8'
          echo '#SBATCH --mem=50gb'
          echo '#SBATCH --time=24:00:00'
          echo '#SBATCH --output=minimap2.'$prefix'_'$sample_ID'.out'
          echo '#SBATCH --error=minimap2.'$prefix'_'$sample_ID'.err'
          echo ''
          echo 'cd $SLURM_SUBMIT_DIR'
          echo ''
          echo 'module load minimap2/2.17-GCC-8.3.0'
          echo ''
          echo 'minimap2 -ax splice:hq -t 8 -uf'
          echo $dir'data/'$hap'/'$genome' \'
          echo $dir'data/reads/'$read_file' \'
          echo '> '$dir'analyses/'$hap'/'$prefix'_'$sample_ID'_minimap2_aln.sam'
        } > $prefix'_'$sample_ID'_minimap2.sh'
        echo -e $species'\t'$prefix'\t'$hap'\t'$dir'\t'$sample_ID'\t'$dir'analyses/'$hap'/'$prefix'_'$sample_ID'_minimap2_aln.sam' >> './minimap2.out'
        sbatch $prefix'_'$sample_ID'_minimap2.sh'
        echo $prefix'_'$sample_ID'_minimap2.sh Submitted!'
      done < $dir'/data/reads/'$read_list
  }
done < $INFILE
