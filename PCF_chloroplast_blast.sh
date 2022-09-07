#!/bin/bash
#SBATCH --job-name=j_BLAST+_multithread
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --time=12:00:00
#SBATCH --output=BLAST+.%j.out
#SBATCH --error=BLAST+.%j.err

cd $SLURM_SUBMIT_DIR

ml BLAST+/2.12.0-gompi-2020b

blastx -num_threads 8 -query /scratch/rdf85141/PPL/plastid_seqs/cp_mt_renamed.fa -out acor_cp_mito_blast.out \
-outfmt 6 -evalue 1e-20 -parse_deflines \
-db /home/rdf85141/genomes/Acor/Aamericanus_586_v1.1.protein_primaryTranscriptOnly.fa
