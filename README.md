PolyCistronFinder (PCF)

Version: 1.0
Author: Rick Field
Email: richard.field@uga.edu
Advisor: Jim Leebens-Mack

Purpose: PolyCistronFinder takes mapped Iso-seq reads and looks for evidence of
polycistronic gene expression, i.e. where at least two predicted gene models are 
overlapped by an Iso-seq read. A common false positive result can arise from "split" 
gene models in the genome annotation. Therefore, PCF can optionally perform local 
or remote BLASTp and analyze the alignments to determine if the gene models likely 
represent subsequences of a "true" gene.

Generalized workflow:
1. Map reads to genomes (minimap2.sh)
2. Sort and filter bam files (samtools_sort.sh, samtools_view_256.sh, samtools_view_2048.sh)
3. Prepare gff3 files for PCF.py (PCF_gff3_prep.sh)
4. Prepare input file for bacth PCF.py runs (PCF_infile_maker.py)
5. Submit multiple PCF.py jobs (PCF.sh)
