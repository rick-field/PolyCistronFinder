## PolyCistronFinder  

**PolyCistronFinder (PCF)** takes mapped Iso-seq reads and looks for evidence of  
polycistronic gene expression, i.e. where at least two predicted gene models are  
overlapped by an Iso-seq read. A common false positive result can arise from "split"  
gene models in the genome annotation. Therefore, PCF can optionally perform local  
or remote BLASTp and analyze the alignments to determine if the gene models likely  
represent subsequences of a "true" gene.  

**Generalized workflow:**  
1. Map reads to genomes  
2. Sort and filter bam files  
3. Prepare gff3 files for PCF.py  
4. Prepare input file for bacth PCF.py runs  
5. Submit multiple PCF.py jobs  


**Usage:**  

```
python PCF.py [options] > <species.out>  
```

**Options:**

```  
-dp --data_path [str]                     Path to data directory. Required.  
-op --output_path [str]                   Path to output directory. Required.  
-bf --bam_file [text file]                Text file with input bam files (first column) and sample prefix (second column). bam file must be sorted and filtered by samtools. Required.  
-bed --target_bed_file [bed file]         Bed file containg gene model predictions in bed format. Must reside in data path. Required.  
-f --filter [int]                         Exclude genes in PPL with less than this number of total reads mapped to the locus. If not provided, value of one read per sample average.  
-ol --read_overlap_percent [float]        Minimum overlap required as a fraction of the read length. Required.  
-p --prefix [str]                         Prefix for output files, probably species or treatment. Required.  
-ow --overwrite_existing_beds [str]       If used, PCF will convert .bam files to .bed format and overwrite existing .bed files.  
-cb --combined_read_bed_output [str]      If used, PCF will output the combined bed file for all reads used as input.  
--reference_pep_fastas [fasta file]       Fasta file for all genes in reference annotation. Must reside in data path. If provided, fasta sequences for each gene within each PPL will be written to a text file.  
--reference_nuc_fastas [fasta file]       Fasta file for all genes in reference annotation. Must reside in data path. If provided, fasta sequences for each gene within each PPL will be written to a text file.  
--name-field [str, int]                   When extracting fasta seqeunces, use the nth ID as sequence name. First argument is the field to use (zero-based index), second argument is the separator in quaotation (i.e. '\t', '.', '|', etc.).  
--local_blast                             Use local NCBI BLAST to align sequences to a protein database.  
--local_db                                Name of local BLAST database to use. Must include full path or have a shell PATH variable.  
--make_local_db_fastas *fasta file*       Fasta file containing peptide sequences from which a BLAST db can be made.  
--remote_blast                            Use NCBI remote BLAST function.  
--remote_db [str]                         Database to use with remote BLAST. default = swissprot  
--evalue [float]                          Expect value cutoff for BLAST. default = 1e-20  
--chrom_ignore [str]                      Chromosome ID to ignore during intersecting.  
```
