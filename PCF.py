# PCF.py
#
# Author: Rick Field
#
# PolyCistronFinder (PCF) takes mapped Iso-seq reads and looks for evidence
# of polycistronic gene expression, i.e. where at least two predicted gene models
# are overlapped by an Iso-seq read. A common false positive result can arise
# from "split" gene models in the genome annotation. Therefore, PCF can
# optionally perform local or remote BLASTp alignments and analyze them to infer
# if the gene models likely represent subsequences of a "true" gene.
#
# Dependencies: pybedtools, pandas, biopython
#
# Example usage:
#
# python PCF.py -dp [path to data dir] -op [path to output dir] -bf [text file
# with bam file names] -bed --target_bed_file [target bed file] -f [int]
# -ol [float] -p [species prefix] --reference_pep_fastas [pep.fa]
# --reference_nuc_fastas [cds.fa] --name-field [e.g. 2 ':'] --evalue [float]
# --chrom_ignore [chrom ID] --local_blast --local_db [local.db]
# --make_local_db_fastas [other_species_pep.fa] --remote_blast
# --remote_db [remote.db] --slim_bed -ow -cb





import os
import sys
import argparse
import datetime
import subprocess
import shutil
import pybedtools
import pandas as pd
from itertools import combinations
from operator import itemgetter
from os.path import exists
from collections import Counter
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast import Record
from Bio.pairwise2 import format_alignment




def starter(bam_file):
    sample_id_dict = {}
    with open(bam_file, "r") as bam_list:
        new_bed_list = []
        for line in bam_list.readlines():
            line = line.strip().split()
            bam_name = line[0]
            read_fasta = line[1]
            sample_ID = line[2]
            print("Starting " + sample_ID + "...")
            sample_id_dict[sample_ID] = read_fasta
            new_bed_name = bam_to_bed(bam_name)
            new_bed_list.append((new_bed_name, sample_ID))

    return new_bed_list, sample_id_dict




def bam_to_bed(bam_name):

    new_bed_name = str(bam_name[:-4] + ".bed")
    if exists(args.data_path + new_bed_name):
        if args.overwrite:
            print(new_bed_name + " exists. Converting %s to bed format and overwrting..." % bam_name)
            bed = pybedtools.BedTool(args.data_path + bam_name).bam_to_bed(split=True)
            bed.saveas(outdir + "/" + new_bed_name)
        elif not args.overwrite:
            print(new_bed_name + " already exists. Skipping conversion...")

    else:
        print("Converting %s to bed format..." % bam_name)
        bed = pybedtools.BedTool(args.data_path + bam_name).bam_to_bed(split=True)
        bed.saveas(args.data_path + "/" + new_bed_name)

    return new_bed_name # String




def bed_maker(bed_list):

    print("Combining bed files...")
    combined_read_bed = []
    for i in bed_list:
        bed_file = args.data_path + i[0]
        sample_ID = str(i[1])
        with open(bed_file, "r") as bed:
            bed = bed.readlines()
            for row in bed:
                row = row.strip('\n').split()
                row.append(sample_ID)
                combined_read_bed.append('\t'.join(row))

    combined_read_bed = list(set(combined_read_bed)) # This step is here to filter duplicate entries in the read file. It was a problem with early Yucca Iso-seq read libraries from JGI.
    combined_read_bed_set = []
    for i in combined_read_bed:
        i = i.strip().split('\t')
        if args.chrom_ignore:
            if i[0] == args.chrom_ignore:
                pass
            else:
                combined_read_bed_set.append(i)
        else:
            combined_read_bed_set.append(i)

    combined_read_bed_dict = {} # Dictionary of reads and their bed coordinates
    read_length_dict = {} # Dictionary of read ids and their lengths
    for read in combined_read_bed_set:
        start = int(read[1])
        stop = int(read[2])
        id = read[3]
        strand = read[5]
        length = stop - start
        try:
            read_length_dict[id] += length
        except KeyError:
            read_length_dict[id] = length

        try:
            if strand == "+":
                combined_read_bed_dict[id].append(read)
                (combined_read_bed_dict[id]).sort(key=itemgetter(1)) # Re-sort list ascending start coordinates
            elif strand == "-":
                combined_read_bed_dict[id].append(read)
                (combined_read_bed_dict[id]).sort(key=itemgetter(1), reverse=True) # Re-sort list descending start coordinates
        except KeyError:
            combined_read_bed_dict[id] = [read]

    gene_bed_dict = {} # Dictinoary of genes and their .bed file lines
    if args.slim_bed:
        gene_bed_name = [i for i in args.target_bed.split(".")[:-1]]
        gene_bed_name.append("gene_line_only.bed")
        gene_bed_name = ".".join(gene_bed_name)
        with open(args.data_path + gene_bed_name, "w") as outfile:
            with open(args.data_path + args.target_bed, "r") as gene_targets_bed:
                gene_targets_bed = gene_targets_bed.readlines()
                for row in gene_targets_bed:
                    if row[0] == "#":
                        pass
                    else:
                        row = row.strip().split()
                        if row[7] == "mRNA":
                            gene = row[3]
                            gene_bed_dict[gene] = row
                            outfile.writelines("\t".join(row) + '\n')
                        else:
                            pass
    else:
        gene_bed_name = args.target_bed
        with open(args.data_path + gene_bed_name, "r") as gene_targets_bed:
            for row in gene_targets_bed:
                row = row.strip().split()
                if row[7] == "mRNA":
                    gene = row[3]
                    gene_bed_dict[gene] = row
                else:
                    pass

    if args.combined_read_bed_ouput:
        print("Writing combined reads bed file...")
        with open(ppl_out_dir + "/" + args.species_prefix + "_PPL_combined_reads_from_all_samples.bed", "w") as filehandle:
            for i in combined_read_bed_set:
                filehandle.write('\t'.join(i) + '\n')
    else:
        pass

    return combined_read_bed_set, combined_read_bed_dict, gene_bed_dict, read_length_dict, gene_bed_name # List, dictionary, dictionary




def intersect_finder(combined_read_bed_set, gene_bed_name): # This function finds the genes that intersect with each read.

    print("Intersecting reads with genes...")
    reads = pybedtools.BedTool(combined_read_bed_set)
    gene_annotation = pybedtools.BedTool(args.data_path + gene_bed_name)
    intersect = reads.intersect(gene_annotation, wo=True, s=True, bed=True, F=args.overlap, nonamecheck=True, split=True, stream="True")
    intersect_list = []
    for row in intersect: # Convert BedTool object to a list
        new_row = []
        for item in row:
            new_row.append(item)
        intersect_list.append(new_row)

    return intersect_list # List




def reads_to_multi_genes(intersect_file): # This function finds all reads that intersect more than one gene model and also makes a dictionary of the gene with the reads that map to them.

    print("Finding reads intersecting multiple gene models...")
    reads_to_genes_dict = {}
    gene_to_reads_dict = {}
    for i in intersect_file:
        read = i[3]
        gene = i[10]
        try:
            reads_to_genes_dict[read].append(gene)
        except KeyError:
            reads_to_genes_dict[read] = [gene]
        try:
            gene_to_reads_dict[gene].append(read)
        except KeyError:
            gene_to_reads_dict[gene] = [read]

    reads_to_multi_genes_dict = {}
    for read in reads_to_genes_dict:
        if len(set(reads_to_genes_dict[read])) > 1: # If more than one gene is intersected by the read, store in dict
            genes_list = [i for i in set(reads_to_genes_dict[read])]
            reads_to_multi_genes_dict[read] = genes_list
        else:
            pass

    return reads_to_multi_genes_dict, gene_to_reads_dict # Dictionary, dictionary




def ppl_curator(reads_to_multi_genes_dict, gene_to_reads_dict, gene_bed_dict, sample_id_dict):
    print("Viewing genes...")
    seen_genes = {} # Dictionary of genes and the genes they've been "seen" with in a given read
    for read in sorted(list(reads_to_multi_genes_dict.keys())):
        genes = reads_to_multi_genes_dict[read]
        for gene in genes:
            for i in genes:
                try:
                    if i not in seen_genes[gene]:
                        seen_genes[gene].append(i)
                    else:
                        pass
                except KeyError:
                    seen_genes[gene] = [i]


    young_ppl_dict = {} # Dictionary that comprehensively associates all other genes related by reads.
                        # e.g. read1: [gene 1, gene 2] + read2: [gene 2, gene 3] = [gene1, gene2, gene 3]
    for gene in sorted(list(seen_genes.keys())):
        young_ppl_dict[gene] = set(seen_genes[gene])
        for g in seen_genes:
            if young_ppl_dict[gene].intersection(set(seen_genes[g])): # Test if at least one gene is common to both sets
                young_ppl_dict[gene] = young_ppl_dict[gene].union(set(seen_genes[g]))
            else:
                pass

    count = 0
    gene_list = list(young_ppl_dict.keys())
    gene_list.sort()
    ### REMOVE DICT ELEMENTS THAT HAVE IDENTICAL VALUES, KEEPING ONLY THE FIRST ELEMENT ENCOUNTERED IN THE LIST ###
    while count < len(gene_list):
        gene = gene_list[count]
        for next_gene in gene_list[count + 1:]:
            try:
                if young_ppl_dict[gene] == young_ppl_dict[next_gene]:
                    del young_ppl_dict[next_gene]
                else:
                    pass
            except KeyError:
                pass
        count += 1

    print("Curating PPL...\n---")
    young_ppl_keys = list(young_ppl_dict.keys())
    young_ppl_count = len(young_ppl_keys)
    teen_ppl_dict = {}
    ppl_read_counts = {}
    count = 1
    for gene in young_ppl_keys:
        reads = gene_to_reads_dict[gene]
        for i in young_ppl_dict[gene]:
            reads = set(gene_to_reads_dict[i]).union(set(reads))
        # Exclude loci that do not pass the filter cutoff
        if not args.filter:
            filter = len(sample_id_dict)
        else:
            filter = args.filter
        if len(reads) >= filter:
            ppl = str(args.species_prefix + "_PPL" + str(count).zfill(5))
            teen_ppl_dict[ppl] = list(young_ppl_dict[gene])
            ppl_read_counts[ppl] = len(reads)
            count += 1
        else:
            pass

    mature_ppl_dict = {} # This is the main
    for ppl in teen_ppl_dict:
        sort_list = []
        first_gene = teen_ppl_dict[ppl][0]
        strand = gene_bed_dict[first_gene][5]
        if strand == "+":
            for gene in teen_ppl_dict[ppl]:
                sort_list.append([gene_bed_dict[gene][3], gene_bed_dict[gene][1]])
        elif strand == "-":
            for gene in teen_ppl_dict[ppl]:
                sort_list.append([gene_bed_dict[gene][3], gene_bed_dict[gene][2]])

        sort_list.sort(key=itemgetter(1))
        mature_ppl_gene_list_sorted = [i[0] for i in sort_list]
        mature_ppl_dict[ppl] = mature_ppl_gene_list_sorted
        mature_ppl_dict[ppl].append(strand)

    print("PPL CURATION SUMMARY\nPPL_ID\tgene_ID\ttotal_reads")
    for ppl in mature_ppl_dict:
        gene = ', '.join(mature_ppl_dict[ppl][:-1])
        line = ppl + '\t' + gene + '\t' + str(ppl_read_counts[ppl])
        print(line)

    print("---")
    print(str(len(young_ppl_dict)) + ' loci were found with at least one read covering at least two gene models ("young" PPL).')
    print(str(len(set(mature_ppl_dict))) + ' loci passed minimum coverage requirement (>' + str(args.filter) + ' reads, "mature" PPL).')
    print("---")

    return mature_ppl_dict, young_ppl_count # Dictionary where PPL ID is the key and the values are the coordinate-sorted gene IDs and the strand; integer count of PPL before filtering




def species_fastas(pep_fasta):

    if exists(fasta_dir):
        os.utime(fasta_dir)
        print(str(fasta_dir) + " exists. Overwriting fasta files in the directory...")
    else:
        os.mkdir(fasta_dir)

    pep_fasta_dict = {}
    new_ids = []
    with open(args.data_path + pep_fasta, "r") as infasta:
        for seq_record in SeqIO.parse(infasta, "fasta"):
            if args.name_field:
                record = str(seq_record.description)
                record = record.split(str(args.name_field[1]))
                id = record[int(args.name_field[0])].split("=")
                id = id[1]
                pep_fasta_dict[id] = str(seq_record.seq)
                new_ids.append([id, seq_record.id])
            else:
                pep_fasta_dict[seq_record.id] = str(seq_record.seq)

    if new_ids:
        with open(fasta_dir + "/" + args.species_prefix + "_peptide_fasta_name_conversions.txt", "w") as filehandle:
            for id in new_ids:
                filehandle.writelines('\t'.join(id))
    else:
        pass

    return pep_fasta_dict




def ppl_fastas(pep_fasta_dict, ppl):

    print("Writing fasta files for " + ppl + "...")
    strand = mature_ppl_dict[ppl][-1]
    with open(fasta_dir + "/" + ppl + ".faa", "w") as filehandle:
        if strand == "+":
            concat_pep_seq = ""
            for gene in mature_ppl_dict[ppl][:-1]:
                filehandle.writelines(">" + gene + '\n' + pep_fasta_dict[gene] + '\n')
                concat_pep_seq += pep_fasta_dict[gene].strip("*")
            filehandle.writelines(">" + ppl + "_concatenated_sequence\n" + concat_pep_seq + '*\n')
        elif strand == "-":
            concat_pep_seq = ""
            for gene in reversed(mature_ppl_dict[ppl][:-1]):
                filehandle.writelines(">" + gene + '\n' + pep_fasta_dict[gene] + '\n')
                concat_pep_seq += pep_fasta_dict[gene].strip("*")
            filehandle.writelines(">" + ppl + "_concatenated_sequence\n" + concat_pep_seq + '*\n')

    return str(ppl + ".faa")




def make_blast_db(protein_fasta):

    print("Building BLAST database from " + protein_fasta + "...")
    blast_db_name = str(protein_fasta + ".blastdb")
    cmd = str("makeblastdb -out " + blast_db_name + " -dbtype prot -in " + protein_fasta)
    print(cmd)
    subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return blast_db_name




def ppl_blast(ppl, blastdb):

    xml_out_path = os.path.join(blast_out_dir, ppl + "_results.xml")
    query = os.path.join(args.output_path, args.species_prefix + "_PPL_fasta_sequences", ppl + ".faa")
    if args.local_blast:
        cmd = str("blastp -out " + xml_out_path + " -outfmt 5 -query " + query + " -db " + blastdb + " -evalue " + str(args.evalue) + " -gapopen " + str(args.gap_open) + " -gapextend " + str(args.gap_extend))
    elif args.remote_blast:
        cmd = str("blastp -out " + xml_out_path + " -outfmt 5 -query " + query + " -db " + blastdb + " -evalue " + str(args.evalue) + " -gapopen " + str(args.gap_open) + " -gapextend " + str(args.gap_extend) + " -remote")

    print(cmd)
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if out.returncode != "0":
        pass
    else:
        print("BLASTp exit code: " + str(out.returncode), out.stderr)

    # TODO: Capture stdout and stderr
    aln_dict = {}
    with open(xml_out_path, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle) # Create a collection of BLAST records
        for blast_record in blast_records:
            aln_dict[blast_record.query] = {}
            for alignment in blast_record.alignments:
                # TODO: create argument for parsing hsp name
                subject_gene_id = ((alignment.title).split(" "))[1]
                aln_dict[blast_record.query][subject_gene_id] = {}
                for hsp in alignment.hsps:
                    aln_dict[blast_record.query][subject_gene_id] = [hsp.expect,
                        hsp.query_start,
                        hsp.query_end,
                        #hsp.match,
                        hsp.sbjct_start,
                        hsp.sbjct_end,
                        hsp.align_length,
                        alignment.length,
                        round(((hsp.sbjct_end - hsp.sbjct_start) / alignment.length)*100, 2)]

    return aln_dict # Dictionary of genes within a PPL and their alignments (HSPs)




def ppl_blast_results_check(aln_dict, ppl):

    print("Evaluating HSPs for each gene in " + ppl + "...")
    # Check if all genes returned BLAST HITS
    empty_count = 0
    genes_with_hsps = []
    empty_genes = []
    for gene in aln_dict:
        if "concatenated" not in str(gene):
            if not aln_dict[gene]:
                empty_count += 1
                empty_genes.append(gene)
            else:
                genes_with_hsps.append(gene)
        else:
            pass

    if empty_count == len(aln_dict) - 1:
        blast_result_status = False
    elif 0 < empty_count < len(aln_dict) - 1:
        blast_result_status = "PARTIAL"
    elif empty_count == 0:
        blast_result_status = "COMPLETE"

    return blast_result_status, genes_with_hsps, empty_genes




def ppl_blast_analyzer(aln_dict, ppl):

    pairwise_gene_combinations = [gene for gene in combinations(aln_dict, 2) if "concatenated" not in str(gene)]
    non_mutual_hsps = []
    mutual_hsp_gene_dict = {}
    for combo in pairwise_gene_combinations:
        subject_1_hsps = set(list(aln_dict[combo[0]].keys()))
        subject_2_hsps = set(list(aln_dict[combo[1]].keys()))
        mutuals = list(subject_1_hsps.intersection(subject_2_hsps))
        non_mutuals = list(subject_1_hsps.symmetric_difference(subject_2_hsps))
        for hsp in mutuals:
            mutual_hsp_gene_dict[hsp] = {}
            for gene in aln_dict:
                if hsp in aln_dict[gene].keys():
                    mutual_hsp_gene_dict[hsp][gene] = aln_dict[gene][hsp]
                else:
                    pass

        non_mutual_hsps.append(hsp for hsp in non_mutuals)

    # If there are mutual best hits, calculate best average escore to determine the best HSP to use.
    best_e_score = 1e-20
    best_hsp = "none"
    aln_sort = []
    genes_with_mutual_hsps = []
    if mutual_hsp_gene_dict:
        for hsp in mutual_hsp_gene_dict:
            e_score = 0
            for gene in mutual_hsp_gene_dict[hsp]:
                e_score += float(mutual_hsp_gene_dict[hsp][gene][0])

            avg_e_score = e_score / len(mutual_hsp_gene_dict[hsp])
            if avg_e_score < best_e_score:
                best_e_score = avg_e_score
                best_hsp = hsp
            else:
                pass

        sorted_alignments = []
        for gene in mutual_hsp_gene_dict[best_hsp]:
            if "concatenated" in gene:
                pass
            else:
                new_list = [gene]
                for i in mutual_hsp_gene_dict[best_hsp][gene]:
                    new_list.append(i)
                sorted_alignments.append(new_list)

            sorted_alignments.sort(key=itemgetter(4))

        for aln in sorted_alignments:
            genes_with_mutual_hsps.append(aln[0])
            query_start = aln[4]
            query_stop = aln[5]
            proportion = aln[-1]
            evalue = aln[1]
            aln_sort.append([ppl.strip('.fa'), best_hsp, aln[0], query_start, query_stop, proportion, evalue])

    elif non_mutual_hsps:
        pass

    aln_sort.sort(key=itemgetter(6))
    return aln_sort, best_hsp, best_e_score, genes_with_mutual_hsps




class PPL:

    def __init__(self, ID):
        self.ID = ID

    def ppl_viewer(self, ID):
        self.genes = [gene for gene in mature_ppl_dict[self.ID][:-1]]
        self.gene_count = len(mature_ppl_dict[self.ID][:-1]) # Total genes in PPL
        self.gene_bed_lines = [gene_bed_dict[gene] for gene in mature_ppl_dict[self.ID][:-1]] # .bed lines for each gene in PPL
        reads = set()
        for gene in self.genes:
            reads = reads.union(set(gene_to_reads_dict[gene]))
        self.reads = [read for read in list(reads)]

    def ppl_blaster(self, ID):
        self.pep_fasta_file = ppl_fastas(pep_fastas, self.ID)
        self.aln_dict = ppl_blast(self.ID, db)
        self.blast_result_status, self.genes_with_hsps, self.empty_genes = ppl_blast_results_check(self.aln_dict, self.ID)
        self.aln_sort, self.best_hsp, self.best_avg_escore, self.genes_with_mutual_hsps = ppl_blast_analyzer(self.aln_dict, self.ID)




def get_parser():
    parser = argparse.ArgumentParser(description='Find putative polycistronic loci using Iso-seq reads.')
    parser.add_argument("-dp", "--data_path", help="Path to data directory. *MUST EXIST*", required=True, dest="data_path")
    parser.add_argument('-op', "--output_path", help="Path to output directory. *MUST EXIST*", required=True, dest="output_path")
    parser.add_argument('-bf', '--bam_file',
        help="Text file with input bam files (first column) and sample prefix (second column). bam file must be sorted and filtered by samtools. *MUST EXIST*",
        required=True, dest="bam_file")
    parser.add_argument("-bed", "--target_bed_file", help="e.g. species.gene.bed. Must reside in data path. *MUST EXIST*", required=True, dest="target_bed")
    parser.add_argument("-f", "--filter", help="Exclude genes in PPL with less than this number of total reads mapped to the locus. If not provided, value of one read per sample average.", required=False, dest="filter", type=int)
    parser.add_argument("-ol", "--read_overlap_percent", help="Minimum overlap required as a fraction of the read length (type=float). *MUST EXIST*", required=True, dest="overlap", type=float)
    parser.add_argument("-p", "--prefix", help="Prefix for output files, probably species or treatment. *MUST EXIST*", required=True, dest="species_prefix")
    parser.add_argument("-ow", "--overwrite_existing_beds", help="If used, PCF will convert .bam files to .bed format and overwrite existing .bed files", required=False,
        dest="overwrite", action="store_true")
    parser.add_argument("-cb", "--combined_read_bed_output", help="If used, PCF will output the combined bed file for all reads used as input.",
        required=False, dest="combined_read_bed_ouput", action="store_true")
    parser.add_argument("--reference_pep_fastas", help="Fasta file for all genes in reference annotation. Must reside in data path. If provided, fasta sequences for each gene within each PPL will be written to a text file.", required=False, dest="ref_pep_fastas")
    parser.add_argument("--reference_nuc_fastas", help="Fasta file for all genes in reference annotation. Must reside in data path. If provided, fasta sequences for each gene within each PPL will be written to a text file.", required=False, dest="ref_nuc_fastas")
    parser.add_argument("--name-field", help="When extracting fasta seqeunces, use the nth ID as sequence name. First argument is the field to use (zero-based index), second argument is the separator in quaotation (i.e. '\t', '.', '|', etc.).",
        dest="name_field", required=False, default=False, nargs=2)
    parser.add_argument("--local_blast", help="Use local NCBI BLAST to align sequences to a protein database.", required=False, action="store_true")
    parser.add_argument("--local_db", help="Name of local BLAST database to use. Must include full path or have a shell PATH variable.", required=False, dest="local_db")
    parser.add_argument("--make_local_db_fastas", help="Fasta file containing peptide sequences from which a BLAST db can be made.", required=False, dest="make_local_db_fastas")
    parser.add_argument("--remote_blast", help="Use NCBI remote BLAST function.", required=False, dest="remote_blast", action="store_true")
    parser.add_argument("--remote_db", help="Database to use with remote BLAST.", required=False, dest="remote_db", default="swissprot")
    parser.add_argument("--evalue", help="Expect value cutoff for BLAST.", required=False, default=1e-20, dest="evalue", type=float)
    parser.add_argument("--gap_open", help="Gap open penalty for BLASTp.", required=False, dest="gap_open", type=int)
    parser.add_argument("--gap_extend", help="Gap extend penalty for BLASTp.", required=False, dest="gap_extend", type=int)
    parser.add_argument("--chrom_ignore", help="Chromosome ID to ignore during intersecting.", required=False, dest="chrom_ignore")
    parser.add_argument("--slim_bed", help="Make a target bed file with only 'gene' features.", required=False, dest="slim_bed")
    args = vars(parser.parse_args())
    return parser



"""RUN SECTION"""
if __name__ == "__main__":

    start_time = datetime.datetime.now()
    print("---\nPCF.py started at " + str(start_time))
    print(' '.join(sys.argv[0:]))
    print("---")

    args = get_parser().parse_args()
    ppl_out_dir = os.path.join(args.output_path, args.species_prefix + "_PPL_output")
    if exists(ppl_out_dir):
        shutil.rmtree(ppl_out_dir)
        os.mkdir(ppl_out_dir)
        print(str(ppl_out_dir) + " exists. Existing files will be deleted...")
    else:
        os.mkdir(ppl_out_dir)

    fasta_dir = os.path.join(args.output_path, args.species_prefix + "_PPL_fasta_sequences")
    if exists(fasta_dir):
        shutil.rmtree(fasta_dir)
        os.mkdir(fasta_dir)
        print(str(fasta_dir) + " exists. Existing files will be deleted...")
    else:
        os.mkdir(fasta_dir)

    new_bed_list, sample_id_dict = starter(args.bam_file)
    combined_read_bed_set, combined_read_bed_dict, gene_bed_dict, read_length_dict, gene_bed_name = bed_maker(new_bed_list)
    intersects = intersect_finder(combined_read_bed_set, gene_bed_name)
    reads_to_multi_genes_dict, gene_to_reads_dict = reads_to_multi_genes(intersects)
    mature_ppl_dict, young_ppl_count = ppl_curator(reads_to_multi_genes_dict, gene_to_reads_dict, gene_bed_dict, sample_id_dict)

    if args.ref_pep_fastas: # Extract peptide sequences for all genes in all PPL, fasta format
        pep_fastas = species_fastas(args.ref_pep_fastas)

    if args.local_blast: # Run local BLAST
        if args.make_local_db_fastas: # Make a BLAST db from peptide sequences
            db = make_blast_db(args.make_local_db_fastas)
        elif args.local_db: # Use a local BLAST db
            db = args.local_db

        blast_out_dir = os.path.join(args.output_path, args.species_prefix + "_PPL_BLAST_output")
        if exists(blast_out_dir):
            os.utime(blast_out_dir)
            print(str(blast_out_dir) + " exists. Existing files will be overwritten...")
        else:
            os.mkdir(blast_out_dir)

        print("Starting local BLASTp searches with " + db + " databse for each mature PPL...")

    elif args.remote_blast: # Run remote BLAST
        db = args.remote_db
        blast_out_dir = os.path.join(args.output_path, args.species_prefix + "_PPL_BLAST_output")
        if exists(blast_out_dir):
            os.utime(blast_out_dir)
            print(str(blast_out_dir) + " exists. Existing files will be overwritten...")
        else:
            os.mkdir(blast_out_dir)

        print("Starting remote BLASTp searches with " + db + " databse for each mature PPL...")




    """PPL CLASS"""

    ppl_summary_dict = {}

    for ppl in mature_ppl_dict:
        ppl = PPL(ppl)
        ppl.ppl_viewer(ppl)
        print("--- " + ppl.ID + " ---")

        if args.local_blast or args.remote_blast:
            ppl.ppl_blaster(ppl)

        if ppl.blast_result_status == "COMPLETE":
            print("COMPLETE BLAST result: " + ppl.ID + " returned BLAST hits for every gene (" + str(len(ppl.genes)) + ").")
            if ppl.aln_sort:
                print(", ".join(ppl.genes_with_mutual_hsps) + " returned mutual HSPs.")
                print("The best HSP aligned to: " + ppl.best_hsp + " with avergae e-score = " + str(ppl.best_avg_escore))
                with open(ppl_out_dir + "/" + args.species_prefix + "_COMPLETE_blast_mutual_hsps.txt", "a") as filehandle:
                    #genes = [i for i in ppl.genes]
                    filehandle.writelines(ppl.ID + '\t' + ", ".join(ppl.genes) + '\n')
            else:
                print("No genes in " + ppl.ID + " returned mutual HSPs.")
                with open(ppl_out_dir + "/" + args.species_prefix + "_COMPLETE_blast_NO_mutual_hsps.txt", "a") as filehandle:
                    #genes = [i for i in ppl.genes]
                    filehandle.writelines(ppl.ID + '\t' + ", ".join(ppl.genes) + '\n')

        elif ppl.blast_result_status == "PARTIAL":
            print("PARTIAL BLAST result: " + ppl.ID + " contains genes that did not return any BLAST hits: " + ", ".join(ppl.empty_genes) + ".")
            if ppl.aln_sort:
                print(", ".join(ppl.genes_with_mutual_hsps) + " returned mutual HSPs.")
                print("The best HSP aligned to: " + ppl.best_hsp + " with avergae e-score = " + str(ppl.best_avg_escore))
                with open(ppl_out_dir + "/" + args.species_prefix + "_PARTIAL_blast_mutual_hsps.txt", "a") as filehandle:
                    #genes = [i for i in ppl.genes]
                    filehandle.writelines(ppl.ID + '\t' + ", ".join(ppl.genes) + '\n')
            else:
                print("No genes in " + ppl.ID + " returned mutual HSPs.")
                with open(ppl_out_dir + "/" + args.species_prefix + "_PARTIAL_blast_NO_mutual_hsps.txt", "a") as filehandle:
                    #genes = [i for i in ppl.genes]
                    filehandle.writelines(ppl.ID + '\t' + ", ".join(ppl.genes) + '\n')

        elif not ppl.blast_result_status:
            print("NO BLAST RESULTS: " + ppl.ID + " did not return any BLAST hits.")
            with open(ppl_out_dir + "/" + args.species_prefix + "_NO_blast_results.txt", "a") as filehandle:
                #genes = [i for i in ppl.genes]
                filehandle.writelines(ppl.ID + '\t' + ", ".join(ppl.genes) + '\n')

        ppl_summary_dict[ppl.ID] = {'genes': ppl.genes,
                                    'gene_count': ppl.gene_count,
                                    'genes_with_hsps': ppl.genes_with_hsps,
                                    'genes_with_mutual_hsps': ppl.genes_with_mutual_hsps,
                                    'best_hsp': ppl.best_hsp,
                                    'best_avg_escore': ppl.best_avg_escore,
                                    'best_hsp_alignments': ppl.aln_sort}


    """OUTPUT SECTION"""

    os.utime(args.output_path)
    with open(ppl_out_dir + "/" + args.species_prefix + "_PPL.bed", "w") as filehandle:
        for ppl in mature_ppl_dict:
            for gene in mature_ppl_dict[ppl][:-1]:
                row = gene_bed_dict[gene]
                row.append(ppl)
                filehandle.writelines('\t'.join(row) + '\n')

    with open(ppl_out_dir + "/" + args.species_prefix + "_PPL_summary.txt", "w") as filehandle:
        gene_count_dict = {}
        for ppl in mature_ppl_dict:
            try:
                gene_count_dict[ppl_summary_dict[ppl]['gene_count']] += 1
            except KeyError:
                gene_count_dict[ppl_summary_dict[ppl]['gene_count']] = 1

        filehandle.writelines("---\nPCF.py run at " + str(start_time) + '\n')
        filehandle.writelines(' '.join(sys.argv[0:]) + '\n')
        filehandle.writelines("---\n")
        filehandle.writelines(str(young_ppl_count) + ' loci were found with at least one read covering at least two gene models ("young" PPL).\n')
        filehandle.writelines(str(len(mature_ppl_dict)) + ' loci passed minimum coverage requirement (>' + str(args.filter) + ' reads, "mature" PPL).\n')
        filehandle.writelines("---\nclass\tcount\tgenes\n")
        total_genes = 0
        for i in sorted(gene_count_dict.keys()):
            filehandle.writelines(str(i) + '\t' + str(gene_count_dict[i]) + '\t' + str(i*gene_count_dict[i]) + '\n')
            total_genes += (i*gene_count_dict[i])

        filehandle.writelines("\t\t---\n\t\t" + str(total_genes))
        filehandle.writelines("\n---\nsample_id\tread_file\n")
        for i in sample_id_dict:
            filehandle.writelines(i + '\t' + sample_id_dict[i] + '\n')

        filehandle.writelines("---\n")
        filehandle.writelines('PPL with mutual HSPs for at least 2 genes\n')
        filehandle.writelines('\t'.join(["PPL_ID", 'subject_gene_ID', 'query_gene_ID', 'query_start', 'query_stop', 'query_cov_%', 'evalue']) + '\n')
        for ppl in ppl_summary_dict:
            if ppl_summary_dict[ppl]['genes_with_mutual_hsps']:
                for alignment in ppl_summary_dict[ppl]['best_hsp_alignments']:
                    alignment = [str(i) for i in alignment]
                    filehandle.writelines('\t'.join(alignment) + '\n')


    ### SORTING BED FILE MATURE PPL READS AND OUTPUTTING ###
    outlist = []
    for ppl in mature_ppl_dict:
        for gene in mature_ppl_dict[ppl][:-1]:
            for reads in gene_to_reads_dict[gene]:
                for i in combined_read_bed_dict[reads]:
                    line = ('\t'.join(i)).strip('\n') + '\t' + ppl + '\n'
                    outlist.append(line)

    outlist = set(outlist)
    cols = ["chr", "start", "stop", "read", "qc", "strand", "sample", "ppl"]
    df = pd.DataFrame(columns = cols)
    index = 0
    for row in outlist:
        row = row.strip().split('\t')
        new_df = pd.DataFrame(
            {"chr":row[0],
            "start":int(row[1]),
            "stop":row[2],
            "read":row[3],
            "qc":row[4],
            "strand":row[5],
            "sample":row[6],
            "ppl":row[7]}, index=[index])
        df = pd.concat([new_df, df], ignore_index=True)
        index += 1
    df.sort_values(by=["chr", "start", "read", "ppl"], ascending=True, kind="stable", inplace=True)#, ignore_index=True)
    df.to_csv(ppl_out_dir + "/" + args.species_prefix + "_PPL_mature_reads.bed", sep='\t', header=False, index=False)

    stop_time = datetime.datetime.now()
    print("---\nPCF.py stopped at " + str(stop_time))
    run_time = stop_time - start_time
    print("Total run time = " + str(run_time))
    print("---")
