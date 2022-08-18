#!/bin/python/
#
# Author: Rick Field
# Email: richard.field@uga.edu
#
# PCF_infile_maker.py produces a file for using when submitting multiple PCF.py jobs. 
# The first input file should have rows with the following fields: species, prefix, haplotype,  
# data directory, sample ID, read file, bam file.
# The second file should have sample ID, gene model annotation (bed format), and fasta file  
# of protein seqeunces for your species.

import sys

bam_list = str(sys.argv[1])
bed_list = str(sys.argv[2])
bam_dict = {}
bed_dict = {}

with open(bam_list, "r") as bams:
    bams = bams.readlines()
    for row in bams:
        row = row.strip().split()
        prefix = row[1]
        try:
            bam_dict[prefix].append(row)
        except KeyError:
            bam_dict[prefix] = [row]

with open(bed_list, "r") as beds:
    beds = beds.readlines()
    for row in beds:
        row = row.strip().split()
        prefix = row[0]
        bed_dict[prefix] = [row[1], row[2]]

for prefix in bam_dict:
    basics = '\t'.join(bam_dict[prefix][0][0:4])
    out_file = prefix + "_" + str(bam_dict[prefix][0][2]) + '_bam_list.txt'
    with open(out_file, "w") as filehandle:
        for i in bam_dict[prefix]:
            bam = i[-1].split('/')
            filehandle.writelines(bam[-1] + '\t' + i[5] + '\t' + i[4] + '\n')
        
        beds = '\t'.join(bed_dict[prefix])
        with open("PCF_samples.txt", "a") as filehandle:
            filehandle.writelines(basics + '\t' + beds + '\t' + out_file + '\n')
