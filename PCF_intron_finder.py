# Author: Rick Field
#
# Purpose: to count introns in each gene in a typical 'gene_exons.gff3'
# formatted file. Prints to std_out. 

import sys



infile = str(sys.argv[1])

with open(infile, "r") as filehandle:
    map_file = {}
    gene_dict = {}
    exon_dict = {}
    for row in filehandle.readlines():
        if row[0] == "#":
            pass
        else:
            row = row.strip().split('\t')
            if row[2] == "gene":
                fields = row[-1].split(";")
                IDs = fields[0].split("=")
                ID = IDs[1]
                gene_dict[ID] = []

            elif row[2] == "mRNA":
                fields = row[-1].split(";")
                longest = fields[2].split("=")
                if longest[1] == "1":
                    parents = fields[-1].split("=")
                    parent = parents[1]
                    IDs = fields[0].split("=")
       	            mRNA_ID = IDs[1]
                    map_file[parent] = mRNA_ID
                else:
                    pass
            elif row[2] == "exon":
                fields = row[-1].split(";")
                IDs = fields[0].split("=")
                full_ID = IDs[1]
                ID = IDs[1].split(".")
                ID = ID[0]
                try:
                    exon_dict[ID].append([full_ID, row[3], row[4]])
                except KeyError:
                    exon_dict[ID] = []
                    exon_dict[ID].append([full_ID, row[3], row[4]])
            else:
                pass 

    for gene in gene_dict:
        mRNA_ID = map_file[gene]
        exons = exon_dict[mRNA_ID]
        exon_count = len(exons) - 1
        exon_list = []
        for exon in exons:
            ex = ", ".join(exon)
            exon_list.append(ex)
        exons = "\t".join(exon_list)
        print(gene + '\t' + str(exon_count) + '\t' + exons)
