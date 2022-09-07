import sys

ppl_file = str(sys.argv[1])
blast_file = str(sys.argv[2])


blast_dict = {}
with open(blast_file, "r") as blast_file:
    for row in blast_file.readlines():
        row = row.strip().split()
        plastid_gene = row[0]
        species_gene = row[1]        
        evalue = row[-2]
        try:
            if evalue < blast_dict[species_gene][1]:
                blast_dict[species_gene] = [plastid_gene, evalue]
            else:
                pass
        except KeyError:
            blast_dict[species_gene] = [plastid_gene, evalue]

with open(ppl_file, "r") as ppl_file:
    for row in ppl_file:
        row = row.strip().split()
        id = row[3]
        try:
            print('\t'.join(row) + '\t' + blast_dict[id][0] + '\t' + blast_dict[id][1])
        except KeyError:
            print('\t'.join(row) + '\tnot_plastid\tNA')
