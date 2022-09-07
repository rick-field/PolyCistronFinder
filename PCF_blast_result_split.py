import sys

infile = str(sys.argv[1])

with open(infile, "r") as id_file:
    for row in id_file.readlines():
        row = row.strip().split()
        genes = row[1:]#.split(", ")
        for gene in genes:
            print(gene.strip(","))
