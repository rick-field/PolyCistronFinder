import sys

infile = str(sys.argv[1])
feature = str(sys.argv[2])

with open(infile, "r") as in_bed:
    if feature == 'mRNA':
        for row in in_bed.readlines():
            row = row.strip().split('\t')
            ids = row[-1].split(";")
            name = ids[1].split("=")
            name = name[1]
            row[3] = name
            print('\t'.join(row))

    elif feature == 'gene':
        for row in in_bed.readlines():
            row = row.strip().split('\t')
            i[3] = i[3].strip(".g")
            print('\t'.join(i))
