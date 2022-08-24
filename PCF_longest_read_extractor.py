import sys
from Bio import SeqIO

inlist = (open(str(sys.argv[1]), "r")).readlines()
read_fasta_list = (open(str(sys.argv[2]), "r")).readlines()

read_ids = []
for row in inlist:
    row = row.strip().split()
    read_ids.append(row[2])

for sample in read_fasta_list:
    sample = sample.strip()
    with open(sample, "r") as read_fasta:
        for seq_record in SeqIO.parse(read_fasta, "fasta"):
            if seq_record.id in read_ids:
                print(">" + str(seq_record.id))
                print(str(seq_record.seq))
            else:
                pass



