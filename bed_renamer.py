import sys

bed = str(sys.argv[1])
bed = (open(bed, "r")).readlines()

for i in bed:
	i = i.strip().split('\t')
	i[3] = i[3].strip(".g")
	print('\t'.join(i))	
