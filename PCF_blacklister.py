import sys



def blacklister(aln_file):
    blacklist_dict = {}
    with open(aln_file, "r") as aln:
        aln = aln.readlines()
        query_id = aln[0][1:].strip()
        for row in aln[1:]:
            if row[0] == ">":
                pass
            else:
                row = row.strip().split()
                chrom = row[0]
                start = int(row[1])
                length = int(row[3])
                stop = start + length
                try:
                    blacklist_dict[chrom].append([start, stop])
                except KeyError:
                    blacklist_dict[chrom] = [[start, stop]]

    return blacklist_dict, query_id



def blacklist_filter(blacklist_dict, ppl_file):
    ppl_dict = {}
    removed_ppl_dict = {}
    with open(ppl_file, "r") as ppl_bed:
        for row in ppl_bed:
            row = row.strip().split()
            ppl = row[-1]
            try:
                ppl_dict[ppl].append(row)
            except KeyError:
                ppl_dict[ppl] = [row]
    
    ppl_ids = [i for i in ppl_dict]
    for ppl in ppl_ids:
        for gene in ppl_dict[ppl]:
            chrom = gene[0]
            start = int(gene[1])
            stop = int(gene[2])
            ppl_interval = set(list(range(start, stop)))
            try:
                for interval in blacklist_dict[chrom]:
                    black_interval = set(list(range(interval[0], interval[1])))
                    if len(ppl_interval.intersection(black_interval)) > 1:
                        del ppl_dict[ppl]
                        try:
                            removed_ppl_dict[ppl].append(gene)
                        except KeyError:
                            removed_ppl_dict[ppl] = [gene]
                    else:
                        pass
            except KeyError:
                pass

    return ppl_dict, removed_ppl_dict



aln_file = str(sys.argv[1])
ppl_file = str(sys.argv[2])
prefix = str(sys.argv[3])
outdir = str(sys.argv[4])

blacklist_regions, query_id = blacklister(aln_file)
final_ppl, removed_ppl = blacklist_filter(blacklist_regions, ppl_file)

for ppl in final_ppl:
    for gene in final_ppl[ppl]:
        print('\t'.join(gene))

with open(outdir + prefix + "_PPL_blacklisted.bed", "w") as filehandle:
    for ppl in removed_ppl:
        for gene in removed_ppl[ppl]:
            print('\t'.join(gene))

with open(outdir + prefix + "_PCF_blacklisted_regions_from_" + query_id + ".bed", "w") as filehandle:
    for region in blacklist_regions:
        for interval in blacklist_regions[region]:
            line = [region, str(interval[0]), str(interval[1]), query_id]
            filehandle.writelines('\t'.join(line) + '\n')
