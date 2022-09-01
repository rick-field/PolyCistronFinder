import sys

def blacklister(samfile):
    blacklist_dict = {}
    with open(samfile, "r") as sam:
        for row in sam:
            if row[0] == "@":
                pass
            else:
                row = row.split('\t')
                cp_name = row[0]
                flag = row[1]
                chrom = row[2]
                start = row[3]
                mq = row[4]
                cigar = row[5]
                length = len(row[9])
                stop = str(int(start) + length)
                if flag == "2048" or flag == "2064":
                    pass
                elif int(mq) < 20:
                    pass
                else:
                    line = [chrom, start, stop, cp_name, mq]
                    try:
                        blacklist_dict[row[0]].append([start, stop])
                    except KeyError:
                        blacklist_dict[row[0]] = [[start, stop]]
                    with open(prefix + "_PCF_blacklisted_regions.bed", "a") as out_blacklist:
                        out_blacklist.writelines('\t'.join(line))

    return blacklist_dict




def blacklist_filter(ppl_bed, blacklist_dict):
    ppl_bed_dict = {}
    with open(ppl_bed, "r") as bed:
        for ppl_row in bed:
            row = ppl_row.strip().split()
            chrom = row[0]
            start = row[1]
            stop = row[2]
            try:
                for interval in blacklist_dict[chrom]:
                    black_range = list(range(int(interval[0]),int(interval[1])))
                    if start in black_range or stop in black_range:
                        with open(prefix + "_PPL_blacklisted.bed", "a") as out_blacklist:
                            out_blacklist.writelines('\t'.join(row))
                    else:
                        print('\t'.join(row))
            except KeyError:
                print('\t'.join(row))




samfile = str(sys.argv[1])
ppl_bed = str(sys.argv[2])
prefix = str(sys.argv[3])

blacklist_dict = blacklister(samfile)
blacklist_filter(ppl_bed, blacklist_dict)
