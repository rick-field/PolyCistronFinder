import sys
from random import sample




def ppl_bed_dict_maker(ppl_file):
    ppl_bed_dict = {}
    with open(ppl_file, "r") as ppl_bed:
        for row in ppl_bed:
            row = row.strip().split('\t')
            ppl = row[-1]
            try:
                ppl_bed_dict[ppl].append(row)
            except KeyError:
                ppl_bed_dict[ppl] = [row]

    return ppl_bed_dict




def genome_bed_list_maker(genome_file):
    genome_bed_list = []
    with open(genome_file, "r") as ppl_bed:
        for row in ppl_bed:
            row = row.strip().split('\t')
            genome_bed_list.append(row)

    return genome_bed_list




def interval_calculator(ppl_bed_dict, genome_bed_list):
    ppl_distances_list = []
    ppl_distances = {}
    for ppl in ppl_bed_dict:
        count = 0
        while count < len(ppl_bed_dict[ppl]) - 1:
            intergenic_distance = int(ppl_bed_dict[ppl][count + 1][1]) - int(ppl_bed_dict[ppl][count][2])
            ppl_distances_list.append(intergenic_distance)
            try:
                ppl_distances[ppl].append(intergenic_distance)
                count += 1
            except KeyError:
                ppl_distances[ppl] = [intergenic_distance]
                count += 1

    total_distance_bp = 0
    total_distance_locus_count = len(ppl_distances_list)
    for i in ppl_distances_list:
        total_distance_bp += i
    avg_ppl_intergenic_distance = total_distance_bp / total_distance_locus_count
    print("Average PPL intergenic distance: " + str(avg_ppl_intergenic_distance))



    genome_distances_list = []
    genome_distances = {}
    rand_samp = sample(list(range(len(genome_bed_list))), len(genome_bed_list))
    i = 0
    count = 0
    while i < len(ppl_distances_list) - 1:
        rand_gene_index = int(rand_samp[i])
        next_gene_index = int(rand_samp[i]) + 1
        # Calculate distance if genes are on the same chromosome and on same strand
        if (genome_bed_list[next_gene_index][5] == genome_bed_list[rand_gene_index][5]) and (genome_bed_list[next_gene_index][0] == genome_bed_list[rand_gene_index][0]):
            intergenic_distance = int(genome_bed_list[next_gene_index][1]) - int(genome_bed_list[rand_gene_index][2])
            genome_distances_list.append(intergenic_distance)
            i += 1
            count += 1

        else:
            i += 1

    total_distance_bp = 0
    total_distance_locus_count = len(genome_distances_list)
    for i in genome_distances_list:
        total_distance_bp += i
    avg_genome_intergenic_distance = total_distance_bp / total_distance_locus_count
    print("Average genome wide intergenic distance: " + str(avg_genome_intergenic_distance))

    return ppl_distances, genome_distances_list


ppl_file = str(sys.argv[1])
genome_file = str(sys.argv[2])
prefix = str(sys.argv[3])

ppl_bed_dict = ppl_bed_dict_maker(ppl_file)
genome_bed_list = genome_bed_list_maker(genome_file)
ppl_distances, genome_distances_list = interval_calculator(ppl_bed_dict, genome_bed_list)

with open(prefix + "_PPL_intergenic_distances.txt", "w") as filehandle:
    for ppl in ppl_distances:
        filehandle.writelines(ppl + '\t' + str(ppl_distances[ppl]))

with open(prefix + "_genome_wide_intergenic_distances.txt", "w") as filehandle:
    for distance in genome_distances_list:
        filehandle.writelines(str(distance) + '\n')
