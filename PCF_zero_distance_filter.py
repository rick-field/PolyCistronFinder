import sys
import subprocess

ppl_bed = str(sys.argv[1])
ppl_bed_basename = ppl_bed.strip(".bed")
prefix = str(sys.argv[2])



ppl_bed_dict = {}
with open(ppl_bed, "r") as ppl_bed:
    for row in ppl_bed.readlines():
        row = row.strip().split('\t')
        ppl = row[-1]
        try:
            ppl_bed_dict[ppl].append(row)
        except KeyError:
            ppl_bed_dict[ppl] = [row]

final_ppl_bed_dict = {}
bad_ppl_dict = {}
ppl_ids = [i for i in ppl_bed_dict]
for ppl in ppl_ids:
    i = 0
    try:
        while i < len(ppl_bed_dict[ppl]) - 1:
            first_gene = ppl_bed_dict[ppl][i]
            next_gene = ppl_bed_dict[ppl][i+1]
            if int(next_gene[1]) - int(first_gene[2]) < 1:
                bad_ppl_dict[ppl] = ppl_bed_dict[ppl]
                del ppl_bed_dict[ppl]
                i += 1
            else:
                i += 1
    except KeyError:
        pass

with open(ppl_bed_basename + "_final.bed", "w") as filehandle:
     for ppl in ppl_bed_dict:
         for row in ppl_bed_dict[ppl]:
             filehandle.writelines('\t'.join(row) + '\n')

with open(ppl_bed_basename + "_final_gene_IDs.txt", "w") as filehandle:
     for ppl in ppl_bed_dict:
         for row in ppl_bed_dict[ppl]:
             filehandle.writelines(row[3] + '\n')

with open(ppl_bed_basename + "_final_ppl_IDs.txt", "w") as filehandle:
     for ppl in ppl_bed_dict:
         filehandle.writelines(ppl + '\n')

with open(ppl_bed_basename + "_zero_intergenic_distance.bed", "w") as filehandle:
     for ppl in	bad_ppl_dict:
       	 for row in bad_ppl_dict[ppl]:
       	     filehandle.writelines('\t'.join(row) + '\n')

with open(ppl_bed_basename + "_zero_intergenic_distance_ppl_IDs.txt", "w") as filehandle:
     for ppl in bad_ppl_dict:
         filehandle.writelines(ppl + '\n')

cmd = "grep -f " + prefix + "_PPL_final_gene_IDs.txt " + prefix + "_COMPLETE_blast_NO_mutual_hsps.txt > " + prefix + "_COMPLETE_blast_NO_mutual_hsps_filtered.txt"
subprocess.run(cmd, shell = True)

cmd = "grep -f " + prefix + "_PPL_final_gene_IDs.txt " + prefix + "_PARTIAL_blast_NO_mutual_hsps.txt > " + prefix + "_PARTIAL_blast_NO_mutual_hsps_filtered.txt"
subprocess.run(cmd, shell = True)

cmd = "grep -f " + prefix + "_PPL_final_gene_IDs.txt " + prefix + "_COMPLETE_blast_mutual_hsps.txt > " + prefix + "_COMPLETE_blast_mutual_hsps_filtered.txt"
subprocess.run(cmd, shell = True)

cmd = "grep -f " + prefix + "_PPL_final_gene_IDs.txt " + prefix + "_PARTIAL_blast_mutual_hsps.txt > " + prefix + "_PARTIAL_blast_mutual_hsps_filtered.txt"
subprocess.run(cmd, shell = True)


cmd = "python PCF_blast_result_split.py " + prefix + "_COMPLETE_blast_mutual_hsps_filtered.txt > " + prefix + "_COMPLETE_blast_mutual_hsps_filtered_gene_IDs.txt"
subprocess.run(cmd, shell = True)

cmd = "python PCF_blast_result_split.py " + prefix + "_COMPLETE_blast_NO_mutual_hsps_filtered.txt > " + prefix + "_COMPLETE_blast_NO_mutual_hsps_filtered_gene_IDs.txt"
subprocess.run(cmd, shell = True)

cmd = "python PCF_blast_result_split.py " + prefix + "_PARTIAL_blast_mutual_hsps_filtered.txt > " + prefix + "_PARTIAL_blast_mutual_hsps_filtered_gene_IDs.txt" 
subprocess.run(cmd, shell = True)

cmd = "python PCF_blast_result_split.py " + prefix + "_PARTIAL_blast_NO_mutual_hsps_filtered.txt > " + prefix + "_PARTIAL_blast_NO_mutual_hsps_filtered_gene_IDs.txt"
subprocess.run(cmd, shell = True)
