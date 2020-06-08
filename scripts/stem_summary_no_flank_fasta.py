#!/usr/bin/env python
from sys import argv

in_fasta = argv[1] # Class I RNA fasta where each RNA starts with the perfectly conserved G and have a 2 nt 3' overhang
stem_len = int(argv[2]) # Length of stem sequence to be analyzed 
out_file = argv[3] # name of output file
stems = {}
total = 0
with open(in_fasta) as fin:
    for line in fin:
        if line[0] == ">":
            total += 1
        else:
            seq = line.strip("\n")
            five_prime_stem  = seq[:stem_len].replace("T", "U")
            three_prime_stem = seq[-2-stem_len:-2][::-1].replace("T", "U")
            comb = five_prime_stem + three_prime_stem
            if comb not in stems:
                stems[comb] = {"count": 1, "5_seq": five_prime_stem, "3_seq": three_prime_stem}
            else:
                stems[comb]["count"] += 1
bps = ["GU", "UG", "AU", "UA", "CG", "GC"]

lines = {x: [] for x in range(stem_len)}
lines["number"] = []
for comb in stems:
    lines["number"].append(str(stems[comb]["count"]))
    for i in range(stem_len-1,-1,-1):
        bp =  stems[comb]["5_seq"][i] + stems[comb]["3_seq"][i]
        if bp in bps:
             bp = stems[comb]["5_seq"][i] + "-" + stems[comb]["3_seq"][i]
        else:
            bp = stems[comb]["5_seq"][i] + " " + stems[comb]["3_seq"][i]
        lines[i].append(bp)

bp_counts = {x: {} for x in range(stem_len)}
for pos in range(stem_len):
    for i in range(len(lines[pos])):
        bp = lines[pos][i]
        count = int(lines["number"][i])
        if bp in bp_counts[pos]:
            bp_counts[pos][bp] += count
        else:
            bp_counts[pos][bp]  = count

with open(out_file, "w") as fout:
    fout.write("total: " + str(total) + "\n\n")
    fout.write("\t" + "\t".join(lines["number"]) + "\n") # print total number of Class I
    for i in range(stem_len-1, -1, -1): # print each stem together with stem count
        fout.write(str(i) + "\t" + "\t".join(lines[i]) + "\n")
    fout.write("\n\n")
    for i in range(stem_len-1, -1, -1):
        fout.write(str(i))
        for x in bp_counts[i]:
            fout.write("\t" + x + " " + str(bp_counts[i][x]) + "\t")
        fout.write("\n")


