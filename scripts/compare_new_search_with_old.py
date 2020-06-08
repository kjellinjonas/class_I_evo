#!/usr/bin/env python
from sys import argv
from subprocess import call

def parse_prev_cand(cands):
        """ Parse text document with candidates identified in previous searches """
        can_index = []
        curated_table = {}
        with open(cands) as fin:
            header = next(fin)
            for line in fin:
                entry = line.rstrip("\n").split("\t")
                name  = entry[0]
                chrom = entry[1]
                start = entry[2]
                stop  = entry[3]
                if int(start) > int(stop):
                    start, stop = stop, start
                if "can" in name:
                    can_index.append(int(name.strip(org + "_can_")))
                curated_table[name] = [chrom, start, stop] + entry[4:]
        return curated_table, can_index


def parse_new_res(new_res, prev_can_index):
    """ Parse genomic coordinate information for candidates identified in
    the search and each is given an index number  """
    can = max(prev_can_index) + 1 
    res = {}
    hit = can
    with open(new_res) as fin:
        for line in fin:
            if line[0] != "#":
                line  = line.split()
                chrom  = line[0]
                start  = line[7]
                stop   = line[8]
                strand = line[9]
                score  = line[14]
                if strand == "-":
                    stop = str(int(stop)-1)
                    start, stop = stop, start
                else:
                    start = str(int(start)-1)
                res[hit] = [chrom, start, stop, strand, score]
                hit += 1
    return res, can


def compare_results(new_res, prev_res, can):
    """ Compare the candidates identified in the latest search with previously identified candidates.
    """
    found = []  # Overlapping candidates will be stored here
    new_can = [] # New candidates will be stored here 
    parsed_res = {} # genomic information for all cands here
    saf = {} # information used to crease saf annotation file here
    # Identify overllaping candidates
    for hit in new_res: 
        for entry in prev_res:
            entry_chrom = prev_res[entry][0]
            entry_start = range(int(prev_res[entry][1])-5,
                                int(prev_res[entry][1])+6)
            entry_stop = range(int(prev_res[entry][2])-5,
                               int(prev_res[entry][2])+6)
            entry_strand = prev_res[entry][3]
            chrom = new_res[hit][0] == entry_chrom
            start = int(new_res[hit][1]) in entry_start
            stop = int(new_res[hit][2]) in entry_stop
            strand = new_res[hit][3] == entry_strand
            if chrom and start and stop and strand:
                parsed_res[hit] = [entry] + new_res[hit][:5] + prev_res[entry][5:]
                saf[hit] = [entry] + res[hit]
                found.append(entry)
    # Identidy new candidates
    for hit in res:
        if hit not in parsed_res:
            entry = org + "_can_" + str(can) 
            new_can.append(hit)
            new = [entry] + new_res[hit]
            parsed_res[hit] = new
            saf[hit] = new
            can += 1
    return new_can, parsed_res, saf, found



# Input files
result_new = argv[1] # CM output file
org = argv[2] # organism identifier, e.g. ddi for D. discoideum
model = argv[3] # CM model version
classIref = "analyses_output/CM_search_output/" + org + "_class_I_candidates.txt" # file with candidates found in previous searches
# Output files
out_saf = result_new.replace(".txt", "_all.saf") # Saf annotation file to be used with featureCounts
out_bed = result_new.replace(".txt", "_all.bed") # bed annotation file with all candidates
out_all_bed = result_new.replace(".txt", "_all25.bed") # bed annotation file with all candidates ≥ 25 in CM
out_new_bed = result_new.replace(".txt", "_new25.bed") # bed annotation file with all new candidates ≥ 25 in CM

# Parse hits from previous searches
curated_table, can_index = parse_prev_cand(classIref)
# Parse hits from the current search
res, can = parse_new_res(result_new, can_index)

new_can, parsed_res, saf, found = compare_results(res, curated_table, can)
    
# Create bed file with all new hits ≥ 25 in CM
with open(out_new_bed, "w") as fout:
    for hit in new_can:
        line = parsed_res[hit][1:4] + [parsed_res[hit][0], parsed_res[hit][-1], parsed_res[hit][-2]]
        if round(float(line[-2])) >= 25:
            fout.write("\t".join(line) + "\n")
    
# Create bed file with all candidates ≥ 25 in CM
with open(out_all_bed, "w") as fout:
    for hit in parsed_res:
        line = parsed_res[hit][1:4] + [parsed_res[hit][0], parsed_res[hit][-1], parsed_res[hit][-2]]
        if round(float(line[-2])) >= 25:
            fout.write("\t".join(line) + "\n")

# Create bed file with all candidates
with open(out_bed, "w") as fout:
    for hit in parsed_res:
        line = parsed_res[hit][1:4] + [parsed_res[hit][0], parsed_res[hit][-1], parsed_res[hit][-2]]
        fout.write("\t".join(line) + "\n")
    
# Create saf annotation file with all candidates for featureCounts
with open(out_saf, "w") as fout:
    for entry in saf:
        fout.write("\t".join(saf[entry][:5]) + "\n")
    
# Update the classIref file with the result from latest search
with open(classIref, "w") as fout:
    fout.write("ID\tchr\tstart\tstop\tstrand\tfound with\n")
    for entry in sorted(curated_table):
        if entry not in found:
            fout.write("\t".join([entry] + curated_table[entry]) + "\n") # Print entry as before
        else:
            fout.write("\t".join([entry] + curated_table[entry]) + ", {}\n".format(model)) # Add that this candidate also was identified in latest search
    for entry in new_can:
        fout.write("\t".join(parsed_res[entry][:-1]) + "\t{}\n".format(model)) # Add newly identified candidates



