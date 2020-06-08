#!/usr/bin/env python
from sys import argv

def parse_published_class_I(published):
    """ Parse name and coordinates for previously
        published class I RNAs """
    pub_dict = {}
    with open(published) as fin:
        for line in fin:
            entry = line.split("\t")
            name  = entry[8].split(";")[2].strip("Name=")
            pub_dict[name] = [entry[0], entry[3], entry[4], entry[6]]
    return pub_dict

def parse_model_v1_res(infernal_out):
    """ Parse class I candidate coordinates and 
        assign an index from infernal out file
         """
    hit = 1
    res = {}
    with open(infernal_out) as fin:
        for line in fin:
            if line[0] != "#":
                line = line.split()
                chrom = line[0]
                start = line[7]
                stop  = line[8]
                strand = line[9]
                score = line[14]
                if strand == "-":
                    stop = str(int(stop)-1)
                    start, stop = stop, start
                else:
                    start = str(int(start)-1)
                res[hit] = [chrom, start, stop, strand, score]
                hit += 1
    return res

def get_pub_ids(parsed_res, parsed_published):
    """ Assign correct ID or assign candidate ID to 
        the parsed result  """
    index = 1
    updated = {}
    pub_not_found = {}
    for hit in parsed_res:
        for entry in parsed_published:
            entry_chrom = parsed_published[entry][0]
            entry_start = range(int(parsed_published[entry][1])-5,
                                int(parsed_published[entry][1])+6)
            entry_stop = range(int(parsed_published[entry][2])-5,
                               int(parsed_published[entry][2])+6)
            entry_strand = parsed_published[entry][3]
            chrom = parsed_res[hit][0] == entry_chrom
            start = int(parsed_res[hit][1]) in entry_start
            stop = int(parsed_res[hit][2]) in entry_stop
            strand = parsed_res[hit][3] == entry_strand
            if chrom and start and stop and strand:
                updated[hit] = [entry] + parsed_res[hit]
    for hit in parsed_res:
        if hit not in updated:
            new = [organism + "_can_" + str(index)] + parsed_res[hit]
            updated[hit] = new
            index += 1
    pub_ids = list(parsed_published.keys())
    for hit in updated:
        if updated[hit][0] in pub_ids:
            pub_ids.remove(updated[hit][0])
    for id in pub_ids:
        pub_not_found[id] = parsed_published[id]
    return updated, pub_not_found

def create_ref_table(path, organism, res, pub_not_found):
    header = ["ID", "chr", "start", "stop", "strand", "found with"]
    with open(path  + organism + "_class_I_candidates.txt", "w") as fout:
        fout.write("\t".join(header) + "\n")
        for entry in res:
            fout.write("\t".join(res[entry][:-1]) + "\tv1\n")
        if organism == "ddi":
            for entry in pub_not_found:
                fout.write("{}\t{}\tnf_v1\n".format(entry, "\t".join(pub_not_found[entry])))


result_v1  = argv[1]
organism   = argv[2]
genome     = argv[3]
classIpath = "analyses_output/CM_search_output/"
published  = parse_published_class_I("annotations/published_class_I.gff")
result     = parse_model_v1_res(result_v1)
parsed_res, pub_not_found = get_pub_ids(result, published)


create_ref_table(classIpath, organism, parsed_res, pub_not_found)

with open(result_v1.replace(".txt", "_all.saf"), "w") as fout:
    for entry in parsed_res:
        fout.write("\t".join(parsed_res[entry][:-1]) + "\n")

out_bed = result_v1.replace(".txt", "_25.bed")
with open(out_bed, "w") as fout:
    for entry in parsed_res:
        score = float(parsed_res[entry][-1])
        name = parsed_res[entry][0]
        if round(score) >= 25:
            if organism == "ddi":
                if "can" not in name: continue
            line = parsed_res[entry][1:4] + [parsed_res[entry][0], parsed_res[entry][-1], parsed_res[entry][4]] 
            fout.write("\t".join(line) + "\n")

with open(result_v1.replace(".txt", "_parsed.txt"), "w") as fout:
    fout.write("\t".join(["name", "chrom", "start", "stop", "strand", "score"]) + "\n")
    for entry in parsed_res:
        fout.write("\t".join(parsed_res[entry]) + "\n")
