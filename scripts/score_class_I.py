#!/usr/bin/env python
from sys import argv

def parse_fimo(motif_table, infernal):
    """ Parse the result from the motif identification upstream of candidates """
    motif = {x: {"DUSE": {"score": -10,
                         "seq": "n/a",
                         "dist": 0,
                         "pval": 1},
                "TGTG": {"score": 0,
                         "seq": "n/a",
                         "dist": 0,
                         "pval": 1},
                "dist": 0} for x in infernal}
    with open(motif_table) as fin:
        next(fin)
        for line in fin:
            motif_id = "DUSE"
            line = line.rstrip("\n").split("\t")
            if len(line) < 2: continue
            if line[5] == "-": continue
            if line[1] == "MEME-2": continue # This motif was not used
            if line[1] == "MEME-3":
                motif_id = "TGTG"
            classI_id = line[2].split("|")[0]
            if float(line[7]) < motif[classI_id][motif_id]["pval"]:
                motif[classI_id][motif_id]["pval"] = float(line[7])
                motif[classI_id][motif_id]["dist"] = 150 - int(line[4])
                motif[classI_id][motif_id]["score"] = float(line[6])
                if motif_id == "TGTG" and motif[classI_id][motif_id]["score"] < 0:
                    motif[classI_id][motif_id]["score"] = 0
                motif[classI_id][motif_id]["seq"] = line[9]
    return motif

def parse_infernal(infernal_out):
    """ Parse the CM search output """
    hit = 1
    res = {}
    with open(infernal_out) as fin:
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
    return res

def parse_candidates(classIref):
    """ Parse previously identified candidates """
    candidates = {}
    with open(classIref) as fin:
        header = next(fin)
        for line in fin:
            entry = line.rstrip("\n").split("\t")
            name = entry[0]
            chrom = entry[1]
            start = entry[2]
            stop = entry[3]
            if int(start) > int(stop):
                start, stop = stop, start
            candidates[name] = [chrom, start, stop] + entry[4:]
    return candidates

def total_score(fimo, infernal):
    """ Calculate the classifier score """
    scores = {x: 0 for x in fimo}
    for can in fimo:
        scores[can] += float(infernal[can][-1]) # Infernal score
        motif_dist = fimo[can]["DUSE"]["dist"] - fimo[can]["TGTG"]["dist"]
        scores[can] += fimo[can]["DUSE"]["score"] # FIMO score if found else -10
        
        # Distance score. +5 if motif dist or DUSE dist in range else -5
        if fimo[can]["DUSE"]["dist"] in range(56, 65):
            scores[can] += 5
        elif motif_dist in range(56, 65):
            scores[can] += 5
        else:
            scores[can] -= 5
    return scores

def parse_names(infernal, classIref):
    name = {}
    for hit in infernal:
        for entry in classIref:
            entry_chrom = classIref[entry][0]
            entry_start = range(int(classIref[entry][1])-5,
                                int(classIref[entry][1])+6)
            entry_stop = range(int(classIref[entry][2])-5,
                               int(classIref[entry][2])+6)
            entry_strand = classIref[entry][3]
            chrom = infernal[hit][0] == entry_chrom
            start = int(infernal[hit][1]) in entry_start
            stop = int(infernal[hit][2]) in entry_stop
            strand = infernal[hit][3] == entry_strand
            if chrom and start and stop and strand:
                name[entry] = infernal[hit]
    return name

infernal_out = argv[1] # out file from CM search
fimo_out = argv[2] # the tab seperated output from FIMO motif identification
org = argv[3] # organism identifier, e.g. ddi for D. discoideum
infernal = parse_infernal(infernal_out)
classIref = parse_candidates("analyses_output/CM_search_output/" + org + "_class_I_candidates.txt")
infernal = parse_names(infernal, classIref)
fimo = parse_fimo(fimo_out, infernal)
scores = total_score(fimo, infernal)
header = ["entry", "total_score", "infernal_score", "DUSE", 
          "duse_score", "duse_dist", "TGTG", "TGTG_score", 
          "TGTG_dist", "inter_dist"]

# Print classifier summary table
with open(infernal_out[infernal_out.rfind("/")+1:].replace(".txt", "_score.txt"), "w") as fout:
    fout.write("\t".join(header) + "\n")
    for entry in infernal:
        line = "\t".join(["{}"] * 10) + "\n"
        fout.write(line.format(entry, scores[entry], infernal[entry][-1],
                               fimo[entry]["DUSE"]["seq"], fimo[entry]["DUSE"]["score"], fimo[entry]["DUSE"]["dist"],
                               fimo[entry]["TGTG"]["seq"], fimo[entry]["TGTG"]["score"], fimo[entry]["TGTG"]["dist"],
                               fimo[entry]["DUSE"]["dist"] - fimo[entry]["TGTG"]["dist"]))