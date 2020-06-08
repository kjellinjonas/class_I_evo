#!/usr/bin/env python
from sys import argv
from Bio import SeqIO, SeqRecord, Seq

def get_coords(in_bed):
    regions = {}
    with open(in_bed) as fin:
        for line in fin:
            details = line.strip("\n").split("\t")
            entry = details[3]
            strand = details[5]
            chrm = details[0]
            if strand == "+":
                low = int(details[1]) - 1 
                high = int(details[2]) + 20
            else:
                low = int(details[1]) - 19
                high = int(details[2]) + 2
            regions[entry] = [chrm, low, high, strand]
    return regions

def get_sequences(regions, indexed_genome):
    sequences = {}
    for entry in regions:
        chrm = regions[entry][0]
        strd = regions[entry][3]
        low  = regions[entry][1]
        high = regions[entry][2]
        header = "{}|{}:{}-{}|{}".format(entry, chrm, low, high, strd)
        seq = indexed_genome[chrm].seq[low:high]
        if strd == "-":
            sequences[header] = str(seq.reverse_complement()).upper()
        else:
            sequences[header] = str(seq).upper()
    return sequences


def write_fasta(seq_dict, outfile):
    with open(outfile, "w") as fout:
        for entry in seq_dict:
            fout.write(">{}\n{}\n".format(entry, seq_dict[entry]))


config_file = argv[1]
with open(config_file) as fin:
    for line in fin:
        line = line.strip("\n").split("\t")
        class_I_bed = line[3]
        genome = SeqIO.index(line[1], "fasta")
        label = line[2]
        regions = get_coords(class_I_bed)
        fasta_out = class_I_bed[class_I_bed.rfind("/")+1:].replace(".bed", "_flank.fa")
        class_I_seqs = get_sequences(regions, genome)
        write_fasta(class_I_seqs, fasta_out)