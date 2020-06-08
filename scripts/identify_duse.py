#!/usr/bin/env python

from sys import argv
from subprocess import call
from Bio import SeqIO, SeqRecord, Seq

def parse_bed(in_file):
    """ Parse genomic information from bed file and adjust to 150 nt upstream of each entry """
    upstream = {}
    downstream = {}
    with open(in_file) as fin:
        for line in fin:
            line = line.split("\t")
            name = line[3]
            chrm = line[0]
            low  = int(line[1])
            high = int(line[2])
            strd = line[5].strip("\n")
            if strd == "+":
                upstream[name] = {"chrm": chrm, "low": low-150, 
                                    "high": low, "strand": strd}
                if upstream[name]["low"] < 0: 
                    print("warning: upstream region for {} extended outside outside of".format(
                          name) +
                          " available genome sequence ({}). adjusted to 0".format(
                            upstream[name]["low"]))
                    upstream[name]["low"] = 0
            else:
                upstream[name] = {"chrm": chrm, "low": high-1, 
                                    "high": high+149, "strand": strd}
    return upstream


def get_sequences(regions, indexed_genome):
    """ Get nucleotide sequences from parsed bed file """
    sequences = {}
    for entry in regions:
        chrm = regions[entry]["chrm"]
        strd = regions[entry]["strand"]
        low  = regions[entry]["low"]
        high = regions[entry]["high"]
        header = "{}|{}:{}-{}|{}".format(entry, chrm, low, high, strd)
        
        seq = indexed_genome[chrm].seq[low:high]
        if strd == "-":
            sequences[header] = str(seq.reverse_complement())
        else:
            sequences[header] = str(seq)
    return sequences


def write_fasta(seq_dict, outfile):
    """ write fasta file with 150 nt upstream sequence of each entry in bed """
    with open(outfile, "w") as fout:
        for entry in seq_dict:
            fout.write(">{}\n{}\n".format(entry, seq_dict[entry]))


def run_fimo(fasta_file, motif):
    """ identify motifs in the upstream sequences """
    outdir = fasta_file.replace(".fa", "_fimo/")
    call("fimo --thresh 1e-2 {} {}".format(
         motif, fasta_file), shell=True)
    return outdir


				



motif = "meme_out/meme.txt"
bed_file = argv[1]
genome   = SeqIO.index(argv[2], "fasta")
fasta_out = bed_file.replace(".bed", "_150bp_upstream.fa")

upstream_regions = parse_bed(bed_file)
upstream_seqs = get_sequences(upstream_regions, genome)
fasta_out = bed_file.replace(".bed", "_150bp_upstream.fa")
write_fasta(upstream_seqs, fasta_out)
motif_folder = run_fimo(fasta_out, motif)

