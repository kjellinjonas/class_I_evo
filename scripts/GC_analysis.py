#!/usr/bin/env python
from sys import argv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.sans-serif'] = "Arial"

def GC_genome(in_file):
    numberGC     = 0.0
    number_bases = 0.0
    with open(in_file) as fin:
        for line in fin:
            if line[0] != ">":
                line = line.strip("\n").upper()
                number_bases += len(line)
                numberGC += line.count("G") + line.count("C")

    return numberGC/number_bases * 100

def GC_class_I(in_file):
    gc = []
    with open(in_file) as fin:
        for line in fin:
            if line[0] != ">":
                line = line.strip("\n").upper()
                gc.append((line.count("G")+ line.count("C"))/float(len(line)) * 100)
    return sum(gc)/len(gc), np.std(gc)



config_file = argv[1] # tab sep file with path to class I fasta in first column and corresponding genome path in the second and label in third
plot_data = {"labels": [], "class_I_gc": [], "class_I_std": [], "genome": []}

# GC analyses of genome and class I seqs for all entries in config file
with open(config_file) as fin:
    for line in fin:
        line = line.strip("\n").split("\t")
        class_I_fasta = line[0]
        genome = line[1]
        label = line[2]
        plot_data["genome"].append(GC_genome(genome))
        class_I_gc, class_I_gc_std = GC_class_I(class_I_fasta)
        plot_data["class_I_gc"].append(class_I_gc)
        plot_data["class_I_std"].append(class_I_gc_std)
        plot_data["labels"].append(label)

# Plot results
ents = np.arange(len(plot_data["labels"]))
width = 0.35
fig, ax = plt.subplots()
class_I_serie = ax.bar(ents, plot_data["class_I_gc"], width, color="r", yerr=plot_data["class_I_std"], ecolor="black")
genome_serie = ax.bar(ents+width, plot_data["genome"], width, color="y")
ax.set_ylabel("GC %")
ax.set_xticks(ents + width / 2)
ax.set_xticklabels(plot_data["labels"], rotation=45, ha="right")
ax.legend((class_I_serie[0], genome_serie[0]), ("Class I", "Genome"), loc="upper left")
plt.subplots_adjust(bottom= 0.2, top = 0.98)
plt.savefig("class_I_gc.pdf")

# Print results as table
with open("class_I_gc.txt", "w") as fout:
    fout.write("Org\tclass_I_gc\tclass_I_std\tgenome_gc\n")
    for i in range(len(plot_data["labels"])):
        fout.write("{}\t{}\t{}\t{}\n".format(plot_data["labels"][i],
                                             plot_data["class_I_gc"][i],
                                             plot_data["class_I_std"][i],
                                             plot_data["genome"][i]))
