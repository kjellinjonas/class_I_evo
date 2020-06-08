#!/usr/bin/env python
from sys import argv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.sans-serif'] = "Arial"

def len_class_I(in_file):
    lens = []
    with open(in_file) as fin:
        for line in fin:
            if line[0] != ">":
                line = line.strip("\n")
                lens.append(len(line))
    return lens



config_file = argv[1] # tab sep file with path to class I fasta in first column and label in third
plot_data = {"labels": [], "class_I_len": [], "class_I_std": []}

# Length analysis of Class I RNAs for all organisms included in config file
with open(config_file) as fin:
    for line in fin:
        line = line.strip("\n").split("\t")
        class_I_fasta = line[0]
        label = line[2]
        class_I_lens = len_class_I(class_I_fasta)
        plot_data["class_I_len"].append(class_I_lens)
        plot_data["labels"].append(label)

# Plot results
ents = np.arange(len(plot_data["labels"]))
fig, ax = plt.subplots()
plt.boxplot(plot_data["class_I_len"], patch_artist=True, sym = "+", boxprops=dict(facecolor="darkgrey", color="black"), flierprops=dict(color="red", markeredgecolor="red"))
ax.set_ylabel("length (nt)")
ax.set_xticklabels(plot_data["labels"], rotation=45, ha="right")
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
plt.subplots_adjust(bottom= 0.25, top = 0.98)
plt.savefig("class_I_len_boxplot.pdf")
