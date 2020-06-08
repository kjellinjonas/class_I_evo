#!/usr/bin/env python
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.sans-serif'] = "Arial"


# identified 11nt motif and pos for all core set class I from meme output 
conserved_motif = {
    "DpuR-25": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-11": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-14": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-19": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-5": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-16": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-9": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-20": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-18": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-23": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-13": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-17": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-3": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-4": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-8": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-7": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-21": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-26": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-1": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-24": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-6": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-12": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-10": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DpuR-2": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dpu"},
    "DdR-62": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-60": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-26": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-49": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-22": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-23C": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-21": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-44": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-25": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-52": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-24B": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-46": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-23B": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-36": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-31": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-47": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-32": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-50": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-41": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-24A": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-33": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-30": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-23A": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-29": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-45": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-59": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-42": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-56": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-35": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "DdR-28": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ddi"},
    "ppa_can_1": {"pos": 6, "motif": "CCTTACAGCAA", "org": "ppa"},
    "dla_can_1": {"pos": 6, "motif": "CCTTACAGCAA", "org": "dla"},
    "DpuR-28": {"pos": 5, "motif": "TCTTACAGCAA", "org": "dpu"},
    "DdR-66": {"pos": 6, "motif": "TCTTACAGCAA", "org": "ddi"},
    "DdR-67": {"pos": 6, "motif": "TCTTACAGCAA", "org": "ddi"},
    "dla_can_199": {"pos": 6, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_3": {"pos": 6, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_4": {"pos": 6, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_9": {"pos": 6, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_11": {"pos": 6, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_190": {"pos": 5, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_27": {"pos": 5, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_26": {"pos": 5, "motif": "TCTTACAGCAA", "org": "dla"},
    "dla_can_7": {"pos": 6, "motif": "TCTTACAGCAA", "org": "dla"},
    "asu_can_3": {"pos": 6, "motif": "TCTTACAGCAA", "org": "asu"},
    "asu_can_26": {"pos": 6, "motif": "TCTTACAGCAA", "org": "asu"},
    "asu_can_2": {"pos": 6, "motif": "TCTTACAGCAA", "org": "asu"},
    "asu_can_5": {"pos": 6, "motif": "TCTTACAGCAA", "org": "asu"},
    "asu_can_36": {"pos": 6, "motif": "TCTTACAGCAA", "org": "asu"},
    "asu_can_4": {"pos": 6, "motif": "TCTTACAGCAA", "org": "asu"},
    "dla_can_87": {"pos": 6, "motif": "CCTTGCAGCAA", "org": "dla"},
    "DdR-34": {"pos": 6, "motif": "CCTCACAGCAA", "org": "ddi"},
    "ppa_can_13": {"pos": 6, "motif": "CCTTACTGCAA", "org": "ppa"},
    "ppa_can_99": {"pos": 5, "motif": "TCTTGCAGCAA", "org": "ppa"},
    "asu_can_20": {"pos": 6, "motif": "CATTACAGCAA", "org": "asu"},
    "DpuR-30": {"pos": 5, "motif": "TCTCACAGCAA", "org": "dpu"},
    "DpuR-29": {"pos": 5, "motif": "TCTCACAGCAA", "org": "dpu"},
    "dla_can_46": {"pos": 6, "motif": "CTTTACAGCAA", "org": "dla"},
    "dfa_can_113": {"pos": 5, "motif": "TCTTACTGCAA", "org": "dfa"},
    "asu_can_59": {"pos": 5, "motif": "CTTTACAGCAA", "org": "asu"},
    "ppa_can_113": {"pos": 5, "motif": "TCATACAGCAA", "org": "ppa"},
    "ppa_can_104": {"pos": 5, "motif": "TCATACAGCAA", "org": "ppa"},
    "DdR-51": {"pos": 6, "motif": "CCCTACAGCAA", "org": "ddi"},
    "DpuR-22": {"pos": 6, "motif": "CCTTACAGCAT", "org": "dpu"},
    "dla_can_6": {"pos": 6, "motif": "CCTTACAACAA", "org": "dla"},
    "dla_can_2": {"pos": 6, "motif": "CCTTACAACAA", "org": "dla"},
    "DpuR-15": {"pos": 6, "motif": "CCTTACAGCTA", "org": "dpu"},
    "ppa_can_97": {"pos": 6, "motif": "CATTGCAGCAA", "org": "ppa"},
    "ppa_can_96": {"pos": 6, "motif": "CATTGCAGCAA", "org": "ppa"},
    "ppa_can_8": {"pos": 6, "motif": "CATTGCAGCAA", "org": "ppa"},
    "ppa_can_4": {"pos": 6, "motif": "CATTGCAGCAA", "org": "ppa"},
    "ppa_can_9": {"pos": 6, "motif": "CATTGCAGCAA", "org": "ppa"},
    "ppa_can_6": {"pos": 6, "motif": "CATTGCAGCAA", "org": "ppa"},
    "ppa_can_7": {"pos": 6, "motif": "CATTGCAGCAA", "org": "ppa"},
    "ppa_can_122": {"pos": 5, "motif": "TCTCGCAGCAA", "org": "ppa"},
    "ppa_can_3": {"pos": 6, "motif": "CATCACAGCAA", "org": "ppa"},
    "ppa_can_5": {"pos": 6, "motif": "CATCACAGCAA", "org": "ppa"},
    "ppa_can_2": {"pos": 6, "motif": "CATCACAGCAA", "org": "ppa"},
    "dla_can_85": {"pos": 6, "motif": "TCTTACAGCTA", "org": "dla"},
    "dla_can_5": {"pos": 6, "motif": "CCTTACATCAA", "org": "dla"},
    "dla_can_88": {"pos": 6, "motif": "CCTTACAGCAC", "org": "dla"},
    "dfa_can_171": {"pos": 5, "motif": "TCTTACAGCTA", "org": "dfa"},
    "asu_can_18": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "asu_can_21": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "asu_can_37": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "asu_can_39": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "asu_can_19": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "asu_can_1": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "asu_can_17": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "asu_can_40": {"pos": 6, "motif": "CATCACAGCAA", "org": "asu"},
    "dfa_can_114": {"pos": 6, "motif": "CCTTCCTGCAA", "org": "dfa"},
    "dfa_can_81": {"pos": 5, "motif": "CCTTCCTGCAA", "org": "dfa"},
    "ppa_can_98": {"pos": 6, "motif": "CTTCACAGCAA", "org": "ppa"},
    "dfa_can_103": {"pos": 5, "motif": "CCTTACTACAA", "org": "dfa"},
    "dfa_can_62": {"pos": 5, "motif": "CCTTACTACAA", "org": "dfa"},
    "DpuR-27": {"pos": 6, "motif": "TATTACAGCAT", "org": "dpu"},
    "dfa_can_77": {"pos": 5, "motif": "TCGTACTGCAA", "org": "dfa"},
    "dfa_can_76": {"pos": 5, "motif": "TCGTACTGCAA", "org": "dfa"},
    "DdR-65": {"pos": 6, "motif": "CTTTACAACAA", "org": "ddi"},
    "asu_can_86": {"pos": 6, "motif": "CTTAACAGCAA", "org": "asu"},
    "ppa_can_10": {"pos": 6, "motif": "CATTATAGCAA", "org": "ppa"},
    "ppa_can_100": {"pos": 6, "motif": "CCTTCCATCAA", "org": "ppa"},
    "DdR-64": {"pos": 6, "motif": "TCATACATCAA", "org": "ddi"},
    "DdR-57": {"pos": 6, "motif": "CCTTTAAGCAA", "org": "ddi"},
    "dfa_can_102": {"pos": 5, "motif": "CCTTCCTACAA", "org": "dfa"},
    "dfa_can_78": {"pos": 6, "motif": "CCTTCCTACAA", "org": "dfa"},
    "dfa_can_79": {"pos": 6, "motif": "CCTTCCTACAA", "org": "dfa"},
    "dfa_can_57": {"pos": 5, "motif": "CCTTCCTACAA", "org": "dfa"},
    "dfa_can_115": {"pos": 6, "motif": "CCTTCCTACAA", "org": "dfa"},
    "dfa_can_88": {"pos": 6, "motif": "CCTTCCTACAA", "org": "dfa"},
    "dla_can_418": {"pos": 6, "motif": "CATTTCAGCAT", "org": "dla"},
    "dfa_can_80": {"pos": 6, "motif": "CCTTTCTACAA", "org": "dfa"},
    "dfa_can_86": {"pos": 6, "motif": "CGCTACTGCAA", "org": "dfa"},
    "asu_can_58": {"pos": 5, "motif": "CATTCCATCAA", "org": "asu"},
    "dfa_can_140": {"pos": 5, "motif": "CATTTCCGCAA", "org": "dfa"},
    "dfa_can_1": {"pos": 5, "motif": "ACTTCCAACAA", "org": "dfa"},
    "dla_can_156": {"pos": 7, "motif": "ATAAACAGCAA", "org": "dla"},
    "dfa_can_84": {"pos": 4, "motif": "CTCGTCAGCAA", "org": "dfa"},
    "dla_can_405": {"pos": 7, "motif": "TTCATCAGCAA", "org": "dla"},
}

# Most expexted accuracy predicted structure from RNAfold
mea = {"ddi": {}, "asu": {}, "dla": {},
       "dpu": {}, "dfa": {}, "ppa": {}}

with open("structure_prediction/curated_core_mea.txt") as fin:
    for line in fin:
        entry = line.strip("\n")
        if entry in conserved_motif:
            seq   = next(fin).strip("\n")
            struc = next(fin).strip("\n")
            org   = conserved_motif[entry]["org"]
            motif = conserved_motif[entry]["motif"].replace("T", "U")
            motif_seq = "-"*(conserved_motif[entry]["pos"]-1)+motif+"-"*(len(seq)-conserved_motif[entry]["pos"]-10)
            motif_struc = struc[conserved_motif[entry]["pos"]-1:conserved_motif[entry]["pos"]+10]
            mea[org][entry] = {"seq": seq, "structure": struc, "motif": motif_seq, "motif_structure": motif_struc}

labels      = ["ddi", "dpu", "dla", "ppa", "asu", "dfa"]
motif_struc = []
full_struc  = []

for org in labels:
    motif_tmp = []
    full_tmp  = []
    stem_seq_tmp = []
    stem_struc_tmp = []
    for entry in mea[org]:
        full = mea[org][entry]["structure"]
        full_tmp.append(float(full.count("(") + full.count(")")) / len(full) * 100)
        motif = mea[org][entry]["motif_structure"]
        motif_tmp.append(float(motif.count("(") + motif.count(")")) / len(motif) * 100)
    motif_struc.append(motif_tmp)
    full_struc.append(full_tmp)

# Plot structure % of full length Class I RNAs
fig, ax = plt.subplots()
plt.boxplot(full_struc, patch_artist=True, sym = "+", boxprops=dict(facecolor="darkgrey", color="black"), flierprops=dict(color="red", markeredgecolor="red"))
ax.set_ylabel("Predicted bp (%)")
ax.set_xticklabels(labels, rotation=45, ha="right")
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
plt.subplots_adjust(bottom= 0.25, top = 0.98)
plt.savefig("structure_percent.pdf")

# Plot struture % of 11 nt sequence motif
fig, ax = plt.subplots()
plt.boxplot(motif_struc, patch_artist=True, sym = "+", boxprops=dict(facecolor="darkgrey", color="black"), flierprops=dict(color="red", markeredgecolor="red"))
ax.set_ylabel("Predicted bp (%)")
ax.set_xticklabels(labels, rotation=45, ha="right")
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
plt.subplots_adjust(bottom= 0.25, top = 0.98)
plt.savefig("motif_structure_percent.pdf")

# Plot structure % per nt in 11 nt sequence motif
fig, ax = plt.subplots()
for org in labels:
    structure = [0.0]*11
    for entry in mea[org]:
        motif = mea[org][entry]["motif_structure"]
        for i in range(len(motif)):
            if motif[i] == "(" or motif[i] == ")":
                structure[i] += 1
    for i in range(len(motif)):
        structure[i] /= len(mea[org])
        structure[i] *= 100
    plt.plot(structure, label=org)

plt.xticks(np.arange(11), ["C", "C", "T", "T", "A", "C", "A", "G", "C", "A", "A"])
plt.legend()
plt.savefig("motif_structure.pdf")





