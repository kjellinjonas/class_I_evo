#!/usr/bin/env python
from refparser import *
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatch
from matplotlib.patches import FancyBboxPatch
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.sans-serif'] = "Arial"

def plot_dist(class_I_dist, org):
    """ Create plot for genomic distribution of CLass I RNAs. Only contigs/chromosome with Class I RNAs are included """
    fig, ax = plt.subplots(figsize=[len(class_I_dist)-1, len(class_I_dist)-2])
    xmin = 0
    xmax = 1
    ymax = len(class_I_dist) * 2.5
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.set_frame_on(False)
    y = 1
    total = 0
    for contig in sorted(class_I_dist, reverse=True):

        if contig == "DDB0232429": # Add indication of chr 2 duplication in Ddi
            rep1 = [(2263132.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"],
                    (3015703.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]]
            rep2 = [(3016083.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"], 
                    (3768654.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]]
            rec1 = mpatch.Rectangle((rep1[0], y-0.1), rep1[1]-rep1[0], 0.2, fc ="lightgray", ec = "black", hatch=r"////", alpha=100.0)
            rec2 = mpatch.Rectangle((rep2[0], y-0.1), rep2[1]-rep2[0], 0.2, fc ="lightgray", ec = "black", hatch=r"\\\\", alpha=100.0)
            ax.add_patch(rec1)
            ax.add_patch(rec2)
            ax.plot([0, rep1[0]], [y, y], lw = 0.5, color="black")
            ax.plot([rep2[1], class_I_dist[contig]["norm_length"]], [y, y], lw = 0.5, color="black")
            for i in range(0, class_I_dist[contig]["length"], 500000):
                xtick = (float(i)/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]
                if xtick > rep1[0] and xtick < rep2[1]: continue
                ax.plot([xtick, xtick], [y+0.05, y-0.05], lw=0.5, color="black")
        else:
            ax.plot([0, class_I_dist[contig]["norm_length"]], [y, y], lw = 0.5, color = "black")
            for i in range(0, class_I_dist[contig]["length"], 500000):
                xtick = (float(i)/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]
                ax.plot([xtick, xtick], [y+0.05, y-0.05], lw=0.5, color="black")
        total += len(class_I_dist[contig]["class_I"])
        ax.text(-0.01, y, contig, ha="right", va="center")
        ax.text(1.1, y, str(len(class_I_dist[contig]["class_I"])), ha="left", va="center")
        ax.plot([0,0], [y+0.1, y-0.1], lw=0.5, color="black")
        ax.plot([class_I_dist[contig]["norm_length"], class_I_dist[contig]["norm_length"]], [y+0.1, y-0.1], lw=0.5, color="black")
        ax.text(0, y+0.12, "0", ha="right", va="bottom", size=7)
        ax.text(class_I_dist[contig]["norm_length"], y+0.12, str(class_I_dist[contig]["length"]), ha="left", va="bottom", size=7)
        if contig == "DDB0232429":
            rep1 = [(2263132.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"],
                    (3015703.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]]
            rep2 = [(3016083.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"], 
                    (3768654.0/class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]]
            rec1 = mpatch.Rectangle((rep1[0], y-0.1), rep1[1]-rep1[0], 0.2, fc ="lightgray", ec = "black", hatch=r"////", alpha=100.0)
            rec2 = mpatch.Rectangle((rep2[0], y-0.1), rep2[1]-rep2[0], 0.2, fc ="lightgray", ec = "black", hatch=r"\\\\", alpha=100.0)
            ax.add_patch(rec1)
            ax.add_patch(rec2)
        sorted_starts = {class_I_dist[contig]["start_stop"][i][0]: {"stop": class_I_dist[contig]["start_stop"][i][1], 
                        "strand": class_I_dist[contig]["strand"][i], "name": class_I_dist[contig]["class_I"][i]} for i in range(len(class_I_dist[contig]["start_stop"]))}
        latest_start = 0.0
        arrow_y = y
        duplicated = {"dpu_can_597": "blue", "dpu_can_598": "blue", "dpu_can_240": "green", "dpu_can_241": "green", "dpu_can_305": "red", "dpu_can_304": "red", 
                      "pvi_can_12": "blue", "pvi_can_13": "blue", "pvi_can_29": "blue", "ppa_can_6": "blue", "ppa_can_96": "blue", "r24A": "blue", 
                      "r24B": "blue", "r23A": "green", "r23B": "green", "r23C": "green", "dpc_can_18": "blue", "dpc_can_16": "blue",
                      "dpc_can_17": "blue", "dpc_can_8": "green", "dpc_can_7": "green"}
        for i in sorted(sorted_starts.keys()):
            start = i
            stop  = sorted_starts[i]["stop"]
            
            if latest_start == 0 and latest_start + 0.01 > start:
                arrow_y = y - 0.25
                latest_start = start
            elif latest_start != 0 and latest_start + 0.01 > start:
                arrow_y -= 0.2
                latest_start = start
            else:
                arrow_y = y - 0.25
                latest_start = start
            if sorted_starts[i]["strand"] == "-":
                a_start = (stop, arrow_y)
                a_stop = (start, arrow_y)
            else:
                a_start = (start, arrow_y)
                a_stop = (stop, arrow_y)
            facecol = "black"
            if sorted_starts[i]["name"] in duplicated:
                facecol = duplicated[sorted_starts[i]["name"]]
            arrow = mpatch.FancyArrowPatch(a_start, a_stop, fc=facecol, ec=facecol, mutation_scale=10)
            ax.add_patch(arrow)
        y += 1.7
    plt.title("{} Class I RNAs (n={})".format(org, total))
    plt.tight_layout()
    plt.savefig(org + "_genomic_dist.pdf", dpi=300, format="pdf")


def collect_class_I(bed, fai):
    """ Parse Class I genomic information """
    class_I_dist = {}
    for class_I in bed.names:
        contig = bed.get_chr(class_I)
        if contig not in class_I_dist:
            c_length = fai.get_length(contig)
            class_I_dist[contig] = {"norm_length": 0, "length": c_length, "class_I": [], "start_stop": [], "strand": []}
        start  = float(bed.get_start(class_I))
        stop   = float(bed.get_stop(class_I))
        strand = bed.get_strand(class_I)
        class_I_dist[contig]["class_I"].append(class_I)
        class_I_dist[contig]["start_stop"].append([start, stop])
        class_I_dist[contig]["strand"].append(strand)
    return class_I_dist

def normalize_contig_lengths(class_I_dist):
    """ Normalize contig/chromosome lengths relative to the longes one """
    class_I_dist = class_I_dist
    lengths = {class_I_dist[contig]["length"]: contig for contig in class_I_dist}
    if len(lengths) != len(class_I_dist): print("something is wrong.")
    by_size = sorted(lengths.keys(), reverse=True)
    for size in by_size:
        contig = lengths[size]
        class_I_dist[contig]["norm_length"] = float(class_I_dist[contig]["length"])/max(by_size)
    return class_I_dist

def normalize_coordinate(class_I_dist):
    """ Adjust Class I RNA genomic coordinates according to the normalized contig lengths """
    class_I_dist = class_I_dist
    for contig in class_I_dist:
        for i in range(len(class_I_dist[contig]["start_stop"])):
            start = class_I_dist[contig]["start_stop"][i][0]
            stop = class_I_dist[contig]["start_stop"][i][1]
            start = (start / class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]
            stop = (stop / class_I_dist[contig]["length"]) * class_I_dist[contig]["norm_length"]
            class_I_dist[contig]["start_stop"][i] = [start, stop]
    return class_I_dist

orgs = {"ddi": {"gff": ddi_gff, "bed": ddi_bed, "id": ddi_id, "fai": ddi_fai},  # more organisms can be added, but then they need to be added also in refparser.py
        "dfi": {"gff": dfi_gff, "bed": dfi_bed, "id": dfi_id, "fai": dfi_fai}}


for org in orgs:
    bed = orgs[org]["bed"]
    fai = orgs[org]["fai"]
    class_I_dist = collect_class_I(bed, fai)
    class_I_dist = normalize_contig_lengths(class_I_dist)
    class_I_dist = normalize_coordinate(class_I_dist)
    plot_dist(class_I_dist, org)




