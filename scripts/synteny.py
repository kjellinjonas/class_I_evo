#!/usr/bin/env python
from sys import argv
from refparser import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.patches as mpatch
from matplotlib.patches import FancyBboxPatch
from matplotlib.pyplot import cm
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.sans-serif'] = "Arial"

def get_xlims(xmin, xmax, class_I, genes):
    xmin = xmin
    xmax = xmax
    for f in flanking_class_I:
        start = bed.get_start(f)
        stop = bed.get_stop(f)
        if start < xmin:
            xmin = start - 20
        if stop > xmax:
            xmax = stop + 20
    for f in flanking_genes:
        if gff.get_name(f) in flanking_class_I: continue
        start = gff.get_start(f)
        stop = gff.get_stop(f)
        if start < xmin:
            xmin = start - 20
        if stop > xmax:
            xmax = stop + 20
    return xmin, xmax


def plot_context(entry, contig, start, stop, flanking_class_I, flanking_genes, ortho_genes, ortho_class_I, org, orgs, flank):
    y = {"ddi": 0, "dfi":0, "dpu":0, "dla":0, "ppa": 0, "asu": 0, "dfa": 0}
    ortho_order = [x for x in orgs]
    fig, ax = plt.subplots(figsize = (8,6))
    xmin = start-1000
    xmax = start+1000
    y_org = y[org]
    contig_y= {x: {} for x in y}
    xmin, xmax = get_xlims(xmin, xmax, flanking_class_I, flanking_genes)
    ax.set_xlim(0,1)
    ax.text(-0, 50, "{}: {}".format(org, entry), size=10)
    number = 0
    tmp_orthos = {}
    #print("parsing and plotting class I ({}) information". format(entry))
    # add Class I to plot
    for f in flanking_class_I:
        contig = bed.get_chr(f)
        start = bed.get_start(f)
        stop = bed.get_stop(f)
        n_start, n_stop, n_max = float(start-xmin), float(stop-xmin), float(xmax-xmin)
        s_start, s_stop, s_max = n_start/n_max, n_stop/n_max, 1
        if bed.get_strand(f) == "-":
            y_tmp = y_org - 3
            al = "top"
            rot = -35
        else:
            y_tmp= y_org  + 2
            al = "bottom"
            rot = 35
        
        rectangle = mpatch.Rectangle((s_start, y_tmp), s_stop-s_start, 1, fc = "black", ec = None)
        ax.add_patch(rectangle)
        ax.text(s_start, y_tmp, f, size = 5, rotation=rot, ha="left", va=al, fontweight='bold')


    #print("Starting to create plot for {}".format(entry)) 
    # Add flanking genes to plot   
    for f in flanking_genes:
        if gff.get_name(f) in flanking_class_I: continue
        start = gff.get_start(f)
        stop = gff.get_stop(f)
        n_start, n_stop, n_max = float(start-xmin), float(stop-xmin), float(xmax-xmin)
        s_start, s_stop, s_max = n_start/n_max, n_stop/n_max, 1
        if gff.get_strand(f) == "-":
            y_tmp = y_org -3
            al = "top"
            rot = -35
        else:
            rot = 35
            y_tmp = y_org +2
            al = "bottom"
        tmp_orthos[f] = s_start
        rectangle = mpatch.Rectangle((s_start, y_tmp), s_stop-s_start, 1, fc = "black", ec = None)
        ax.add_patch(rectangle)
        col = "black"
        ax.text(s_start, y_tmp, gff.get_name(f).strip("\n"), size = 5, rotation=rot, color =col, ha="left", va=al)
    x_ticks = np.arange(0, 1, 1000.0/n_max)
    for coord in x_ticks:
        ax.plot([coord, coord], [y[org]+1,y[org]-1], lw = 0.5, color ="darkgray")

    entries = 50 # Used to adjust y coordinates in plot
    written = [] # Information about plotted chromsomes/contigs stored here
    for o in ortho_order:
        ortho_context = {}
        tmp_gff = orgs[o]["gff"]
        tmp_bed = orgs[o]["bed"]
        start_pos={}
        class_I_starts = {}
        for ngene in ortho_class_I[o]:
            if len(ortho_class_I[o][ngene]) > 0:
                number += 1
                if o not in ortho_context: ortho_context[o] = {}
                for cI in ortho_class_I[o][ngene]:
                    contig = tmp_bed.get_chr(cI)
                    if contig not in class_I_starts: class_I_starts[contig] = []
                    start = tmp_bed.get_start(cI)
                    if contig not in ortho_context[o]: ortho_context[o][contig]= {"gene": [], "class_I": []}
                    if cI not in ortho_context[o][contig]["class_I"]: ortho_context[o][contig]["class_I"].append(cI)
                    if contig not in start_pos: start_pos[contig] = {}
                    start_pos[contig][start] = cI
                    class_I_starts[contig].append(start)
        for fgene in ortho_genes[o]:
            if o not in ortho_context: continue # only care about orthos if there is a potentially conserved class I
            for ngene in ortho_genes[o][fgene]:
                contig = tmp_gff.get_chr(ngene)
                start  = tmp_gff.get_start(ngene)
                if contig in ortho_context[o]:
                    if ngene not in ortho_context[o][contig]["gene"]:
                        ortho_context[o][contig]["gene"].append(ngene)
                        start_pos[contig][start] = ngene
        
        for contig in start_pos:
            tmp_y = y[org] - entries
            entries += 50
            start = min(start_pos[contig].keys())
            stop = tmp_gff.get_stop(start_pos[contig][max(start_pos[contig].keys())])
            if "not found" in str(stop):
                stop = tmp_bed.get_stop(start_pos[contig][max(start_pos[contig].keys())])
            orig_flank = flank
            if min(class_I_starts[contig]) -orig_flank > start:
                start = min(class_I_starts[contig]) -orig_flank
            if max(class_I_starts[contig]) + orig_flank < stop:
                stop = max(class_I_starts[contig]) + orig_flank


            if o == org: 
                if bed.get_chr(entry) == contig:
                    if start in range(xmin, xmax) and stop in range(xmin, xmax):
                        start, stop = xmin, xmax
            all_ortho_context = tmp_gff.get_region(contig, start, stop)
            for ogene in all_ortho_context:
                if tmp_gff.get_start(ogene) not in start_pos[contig]:
                    start_pos[contig][tmp_gff.get_start(ogene)] = ogene
            check_pos = list(start_pos[contig].keys())
            for s in check_pos:
                if s < start or s > stop: del start_pos[contig][s]

            if min(start_pos[contig]) in class_I_starts[contig]:
                check_for_genes = {x: 0 for x in tmp_gff.get_region(contig, 0, min(start_pos[contig]))}
                latest_start = 0
                latest_gene = "n/a"
                for check in check_for_genes:
                     if tmp_gff.get_start(check) > latest_start:
                        #print("adjusting x start for " + entry)
                        latest_start = tmp_gff.get_start(check)
                        latest_gene = check
                if latest_start > 0: start_pos[contig][latest_start] = latest_gene
            if max(start_pos[contig]) in class_I_starts[contig]:
                check_for_genes = {x: 0 for x in tmp_gff.get_region(contig, max(start_pos[contig]), orgs[o]["fai"].get_length(contig))}
                latest_stop = orgs[o]["fai"].get_length(contig)
                latest_gene = "n/a"
                for check in check_for_genes:
                     if tmp_gff.get_start(check) < latest_stop:
                        #print("adjusting x stop for " + entry)
                        latest_stop = tmp_gff.get_start(check)
                        latest_gene = check
                if latest_stop < orgs[o]["fai"].get_length(contig): start_pos[contig][latest_stop] = latest_gene

            i = 0
            for start in sorted(start_pos[contig].keys()):
                tmp_min_x = min(start_pos[contig])
                if start_pos[contig][max(start_pos[contig].keys())] in tmp_gff.parsed:
                    tmp_max_x = tmp_gff.get_stop(start_pos[contig][max(start_pos[contig].keys())])
                else:
                    tmp_max_x = tmp_bed.get_stop(start_pos[contig][max(start_pos[contig].keys())])
                ogene = start_pos[contig][start] # orthologous gene
                if ogene in tmp_gff.parsed:
                    stop = tmp_gff.get_stop(ogene)
                    col = "lightgray"
                    found = 1
                    fw = "normal"

                    if found: col = "black"
                    

                    if tmp_gff.get_strand(ogene) == "-":
                        a = "top"
                        change_y = tmp_y-3
                        rot = -35
                    else:
                        rot = +35
                        a = "bottom"
                        change_y = tmp_y+2
                    
                else:
                    stop = tmp_bed.get_stop(ogene)
                    fw = "bold"
                    col = "black"
                   
                    
                    a = "bottom"
                    color = "black"
                    if tmp_bed.get_strand(ogene) == "-":
                        change_y = tmp_y - 3
                        a = "top"
                        rot = -35
                    else:
                        rot = 35
                        change_y = tmp_y +2
                        a = "bottom"
                
                n_start, n_stop, n_max = float(start-tmp_min_x), float(stop-tmp_min_x), float(tmp_max_x-tmp_min_x)
                s_start, s_stop, s_max = n_start/n_max, n_stop/n_max, 1
                x_ticks = np.arange(0, 1, 1000.0/n_max)
                org_chrom = o + ": " + contig
                if org_chrom not in written:
                    ax.text(-0.2, tmp_y, o + ": " + contig, va="center", size=5)
                    written.append(org_chrom)
                    ax.plot([0, 1], [tmp_y, tmp_y], lw=0.5, color="darkgray")
                    for coord in x_ticks:
                        ax.plot([coord, coord], [tmp_y+1,tmp_y-1], lw = 0.5, color ="darkgray")

                rectangle = mpatch.Rectangle((s_start, change_y), s_stop-s_start, 1, fc = "black", ec = None)
                ax.add_patch(rectangle)
                ax.text(s_start, change_y, ogene, size = 5, color = col, rotation=rot, ha="left", va=a, fontweight = fw)
                for f in ortho_genes[o]:
                    if ogene in ortho_genes[o][f]:
                        connection_x = tmp_orthos[f]
                        if gff.get_strand(f) == "-":
                            connection_y = y_org - 3
                        else:
                            connection_y = y_org +2
                        ax.plot([connection_x, s_start], [connection_y, change_y], lw =0.5, color="lightgray", linestyle="dashed")
                i +=1


    ax.text(-0.2, y[org], org + ": " + contig, va="center", size=5)
    ax.plot([0, 1], [y[org], y[org]], lw=0.5, color="darkgray")

    fig.patch.set_visible(False)
    ax.axis('off')
    ax.set_ylim(y[org]-entries,50)
    ax.set_xlim(0, 1)
    plt.tight_layout()
    if number >= 1: 
        plt.savefig(entry + "_" + str(flank) + "_".join(list(orgs.keys())) + ".pdf", dpi=300, format="pdf")
    plt.close()
    

def check_org(input_org, include):
    supported_orgs = {"ddi": {"gff": ddi_gff, "bed": ddi_bed, "id": ddi_id, "fai": ddi_fai}, # more organisms can be added if added in refparser.py 
                      "dfi": {"gff": dfi_gff, "bed": dfi_bed, "id": dfi_id, "fai": dfi_fai},}
    if input_org not in supported_orgs:
        print("{} not supported..".format(input_org))
        print("supported organisms: ".format(" ".join(supported_orgs.keys())))
        
    if "all" in include:
        to_include = supported_orgs
    else:
        include += ","+input_org
        to_include = {x: supported_orgs[x] for x in include.split(",")}
    return to_include


def get_orthos(genes, org, org_dict):
    ortho_dict = {x: {} for x in org_dict}
    #ortho_dict = {x: {} for x in org_dict if x != org}
    for i in ortho_dict:
        if len(ortho_dict) > 1 and i == org: continue
        gff, id_conv = org_dict[i]["gff"], org_dict[i]["id"]
        for gene in genes:
            name     = org_dict[org]["gff"].get_name(gene)
            ortho_id = org_dict[org]["id"].get_orthoID(name)
            if ortho_id in orthologues.orthos[org]:
                ortho_genes = [x for x in orthologues.orthos[org][ortho_id][i]]
            else:
                ortho_genes = []

            ortho_dict[i][gene] = []
            for x in ortho_genes:
                if x != "":
                    orthoID = org_dict[i]["id"].get_ID(x)
                    if gene not in ortho_dict[i]:
                        ortho_dict[i][gene] = [orthoID]
                    else:
                        ortho_dict[i][gene].append(orthoID)
    return ortho_dict 


def get_regions_orthos(ortho_dict):
    ortho_class_I = {i:{} for i in ortho_dict}
    for i in ortho_dict:
        gff = org_dict[i]["gff"]
        bed = org_dict[i]["bed"]
        for gene in ortho_dict[i]:
            ortho = ortho_dict[i][gene]
            contig, start, stop = gff.get_chr(entry), gff.get_start(entry), gff.get_stop(entry)



def nearby_class_I(ortho_dict, org_dict, flank):
    ortho_class_I = {i:{} for i in ortho_dict}
    for i in ortho_dict:
        gff = org_dict[i]["gff"]
        bed = org_dict[i]["bed"]
        for gene in ortho_dict[i]:
                if len(ortho_dict[i][gene]) > 0:
                    for ortho in ortho_dict[i][gene]:
                        if gene not in ortho_class_I[i]: ortho_class_I[i][gene] = []
                        contig, start, stop = gff.get_chr(ortho), gff.get_start(ortho), gff.get_stop(ortho)
                        if "not found" in contig:
                            #print("Warning. {} not found. skipping...".format(ortho))
                            continue
                        flanking_class_I    = bed.get_region(contig, start-flank, stop+flank)
                        ortho_class_I[i][gene] += flanking_class_I
    return ortho_class_I


organism = argv[1] 
flanking = int(argv[2])
to_be_incl = argv[3] # "all" / "org" / "org1,org2,org3"
orgs     = check_org(organism, to_be_incl)
bed      = orgs[organism]["bed"]
gff      = orgs[organism]["gff"]


#print("Checking for synteny.")
for entry in bed.parsed:
    #print("parsing {} genomic information".format(entry))
    contig, start, stop = bed.get_chr(entry), bed.get_start(entry), bed.get_stop(entry)
    #print("Getting flanking class I and genes ({}bp)".format(flanking))
    flanking_class_I    = bed.get_region(contig, start-flanking, stop+flanking)
    flanking_genes      = gff.get_region(contig, start-flanking, stop+flanking)
    #print("looking for orthologues of flanking genes..")
    orthos              = get_orthos(flanking_genes, organism, orgs)
    #print("looking for class I candidates nearby orthologous genes ({}bp)".format(flanking))
    ortho_class_I       = nearby_class_I(orthos, orgs, flanking)

    plot_context(entry, contig, start, stop, flanking_class_I, flanking_genes, orthos, ortho_class_I, organism, orgs, flanking)
        

