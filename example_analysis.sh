#!/bin/bash

# requirements python 3.x, subprocess, Biopython, numpy, matplotlib, FIMO

# Create output folders 
mkdir analyses_output analyses_output/CM_search_output analyses_output/stem_mutations analyses_output/plots analyses_output/genomic_distribution \
analyses_output/synteny analyses_output/fasta analyses_output/annotations analyses_output/classifier analyses_output/CM_search_output/ddi_v1 \
analyses_output/CM_search_output/ddi_v2 analyses_output/CM_search_output/ddi_v3 analyses_output/CM_search_output/ddi_v4 \
analyses_output/CM_search_output/ddi_v5 analyses_output/CM_search_output/ddi_v6


# Parse CM search result
echo "Example of parsing of results for D. discoideum CM searches"
echo "Parsing model v1 CM output"
python scripts/parse_first_CM_search.py CM_output/ddi_v1.txt ddi fasta/dicty_chromosomal.fa
mv CM_output/ddi_v1_parsed.txt CM_output/ddi_v1_25.bed CM_output/ddi_v1_all.saf analyses_output/CM_search_output/ddi_v1/

echo "Parsing model v2 CM output"
python scripts/compare_new_search_with_old.py CM_output/ddi_v2.txt ddi v2
mv CM_output/*saf CM_output/*bed analyses_output/CM_search_output/ddi_v2/

echo "Parsing model v3 CM output"
python scripts/compare_new_search_with_old.py CM_output/ddi_v3.txt ddi v3
mv CM_output/*saf CM_output/*bed analyses_output/CM_search_output/ddi_v3/

echo "Parsing model v4 CM output"
python scripts/compare_new_search_with_old.py CM_output/ddi_v4.txt ddi v4
mv CM_output/*saf CM_output/*bed analyses_output/CM_search_output/ddi_v4/

echo "Parsing model v5 CM output"
python scripts/compare_new_search_with_old.py CM_output/ddi_v5.txt ddi v5
mv CM_output/*saf CM_output/*bed analyses_output/CM_search_output/ddi_v5/

echo "Parsing model v6 CM output"
python scripts/compare_new_search_with_old.py CM_output/ddi_v6.txt ddi v6
mv CM_output/*saf CM_output/*bed analyses_output/CM_search_output/ddi_v6/

echo "Parsing last CM search (v6 increased sensititvity)"
python scripts/compare_new_search_with_old.py CM_output/ddi_maxSens_v6.txt ddi v6Sens

echo "Identify promotor motifs upstream of CLass I candidates"
python scripts/identify_duse.py CM_output/ddi_maxSens_v6_all.bed fasta/dicty_chromosomal.fa

echo "Running Class I classifier"
python scripts/score_class_I.py CM_output/ddi_maxSens_v6.txt fimo_out/fimo.tsv ddi

mv fimo_out analyses_output/classifier/
mv *score.txt analyses_output/classifier/
mv CM_output/*saf CM_output/*bed analyses_output/annotations/
mv CM_output/*fa analyses_output/fasta/

echo "GC and length analyses of Class I RNAs"
python scripts/GC_analysis.py config.txt

python scripts/plot_class_I_lens.py config.txt

mv class_I_len_boxplot.pdf class_I_gc.txt class_I_gc.pdf analyses_output/plots/
echo "Investigating shared synteny between D. discoideum and D. firmibasis"
python scripts/synteny.py ddi 10000 dfi
mv *pdf analyses_output/synteny/

echo "Plotting genomic distribution of Class I RNA genes"
python scripts/genomic_distribution.py 
mv *pdf analyses_output/genomic_distribution/

echo "Plotting structure % of full length Class I RNAs and 11 nt motif"
python scripts/structure_motif.py
mv *pdf analyses_output/plots/

echo "Checking for compensatory mutation in stem structure"
python scripts/stem_summary_no_flank_fasta.py fasta/ddi_curated.fa 6 stem_mutations_ddi.txt
python scripts/stem_summary_no_flank_fasta.py fasta/dfi_curated.fa 6 stem_mutations_dfi.txt
mv stem_mutations* analyses_output/stem_mutations/

echo "Parsing Class I RNA fasta with flanking sequences"
python scripts/get_fasta_w_flank.py config.txt

echo "Removing variable region from Class I RNAs"
python scripts/filter_fastas.py ddi_curated_flank.fa
mv *fa analyses_output/fasta

echo "Results can be found in analyses_output"