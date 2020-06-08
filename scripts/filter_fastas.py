#!/usr/bin/env python
from sys import argv

low  = 100
high = 0

filtered = {}

in_fasta = argv[1]
print(in_fasta)
with open(in_fasta) as fin:
	for line in fin:
		if line[0] == ">":
			seq = next(fin)
			fSeq = seq[:20] + "N"*10 + seq[-30:]
			removed = len(seq) - len(fSeq)
			if removed >= high: high = removed
			if removed <= low: low = removed
			filtered[line] = fSeq

print("Length range of removed variable region.")
print("min: " + str(low))
print("max: " + str(high))
with open(in_fasta.replace("fa", "filtered.fa"), "w") as fout:
	for entry in filtered:
		fout.write(entry)
		fout.write(filtered[entry])
