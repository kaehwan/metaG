#!/usr/bin/env python3

""" 
USAGE: ${SOFT}/rm_short_contigs.py ${out}/metaspades/scaffolds.fasta $min_len > ${out}/final_assembly.fasta
"""
import sys

# example metaspades contig/scaffolds: ">NODE_1_length_440246_cov_10.775182"
for line in open(sys.argv[1]):
    if not line.startswith(">"):
        print(line.strip())
    else:
        if int(line.strip().split("_")[3]) < int(sys.argv[2]):
            break
        else:
            print(line.strip())
