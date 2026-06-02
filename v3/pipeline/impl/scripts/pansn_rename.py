#!/usr/bin/env python3
"""
PanSN rename: rewrite FASTA headers to SAMPLE#1#CONTIG format.
Usage: pansn_rename.py <fasta_path> <sample_name>
Writes renamed FASTA to stdout.
"""
import sys

if len(sys.argv) != 3:
    print("Usage: pansn_rename.py <fasta> <sample>", file=sys.stderr)
    sys.exit(1)

in_fa, sample = sys.argv[1], sys.argv[2]
for line in open(in_fa):
    if line.startswith('>'):
        contig = line[1:].strip().split()[0]
        print(f'>{sample}#1#{contig}')
    else:
        print(line, end='')
