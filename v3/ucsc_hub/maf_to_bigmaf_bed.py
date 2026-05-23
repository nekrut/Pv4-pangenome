#!/usr/bin/env python3
"""
Convert a preprocessed MAF to BED3+1 directly (bypass mafToBigMaf's overlap check).
The 4th field is the full MAF block text with newlines replaced by ';'.

Usage: maf_to_bigmaf_bed.py <ref_acc> <input.maf> <output.bed>
"""
import sys

ref_acc, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]


def species_of(seq_name):
    parts = seq_name.split('.')
    if len(parts) >= 2 and parts[0].startswith(('GCA_', 'GCF_')):
        return parts[0] + '.' + parts[1]
    return parts[0]


def emit_block(out, block):
    """Block is a list of MAF lines (including the 'a ' line). Find the ref s-line, then emit BED3+1."""
    ref_chrom = None
    ref_start = None
    ref_size = None
    for line in block:
        if line.startswith('s '):
            parts = line.split()
            sname = parts[1]
            if species_of(sname) == ref_acc:
                # Strip species prefix to match chrom-sizes file naming
                ref_chrom = sname[len(ref_acc) + 1:] if sname.startswith(ref_acc + '.') else sname
                ref_start = int(parts[2])
                ref_size = int(parts[3])
                break
    if ref_chrom is None:
        return
    # MAF block text: join lines with ';' (UCSC bigMaf convention)
    # Strip trailing newlines, trim leading whitespace per line
    text_lines = [l.rstrip('\n') for l in block]
    block_text = ';'.join(text_lines)
    out.write(f"{ref_chrom}\t{ref_start}\t{ref_start + ref_size}\t{block_text}\n")


with open(src) as fh, open(dst, 'w') as out:
    block = []
    in_body = False
    n_emit = 0
    for line in fh:
        if not in_body:
            if line.startswith('##maf'):
                in_body = True
            continue
        if line.startswith('a '):
            if block:
                emit_block(out, block)
                n_emit += 1
            block = [line]
        elif line.strip() == '':
            if block:
                emit_block(out, block)
                n_emit += 1
                block = []
        else:
            if block:
                block.append(line)
    if block:
        emit_block(out, block)
        n_emit += 1

print(f"  Emitted {n_emit} BED3+1 records to {dst}", flush=True)
