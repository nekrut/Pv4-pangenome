#!/usr/bin/env python3
"""
Phase E — rbest chain overlap edges.

For each rbest chain file, identify gene pairs (strain_a.gene_a ↔ strain_b.gene_b)
where the chain overlaps both gene bodies by >= min_overlap (reciprocal).

Output TSV columns: strain_a, gene_a, strain_b, gene_b, overlap_a, overlap_b

Requires: Python 3.9+
"""

import argparse
import csv
import glob
import sys
from collections import defaultdict
from pathlib import Path


def parse_chain_header(line):
    """Parse a UCSC chain header line into a dict."""
    parts = line.strip().split()
    if len(parts) < 13 or parts[0] != 'chain':
        return None
    return {
        'score':   int(parts[1]),
        'tName':   parts[2],
        'tSize':   int(parts[3]),
        'tStrand': parts[4],
        'tStart':  int(parts[5]),
        'tEnd':    int(parts[6]),
        'qName':   parts[7],
        'qSize':   int(parts[8]),
        'qStrand': parts[9],
        'qStart':  int(parts[10]),
        'qEnd':    int(parts[11]),
        'id':      parts[12],
    }


def load_bed_genes(pattern):
    """Load genes from all BED files matching pattern.
    Returns dict strain -> list of (chrom, start, end, gene_id).
    """
    genes: dict = defaultdict(list)
    for bed_path in glob.glob(pattern):
        strain = Path(bed_path).stem   # filename without extension = strain name
        with open(bed_path) as f:
            for ln in f:
                if ln.startswith('#') or not ln.strip():
                    continue
                parts = ln.rstrip('\n').split('\t')
                if len(parts) < 4:
                    continue
                chrom, start, end, gid = parts[0], int(parts[1]), int(parts[2]), parts[3]
                genes[strain].append((chrom, start, end, gid))
    return genes


def reciprocal_overlap(cs, ce, gs, ge):
    """Overlap fraction between chain segment [cs,ce) and gene [gs,ge)."""
    ov = max(0, min(ce, ge) - max(cs, gs))
    if ov == 0:
        return 0.0
    gene_len = max(1, ge - gs)
    return ov / gene_len


def genes_overlapping(genes_list, chrom, start, end, min_ov):
    """Return list of (gene_id, overlap) for genes with overlap >= min_ov."""
    hits = []
    for g_chrom, g_start, g_end, g_id in genes_list:
        if g_chrom != chrom:
            continue
        ov = reciprocal_overlap(start, end, g_start, g_end)
        if ov >= min_ov:
            hits.append((g_id, ov))
    return hits


def main():
    ap = argparse.ArgumentParser(
        description="Phase E: extract gene-pair edges from rbest chain files")
    ap.add_argument('--chains', required=True,
                    help='Glob pattern for rbest chain files (quote it)')
    ap.add_argument('--annotations', required=True,
                    help='Glob pattern for per-strain BED files (quote it)')
    ap.add_argument('--strains', required=True,
                    help='Space-separated strain list')
    ap.add_argument('--min_overlap', type=float, default=0.90)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    all_strains = args.strains.split()
    genes_by_strain = load_bed_genes(args.annotations)
    print(f"Loaded genes for strains: {list(genes_by_strain.keys())}", file=sys.stderr)

    chain_files = glob.glob(args.chains)
    print(f"Processing {len(chain_files)} rbest chain files...", file=sys.stderr)

    edges = []
    for chain_path in chain_files:
        # filename: {A}.{B}.rbest.chain
        stem = Path(chain_path).stem.replace('.rbest', '')
        parts = stem.split('.')
        if len(parts) < 2:
            continue
        strain_a, strain_b = parts[0], parts[1]

        genes_a = genes_by_strain.get(strain_a, [])
        genes_b = genes_by_strain.get(strain_b, [])

        with open(chain_path) as f:
            cur_header = None
            for ln in f:
                ln = ln.rstrip('\n')
                if ln.startswith('chain'):
                    cur_header = parse_chain_header(ln)
                elif cur_header and ln:
                    pass   # chain data line; we use the header coords for gene overlap

            # Use header-level coords for gene-pair matching
            # (block-level would be more precise but expensive for this stub)

        with open(chain_path) as f:
            for ln in f:
                ln = ln.rstrip('\n')
                if not ln.startswith('chain'):
                    continue
                h = parse_chain_header(ln)
                if not h:
                    continue
                t_hits = genes_overlapping(
                    genes_a, h['tName'], h['tStart'], h['tEnd'], args.min_overlap)
                q_hits = genes_overlapping(
                    genes_b, h['qName'], h['qStart'], h['qEnd'], args.min_overlap)
                for g_a, ov_a in t_hits:
                    for g_b, ov_b in q_hits:
                        edges.append({
                            'strain_a': strain_a, 'gene_a': g_a,
                            'strain_b': strain_b, 'gene_b': g_b,
                            'overlap_a': f'{ov_a:.3f}',
                            'overlap_b': f'{ov_b:.3f}',
                        })

    print(f"  {len(edges)} rbest edges found", file=sys.stderr)
    fields = ['strain_a', 'gene_a', 'strain_b', 'gene_b', 'overlap_a', 'overlap_b']
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(edges)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
