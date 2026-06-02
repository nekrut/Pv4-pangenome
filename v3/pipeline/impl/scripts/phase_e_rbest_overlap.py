#!/usr/bin/env python3
"""
Phase E — rbest chain overlap edges (block-level projection).

For each rbest chain, project every target gene's interval through the chain's
aligned BLOCKS into query coordinates, then connect it to the single query gene
that best overlaps the projected location. An edge is emitted only when the
target gene is well-covered by aligned blocks (overlap_a >= min_overlap) AND its
best query gene reciprocally overlaps the projection (overlap_b >= min_overlap).

This replaces the earlier stub that took the cartesian product of all genes under
a chain's whole-chromosome header span (which produced millions of spurious edges).

Output TSV columns: strain_a, gene_a, strain_b, gene_b, overlap_a, overlap_b
Requires: Python 3.9+
"""

import argparse
import csv
import glob
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path


def parse_chain_header(line):
    p = line.strip().split()
    if len(p) < 13 or p[0] != 'chain':
        return None
    return {'tName': p[2], 'tSize': int(p[3]), 'tStart': int(p[5]), 'tEnd': int(p[6]),
            'qName': p[7], 'qSize': int(p[8]), 'qStrand': p[9],
            'qStart': int(p[10]), 'qEnd': int(p[11])}


def load_bed_genes(pattern):
    genes = defaultdict(list)
    for bed_path in glob.glob(pattern):
        strain = Path(bed_path).stem
        with open(bed_path) as f:
            for ln in f:
                if ln.startswith('#') or not ln.strip():
                    continue
                parts = ln.rstrip('\n').split('\t')
                if len(parts) < 4:
                    continue
                genes[strain].append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return genes


def index_by_chrom(gene_list):
    """chrom -> (sorted_starts, sorted list of (start,end,gid)) for bisect lookup."""
    by = defaultdict(list)
    for chrom, s, e, gid in gene_list:
        by[chrom].append((s, e, gid))
    idx = {}
    for chrom, lst in by.items():
        lst.sort()
        idx[chrom] = ([g[0] for g in lst], lst)
    return idx


def iter_chains(path):
    """Yield (header, blocks) where blocks are (t0,t1,qf0,qf1,strand) in forward query coords."""
    h = None
    blocks = []
    t = q = 0
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip('\n')
            if ln.startswith('chain'):
                if h:
                    yield h, blocks
                h = parse_chain_header(ln)
                blocks = []
                if not h:
                    continue
                t = h['tStart']
                q = h['qStart']
            elif h and ln.strip():
                f3 = ln.split()
                size = int(f3[0])
                if h['qStrand'] == '+':
                    qf0, qf1 = q, q + size
                else:
                    qf1 = h['qSize'] - q
                    qf0 = h['qSize'] - (q + size)
                blocks.append((t, t + size, qf0, qf1, h['qStrand']))
                if len(f3) >= 3:
                    t += size + int(f3[1])
                    q += size + int(f3[2])
        if h:
            yield h, blocks


def project_gene(gs, ge, blocks):
    """Project target gene [gs,ge) through blocks -> (aligned_len, qmin, qmax)."""
    aligned = 0
    qmin = qmax = None
    for t0, t1, qf0, qf1, strand in blocks:
        lo, hi = max(gs, t0), min(ge, t1)
        if hi <= lo:
            continue
        aligned += hi - lo
        if strand == '+':
            a, b = qf0 + (lo - t0), qf0 + (hi - t0)
        else:
            a, b = qf1 - (hi - t0), qf1 - (lo - t0)
        lo_q, hi_q = (a, b) if a <= b else (b, a)
        qmin = lo_q if qmin is None else min(qmin, lo_q)
        qmax = hi_q if qmax is None else max(qmax, hi_q)
    return aligned, qmin, qmax


def best_query_gene(qmin, qmax, qchrom_idx):
    """Best (gid, overlap_b) among query genes overlapping [qmin,qmax)."""
    if qchrom_idx is None or qmin is None:
        return None
    starts, lst = qchrom_idx
    best = None
    hi = bisect_right(starts, qmax)
    for i in range(hi):
        gs, ge, gid = lst[i]
        if ge <= qmin:
            continue
        ov = min(qmax, ge) - max(qmin, gs)
        if ov <= 0:
            continue
        frac = ov / max(1, ge - gs)
        if best is None or frac > best[1]:
            best = (gid, frac)
    return best


def main():
    ap = argparse.ArgumentParser(description="Phase E: rbest chain block-level gene edges")
    ap.add_argument('--chains', required=True)
    ap.add_argument('--annotations', required=True)
    ap.add_argument('--strains', required=True)
    ap.add_argument('--min_overlap', type=float, default=0.90)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    genes_by_strain = load_bed_genes(args.annotations)
    print(f"Loaded genes for strains: {list(genes_by_strain.keys())}", file=sys.stderr)
    chain_files = glob.glob(args.chains)
    print(f"Processing {len(chain_files)} rbest chain files...", file=sys.stderr)

    edges = []
    seen = set()
    for chain_path in chain_files:
        stem = Path(chain_path).stem.replace('.rbest', '')
        parts = stem.split('.')
        if len(parts) < 2:
            continue
        strain_a, strain_b = parts[0], parts[1]
        q_idx = index_by_chrom(genes_by_strain.get(strain_b, []))
        a_by_chrom = defaultdict(list)
        for g in genes_by_strain.get(strain_a, []):
            a_by_chrom[g[0]].append(g)

        for h, blocks in iter_chains(chain_path):
            if not blocks:
                continue
            for chrom, s, e, gid in a_by_chrom.get(h['tName'], []):
                if e <= h['tStart'] or s >= h['tEnd']:
                    continue
                aligned, qmin, qmax = project_gene(s, e, blocks)
                ov_a = aligned / max(1, e - s)
                if ov_a < args.min_overlap or qmin is None:
                    continue
                bq = best_query_gene(qmin, qmax, q_idx.get(h['qName']))
                if not bq or bq[1] < args.min_overlap:
                    continue
                key = (strain_a, gid, strain_b, bq[0])
                if key in seen:
                    continue
                seen.add(key)
                edges.append({'strain_a': strain_a, 'gene_a': gid,
                              'strain_b': strain_b, 'gene_b': bq[0],
                              'overlap_a': f'{ov_a:.3f}', 'overlap_b': f'{bq[1]:.3f}'})

    print(f"  {len(edges)} rbest edges found", file=sys.stderr)
    fields = ['strain_a', 'gene_a', 'strain_b', 'gene_b', 'overlap_a', 'overlap_b']
    with open(args.output, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(edges)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
