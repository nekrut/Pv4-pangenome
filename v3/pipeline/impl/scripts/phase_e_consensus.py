#!/usr/bin/env python3
"""
Phase E — consensus ortholog table for the Pk v1 pipeline.

Builds a multigraph from per-anchor Phase C.4 classification.tsv files:
  - Each (strain, gene_id) pair is a graph node.
  - Edges connect orthologous genes; weight by source/intactness.
  - Connected components → orthogroups.
  - Labels: CORE-1:1 | CORE-VAR | FAMILY | PARTIAL | LINEAGE-SPECIFIC

Output: work/03_consensus/ortholog_table.tsv
  orthogroup_id, label, n_strains, max_copies, {strain columns...}

Ported and parameterized from
  /media/anton/data/sandbox/Pv4/v3/scripts/phase_e_consensus.py
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


WEIGHTS = {
    ('cesar2', 'I'):  1.00,
    ('cesar2', 'PI'): 0.70,
    ('cesar2', 'UL'): 0.40,
    ('cesar2', 'PG'): 0.40,
    ('cesar2', 'L'):  0.00,
    ('cesar2', 'M'):  0.00,
    ('cesar2', 'FI'): 0.20,
    ('liftoff', 'I'): 0.95,
}


def edge_weight(source, intactness):
    if source == 'liftoff':
        return 0.95
    if source == 'none':
        return 0.0
    return WEIGHTS.get((source, intactness), 0.10)


def normalize_gene_id(gid):
    if not gid or gid == 'None':
        return gid
    # _tN PlasmoDB transcript suffix
    m = re.match(r'^(.+)_t\d+$', gid)
    if m:
        return m.group(1)
    # .N transcript suffix
    m = re.match(r'^(.+)\.\d+$', gid)
    if m:
        return m.group(1)
    # _N Liftoff extra-copy suffix (small integer)
    m = re.match(r'^(.+)_(\d+)$', gid)
    if m and len(m.group(2)) <= 2 and not m.group(1).endswith('_'):
        return m.group(1)
    return gid


class UnionFind:
    def __init__(self):
        self.parent: dict = {}

    def find(self, x):
        while self.parent.get(x, x) != x:
            self.parent[x] = self.parent.get(self.parent[x], self.parent[x])
            x = self.parent[x]
        self.parent.setdefault(x, x)
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def load_classifications(liftoff_dir, anchors, all_strains):
    """Yield (anchor, query, ref_gene_canonical, q_gene_canonical, source, intactness, weight)."""
    base = Path(liftoff_dir)
    for anchor in anchors:
        sub = base / f'{anchor}-as-ref'
        if not sub.exists():
            continue
        for q in all_strains:
            if q == anchor:
                continue
            cls_path = sub / f'{q}.classification.tsv'
            if not cls_path.exists():
                continue
            with open(cls_path) as fh:
                r = csv.DictReader(fh, delimiter='\t')
                for row in r:
                    rg = row.get('reference_gene_id', '')
                    qg = row.get('query_gene_id', '')
                    src = row.get('source', '')
                    intact = row.get('intactness', '')
                    if not qg or qg == 'None':
                        continue
                    w = edge_weight(src, intact)
                    if w == 0:
                        continue
                    yield anchor, q, rg, qg, src, intact, w


def reciprocal_overlap(a, b):
    s1, e1, _ = a
    s2, e2, _ = b
    ov = max(0, min(e1, e2) - max(s1, s2))
    if ov == 0:
        return 0.0
    return min(ov / max(1, e1 - s1), ov / max(1, e2 - s2))


def collapse_positions(gene_id_list, strain, node_pos):
    """Group gene IDs by overlapping query position; return list of position-clusters."""
    recs = [(gid, node_pos.get(f'{strain}#{gid}')) for gid in gene_id_list]
    positioned = [(p, gid) for gid, p in recs if p is not None]
    no_pos = [gid for gid, p in recs if p is None]
    if not positioned:
        return [[g] for g in dict.fromkeys(no_pos)]
    by_chr: dict = defaultdict(list)
    for (chrom, s, e), gid in positioned:
        by_chr[chrom].append((s, e, gid))
    clusters = []
    for chrom, lst in by_chr.items():
        lst.sort()
        used = [False] * len(lst)
        for i in range(len(lst)):
            if used[i]:
                continue
            s_i, e_i, g_i = lst[i]
            grp = [g_i]
            used[i] = True
            for j in range(i + 1, len(lst)):
                if used[j]:
                    continue
                s_j, e_j, g_j = lst[j]
                if s_j > e_i:
                    break
                ov = max(0, min(e_i, e_j) - max(s_i, s_j))
                if ov and min(ov / (e_i - s_i + 1), ov / (e_j - s_j + 1)) >= 0.2:
                    grp.append(g_j)
                    used[j] = True
            clusters.append(grp)
    for g in dict.fromkeys(no_pos):
        clusters.append([g])
    return clusters


def main():
    ap = argparse.ArgumentParser(description="Phase E: consensus ortholog table")
    ap.add_argument('--liftoff_dir', required=True, help='work/02d_merged/')
    ap.add_argument('--rbest', required=True, help='work/03_consensus/rbest_edges.tsv')
    ap.add_argument('--graph', required=True, help='work/03_consensus/graph_edges.tsv')
    ap.add_argument('--anchors', required=True, help='Space-separated anchor strain list')
    ap.add_argument('--strains', required=True, help='Space-separated all-strain list')
    ap.add_argument('--ref', required=True, help='Reference strain name')
    ap.add_argument('--output', required=True, help='Output TSV path')
    args = ap.parse_args()

    anchors = args.anchors.split()
    all_strains = args.strains.split()

    print('Loading per-anchor classifications...', flush=True)
    raw = list(load_classifications(args.liftoff_dir, anchors, all_strains))
    print(f'  raw edges: {len(raw)}')

    uf = UnionFind()
    pos_records: dict = defaultdict(list)
    node_pos: dict = {}

    # Seed positions from classification.tsv query_chrom/start/end columns
    base = Path(args.liftoff_dir)
    for anchor in anchors:
        sub = base / f'{anchor}-as-ref'
        if not sub.exists():
            continue
        for q in all_strains:
            if q == anchor:
                continue
            cls_path = sub / f'{q}.classification.tsv'
            if not cls_path.exists():
                continue
            with open(cls_path) as fh:
                r = csv.DictReader(fh, delimiter='\t')
                for row in r:
                    rg = normalize_gene_id(row.get('reference_gene_id', ''))
                    qg = normalize_gene_id(row.get('query_gene_id', ''))
                    src = row.get('source', '')
                    intact = row.get('intactness', '')
                    if not qg or qg == 'None':
                        continue
                    if src == 'none' or edge_weight(src, intact) == 0:
                        continue
                    a_node = f'{anchor}#{rg}'
                    q_node = f'{q}#{qg}'
                    uf.union(a_node, q_node)
                    chrom = row.get('query_chrom', '')
                    start_s = row.get('query_start', '')
                    end_s = row.get('query_end', '')
                    if chrom and start_s and end_s:
                        try:
                            s_i, e_i = int(start_s), int(end_s)
                            pos_records[(q, chrom)].append((s_i, e_i, q_node))
                            node_pos[q_node] = (chrom, s_i, e_i)
                        except ValueError:
                            pass

    # Interval-based aliasing (90% reciprocal overlap → same physical gene)
    aliases_merged = 0
    for key, recs in pos_records.items():
        recs.sort()
        for i in range(len(recs)):
            si, ei, ni = recs[i]
            for j in range(i + 1, len(recs)):
                sj, ej, nj = recs[j]
                if sj > ei:
                    break
                if uf.find(ni) == uf.find(nj):
                    continue
                if reciprocal_overlap(recs[i], recs[j]) >= 0.9:
                    uf.union(ni, nj)
                    aliases_merged += 1
    print(f'  position aliases merged: {aliases_merged}')

    # Incorporate rbest chain edges
    rbest_edges_added = 0
    if Path(args.rbest).exists():
        with open(args.rbest) as fh:
            r = csv.DictReader(fh, delimiter='\t')
            for row in r:
                sa = row.get('strain_a', '')
                ga = normalize_gene_id(row.get('gene_a', ''))
                sb = row.get('strain_b', '')
                gb = normalize_gene_id(row.get('gene_b', ''))
                if sa and ga and sb and gb:
                    uf.union(f'{sa}#{ga}', f'{sb}#{gb}')
                    rbest_edges_added += 1
    print(f'  rbest chain edges: {rbest_edges_added}')

    # Incorporate graph co-membership edges
    graph_edges_added = 0
    if Path(args.graph).exists():
        with open(args.graph) as fh:
            r = csv.DictReader(fh, delimiter='\t')
            for row in r:
                sa = row.get('strain_a', '')
                ga = normalize_gene_id(row.get('gene_a', ''))
                sb = row.get('strain_b', '')
                gb = normalize_gene_id(row.get('gene_b', ''))
                if sa and ga and sb and gb:
                    uf.union(f'{sa}#{ga}', f'{sb}#{gb}')
                    graph_edges_added += 1
    print(f'  graph co-membership edges: {graph_edges_added}')

    # Connected components
    comps: dict = defaultdict(set)
    for node in uf.parent:
        comps[uf.find(node)].add(node)
    print(f'  orthogroup count (all connected components): {len(comps)}')

    N_ALL = len(all_strains)
    rows_out = []
    for cid, nodes in comps.items():
        per_strain: dict = defaultdict(list)
        for n in nodes:
            parts = n.split('#', 1)
            if len(parts) == 2:
                per_strain[parts[0]].append(parts[1])
        if len(per_strain) < 2:
            continue
        strain_clusters = {s: collapse_positions(gs, s, node_pos) for s, gs in per_strain.items()}
        present_strains = [s for s in all_strains if s in per_strain]
        n_strains = len(present_strains)
        max_copies = max(len(c) for c in strain_clusters.values())
        if n_strains == N_ALL and max_copies == 1:
            label = 'CORE-1:1'
        elif n_strains == N_ALL and max_copies >= 2:
            label = 'CORE-VAR'
        elif max_copies >= 3:
            label = 'FAMILY'
        elif n_strains <= 2:
            label = 'LINEAGE-SPECIFIC'
        else:
            label = 'PARTIAL'
        row = {
            'orthogroup_id': f'OG{len(rows_out)+1:06d}',
            'label': label,
            'n_strains': n_strains,
            'max_copies': max_copies,
        }
        for s in all_strains:
            if s in strain_clusters:
                row[s] = ','.join('|'.join(c) for c in strain_clusters[s])
            else:
                row[s] = '-'
        rows_out.append(row)

    print(f'  multi-strain orthogroups: {len(rows_out)}')
    labels = Counter(r['label'] for r in rows_out)
    for k, n in labels.most_common():
        print(f'    {k}: {n}')

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ['orthogroup_id', 'label', 'n_strains', 'max_copies'] + all_strains
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows_out)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
