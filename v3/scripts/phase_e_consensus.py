#!/usr/bin/env python3
"""
Phase E — consensus ortholog table.

Build the multigraph from per-anchor Phase D classification.tsv files:
  - Each strain×gene_id is a node.
  - Edges connect orthologous gene_ids; weight by source/intactness:
       CESAR2-I  → 1.00
       Liftoff   → 0.95
       CESAR2-PI → 0.70
       UL/PG     → 0.40
       L/M       → 0
  - Connected components after thresholding give orthogroups.
  - Label each by 8-strain occupancy:
       CORE-1:1            : exactly 1 gene per strain in all 8 strains
       CORE-VAR            : present in all 8 but multi-copy in some
       FAMILY              : multi-copy in ≥4 strains (or contains an R8/family gene)
       LINEAGE-SPECIFIC    : present in 1-2 strains
       LOST                : present elsewhere but absent in ≥1 strain

Output: work/03_consensus/ortholog_table.tsv with columns:
  orthogroup_id  label  PvP01  Sal-I  PvW1  PAM  PvSY56  PvT01  PvC01  MHC087  provenance
"""
import csv
import os
from pathlib import Path
from collections import defaultdict


ANCHORS = ['PvP01', 'PvW1', 'PAM', 'PvSY56']
QUERIES_PER_ANCHOR = {
    'PvP01':  ['Sal-I', 'PvT01', 'PvC01', 'PvW1', 'PAM', 'PvSY56', 'MHC087'],
    'PvW1':   ['PvP01', 'Sal-I', 'PvT01', 'PvC01', 'PAM', 'PvSY56', 'MHC087'],
    'PAM':    ['PvP01', 'Sal-I', 'PvT01', 'PvC01', 'PvW1', 'PvSY56', 'MHC087'],
    'PvSY56': ['PvP01', 'Sal-I', 'PvT01', 'PvC01', 'PvW1', 'PAM', 'MHC087'],
}
ALL_STRAINS = ['PvP01', 'Sal-I', 'PvW1', 'PAM', 'PvSY56', 'PvT01', 'PvC01', 'MHC087']

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
    """Strip transcript/extra-copy suffixes so the same gene is one ID across anchors.

    Liftoff convention: 'PVP01_xxx.1' → 'PVP01_xxx'  (.N)
    PlasmoDB PvW1/PAM:  'PVW1_xxx_t1' → 'PVW1_xxx'    (_tN)
    Liftoff extra copy: 'PVP01_xxx_1' → 'PVP01_xxx'   (_N for small N)
    """
    import re
    if not gid or gid == 'None':
        return gid
    # Try the _tN suffix first (PlasmoDB transcripts)
    m = re.match(r'^(.+)_t\d+$', gid)
    if m:
        return m.group(1)
    # Then .N
    m = re.match(r'^(.+)\.\d+$', gid)
    if m:
        return m.group(1)
    # Then _N (small number — Liftoff extra-copy)
    m = re.match(r'^(.+)_(\d+)$', gid)
    if m and len(m.group(2)) <= 2 and not m.group(1).endswith('_'):
        return m.group(1)
    return gid


def load_classifications():
    """Yield (anchor, query, ref_gene, q_gene, source, intactness, weight)."""
    base = Path('work/02d_merged')
    for anchor in ANCHORS:
        sub = base / f'{anchor}-as-ref'
        if not sub.exists():
            continue
        for q in QUERIES_PER_ANCHOR[anchor]:
            cls_path = sub / f'{q}.classification.tsv'
            if not cls_path.exists():
                continue
            with open(cls_path) as fh:
                r = csv.DictReader(fh, delimiter='\t')
                for row in r:
                    rg = row['reference_gene_id']
                    qg = row['query_gene_id']
                    src = row['source']
                    intact = row['intactness']
                    if not qg or qg == 'None':
                        continue
                    w = edge_weight(src, intact)
                    if w == 0:
                        continue
                    yield anchor, q, rg, qg, src, intact, w


class UnionFind:
    def __init__(self):
        self.parent = {}
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


def main():
    print('Loading per-anchor classifications...', flush=True)
    edges = []  # (anchor, query, ref_gene_canonical, q_gene_canonical, weight, src)
    raw = list(load_classifications())
    print(f'  raw edges: {len(raw)}')

    # Build PvP01-canonical gene ID for each (anchor, gene_id)
    # The classification.tsv rows from anchor=A list ref_gene_id in A's gene-ID space.
    # To anchor everything in PvP01-space, we'll add (anchor#gene_id) tokens to the union-find,
    # then merge across anchors via the same query gene IDs.
    # Node format: 'strain#gene_id'
    # For 'PvP01' (anchor or query), strain prefix = 'PvP01'
    uf = UnionFind()
    weight_accum = defaultdict(float)
    src_accum = defaultdict(set)

    # Position-based aliasing: collect (strain, chrom, start, end) → set of nodes.
    # Then merge nodes that have overlapping intervals on the same strain+chrom.
    pos_records = defaultdict(list)  # (strain, chrom) → list of (start, end, node)
    node_pos = {}  # node → (chrom, start, end), used at labeling time

    # Pre-seed node_pos and pos_records from each strain's native BED so that
    # ref_gene_ids in anchor rows (which classification.tsv does NOT give us coords for)
    # have known positions in their home genome.
    native_bed = {
        'PvP01':  'inputs/annotations/PvP01.bed',
        'PvW1':   'inputs/annotations/PvW1.bed',
        'PAM':    'inputs/annotations/PAM.bed',
        'PvSY56': 'inputs/annotations/PvSY56.bed',
    }
    for strain, bed in native_bed.items():
        if not Path(bed).exists():
            continue
        with open(bed) as fh:
            for ln in fh:
                f = ln.rstrip('\n').split('\t')
                if len(f) < 4: continue
                chrom, s, e, gid = f[0], int(f[1]), int(f[2]), f[3]
                gid = normalize_gene_id(gid)
                node = f'{strain}#{gid}'
                node_pos[node] = (chrom, s, e)
                pos_records[(strain, chrom)].append((s, e, node))

    # Read row coords from classification.tsv (need them again — repeat the loop with full rows)
    base = Path('work/02d_merged')
    for anchor in ANCHORS:
        sub = base / f'{anchor}-as-ref'
        if not sub.exists():
            continue
        for q in QUERIES_PER_ANCHOR[anchor]:
            cls_path = sub / f'{q}.classification.tsv'
            if not cls_path.exists():
                continue
            with open(cls_path) as fh:
                r = csv.DictReader(fh, delimiter='\t')
                for row in r:
                    rg = normalize_gene_id(row['reference_gene_id'])
                    qg = normalize_gene_id(row['query_gene_id'])
                    if not qg or qg == 'None':
                        continue
                    src = row['source']; intact = row['intactness']
                    if src == 'none' or edge_weight(src, intact) == 0:
                        continue
                    a_node = f'{anchor}#{rg}'
                    q_node = f'{q}#{qg}'
                    uf.union(a_node, q_node)
                    # Record query position for later interval-based merge
                    chrom = row.get('query_chrom', '')
                    start = row.get('query_start', '')
                    end = row.get('query_end', '')
                    if chrom and start and end:
                        try:
                            s_i, e_i = int(start), int(end)
                            pos_records[(q, chrom)].append((s_i, e_i, q_node))
                            node_pos[q_node] = (chrom, s_i, e_i)
                        except ValueError:
                            pass

    # Interval-based aliasing: two nodes on same (strain, chrom) with >=50% reciprocal overlap
    # are the same physical gene under different ID conventions → merge.
    def reciprocal_overlap(a, b):
        s1, e1, _ = a; s2, e2, _ = b
        ov = max(0, min(e1, e2) - max(s1, s2))
        if ov == 0: return 0
        return min(ov / max(1, e1 - s1), ov / max(1, e2 - s2))

    aliases_merged = 0
    for key, recs in pos_records.items():
        # Sort by start; sweep with a small look-back window
        recs.sort()
        for i in range(len(recs)):
            si, ei, ni = recs[i]
            for j in range(i + 1, len(recs)):
                sj, ej, nj = recs[j]
                if sj > ei:  # past the right edge — sorted by start so no further overlap
                    break
                if uf.find(ni) == uf.find(nj):
                    continue  # already merged
                if reciprocal_overlap(recs[i], recs[j]) >= 0.9:
                    uf.union(ni, nj)
                    aliases_merged += 1
    print(f'  position aliases merged: {aliases_merged}')

    # Legacy edge-weight tracking from the original raw list
    for anchor, q, ref_gene, q_gene, src, intact, w in raw:
        ref_gene = normalize_gene_id(ref_gene); q_gene = normalize_gene_id(q_gene)
        a_node = f'{anchor}#{ref_gene}'; q_node = f'{q}#{q_gene}'
        e = tuple(sorted([a_node, q_node]))
        weight_accum[e] += w
        src_accum[e].add(f'{src}:{intact}({anchor}→{q})')

    # Connected components
    comps = defaultdict(set)
    for node in uf.parent:
        comps[uf.find(node)].add(node)
    print(f'  orthogroup count (connected components, all edges): {len(comps)}')

    # For each component, build the per-strain occupancy.
    # Count COPIES by distinct physical position, not distinct gene_id —
    # otherwise cross-anchor aliases (PVPAM_xxx vs PVP01_xxx for the same physical
    # gene) inflate the copy count.
    def collapse_positions(gene_id_list, strain):
        """Group gene IDs by overlapping query position; return list of position-clusters."""
        recs = [(gid, node_pos.get(f'{strain}#{gid}')) for gid in gene_id_list]
        # Genes without position info → each unique gid is its own cluster
        positioned = [(p, gid) for gid, p in recs if p is not None]
        no_pos    = [gid for gid, p in recs if p is None]
        if not positioned:
            return [[g] for g in dict.fromkeys(no_pos)]
        # Group by chrom, then sweep with ≥90% reciprocal overlap
        from collections import defaultdict as _dd
        by_chr = _dd(list)
        for (chrom, s, e), gid in positioned:
            by_chr[chrom].append((s, e, gid))
        clusters = []
        for chrom, lst in by_chr.items():
            lst.sort()
            used = [False] * len(lst)
            for i in range(len(lst)):
                if used[i]: continue
                s_i, e_i, g_i = lst[i]
                grp = [g_i]; used[i] = True
                for j in range(i + 1, len(lst)):
                    if used[j]: continue
                    s_j, e_j, g_j = lst[j]
                    if s_j > e_i: break
                    ov = max(0, min(e_i, e_j) - max(s_i, s_j))
                    if ov and min(ov/(e_i-s_i+1), ov/(e_j-s_j+1)) >= 0.2:
                        grp.append(g_j); used[j] = True
                clusters.append(grp)
        # Append no_pos as singleton clusters
        for g in dict.fromkeys(no_pos):
            clusters.append([g])
        return clusters

    rows = []
    for cid, nodes in comps.items():
        per_strain = defaultdict(list)
        for n in nodes:
            strain, gid = n.split('#', 1)
            per_strain[strain].append(gid)
        if len(per_strain) < 2:
            continue
        # Collapse per strain by physical position
        strain_clusters = {s: collapse_positions(gs, s) for s, gs in per_strain.items()}
        present_strains = [s for s in ALL_STRAINS if s in per_strain]
        n_strains = len(present_strains)
        max_copies = max(len(c) for c in strain_clusters.values())
        n_single = sum(1 for v in per_strain.values() if len(v) == 1)
        if n_strains == 8 and max_copies == 1:
            label = 'CORE-1:1'
        elif n_strains == 8 and max_copies >= 2:
            label = 'CORE-VAR'
        elif max_copies >= 3:
            label = 'FAMILY'
        elif n_strains <= 2:
            label = 'LINEAGE-SPECIFIC'
        else:
            label = 'PARTIAL'
        row = {
            'orthogroup_id': f'OG{len(rows)+1:06d}',
            'label': label,
            'n_strains': n_strains,
            'max_copies': max_copies,
        }
        for s in ALL_STRAINS:
            if s in strain_clusters:
                # One entry per physical-position cluster, IDs joined with '|' within a cluster
                row[s] = ','.join('|'.join(c) for c in strain_clusters[s])
            else:
                row[s] = '-'
        rows.append(row)

    print(f'  multi-strain orthogroups: {len(rows)}')

    # Stats
    from collections import Counter
    labels = Counter(r['label'] for r in rows)
    for k, n in labels.most_common():
        print(f'    {k}: {n}')

    out = Path('work/03_consensus/ortholog_table.tsv')
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ['orthogroup_id', 'label', 'n_strains', 'max_copies'] + ALL_STRAINS
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
