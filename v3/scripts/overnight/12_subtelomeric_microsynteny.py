#!/usr/bin/env python3
"""Subtelomeric microsynteny plot: for each PvP01 chromosome's 5' and 3' end
(first/last 300 kb), draw a multi-strain ribbon plot of gene order colored by
variant-antigen family.
"""
import csv
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
import os
os.chdir('/media/anton/data/sandbox/Pv4/v3')

OUT_DIR = Path('writeup/microsynteny')
OUT_DIR.mkdir(exist_ok=True, parents=True)

SUBTEL_BP = 300_000
STRAINS = ['PvP01', 'Sal-I', 'PvW1', 'PAM', 'PvSY56', 'PvT01', 'PvC01', 'MHC087']

# Per-strain BED files (or merged annotation GFF if BED missing)
strain_beds = {
    'PvP01':  'inputs/annotations/PvP01.bed',
    'PvW1':   'inputs/annotations/PvW1.bed',
    'PAM':    'inputs/annotations/PAM.bed',
    'PvSY56': 'inputs/annotations/PvSY56.bed',
    # Sal-I, PvT01, PvC01, MHC087 — use Liftoff-projection annotation (PvP01 IDs)
}
# Fallback: extract per-strain gene intervals from merged GFFs
strain_gffs = {
    'Sal-I':  'work/02d_merged/PvP01-as-ref/Sal-I.annotation.gff3',
    'PvT01':  'work/02d_merged/PvP01-as-ref/PvT01.annotation.gff3',
    'PvC01':  'work/02d_merged/PvP01-as-ref/PvC01.annotation.gff3',
    'MHC087': 'work/02d_merged/PvP01-as-ref/MHC087.annotation.gff3',
}

def parse_gff_genes(path):
    """Yield (chrom, start, end, gene_id, strand)."""
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip(): continue
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('protein_coding_gene', 'ncRNA_gene', 'pseudogene'): continue
            m = re.search(r'ID=([^;]+)', f[8])
            if not m: continue
            gid = m.group(1).split(':', 1)[-1]
            yield f[0], int(f[3]) - 1, int(f[4]), gid, f[6]

# Load family labels (PvP01-strain-only)
gene_family = {}
with open('work/05_families/family_table.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        gene_family[(row['strain'], row['gene_id'])] = row['family']

FAMILY_COLORS = {
    'PIR':         '#d62728',  # red
    'PHIST':       '#9467bd',  # purple
    'Pv-fam':      '#ff7f0e',  # orange
    'MSP':         '#1f77b4',  # blue
    'DBP':         '#17becf',  # cyan
    'EBA':         '#bcbd22',  # olive
    'RBP':         '#e377c2',  # pink
    'AMA':         '#8c564b',  # brown
    'RAP':         '#2ca02c',  # green
    'SERA':        '#9edae5',
    'TRAg':        '#f7b6d2',
    'STP1':        '#ffbb78',
    'RESA':        '#c5b0d5',
    'tRNA':        '#dddddd',
    'rRNA':        '#dddddd',
    'ncRNA':       '#dddddd',
    'other':       '#cccccc',
    'conserved':   '#aaaaaa',
    'hypothetical': '#bbbbbb',
    'unannotated': '#eeeeee',
}

def family_color(family):
    return FAMILY_COLORS.get(family, '#cccccc')

# Per-strain genes per chromosome
strain_genes = defaultdict(lambda: defaultdict(list))  # strain → chrom → [(start, end, gid, strand, family)]
print("Loading per-strain genes...")
for strain in STRAINS:
    if strain in strain_beds and Path(strain_beds[strain]).exists():
        with open(strain_beds[strain]) as fh:
            for ln in fh:
                f = ln.rstrip('\n').split('\t')
                if len(f) < 4: continue
                chrom, s, e, gid = f[0], int(f[1]), int(f[2]), f[3]
                fam = gene_family.get((strain, gid), 'unannotated')
                strain_genes[strain][chrom].append((s, e, gid, f[5] if len(f)>=6 else '+', fam))
    elif strain in strain_gffs and Path(strain_gffs[strain]).exists():
        for chrom, s, e, gid, st in parse_gff_genes(strain_gffs[strain]):
            fam = gene_family.get((strain, gid), 'unannotated')
            strain_genes[strain][chrom].append((s, e, gid, st, fam))
    print(f"  {strain}: {sum(len(v) for v in strain_genes[strain].values())} genes")

# Get PvP01 chromosomes (LT635xxx 14 main chrs + LT635626.1 API + LT635627.1 MIT)
PVP01_CHRS_FULL = {
    'LT635612.2': ('chr1',  1021664), 'LT635613.2': ('chr2',   956327),
    'LT635614.2': ('chr3',   896704), 'LT635615.1': ('chr4',  1012024),
    'LT635616.2': ('chr5',  1524814), 'LT635617.2': ('chr6',  1042791),
    'LT635618.2': ('chr7',  1652210), 'LT635619.1': ('chr8',  1761288),
    'LT635620.2': ('chr9',  2237066), 'LT635621.2': ('chr10', 1548844),
    'LT635622.1': ('chr11', 2131221), 'LT635623.1': ('chr12', 3182763),
    'LT635624.2': ('chr13', 2093556), 'LT635625.2': ('chr14', 3153402),
}

# Cross-strain orthology from Phase E ortholog_table.tsv
ortho = {}  # (strain, gene_id) → orthogroup_id
print("Loading orthogroups...")
with open('work/03_consensus/ortholog_table.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        og = row['orthogroup_id']
        for strain in STRAINS:
            cell = row.get(strain, '-')
            if cell == '-': continue
            for gid in re.split(r'[,|]', cell):
                if gid:
                    ortho[(strain, gid)] = og
print(f"  loaded {len(ortho)} (strain, gene) → orthogroup mappings")

# For each PvP01 chromosome end, draw a microsynteny plot
def draw_subtel(pv_chrom_acc, end):
    """end: 'L' (5'/start) or 'R' (3'/end). Plot last 300 kb of all 8 strains'
    contigs that have orthologs in this PvP01 region."""
    chr_name, chr_len = PVP01_CHRS_FULL[pv_chrom_acc]
    if end == 'L':
        pv_window = (0, SUBTEL_BP)
        title_end = "5' subtelomere"
    else:
        pv_window = (chr_len - SUBTEL_BP, chr_len)
        title_end = "3' subtelomere"

    # PvP01 genes in window
    pv_genes_in_window = [g for g in strain_genes['PvP01'][pv_chrom_acc]
                          if g[0] >= pv_window[0] and g[1] <= pv_window[1]]
    if not pv_genes_in_window:
        return None

    # For each non-PvP01 strain, collect ortholog genes (anywhere)
    other_strain_genes = {}
    pv_ogs = set()
    for g in pv_genes_in_window:
        og = ortho.get(('PvP01', g[2]))
        if og: pv_ogs.add(og)

    for strain in STRAINS[1:]:
        # Find genes in this strain whose orthogroup is in pv_ogs
        hit_genes = []
        for chrom, genes in strain_genes[strain].items():
            for s, e, gid, st, fam in genes:
                og = ortho.get((strain, gid))
                if og and og in pv_ogs:
                    hit_genes.append((chrom, s, e, gid, st, fam, og))
        # Sort by chrom, position
        hit_genes.sort()
        other_strain_genes[strain] = hit_genes

    # Plotting
    n_strains = len(STRAINS)
    fig, ax = plt.subplots(figsize=(18, 1.2 * n_strains + 1))
    Y = {s: n_strains - 1 - i for i, s in enumerate(STRAINS)}
    BAR_H = 0.5
    BAR_Y_OFFSET = -BAR_H / 2

    # Draw PvP01 genes on its track
    for s, e, gid, st, fam in pv_genes_in_window:
        # Convert to local coord (subtract pv_window[0])
        x0 = s - pv_window[0]
        x1 = e - pv_window[0]
        color = family_color(fam)
        ax.add_patch(Rectangle((x0, Y['PvP01'] + BAR_Y_OFFSET), x1 - x0, BAR_H,
                               facecolor=color, edgecolor='black', linewidth=0.3))

    # Draw non-PvP01 strain genes. Layout: pack all hit genes into a single
    # horizontal track in their original gene-order; just place them sequentially.
    for strain in STRAINS[1:]:
        hits = other_strain_genes[strain]
        # Pack: use the actual genomic position within the strain's hit-genome region;
        # normalize so the first and last hit map to [0, SUBTEL_BP*1.2]
        if not hits: continue
        # Group by contig
        hits_by_chrom = defaultdict(list)
        for h in hits:
            hits_by_chrom[h[0]].append(h)
        # If hits span multiple contigs, just render them sequentially on the track
        x_cursor = 0
        for chrom, chrom_hits in sorted(hits_by_chrom.items(), key=lambda x: -len(x[1])):
            # Sort by position
            chrom_hits.sort(key=lambda h: h[1])
            for c, s, e, gid, st, fam, og in chrom_hits:
                # Use a fixed gene-width of ~3kb for visibility
                w = max(1500, e - s)
                color = family_color(fam)
                ax.add_patch(Rectangle((x_cursor, Y[strain] + BAR_Y_OFFSET),
                                       w, BAR_H,
                                       facecolor=color, edgecolor='black', linewidth=0.3))
                # Connector: ribbon from PvP01 to here
                # Find PvP01 gene with same OG
                for ps, pe, pgid, pst, pfam in pv_genes_in_window:
                    if ortho.get(('PvP01', pgid)) == og:
                        px = (ps + pe)/2 - pv_window[0]
                        ax.plot([px, x_cursor + w/2], [Y['PvP01'] - BAR_H/2, Y[strain] + BAR_H/2],
                                color=color, linewidth=0.3, alpha=0.5, zorder=0)
                        break
                x_cursor += w + 500
            x_cursor += 5000  # gap between contigs

    # Y-axis labels
    for s in STRAINS:
        ax.text(-5000, Y[s], s, ha='right', va='center', fontsize=10, fontweight='bold')

    # X-axis
    ax.set_xlim(-50000, max(SUBTEL_BP, max((len(other_strain_genes[s]) * 5500 for s in STRAINS[1:] if other_strain_genes[s]), default=SUBTEL_BP)) + 10000)
    ax.set_ylim(-0.7, n_strains - 0.3)
    ax.set_yticks([])
    ax.set_xlabel(f"position (bp) — PvP01 track: real coords from {pv_window[0]} to {pv_window[1]}; "
                  "other strains: packed by orthology-hit order", fontsize=9)
    for sp in ('top','left','right'):
        ax.spines[sp].set_visible(False)
    ax.set_title(f"PvP01 {chr_name} {title_end} ({pv_window[0]:,}-{pv_window[1]:,} bp) — "
                 f"{len(pv_genes_in_window)} PvP01 genes, "
                 f"{len(pv_ogs)} orthogroups; cross-strain orthologs colored by family",
                 fontsize=10)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = []
    seen_fams = set()
    for s, e, gid, st, fam in pv_genes_in_window:
        if fam not in seen_fams and fam != 'unannotated':
            seen_fams.add(fam)
    # Add families found in non-PvP01 strains too
    for strain in STRAINS[1:]:
        for c, s, e, gid, st, fam, og in other_strain_genes[strain]:
            if fam not in seen_fams and fam not in ('unannotated', 'other', 'conserved', 'hypothetical'):
                seen_fams.add(fam)
    for fam in sorted(seen_fams):
        legend_elements.append(Patch(facecolor=family_color(fam), edgecolor='black', label=fam))
    if legend_elements:
        ax.legend(handles=legend_elements, loc='lower right', ncol=4, fontsize=8, frameon=False)

    plt.tight_layout()
    out_path = OUT_DIR / f'{chr_name}_{end}.png'
    plt.savefig(out_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    return out_path


# Build all 14 chrs × 2 ends = 28 plots
print("\nRendering subtelomeric microsynteny plots...")
n_done = 0
for pv_acc in sorted(PVP01_CHRS_FULL.keys()):
    for end in ('L', 'R'):
        out = draw_subtel(pv_acc, end)
        if out:
            print(f"  ✓ {out}")
            n_done += 1
print(f"\nWrote {n_done} subtelomeric microsynteny plots to {OUT_DIR}")
