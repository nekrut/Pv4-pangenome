#!/usr/bin/env python3
"""
Synvisio-style pairwise synteny plot for the v3 pangenome.

For each of 7 query strains, render a stacked panel:
  TOP track: PvP01 reference chromosomes (16: chr1-14 + API + MIT)
  BOTTOM track: query strain's contigs, laid out by their dominant PvP01 chrom
  RIBBONS: bezier curves connecting chain blocks; twisted for inversions

Chain source: work/01_chains/${T_ACC}.${Q_ACC}.cleaned.chain (Phase B output).
Filters chains by minimum score to keep figure legible.
"""
import argparse
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch, Patch
from matplotlib.path import Path as MPath

# PvP01 reference chromosomes — GenBank accessions
PVP01_CHRS = {
    'LT635612.2': ('chr1',  1021664),
    'LT635613.2': ('chr2',   956327),
    'LT635614.2': ('chr3',   896704),
    'LT635615.1': ('chr4',  1012024),
    'LT635616.2': ('chr5',  1524814),
    'LT635617.2': ('chr6',  1042791),
    'LT635618.2': ('chr7',  1652210),
    'LT635619.1': ('chr8',  1761288),
    'LT635620.2': ('chr9',  2237066),
    'LT635621.2': ('chr10', 1548844),
    'LT635622.1': ('chr11', 2131221),
    'LT635623.1': ('chr12', 3182763),
    'LT635624.2': ('chr13', 2093556),
    'LT635625.2': ('chr14', 3153402),
    'LT635626.1': ('API',     29582),
    'LT635627.1': ('MIT',      5989),
}
PVP01_ORDER = ['LT635612.2','LT635613.2','LT635614.2','LT635615.1','LT635616.2',
               'LT635617.2','LT635618.2','LT635619.1','LT635620.2','LT635621.2',
               'LT635622.1','LT635623.1','LT635624.2','LT635625.2','LT635626.1','LT635627.1']

# Distinct colors for 16 chromosomes (perceptually-separated qualitative palette)
CHR_COLORS = [
    '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f',
    '#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c49c94',
]
CHR_COLOR = {c: CHR_COLORS[i] for i, c in enumerate(PVP01_ORDER)}

QUERIES = {
    'Sal-I':  'GCA_000002415.2',
    'PvT01':  'GCA_900093545.1',
    'PvC01':  'GCA_900093535.1',
    'PvW1':   'GCA_914969965.1',
    'PAM':    'GCA_949152365.1',
    'PvSY56': 'GCA_003402215.1',
    'MHC087': 'GCA_040114635.1',
}
QUERIES_ORDER = ['PvW1', 'PAM', 'Sal-I', 'PvT01', 'PvC01', 'PvSY56', 'MHC087']

PVP01_ACC = 'GCA_900093555.2'


def parse_chains(path, min_score=10000):
    """Yield (score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd)."""
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith('chain '):
                continue
            f = ln.split()
            if len(f) < 13:
                continue
            score = int(f[1])
            if score < min_score:
                continue
            yield {
                'score': score,
                'tName': f[2], 'tSize': int(f[3]), 'tStrand': f[4],
                'tStart': int(f[5]), 'tEnd': int(f[6]),
                'qName': f[7], 'qSize': int(f[8]), 'qStrand': f[9],
                'qStart': int(f[10]), 'qEnd': int(f[11]),
            }


def get_query_contig_layout(chains, top_n_contigs=None):
    """Assign each query contig to its dominant PvP01 chrom, return ordered layout.

    Returns:
        order: list of (qName, qSize, dominant_PvP01_chrom_acc, sum_match_bp)
    """
    contig_dominance = defaultdict(lambda: defaultdict(int))  # qName → {tName → bp}
    contig_size = {}
    for c in chains:
        if c['tName'] not in CHR_COLOR:
            continue
        match_bp = c['tEnd'] - c['tStart']
        contig_dominance[c['qName']][c['tName']] += match_bp
        contig_size[c['qName']] = c['qSize']

    order = []
    for qname, sizes in contig_dominance.items():
        dom = max(sizes, key=sizes.get)
        order.append((qname, contig_size[qname], dom, sizes[dom]))

    # Sort: first by PvP01 chrom order, then by dominance (largest first)
    order.sort(key=lambda x: (PVP01_ORDER.index(x[2]), -x[3]))
    if top_n_contigs:
        order = order[:top_n_contigs]
    return order


def bezier_ribbon(x_top_l, x_top_r, y_top, x_bot_l, x_bot_r, y_bot, twist=False):
    if twist:
        x_bot_l, x_bot_r = x_bot_r, x_bot_l
    mid = (y_top + y_bot) / 2.0
    verts = [
        (x_top_l, y_top),
        (x_top_l, mid), (x_bot_l, mid), (x_bot_l, y_bot),
        (x_bot_r, y_bot),
        (x_bot_r, mid), (x_top_r, mid), (x_top_r, y_top),
        (x_top_l, y_top),
    ]
    codes = [MPath.MOVETO,
             MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
             MPath.LINETO,
             MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
             MPath.CLOSEPOLY]
    return MPath(verts, codes)


def render_pairwise(query_name, chain_path, out_path, min_score, top_n_contigs):
    chains = list(parse_chains(chain_path, min_score=min_score))
    if not chains:
        print(f"  [{query_name}] no chains above min_score={min_score}")
        return
    # Filter chains targeting PvP01 chroms only
    chains = [c for c in chains if c['tName'] in CHR_COLOR]
    # Layout
    layout = get_query_contig_layout(chains, top_n_contigs=top_n_contigs)
    chosen_qnames = {q[0] for q in layout}
    chains = [c for c in chains if c['qName'] in chosen_qnames]
    print(f"  [{query_name}] chains={len(chains)} contigs={len(layout)}")

    # Build PvP01 x-coordinate layout
    gap = 100_000
    t_offsets = {}
    x = 0
    for acc in PVP01_ORDER:
        t_offsets[acc] = x
        x += PVP01_CHRS[acc][1] + gap
    t_total = x - gap

    # Build query x-coordinate layout (proportional sizes with gaps within each
    # dominant PvP01 chrom)
    q_offsets = {}
    # Group query contigs by their dominant chrom and place them under that chrom
    dom_groups = defaultdict(list)
    for qname, qsize, dom_acc, _wt in layout:
        dom_groups[dom_acc].append((qname, qsize))
    for acc, ctgs in dom_groups.items():
        # Total contig length in this group
        total = sum(s for _, s in ctgs)
        target = PVP01_CHRS[acc][1]
        # Scale factor — fit contigs into target's width
        scale = target / total if total > target else 1.0
        # Place at t_offsets[acc]
        x = t_offsets[acc]
        small_gap = 5000
        for qname, qsize in ctgs:
            q_offsets[qname] = (x, qsize * scale)
            x += qsize * scale + small_gap

    # Plot
    fig, ax = plt.subplots(figsize=(20, 6))
    ROW_H = 100000
    y_top = 700000
    y_bot = 0

    # Top track: PvP01 chroms
    for acc in PVP01_ORDER:
        size = PVP01_CHRS[acc][1]
        x0 = t_offsets[acc]
        col = CHR_COLOR[acc]
        ax.add_patch(Rectangle((x0, y_top), size, ROW_H, facecolor=col,
                               edgecolor='black', linewidth=0.5, alpha=0.9))
        # Label
        ax.text(x0 + size/2, y_top + ROW_H + 30000, PVP01_CHRS[acc][0],
                ha='center', va='bottom', fontsize=8, color=col, fontweight='bold')

    # Bottom track: query contigs
    for qname, qsize, dom_acc, _wt in layout:
        x0, w = q_offsets[qname]
        col = CHR_COLOR[dom_acc]
        ax.add_patch(Rectangle((x0, y_bot), w, ROW_H, facecolor=col,
                               edgecolor='black', linewidth=0.3, alpha=0.7))

    # Ribbons
    for c in chains:
        if c['qName'] not in q_offsets:
            continue
        t_x0 = t_offsets[c['tName']] + c['tStart']
        t_x1 = t_offsets[c['tName']] + c['tEnd']
        q_off, q_w = q_offsets[c['qName']]
        q_scale = q_w / c['qSize']
        if c['qStrand'] == '-':
            # qStart/qEnd are on +strand coords already in chain spec; - flips
            q_x0 = q_off + (c['qSize'] - c['qEnd']) * q_scale
            q_x1 = q_off + (c['qSize'] - c['qStart']) * q_scale
            twist = True
        else:
            q_x0 = q_off + c['qStart'] * q_scale
            q_x1 = q_off + c['qEnd'] * q_scale
            twist = False
        col = CHR_COLOR[c['tName']]
        path = bezier_ribbon(t_x0, t_x1, y_top,
                             q_x0, q_x1, y_bot + ROW_H,
                             twist=twist)
        # Alpha scaled by chain score (smaller chains more transparent)
        alpha = 0.6 if c['score'] > 1000000 else (0.4 if c['score'] > 100000 else 0.2)
        ax.add_patch(PathPatch(path, facecolor=col, edgecolor='none', alpha=alpha))

    # Labels
    ax.text(-2000000, y_top + ROW_H/2, 'PvP01', ha='right', va='center',
            fontsize=12, fontweight='bold')
    ax.text(-2000000, y_bot + ROW_H/2, query_name, ha='right', va='center',
            fontsize=12, fontweight='bold')

    ax.set_xlim(-3500000, t_total + 500000)
    ax.set_ylim(-200000, y_top + ROW_H + 200000)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f'PvP01 → {query_name} synteny (Phase B chains; min score {min_score})',
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chain-dir', default='work/01_chains')
    ap.add_argument('--out-dir', default='writeup/synteny')
    ap.add_argument('--min-score', type=int, default=10000)
    ap.add_argument('--top-contigs', type=int, default=50)
    args = ap.parse_args()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    for qname in QUERIES_ORDER:
        qacc = QUERIES[qname]
        chain = Path(args.chain_dir) / f'{PVP01_ACC}.{qacc}.cleaned.chain'
        if not chain.exists():
            print(f'  [{qname}] missing chain {chain}')
            continue
        out = outdir / f'synteny_PvP01_to_{qname}.png'
        render_pairwise(qname, chain, out, args.min_score, args.top_contigs)


if __name__ == '__main__':
    main()
