#!/usr/bin/env python3
"""
Synvisio-style synteny plots from PGGB-graph-derived blocks (output of
scripts/graph_synteny_blocks.py). Same rendering style as
scripts/synvisio_pangenome.py, but the input is graph blocks (BEDPE-like)
not pairwise chains.
"""
import argparse
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MPath

PVP01_CHRS = {
    'LT635612.2': ('chr1',  1021664), 'LT635613.2': ('chr2',   956327),
    'LT635614.2': ('chr3',   896704), 'LT635615.1': ('chr4',  1012024),
    'LT635616.2': ('chr5',  1524814), 'LT635617.2': ('chr6',  1042791),
    'LT635618.2': ('chr7',  1652210), 'LT635619.1': ('chr8',  1761288),
    'LT635620.2': ('chr9',  2237066), 'LT635621.2': ('chr10', 1548844),
    'LT635622.1': ('chr11', 2131221), 'LT635623.1': ('chr12', 3182763),
    'LT635624.2': ('chr13', 2093556), 'LT635625.2': ('chr14', 3153402),
    'LT635626.1': ('API',     29582), 'LT635627.1': ('MIT',      5989),
}
PVP01_ORDER = ['LT635612.2','LT635613.2','LT635614.2','LT635615.1','LT635616.2',
               'LT635617.2','LT635618.2','LT635619.1','LT635620.2','LT635621.2',
               'LT635622.1','LT635623.1','LT635624.2','LT635625.2','LT635626.1','LT635627.1']
CHR_COLORS = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f',
              '#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#c49c94']
CHR_COLOR = {c: CHR_COLORS[i] for i, c in enumerate(PVP01_ORDER)}

QUERIES = {
    'Sal-I':  'GCA_000002415.2', 'PvT01':  'GCA_900093545.1',
    'PvC01':  'GCA_900093535.1', 'PvW1':   'GCA_914969965.1',
    'PAM':    'GCA_949152365.1', 'PvSY56': 'GCA_003402215.1',
    'MHC087': 'GCA_040114635.1',
}


def parse_blocks(path):
    """Yield block dicts from BEDPE TSV."""
    with open(path) as fh:
        next(fh)  # header
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            yield {
                'tName': f[0], 'tStart': int(f[1]), 'tEnd': int(f[2]),
                'qName': f[3], 'qStart': int(f[4]), 'qEnd': int(f[5]),
                'strand': f[6], 'score': int(f[7]),
            }


def extract_short_id(panSN_name):
    """GCA_900093555.2#1#LT635625.2#0 -> LT635625.2 (3rd component)."""
    parts = panSN_name.split('#')
    return parts[2] if len(parts) >= 3 else panSN_name


def bezier_ribbon(x_top_l, x_top_r, y_top, x_bot_l, x_bot_r, y_bot, twist=False):
    if twist:
        x_bot_l, x_bot_r = x_bot_r, x_bot_l
    mid = (y_top + y_bot) / 2.0
    verts = [
        (x_top_l, y_top), (x_top_l, mid), (x_bot_l, mid), (x_bot_l, y_bot),
        (x_bot_r, y_bot), (x_bot_r, mid), (x_top_r, mid), (x_top_r, y_top),
        (x_top_l, y_top),
    ]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
             MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
    return MPath(verts, codes)


def get_query_layout(blocks_full, top_n=50):
    """Per query contig: dominant target chrom + total bp + intrinsic size."""
    contig_dom = defaultdict(lambda: defaultdict(int))  # qShortId → tShortId → bp
    contig_max = defaultdict(int)  # qShortId → max(qEnd) (approximate contig size)
    for b in blocks_full:
        t_acc = extract_short_id(b['tName'])
        q_acc = extract_short_id(b['qName'])
        if t_acc not in CHR_COLOR:
            continue
        contig_dom[q_acc][t_acc] += b['tEnd'] - b['tStart']
        contig_max[q_acc] = max(contig_max[q_acc], b['qEnd'])
    order = []
    for q, doms in contig_dom.items():
        dom = max(doms, key=doms.get)
        order.append((q, contig_max[q], dom, doms[dom]))
    order.sort(key=lambda x: (PVP01_ORDER.index(x[2]), -x[3]))
    if top_n:
        order = order[:top_n]
    return order


def render(query_name, blocks_path, out_path, top_n=50):
    blocks = [b for b in parse_blocks(blocks_path)
              if extract_short_id(b['tName']) in CHR_COLOR]
    if not blocks:
        print(f"  [{query_name}] no blocks")
        return
    layout = get_query_layout(blocks, top_n=top_n)
    chosen_q = {q[0] for q in layout}
    blocks = [b for b in blocks if extract_short_id(b['qName']) in chosen_q]
    print(f"  [{query_name}] blocks={len(blocks)} q_contigs={len(layout)}")

    gap = 100_000
    t_offsets = {}
    x = 0
    for acc in PVP01_ORDER:
        t_offsets[acc] = x
        x += PVP01_CHRS[acc][1] + gap
    t_total = x - gap

    dom_groups = defaultdict(list)
    for q, qsize, dom, _wt in layout:
        dom_groups[dom].append((q, qsize))

    q_offsets = {}
    for acc, ctgs in dom_groups.items():
        total = sum(s for _, s in ctgs)
        target = PVP01_CHRS[acc][1]
        scale = target / total if total > target else 1.0
        x = t_offsets[acc]
        small_gap = 5000
        for q, qsize in ctgs:
            q_offsets[q] = (x, qsize * scale)
            x += qsize * scale + small_gap

    fig, ax = plt.subplots(figsize=(20, 6))
    ROW_H = 100000
    y_top = 700000
    y_bot = 0

    for acc in PVP01_ORDER:
        size = PVP01_CHRS[acc][1]
        x0 = t_offsets[acc]
        col = CHR_COLOR[acc]
        ax.add_patch(Rectangle((x0, y_top), size, ROW_H, facecolor=col,
                               edgecolor='black', linewidth=0.5, alpha=0.9))
        ax.text(x0 + size/2, y_top + ROW_H + 30000, PVP01_CHRS[acc][0],
                ha='center', va='bottom', fontsize=8, color=col, fontweight='bold')

    for q, qsize, dom, _wt in layout:
        x0, w = q_offsets[q]
        col = CHR_COLOR[dom]
        ax.add_patch(Rectangle((x0, y_bot), w, ROW_H, facecolor=col,
                               edgecolor='black', linewidth=0.3, alpha=0.7))

    for b in blocks:
        t_acc = extract_short_id(b['tName'])
        q_acc = extract_short_id(b['qName'])
        if q_acc not in q_offsets:
            continue
        t_x0 = t_offsets[t_acc] + b['tStart']
        t_x1 = t_offsets[t_acc] + b['tEnd']
        q_off, q_w = q_offsets[q_acc]
        # qStart/qEnd are in cumulative path bp (per node sum). Use contig-relative for proportional rendering:
        # since contig_max approximates contig size, use those as relative coords
        q_scale = q_w / max(1, q_offsets[q_acc][1] if q_offsets[q_acc][1] > 0 else 1)
        # contig_max as approximation:
        q_max = max(b['qEnd'], 1)
        # Project linearly
        q_x0 = q_off + b['qStart'] / max(b['qEnd'], q_max) * q_w
        q_x1 = q_off + b['qEnd']   / max(b['qEnd'], q_max) * q_w
        twist = (b['strand'] == '-')
        col = CHR_COLOR[t_acc]
        path = bezier_ribbon(t_x0, t_x1, y_top, q_x0, q_x1, y_bot + ROW_H, twist=twist)
        score = b['score']
        alpha = 0.65 if score > 100000 else (0.4 if score > 10000 else 0.2)
        ax.add_patch(PathPatch(path, facecolor=col, edgecolor='none', alpha=alpha))

    ax.text(-2000000, y_top + ROW_H/2, 'PvP01', ha='right', va='center',
            fontsize=12, fontweight='bold')
    ax.text(-2000000, y_bot + ROW_H/2, query_name, ha='right', va='center',
            fontsize=12, fontweight='bold')
    ax.set_xlim(-3500000, t_total + 500000)
    ax.set_ylim(-200000, y_top + ROW_H + 200000)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f'PvP01 → {query_name} synteny (PGGB graph blocks; node-shared ≥1kb)',
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--blocks-dir', default='writeup/synteny_graph')
    ap.add_argument('--out-dir', default='writeup/synteny_graph')
    ap.add_argument('--top-contigs', type=int, default=50)
    args = ap.parse_args()
    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for qname in ['PvW1','PAM','Sal-I','PvT01','PvC01','PvSY56','MHC087']:
        blocks = Path(args.blocks_dir) / f'PvP01_to_{qname}.blocks.tsv'
        if not blocks.exists():
            print(f"  [{qname}] missing {blocks}")
            continue
        out = outdir / f'graph_synteny_PvP01_to_{qname}.png'
        render(qname, blocks, out, top_n=args.top_contigs)


if __name__ == '__main__':
    main()
