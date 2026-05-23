#!/usr/bin/env python3
"""
Synvisio-style synteny plots directly from wfmash PAF (Path A1 output).

Same rendering style as synvisio_pangenome.py and synvisio_graph.py, but the
input is the raw wfmash PAF — the same aligner PGGB uses internally to build
the graph. No chain net-filtering, no GFA P-line parsing.

PAF column convention (wfmash output for our PvP01-as-target invocation):
  col 1 = query   = OTHER strain contig
  col 6 = target  = PvP01 chromosome
  col 5 = strand
  col 11 = block length (post-CIGAR)
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


def parse_paf(path, min_block=1000):
    """Yield blocks: {qName, qLen, qStart, qEnd, strand, tName, tLen, tStart, tEnd, blockLen}.
    Note: wfmash invocation had PvP01 as target. So tName ∈ PvP01 chroms, qName is OTHER strain contig.
    """
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 12:
                continue
            blockLen = int(f[10])
            if blockLen < min_block:
                continue
            yield {
                'qName': f[0], 'qLen': int(f[1]), 'qStart': int(f[2]), 'qEnd': int(f[3]),
                'strand': f[4],
                'tName': f[5], 'tLen': int(f[6]), 'tStart': int(f[7]), 'tEnd': int(f[8]),
                'blockLen': blockLen,
            }


def get_query_layout(blocks, top_n=50):
    """Assign each query contig to its dominant PvP01 chrom, sort by dom chrom."""
    contig_dom = defaultdict(lambda: defaultdict(int))
    contig_size = {}
    for b in blocks:
        if b['tName'] not in CHR_COLOR:
            continue
        contig_dom[b['qName']][b['tName']] += b['blockLen']
        contig_size[b['qName']] = b['qLen']
    order = []
    for q, doms in contig_dom.items():
        dom = max(doms, key=doms.get)
        order.append((q, contig_size[q], dom, doms[dom]))
    order.sort(key=lambda x: (PVP01_ORDER.index(x[2]), -x[3]))
    if top_n:
        order = order[:top_n]
    return order


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


def render(query_name, paf_path, out_path, top_n=50, min_block=1000):
    blocks = [b for b in parse_paf(paf_path, min_block=min_block)
              if b['tName'] in CHR_COLOR]
    if not blocks:
        print(f"  [{query_name}] no PAF blocks ≥{min_block} bp")
        return
    layout = get_query_layout(blocks, top_n=top_n)
    chosen = {q[0] for q in layout}
    blocks = [b for b in blocks if b['qName'] in chosen]
    print(f"  [{query_name}] PAF blocks={len(blocks)} q_contigs={len(layout)}")

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
            q_offsets[q] = (x, qsize * scale, qsize)  # (offset, drawn_width, real_size)
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
        x0, w, _ = q_offsets[q]
        col = CHR_COLOR[dom]
        ax.add_patch(Rectangle((x0, y_bot), w, ROW_H, facecolor=col,
                               edgecolor='black', linewidth=0.3, alpha=0.7))

    for b in blocks:
        if b['qName'] not in q_offsets:
            continue
        t_x0 = t_offsets[b['tName']] + b['tStart']
        t_x1 = t_offsets[b['tName']] + b['tEnd']
        q_off, q_w, q_real = q_offsets[b['qName']]
        q_scale = q_w / q_real
        if b['strand'] == '-':
            q_x0 = q_off + (q_real - b['qEnd']) * q_scale
            q_x1 = q_off + (q_real - b['qStart']) * q_scale
            twist = True
        else:
            q_x0 = q_off + b['qStart'] * q_scale
            q_x1 = q_off + b['qEnd'] * q_scale
            twist = False
        col = CHR_COLOR[b['tName']]
        path = bezier_ribbon(t_x0, t_x1, y_top, q_x0, q_x1, y_bot + ROW_H, twist=twist)
        # Alpha by block size
        blen = b['blockLen']
        alpha = 0.65 if blen > 100000 else (0.4 if blen > 10000 else 0.2)
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
    ax.set_title(f'PvP01 → {query_name} synteny (wfmash PAF; -s 5000 -p 90 -n 1; min block ≥{min_block} bp)',
                 fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paf-dir', default='projection/A1_wfmash')
    ap.add_argument('--out-dir', default='writeup/synteny_paf')
    ap.add_argument('--top-contigs', type=int, default=50)
    ap.add_argument('--min-block', type=int, default=1000)
    ap.add_argument('--use-filtered', action='store_true', default=False,
                    help='Use PvP01_vs_*.filtered.paf if present (our v3 awk filter)')
    args = ap.parse_args()

    QUERIES = {
        'Sal-I':  'GCA_000002415.2', 'PvT01':  'GCA_900093545.1',
        'PvC01':  'GCA_900093535.1', 'PvW1':   'GCA_914969965.1',
        'PAM':    'GCA_949152365.1', 'PvSY56': 'GCA_003402215.1',
        'MHC087': 'GCA_040114635.1',
    }
    QUERIES_ORDER = ['PvW1', 'PAM', 'Sal-I', 'PvT01', 'PvC01', 'PvSY56', 'MHC087']

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    for qname in QUERIES_ORDER:
        qacc = QUERIES[qname]
        if args.use_filtered:
            paf = Path(args.paf_dir) / f'PvP01_vs_{qacc}.filtered.paf'
        else:
            paf = Path(args.paf_dir) / f'PvP01_vs_{qacc}.paf'
        if not paf.exists():
            print(f"  [{qname}] missing PAF {paf}")
            continue
        out = outdir / f'paf_synteny_PvP01_to_{qname}.png'
        render(qname, paf, out, top_n=args.top_contigs, min_block=args.min_block)


if __name__ == '__main__':
    main()
