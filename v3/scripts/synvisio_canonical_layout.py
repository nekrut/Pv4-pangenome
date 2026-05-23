#!/usr/bin/env python3
"""
Render PAM synteny from any of the 5 source formats using a CANONICAL PAM
contig layout (computed from the union of all sources). This way the same
PAM contig sits at the same x-position across all panels.
"""
import argparse, sys, json, re
from pathlib import Path
from collections import defaultdict
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


def strip_pansn(name):
    # GCA_xxx#1#CONTIG#0 → CONTIG  (or pass through if not PanSN)
    parts = name.split('#')
    return parts[2] if len(parts) >= 3 else name


def parse_chain(path):
    """Yield blocks: {tName, tStart, tEnd, qName, qStart, qEnd, strand, score}."""
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith('chain '): continue
            f = ln.split()
            if len(f) < 13: continue
            yield dict(tName=f[2], tStart=int(f[5]), tEnd=int(f[6]),
                       qName=f[7], qSize=int(f[8]), qStart=int(f[10]), qEnd=int(f[11]),
                       strand=f[9], score=int(f[1]))


def parse_paf(path, min_block=1000):
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 12: continue
            bl = int(f[10])
            if bl < min_block: continue
            yield dict(tName=f[5], tStart=int(f[7]), tEnd=int(f[8]),
                       qName=f[0], qSize=int(f[1]), qStart=int(f[2]), qEnd=int(f[3]),
                       strand=f[4], score=bl)


def parse_graph(path, min_block=1000):
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            sc = int(f[7])
            if sc < min_block: continue
            yield dict(tName=strip_pansn(f[0]), tStart=int(f[1]), tEnd=int(f[2]),
                       qName=strip_pansn(f[3]), qSize=0,  # unknown — fill from contig_size if present
                       qStart=int(f[4]), qEnd=int(f[5]),
                       strand=f[6], score=sc)


SOURCES = {
    '1_kegalign_default': ('chain', 'work/01_chains/GCA_900093555.2.GCA_949152365.1.cleaned.chain'),
    '2_kegalign_tuned':   ('chain', 'projection/A1_wfmash/synteny_rerun/PvP01_to_PAM.tuned.cleaned.chain'),
    '3_wfmash_n1':        ('paf',   'projection/A1_wfmash/PvP01_vs_GCA_949152365.1.paf'),
    '4_wfmash_pggblike':  ('paf',   'projection/A1_wfmash/synteny_rerun/PvP01_vs_PAM.pggblike.paf'),
    '5_graph_blocks':     ('graph', 'writeup/synteny_graph/PvP01_to_PAM.blocks.tsv'),
    '6_graph_blocks_2way': ('graph', 'writeup/synteny_2way/PvP01_to_PAM_2way.blocks.tsv'),
}


def normalize_qname(qname, source_kind):
    """Make qName comparable across sources. PAF / chain use raw FASTA names,
    graph blocks use the PanSN's contig segment."""
    return strip_pansn(qname)


def build_canonical_layout(min_block=1000):
    """Aggregate dominance from ALL sources; one layout used by all renderers."""
    contig_dom = defaultdict(lambda: defaultdict(int))  # qname → {tName → bp}
    contig_size = {}
    for src_name, (kind, path) in SOURCES.items():
        if kind == 'chain':
            blocks = parse_chain(path)
        elif kind == 'paf':
            blocks = parse_paf(path, min_block=min_block)
        elif kind == 'graph':
            blocks = parse_graph(path, min_block=min_block)
        for b in blocks:
            t = b['tName']
            if t not in CHR_COLOR: continue
            q = normalize_qname(b['qName'], kind)
            contig_dom[q][t] += b['tEnd'] - b['tStart']
            if b.get('qSize', 0) > 0:
                contig_size[q] = max(contig_size.get(q, 0), b['qSize'])
            else:
                contig_size[q] = max(contig_size.get(q, 0), b['qEnd'])
    layout = []
    for q, doms in contig_dom.items():
        dom = max(doms, key=doms.get)
        layout.append((q, contig_size.get(q, 0), dom, doms[dom]))
    layout.sort(key=lambda x: (PVP01_ORDER.index(x[2]), -x[3]))
    return layout


def make_q_offsets(layout, t_offsets):
    """Pack contigs under their dominant target chrom."""
    dom_groups = defaultdict(list)
    for q, qsize, dom, _wt in layout:
        dom_groups[dom].append((q, qsize))
    q_offsets = {}
    for acc, ctgs in dom_groups.items():
        total = sum(s for _, s in ctgs)
        target = PVP01_CHRS[acc][1]
        scale = target / total if total > target else 1.0
        x = t_offsets[acc]; small_gap = 5000
        for q, qsize in ctgs:
            q_offsets[q] = (x, qsize * scale, qsize)
            x += qsize * scale + small_gap
    return q_offsets


def bezier_ribbon(x0L, x0R, yT, x1L, x1R, yB, twist=False):
    if twist:
        x1L, x1R = x1R, x1L
    mid = (yT + yB) / 2
    verts = [(x0L,yT),(x0L,mid),(x1L,mid),(x1L,yB),(x1R,yB),(x1R,mid),(x0R,mid),(x0R,yT),(x0L,yT)]
    codes = [MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,
             MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
    return MPath(verts, codes)


def render_one(source_key, kind, path, t_offsets, q_offsets, out_png, title):
    if kind == 'chain':
        blocks = list(parse_chain(path))
    elif kind == 'paf':
        blocks = list(parse_paf(path))
    elif kind == 'graph':
        blocks = list(parse_graph(path))
    # Filter to canonical layout's contigs
    valid_q = set(q_offsets)
    blocks = [b for b in blocks
              if b['tName'] in CHR_COLOR and normalize_qname(b['qName'], kind) in valid_q]

    fig, ax = plt.subplots(figsize=(20, 6))
    ROW_H = 100000; yT = 700000; yB = 0
    t_total = max(t_offsets.values()) + max(PVP01_CHRS[a][1] for a in PVP01_ORDER)
    # Top track
    for acc in PVP01_ORDER:
        size = PVP01_CHRS[acc][1]; x0 = t_offsets[acc]; col = CHR_COLOR[acc]
        ax.add_patch(Rectangle((x0, yT), size, ROW_H, facecolor=col, edgecolor='black',
                               linewidth=0.5, alpha=0.9))
        ax.text(x0 + size/2, yT + ROW_H + 30000, PVP01_CHRS[acc][0],
                ha='center', va='bottom', fontsize=8, color=col, fontweight='bold')
    # Bottom contigs in canonical order
    for q, (x0, w, _real) in q_offsets.items():
        # color by canonical dominant chrom; compute from layout — pass via attr or recompute
        # quick recompute: find which chrom this contig was assigned to
        pass
    # Use the layout's dom assignment from x0 mapping (find chrom by t_offsets bucket)
    sorted_offsets = sorted(t_offsets.items(), key=lambda x: x[1])
    for q, (x0, w, real) in q_offsets.items():
        # find dom chrom by checking which target chrom interval x0 falls into
        dom = None
        for acc in PVP01_ORDER:
            if x0 >= t_offsets[acc] and x0 < t_offsets[acc] + PVP01_CHRS[acc][1] + 200000:
                dom = acc
        col = CHR_COLOR.get(dom, 'gray')
        ax.add_patch(Rectangle((x0, yB), w, ROW_H, facecolor=col, edgecolor='black',
                               linewidth=0.3, alpha=0.7))

    # Ribbons
    for b in blocks:
        q = normalize_qname(b['qName'], kind)
        if q not in q_offsets: continue
        t_x0 = t_offsets[b['tName']] + b['tStart']
        t_x1 = t_offsets[b['tName']] + b['tEnd']
        q_off, q_w, q_real = q_offsets[q]
        q_scale = q_w / max(q_real, 1)
        if b['strand'] == '-':
            q_x0 = q_off + (q_real - b['qEnd']) * q_scale
            q_x1 = q_off + (q_real - b['qStart']) * q_scale
            twist = True
        else:
            q_x0 = q_off + b['qStart'] * q_scale
            q_x1 = q_off + b['qEnd'] * q_scale
            twist = False
        col = CHR_COLOR[b['tName']]
        score = b['score']
        alpha = 0.65 if score > 100000 else (0.4 if score > 10000 else 0.2)
        ax.add_patch(PathPatch(bezier_ribbon(t_x0, t_x1, yT, q_x0, q_x1, yB + ROW_H, twist),
                               facecolor=col, edgecolor='none', alpha=alpha))

    ax.text(-2000000, yT + ROW_H/2, 'PvP01', ha='right', va='center', fontsize=12, fontweight='bold')
    ax.text(-2000000, yB + ROW_H/2, 'PAM',   ha='right', va='center', fontsize=12, fontweight='bold')
    ax.set_xlim(-3500000, t_total + 500000)
    ax.set_ylim(-200000, yT + ROW_H + 200000)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title, fontsize=13)
    plt.tight_layout(); plt.savefig(out_png, dpi=140, bbox_inches='tight'); plt.close(fig)


def main():
    # 1) Canonical layout from union of all sources
    layout = build_canonical_layout(min_block=1000)
    print(f'canonical layout: {len(layout)} PAM contigs')
    # t_offsets
    gap = 100_000
    t_offsets = {}; x = 0
    for acc in PVP01_ORDER:
        t_offsets[acc] = x
        x += PVP01_CHRS[acc][1] + gap
    q_offsets = make_q_offsets(layout, t_offsets)
    # 2) Render each source
    titles = {
        '1_kegalign_default': 'KegAlign chains, default (HoxD70, hspthresh 3000, seed 12of19)',
        '2_kegalign_tuned':   'KegAlign chains, TUNED (matrix +100/-100, hspthresh 4500, seed 14of22)',
        '3_wfmash_n1':        'wfmash PAF (A1 params: -n 1, best-per-segment)',
        '4_wfmash_pggblike':  'wfmash PAF (PGGB-build params: -n 8, multi-mapping)',
        '5_graph_blocks':     'PGGB graph blocks, 8-way (all strains)',
        '6_graph_blocks_2way': 'PGGB graph blocks, 2-way (PvP01 + PAM only)',
    }
    outdir = Path('writeup/synteny_canonical')
    outdir.mkdir(parents=True, exist_ok=True)
    for k, (kind, path) in SOURCES.items():
        render_one(k, kind, path, t_offsets, q_offsets,
                   str(outdir / f'PvP01_to_PAM__{k}.png'),
                   f'{titles[k]}: PvP01 → PAM')
        print(f'  wrote {outdir / (k + ".png")}')


if __name__ == '__main__':
    main()
