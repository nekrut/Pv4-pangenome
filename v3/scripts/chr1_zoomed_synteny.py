#!/usr/bin/env python3
"""
Render chr1-zoomed synteny panels for the 5 sources, then combine with the
odgi viz panels into one chr1-focused image.

For each source:
  - x-axis = PvP01 chr 1 only (LT635612.2, 1,021,664 bp)
  - bottom track = PAM contigs that have any block touching chr1
  - PAM contigs sorted by their chr1 alignment weight (largest first)
  - blocks targeting OTHER chromosomes on the same PAM contig are still drawn
    as truncated ribbons reaching the edge — visualizes that the PAM contig
    has cargo from other PvP01 chromosomes too
"""
import sys
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path as MPath

sys.path.insert(0, 'scripts')
from synvisio_canonical_layout import (
    PVP01_CHRS, CHR_COLOR, SOURCES, strip_pansn,
    parse_chain, parse_paf, parse_graph, normalize_qname)
from synteny_with_labels import chain_with_density

CHR1 = 'LT635612.2'
CHR1_LEN = PVP01_CHRS[CHR1][1]


def get_blocks(kind, path):
    if kind == 'chain':
        return chain_with_density(path, min_density=0.5)
    elif kind == 'paf':
        return list(parse_paf(path))
    elif kind == 'graph':
        return list(parse_graph(path))


def bezier_ribbon(x0L, x0R, yT, x1L, x1R, yB, twist=False):
    if twist: x1L, x1R = x1R, x1L
    mid = (yT + yB) / 2
    verts = [(x0L,yT),(x0L,mid),(x1L,mid),(x1L,yB),(x1R,yB),(x1R,mid),(x0R,mid),(x0R,yT),(x0L,yT)]
    codes = [MPath.MOVETO,MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.LINETO,
             MPath.CURVE4,MPath.CURVE4,MPath.CURVE4,MPath.CLOSEPOLY]
    return MPath(verts, codes)


def render_chr1(source_key, kind, path, out_png, title):
    all_blocks = get_blocks(kind, path)
    # PAM contigs with any chr1 alignment, sorted by chr1 alignment weight
    chr1_weight = defaultdict(int)
    contig_size = {}
    for b in all_blocks:
        if b['tName'] == CHR1:
            qn = normalize_qname(b['qName'], kind)
            chr1_weight[qn] += b.get('aligned_bp', b['tEnd'] - b['tStart'])
            contig_size[qn] = max(contig_size.get(qn, 0), b.get('qSize', b['qEnd']))
        # Also record any PAM contigs we see (for finding their size)
        qn = normalize_qname(b['qName'], kind)
        if qn in contig_size:
            contig_size[qn] = max(contig_size.get(qn, 0), b.get('qSize', b['qEnd']))

    sorted_contigs = sorted(chr1_weight.keys(), key=lambda q: -chr1_weight[q])
    if not sorted_contigs:
        print(f'  [{source_key}] no chr1 blocks')
        return

    # Layout: chr1 is wide; below it, just-touching PAM contigs side-by-side
    fig, ax = plt.subplots(figsize=(20, 5))
    ROW_H = 100000; yT = 700000; yB = 0

    # Top track: just PvP01 chr1
    ax.add_patch(Rectangle((0, yT), CHR1_LEN, ROW_H, facecolor=CHR_COLOR[CHR1],
                           edgecolor='black', linewidth=0.5, alpha=0.9))
    ax.text(CHR1_LEN/2, yT + ROW_H + 30000, 'chr1 (LT635612.2)',
            ha='center', va='bottom', fontsize=11,
            color=CHR_COLOR[CHR1], fontweight='bold')

    # Bottom track: PAM contigs that touch chr1
    # Pack into chr1 width; each contig sized proportional to its chr1-aligned bp
    total_chr1_bp = sum(chr1_weight.values())
    pad_between = max(1000, CHR1_LEN * 0.005)
    available_w = CHR1_LEN - pad_between * (len(sorted_contigs) - 1)
    q_offsets = {}
    x = 0
    for q in sorted_contigs:
        w = max(5000, available_w * (chr1_weight[q] / total_chr1_bp))
        q_offsets[q] = (x, w, contig_size[q])
        x += w + pad_between

    # Draw PAM contig rectangles with labels
    for q in sorted_contigs:
        x0, w, real = q_offsets[q]
        # color by dominant target of this contig (could be chr1 or another)
        # compute contig's dominant target chr
        dom = CHR1  # all here touch chr1; let's also note if their dominant elsewhere
        # check dominant chr from all_blocks for this contig
        c_t_bp = defaultdict(int)
        for b in all_blocks:
            if normalize_qname(b['qName'], kind) == q and b['tName'] in CHR_COLOR:
                c_t_bp[b['tName']] += b.get('aligned_bp', b['tEnd'] - b['tStart'])
        if c_t_bp:
            dom = max(c_t_bp, key=c_t_bp.get)
        col = CHR_COLOR[dom]
        ax.add_patch(Rectangle((x0, yB), w, ROW_H, facecolor=col, edgecolor='black',
                               linewidth=0.4, alpha=0.7))
        short = q if len(q) <= 14 else q[:7] + '…' + q[-6:]
        label = f'{short} ({real/1000:.0f}kb)\ndom={PVP01_CHRS[dom][0]}'
        ax.text(x0 + w/2, yB - 50000, label,
                ha='center', va='top', fontsize=7, color='#222',
                rotation=-30, rotation_mode='anchor')

    # Draw ribbons for chr1 blocks; for blocks targeting OTHER chrs on the same
    # PAM contig, draw "dangling" ribbons that exit the top track at the chr1
    # boundary
    for b in all_blocks:
        qn = normalize_qname(b['qName'], kind)
        if qn not in q_offsets:
            continue
        q_off, q_w, q_real = q_offsets[qn]
        q_scale = q_w / max(q_real, 1)
        if b['strand'] == '-':
            q_x0 = q_off + (q_real - b['qEnd']) * q_scale
            q_x1 = q_off + (q_real - b['qStart']) * q_scale
            twist = True
        else:
            q_x0 = q_off + b['qStart'] * q_scale
            q_x1 = q_off + b['qEnd'] * q_scale
            twist = False
        # If target is chr1, ribbon goes from chr1 position to query position
        if b['tName'] == CHR1:
            t_x0 = b['tStart']
            t_x1 = b['tEnd']
            col = CHR_COLOR[CHR1]
            score = b.get('score', b['tEnd'] - b['tStart'])
            alpha = 0.65 if score > 100000 else (0.4 if score > 10000 else 0.25)
            ax.add_patch(PathPatch(bezier_ribbon(t_x0, t_x1, yT, q_x0, q_x1, yB + ROW_H, twist),
                                   facecolor=col, edgecolor='none', alpha=alpha))
        elif b['tName'] in CHR_COLOR:
            # Block targets OTHER chr; draw dangling ribbon to the top edge
            # showing this contig also has cargo from elsewhere
            other_col = CHR_COLOR[b['tName']]
            # Truncated: top endpoints at the chr1 boundaries (either side)
            t_dummy_x = -200000 if q_off < CHR1_LEN/2 else CHR1_LEN + 200000
            # Use a thin "exiting" indicator
            ax.add_patch(Rectangle((q_x0, yB + ROW_H + 5000),
                                   q_x1 - q_x0, 25000,
                                   facecolor=other_col, edgecolor='none', alpha=0.6))

    ax.text(-150000, yT + ROW_H/2, 'PvP01', ha='right', va='center',
            fontsize=12, fontweight='bold')
    ax.text(-150000, yB + ROW_H/2, 'PAM',   ha='right', va='center',
            fontsize=12, fontweight='bold')
    ax.set_xlim(-300000, CHR1_LEN + 100000)
    ax.set_ylim(-650000, yT + ROW_H + 200000)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title, fontsize=12)
    plt.tight_layout(); plt.savefig(out_png, dpi=140, bbox_inches='tight'); plt.close(fig)
    print(f'  wrote {out_png}')


def main():
    outdir = Path('writeup/synteny_3way/chr1_zoomed_panels')
    outdir.mkdir(parents=True, exist_ok=True)
    titles = {
        '1_kegalign_default': 'KegAlign chains, default — density ≥50%',
        '2_kegalign_tuned':   'KegAlign chains, TUNED — density ≥50%',
        '3_wfmash_n1':        'wfmash PAF (A1: -n 1)',
        '4_wfmash_pggblike':  'wfmash PAF (PGGB: -n 8)',
        '5_graph_blocks':     'PGGB graph blocks (8-way)',
        '6_graph_blocks_2way': 'PGGB graph blocks (2-way, PvP01+PAM)',
    }
    SOURCES_EXT = dict(SOURCES)
    SOURCES_EXT['6_graph_blocks_2way'] = ('graph', 'writeup/synteny_2way/PvP01_to_PAM_2way.blocks.tsv')
    for k, (kind, path) in SOURCES_EXT.items():
        render_chr1(k, kind, path, str(outdir / f'{k}.png'),
                    f'PvP01 chr1 → PAM: {titles[k]}')


if __name__ == '__main__':
    main()
