#!/usr/bin/env python3
"""
Dot plot of PvP01 chr 1 + chr 2 vs PAM, from KegAlign tuned chains.
Each chain block (size, dt, dq) contributes one matched segment as a line.
"""
import sys
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

CHR1 = ('LT635612.2', 'chr1', 1_021_664)
CHR2 = ('LT635613.2', 'chr2',   956_327)
TARGETS = [CHR1, CHR2]
CHAIN_PATH = 'projection/A1_wfmash/synteny_rerun/PvP01_to_PAM.tuned.cleaned.chain'


def parse_chains_with_blocks(path, targets, min_density=0.5):
    target_names = {t[0] for t in targets}
    chains = []
    cur_header = None
    cur_blocks = []  # list of (size, dt, dq); last entry just (size,)
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip('\n')
            if ln.startswith('chain '):
                if cur_header:
                    chains.append((cur_header, cur_blocks))
                f = ln.split()
                cur_header = dict(score=int(f[1]), tName=f[2], tSize=int(f[3]),
                                  tStrand=f[4], tStart=int(f[5]), tEnd=int(f[6]),
                                  qName=f[7], qSize=int(f[8]), qStrand=f[9],
                                  qStart=int(f[10]), qEnd=int(f[11]))
                cur_blocks = []
            elif ln.strip() and cur_header is not None:
                parts = ln.split('\t')
                if len(parts) == 1 and parts[0].isdigit():
                    cur_blocks.append((int(parts[0]),))
                elif len(parts) == 3:
                    cur_blocks.append((int(parts[0]), int(parts[1]), int(parts[2])))
        if cur_header:
            chains.append((cur_header, cur_blocks))

    out = []
    for header, blocks in chains:
        if header['tName'] not in target_names:
            continue
        aligned = sum(b[0] for b in blocks)
        span = header['tEnd'] - header['tStart']
        density = aligned / span if span else 0
        if density < min_density:
            continue
        out.append((header, blocks, density))
    return out


def main():
    chains = parse_chains_with_blocks(CHAIN_PATH, TARGETS, min_density=0.5)
    print(f'kept {len(chains)} chains on chr1/chr2 with density ≥50%')

    # Identify PAM contigs that hit chr1 or chr2 (dominant target by aligned_bp)
    contig_target_bp = defaultdict(lambda: defaultdict(int))
    contig_size = {}
    for header, blocks, density in chains:
        contig_target_bp[header['qName']][header['tName']] += sum(b[0] for b in blocks)
        contig_size[header['qName']] = header['qSize']

    contigs = []
    for q, tdict in contig_target_bp.items():
        dom = max(tdict, key=tdict.get)
        contigs.append((q, contig_size[q], dom, tdict[dom]))
    contigs.sort(key=lambda x: ([t[0] for t in TARGETS].index(x[2]), -x[3]))
    print(f'PAM contigs: {len(contigs)}')

    # X axis: chr1 + gap + chr2
    gap = 100_000
    t_offsets = {CHR1[0]: 0, CHR2[0]: CHR1[2] + gap}
    x_total = CHR1[2] + gap + CHR2[2]

    # Y axis: PAM contigs stacked
    y_offsets = {}
    y_total = 0
    y_gap = 50_000
    for q, qsize, dom, _wt in contigs:
        y_offsets[q] = y_total
        y_total += qsize + y_gap

    fig, ax = plt.subplots(figsize=(14, 14))
    chr_colors = {CHR1[0]: '#1f77b4', CHR2[0]: '#ff7f0e'}

    # Draw chromosome boundaries on x axis
    for acc, name, size in TARGETS:
        x0 = t_offsets[acc]
        ax.axvspan(x0, x0 + size, ymin=0, ymax=1, color=chr_colors[acc], alpha=0.05)
        ax.text(x0 + size/2, -y_total*0.025, name,
                ha='center', va='top', fontsize=14, fontweight='bold',
                color=chr_colors[acc])

    # Draw contig boundaries on y axis
    for q, qsize, dom, _wt in contigs:
        y0 = y_offsets[q]
        ax.axhspan(y0, y0 + qsize, xmin=0, xmax=1, color=chr_colors[dom], alpha=0.04)
        short = q if len(q) <= 14 else q[:7] + '…' + q[-6:]
        ax.text(-x_total*0.005, y0 + qsize/2, f'{short} ({qsize/1000:.0f}kb)',
                ha='right', va='center', fontsize=7, color='#333')

    # Draw block lines per chain
    for header, blocks, density in chains:
        t_x_base = t_offsets[header['tName']]
        y_base = y_offsets[header['qName']]
        t_pos = header['tStart']
        # qStart/qEnd: qStart < qEnd is true even when strand=-; for - strand,
        # walk query from qEnd backward
        q_pos_init = header['qStart'] if header['qStrand'] == '+' else (header['qSize'] - header['qStart'])
        q_pos = q_pos_init
        col = chr_colors[header['tName']]
        for b in blocks:
            size = b[0]
            # Draw matching segment
            x0 = t_x_base + t_pos
            x1 = x0 + size
            if header['qStrand'] == '+':
                y0 = y_base + q_pos
                y1 = y_base + q_pos + size
            else:
                # query strand = - → segment goes the other way
                y0 = y_base + (header['qSize'] - q_pos)
                y1 = y_base + (header['qSize'] - (q_pos + size))
            ax.plot([x0, x1], [y0, y1], color=col, linewidth=0.6, alpha=0.7,
                    solid_capstyle='butt')
            # Advance positions
            if len(b) == 3:
                size, dt, dq = b
                t_pos += size + dt
                if header['qStrand'] == '+':
                    q_pos += size + dq
                else:
                    q_pos += size + dq  # still increment by aligned + query gap
            else:
                t_pos += size
                if header['qStrand'] == '+':
                    q_pos += size
                else:
                    q_pos += size

    # Vertical line between chr1 and chr2
    ax.axvline(CHR1[2] + gap/2, color='#888', linestyle='--', linewidth=0.8)

    ax.set_xlim(-x_total*0.05, x_total*1.01)
    ax.set_ylim(-y_total*0.05, y_total*1.01)
    ax.set_xlabel('PvP01 chromosome (bp)', fontsize=12)
    ax.set_ylabel('PAM contigs (concatenated, bp)', fontsize=12)
    ax.set_title('Dot plot: PvP01 chr1 + chr2 vs PAM (KegAlign tuned chains, density ≥50%)\n'
                 'Diagonal = colinear (+); anti-diagonal = inverted (−); colored by PvP01 chrom',
                 fontsize=12)
    plt.tight_layout()
    out = Path('writeup/synteny_3way/dotplot_PvP01_chr1_chr2_vs_PAM.png')
    out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
