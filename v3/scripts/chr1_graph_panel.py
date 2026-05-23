#!/usr/bin/env python3
"""
Custom graph-coverage panel for PvP01 chr1, bp-aligned with the synteny ribbons.

Reads the graph_synteny_blocks.py TSV (which has tStart/tEnd as PvP01 chr1 bp
coordinates), draws:
  - Top row: PvP01 chr1 reference bar (always covered)
  - Bottom: 2-strain "depth" track — for each chr1 position, color based on
    whether 1 (strain-specific bubble) or 2 strains traverse the node
  - One row per PAM contig that hits chr1 (sorted by chr1 alignment weight),
    with bp-aligned bars showing where they cover chr1
"""
import sys
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

CHR1 = 'LT635612.2'
CHR1_LEN = 1_021_664
GRAPH_2WAY_TSV = 'writeup/synteny_2way/PvP01_to_PAM_2way_chr1.blocks.tsv'


def parse_chr1_blocks(path):
    """Yield blocks: {tStart, tEnd, qName, qStart, qEnd, strand, score}.
    PanSN-strip qName."""
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            t_acc = f[0].split('#')[2] if '#' in f[0] else f[0]
            if t_acc != CHR1:
                continue
            q_acc = f[3].split('#')[2] if '#' in f[3] else f[3]
            yield dict(tStart=int(f[1]), tEnd=int(f[2]),
                       qName=q_acc, qStart=int(f[4]), qEnd=int(f[5]),
                       strand=f[6], score=int(f[7]))


def main():
    blocks = list(parse_chr1_blocks(GRAPH_2WAY_TSV))
    print(f'  {len(blocks)} chr1 blocks (PAM only, 2-way graph)')
    # Sort PAM contigs by total chr1 alignment weight
    pam_weight = defaultdict(int)
    for b in blocks:
        pam_weight[b['qName']] += b['tEnd'] - b['tStart']
    contigs = sorted(pam_weight, key=lambda q: -pam_weight[q])
    print(f'  {len(contigs)} PAM contigs touching chr1')

    # Build per-bp coverage (in 1kb bins for speed)
    BIN = 1000
    n_bins = CHR1_LEN // BIN + 1
    # PvP01 always covers entire chr1
    cov_pvp01 = [1] * n_bins
    cov_pam = [0] * n_bins
    pam_by_contig = {q: [0] * n_bins for q in contigs}
    for b in blocks:
        s = b['tStart'] // BIN
        e = b['tEnd'] // BIN + 1
        for i in range(s, min(e, n_bins)):
            cov_pam[i] = 1
            pam_by_contig[b['qName']][i] = 1

    # Plot
    n_rows = 2 + 2 + len(contigs)  # PvP01 ref + spacer + 2-strain depth + per-contig rows
    fig_h = max(3, 0.4 * n_rows + 1.5)
    fig, ax = plt.subplots(figsize=(20, fig_h))
    PVP01_COLOR = '#1f77b4'
    PAM_COLOR = '#ff7f0e'
    SHARED_COLOR = '#666666'
    BUBBLE_COLOR = '#ff7f0e'  # PAM-only bubble

    y = 0
    row_h = 0.8

    # Row: PvP01 chr1 (always full coverage)
    ax.add_patch(Rectangle((0, y), CHR1_LEN, row_h, facecolor=PVP01_COLOR,
                           edgecolor='black', linewidth=0.5))
    ax.text(-25000, y + row_h/2, 'PvP01 chr1', ha='right', va='center',
            fontsize=10, fontweight='bold', color=PVP01_COLOR)
    y += row_h + 0.4

    # Row: 2-strain depth (PvP01 + any PAM contig)
    for i, c in enumerate(cov_pam):
        x_start = i * BIN
        x_end = min((i + 1) * BIN, CHR1_LEN)
        if c >= 1:
            ax.add_patch(Rectangle((x_start, y), x_end - x_start, row_h,
                                   facecolor=SHARED_COLOR, edgecolor='none'))
        else:
            ax.add_patch(Rectangle((x_start, y), x_end - x_start, row_h,
                                   facecolor='#eee', edgecolor='none'))
    ax.text(-25000, y + row_h/2, 'Depth (any PAM)', ha='right', va='center',
            fontsize=10, fontweight='bold', color=SHARED_COLOR)
    y += row_h + 0.4

    # Spacer / divider
    ax.axhline(y - 0.1, color='#888', linestyle='--', linewidth=0.5)
    ax.text(-25000, y - 0.05, 'per-contig:', ha='right', va='top',
            fontsize=9, color='#666', style='italic')
    y += 0.3

    # Per-PAM-contig rows
    for q in contigs:
        cov = pam_by_contig[q]
        for i, c in enumerate(cov):
            if c < 1:
                continue
            x_start = i * BIN
            x_end = min((i + 1) * BIN, CHR1_LEN)
            ax.add_patch(Rectangle((x_start, y), x_end - x_start, row_h,
                                   facecolor=PAM_COLOR, edgecolor='none'))
        short = q if len(q) <= 14 else q[:7] + '…' + q[-6:]
        wt_kb = pam_weight[q] / 1000
        ax.text(-25000, y + row_h/2, f'{short} ({wt_kb:.0f}kb)',
                ha='right', va='center', fontsize=8, color='#333')
        y += row_h + 0.2

    ax.set_xlim(-180000, CHR1_LEN + 50000)
    ax.set_ylim(y + 0.5, -0.5)  # invert so PvP01 is on top
    ax.set_yticks([]); ax.set_xticks([0, 250000, 500000, 750000, 1_000_000])
    ax.set_xticklabels(['0', '250 kb', '500 kb', '750 kb', '1 Mb'])
    ax.set_xlabel('PvP01 chr1 position (bp)', fontsize=11)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.set_title('PvP01 chr1: per-contig graph coverage from 2-way PGGB '
                 '(rows = PAM contigs; orange = PAM covers this 1 kb bin; '
                 'gray = PvP01 only)',
                 fontsize=11)
    plt.tight_layout()
    out = Path('writeup/graph_viz/chr1_2way_coverage_bp.png')
    plt.savefig(out, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
