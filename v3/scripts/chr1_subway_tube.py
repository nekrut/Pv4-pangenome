#!/usr/bin/env python3
"""
'Subway tube' visualization of PvP01 chr1 with the 2-way PGGB graph:
two horizontal tube-tracks (PvP01 + PAM) that merge into a single fused
backbone where both strains traverse the same nodes, and diverge into
parallel rails where one strain has a bubble (insertion/deletion).
"""
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch

CHR1 = 'LT635612.2'
CHR1_LEN = 1_021_664
GRAPH_2WAY_TSV = 'writeup/synteny_2way/PvP01_to_PAM_2way_chr1.blocks.tsv'


def parse_chr1_blocks(path):
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            t_acc = f[0].split('#')[2] if '#' in f[0] else f[0]
            if t_acc != CHR1:
                continue
            q_acc = f[3].split('#')[2] if '#' in f[3] else f[3]
            yield dict(tStart=int(f[1]), tEnd=int(f[2]), qName=q_acc,
                       qStart=int(f[4]), qEnd=int(f[5]), strand=f[6],
                       score=int(f[7]))


def main():
    blocks = list(parse_chr1_blocks(GRAPH_2WAY_TSV))
    # union of all PAM coverage at 100-bp resolution
    BIN = 100
    n_bins = CHR1_LEN // BIN + 1
    cov_pam = [0] * n_bins
    for b in blocks:
        s = b['tStart'] // BIN
        e = b['tEnd'] // BIN + 1
        for i in range(s, min(e, n_bins)):
            cov_pam[i] = 1

    # Merge consecutive equal-state bins into segments
    segments = []  # list of (start_bp, end_bp, both_strains)
    if not cov_pam:
        return
    cur = cov_pam[0]; seg_start = 0
    for i in range(1, n_bins):
        if cov_pam[i] != cur:
            segments.append((seg_start * BIN, i * BIN, cur))
            seg_start = i
            cur = cov_pam[i]
    segments.append((seg_start * BIN, min(n_bins * BIN, CHR1_LEN), cur))
    print(f'  {len(segments)} alternating-coverage segments')

    # Render
    fig, ax = plt.subplots(figsize=(20, 4))
    PVP01_COLOR = '#1f77b4'
    PAM_COLOR   = '#ff7f0e'
    SHARED_COLOR = '#888'
    TUBE_W = 0.45     # half-thickness of each strain's tube
    MERGE_Y = 0      # y where both tubes merge
    DIVERGE_OFFSET = 1.2  # how far the strains separate at a bubble

    def draw_tube(ax, x0, x1, y, color, thickness=TUBE_W*2, edgecolor='black'):
        ax.add_patch(FancyBboxPatch((x0, y - thickness/2), x1 - x0, thickness,
                                    boxstyle='round,pad=0,rounding_size=0',
                                    facecolor=color, edgecolor=edgecolor,
                                    linewidth=0.5))

    # Iterate segments + draw two tubes
    pv_y = MERGE_Y
    pam_y = MERGE_Y
    last_x = 0
    for s, e, both in segments:
        if both == 1:
            # Both strains share — merged tube
            if pv_y != MERGE_Y or pam_y != MERGE_Y:
                # Draw converging diagonals from previous position
                _draw_bezier(ax, last_x, pv_y, s, MERGE_Y, PVP01_COLOR)
                _draw_bezier(ax, last_x, pam_y, s, MERGE_Y, PAM_COLOR)
                pv_y = MERGE_Y
                pam_y = MERGE_Y
            # Draw fused backbone — wider single tube (no edge color so two strains look fused)
            draw_tube(ax, s, e, MERGE_Y, SHARED_COLOR, thickness=TUBE_W*2.5)
        else:
            # Diverge: PAM doesn't cover; PvP01 keeps going, PAM has nothing
            if pv_y != MERGE_Y or pam_y != -DIVERGE_OFFSET:
                _draw_bezier(ax, last_x, pv_y, s, MERGE_Y, PVP01_COLOR)
                _draw_bezier(ax, last_x, pam_y, s, -DIVERGE_OFFSET, PAM_COLOR)
                pv_y = MERGE_Y
                pam_y = -DIVERGE_OFFSET
            # PvP01 continues
            draw_tube(ax, s, e, MERGE_Y, PVP01_COLOR)
            # PAM is missing — show dashed gap-line
            ax.plot([s, e], [-DIVERGE_OFFSET, -DIVERGE_OFFSET],
                    color=PAM_COLOR, linewidth=1.5, linestyle=(0, (4, 3)),
                    alpha=0.6)
        last_x = e

    # Legend
    ax.text(-180000, MERGE_Y, 'PvP01 / shared', ha='right', va='center',
            fontsize=10, color=PVP01_COLOR, fontweight='bold')
    ax.text(-180000, -DIVERGE_OFFSET, 'PAM (when not shared)', ha='right', va='center',
            fontsize=10, color=PAM_COLOR, fontweight='bold')
    # x-axis
    ax.set_xticks([0, 250_000, 500_000, 750_000, 1_000_000])
    ax.set_xticklabels(['0', '250 kb', '500 kb', '750 kb', '1 Mb'])
    ax.set_xlim(-300000, CHR1_LEN + 50000)
    ax.set_ylim(-2.5, 1.5)
    ax.set_yticks([])
    ax.set_xlabel('PvP01 chr1 position (bp)', fontsize=11)
    for sp in ('top', 'left', 'right'):
        ax.spines[sp].set_visible(False)
    ax.set_title('Subway map of PvP01 chr1 in the 2-way PGGB graph — '
                 'gray tube = shared backbone, blue = PvP01-only segments, '
                 'orange dashes = PAM gaps relative to PvP01',
                 fontsize=11)
    plt.tight_layout()
    out = Path('writeup/graph_viz/chr1_2way_subway.png')
    plt.savefig(out, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out}')


def _draw_bezier(ax, x0, y0, x1, y1, color):
    """Smooth diverge/converge curve between two y levels."""
    if x0 == x1:
        return
    mid_x = (x0 + x1) / 2
    verts = [(x0, y0), (mid_x, y0), (mid_x, y1), (x1, y1)]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    p = MPath(verts, codes)
    ax.add_patch(PathPatch(p, edgecolor=color, fill=False, linewidth=2, alpha=0.7))


if __name__ == '__main__':
    main()
