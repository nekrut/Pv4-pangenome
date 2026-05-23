#!/usr/bin/env python3
"""
Base-resolution subway map of a chr1 region.

Parse the GFA subgraph; lay out nodes by their PvP01 chr1 position; for each
node show its sequence; draw paths as colored lines threading through nodes.
Where PvP01 and PAM share a node → lines overlap. Where they differ → lines
branch into separate parallel nodes (SNP) or one path skips a node (indel).
"""
import sys
from pathlib import Path
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D

GFA = 'cactus_2way/chr1_zoom/chr1_4200_4600.gfa'
PVP01_PATH = 'GCA_900093555.2#1#LT635612.2#0:4200-4600'
CHR1_START = 4200  # offset for x-axis labels


def parse_gfa(path):
    nodes = {}  # nid (str) → seq
    edges = set()  # (a, ori, b, ori)
    paths = {}  # name → list of (nid, ori)
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if not f: continue
            if f[0] == 'S':
                nodes[f[1]] = f[2]
            elif f[0] == 'L':
                edges.add((f[1], f[2], f[3], f[4]))
            elif f[0] == 'P':
                name = f[1]
                steps = []
                for tok in f[2].split(','):
                    if not tok: continue
                    steps.append((tok[:-1], tok[-1]))
                paths[name] = steps
    return nodes, edges, paths


def main():
    nodes, edges, paths = parse_gfa(GFA)
    print(f'  {len(nodes)} nodes, {len(edges)} edges, {len(paths)} paths')

    # PvP01 traversal: compute x-position (bp) for each PvP01 node
    pv_steps = paths.get(PVP01_PATH, [])
    pv_node_pos = {}    # node_id → (start_bp, end_bp, orientation)
    pos = CHR1_START
    for nid, ori in pv_steps:
        seq = nodes[nid]
        pv_node_pos[nid] = (pos, pos + len(seq), ori)
        pos += len(seq)
    pv_end = pos
    print(f'  PvP01 path covers {CHR1_START}-{pv_end} bp ({pv_end-CHR1_START} bp)')

    # PAM paths
    pam_paths = {n: s for n, s in paths.items() if 'GCA_949152365.1' in n}

    # For non-PvP01 nodes, assign x by interpolating from flanking PvP01 nodes
    # found in the PAM path
    def x_of(node_id, fallback_pos=None):
        if node_id in pv_node_pos:
            return pv_node_pos[node_id][0]
        return fallback_pos

    # Per PAM path: compute x-positions of its nodes
    pam_node_x = defaultdict(dict)  # pam_path → node → x
    for pname, steps in pam_paths.items():
        # Walk; for nodes in PvP01, use that x; for others, interpolate between flanks
        last_anchor_x = CHR1_START
        # First pass: find anchors (PvP01-shared nodes)
        anchor_indices = [(i, pv_node_pos[nid][0]) for i, (nid, _) in enumerate(steps) if nid in pv_node_pos]
        for i, (nid, ori) in enumerate(steps):
            if nid in pv_node_pos:
                pam_node_x[pname][nid] = pv_node_pos[nid][0]
            else:
                # Find nearest anchors before and after
                prev_anchor_x = next((a[1] for a in reversed(anchor_indices) if a[0] < i), CHR1_START)
                next_anchor_x = next((a[1] for a in anchor_indices if a[0] > i), pv_end)
                pam_node_x[pname][nid] = (prev_anchor_x + next_anchor_x) / 2

    # Plotting
    fig, ax = plt.subplots(figsize=(28, 6))
    Y_PVP01 = 0
    Y_PAM_OFFSETS = {p: -2 - i*2 for i, p in enumerate(pam_paths)}
    NODE_H = 0.7
    NODE_PAD = 0.05

    base_colors = {'A': '#3c8bd0', 'C': '#e67e22', 'G': '#27ae60', 'T': '#c0392b', 'N': '#888'}

    # Draw PvP01 nodes
    for nid, ori in pv_steps:
        seq = nodes[nid]
        x0, x1, _ = pv_node_pos[nid]
        for i, base in enumerate(seq):
            col = base_colors.get(base.upper(), '#999')
            ax.add_patch(Rectangle((x0 + i, Y_PVP01 - NODE_H/2), 1 - NODE_PAD, NODE_H,
                                   facecolor=col, edgecolor='black', linewidth=0.3, alpha=0.85))
            if len(seq) <= 10 or i % 2 == 0:
                ax.text(x0 + i + 0.5, Y_PVP01, base.upper(),
                        ha='center', va='center', fontsize=7, color='white',
                        fontweight='bold')

    # Draw PAM nodes (those not in PvP01) — at offset rows
    for pname, steps in pam_paths.items():
        py = Y_PAM_OFFSETS[pname]
        # Per-node draw
        for nid, ori in steps:
            x = pam_node_x[pname][nid]
            seq = nodes[nid]
            if nid in pv_node_pos:
                # shared node — draw a thin connector tie (vertical line from pv row to pam row at this x)
                # but only one tie per shared node
                pass
            else:
                # PAM-specific node (bubble alt allele or insertion)
                for i, base in enumerate(seq):
                    col = base_colors.get(base.upper(), '#999')
                    ax.add_patch(Rectangle((x + i - 0.5, py - NODE_H/2), 1 - NODE_PAD, NODE_H,
                                           facecolor=col, edgecolor='black', linewidth=0.3, alpha=0.85))
                    ax.text(x + i, py, base.upper(),
                            ha='center', va='center', fontsize=7, color='white',
                            fontweight='bold')

    # Draw path traversal as connecting lines
    # PvP01 backbone line (top)
    pv_xs = []
    for nid, _ in pv_steps:
        x0, x1, _ = pv_node_pos[nid]
        pv_xs.append((x0 + x1) / 2)
    ax.plot(pv_xs, [Y_PVP01]*len(pv_xs), color='#1f4e79', linewidth=1.5, alpha=0.4, zorder=0)
    ax.text(CHR1_START - 5, Y_PVP01, 'PvP01', ha='right', va='center',
            fontsize=11, fontweight='bold', color='#1f4e79')

    # PAM path lines: thread through nodes (shared on top row, bubble nodes on PAM row)
    for pname, steps in pam_paths.items():
        py = Y_PAM_OFFSETS[pname]
        xs = []
        ys = []
        for nid, _ in steps:
            x = pam_node_x[pname][nid]
            seq = nodes[nid]
            y_here = Y_PVP01 if nid in pv_node_pos else py
            xs.append(x + len(seq)/2)
            ys.append(y_here)
        ax.plot(xs, ys, color='#cc6600', linewidth=2, alpha=0.6, zorder=1)
        # Path label
        short = pname.split('#')[2]
        if len(short) > 14:
            short = short[:7] + '…' + short[-6:]
        ax.text(CHR1_START - 5, py, f'PAM ({short})', ha='right', va='center',
                fontsize=10, fontweight='bold', color='#cc6600')

    # x-axis
    ax.set_xlim(CHR1_START - 30, pv_end + 10)
    ax.set_ylim(min(Y_PAM_OFFSETS.values()) - 1.5, Y_PVP01 + 1.0)
    xticks = list(range(CHR1_START, pv_end + 1, 50))
    ax.set_xticks(xticks)
    ax.set_xticklabels([f'{t}' for t in xticks], fontsize=8)
    ax.set_yticks([])
    ax.set_xlabel('PvP01 chr1 position (bp)', fontsize=11)
    for sp in ('top', 'left', 'right'):
        ax.spines[sp].set_visible(False)

    # Legend for bases
    legend_handles = [Line2D([0], [0], marker='s', color='w',
                              markerfacecolor=base_colors[b], markersize=10,
                              label=b) for b in ['A','C','G','T']]
    ax.legend(handles=legend_handles, loc='lower right', ncol=4,
              frameon=False, fontsize=10)

    ax.set_title(f'PGGB 2-way graph — base-resolution subway map for PvP01 chr1:{CHR1_START}-{pv_end} '
                 '(PvP01 row top, PAM contigs below; nodes labeled by base; '
                 'bubbles where the strains differ)',
                 fontsize=11)
    plt.tight_layout()
    out = Path('writeup/graph_viz/chr1_zoom_subway.png')
    plt.savefig(out, dpi=140, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
