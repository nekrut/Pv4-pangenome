#!/usr/bin/env python3
"""
Re-render the 5 canonical-layout panels for PvP01 → PAM with:
  - Sparse-chain filter (only blocks with ≥50% alignment density survive in chain tracks)
  - PAM scaffold labels below each contig rectangle
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
    PVP01_CHRS, PVP01_ORDER, CHR_COLOR, SOURCES, strip_pansn,
    parse_chain, parse_paf, parse_graph, normalize_qname, make_q_offsets,
    bezier_ribbon, build_canonical_layout)


def chain_with_density(path, min_density=0.5):
    """parse_chain + compute aligned_bp; filter sparse chains."""
    chains = []
    cur = None
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip('\n')
            if ln.startswith('chain '):
                if cur:
                    chains.append(cur)
                f = ln.split()
                cur = dict(score=int(f[1]), tName=f[2], tStart=int(f[5]), tEnd=int(f[6]),
                           qName=f[7], qSize=int(f[8]), qStart=int(f[10]), qEnd=int(f[11]),
                           strand=f[9], aligned_bp=0)
            elif ln.strip() and cur is not None:
                parts = ln.split('\t')
                if parts and parts[0].isdigit():
                    cur['aligned_bp'] += int(parts[0])
        if cur:
            chains.append(cur)
    for c in chains:
        t_span = c['tEnd'] - c['tStart']
        c['density'] = c['aligned_bp'] / t_span if t_span > 0 else 0
    return [c for c in chains if c['density'] >= min_density]


def render(source_key, kind, path, t_offsets, q_offsets, q_dom, out_png, title,
           min_chain_density=0.5):
    if kind == 'chain':
        blocks = chain_with_density(path, min_density=min_chain_density)
    elif kind == 'paf':
        blocks = list(parse_paf(path))
    elif kind == 'graph':
        blocks = list(parse_graph(path))
    valid_q = set(q_offsets)
    blocks = [b for b in blocks
              if b['tName'] in CHR_COLOR and normalize_qname(b['qName'], kind) in valid_q]

    fig, ax = plt.subplots(figsize=(20, 7.2))
    ROW_H = 100000; yT = 700000; yB = 0
    t_total = max(t_offsets.values()) + max(PVP01_CHRS[a][1] for a in PVP01_ORDER)
    for acc in PVP01_ORDER:
        size = PVP01_CHRS[acc][1]; x0 = t_offsets[acc]; col = CHR_COLOR[acc]
        ax.add_patch(Rectangle((x0, yT), size, ROW_H, facecolor=col, edgecolor='black',
                               linewidth=0.5, alpha=0.9))
        ax.text(x0 + size/2, yT + ROW_H + 30000, PVP01_CHRS[acc][0],
                ha='center', va='bottom', fontsize=9, color=col, fontweight='bold')
    # Bottom contigs + labels
    for q, (x0, w, real) in q_offsets.items():
        dom = q_dom.get(q)
        col = CHR_COLOR.get(dom, 'gray')
        ax.add_patch(Rectangle((x0, yB), w, ROW_H, facecolor=col, edgecolor='black',
                               linewidth=0.3, alpha=0.7))
        # Label: scaffold name + size in kb, rotated 90° for compactness
        n_kb = real / 1000
        # short scaffold id: CASCJQ010000001.1 → CASCJQ...001
        short = q
        if len(short) > 12:
            short = short[:7] + '…' + short[-6:]
        label = f'{short} ({n_kb:.0f}kb)'
        ax.text(x0 + w/2, yB - 40000, label,
                ha='center', va='top', fontsize=7, color='#333',
                rotation=-30, rotation_mode='anchor')
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
    ax.set_ylim(-700000, yT + ROW_H + 200000)
    ax.set_yticks([]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title, fontsize=13)
    plt.tight_layout(); plt.savefig(out_png, dpi=140, bbox_inches='tight'); plt.close(fig)


def main():
    layout = build_canonical_layout(min_block=1000)
    print(f'canonical layout: {len(layout)} PAM contigs')
    gap = 100_000
    t_offsets = {}; x = 0
    for acc in PVP01_ORDER:
        t_offsets[acc] = x; x += PVP01_CHRS[acc][1] + gap
    q_offsets = make_q_offsets(layout, t_offsets)
    q_dom = {q: dom for q, _, dom, _ in layout}

    titles = {
        '1_kegalign_default': 'KegAlign chains, default (HoxD70, hspthresh 3000, seed 12of19) — density ≥50%',
        '2_kegalign_tuned':   'KegAlign chains, TUNED (matrix +100/-100, hspthresh 4500, seed 14of22) — density ≥50%',
        '3_wfmash_n1':        'wfmash PAF (A1 params: -n 1, best-per-segment)',
        '4_wfmash_pggblike':  'wfmash PAF (PGGB-build params: -n 8, multi-mapping)',
        '5_graph_blocks':     'PGGB graph blocks (shared-node runs)',
    }
    outdir = Path('writeup/synteny_canonical_filtered')
    outdir.mkdir(parents=True, exist_ok=True)
    for k, (kind, path) in SOURCES.items():
        out = outdir / f'PvP01_to_PAM__{k}.png'
        render(k, kind, path, t_offsets, q_offsets, q_dom, str(out),
               f'{titles[k]}: PvP01 → PAM',
               min_chain_density=0.5)
        print(f'  wrote {out}')


if __name__ == '__main__':
    main()
