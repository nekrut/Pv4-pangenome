#!/usr/bin/env python3
"""
Make a self-contained SVG showing 5 synteny tracks on PvP01 chr 1.
SVG is vector → infinitely zoomable in any browser. No server needed.
"""
from pathlib import Path
import html

CHR_LEN = 1_021_664  # LT635612.2

tracks = [
    ('1. KegAlign default (HoxD70, hspthresh 3000)',                'writeup/browser/tracks/1_kegalign_default.bed', '#1f77b4'),
    ('2. KegAlign tuned (custom +100/-100, hspthresh 4500, 14of22)', 'writeup/browser/tracks/2_kegalign_tuned.bed',   '#ff7f0e'),
    ('3. wfmash (-n 1, A1 chain-gen)',                              'writeup/browser/tracks/3_wfmash_n1.bed',        '#2ca02c'),
    ('4. wfmash (-n 8, PGGB build)',                                'writeup/browser/tracks/4_wfmash_pggblike.bed',  '#d62728'),
    ('5. PGGB graph blocks (shared-node runs)',                      'writeup/browser/tracks/5_graph_blocks.bed',    '#9467bd'),
]

# Drawing
WIDTH = 1800
LEFT = 240; RIGHT = 60
TRACK_W = WIDTH - LEFT - RIGHT
TRACK_H = 60
TRACK_GAP = 30
TOP = 100
N_TRACKS = len(tracks)
HEIGHT = TOP + N_TRACKS * (TRACK_H + TRACK_GAP) + 100


def x_of(bp):
    return LEFT + bp / CHR_LEN * TRACK_W


svg_parts = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
    f'viewBox="0 0 {WIDTH} {HEIGHT}" style="background:#fafafa;font-family:ui-sans-serif,sans-serif">',
    f'<text x="{LEFT}" y="40" font-size="20" font-weight="bold">PvP01 chr 1 (LT635612.2, 1,021,664 bp) → PAM synteny — 5 sources</text>',
    f'<text x="{LEFT}" y="62" font-size="12" fill="#555">Each rectangle = one alignment block ≥1 kb. Color by track. Stroke for − strand.</text>',
]

# X-axis ticks
for tick in range(0, CHR_LEN + 1, 100_000):
    x = x_of(tick)
    label_mb = f'{tick/1_000_000:.1f} Mb'
    svg_parts.append(f'<line x1="{x:.1f}" y1="{TOP-12}" x2="{x:.1f}" y2="{TOP-2}" stroke="#888" stroke-width="1"/>')
    svg_parts.append(f'<text x="{x:.1f}" y="{TOP-18}" font-size="11" text-anchor="middle" fill="#666">{label_mb}</text>')

# Draw each track
for i, (label, bed_path, color) in enumerate(tracks):
    y0 = TOP + i * (TRACK_H + TRACK_GAP)
    # Track baseline
    svg_parts.append(f'<line x1="{LEFT}" y1="{y0+TRACK_H/2}" x2="{x_of(CHR_LEN):.1f}" y2="{y0+TRACK_H/2}" stroke="#ddd" stroke-width="2"/>')
    # Label
    svg_parts.append(f'<text x="{LEFT-10}" y="{y0+TRACK_H/2+4}" font-size="13" font-weight="bold" '
                     f'text-anchor="end" fill="{color}">{html.escape(label)}</text>')
    # Boxes
    if Path(bed_path).exists():
        with open(bed_path) as fh:
            for ln in fh:
                f = ln.rstrip('\n').split('\t')
                if len(f) < 6: continue
                start, end = int(f[1]), int(f[2])
                name, score_s, strand = f[3], f[4], f[5]
                x0 = x_of(start)
                w  = max(1.0, x_of(end) - x_of(start))
                # Score → opacity 0.3-1.0
                try:
                    sc = int(score_s)
                except ValueError:
                    sc = 500
                opac = 0.3 + (sc / 1000) * 0.7
                stroke = '#000' if strand == '-' else 'none'
                svg_parts.append(
                    f'<rect x="{x0:.2f}" y="{y0+8}" width="{w:.2f}" height="{TRACK_H-16}" '
                    f'fill="{color}" fill-opacity="{opac:.2f}" stroke="{stroke}" stroke-width="0.8">'
                    f'<title>{html.escape(name)} {start}-{end} strand={strand}</title></rect>'
                )

svg_parts.append('</svg>')

out = Path('writeup/browser/chr1_tracks.svg')
out.write_text('\n'.join(svg_parts))
print(f'wrote {out} ({out.stat().st_size} bytes)')

# Also build a tiny HTML wrapper for convenience
html_out = Path('writeup/browser/chr1_tracks.html')
html_out.write_text(f'''<!doctype html>
<html><head><meta charset="utf-8"><title>PvP01 chr 1 synteny tracks</title>
<style>body {{ margin:0; padding:16px; background:#fff; font-family:ui-sans-serif,sans-serif; }} svg {{ max-width:100%; height:auto; }}</style>
</head><body>
<p style="font-size:13px;color:#555;margin:0 0 8px 0">
Hover blocks for query coordinates. Zoom via Ctrl/Cmd+wheel or browser zoom. <a href="index.html">IGV.js version</a>.
</p>
{out.read_text()}
</body></html>''')
print(f'wrote {html_out}')
