#!/usr/bin/env python3
"""Stack chr1-focused panels into writeup/synteny_3way/chr1_PvP01_PAM_graph.png.

Order:
  1. KegAlign chains, default
  2. KegAlign chains, TUNED
  3. wfmash PAF (PGGB params, -n 8)
  4. PGGB graph blocks (8-way)
  5. PGGB graph blocks (2-way, PvP01+PAM)         ← new ribbon panel
  6. PGGB 2-way coverage on PvP01 chr1, bp-aligned
  7. Subway map (PvP01/shared vs PAM-when-not-shared)
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None

PANELS = [
    ('writeup/synteny_3way/chr1_zoomed_panels/1_kegalign_default.png',
     'KegAlign chains, default (density ≥50%)'),
    ('writeup/synteny_3way/chr1_zoomed_panels/2_kegalign_tuned.png',
     'KegAlign chains, TUNED (density ≥50%)'),
    ('writeup/synteny_3way/chr1_zoomed_panels/4_wfmash_pggblike.png',
     'wfmash PAF (PGGB params: -n 8)'),
    ('writeup/synteny_3way/chr1_zoomed_panels/5_graph_blocks.png',
     'PGGB graph blocks (8-way, all strains)'),
    ('writeup/synteny_3way/chr1_zoomed_panels/6_graph_blocks_2way.png',
     'PGGB graph blocks (2-way, PvP01 + PAM only)'),
    ('writeup/graph_viz/chr1_2way_coverage_bp.png',
     'PGGB 2-way graph: per-contig coverage on PvP01 chr1, bp-aligned'),
    ('writeup/graph_viz/chr1_2way_subway.png',
     'Subway map: gray = shared backbone, blue = PvP01-only, dashed orange = PAM gap'),
]

try:
    title_font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 32)
except OSError:
    title_font = ImageFont.load_default()

BANNER_H = 60
PAD = 20

imgs = []
for p, title in PANELS:
    im = Image.open(p).convert('RGB')
    imgs.append((im, title))

max_w = max(im.width for im, _ in imgs)
# Scale all images to the same width (no upscaling)
scaled = []
for im, title in imgs:
    if im.width != max_w:
        new_h = int(im.height * (max_w / im.width))
        im = im.resize((max_w, new_h), Image.LANCZOS)
    scaled.append((im, title))

total_h = sum(im.height for im, _ in scaled) + (BANNER_H + PAD) * len(scaled) + PAD
out = Image.new('RGB', (max_w, total_h), 'white')
draw = ImageDraw.Draw(out)
y = PAD
for im, title in scaled:
    draw.text((PAD, y), title, fill='black', font=title_font)
    y += BANNER_H
    out.paste(im, (0, y))
    y += im.height + PAD

out_path = Path('writeup/synteny_3way/chr1_PvP01_PAM_graph.png')
out.save(out_path, dpi=(150, 150))
print(f'wrote {out_path} {out.size}')
