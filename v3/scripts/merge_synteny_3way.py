#!/usr/bin/env python3
"""Stack chain/PAF/graph synteny plots into one PNG per pair."""
import sys
from pathlib import Path
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    sys.exit("pip install Pillow")

strains = ['PvW1', 'PAM', 'Sal-I', 'PvT01', 'PvC01', 'PvSY56', 'MHC087']
chain_dir = Path('writeup/synteny')
paf_dir = Path('writeup/synteny_paf')
graph_dir = Path('writeup/synteny_graph')
out_dir = Path('writeup/synteny_3way')
out_dir.mkdir(parents=True, exist_ok=True)

try:
    font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 28)
except OSError:
    font = ImageFont.load_default()

n_done = 0
n_partial = 0
for s in strains:
    p_chain = chain_dir / f'synteny_PvP01_to_{s}.png'
    p_paf = paf_dir / f'paf_synteny_PvP01_to_{s}.png'
    p_graph = graph_dir / f'graph_synteny_PvP01_to_{s}.png'
    available = [(p, lbl) for p, lbl in
                 [(p_chain, 'KegAlign chains (Phase B)'),
                  (p_paf, 'wfmash PAF (raw A1)'),
                  (p_graph, 'PGGB graph blocks (shared-node runs)')]
                 if p.exists()]
    if not available:
        print(f'  skip {s} — no images')
        continue
    # Open all and align widths to the widest
    images = [(Image.open(p).convert('RGB'), lbl) for p, lbl in available]
    target_w = max(im.width for im, _ in images)
    images = [(im if im.width == target_w
               else im.resize((target_w, int(im.height * target_w / im.width)),
                              Image.LANCZOS), lbl)
              for im, lbl in images]
    banner_h = 50
    pad = 20
    canvas_h = sum(im.height + banner_h + pad for im, _ in images) - pad
    canvas = Image.new('RGB', (target_w, canvas_h), 'white')
    draw = ImageDraw.Draw(canvas)
    y = 0
    for im, lbl in images:
        draw.text((20, y + 10), f'{lbl}: PvP01 → {s}', fill='black', font=font)
        y += banner_h
        canvas.paste(im, (0, y))
        y += im.height + pad
    out = out_dir / f'PvP01_to_{s}_3way.png'
    canvas.save(out, 'PNG', optimize=True)
    if len(images) == 3:
        print(f'  wrote {out}  (3-way)')
        n_done += 1
    else:
        print(f'  wrote {out}  (PARTIAL — only {len(images)} sources)')
        n_partial += 1

print(f'\n{n_done} complete 3-way + {n_partial} partial')
