#!/usr/bin/env python3
"""Stack chain-based and PAF-based synteny plots into one PNG per pair."""
import sys
from pathlib import Path
try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    sys.exit("pip install Pillow")

strains = ['PvW1', 'PAM', 'Sal-I', 'PvT01', 'PvC01', 'PvSY56', 'MHC087']
chain_dir = Path('writeup/synteny')
paf_dir = Path('writeup/synteny_paf')
out_dir = Path('writeup/synteny_compare')
out_dir.mkdir(parents=True, exist_ok=True)

try:
    font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 28)
except OSError:
    font = ImageFont.load_default()

for s in strains:
    p_chain = chain_dir / f'synteny_PvP01_to_{s}.png'
    p_paf = paf_dir / f'paf_synteny_PvP01_to_{s}.png'
    if not p_chain.exists() or not p_paf.exists():
        print(f'  skip {s}')
        continue
    im_c = Image.open(p_chain).convert('RGB')
    im_p = Image.open(p_paf).convert('RGB')
    # Scale paf image to chain width if different
    if im_p.width != im_c.width:
        ratio = im_c.width / im_p.width
        im_p = im_p.resize((im_c.width, int(im_p.height * ratio)), Image.LANCZOS)
    # Make banner per panel
    banner_h = 50
    canvas = Image.new('RGB', (im_c.width, banner_h*2 + im_c.height + im_p.height + 20), 'white')
    draw = ImageDraw.Draw(canvas)
    y = 0
    draw.text((20, y + 10), f'TOP — KegAlign chains (Phase B): PvP01 → {s}', fill='black', font=font)
    y += banner_h
    canvas.paste(im_c, (0, y))
    y += im_c.height + 20
    draw.text((20, y + 10), f'BOTTOM — wfmash PAF (raw A1): PvP01 → {s}', fill='black', font=font)
    y += banner_h
    canvas.paste(im_p, (0, y))
    out = out_dir / f'PvP01_to_{s}_compare.png'
    canvas.save(out, 'PNG', optimize=True)
    print(f'  wrote {out}')
