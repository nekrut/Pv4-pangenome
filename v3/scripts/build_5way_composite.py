#!/usr/bin/env python3
"""Stack the 6 PvP01 → PAM whole-genome ribbon panels into PvP01_to_PAM_5way.png
(renamed-but-still-called-5way for continuity; now 6 panels)."""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
Image.MAX_IMAGE_PIXELS = None

PANELS = [
    ('writeup/synteny_canonical/PvP01_to_PAM__1_kegalign_default.png',
     'KegAlign chains, default (HoxD70, hspthresh 3000, seed 12of19)'),
    ('writeup/synteny_canonical/PvP01_to_PAM__2_kegalign_tuned.png',
     'KegAlign chains, TUNED (matrix +100/-100, hspthresh 4500, seed 14of22)'),
    ('writeup/synteny_canonical/PvP01_to_PAM__4_wfmash_pggblike.png',
     'wfmash PAF (PGGB-build params: -n 8, multi-mapping)'),
    ('writeup/synteny_canonical/PvP01_to_PAM__6_graph_blocks_2way.png',
     'PGGB graph blocks, 2-way (PvP01 + PAM only)'),
]

try:
    title_font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf', 32)
except OSError:
    title_font = ImageFont.load_default()

BANNER_H = 60
PAD = 20

imgs = [(Image.open(p).convert('RGB'), t) for p, t in PANELS]
max_w = max(im.width for im, _ in imgs)
scaled = []
for im, t in imgs:
    if im.width != max_w:
        new_h = int(im.height * (max_w / im.width))
        im = im.resize((max_w, new_h), Image.LANCZOS)
    scaled.append((im, t))

total_h = sum(im.height for im, _ in scaled) + (BANNER_H + PAD) * len(scaled) + PAD
out = Image.new('RGB', (max_w, total_h), 'white')
draw = ImageDraw.Draw(out)
y = PAD
for im, t in scaled:
    draw.text((PAD, y), t, fill='black', font=title_font)
    y += BANNER_H
    out.paste(im, (0, y))
    y += im.height + PAD

p = Path('writeup/synteny_3way/PvP01_to_PAM_5way.png')
out.save(p, dpi=(150, 150))
print(f'wrote {p} {out.size}')
