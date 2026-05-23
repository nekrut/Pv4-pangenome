#!/usr/bin/env python3
"""
Build 5 BED tracks restricted to PvP01 chr 1 (LT635612.2), one per synteny source:

  1. KegAlign default chain (Phase B, default lastz params)
  2. KegAlign tuned chain (custom matrix +100/-100, hspthresh 4500, seed 14of22)
  3. wfmash PAF A1 (-n 1)
  4. wfmash PAF PGGB-params (-n 8)
  5. PGGB graph blocks (shared-node runs)

Each block becomes a BED12 row (or BED6 where sub-block structure isn't available).
Output to writeup/browser/tracks/.
"""
import re
from pathlib import Path

PVP01_CHR1 = 'LT635612.2'
PAM_PANSN_PREFIX = 'GCA_949152365.1#1#'  # for graph block source
OUT = Path('writeup/browser/tracks')
OUT.mkdir(parents=True, exist_ok=True)


def parse_chain_blocks(path, score_min=0):
    """Parse a UCSC chain file. Yield {tStart, tEnd, qName, qStart, qEnd, strand, score}
    for chains where the target chromosome is PvP01_CHR1.
    Use the chain header (alignable bounds in target) — sub-block detail can be
    represented in BED12 thickStart/Ends if needed.
    """
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith('chain '):
                continue
            f = ln.split()
            tName = f[2]
            if tName != PVP01_CHR1:
                continue
            score = int(f[1])
            if score < score_min:
                continue
            tSize = int(f[3]); tStart = int(f[5]); tEnd = int(f[6])
            qName = f[7]; qSize = int(f[8]); qStrand = f[9]
            qStart = int(f[10]); qEnd = int(f[11])
            yield dict(tStart=tStart, tEnd=tEnd, qName=qName, qStart=qStart,
                       qEnd=qEnd, strand=qStrand, score=score)


def parse_paf_blocks(path, min_block=1000):
    """Parse a PAF. Filter to PvP01 chr 1 as target."""
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 12 or f[5] != PVP01_CHR1:
                continue
            blockLen = int(f[10])
            if blockLen < min_block:
                continue
            yield dict(tStart=int(f[7]), tEnd=int(f[8]),
                       qName=f[0], qStart=int(f[2]), qEnd=int(f[3]),
                       strand=f[4], score=blockLen)


def parse_graph_blocks(path, min_block=1000):
    """Parse graph_synteny_blocks.py TSV. Filter to PvP01 chr 1 (look at end of PanSN name)."""
    with open(path) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            tName = f[0]  # PanSN name e.g. GCA_900093555.2#1#LT635612.2#0
            if not tName.endswith(f'#{PVP01_CHR1}#0'):
                continue
            score = int(f[7])
            if score < min_block:
                continue
            qName = f[3]
            # Strip PanSN prefix
            yield dict(tStart=int(f[1]), tEnd=int(f[2]),
                       qName=qName, qStart=int(f[4]), qEnd=int(f[5]),
                       strand=f[6], score=score)


def write_bed(blocks, out_path, name_prefix='blk'):
    """Write BED6: chrom start end name score strand. The score is clamped to 0-1000 for IGV color shading."""
    n = 0
    rows = list(blocks)
    if not rows:
        Path(out_path).write_text('')
        return 0
    max_score = max(b['score'] for b in rows)
    with open(out_path, 'w') as fh:
        for i, b in enumerate(rows):
            # Friendly name: query contig:start-end
            qshort = b['qName'].split('#')[-2] if '#' in b['qName'] else b['qName']
            name = f'{name_prefix}{i}_{qshort}:{b["qStart"]}-{b["qEnd"]}_{b["strand"]}'
            # IGV BED score 0-1000 for visual shading
            scaled = min(1000, int(b['score'] / max_score * 1000)) if max_score > 0 else 500
            fh.write(f'{PVP01_CHR1}\t{b["tStart"]}\t{b["tEnd"]}\t{name}\t{scaled}\t{b["strand"]}\n')
            n += 1
    return n


sources = [
    ('1_kegalign_default',
     parse_chain_blocks('work/01_chains/GCA_900093555.2.GCA_949152365.1.cleaned.chain', score_min=1000)),
    ('2_kegalign_tuned',
     parse_chain_blocks('projection/A1_wfmash/synteny_rerun/PvP01_to_PAM.tuned.cleaned.chain', score_min=1000)),
    ('3_wfmash_n1',
     parse_paf_blocks('projection/A1_wfmash/PvP01_vs_GCA_949152365.1.paf', min_block=1000)),
    ('4_wfmash_pggblike',
     parse_paf_blocks('projection/A1_wfmash/synteny_rerun/PvP01_vs_PAM.pggblike.paf', min_block=1000)),
    ('5_graph_blocks',
     parse_graph_blocks('writeup/synteny_graph/PvP01_to_PAM.blocks.tsv', min_block=1000)),
]

print(f'tracks for {PVP01_CHR1} (PvP01 chr 1):\n')
for name, blocks in sources:
    n = write_bed(blocks, OUT / f'{name}.bed')
    print(f'  {name}: {n} blocks → {OUT / (name + ".bed")}')
