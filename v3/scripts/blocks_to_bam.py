#!/usr/bin/env python3
"""
Convert each synteny source to BAM aligned to PvP01 chr 1 (LT635612.2):

  PAF  → SAM (using cg:Z: CIGAR tag directly)
  chain → SAM (CIGAR built from chain block triples [size, dt, dq])
  graph blocks → SAM (CIGAR = single M op since no per-base alignment available)

Each block becomes one SAM record. Output BAMs are sorted + indexed.
"""
import re
import subprocess
from pathlib import Path

CHR = 'LT635612.2'
CHR_LEN = 1_021_664
OUT = Path('writeup/browser/tracks_bam')
OUT.mkdir(parents=True, exist_ok=True)

SAM_HEADER = f'@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{CHR}\tLN:{CHR_LEN}\n'


def cigar_from_paf(cg_z):
    """PAF cg:Z: encoding uses M/I/D directly (or =/X for match/mismatch)."""
    return cg_z


def write_sam_paf(paf_path, sam_path, pg='wfmash'):
    """One record per PAF line, target=PvP01_CHR1."""
    with open(paf_path) as fin, open(sam_path, 'w') as fout:
        fout.write(SAM_HEADER)
        fout.write(f'@PG\tID:{pg}\tPN:{pg}\n')
        n = 0
        for ln in fin:
            f = ln.rstrip('\n').split('\t')
            if len(f) < 12 or f[5] != CHR:
                continue
            qName = f[0]; qStart = int(f[2]); qEnd = int(f[3])
            strand = f[4]; tStart = int(f[7]); tEnd = int(f[8])
            cigar = None
            for tag in f[12:]:
                if tag.startswith('cg:Z:'):
                    cigar = tag[5:]; break
            if not cigar:
                continue
            flag = 16 if strand == '-' else 0
            fout.write(
                f'{qName}_{qStart}_{qEnd}\t{flag}\t{CHR}\t{tStart+1}\t60\t{cigar}\t*\t0\t0\t*\t*\n'
            )
            n += 1
    return n


def cigar_from_chain_blocks(chain_lines, t_strand, q_strand):
    """Chain body: sequence of lines '{size}\\t{dt}\\t{dq}' (last line just '{size}').
    Build CIGAR for the chain as M for match block, D for target-only gap, I for query-only gap."""
    cigar = []
    for raw in chain_lines:
        parts = raw.strip().split('\t')
        if not parts or not parts[0]:
            continue
        size = int(parts[0])
        if size > 0:
            cigar.append(f'{size}M')
        if len(parts) == 3:
            dt = int(parts[1]); dq = int(parts[2])
            # dt = bases in target only (=> D from query perspective)
            # dq = bases in query only (=> I from query perspective)
            if dt > 0: cigar.append(f'{dt}D')
            if dq > 0: cigar.append(f'{dq}I')
    return ''.join(cigar)


def write_sam_chain(chain_path, sam_path, pg='kegalign_chain'):
    """One record per chain (in target=CHR). CIGAR built from chain body."""
    with open(chain_path) as fin, open(sam_path, 'w') as fout:
        fout.write(SAM_HEADER)
        fout.write(f'@PG\tID:{pg}\tPN:{pg}\n')
        in_chain = False
        cur_header = None
        cur_body = []
        n = 0
        def flush():
            nonlocal n, cur_header, cur_body
            if cur_header is None:
                return
            score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, cid = cur_header
            if tName == CHR:
                cigar = cigar_from_chain_blocks(cur_body, tStrand, qStrand)
                flag = 16 if qStrand == '-' else 0
                qname_full = f'{qName}_{qStart}_{qEnd}'
                fout.write(f'{qname_full}\t{flag}\t{CHR}\t{tStart+1}\t60\t{cigar}\t*\t0\t0\t*\t*\n')
                n += 1
            cur_header = None
            cur_body = []
        for ln in fin:
            ln = ln.rstrip('\n')
            if ln.startswith('chain '):
                flush()
                f = ln.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                cur_header = (int(f[1]), f[2], int(f[3]), f[4], int(f[5]), int(f[6]),
                              f[7], int(f[8]), f[9], int(f[10]), int(f[11]),
                              f[12] if len(f) > 12 else '0')
                cur_body = []
            elif ln.strip() and cur_header is not None:
                cur_body.append(ln)
        flush()
    return n


def write_sam_graph(blocks_tsv, sam_path, pg='graph_blocks'):
    """One record per graph block. No per-base alignment → CIGAR is one big M."""
    with open(blocks_tsv) as fin, open(sam_path, 'w') as fout:
        fout.write(SAM_HEADER)
        fout.write(f'@PG\tID:{pg}\tPN:{pg}\n')
        next(fin)  # header
        n = 0
        for ln in fin:
            f = ln.rstrip('\n').split('\t')
            tName = f[0]
            if not tName.endswith(f'#{CHR}#0'):
                continue
            tStart, tEnd = int(f[1]), int(f[2])
            qName = f[3]
            qStart, qEnd = int(f[4]), int(f[5])
            strand = f[6]
            flag = 16 if strand == '-' else 0
            cigar = f'{tEnd - tStart}M'
            qname_full = f'{qName.split("#")[-2]}_{qStart}_{qEnd}'
            fout.write(f'{qname_full}\t{flag}\t{CHR}\t{tStart+1}\t60\t{cigar}\t*\t0\t0\t*\t*\n')
            n += 1
    return n


def sort_and_index(sam_path, bam_path):
    """samtools sort + index, all via docker."""
    cmd = (f"docker run --rm -v {Path('.').absolute()}:/v3 "
           "quay.io/biocontainers/samtools:1.20--h50ea8bc_0 bash -c "
           f"'samtools sort -O bam /v3/{sam_path} > /v3/{bam_path} && "
           f"samtools index /v3/{bam_path}'")
    subprocess.run(cmd, shell=True, check=True, capture_output=True)


sources = [
    ('1_kegalign_default',  'chain', 'work/01_chains/GCA_900093555.2.GCA_949152365.1.cleaned.chain'),
    ('2_kegalign_tuned',    'chain', 'projection/A1_wfmash/synteny_rerun/PvP01_to_PAM.tuned.cleaned.chain'),
    ('3_wfmash_n1',         'paf',   'projection/A1_wfmash/PvP01_vs_GCA_949152365.1.paf'),
    ('4_wfmash_pggblike',   'paf',   'projection/A1_wfmash/synteny_rerun/PvP01_vs_PAM.pggblike.paf'),
    ('5_graph_blocks',      'graph', 'writeup/synteny_graph/PvP01_to_PAM.blocks.tsv'),
]

for name, kind, path in sources:
    sam = OUT / f'{name}.sam'
    bam = OUT / f'{name}.bam'
    if kind == 'paf':
        n = write_sam_paf(path, sam, pg=name)
    elif kind == 'chain':
        n = write_sam_chain(path, sam, pg=name)
    elif kind == 'graph':
        n = write_sam_graph(path, sam, pg=name)
    print(f'  {name}: {n} records → {sam}')
    sort_and_index(sam, bam)
    sam.unlink()
    print(f'    → {bam} (+ .bai)')
