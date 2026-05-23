#!/usr/bin/env python3
"""
Convert PAF (with cg:Z:CIGAR tags) → UCSC chain format.

Chain format spec: https://genome.ucsc.edu/goldenPath/help/chain.html

PAF columns: qName qLen qStart qEnd strand tName tLen tStart tEnd matches blockLen mapq [tags]
Tag cg:Z: holds the CIGAR (M/I/D/=/X/N/S/H/P).
We emit:
  chain score tName tSize + tStart tEnd qName qSize qStrand qStart qEnd id
  size [dt dq]
  ...
"""
import sys, re

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')

def parse_paf(line):
    f = line.rstrip('\n').split('\t')
    pa = {
        'qName': f[0], 'qLen': int(f[1]), 'qStart': int(f[2]), 'qEnd': int(f[3]),
        'strand': f[4],
        'tName': f[5], 'tLen': int(f[6]), 'tStart': int(f[7]), 'tEnd': int(f[8]),
        'matches': int(f[9]), 'block': int(f[10]),
    }
    cigar = None
    for tag in f[12:]:
        if tag.startswith('cg:Z:'):
            cigar = tag[5:]
            break
    pa['cigar'] = cigar
    return pa

def emit_chain(pa, chain_id):
    ops = CIGAR_RE.findall(pa['cigar']) if pa['cigar'] else []
    # accumulate matching-block sizes interleaved with insert/delete gaps
    blocks = []        # list of (size, dt, dq); last entry: (size,)
    cur_match = 0
    pending_dt = 0
    pending_dq = 0
    for n_str, op in ops:
        n = int(n_str)
        if op in ('=', 'M', 'X'):
            if pending_dt or pending_dq:
                # flush previous matching block + its trailing gap
                blocks.append((cur_match, pending_dt, pending_dq))
                cur_match = 0
                pending_dt = pending_dq = 0
            cur_match += n
        elif op == 'I':
            # insertion in query relative to target → dq (query gap) on chain reads as "query insert"
            pending_dq += n
        elif op == 'D':
            pending_dt += n
        elif op in ('N', 'S', 'H', 'P'):
            pass
    blocks.append((cur_match,))
    score = pa['matches']
    qStart, qEnd = pa['qStart'], pa['qEnd']
    if pa['strand'] == '-':
        qStart, qEnd = pa['qLen'] - pa['qEnd'], pa['qLen'] - pa['qStart']
    out = []
    out.append(f"chain {score} {pa['tName']} {pa['tLen']} + {pa['tStart']} {pa['tEnd']} "
               f"{pa['qName']} {pa['qLen']} {pa['strand']} {qStart} {qEnd} {chain_id}")
    for blk in blocks[:-1]:
        out.append(f"{blk[0]}\t{blk[1]}\t{blk[2]}")
    out.append(f"{blocks[-1][0]}")
    out.append("")
    return "\n".join(out)

def main():
    if len(sys.argv) < 2:
        print("usage: paf2chain.py input.paf > output.chain", file=sys.stderr)
        sys.exit(1)
    n_in = n_out = 0
    next_id = 1
    with open(sys.argv[1]) as fh:
        for line in fh:
            if not line.strip(): continue
            n_in += 1
            pa = parse_paf(line)
            if pa['cigar'] is None:
                continue
            print(emit_chain(pa, next_id))
            next_id += 1
            n_out += 1
    print(f"paf2chain: {n_out}/{n_in} rows had CIGAR and were converted", file=sys.stderr)

if __name__ == '__main__':
    main()
