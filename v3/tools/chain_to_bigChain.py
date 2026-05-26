#!/usr/bin/env python3
"""Convert UCSC chain (gz or plain) to bigChain.bed + bigLink.bed.

Walks chain blocks to emit the link rows in chain coordinates,
mirroring `hgLoadChain -noBin` output structure.
"""
import gzip
import sys
from pathlib import Path


def open_text(p):
    p = Path(p)
    if p.suffix == ".gz":
        return gzip.open(p, "rt")
    return open(p, "r")


def convert(chain_path, big_bed, big_link):
    with open_text(chain_path) as fh, open(big_bed, "w") as out_bed, open(big_link, "w") as out_link:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("chain "):
                parts = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                score = int(parts[1])
                tName = parts[2]
                tSize = int(parts[3])
                tStrand = parts[4]
                tStart = int(parts[5])
                tEnd = int(parts[6])
                qName = parts[7]
                qSize = int(parts[8])
                qStrand = parts[9]
                qStart = int(parts[10])
                qEnd = int(parts[11])
                chain_id = parts[12]
                if tStrand != "+":
                    sys.stderr.write(f"Skipping chain {chain_id}: tStrand={tStrand} not supported\n")
                    skip_chain = True
                    continue
                skip_chain = False
                # bigChain row: tName tStart tEnd id 1000 qStrand tSize qName qSize qStart qEnd chainScore
                out_bed.write("\t".join(str(x) for x in [
                    tName, tStart, tEnd, chain_id, 1000, qStrand,
                    tSize, qName, qSize, qStart, qEnd, score
                ]) + "\n")
                t_cur = tStart
                q_cur = qStart
            else:
                if skip_chain:
                    continue
                parts = line.split()
                if len(parts) == 3:
                    size, dt, dq = int(parts[0]), int(parts[1]), int(parts[2])
                elif len(parts) == 1:
                    size, dt, dq = int(parts[0]), 0, 0
                else:
                    continue
                # bigLink row: tName tStart tEnd chain_id qStart
                out_link.write("\t".join(str(x) for x in [
                    tName, t_cur, t_cur + size, chain_id, q_cur
                ]) + "\n")
                t_cur += size + dt
                q_cur += size + dq


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("usage: chain_to_bigChain.py in.chain[.gz] out.bigChain.bed out.bigLink.bed\n")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2], sys.argv[3])
