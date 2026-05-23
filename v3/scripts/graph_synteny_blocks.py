#!/usr/bin/env python3
"""
Extract pairwise synteny blocks directly from a PGGB GFA.

For each (target, query) path pair, find contiguous runs of shared nodes.
This is what odgi untangle does, computed from the GFA directly because
odgi 0.8.2 and 0.9.4 untangle both crash on the v2 graph (SDSL int_vector
assertion).

Output: chain-like BEDPE per pair: tName tStart tEnd qName qStart qEnd strand score
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


def parse_gfa_paths(gfa_path, prefix_target='GCA_900093555.2', prefixes_query=None):
    """Single-pass parse of GFA:
       node_lens: dict node_id (int) → length
       paths_target: dict path_name → list[(node_id, '+'/'-', cum_start_bp)]
       paths_query: dict path_name → list[(node_id, '+'/'-', cum_start_bp)]
    """
    node_lens = {}
    paths_target = {}
    paths_query = {}
    if prefixes_query is None:
        prefixes_query = []
    print(f"  parsing GFA: {gfa_path}", file=sys.stderr)
    with open(gfa_path) as fh:
        for ln in fh:
            if ln.startswith('S\t'):
                f = ln.split('\t', 3)
                node_id = int(f[1])
                seq_len = len(f[2])
                node_lens[node_id] = seq_len
            elif ln.startswith('P\t'):
                f = ln.rstrip('\n').split('\t')
                if len(f) < 3:
                    continue
                pname = f[1]
                # Match prefix
                is_target = pname.startswith(prefix_target + '#')
                is_query = any(pname.startswith(pq + '#') for pq in prefixes_query)
                if not (is_target or is_query):
                    continue
                # Parse node sequence: '1+,2-,3+,...'
                steps = []
                cum = 0
                for tok in f[2].split(','):
                    if not tok:
                        continue
                    strand = tok[-1]
                    nid = int(tok[:-1])
                    nl = node_lens.get(nid, 0)
                    steps.append((nid, strand, cum))
                    cum += nl
                if is_target:
                    paths_target[pname] = steps
                if is_query:
                    paths_query[pname] = steps
    print(f"  parsed {len(node_lens)} nodes, "
          f"{len(paths_target)} target paths, "
          f"{len(paths_query)} query paths", file=sys.stderr)
    return node_lens, paths_target, paths_query


def derive_blocks(target_path_steps, query_path_steps_dict, node_lens, min_block_bp=1000):
    """For one target path, find blocks of contiguous shared nodes with each
    query path.

    Returns list of dicts: {qPath, tStart, tEnd, qStart, qEnd, strand, score}
    """
    # Build node → (q_path_name, cum_position, strand) lookup. A node may
    # appear multiple times in a path (loops/repeats) — we store all instances.
    node_to_q = defaultdict(list)
    for qname, steps in query_path_steps_dict.items():
        for nid, strand, cum in steps:
            node_to_q[nid].append((qname, cum, strand))

    blocks = []
    # State: current block being extended
    current = {}  # (qname, q_strand) → {tStart, tEnd, qStart, qEnd, last_q_cum, last_t_cum, score}
    for nid, t_strand, t_cum in target_path_steps:
        nl = node_lens.get(nid, 0)
        if nl == 0:
            continue
        if nid not in node_to_q:
            # Flush all current blocks if no match — they'll be re-evaluated below
            continue
        # For each query instance of this node, extend or start a block
        matched_keys = set()
        for qname, q_cum, q_strand in node_to_q[nid]:
            # Relative strand
            rel = '+' if t_strand == q_strand else '-'
            key = (qname, rel)
            matched_keys.add(key)
            if key in current:
                # Continuation candidate: q_cum must be near last_q_cum (within ~node-length tolerance)
                blk = current[key]
                expected = blk['last_q_cum_end']
                if rel == '+':
                    delta = q_cum - expected
                else:
                    # For - strand the query path is being traversed in reverse
                    delta = blk['last_q_cum_start'] - (q_cum + nl)
                # Allow gap up to a few node-lengths (graph bubbles)
                if abs(delta) < 5000:
                    blk['tEnd'] = t_cum + nl
                    blk['score'] += nl
                    if rel == '+':
                        blk['qEnd'] = max(blk['qEnd'], q_cum + nl)
                        blk['last_q_cum_end'] = q_cum + nl
                    else:
                        blk['qStart'] = min(blk['qStart'], q_cum)
                        blk['last_q_cum_start'] = q_cum
                    continue
                # Gap too large — close previous, start new
                if blk['score'] >= min_block_bp:
                    blocks.append({**blk, 'qName': qname})
                del current[key]
            # Start new block
            current[key] = {
                'tStart': t_cum,
                'tEnd': t_cum + nl,
                'qStart': q_cum,
                'qEnd': q_cum + nl,
                'strand': rel,
                'score': nl,
                'last_q_cum_end': q_cum + nl,
                'last_q_cum_start': q_cum,
            }
        # Close blocks for query paths that didn't match this node
        for key in list(current):
            if key not in matched_keys:
                blk = current[key]
                if blk['score'] >= min_block_bp:
                    blocks.append({**blk, 'qName': key[0]})
                del current[key]
    # Final flush
    for key, blk in current.items():
        if blk['score'] >= min_block_bp:
            blocks.append({**blk, 'qName': key[0]})
    return blocks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gfa', required=True)
    ap.add_argument('--target-prefix', default='GCA_900093555.2')
    ap.add_argument('--query-prefix', action='append', required=True,
                    help='Repeatable, e.g. --query-prefix GCA_914969965.1')
    ap.add_argument('--min-block-bp', type=int, default=1000)
    ap.add_argument('--out-tsv', required=True)
    args = ap.parse_args()

    node_lens, paths_t, paths_q = parse_gfa_paths(
        args.gfa, args.target_prefix, args.query_prefix
    )

    n_blocks_total = 0
    with open(args.out_tsv, 'w') as fout:
        fout.write('tName\ttStart\ttEnd\tqName\tqStart\tqEnd\tstrand\tscore\n')
        for tname, steps in paths_t.items():
            print(f"  target {tname} (steps={len(steps)})", file=sys.stderr)
            blocks = derive_blocks(steps, paths_q, node_lens, args.min_block_bp)
            for b in blocks:
                fout.write(f"{tname}\t{b['tStart']}\t{b['tEnd']}\t{b['qName']}\t"
                           f"{b['qStart']}\t{b['qEnd']}\t{b['strand']}\t{b['score']}\n")
            n_blocks_total += len(blocks)
            print(f"    blocks={len(blocks)}", file=sys.stderr)
    print(f"\nDone: {n_blocks_total} blocks → {args.out_tsv}", file=sys.stderr)


if __name__ == '__main__':
    main()
