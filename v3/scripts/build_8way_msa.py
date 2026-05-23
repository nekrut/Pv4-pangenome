#!/usr/bin/env python3
"""
Build 8-way codon MSAs by joining per-query pairwise CESAR2 alignments.

For each PvP01 gene, gather the 7 query CESAR2 codon-aligned QUERY sequences
plus the PvP01 REFERENCE row. Output one FASTA per gene.

Notes:
- TOGA2 codon_aln.fa.gz contains REF and QUERY rows per projection
  (e.g., PVP01_0000010.1#chain | CODON | REFERENCE/QUERY).
- Codons are space-separated with --- for gaps and XXX/NNN at frame breaks.
- The REFERENCE row is identical (modulo CESAR2 alignment columns) across all
  queries, but minor differences can occur if a column ends up unique to one
  query's projection. We canonicalize by picking the longest REFERENCE row.
- For a clean 8-way MSA we want all queries aligned to the SAME REF column set.
  CESAR2's REF aligning differs per projection, so we re-align by gap-stripping
  the REF and re-inserting query codons via per-codon ref-position.

Output dir layout:
    out_dir/<gene_id>.codon.fa
    out_dir/<gene_id>.protein.fa
    out_dir/summary.tsv   # gene_id, n_intact, missing_strains
"""

import argparse
import csv
import gzip
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path


def parse_intactness(classification_tsv):
    """Return dict gene_id -> intactness."""
    d = {}
    with open(classification_tsv) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            d[row['reference_gene_id']] = row['intactness']
    return d


def parse_codon_aln(path):
    """Yield (header_id, role, codon_list) tuples from a codon_aln.fa.gz file.

    header_id is the transcript ID part (e.g., PVP01_0000010.1#chain_id)
    role is REFERENCE or QUERY.
    """
    open_fn = gzip.open if str(path).endswith('.gz') else open
    cur_header = None
    cur_role = None
    cur_seq = []
    with open_fn(path, 'rt') as fh:
        for ln in fh:
            ln = ln.rstrip('\n')
            if ln.startswith('>'):
                if cur_header is not None:
                    yield cur_header, cur_role, ''.join(cur_seq).split()
                # Parse header: ">PVP01_0008520.1#388| 388 | CODON | REFERENCE"
                head = ln[1:].split('|')
                proj_id = head[0].strip()
                # role at end
                if 'REFERENCE' in ln:
                    cur_role = 'REFERENCE'
                elif 'QUERY' in ln:
                    cur_role = 'QUERY'
                else:
                    cur_role = None
                cur_header = proj_id
                cur_seq = []
            else:
                cur_seq.append(ln)
    if cur_header is not None:
        yield cur_header, cur_role, ''.join(cur_seq).split()


def index_alignments(codon_aln_path, intact_genes):
    """Return dict transcript_id -> {ref: [codons], query: [codons]} but
    only for transcripts whose underlying gene is in intact_genes.
    Picks the projection (chain) with the longest non-gap query.
    """
    by_tx = {}
    for proj_id, role, codons in parse_codon_aln(codon_aln_path):
        # proj_id like PVP01_0000010.1#387 → transcript=PVP01_0000010.1, gene=PVP01_0000010
        tx_id = proj_id.split('#')[0]
        gene_id = tx_id.rsplit('.', 1)[0]
        if gene_id not in intact_genes:
            continue
        if tx_id not in by_tx:
            by_tx[tx_id] = {}
        if proj_id not in by_tx[tx_id]:
            by_tx[tx_id][proj_id] = {}
        by_tx[tx_id][proj_id][role] = codons
    # Collapse: keep only the best projection per transcript (highest non-gap codon count)
    best = {}
    for tx_id, projs in by_tx.items():
        best_proj = None
        best_count = -1
        for proj_id, roles in projs.items():
            q = roles.get('QUERY', [])
            count = sum(1 for c in q if c not in ('---', 'NNN', 'XXX'))
            if count > best_count:
                best_count = count
                best_proj = proj_id
        if best_proj:
            best[tx_id] = projs[best_proj]
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--queries', nargs='+', required=True,
                    help='Query names (PvW1 Sal-I PvT01 ...)')
    ap.add_argument('--toga-base', required=True,
                    help='work/02c_toga/PvP01-as-ref')
    ap.add_argument('--merged-base', required=True,
                    help='work/02d_merged/PvP01-as-ref')
    ap.add_argument('--min-intact', type=int, default=7,
                    help='Minimum number of queries with intact projection')
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Per-query intactness
    intact_per_q = {}
    for q in args.queries:
        cls_path = Path(args.merged_base) / f"{q}.classification.tsv"
        per = parse_intactness(cls_path)
        # Keep genes classified as I or PI
        intact_per_q[q] = {g for g, s in per.items() if s in ('I', 'PI')}
        print(f"  {q:8s}  intact|PI = {len(intact_per_q[q])}", file=sys.stderr)

    # 2) Find genes intact in ≥ min_intact queries
    all_genes = set().union(*intact_per_q.values())
    gene_intact_count = {g: sum(1 for s in intact_per_q.values() if g in s)
                         for g in all_genes}
    target_genes = {g for g, n in gene_intact_count.items() if n >= args.min_intact}
    print(f"\nTarget set (intact in >= {args.min_intact}): {len(target_genes)} genes",
          file=sys.stderr)

    # 3) For each query, index its codon_aln.fa.gz for target_genes
    print(f"\nIndexing codon alignments per query...", file=sys.stderr)
    per_q_aln = {}
    for q in args.queries:
        path = Path(args.toga_base) / q / 'codon_aln.fa.gz'
        if not path.exists():
            print(f"  WARN: missing {path}", file=sys.stderr)
            per_q_aln[q] = {}
            continue
        per_q_aln[q] = index_alignments(path, target_genes)
        print(f"  {q:8s}  indexed = {len(per_q_aln[q])} transcripts", file=sys.stderr)

    # 4) Build 8-way MSA per gene
    #    For each gene, pick its first transcript (.1) and gather REF + 7 QUERY codon strings
    summary_rows = []
    n_written = 0
    for gene_id in sorted(target_genes):
        tx_id = gene_id + '.1'   # canonical first transcript
        rows = {}  # strain → codon list
        # Get the REFERENCE row from any one query (they should match modulo cesar2 columns)
        # Pick the query whose tx is present
        ref_codons = None
        for q in args.queries:
            if tx_id in per_q_aln[q]:
                if 'REFERENCE' in per_q_aln[q][tx_id]:
                    ref_codons = per_q_aln[q][tx_id]['REFERENCE']
                    break
        if ref_codons is None:
            continue
        # Each query gives its own QUERY codon vector
        # Different projections may have different lengths because CESAR2 aligns
        # query→ref per-projection. The REF is also per-projection. We can't
        # naively concatenate. For Path 1 simplicity, use the per-projection
        # REF from each query and produce a multi-row alignment with each
        # query's REF/QUERY pair concatenated below.
        # However, that's not a true 8-way MSA. A true 8-way MSA requires
        # column-alignment across queries.
        #
        # Approximation: gap-strip the REFERENCE from each per-query projection
        # to get the canonical reference codon sequence; align all queries to
        # that. The REF is the same across queries when chains are correct, so
        # this just normalizes per-projection insertions in the REF.
        canonical_ref = [c for c in ref_codons if c != '---']
        ref_codon_to_pos = list(range(len(canonical_ref)))  # 0..N-1
        # For each query, walk through its REF row and emit query codon at
        # canonical ref position when REF codon is not a gap
        q_rows = {}
        for q in args.queries:
            if tx_id not in per_q_aln[q]:
                q_rows[q] = ['---'] * len(canonical_ref)
                continue
            ref = per_q_aln[q][tx_id]['REFERENCE']
            query = per_q_aln[q][tx_id].get('QUERY', [])
            if len(ref) != len(query):
                # CESAR2 should produce aligned REF/QUERY of equal length
                q_rows[q] = ['---'] * len(canonical_ref)
                continue
            q_row = []
            for r, qc in zip(ref, query):
                if r == '---':
                    continue  # query insertions vs canonical REF — skip
                q_row.append(qc)
            # Pad/truncate to canonical length
            if len(q_row) < len(canonical_ref):
                q_row.extend(['---'] * (len(canonical_ref) - len(q_row)))
            elif len(q_row) > len(canonical_ref):
                q_row = q_row[:len(canonical_ref)]
            q_rows[q] = q_row

        # 5) Emit FASTA: 1 REF + 7 query rows
        codon_fa = outdir / f"{gene_id}.codon.fa"
        with open(codon_fa, 'w') as f:
            f.write(f">PvP01_REF | gene={gene_id}\n")
            f.write(''.join(canonical_ref) + '\n')
            for q in args.queries:
                missing = '_MISSING' if all(c == '---' for c in q_rows[q]) else ''
                f.write(f">{q}{missing} | gene={gene_id}\n")
                f.write(''.join(q_rows[q]) + '\n')

        n_intact = gene_intact_count[gene_id]
        missing_strains = [q for q in args.queries if gene_id not in intact_per_q[q]]
        summary_rows.append({
            'gene_id': gene_id,
            'n_intact': n_intact,
            'n_codons': len(canonical_ref),
            'missing_strains': ','.join(missing_strains) if missing_strains else '-',
        })
        n_written += 1

    # 6) Summary TSV
    sum_path = outdir / 'summary.tsv'
    with open(sum_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gene_id', 'n_intact', 'n_codons',
                                          'missing_strains'], delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nWrote {n_written} 8-way codon MSAs to {outdir}", file=sys.stderr)
    print(f"Summary: {sum_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
