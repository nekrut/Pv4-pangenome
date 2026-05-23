#!/usr/bin/env python3
"""
Phase C.4: merge Liftoff-clean + TOGA2/CESAR2 outputs into unified annotation
per query.

Rules:
  - Genes in liftoff_clean.gff3 → source=liftoff, intactness=I.
  - Genes in TOGA2 query_annotation.bed → source=cesar2, intactness from
    loss_summary.tsv (I/PI/UL/L/PG/M).
  - Genes flagged by triage but absent from TOGA2 output → source=cesar2,
    intactness=M (or omit, with warning).

Outputs per query:
  - {Q}.annotation.gff3 — unified GFF (Liftoff clean entries + CESAR2 entries)
  - {Q}.classification.tsv — flat per-gene table:
      reference_gene_id, query_gene_id, source, intactness, query_chrom,
      query_start, query_end, query_strand, orthology_class
"""

import argparse
import csv
import sys
from pathlib import Path
from collections import defaultdict


def parse_gff_attributes(s):
    d = {}
    for kv in s.strip().rstrip(';').split(';'):
        kv = kv.strip()
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k.strip()] = v.strip()
    return d


def load_liftoff_clean(gff_path):
    """Return dict reference_gene_id -> list of GFF lines (gene + children).
    Assumes GFF order: each gene's child features (mRNA, exon, CDS) immediately
    follow the gene line until the next gene.
    """
    genes = {}      # ref_id → [lines]
    current_ref_id = None
    if not Path(gff_path).exists():
        return genes
    GENE_TYPES = {'gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'}
    with open(gff_path) as f:
        for ln in f:
            if ln.startswith('#') or not ln.strip():
                continue
            fields = ln.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2]
            if ftype in GENE_TYPES:
                attrs = parse_gff_attributes(fields[8])
                gid = attrs.get('ID', '')
                ref_id = gid
                if '_' in gid:
                    parts = gid.rsplit('_', 1)
                    if len(parts[1]) <= 2 and parts[1].isdigit() and not parts[0].endswith('_'):
                        ref_id = parts[0]
                current_ref_id = ref_id
                genes.setdefault(ref_id, []).append(ln)
            elif current_ref_id is not None:
                genes[current_ref_id].append(ln)
    return genes


def load_toga2_loss_summary(path):
    """Return dict projection_id -> status (e.g., 'PVP01_0000020.1#15' -> 'I')."""
    if not Path(path).exists():
        return {}
    status = {}
    with open(path) as f:
        for ln in f:
            ln = ln.rstrip('\n')
            if not ln or ln.startswith('level'):
                continue
            fields = ln.split('\t')
            if len(fields) >= 3 and fields[0] == 'PROJECTION':
                status[fields[1]] = fields[2]
    return status


def load_toga2_orthology(path):
    """Return dict reference_gene_id -> {class, query_gene_id, query_transcript_id}."""
    if not Path(path).exists():
        return {}
    out = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter='\t')
        for row in r:
            t_gene = row.get('t_gene')
            if not t_gene:
                continue
            out.setdefault(t_gene, []).append({
                'q_gene': row.get('q_gene'),
                'q_tx': row.get('q_transcript'),
                'class': row.get('orthology_class'),
                't_tx': row.get('t_transcript'),
            })
    return out


def load_toga2_query_bed(path):
    """Return dict query_gene_id -> BED12 line (raw, query coords)."""
    if not Path(path).exists():
        return {}
    out = {}
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith('#'):
                continue
            fields = ln.rstrip('\n').split('\t')
            if len(fields) >= 4:
                out[fields[3]] = ln
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--query', required=True, help='Query name (PvW1, Sal-I, etc.)')
    ap.add_argument('--triage-dir', required=True, help='work/02b_triage/PvP01-as-ref/<Q>/')
    ap.add_argument('--toga-dir', required=True, help='work/02c_toga/PvP01-as-ref/<Q>/')
    ap.add_argument('--out-dir', required=True, help='work/02d_merged/PvP01-as-ref/')
    ap.add_argument('--ref-bed', help='inputs/annotations/PvP01.bed (gene-level)', required=True)
    args = ap.parse_args()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Load Liftoff-clean GFF
    liftoff_clean = load_liftoff_clean(f"{args.triage_dir}/liftoff_clean.gff3")

    # 2) Load TOGA2 outputs
    loss = load_toga2_loss_summary(f"{args.toga_dir}/loss_summary.tsv")
    ortho = load_toga2_orthology(f"{args.toga_dir}/orthology_classification.tsv")
    query_bed = load_toga2_query_bed(f"{args.toga_dir}/query_annotation.bed")
    query_genes_bed = load_toga2_query_bed(f"{args.toga_dir}/query_genes.bed")

    # 3) Reference BED gene set
    ref_genes = set()
    with open(args.ref_bed) as f:
        for ln in f:
            fields = ln.rstrip('\n').split('\t')
            if len(fields) >= 4:
                ref_genes.add(fields[3])

    # 4) Classification table rows
    rows = []
    seen_ref = set()

    # 4a) Liftoff clean
    for ref_id, lines in liftoff_clean.items():
        seen_ref.add(ref_id)
        # First line is gene; get chrom/start/end/strand
        gene_line = lines[0]
        fields = gene_line.rstrip('\n').split('\t')
        if len(fields) < 9:
            continue
        chrom, start, end, strand = fields[0], fields[3], fields[4], fields[6]
        attrs = parse_gff_attributes(fields[8])
        rows.append({
            'reference_gene_id': ref_id,
            'query_gene_id': attrs.get('ID', ''),
            'source': 'liftoff',
            'intactness': 'I',
            'query_chrom': chrom,
            'query_start': start,
            'query_end': end,
            'query_strand': strand,
            'orthology_class': 'liftoff_clean',
        })

    # 4b) TOGA2 / CESAR2
    # Walk orthology_classification.tsv for genes with q_gene != "None"
    for ref_gene_id, projections in ortho.items():
        for p in projections:
            q_gene = p['q_gene']
            if q_gene in (None, '', 'None'):
                # Lost in query, no projection
                rows.append({
                    'reference_gene_id': ref_gene_id,
                    'query_gene_id': '',
                    'source': 'cesar2',
                    'intactness': 'L',
                    'query_chrom': '',
                    'query_start': '',
                    'query_end': '',
                    'query_strand': '',
                    'orthology_class': p.get('class', 'one2zero'),
                })
                seen_ref.add(ref_gene_id)
                continue
            # Look up status — use t_transcript as proxy (e.g., PVP01_0000020.1#15)
            # Loss summary keys are projection IDs (PVP01_0000020.1#chainID)
            tx = p.get('t_tx', '')
            status_keys = [k for k in loss if k.startswith(tx + '#')]
            status = '?'
            if status_keys:
                # Pick the first matching status; usually one per transcript
                status = loss[status_keys[0]]
            # Get query BED coords
            bed_line = query_bed.get(q_gene, '') or query_genes_bed.get(q_gene, '')
            if bed_line:
                fields = bed_line.rstrip('\n').split('\t')
                rows.append({
                    'reference_gene_id': ref_gene_id,
                    'query_gene_id': q_gene,
                    'source': 'cesar2',
                    'intactness': status,
                    'query_chrom': fields[0] if len(fields) > 0 else '',
                    'query_start': fields[1] if len(fields) > 1 else '',
                    'query_end': fields[2] if len(fields) > 2 else '',
                    'query_strand': fields[5] if len(fields) > 5 else '',
                    'orthology_class': p.get('class', ''),
                })
            else:
                rows.append({
                    'reference_gene_id': ref_gene_id,
                    'query_gene_id': q_gene,
                    'source': 'cesar2',
                    'intactness': status,
                    'query_chrom': '',
                    'query_start': '',
                    'query_end': '',
                    'query_strand': '',
                    'orthology_class': p.get('class', ''),
                })
            seen_ref.add(ref_gene_id)

    # 4c) Reference genes neither in liftoff-clean nor TOGA2 → mark as missing/unprojected
    missing = ref_genes - seen_ref
    for ref_id in missing:
        rows.append({
            'reference_gene_id': ref_id,
            'query_gene_id': '',
            'source': 'none',
            'intactness': 'M',
            'query_chrom': '',
            'query_start': '',
            'query_end': '',
            'query_strand': '',
            'orthology_class': 'unprojected',
        })

    # 5) Write classification.tsv
    out_cls = outdir / f"{args.query}.classification.tsv"
    fieldnames = ['reference_gene_id', 'query_gene_id', 'source', 'intactness',
                  'query_chrom', 'query_start', 'query_end', 'query_strand',
                  'orthology_class']
    with open(out_cls, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for row in rows:
            w.writerow(row)

    # 6) Write merged GFF (liftoff GFF entries + TOGA2 BED12 converted to GFF lines)
    out_gff = outdir / f"{args.query}.annotation.gff3"
    n_lo = 0
    n_cs = 0
    with open(out_gff, 'w') as f:
        f.write('##gff-version 3\n')
        f.write(f'# Phase C.4 merged annotation for {args.query}\n')
        f.write('# source=liftoff means triage-clean Liftoff projection; intactness=I assumed\n')
        f.write('# source=cesar2 means TOGA2/CESAR2 fallback; intactness from TOGA2 loss_summary\n')
        # Liftoff clean GFF entries with source tagged
        for ref_id, lines in liftoff_clean.items():
            for ln in lines:
                fields = ln.rstrip('\n').split('\t')
                if len(fields) >= 9:
                    # Append source=liftoff to attrs
                    fields[8] = fields[8].rstrip(';') + ';source=liftoff;intactness=I'
                    f.write('\t'.join(fields) + '\n')
                    if fields[2] in ('gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'):
                        n_lo += 1
        # Build q_gene → (ref_id, orth_class, t_tx) lookup once
        q_to_ref = {}
        for rg, projs in ortho.items():
            for p in projs:
                qg = p.get('q_gene')
                if qg and qg != 'None':
                    q_to_ref[qg] = (rg, p.get('class', ''), p.get('t_tx', ''))
        # Build t_tx prefix → status lookup
        tx_to_status = {}
        for k, v in loss.items():
            if '#' in k:
                tx_prefix = k.rsplit('#', 1)[0]
                if tx_prefix not in tx_to_status:
                    tx_to_status[tx_prefix] = v
        # TOGA2 query_genes.bed lines as GFF
        for q_gene, bed_ln in query_genes_bed.items():
            fields = bed_ln.rstrip('\n').split('\t')
            if len(fields) < 6:
                continue
            chrom, start, end, name, score, strand = fields[:6]
            ref_gene_id, orth_class, t_tx = q_to_ref.get(q_gene, ('', '', ''))
            intactness = tx_to_status.get(t_tx, '?')
            attrs = f"ID={name};reference_gene_id={ref_gene_id};source=cesar2;intactness={intactness};orthology_class={orth_class}"
            f.write(f"{chrom}\tTOGA2\tprotein_coding_gene\t{int(start)+1}\t{end}\t.\t{strand}\t.\t{attrs}\n")
            n_cs += 1

    # 7) Summary stats
    from collections import Counter
    src_counts = Counter(r['source'] for r in rows)
    intact_counts = Counter(r['intactness'] for r in rows)
    print(f"[{args.query}] {len(rows)} rows; sources={dict(src_counts)} intactness={dict(intact_counts)}")
    print(f"  GFF: {n_lo} liftoff genes + {n_cs} CESAR2 genes → {out_gff}")
    print(f"  Classification: {out_cls}")


if __name__ == '__main__':
    main()
