#!/usr/bin/env python3
"""
8-way codon MSA builder, v2.

For each PvP01 gene intact in >= K queries:
  1. Extract PvP01 reference CDS from inputs/annotations/PvP01.genbank.gff3 + PvP01.fa
  2. For each of 7 queries, extract query CDS:
       - If gene has Liftoff projection in {Q}.lifted.gff3 → use Liftoff coords + query FASTA
       - Else if gene has CESAR2 nucleotide.fa.gz entry → use that
       - Else mark as missing (gap row)
  3. Translate all 8 CDS to protein
  4. MAFFT-LINSI protein MSA
  5. pal2nal codon back-translation → 8-way codon MSA

Output: out_dir/{gene_id}.codon.aln.fa + .protein.aln.fa + summary.tsv

Note: this extracts via pyfaidx + GFF parsing rather than gffread,
avoiding docker per-gene calls.
"""

import argparse
import csv
import gzip
import sys
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    from pyfaidx import Fasta
except ImportError:
    sys.exit("pyfaidx required")

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def translate(cds):
    """Translate CDS to protein (X for ambiguous, * for stop)."""
    aa = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3].upper()
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def parse_gff_cds(gff_path, target_genes):
    """Return dict gene_id -> [(chrom, start, end, strand, frame), ...] (CDS segments)."""
    out = defaultdict(list)
    if not Path(gff_path).exists():
        return out
    # First pass: build transcript→gene map
    tx_to_gene = {}
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype in ('mRNA','transcript','pseudogenic_transcript'):
                attrs = {kv.split('=',1)[0]: kv.split('=',1)[1]
                         for kv in f[8].split(';') if '=' in kv}
                tx_id = attrs.get('ID')
                parent = attrs.get('Parent')
                if tx_id and parent:
                    tx_to_gene[tx_id] = parent
    # Second pass: gather CDS
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            attrs = {kv.split('=',1)[0]: kv.split('=',1)[1]
                     for kv in f[8].split(';') if '=' in kv}
            parent = attrs.get('Parent', '')
            # Parent may be transcript or gene
            # Liftoff: Parent=PVP01_xxx.1 (transcript)
            # Strip "_N" extra-copy suffix
            gene_id = tx_to_gene.get(parent, parent.rsplit('.', 1)[0])
            # Normalize Liftoff extra-copy suffix on gene_id
            if '_' in gene_id:
                parts = gene_id.rsplit('_', 1)
                if len(parts[1]) <= 2 and parts[1].isdigit() and not parts[0].endswith('_'):
                    gene_id = parts[0]
            if target_genes is not None and gene_id not in target_genes:
                continue
            chrom = f[0]
            start = int(f[3])
            end = int(f[4])
            strand = f[6]
            phase = int(f[7]) if f[7] != '.' else 0
            out[gene_id].append((chrom, start, end, strand, phase, parent))
    return out


def extract_cds(fasta, segments):
    """Extract CDS sequence given list of (chrom, start, end, strand, phase, parent).
    Returns the concatenated CDS, oriented in transcription direction.
    """
    if not segments:
        return ''
    # All segments share strand
    strand = segments[0][3]
    parent = segments[0][5]
    # Filter to single transcript (e.g., .1)
    same_tx = [s for s in segments if s[5] == parent]
    if not same_tx:
        same_tx = segments
    # Sort: + strand → start ascending; - strand → start descending
    same_tx.sort(key=lambda x: x[1], reverse=(strand == '-'))
    parts = []
    for chrom, start, end, strd, phase, _p in same_tx:
        try:
            seq = str(fasta[chrom][start-1:end]).upper()
        except (KeyError, IndexError):
            continue
        if strd == '-':
            seq = revcomp(seq)
        parts.append(seq)
    return ''.join(parts)


def run_mafft(input_fa, output_fa, threads=2):
    """Run MAFFT --auto on input_fa, write to output_fa."""
    cmd = ['docker', 'run', '--rm', '-v', f'{Path(input_fa).parent.absolute()}:/d',
           '-w', '/d', 'quay.io/biocontainers/mafft:7.313--0',
           'mafft', '--auto', '--thread', str(threads),
           f'/d/{Path(input_fa).name}']
    with open(output_fa, 'w') as out_fh:
        r = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.DEVNULL)
    return r.returncode == 0


def run_pal2nal(prot_aln, nucl_fa, codon_aln):
    """Run pal2nal protein-aln nucl → codon-aln."""
    cmd = ['docker', 'run', '--rm', '-v', f'{Path(prot_aln).parent.absolute()}:/d',
           '-w', '/d', 'quay.io/biocontainers/pal2nal:14.1--pl5321hdfd78af_3',
           'pal2nal.pl', f'/d/{Path(prot_aln).name}', f'/d/{Path(nucl_fa).name}',
           '-output', 'fasta', '-codontable', '1']
    with open(codon_aln, 'w') as out_fh:
        r = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.DEVNULL)
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--queries', nargs='+', required=True)
    ap.add_argument('--ref-gff', required=True)
    ap.add_argument('--ref-fasta', required=True)
    ap.add_argument('--query-gff', nargs='+', required=True,
                    help='In order matching --queries: per-query annotation GFF')
    ap.add_argument('--query-fasta', nargs='+', required=True)
    ap.add_argument('--merged-base', required=True)
    ap.add_argument('--min-intact', type=int, default=7)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--max-genes', type=int, default=0,
                    help='Limit to first N genes (for testing); 0 = all')
    ap.add_argument('--threads-per-gene', type=int, default=2)
    ap.add_argument('--shard', type=int, default=0,
                    help='Process only genes where gene_index %% num_shards == shard')
    ap.add_argument('--num-shards', type=int, default=1)
    args = ap.parse_args()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) Per-query intactness
    intact_per_q = {}
    for q in args.queries:
        cls_path = Path(args.merged_base) / f"{q}.classification.tsv"
        s = set()
        with open(cls_path) as fh:
            r = csv.DictReader(fh, delimiter='\t')
            for row in r:
                if row['intactness'] in ('I', 'PI'):
                    s.add(row['reference_gene_id'])
        intact_per_q[q] = s
    all_genes = set().union(*intact_per_q.values())
    gene_count = {g: sum(1 for s in intact_per_q.values() if g in s) for g in all_genes}
    target_genes = {g for g, n in gene_count.items() if n >= args.min_intact}
    print(f"Target: {len(target_genes)} genes intact in >= {args.min_intact} queries",
          file=sys.stderr)
    if args.max_genes > 0:
        target_genes = set(sorted(target_genes)[:args.max_genes])
        print(f"  (truncated to first {len(target_genes)} for testing)", file=sys.stderr)

    # 2) Load reference GFF + FASTA, extract CDS per gene
    print(f"Loading reference {args.ref_gff} ...", file=sys.stderr)
    ref_cds_map = parse_gff_cds(args.ref_gff, target_genes)
    ref_fa = Fasta(args.ref_fasta)
    print(f"  ref CDS for {len(ref_cds_map)} genes", file=sys.stderr)

    # 3) Per-query GFF + FASTA: extract CDS for each gene from merged annotation
    print(f"Loading per-query GFFs ...", file=sys.stderr)
    query_cds_maps = {}
    query_fastas = {}
    for q, gff, fa in zip(args.queries, args.query_gff, args.query_fasta):
        cds_map = parse_gff_cds(gff, target_genes)
        query_cds_maps[q] = cds_map
        query_fastas[q] = Fasta(fa)
        print(f"  {q:8s}  {len(cds_map)} genes", file=sys.stderr)

    # 4) For each gene, extract 8 CDS, translate, write protein FASTA
    print(f"Building alignments ...", file=sys.stderr)
    summary_rows = []
    n_ok = 0
    n_skip = 0
    tmp_dir = outdir / f'tmp_shard{args.shard}'
    tmp_dir.mkdir(exist_ok=True)
    sorted_genes = sorted(target_genes)
    # Skip genes already done (by .codon.aln.fa presence) — supports resume
    done_set = {p.stem.replace('.codon.aln', '') for p in outdir.glob('*.codon.aln.fa')}
    for i, gene_id in enumerate(sorted_genes):
        if args.num_shards > 1 and (i % args.num_shards) != args.shard:
            continue
        if gene_id in done_set:
            continue
        if (i + 1) % 100 == 0:
            print(f"  ... {i+1} / {len(target_genes)} (ok={n_ok}, skip={n_skip})",
                  file=sys.stderr)
        ref_segs = ref_cds_map.get(gene_id, [])
        ref_cds = extract_cds(ref_fa, ref_segs) if ref_segs else ''
        if not ref_cds or len(ref_cds) % 3 != 0:
            n_skip += 1
            continue
        ref_prot = translate(ref_cds).rstrip('*')
        if '*' in ref_prot:  # internal stop
            n_skip += 1
            continue
        seqs_cds = {'PvP01_REF': ref_cds}
        seqs_prot = {'PvP01_REF': ref_prot}
        missing = []
        for q in args.queries:
            q_segs = query_cds_maps[q].get(gene_id, [])
            q_cds = extract_cds(query_fastas[q], q_segs) if q_segs else ''
            if not q_cds:
                missing.append(q)
                continue
            # Truncate to nearest codon
            q_cds = q_cds[:(len(q_cds) // 3) * 3]
            q_prot = translate(q_cds).rstrip('*')
            seqs_cds[q] = q_cds
            seqs_prot[q] = q_prot
        # Need at least min_intact queries
        if len(seqs_cds) - 1 < args.min_intact:
            n_skip += 1
            continue

        # Write nucl + protein FASTAs in tmp_dir (keep all in same dir for docker mount)
        nucl_fa = tmp_dir / f"{gene_id}.nucl.fa"
        prot_fa = tmp_dir / f"{gene_id}.prot.fa"
        prot_aln = tmp_dir / f"{gene_id}.prot.aln.fa"
        codon_aln_tmp = tmp_dir / f"{gene_id}.codon.aln.fa"

        with open(nucl_fa, 'w') as fh:
            for k, v in seqs_cds.items():
                fh.write(f">{k}\n{v}\n")
        with open(prot_fa, 'w') as fh:
            for k, v in seqs_prot.items():
                fh.write(f">{k}\n{v}\n")
        ok_mafft = run_mafft(prot_fa, prot_aln, threads=args.threads_per_gene)
        if not ok_mafft:
            n_skip += 1
            continue
        # pal2nal in tmp_dir (prot_aln + nucl_fa both there)
        ok_pal = run_pal2nal(prot_aln, nucl_fa, codon_aln_tmp)
        if not ok_pal or not codon_aln_tmp.exists() or codon_aln_tmp.stat().st_size == 0:
            n_skip += 1
            codon_aln_tmp.unlink(missing_ok=True)
            continue
        # Move final outputs to outdir
        codon_aln = outdir / f"{gene_id}.codon.aln.fa"
        protein_aln_out = outdir / f"{gene_id}.protein.aln.fa"
        prot_aln.rename(protein_aln_out)
        codon_aln_tmp.rename(codon_aln)
        # Clean nucl/prot inputs
        nucl_fa.unlink(missing_ok=True)
        prot_fa.unlink(missing_ok=True)
        n_ok += 1
        summary_rows.append({
            'gene_id': gene_id,
            'n_strains': len(seqs_cds),
            'missing_strains': ','.join(missing) if missing else '-',
            'ref_cds_len': len(ref_cds),
            'ref_prot_len': len(ref_prot),
        })

    # Clean tmp
    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)

    # Summary
    sum_path = outdir / f'summary_shard{args.shard}.tsv'
    with open(sum_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['gene_id','n_strains','missing_strains',
                                            'ref_cds_len','ref_prot_len'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nDone: {n_ok} MSAs built, {n_skip} skipped", file=sys.stderr)
    print(f"Output: {outdir}", file=sys.stderr)
    print(f"Summary: {sum_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
