#!/usr/bin/env python3
"""
build_msa.py — codon + protein MSA builder for the Pk v1 pipeline.

For each orthogroup in ortholog_table.tsv with >= min_intact strains intact:
  1. Extract per-strain CDS via GFF + FASTA (pyfaidx, no docker per-gene).
  2. Translate to protein.
  3. Protein MSA via MAFFT-LINSI (via cmd wrapper).
  4. Codon back-translation via pal2nal.
  5. Write {gene}.codon.aln.fa + {gene}.protein.aln.fa

Supports --shard / --num-shards for parallelism across multiple processes.
Idempotent: skips genes whose .codon.aln.fa already exists in out_dir.

Parameterized from build_8way_msa_v2.py
  (/media/anton/data/sandbox/Pv4/v3/scripts/build_8way_msa_v2.py)
to accept dynamic strain set and min_intact rather than hard-coded 8-way.
"""

import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

try:
    from pyfaidx import Fasta
except ImportError:
    sys.exit("pyfaidx required — run inside the bcftools container or install pyfaidx")

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(cds: str) -> str:
    aa = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i + 3].upper()
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def strip_internal_stops(cds: str) -> str:
    """Replace in-frame internal stop codons (not the final codon) with NNN."""
    if len(cds) < 3:
        return cds
    codons = [cds[i:i + 3] for i in range(0, len(cds), 3)]
    n_full = len(cds) // 3
    fixed = []
    for i, codon in enumerate(codons[:n_full - 1]):  # all but last full codon
        fixed.append('NNN' if codon.upper() in STOP_CODONS else codon)
    fixed.extend(codons[n_full - 1:])  # last codon + any remainder
    return ''.join(fixed)


def parse_gff_cds(gff_path, target_genes=None):
    """Return dict gene_id -> [(chrom, start, end, strand, phase, parent)]."""
    out: dict = defaultdict(list)
    if not Path(gff_path).exists():
        return out
    tx_to_gene: dict = {}
    GENE_TYPES = {'gene', 'protein_coding_gene', 'ncRNA_gene', 'pseudogene'}
    TX_TYPES = {'mRNA', 'transcript', 'pseudogenic_transcript'}
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype in TX_TYPES:
                attrs = {kv.split('=', 1)[0]: kv.split('=', 1)[1]
                         for kv in f[8].split(';') if '=' in kv}
                tx_id = attrs.get('ID')
                parent = attrs.get('Parent')
                if tx_id and parent:
                    tx_to_gene[tx_id] = parent
    with open(gff_path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'CDS':
                continue
            attrs = {kv.split('=', 1)[0]: kv.split('=', 1)[1]
                     for kv in f[8].split(';') if '=' in kv}
            parent = attrs.get('Parent', '')
            gene_id = tx_to_gene.get(parent, parent.rsplit('.', 1)[0])
            # Normalize Liftoff extra-copy suffix
            m = re.match(r'^(.+)_(\d+)$', gene_id)
            if m and len(m.group(2)) <= 2 and not m.group(1).endswith('_'):
                gene_id = m.group(1)
            if target_genes is not None and gene_id not in target_genes:
                continue
            chrom = f[0]
            start = int(f[3])
            end = int(f[4])
            strand = f[6]
            phase = int(f[7]) if f[7] != '.' else 0
            out[gene_id].append((chrom, start, end, strand, phase, parent))
    return out


def extract_cds(fasta, segments) -> str:
    if not segments:
        return ''
    strand = segments[0][3]
    parent = segments[0][5]
    same_tx = [s for s in segments if s[5] == parent] or segments
    same_tx.sort(key=lambda x: x[1], reverse=(strand == '-'))
    parts = []
    for chrom, start, end, strd, phase, _p in same_tx:
        try:
            seq = str(fasta[chrom][start - 1:end]).upper()
        except (KeyError, IndexError):
            continue
        if strd == '-':
            seq = revcomp(seq)
        parts.append(seq)
    return ''.join(parts)


def run_mafft(input_fa: Path, output_fa: Path, threads: int = 2,
              work_dir: str = '') -> bool:
    """Run MAFFT --localpair --maxiterate 1000 via run_in_container.sh."""
    wrapper = Path(work_dir) / 'pipeline' / 'lib' / 'run_in_container.sh'
    if not wrapper.exists():
        wrapper = Path(work_dir) / 'pipeline' / 'lib' / 'run_in_apptainer.sh'
    cmd = [str(wrapper), 'mafft',
           '--localpair', '--maxiterate', '1000',
           '--thread', str(threads),
           str(input_fa)]
    with open(output_fa, 'w') as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL)
    return r.returncode == 0


def run_pal2nal(prot_aln: Path, nucl_fa: Path, codon_aln: Path,
                work_dir: str = '') -> bool:
    """Run pal2nal via run_in_container.sh."""
    wrapper = Path(work_dir) / 'pipeline' / 'lib' / 'run_in_container.sh'
    if not wrapper.exists():
        wrapper = Path(work_dir) / 'pipeline' / 'lib' / 'run_in_apptainer.sh'
    cmd = [str(wrapper), 'pal2nal',
           'pal2nal.pl',
           str(prot_aln), str(nucl_fa),
           '-output', 'fasta', '-codontable', '1']
    with open(codon_aln, 'w') as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL)
    return r.returncode == 0


def load_ortho_table(ortho_path, all_strains, min_intact, ref_strain):
    """Return set of gene IDs (reference strain) that have >= min_intact intact strains."""
    target_genes: set = set()
    with open(ortho_path) as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            ref_val = row.get(ref_strain, '-')
            if ref_val == '-':
                continue
            # Count non-dash, non-empty strains
            n_present = sum(1 for s in all_strains if row.get(s, '-') not in ('-', ''))
            if n_present >= min_intact:
                # ref_val may be 'gene_id' or 'gene_id|alias'
                for gid in ref_val.split(','):
                    gene_id = gid.split('|')[0].strip()
                    if gene_id:
                        target_genes.add(gene_id)
    return target_genes


def main():
    ap = argparse.ArgumentParser(description="Build codon + protein MSAs")
    ap.add_argument('--ortho', required=True, help='work/03_consensus/ortholog_table.tsv')
    ap.add_argument('--strains', required=True, help='Space-separated all-strain list')
    ap.add_argument('--ref', required=True, help='Reference strain name')
    ap.add_argument('--ref-gff', required=True, help='inputs/annotations/{REF}.fixed.gff3')
    ap.add_argument('--ref-fasta', required=True, help='genomes/softmasked/{REF}.fa')
    ap.add_argument('--queries', nargs='+', required=True,
                    help='Query strain names in order matching --query-gff and --query-fasta')
    ap.add_argument('--query-gff', nargs='+', required=True,
                    help='Per-query merged annotation GFF paths')
    ap.add_argument('--query-fasta', nargs='+', required=True,
                    help='Per-query softmasked FASTA paths')
    ap.add_argument('--merged-base', required=True,
                    help='work/02d_merged/{REF}-as-ref/ directory')
    ap.add_argument('--min-intact', type=int, default=6)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--threads-per-gene', type=int, default=2)
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--num-shards', type=int, default=1)
    ap.add_argument('--max-genes', type=int, default=0,
                    help='Limit to first N genes (0 = all; for testing)')
    ap.add_argument('--work-dir', default='',
                    help='WORK root for locating run_in_container.sh wrapper')
    args = ap.parse_args()

    all_strains = args.strains.split()
    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Determine WORK dir for wrapper discovery
    work_dir = args.work_dir or str(Path(args.out_dir).parents[3])

    # 1) Target gene set
    target_genes = load_ortho_table(args.ortho, all_strains, args.min_intact, args.ref)
    print(f"Target: {len(target_genes)} genes intact in >= {args.min_intact} strains",
          file=sys.stderr)
    if args.max_genes > 0:
        target_genes = set(sorted(target_genes)[:args.max_genes])
        print(f"  (truncated to first {len(target_genes)} for testing)", file=sys.stderr)

    # 2) Load reference GFF + FASTA
    print(f"Loading reference GFF: {args.ref_gff}", file=sys.stderr)
    ref_cds_map = parse_gff_cds(args.ref_gff, target_genes)
    ref_fa = Fasta(args.ref_fasta)

    # 3) Per-query GFF + FASTA
    query_cds_maps: dict = {}
    query_fastas: dict = {}
    for q, gff, fa_path in zip(args.queries, args.query_gff, args.query_fasta):
        query_cds_maps[q] = parse_gff_cds(gff, target_genes)
        query_fastas[q] = Fasta(fa_path)
        print(f"  {q:12s}  {len(query_cds_maps[q])} target genes found in GFF",
              file=sys.stderr)

    # 4) Build MSAs
    tmp_dir = outdir / f'tmp_shard{args.shard}'
    tmp_dir.mkdir(exist_ok=True)
    done_set = {p.stem.replace('.codon.aln', '') for p in outdir.glob('*.codon.aln.fa')}
    sorted_genes = sorted(target_genes)
    n_ok = n_skip = 0
    summary_rows = []

    for i, gene_id in enumerate(sorted_genes):
        if args.num_shards > 1 and (i % args.num_shards) != args.shard:
            continue
        if gene_id in done_set:
            continue
        if (i + 1) % 100 == 0:
            print(f"  ... {i+1}/{len(sorted_genes)} (ok={n_ok}, skip={n_skip})",
                  file=sys.stderr)

        ref_segs = ref_cds_map.get(gene_id, [])
        ref_cds = extract_cds(ref_fa, ref_segs) if ref_segs else ''
        if not ref_cds or len(ref_cds) % 3 != 0:
            n_skip += 1
            continue
        ref_prot = translate(ref_cds).rstrip('*')
        if '*' in ref_prot:   # internal stop in reference
            n_skip += 1
            continue

        seqs_cds = {args.ref: ref_cds}
        seqs_prot = {args.ref: ref_prot}
        missing = []

        for q in args.queries:
            q_segs = query_cds_maps[q].get(gene_id, [])
            q_cds = extract_cds(query_fastas[q], q_segs) if q_segs else ''
            if not q_cds:
                missing.append(q)
                continue
            # Truncate to nearest codon
            q_cds = q_cds[:(len(q_cds) // 3) * 3]
            # Strip internal stops (F.2 lesson)
            q_cds = strip_internal_stops(q_cds)
            q_prot = translate(q_cds).rstrip('*')
            seqs_cds[q] = q_cds
            seqs_prot[q] = q_prot

        if len(seqs_cds) - 1 < args.min_intact:
            n_skip += 1
            continue

        # Write to tmp_dir for docker volume mount
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

        if not run_mafft(prot_fa, prot_aln, threads=args.threads_per_gene,
                         work_dir=work_dir):
            n_skip += 1
            continue

        if not run_pal2nal(prot_aln, nucl_fa, codon_aln_tmp, work_dir=work_dir):
            n_skip += 1
            codon_aln_tmp.unlink(missing_ok=True)
            continue

        if not codon_aln_tmp.exists() or codon_aln_tmp.stat().st_size == 0:
            n_skip += 1
            codon_aln_tmp.unlink(missing_ok=True)
            continue

        # Validate codon MSA length = 3 × protein MSA length
        codon_lines = codon_aln_tmp.read_text().splitlines()
        prot_lines = prot_aln.read_text().splitlines()
        codon_seq = next((l for l in codon_lines if not l.startswith('>')), '')
        prot_seq = next((l for l in prot_lines if not l.startswith('>')), '')
        if codon_seq and prot_seq and len(codon_seq) != 3 * len(prot_seq):
            n_skip += 1
            codon_aln_tmp.unlink(missing_ok=True)
            continue

        # Move to final output
        codon_aln = outdir / f"{gene_id}.codon.aln.fa"
        prot_aln_out = outdir / f"{gene_id}.protein.aln.fa"
        prot_aln.rename(prot_aln_out)
        codon_aln_tmp.rename(codon_aln)
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

    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)

    sum_path = outdir / f'summary_shard{args.shard}.tsv'
    with open(sum_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh,
                           fieldnames=['gene_id', 'n_strains', 'missing_strains',
                                       'ref_cds_len', 'ref_prot_len'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)

    print(f"\nDone shard {args.shard}/{args.num_shards}: "
          f"{n_ok} MSAs built, {n_skip} skipped", file=sys.stderr)


if __name__ == '__main__':
    main()
