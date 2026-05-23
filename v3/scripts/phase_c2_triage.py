#!/usr/bin/env python3
"""
Phase C.2 triage script for the P. vivax orthology pipeline (v3.1).

For each gene projected by Liftoff, decide whether the projection is
clean enough to accept as-is or whether it needs to be re-projected by
CESAR2 (the slow but rigorous codon-aware fallback in Phase C.3).

Eight triage rules are applied. A gene is flagged for CESAR2 fallback
if ANY rule fires:

    R1 — frame disruption (valid_ORFs=0, CDS length not divisible by 3,
         or internal stop codon)
    R2 — sequence identity below threshold (default 0.95 for core genes,
         0.85 for known family genes)
    R3 — reference CDS coverage below threshold (default 0.90)
    R4 — copy-number variation (extra_copy_number > 0 for a non-family gene)
    R5 — partial mapping flagged by Liftoff
    R6 — splice-site disruption (non-canonical GT/AG or AT/AC dinucleotides)
    R7 — subtelomeric location (within flank_bp of chromosome end)
    R8 — known family membership (force CESAR2 regardless of clean Liftoff)

Outputs:
    triage.tsv            per-gene decision and rules triggered
    needs_cesar2.bed      BED of flagged genes (in reference coords) for TOGA2
    liftoff_clean.gff3    Liftoff GFF filtered to non-flagged genes
    summary.json          counts and fallback rate

Requires: Python 3.9+, pyfaidx
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    from pyfaidx import Fasta
except ImportError:
    sys.exit("Missing dependency: pyfaidx. Install with `pip install pyfaidx` "
             "or `conda install -c bioconda pyfaidx`.")

LOG = logging.getLogger("triage")

# Plasmodium uses the standard genetic code (NCBI translation table 1).
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Canonical splice-site dinucleotide pairs: (donor, acceptor).
# GT-AG is the major spliceosome; AT-AC is the minor (U12) spliceosome.
CANONICAL_SPLICE = {("GT", "AG"), ("AT", "AC")}

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# Liftoff appends _<integer> to extra-copy IDs.
EXTRA_COPY_RE = re.compile(r"^(.+?)_(\d+)$")


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def parse_gff_attributes(attr_str: str) -> dict:
    """Parse the GFF3 attributes column (col 9) into a dict."""
    d = {}
    for kv in attr_str.strip().rstrip(';').split(';'):
        kv = kv.strip()
        if '=' in kv:
            k, v = kv.split('=', 1)
            d[k.strip()] = v.strip()
    return d


class GeneRecord:
    """In-memory representation of a Liftoff-projected gene."""
    __slots__ = ('gene_id', 'reference_id', 'chrom', 'start', 'end',
                 'strand', 'attrs', 'transcripts')

    def __init__(self, gene_id, chrom, start, end, strand, attrs):
        self.gene_id = gene_id
        self.reference_id = normalize_gene_id(gene_id)
        self.chrom = chrom
        self.start = start   # 1-based GFF3 coordinates
        self.end = end
        self.strand = strand
        self.attrs = attrs
        self.transcripts = []


def normalize_gene_id(gid: str) -> str:
    """Strip Liftoff extra-copy suffix (_1, _2, ...) to recover the
    reference gene ID."""
    m = EXTRA_COPY_RE.match(gid)
    if not m:
        return gid
    core, suffix = m.group(1), m.group(2)
    if len(suffix) <= 2 and not core.endswith('_'):
        return core
    return gid


def parse_liftoff_gff(path: Path) -> list:
    """Parse a Liftoff GFF3 file into a list of GeneRecord objects."""
    genes = {}
    tx_to_gene = {}
    tx_features = defaultdict(lambda: {"exon": [], "CDS": []})

    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, _src, ftype, start, end, _score, strand, _phase, attr_str = parts[:9]
            attrs = parse_gff_attributes(attr_str)
            start, end = int(start), int(end)

            if ftype in ("gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"):
                gid = attrs.get("ID")
                if not gid:
                    continue
                genes[gid] = GeneRecord(gid, chrom, start, end, strand, attrs)
            elif ftype in ("mRNA", "transcript", "pseudogenic_transcript",
                           "tRNA", "rRNA", "ncRNA", "snoRNA", "snRNA"):
                tx_id = attrs.get("ID")
                parent = attrs.get("Parent")
                if tx_id and parent:
                    tx_to_gene[tx_id] = parent
            elif ftype in ("exon", "CDS"):
                parent = attrs.get("Parent")
                if parent:
                    tx_features[parent][ftype].append((start, end))

    for tx_id, gid in tx_to_gene.items():
        if gid in genes:
            exons = sorted(tx_features[tx_id]["exon"])
            cdss = sorted(tx_features[tx_id]["CDS"])
            genes[gid].transcripts.append((tx_id, exons, cdss))

    return list(genes.values())


# ---------------------------------------------------------------------------
# Sequence extraction and codon checks
# ---------------------------------------------------------------------------

def extract_sequence(fa, chrom, start, end, strand) -> str:
    """Extract sequence from 1-based inclusive GFF coordinates."""
    seq = str(fa[chrom][start - 1:end]).upper()
    if strand == '-':
        seq = seq.translate(COMPLEMENT)[::-1]
    return seq


def extract_cds(fa, chrom, cds_segments, strand) -> str:
    """Concatenate CDS segments in transcription order."""
    segs = sorted(cds_segments, reverse=(strand == '-'))
    return ''.join(extract_sequence(fa, chrom, s, e, strand) for s, e in segs)


def has_internal_stop(cds_nt: str) -> bool:
    """True if any in-frame stop codon occurs before the final codon."""
    if len(cds_nt) % 3 != 0 or len(cds_nt) < 6:
        return False
    n_codons = len(cds_nt) // 3
    for i in range(n_codons - 1):
        if cds_nt[i * 3:(i + 1) * 3] in STOP_CODONS:
            return True
    return False


def get_splice_sites(fa, chrom, exons, strand) -> list:
    """Extract (donor, acceptor) dinucleotides for each intron."""
    if len(exons) < 2:
        return []

    exons_sorted = sorted(exons, reverse=(strand == '-'))
    sites = []
    for i in range(len(exons_sorted) - 1):
        if strand == '+':
            intron_start = exons_sorted[i][1] + 1
            intron_end = exons_sorted[i + 1][0] - 1
            if intron_end < intron_start + 3:
                continue
            donor = extract_sequence(fa, chrom, intron_start, intron_start + 1, '+')
            acceptor = extract_sequence(fa, chrom, intron_end - 1, intron_end, '+')
        else:
            intron_start = exons_sorted[i + 1][1] + 1
            intron_end = exons_sorted[i][0] - 1
            if intron_end < intron_start + 3:
                continue
            donor = extract_sequence(fa, chrom, intron_end - 1, intron_end, '-')
            acceptor = extract_sequence(fa, chrom, intron_start, intron_start + 1, '-')
        sites.append((donor, acceptor))
    return sites


def is_subtelomeric(chrom, start, end, chrom_sizes, flank_bp) -> bool:
    size = chrom_sizes.get(chrom)
    if not size:
        return False
    return start < flank_bp or end > size - flank_bp


# ---------------------------------------------------------------------------
# Auxiliary file readers
# ---------------------------------------------------------------------------

def read_family_list(path) -> dict:
    """Read TSV (gene_id <tab> family_name)."""
    fams = {}
    if not path:
        return fams
    p = Path(path)
    if not p.exists():
        LOG.warning(f"Family list not found at {p}; rule R8 will not fire")
        return fams
    with open(p) as f:
        for row in csv.reader(f, delimiter='\t'):
            if not row or row[0].startswith('#') or row[0] == 'gene_id':
                continue
            if len(row) >= 2:
                fams[row[0]] = row[1]
    return fams


def read_reference_bed(path) -> dict:
    """Read reference BED into dict of gene_id -> full BED line."""
    bed_lines = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 4:
                bed_lines[parts[3]] = line
    return bed_lines


# ---------------------------------------------------------------------------
# Core triage logic
# ---------------------------------------------------------------------------

def triage_gene(gene: GeneRecord, fa, chrom_sizes, family_list, args):
    """Apply all 8 rules to a single gene. Returns (triggers, is_family, family_name)."""
    triggers = []
    family_membership = family_list.get(gene.reference_id) or family_list.get(gene.gene_id)
    is_family = family_membership is not None

    # R8 — known family membership
    if is_family:
        triggers.append("R8_family")

    # R1 — frame disruption
    valid_orfs_attr = gene.attrs.get("valid_ORFs", gene.attrs.get("valid_ORF", ""))
    if valid_orfs_attr in {"0", "False", "false"}:
        triggers.append("R1a_valid_ORF_flag")

    frame_bad = False
    internal_stop = False
    for _tx_id, _exons, cdss in gene.transcripts:
        if not cdss:
            continue
        try:
            cds_nt = extract_cds(fa, gene.chrom, cdss, gene.strand)
        except (KeyError, ValueError) as e:
            LOG.debug(f"CDS extract failed for {gene.gene_id}: {e}")
            continue
        if len(cds_nt) % 3 != 0:
            frame_bad = True
        if has_internal_stop(cds_nt):
            internal_stop = True
    if frame_bad and not any(t.startswith("R1a") for t in triggers):
        triggers.append("R1b_cds_length")
    if internal_stop:
        triggers.append("R1c_internal_stop")

    # R2 — identity threshold
    try:
        seq_id = float(gene.attrs.get("sequence_ID", "1.0"))
    except ValueError:
        seq_id = 1.0
    id_min = args.family_identity_min if is_family else args.core_identity_min
    if seq_id < id_min:
        triggers.append(f"R2_identity_{seq_id:.3f}")

    # R3 — coverage threshold
    try:
        coverage = float(gene.attrs.get("coverage", "1.0"))
    except ValueError:
        coverage = 1.0
    if coverage < args.core_coverage_min:
        triggers.append(f"R3_coverage_{coverage:.3f}")

    # R4 — copy-number variation
    try:
        extra_copies = int(gene.attrs.get("extra_copy_number", "0"))
    except ValueError:
        extra_copies = 0
    if extra_copies > 0 and not is_family:
        triggers.append(f"R4_extra_copies_{extra_copies}")

    # R5 — partial mapping
    if gene.attrs.get("partial_mapping", "").lower() == "true":
        triggers.append("R5_partial")

    # R6 — splice-site disruption
    splice_bad = None
    for _tx_id, exons, _cdss in gene.transcripts:
        if len(exons) < 2:
            continue
        try:
            sites = get_splice_sites(fa, gene.chrom, exons, gene.strand)
        except (KeyError, ValueError):
            continue
        for donor, acceptor in sites:
            if (donor, acceptor) not in CANONICAL_SPLICE:
                splice_bad = (donor, acceptor)
                break
        if splice_bad:
            break
    if splice_bad:
        triggers.append(f"R6_splice_{splice_bad[0]}_{splice_bad[1]}")

    # R7 — subtelomeric location
    if is_subtelomeric(gene.chrom, gene.start, gene.end,
                       chrom_sizes, args.subtelomere_bp):
        triggers.append("R7_subtelomeric")

    return triggers, is_family, family_membership


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Phase C.2 triage: decide which Liftoff-projected genes "
                    "need CESAR2 fallback in Phase C.3.")
    parser.add_argument("--liftoff-gff", required=True)
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--reference-bed", required=True)
    parser.add_argument("--family-list", default=None)
    parser.add_argument("--subtelomere-bp", type=int, default=100_000)
    parser.add_argument("--core-identity-min", type=float, default=0.95)
    parser.add_argument("--core-coverage-min", type=float, default=0.90)
    parser.add_argument("--family-identity-min", type=float, default=0.85)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--query-name", required=True)
    parser.add_argument("--log-level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s %(levelname)s %(message)s")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    LOG.info(f"Indexing query FASTA: {args.query_fasta}")
    fa = Fasta(args.query_fasta)
    chrom_sizes = {name: len(fa[name]) for name in fa.keys()}
    LOG.info(f"  {len(chrom_sizes)} sequences, "
             f"total {sum(chrom_sizes.values())/1e6:.1f} Mb")

    family_list = read_family_list(args.family_list)
    LOG.info(f"Loaded {len(family_list)} family-gene assignments")

    LOG.info(f"Parsing Liftoff GFF: {args.liftoff_gff}")
    genes = parse_liftoff_gff(args.liftoff_gff)
    LOG.info(f"  parsed {len(genes)} genes")

    ref_bed_lines = read_reference_bed(args.reference_bed)
    LOG.info(f"  {len(ref_bed_lines)} reference BED entries")

    triage_rows = []
    flagged_ref_ids = set()
    clean_gene_ids = set()
    rule_counter = defaultdict(int)

    for gene in genes:
        triggers, is_family, fam = triage_gene(
            gene, fa, chrom_sizes, family_list, args
        )
        if triggers:
            flagged_ref_ids.add(gene.reference_id)
            for t in triggers:
                key = '_'.join(t.split('_', 2)[:2]) if '_' in t else t
                rule_counter[key] += 1
        else:
            clean_gene_ids.add(gene.gene_id)

        triage_rows.append({
            "gene_id": gene.gene_id,
            "reference_id": gene.reference_id,
            "chrom": gene.chrom,
            "start": gene.start,
            "end": gene.end,
            "strand": gene.strand,
            "is_family": is_family,
            "family": fam or "",
            "sequence_ID": gene.attrs.get("sequence_ID", ""),
            "coverage": gene.attrs.get("coverage", ""),
            "extra_copy_number": gene.attrs.get("extra_copy_number", "0"),
            "valid_ORFs": gene.attrs.get("valid_ORFs",
                                         gene.attrs.get("valid_ORF", "")),
            "decision": "CESAR2_FALLBACK" if triggers else "LIFTOFF_OK",
            "rules_triggered": ",".join(triggers),
        })

    # Write triage.tsv
    triage_path = outdir / "triage.tsv"
    with open(triage_path, 'w', newline='') as fout:
        writer = csv.DictWriter(fout, fieldnames=list(triage_rows[0].keys()),
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(triage_rows)

    # Write needs_cesar2.bed
    cesar2_bed_path = outdir / "needs_cesar2.bed"
    n_written = n_missing = 0
    with open(cesar2_bed_path, 'w') as fout:
        for ref_id in sorted(flagged_ref_ids):
            if ref_id in ref_bed_lines:
                fout.write(ref_bed_lines[ref_id])
                n_written += 1
            else:
                n_missing += 1
    LOG.info(f"BED: {n_written} written, {n_missing} flagged but not in reference BED")

    # Write liftoff_clean.gff3
    clean_gff_path = outdir / "liftoff_clean.gff3"
    current_gene_id = None
    write_current = True
    n_lines_kept = 0
    with open(args.liftoff_gff) as fin, open(clean_gff_path, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#'):
                fout.write(line)
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            ftype = parts[2]
            attrs = parse_gff_attributes(parts[8])
            if ftype in ("gene", "protein_coding_gene", "ncRNA_gene", "pseudogene"):
                current_gene_id = attrs.get("ID")
                write_current = current_gene_id in clean_gene_ids
            if write_current:
                fout.write(line)
                n_lines_kept += 1

    # Summary
    summary = {
        "query": args.query_name,
        "total_genes": len(genes),
        "liftoff_clean": len(clean_gene_ids),
        "needs_cesar2": len(flagged_ref_ids),
        "needs_cesar2_in_bed": n_written,
        "needs_cesar2_missing_from_bed": n_missing,
        "fallback_rate": (len(flagged_ref_ids) / len(genes)) if genes else 0.0,
        "rule_counts": dict(rule_counter),
        "thresholds": {
            "core_identity_min": args.core_identity_min,
            "core_coverage_min": args.core_coverage_min,
            "family_identity_min": args.family_identity_min,
            "subtelomere_bp": args.subtelomere_bp,
        },
    }
    with open(outdir / "summary.json", 'w') as fout:
        json.dump(summary, fout, indent=2, sort_keys=True)

    LOG.info(f"Triage complete for {args.query_name}: "
             f"{summary['needs_cesar2']}/{summary['total_genes']} "
             f"({summary['fallback_rate']:.1%}) flagged for CESAR2 fallback")

    if summary['fallback_rate'] > 0.50:
        LOG.warning("Fallback rate exceeds 50% — triage may be too strict.")
    if summary['fallback_rate'] < 0.02 and len(family_list) > 0:
        LOG.warning("Fallback rate below 2% despite family list — check coverage.")


if __name__ == "__main__":
    main()
