# Orthology Reconstruction Across 8 *Plasmodium vivax* References
## v3.1 — Hybrid Liftoff + TOGA2/CESAR2 Pipeline with Embedded Triage Script

**Plan for Claude Code execution.** Multi-phase pipeline to compute orthologs and generate codon-aware MSAs across eight *P. vivax* reference assemblies, starting from a PGGB pangenome graph and all-vs-all lastz `.axt` alignments that already exist on disk.

**Change from v3.** Phase C.2 now includes the full triage script (`scripts/phase_c2_triage.py`) inline, along with invocation examples and test fixtures. All other phases unchanged.

---

## 0. Assemblies and annotation tiers

| accession | strain | Tier | PlasmoDB ID | annotation source |
|---|---|---|---|---|
| GCA_000002415.2 | Sal-I  | 1 (anchor-eligible but fragmented) | `PvivaxSal1` | PlasmoDB-68 |
| GCA_900093555.2 | PvP01  | 1 (primary anchor) | `PvivaxP01` | PlasmoDB-68 |
| GCA_914969965.1 | PvW1   | 1 | `PvivaxPvW1` | PlasmoDB-68 |
| GCA_949152365.1 | PAM    | 1 | `PvivaxPAM` | PlasmoDB-68 |
| GCA_003402215.1 | PvSY56 | 1 | `PvivaxPvSY56` | PlasmoDB-68 |
| GCA_900093545.1 | PvT01  | 2 | — | NCBI Datasets |
| GCA_900093535.1 | PvC01  | 2 | — | NCBI Datasets |
| GCA_040114635.1 | MHC087 | 2 | — | NCBI Datasets |

**Tier policy.** Tier 1 PlasmoDB-annotated → eligible as anchor (Sal-I excluded by default for fragmentation). Tier 2 NCBI-only → query-only, annotations derived by projection. MHC087 needs Phase A characterization before downstream use.

---

## 1. Goals and definition of success

**Primary goal.** Unified per-gene table classifying each gene as CORE-1:1 / CORE-VAR / FAMILY / LINEAGE-SPECIFIC / LOST across all 8 assemblies.

**Secondary goal.** Codon-aware MSAs for CORE-1:1 and CORE-VAR groups, ready for HyPhy and PAML.

**Tertiary goal.** Subfamily-resolved phylogenetic placement plus per-subfamily codon MSAs for FAMILY groups.

**Success criteria.**
- ≥ 90% of BUSCO `plasmodium_odb10` single-copy markers as CORE-1:1.
- ≥ 80% concordance with PlasmoDB OrthoMCL legacy calls on Tier-1 genomes.
- ≥ 70% concordance of pir/vir subfamily calls with published A–L scheme.
- All emitted codon MSAs in-frame, no internal stops, divisible by 3.
- **v3 audit check.** Spot-check 100 random Liftoff-only `I` calls by running CESAR2 on them; require ≥ 98% agreement. Halt if lower.

---

## 2. Inputs the pipeline expects on disk

```
inputs/
├── assemblies/
│   ├── PvP01.fa         # GCA_900093555.2
│   ├── PvW1.fa          # GCA_914969965.1
│   ├── PAM.fa           # GCA_949152365.1
│   ├── PvSY56.fa        # GCA_003402215.1
│   ├── Sal-I.fa         # GCA_000002415.2
│   ├── PvT01.fa         # GCA_900093545.1
│   ├── PvC01.fa         # GCA_900093535.1
│   └── MHC087.fa        # GCA_040114635.1
├── annotations/
│   ├── plasmodb-68/{PvP01,PvW1,PAM,PvSY56,Sal-I}.{gff3,proteins.fa}
│   ├── plasmodb-68/PvP01.family_list.tsv      # known family genes for triage rule R8
│   └── ncbi-datasets/{PvT01,PvC01,MHC087}.{gff3,proteins.fa}
├── lastz/
│   └── axt/                                    # 56 directed files
└── pggb/
    ├── pv.gfa
    ├── pv.og
    └── pv.smooth.vcf.gz
```

**`PvP01.family_list.tsv` schema.** TSV with header `gene_id\tfamily` (e.g. `PVP01_1248200\tpir`). Build from the PlasmoDB-68 functional annotation by selecting genes annotated as pir, vir, PHIST, MSP3, SERA, or PfRH-like.

If any input is missing, halt at Phase A.

---

## 3. Environment

```yaml
name: pv-ortho
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - liftoff>=1.6.3
  - minimap2>=2.28
  - pyfaidx>=0.7
  - ucsc-axtchain
  - ucsc-chainnet
  - ucsc-chaincleaner
  - ucsc-axttomaf
  - ucsc-mafgene
  - ucsc-chainswap
  - ucsc-chainstats
  - ucsc-fatotwobit
  - ucsc-gff3togenepred
  - ucsc-genepredtobed
  - multiz
  - lastz>=1.04
  - mafft
  - macse=2.07
  - hmmer>=3.4
  - diamond>=2.1
  - mcl
  - iqtree>=2.3
  - trimal
  - hyphy>=2.5.62
  - paml>=4.10
  - busco>=5.8
  - odgi>=0.9
  - vg>=1.55
  - bedtools
  - samtools
  - seqkit
  - gffread
  - gffcompare
  - agat
  - ncbi-datasets-cli
  - biopython
  - pandas
  - pyranges
  - networkx
  - mash
  - nextflow
```

TOGA2 and chainCleaner/RepeatFiller installed separately per Hiller lab instructions.

---

## 4. Working directory and repo layout

```
pv-ortho-pipeline/
├── scripts/
│   ├── phase_c2_triage.py           # the triage script (see Phase C.2 below)
│   ├── phase_c4_merge.py
│   ├── phase_e_consensus.py
│   └── ...
├── environment.yml
├── inputs/                           # symlink to data location
└── work/
    ├── 00_inventory/
    ├── 01_chains/
    ├── 02a_liftoff/
    ├── 02b_triage/
    ├── 02c_toga/
    ├── 02d_merged/
    ├── 03_consensus/
    ├── 04_pggb_xcheck/
    ├── 05_families/
    ├── 06_msa/
    ├── 07_msa_qc/
    ├── 08_selection_ready/
    └── logs/
        └── checkpoints/
```

---

## Phase A — Inventory, validation, harmonization

Unchanged from v3. FASTA stats, axt matrix check, BUSCO baseline, annotation harmonization via AGAT, MHC087 mash probe, PGGB graph validation.

**Additional Phase A step for v3.1.** Build `PvP01.family_list.tsv` by parsing the PlasmoDB-68 GFF and selecting genes with `product` or `description` matching the family regex: `(pir|vir|PHIST|MSP3|SERA|PfRH|RH[0-9]|EBP|MSP\d)`. Save to `inputs/annotations/plasmodb-68/PvP01.family_list.tsv`.

**Checkpoint.** `00_inventory.done`.

---

## Phase B — Chains from axt files

Unchanged. 56 directed cleaned chains + 28 reciprocal-best chains.

**Checkpoint.** `01_chains.done`.

---

## Phase C — Hybrid projection from PvP01 (primary anchor)

### C.1 — Liftoff fast pass

For each query Q ∈ {PvW1, PAM, PvSY56, Sal-I, PvT01, PvC01, MHC087}:

```bash
liftoff \
    -g inputs/annotations/plasmodb-68/PvP01.gff3 \
    -o work/02a_liftoff/PvP01-as-ref/${Q}/${Q}.lifted.gff3 \
    -u work/02a_liftoff/PvP01-as-ref/${Q}/${Q}.unmapped.txt \
    -dir work/02a_liftoff/PvP01-as-ref/${Q}/intermediate \
    -copies \
    -sc 0.90 \
    -d 5 \
    -flank 0.1 \
    -polish \
    -p 8 \
    inputs/assemblies/${Q}.fa \
    inputs/assemblies/PvP01.fa
```

Output GFF includes per-gene attributes: `sequence_ID`, `coverage`, `extra_copy_number`, `valid_ORFs`, `partial_mapping`.

### C.2 — Triage

The triage script applies eight rules and routes ~10–20% of genes to CESAR2 fallback while letting the clean majority be served by Liftoff.

**Rules.**

| Rule | Trigger condition |
|---|---|
| R1a | `valid_ORFs=0` or `valid_ORF=False` from Liftoff polish step |
| R1b | CDS length not divisible by 3 (any transcript) |
| R1c | Internal stop codon in CDS (any transcript) |
| R2  | `sequence_ID` < 0.95 for core genes (0.85 for known family genes) |
| R3  | `coverage` < 0.90 of reference CDS |
| R4  | `extra_copy_number` > 0 for a non-family gene |
| R5  | `partial_mapping=True` from Liftoff |
| R6  | Non-canonical splice-site dinucleotide (not GT/AG or AT/AC) at any intron boundary |
| R7  | Gene falls within ±100 kb of a chromosome end (subtelomeric) |
| R8  | Gene is in the family list (force CESAR2 regardless of clean Liftoff) |

**Invocation, per query.**

```bash
python scripts/phase_c2_triage.py \
    --liftoff-gff   work/02a_liftoff/PvP01-as-ref/${Q}/${Q}.lifted.gff3 \
    --query-fasta   inputs/assemblies/${Q}.fa \
    --reference-bed inputs/annotations/PvP01.bed \
    --family-list   inputs/annotations/plasmodb-68/PvP01.family_list.tsv \
    --output-dir    work/02b_triage/PvP01-as-ref/${Q}/ \
    --query-name    ${Q}
```

**Outputs.**
- `triage.tsv` — per-gene decision and triggered rules
- `needs_cesar2.bed` — BED of flagged reference genes for Phase C.3
- `liftoff_clean.gff3` — Liftoff GFF filtered to non-flagged genes
- `summary.json` — counts, fallback rate, per-rule statistics

**Sanity expectations on `summary.json`.**
- Fallback rate 10–25% for Tier-1 close-relative queries (PvW1, PAM, PvSY56)
- Fallback rate 25–40% for Sal-I (fragmentation drives R5)
- Fallback rate 10–25% for PvT01, PvC01
- MHC087 fallback rate depends on Phase A characterization
- If fallback rate < 2% with a populated family list, R8 isn't firing — investigate

**The script itself follows below for reference. Save as `scripts/phase_c2_triage.py`.**

```python
#!/usr/bin/env python3
"""
Phase C.2 triage script for the P. vivax orthology pipeline (v3).

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


def parse_liftoff_gff(path: Path) -> list[GeneRecord]:
    """Parse a Liftoff GFF3 file into a list of GeneRecord objects."""
    genes: dict[str, GeneRecord] = {}
    tx_to_gene: dict[str, str] = {}
    tx_features: dict[str, dict[str, list[tuple[int, int]]]] = defaultdict(
        lambda: {"exon": [], "CDS": []}
    )

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

            if ftype == "gene":
                gid = attrs.get("ID")
                if not gid:
                    continue
                genes[gid] = GeneRecord(gid, chrom, start, end, strand, attrs)
            elif ftype in ("mRNA", "transcript"):
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


def get_splice_sites(fa, chrom, exons, strand) -> list[tuple[str, str]]:
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

def read_family_list(path) -> dict[str, str]:
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
            if not row or row[0].startswith('#'):
                continue
            if len(row) >= 2:
                fams[row[0]] = row[1]
    return fams


def read_reference_bed(path) -> dict[str, str]:
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
    triggers: list[str] = []
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
    flagged_ref_ids: set[str] = set()
    clean_gene_ids: set[str] = set()
    rule_counter: dict[str, int] = defaultdict(int)

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
            if ftype == "gene":
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
```

**Test fixtures for the triage script.** Before the first full pipeline run, exercise the triage script against synthetic fixtures, one per rule:

```
tests/triage/
├── fixture_R1a_valid_ORF.gff3        # Liftoff valid_ORFs=0
├── fixture_R1b_length.gff3           # CDS length not divisible by 3
├── fixture_R1c_internal_stop.gff3    # TAA in middle of CDS
├── fixture_R2_low_identity.gff3      # sequence_ID=0.91
├── fixture_R3_low_coverage.gff3      # coverage=0.82
├── fixture_R4_extra_copies.gff3      # extra_copy_number=2, non-family
├── fixture_R5_partial.gff3           # partial_mapping=True
├── fixture_R6_bad_splice.gff3        # AT-AG (non-canonical)
├── fixture_R7_subtelomeric.gff3      # gene at chrom_start + 50kb
├── fixture_R8_family.gff3            # gene in family_list, clean elsewhere
├── fixture_clean.gff3                # passes all rules
└── tiny_query.fa                     # minimal FASTA for the fixtures above
```

Each fixture should trigger exactly the named rule (plus possibly R8 if the gene is in the family list). A small `pytest` harness verifies this. Failure on any fixture blocks the pipeline.

### C.3 — TOGA2/CESAR2 fallback

```bash
toga2.py \
    --target_2bit inputs/assemblies/PvP01.2bit \
    --query_2bit inputs/assemblies/${Q}.2bit \
    --chain work/01_chains/PvP01.${Q}.cleaned.chain \
    --reference_annotation work/02b_triage/PvP01-as-ref/${Q}/needs_cesar2.bed \
    --isoforms inputs/annotations/PvP01.isoforms.tsv \
    --output work/02c_toga/PvP01-as-ref/${Q}/ \
    --nextflow_config local
```

### C.4 — Merge Liftoff and CESAR2 outputs

Per query, produce a unified projected annotation in `work/02d_merged/PvP01-as-ref/${Q}.annotation.gff3` plus a flat classification table `${Q}.classification.tsv`. Rules:
- Genes in `liftoff_clean.gff3` → source=liftoff, intactness=I.
- Genes in CESAR2 output → source=cesar2, intactness from TOGA2 (I/PI/UL/L/PG/M).
- Genes flagged by triage but absent from CESAR2 output (rare) → source=cesar2, intactness=M, retain Liftoff coordinates with a warning.

**Validation.** Per-query: gene count ≈ PvP01 ± 5%; BUSCO single-copy markers ≥ 95% I (≥ 85% for Sal-I). Run the audit: spot-check 100 random Liftoff-only `I` calls by forcing CESAR2 on them; require ≥ 98% agreement. Halt below 98%.

**Checkpoint.** `02_projection_primary.done`.

---

## Phase D — Hybrid projection with reciprocal anchors

Repeat C.1–C.4 with each of PvW1, PAM, PvSY56 as anchor. Same triage script, just different `--liftoff-gff` and `--reference-bed` inputs per anchor.

**Skip option.** With Liftoff doing the bulk, Phase D drops to ~15–30 CPU-hours. Run by default.

**Checkpoint.** `02_projection_reciprocal.done`.

---

## Phase E — Consensus ortholog table

Unchanged from v3. Multigraph reconciliation across anchors, source-factor weighted edges (1.0 CESAR2-I, 0.95 Liftoff-clean, 0.7 CESAR2-PI, 0.4 UL/PG, 0 L/M), connected components, conflict resolution by weighted majority. Per-group provenance tracked.

**Output.** `work/03_consensus/ortholog_table.tsv`.

**Checkpoint.** `03_consensus.done`.

---

## Phase F — PGGB cross-validation

Unchanged. `odgi extract` per CORE-1:1 group, confirm 8-path traversal, flag disagreements.

**Checkpoint.** `04_pggb_xcheck.done`.

---

## Phase G — Gene family handling

Unchanged from v3. Family candidates from Phase E + force-routed R8 genes. HMM-based subfamily ID + all-vs-all DIAMOND + synteny-weighted MCL + phylogenetic placement.

**Checkpoint.** `05_families.done`.

---

## Phase H — Codon-aware MSA generation

Unchanged from v3. Three routes: Liftoff-projected → MAFFT-on-protein then MACSE back-translation; CESAR2-projected → CESAR2 codon alignment refined by MACSE; family clusters → direct MACSE.

**Checkpoint.** `06_msa.done`.

---

## Phase I — MSA QC

Unchanged. HmmCleaner, trimAl, frame/coverage/composition checks, 8-way completeness.

**Checkpoint.** `07_msa_qc.done`.

---

## Phase J — Selection analysis preparation

Unchanged. IQ-TREE with `-st CODON1`, HyPhy bundles, codeml control files with F3x4 default.

**Checkpoint.** `08_selection_ready.done`.

---

## Decisions deferred to operator

1. MHC087 identity and provenance.
2. lastz scoring matrix in original axt.
3. Sal-I as anchor — default no.
4. Phase D scope — all 3 reciprocal, subset, or skip.
5. Annotation precedence for Tier 2 — projection vs native NCBI.
6. pir A–L renaming.
7. pir outgroup (*P. cynomolgi*).
8. Selection test partitions.
9. Triage threshold tuning if the v3 audit (≥ 98% Liftoff/CESAR2 agreement) fails.

---

## Logging, checkpointing, resumability

Per-phase logs and structured `summary.json`. Checkpoint markers carry timestamp, git SHA, MD5 manifest. Resume on input MD5 change or `--force <phase>`. Triage provenance preserved per-gene through the entire pipeline.

---

## Expected runtime and resource footprint

| Phase | CPU-hours (v3.1) | v2 ref | Peak RAM |
|---|---|---|---|
| A | 2–4 | 2–4 | 8 GB |
| B | 12–25 | 12–25 | 8 GB |
| C.1 Liftoff | 1–2 | — | 8 GB |
| C.2 Triage | < 0.5 | — | 4 GB |
| C.3 CESAR2 | 4–8 | — | 16 GB |
| C.4 Merge | < 0.5 | — | 4 GB |
| **C total** | **6–11** | 30–55 | |
| D | 15–30 | 90–165 | 16 GB |
| E | 1 | 1 | 8 GB |
| F | 2 | 2 | 8 GB |
| G | 16–40 | 16–40 | 32 GB |
| H | 8–20 | 8–20 | 16 GB |
| I | 2 | 2 | 4 GB |
| J | 1 | 1 | 4 GB |

**Total v3.1: 65–135 CPU-hours** with full reciprocal anchoring (~3× faster than v2). 2–5 wall-hours on a 32-core node.

---

## Tests before first full run

1. Mini-multigenome (200 kb chr5 region, all 8 references). A–J in < 15 min.
2. Synthetic pir cluster (30 genes, 5 lost, 8 duplicated). Phase G recovers duplications.
3. Frame-shift fixture. Triage rule R1 fires; CESAR2 produces correct alignment.
4. ID-mapping fixture. Phase A harmonization produces consistent IDs.
5. MHC087 outlier fixture (mash distance 0.08). Phase A halts.
6. Liftoff-vs-CESAR2 agreement audit. 100 clean genes; ≥ 99% Liftoff-only, full agreement on forced CESAR2.
7. Triage rule coverage. Each fixture in `tests/triage/` fires exactly its named rule.

---

## Citations

- **Liftoff** — Shumate & Salzberg 2021, *Bioinformatics*
- TOGA2 / CESAR2 / chainCleaner / RepeatFiller — Hiller lab
- PGGB — Garrison et al. 2023
- MACSE v2 — Ranwez et al. 2018
- OrthoFinder methodology — Emms & Kelly 2019; v3 preprint 2025
- HyPhy / aBSREL / BUSTED / RELAX — Kosakovsky Pond lab
- PvP01 reference and pir family redefinition — Auburn et al. 2016; Lopez et al. 2013
- PAM reference — De Meulenaere et al. 2023
- PvSY56 reference — primary publication (operator to add)
- BUSCO `plasmodium_odb10` — Manni et al. 2021
- PlasmoDB-68 release — Amos et al. 2022
- pyfaidx — Shirley et al. 2015