# MalariaGEN cohort VCF projection (Path A2 KegAlign)

The MalariaGEN Pv4 cohort — 1,895 samples genotyped against PvP01 — is the population-level variation layer we want to display on every other strain assembly. Three projection paths were tested in v3 (A1 wfmash chains, A2 KegAlign chains, B graph-native). **A2 won the drug-resistance QC on all 19 genes** and is the canonical projection. This document covers A2 only.

**Related documents.** [MULTIZ.md](MULTIZ.md) covers the upstream KegAlign + chain pipeline that produces the cleaned chains we use for liftover (Steps 2–3 below depend on those chains); [PANGENOME.md → How the graph fed downstream analyses](PANGENOME.md#how-the-graph-fed-downstream-analyses) covers the non-canonical Path B graph-native branch (mentioned as a comparator in the cross-validation section here); [ORTHOLOGY.md](ORTHOLOGY.md) covers the consensus ortholog table that maps drug-resistance gene IDs between strains for the QC step.

## How we did it

### Source data

MalariaGEN Pv4, per-chromosome cohort VCFs, in PvP01 coordinates with PlasmoDB-style chromosome names (`PvP01_01_v1`…`PvP01_14_v1`, `PvP01_API_v1`, `PvP01_MIT_v1`):

```
/media/anton/scratch/malariagen_pv4/Pv4_PvP01_01_v1.vcf.gz
...
/media/anton/scratch/malariagen_pv4/Pv4_PvP01_14_v1.vcf.gz
/media/anton/scratch/malariagen_pv4/Pv4_PvP01_API_v1.vcf.gz
/media/anton/scratch/malariagen_pv4/Pv4_PvP01_MIT_v1.vcf.gz
```

16 files, ~24 GB total, 1,895 sample columns each.

### Step 1 — Chromosome rename to GenBank style

PvP01 assembly FASTA uses GenBank `LT635xxx.N` accessions; MalariaGEN uses PlasmoDB `PvP01_NN_v1` style. Rename so the chain pipeline's source coords match:

```
bcftools annotate \
  --rename-chrs inputs/annotations/PvP01_plasmodb_to_genbank.tsv \
  /scratch/malariagen_pv4/Pv4_PvP01_{NN}_v1.vcf.gz \
  -Oz -o projection/A1_wfmash/mg_renamed/Pv4_{NN}.vcf.gz
bcftools index projection/A1_wfmash/mg_renamed/Pv4_{NN}.vcf.gz
```

The rename map is a 2-column TSV (`PvP01_01_v1<TAB>LT635612.2`, etc.) at `inputs/annotations/PvP01_plasmodb_to_genbank.tsv`. Coverage must be total — if a contig in the source VCF is missing from the map, bcftools errors.

### Step 2 — CrossMap projection per chromosome per target

For each of the 7 non-PvP01 strain assemblies (`GCA_000002415.2` Sal-I through `GCA_040114635.1` MHC087), and each of the 16 chromosomes:

```
CrossMap vcf \
  work/01_chains/PvP01.{strain}.cleaned.chain \
  projection/A1_wfmash/mg_renamed/Pv4_{NN}.vcf.gz \
  genomes/softmasked/{GCA}.fa \
  /dev/stdout |
  bcftools sort -T /media/anton/scratch/path_a2_sort_tmp -Oz \
    -o projection/A2_lastz/per_chr/{GCA}/Pv4_{NN}_on_{GCA}.vcf.gz
bcftools index projection/A2_lastz/per_chr/{GCA}/Pv4_{NN}_on_{GCA}.vcf.gz
```

Three things to get right:

1. **`bcftools sort` is not optional**. CrossMap interleaves chromosomes when chain segments cross contigs. Without sort the index step fails.
2. **The sort temp dir must be on fast disk** (NVMe scratch, ≥ 100 GB free). Default `/tmp` (often tmpfs) runs out of RAM on 25 GB sorts.
3. **The chain MUST be the cleaned chain** built by the KegAlign + chain pipeline ([MULTIZ.md → How we built them](MULTIZ.md#how-we-built-them), "Chain pipeline" subsection), not the AXT directly, not the rbest chain. The cleaned chain is the canonical UCSC chain build (axtChain → chainSort → chainPreNet → chainNet → netChainSubset → chainStitchId). rbest is overly conservative for variant projection — drops ~40 % of liftable variants.

### Step 3 — Concat per-chr → cohort

```
bcftools concat -a projection/A2_lastz/per_chr/{GCA}/*.vcf.gz \
  -Oz -o projection/A2_lastz/Pv4_cohort_on_{GCA}.vcf.gz
bcftools index projection/A2_lastz/Pv4_cohort_on_{GCA}.vcf.gz
```

`-a` allows duplicate records (some inversion-crossing chain segments can lift one source variant to two target positions; we keep both).

Wall time: ~30 min per target × 7 targets = ~3.5 hours sequential, ~1 hour at `xargs -P 4`.

The driver script is `scripts/overnight/07b_a2_redo.sh`.

## Outputs

| File                                                                                                   |         Count | Size each    | Total  | Where                                  |
| ------------------------------------------------------------------------------------------------------ | ------------: | ------------ | ------ | -------------------------------------- |
| `projection/A1_wfmash/mg_renamed/Pv4_{NN}.vcf.gz` (+ `.csi`) — renamed MalariaGEN source               |            16 | ~1.5 GB      | 24 GB  | **Dropbox** — `Pv4_v3/mg_renamed/`     |
| `projection/A2_lastz/per_chr/{GCA}/Pv4_{NN}_on_{GCA}.vcf.gz` — per-chr intermediates                   |           112 | ~200 MB–2 GB | ~80 GB | Not preserved (deleted after concat)   |
| `projection/A2_lastz/Pv4_cohort_on_{GCA}.vcf.gz` (+ `.csi`) — **the canonical 7 cohort VCFs**          |             7 | 4.5–23 GB    | 148 GB | **Dropbox** — `Pv4_v3/A2_lastz/`       |
| `projection/A2_lastz/chain/` — copies of the chains used                                               | 8 small files | <1 MB        | 8 MB   | **Dropbox** — `Pv4_v3/A2_lastz_chain/` |
| `work/09_projection_compare/{GCA}/A2.sites.tsv.gz` — per-target site lists for cross-method comparison |             7 | 10–13 MB     | 80 MB  | **Git**                                |

The chains themselves (`work/01_chains/PvP01.{strain}.cleaned.chain.gz`) live in Git and are documented in `MULTIZ.md`. The cohort VCFs and the renamed source VCFs are the heavy outputs — both on Dropbox.

## How these projected VCFs fed downstream analyses

**Drug-resistance QC across strains.** The lifted-cohort VCFs were the input to the per-gene drug-resistance check on 19 known *P. vivax* resistance loci (dhfr 4 codons, dhps 5 SNPs, mdr1 3 SNPs, plus pvcrt and pvmrp1 markers). For each of the 19 genes, on each of the 7 non-PvP01 references, we ran a site-coverage + allele-concordance score. **A2 won all 19**: median per-gene score 0.72, vs Path A1 wfmash at 0.13 and Path B graph-native at 0.45. Driver script: `scripts/overnight/09_phase5_drug_resistance_qc.py`. Output table at `work/09_projection_compare/drug_resistance_qc.tsv` (138 rows = 19 genes × 7 targets + summary). The per-strain gene IDs come from the consensus ortholog table — see [ORTHOLOGY.md → Outputs](ORTHOLOGY.md#outputs).

**3-way intersection (A1 ∩ A2 ∩ B) cross-validation.** All three paths produced the same downstream outputs, intersected per-site for cross-validation. The per-target site TSVs at `work/09_projection_compare/{GCA}/{A1,A2,B}.sites.tsv.gz` (21 files) feed the intersection table at `work/09_projection_compare/intersection.tsv`. A2 sites had the highest agreement with the canonical PvP01 sites (typical recovery 92–98 %). Path B (graph-native) details: [PANGENOME.md → How the graph fed downstream analyses](PANGENOME.md#how-the-graph-fed-downstream-analyses), final paragraph.

**Future — BRC population-genomic tracks.** The A2 cohort VCFs are the input for computing iHS, FST, π, Tajima's D, and per-country allele frequencies on each non-PvP01 reference. These would land as BigWig and BigBed tracks under the `brc_pangenome_popgen_pv_v1` composite in the UCSC hub. The brc-analytics issue (#1279) lays out the data model. v3 didn't compute these — flagged for a follow-up workflow.

The 8th projection (PvP01 onto PvP01) is the source itself — the 16 `mg_renamed/` per-chr files are the canonical-coords cohort, just chrname-translated.

## Re-running on a different cohort

If you have a per-chromosome multi-sample cohort VCF in the reference's coords and a target set of related assemblies, the recipe carries straight over. Walk-through (assumes soft-masking and the KegAlign + chain pipeline are done — see [MULTIZ.md → Re-running on a different species](MULTIZ.md#re-running-on-a-different-species)):

1. Drop the per-chromosome source VCFs into `inputs/cohort_vcf/` and set `COHORT_VCF_DIR` + `COHORT_CHROM_GLOB` in `pipeline/species.conf`.
2. Build the chromosome rename map: `inputs/annotations/{ref}_plasmodb_to_genbank.tsv` (or any equivalent — 2-col TSV mapping source chrnames to assembly-FASTA chrnames). The helper `pipeline/lib/fix_gff_chroms.sh` builds one automatically from the GFF + FASTA header diff if your annotation and assembly use different naming.
3. `bash pipeline/11_project_vcf.sh` — this runs Steps 1, 2, 3 above for every target. Wall: ~30 min per target on a 32-core box with NVMe scratch.

Three parameters worth tuning by scale:

- For huge cohorts (≥ 10,000 samples), the per-chr intermediates balloon. Use `bcftools view -G` to drop sample columns before projection if you only need site-level data, then rejoin samples after.
- For inter-species cohorts, drop `-linearGap=loose` → `-linearGap=medium` in the chain build (Phase C) so CrossMap doesn't try to lift over fragmented chain segments.
- The sort temp dir size scales linearly with per-chr VCF size. For Pv4-scale data, allow 30 GB scratch per target running concurrently.

The *P. knowlesi* scaffold at `Pk/v1/pipeline/11_project_vcf.sh` is the parameterized version of this. A `COHORT_VCF_DIR=TODO_NO_COHORT_VCF_AVAILABLE` placeholder is in the Pk species.conf because no community Pk cohort exists yet — the recipe is ready when one does.

## Galaxy tool-wrap list

Path A2 is the simplest of the three paths to surface in Galaxy because all the tools already exist. No new wrappers needed; the work is workflow assembly.

| Tool                                           | Galaxy state | Note                                                             |
| ---------------------------------------------- | ------------ | ---------------------------------------------------------------- |
| **bcftools annotate** (`--rename-chrs`)        | ✅ exists     | IUC; supports the rename-chrs flag                               |
| **CrossMap vcf**                               | ✅ exists     | IUC; the relevant subcommand is `crossmap_vcf`                   |
| **bcftools sort**                              | ✅ exists     | Make sure the wrapper exposes `-T` (temp-dir) for large inputs   |
| **bcftools concat**                            | ✅ exists     | IUC; need `-a` flag for allowing duplicates                      |
| **bcftools index** (CSI)                       | ✅ exists     |                                                                  |
| **Per-target per-chr collection split + join** | ✅ built-in   | Galaxy's collection rules handle the `{target} × {chrom}` matrix |

Workflow shape: take a `list:list` collection (outer = target assembly, inner = per-chr VCF), iterate `crossmap_vcf` over it, then collapse the inner collection with `bcftools concat`. The chain files come in as a parallel `list` collection (one chain per target). The chromosome rename map is a single auxiliary input.

Suggested workflow name: `malariagen-cohort-project-via-chain` — generic enough to work on any cohort + any chain-based liftover. Category: `POPULATION_GENOMICS` (per the brc-analytics#1279 PR plan).

The 3-way cross-validation step (A1 ∩ A2 ∩ B) is a separate add-on workflow if it's wanted; the per-method site TSVs are produced by a small `bcftools query -f '%CHROM\t%POS\n'` per cohort VCF, which is trivially wrappable.
