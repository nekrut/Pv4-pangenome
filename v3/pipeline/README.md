# Pangenome + selection-scan pipeline — overview

Reproduces the 27 essential outputs (`*`-marked in `writeup/OUTLINE.md`) from the Pv4 v3 *P. vivax* analysis. Designed to run on any apicomplexan-scale species panel (5–15 haploid assemblies, ~25 Mb each, gene-dense). Two execution plans:

- **`LOCAL.md`** — Claude Code / bash driver scripts on a single workstation (GPU optional). What v3 actually used.
- **`GALAXY.md`** — Galaxy workflow (`.ga`) with existing + to-be-wrapped tools. For multi-user / cloud / reproducibility.

## Inputs (assumed identical structure to v3)

| Input | Format | Count | Where it comes from |
|---|---|---|---|
| Haploid genome assemblies | `.fa` (one contig per chromosome ideal) | N (8 in v3) | NCBI Datasets, PlasmoDB, EuPathDB |
| Per-strain annotation | `.gff3` (PlasmoDB or NCBI style, gene/mRNA/CDS) | N | Same as above |
| Per-strain protein FASTA | `.fa` (one record per gene) | N | PlasmoDB AnnotatedProteins, or `gffread -y` |
| Variant cohort VCF | per-chromosome `.vcf.gz` (multi-sample) | 1 set | MalariaGEN / Pf3K / equivalent population panel |
| Reference strain choice | string (one of the N strains) | 1 | The strain whose coordinates the cohort VCF uses |

## Outputs (the 27 essentials)

Files produced by the pipeline, grouped by phase:

| Phase | Tool(s) | Output | Files |
|---|---|---|---|
| A. Inventory | mash, BUSCO | `00_inventory/mash/` | NxN matrix |
| B. Mask | longdust + sdust + bedtools | `genomes/softmasked/{S}.fa` | N |
| C. Pairwise align | KegAlign (GPU) | `A2_kegalign/axt/{S1}__vs__{S2}.axt` | N×(N-1)/2 pairs × 2 dirs |
| C. Chain pipeline | axtChain → chainSort → chainPreNet → chainNet → netChainSubset → chainStitchId; chainSwap+chainNet for rbest | `work/01_chains/{src}.{tgt}.cleaned.chain` + `.rbest.chain` | N×(N-1)/2 × 2 dirs |
| C. Annotation projection | Liftoff + TOGA2/CESAR2 merge | `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` + `.classification.tsv` | (N-1) per anchor × 4 anchors |
| D. PGGB graph | wfmash + seqwish + smoothxg + odgi | `inputs/pggb/pv.{og,gfa}` | 1 |
| E. Consensus orthology | Phase E union-find on chains + annotations | `work/03_consensus/ortholog_table.tsv` | 1 (~7,500 rows) |
| F. MSAs | gffread + MAFFT-LINSI + pal2nal | `work/06_msa/core_v3/{gene}.{codon,protein}.aln.fa` (strict) + `core_relaxed/` (relaxed) | 2 × (~1,500 + ~4,200) |
| G. ML trees | IQ-TREE3 `-m MFP -B 1000` | `work/06_msa/{set}_trees/{gene}/{gene}.treefile` | ~5,800 |
| H. HyPhy bulk | BUSTED (per gene) | `work/06_msa/{set}_hyphy/{gene}/busted.json` | ~5,800 |
| I. Multiz | axtToMaf + multiz | `work/07_multiz/{hinge}/{hinge}.multiz.maf` | N (one per strain as hinge) |
| J. Cohort projection | CrossMap + bcftools sort/concat | `projection/A2_lastz/Pv4_cohort_on_{S}.vcf.gz` | N-1 |
| J. Renamed cohort | bcftools annotate --rename-chrs | `projection/A1_wfmash/mg_renamed/Pv4_{chr}.vcf.gz` | per-chr (16 in v3) |

## Quick adapt-to-new-species checklist

1. Drop the new assemblies + annotations into `inputs/assemblies/` + `inputs/annotations/`
2. Edit one config file: `pipeline/species.conf` (strain list, reference strain, cohort VCF path, chromosome rename map)
3. Run `pipeline/run_all.sh` (local) or import `pipeline/galaxy_workflow.ga` (Galaxy)

Wall-time estimate on a 32-core / 1 GPU workstation, gene-dense ~25 Mb genomes: **~24 hours end-to-end**.

## What's NOT in the pipeline (and should be added before reuse)

- **GARD recombination pre-screen** before HyPhy BUSTED — top hits in v3 weren't GARD-filtered
- **BUSTED-MH on bulk** — single-hit BUSTED inflated by MNS events (see Pvs230 case in v3)
- **Outgroup injection** for rooted selection — v3 has no outgroup; *P. cynomolgi* / *P. knowlesi* needed for rooted dN/dS
- **MACSE in addition to MAFFT** — handles frameshifts + pseudogenes (relevant for variant antigen families)

These were on the original plan but cut for v3 to ship. The `*`-marked outputs do not depend on them.

## Phase dependency graph

```
       A.inventory
            |
       B.mask ──────┐
            |       |
       C.align (GPU)|
            |       |
       C.chain      |
       /  |  \      |
  C.proj  D.graph   |
       \  /         |
       E.ortho      |
            |       |
   F.MSA(strict + relaxed) ── G.trees ── H.HyPhy
                  |
            I.multiz (uses C.align AXTs)
                  
       C.chain ── J.cohortVCF (CrossMap)
```

A → B → {C, D in parallel} → E → F → {G, H in parallel}; I depends only on C; J depends on C.chain + variant cohort. Phases F-J can be sharded.

## Cost estimate per species

- Compute: 24 hrs × 32 cores ≈ 770 CPU-hrs ≈ $30 spot / $80 on-demand on AWS
- GPU: 6 hrs × A10G/A100 for KegAlign ≈ $5–15
- Storage: ~50 GB output (gzipped), peak ~450 GB intermediates
- Manual labor: 30 min config + ~2 hrs QC at the end

## Test-on-new-species smoke test

`pipeline/smoke_test.sh` runs all 11 phases on a **single chromosome** of 3 strains. Wall time ~30 min on the same hardware. Validates that the config + tooling are correct before committing to the full N-strain run.

See `LOCAL.md` and `GALAXY.md` for execution details.
