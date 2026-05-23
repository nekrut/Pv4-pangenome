# Pv4 — *Plasmodium vivax* 8-strain pangenome

Reproducible artifacts for the Pv4 v3 analysis: pangenome graph, ortholog tables, codon MSAs, ML trees, HyPhy selection scans, cohort VCF projections, and a UCSC track hub.

## Contents

- **`v3/`** — main analysis: 8-strain pangenome (PGGB), 27 starred-essential outputs documented in `v3/writeup/OUTLINE.md`. See `v3/pipeline/LOCAL.md` for the executable recipe.
- **`v3/ucsc_hub/`** — UCSC track-hub manifests + small BBs. Large `.multiz.maf.bb` files (~4 GB total) are on Dropbox.
- **`v2/pggb_out/`** — 8-way PGGB graph (`*.smooth.fix.gfa.gz` only — 75 MB). `.og` (685 MB) and `.og.lay` (109 MB) on Dropbox.
- **`PlanA.md`, `PlanB.md`, `Pv4_samples.txt`, `ena_manifest.tsv`** — initial study design and sample inventory.

## Large files on Dropbox

342 GB of analysis intermediates + final products (cohort VCFs, AXTs, multi-way MAFs, soft-masked FASTAs, the PGGB `.og` graph, etc.). See `v3/writeup/LARGE_FILES_DROPBOX.md` for the full manifest + MD5 checksums + per-file recreate recipes.

Folder share URL: https://www.dropbox.com/scl/fo/gx1mta4adubja4bsxxmgm/AOni3YRX8TS1E-saUblX-eo?rlkey=0ksw8a5hhkxhy3sljqyycr7bu&dl=0

## Pipeline reproduction

- `v3/pipeline/LOCAL.md` — bash + container-based execution, 11 phases (A through K).
- `v3/pipeline/GALAXY.md` — Galaxy workflow port plan.
- `v3/pipeline/README.md` — overview + swap-to-new-species checklist.

A second-species validation on *P. knowlesi* is in `Pk/` (pipeline scaffold, 4,041 LoC, awaiting input fetch).

## BRC-analytics integration

This dataset is the first BRC pangenome resource — see https://github.com/galaxyproject/brc-analytics/issues/1279.
