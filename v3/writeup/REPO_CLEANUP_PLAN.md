# Repo cleanup plan — prepare v3/ for GitHub commit

GitHub limits: 100 MB hard / 50 MB warning per file; recommended <1 GB repo total.
v3/ is currently 475 GB. Need to triage into 3 buckets: **KEEP** (commit, gzip if text), **DROPBOX** (move out, list externally), **DELETE** (intermediate, recreate from scripts).

## Summary of triage decisions

| Bucket | Files | Size | Action |
|:---|---:|---:|:---|
| **DROPBOX-HARD** (>100 MB each) | 249 | 460 GB | Move to `/media/anton/scratch/Pv4_dropbox_staging/` for upload |
| **DROPBOX-WARN** (50–100 MB) | 54 | 4.4 GB | Gzip first; commit if <50 MB after, otherwise move to staging |
| **GZIP** (5–50 MB) | 318 | 4.9 GB | Gzip in place; commit gzipped version |
| **<5 MB** | thousands | <500 MB | Commit as-is (scripts, configs, small TSVs) |

## Dropbox-bound large files

By directory, with recreation recipe:

### `projection/A1_wfmash/Pv4_cohort_on_*.vcf.gz` (~92 GB, 7 files)
A1 cohort VCFs from wfmash chains + CrossMap.
**Recreate**: re-run wfmash → paftools.js → CrossMap pipeline (~3 hrs). Driver script: `scripts/overnight/...` (none preserved as a single driver — was run interactively).

### `projection/A2_lastz/Pv4_cohort_on_*.vcf.gz` (148 GB, 7 files)
A2 cohort VCFs from KegAlign chains + CrossMap + bcftools sort.
**Recreate**: `scripts/overnight/07b_a2_redo.sh` 6 (parallel) — needs KegAlign chains under `work/01_chains/` (kept in git, gzipped) + `genomes/softmasked/` (delete from git, recreate via `v2/mask_one.sh`).

### `projection/B_graph/` (164 GB total, 112 large files)
- `cohorts/{GCA}/Pv4_*.vcf.gz` (~82 GB, 112 files) — per-chr intermediates, **DELETE** (regenerable)
- `Pv4_cohort_on_{GCA}.vcf.gz` (~82 GB, 7 files) — final cohort VCFs, **DROPBOX**
- `pangenome_on_{GCA}.vcf` (~2.5 GB, 7 files) — vg deconstruct output, **DROPBOX**
- `sites/mg_on_{GCA}.tsv` (1.3 GB) — odgi position lookups, **GZIP+DROPBOX** (still big after gzip)
**Recreate**: `scripts/overnight/10c_path_b_final.sh` — needs 8-way PGGB graph (symlinked to v2) + MalariaGEN per-chr VCFs (kept in mg_renamed/, gzipped — move to dropbox).

### `projection/A2_kegalign/axt/` (20 GB, 28 files + others)
Raw KegAlign AXT alignment outputs.
**Recreate**: GPU KegAlign run via `scripts/run_kegalign.sh` on the 8 softmasked assemblies. Wall time ~6 hrs on the GPU node.

### `projection/A1_wfmash/mg_renamed/Pv4_{01..14,API,MIT}.vcf.gz` (~25 GB)
MalariaGEN Pv4 with CHROMs renamed from `PvP01_NN_v1` → `LT635xxx`.
**Recreate**: `bcftools annotate --rename-chrs inputs/annotations/PvP01_plasmodb_to_genbank.tsv` on the source `/media/anton/scratch/malariagen_pv4/Pv4_PvP01_*.vcf.gz` (preserved on scratch).

### `work/07_multiz/*/{hinge}.multiz.maf` and `{hinge}_vs_{other}.maf` (37 GB, 52 large files)
- 8 final `{hinge}.multiz.maf` files (~3 GB each, **DROPBOX**)
- 56 pairwise `{hinge}_vs_{other}.maf` intermediates (~600 MB–1.5 GB each, **DELETE**)
**Recreate**: `scripts/overnight/01_multiz_all.sh` — needs A2 AXTs (also in dropbox).

### `inputs/pggb/pv.{og,gfa}` (1 GB total, symlinks to v2)
Symlinks — replace with text files containing the v2 source path. Actual graph stays in v2/.

### `panagram_2way/index/` (737 MB)
Panagram k-mer index. **DROPBOX**.
**Recreate**: see `panagram_2way/NOTES.md` for full recipe.

### `cactus_2way/pggb_out/pggb_2way.fa.gz.*.smooth.fix.{og,gfa}` (~140 MB)
2-way PvP01+PAM PGGB graph. **DROPBOX**.
**Recreate**: `pggb -i cactus_2way/pggb_in/pggb_2way.fa.gz -n 2 -s 5000 -p 90 -t 16 -o cactus_2way/pggb_out/`.

### `writeup/browser/` (~800 MB)
IGV.js + BAM tracks. **DROPBOX**.

## Intermediate deletions (and recreate recipes)

### `work/02a_liftoff/` (1.9 GB)
Raw per-strain Liftoff projections (pre-triage, pre-CESAR2 merge).
**Recreate**: `scripts/run_phase_c1_liftoff.sh` (~2 hrs per anchor).
**Why delete**: superseded by merged `work/02d_merged/{anchor}-as-ref/*.annotation.gff3`.

### `work/02b_triage/` (172 MB)
Triage 8-rule classifications, pre-CESAR2 rescue.
**Recreate**: `scripts/phase_c2_triage.py`.
**Why delete**: superseded by merged outputs.

### `work/02c_toga/*-as-ref/` (398 MB)
TOGA2/CESAR2 fallback workdirs (BED12, projections, intermediate logs).
**Recreate**: `scripts/run_phase_c3_toga2.sh` per anchor. ~6 hrs per anchor.
**Why delete**: superseded by merged outputs; only `classification.tsv` + `annotation.gff3` retained in `work/02d_merged/`.

### `work/06_msa/8way_{strict,relaxed}{,_v2}/` (~100 MB)
Earlier MSA prototypes.
**Recreate**: `scripts/build_8way_msa_v2.py` with appropriate flags.
**Why delete**: superseded by `core_v3/` and `core_relaxed/`.

### `work/08_genespace/` (207 MB)
GENESPACE failed attempts (4 staging directories from header debugging).
**Recreate**: skip — substituted with OrthoFinder3 (see `work/08_orthofinder/`).
**Why delete**: not used.

### `intermediate_files/`, `TOGA2_ref_annotation_*/`, `toga2_run_*/`
Explicit intermediate dirs left over from runs.
**Recreate**: re-run TOGA2.

### `/media/anton/scratch/A2_lastz_v2/` (270 GB)
A2 CrossMap per-chr intermediates on scratch.
**Recreate**: `scripts/overnight/07b_a2_redo.sh`.
**Why delete**: outputs were concatenated into final cohort VCFs; per-chr inputs no longer needed.

### `panagram_2way/{venv,panagram}/` (~1 GB)
Python venv + git clone of panagram source.
**Recreate**: `cd panagram_2way && python3 -m venv venv && source venv/bin/activate && git clone --recursive https://github.com/kjenike/panagram && cd panagram && pip install .`
**Why delete**: standard pip install, not project-specific.

### Logs
`logs/overnight/{relaxed_iqtree,hyphy_priority,hyphy_bulk,va_hyphy,path_b_final}_*.log` (many small)
**Recreate**: re-run respective driver scripts; logs are by-products.
**Why delete**: not load-bearing.

## What stays in git after triage

| Path | Why kept |
|:---|:---|
| `scripts/` | All driver scripts (recreate everything) |
| `writeup/` (sans `browser/`) | Reports, blog posts, OUTLINE, plans, plots |
| `inputs/annotations/` | Source GFFs / BEDs / proteomes (gzip large ones) |
| `inputs/pggb/pv.{og,gfa}` | Replace symlinks with text pointers to v2 path |
| `work/00_inventory/` | Mash + BUSCO + AGAT outputs (small) |
| `work/01_chains/*.{cleaned,rbest}.chain` | Gzipped chain files (small enough) |
| `work/02d_merged/{anchor}-as-ref/*.{annotation.gff3,classification.tsv}` | Merged annotation outputs (the load-bearing C/D outputs) |
| `work/03_consensus/ortholog_table.tsv` | Phase E consensus orthogroup table |
| `work/04_graph_validation/` | Phase F PAV table + BED |
| `work/05_families/` | Phase G family table |
| `work/05_priorities/` | Phase J priority gene list |
| `work/06_msa/{core_v3,core_relaxed,core_v3_clean,core_relaxed_clean,core_v3_trees,core_relaxed_trees}/` | MSAs + trees (gzip individual files) |
| `work/06_msa/core_v3_hyphy/`, `core_relaxed_hyphy/`, `core_relaxed_hyphy_va/` | HyPhy JSON results |
| `work/08_orthofinder/results/Results_May21/Orthogroups/` | Orthogroups table (gzip) |
| `work/09_projection_compare/` | Intersection table + per-method site TSVs (gzip site TSVs) |
| `tools/` | Custom build helpers |
| Top-level `.md` files | OUTLINE, plans, reports |
| `.gitignore` | Excludes the dropbox + scratch + intermediate dirs |

## Execution order

1. **Build manifests** (this file + `LARGE_FILES_DROPBOX.tsv` + `RECREATE_INTERMEDIATES.md`) — done
2. **Gzip text/data files** in the KEEP set — non-destructive
3. **Stage dropbox files** to `/media/anton/scratch/Pv4_dropbox_staging/` — `mv`, not `cp`
4. **Delete intermediates** — destructive, do last
5. **Write `.gitignore`** to prevent re-tracking the moved/deleted paths
6. **Final size check** — repo should be <1 GB total
