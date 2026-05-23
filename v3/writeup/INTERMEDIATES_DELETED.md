# Intermediate directories deleted during repo cleanup

These directories existed during the analysis but were removed before git commit because they:
1. Are intermediate state superseded by downstream outputs
2. Are easily regenerable from the preserved scripts
3. Contain non-project tool data (BUSCO lineages, Python venvs, etc.)

Below: each deleted path with its recreation recipe.

## Phase C — Liftoff/TOGA2 intermediates

### `work/02a_liftoff/` (1.9 GB)
Raw per-strain Liftoff projections before triage.
**Superseded by**: `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` (preserved as `*_archive.tar.gz`).
**Recreate**: `scripts/run_phase_c1_liftoff.sh` — ~2 hrs per anchor.

### `work/02b_triage/` (172 MB)
Triage 8-rule classifications pre-CESAR2 rescue.
**Superseded by**: merged GFF3 + `*.classification.tsv` in 02d_merged archives.
**Recreate**: `scripts/phase_c2_triage.py`. Inputs: 02a_liftoff outputs + per-strain GFF3s in `inputs/annotations/`.
**Note**: Removal blocked by root-owned files (docker artifacts). See `writeup/NEEDS_SUDO_CLEANUP.md`.

### `work/02c_toga/*-as-ref/` (398 MB combined)
TOGA2/CESAR2 fallback workdirs (BED12, intermediate logs).
**Superseded by**: merged outputs.
**Recreate**: `scripts/run_phase_c3_toga2.sh` per anchor — ~6 hrs.
**Note**: Removal partially blocked by root-owned files. See `writeup/NEEDS_SUDO_CLEANUP.md`.

## Phase E/F failed attempts

### `work/08_genespace/` (207 MB)
GENESPACE failed runs (4 header-debugging attempts).
**Reason**: `parse_annotations` matched 0 genes from non-standard GFFs (filed jtlovell/GENESPACE#206).
**Substituted by**: OrthoFinder3 in `work/08_orthofinder/` (preserved as `_archive.tar.gz`).
**Note**: Removal partially blocked by root-owned files. See `writeup/NEEDS_SUDO_CLEANUP.md`.

## Phase G MSA prototypes

### `work/06_msa/8way_{strict,relaxed}{,_v2}/` (~110 MB)
Earlier MSA prototypes (4 sub-dirs).
**Superseded by**: `core_v3/` (strict, min_intact=7) and `core_relaxed/` (min_intact=5), preserved as `*_archive.tar.gz`.
**Recreate**: `scripts/build_8way_msa_v2.py` with appropriate `--min_intact` flag.

## Tool downloads + venvs

### `work/00_inventory/busco/busco_downloads/` (~1.4 GB)
BUSCO lineage HMM downloads (plasmodium_odb10).
**Recreate**: BUSCO redownloads automatically on next run with `--download` flag.

### `panagram_2way/{venv,panagram}/` (~1.8 GB)
Python venv + git clone of panagram source.
**Recreate**:
```bash
cd panagram_2way
python3 -m venv venv && source venv/bin/activate
git clone --recursive https://github.com/kjenike/panagram
cd panagram && pip install .
```

## TOGA2 internal workdirs

### `intermediate_files/`, `TOGA2_ref_annotation_*/`, `toga2_run_*/` (~70 MB combined)
TOGA2 leftover scratch dirs.
**Recreate**: re-run TOGA2 — these are regenerated automatically.
**Note**: Removal blocked by root-owned files. See `writeup/NEEDS_SUDO_CLEANUP.md`.

## Path B intermediates

### `projection/B_graph/cohorts/{GCA_*}/Pv4_*.vcf.gz` (~82 GB total)
Per-chromosome cohort intermediates (16 chromosomes × 7 targets = 112 files).
**Superseded by**: `Pv4_cohort_on_{GCA}.vcf.gz` (concatenated; moved to Dropbox).
**Recreate**: `scripts/overnight/10c_path_b_final.sh`.

### `projection/B_graph/pangenome_on_{GCA}.vcf` (~2.5 GB total)
Raw vg deconstruct output before site-filtering.
**Recreate**: `vg deconstruct -P {GCA}_path v2/vg_idx/pv.gbz`.

### `projection/B_graph/sites/mg_on_{GCA}.tsv` (~9.5 GB total)
odgi position lookups for MalariaGEN sites.
**Recreate**: `odgi position -i v2/pggb_out/...og -b mg_sites.bed -r {GCA}_path > sites/mg_on_{GCA}.tsv`.

## Logs

### `logs/overnight/*.log` (50+ files)
Driver script run logs.
**Recreate**: re-run the corresponding script.


## Newer deletions (this cleanup pass)

### `projection/A1_wfmash/synteny_rerun/` (~365 MB → staged)
Earlier exploratory wfmash run (tuned axt + PAF + cleaned chain).
**Reason for staging**: superseded by canonical A2 KegAlign chains in `work/01_chains/`.

### `projection/A1_wfmash/test_MIT_on_Sal1.vcf` (~14 MB)
One-off chrM test projection.
**Recreate**: re-run A1 pipeline on just MalariaGEN MIT VCF.

### `projection/A1_wfmash/MalariaGEN_MIT_renamed.vcf.gz`, `MalariaGEN_chr05_renamed.vcf.gz`
Single-chr early-prototype renamed VCFs (precursor to mg_renamed/).
**chr05 (1.6 GB)**: staged to Dropbox alongside mg_renamed/.
**MIT (9 MB)**: deleted (recreatable from `Pv4_PvP01_API_v1.vcf.gz` + rename-chrs).

### `projection/A2_kegalign/work/` (~328 MB)
KegAlign intermediate workdir (per-pair scratch).
**Recreate**: re-run `scripts/run_kegalign.sh`.

### `projection/A2_kegalign/2bit/`, `chain/`, `lifted/`, `work2/` (small + empty)
KegAlign intermediate format conversions.
**Recreate**: `faToTwoBit` on each genomes/softmasked/*.fa.

### `projection/A2_lastz/chain/` (~8 MB)
Duplicate chain files (canonical copies live in `work/01_chains/`).
**Recreate**: copy from `work/01_chains/*.cleaned.chain.gz` + gunzip.

### `projection/B_graph/cohorts/{GCA_*}/Pv4_*.vcf.gz` (~82 GB)
Per-chr per-target Path B intermediates (16 chr × 7 targets).
**Recreate**: `scripts/overnight/10c_path_b_final.sh` — the per-chr build phase.

### `projection/B_graph/pangenome_on_{GCA}.vcf` (~2.5 GB)
Raw vg deconstruct output.
**Recreate**: `vg deconstruct -P {GCA}_path v2/vg_idx/pv.gbz`.

### `projection/B_graph/sites/mg_on_{GCA}.tsv` (~1.3 GB)
odgi position lookups for MalariaGEN sites.
**Recreate**: `odgi position -i v2/pggb_out/...og -b mg_sites.bed -r {GCA}_path > mg_on_{GCA}.tsv`.

### `projection/B_graph/ref_gbz/{GCA}.gbz` (~614 MB → staged)
Per-target GBZ subgraphs.
**Recreate**: `vg convert -g v2/pggb_out/...gfa -P {GCA} > {GCA}.gbz`.

### `inputs/assemblies/*.fa.mmi` (~750 MB)
minimap2 indexes for the 8 assemblies.
**Recreate**: `minimap2 -d {GCA}.fa.mmi {GCA}.fa` (instant).

### `work/00_inventory/busco/busco_downloads/` (~1.4 GB)
BUSCO lineage HMM downloads (plasmodium_odb10).
**Recreate**: BUSCO auto-downloads on next run.

### `cactus_2way/{chr1,chr1_only,input}/` (~261 MB → staged)
2-way PvP01+PAM PGGB outputs for chr1.
**Recreate**: `pggb -i cactus_2way/pggb_in/pggb_2way.fa.gz` (~30 min).

