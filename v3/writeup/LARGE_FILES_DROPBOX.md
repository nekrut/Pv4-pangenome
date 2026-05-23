# Large files for Dropbox upload

Pv4/v3 analysis intermediates and outputs too large for GitHub. Each entry: path, size, essentiality, and recreate recipe.

**Dropbox folder (public, world-readable)**: https://www.dropbox.com/scl/fo/gx1mta4adubja4bsxxmgm/AOni3YRX8TS1E-saUblX-eo?rlkey=0ksw8a5hhkxhy3sljqyycr7bu&dl=0

**Local staging** (cleanup mirror): `/media/anton/scratch/Pv4_dropbox_staging/`
**Total uploaded**: 339 GB across 369 files (+757 MB pggb_8way graph, added after main upload)
**MD5 manifest**: `writeup/LARGE_FILES_DROPBOX.tsv` (path, size, md5, essential, category)

`⭐` marks files that appear in the `*`-marked essentials list in OUTLINE.md — preserve carefully.

## Download

For a single file, append `?dl=1` to its Dropbox URL to force direct download (no preview interstitial). For the whole archive, use rclone:

```bash
rclone copy dropbox-public:Pv4_v3/ ./Pv4_v3/ --transfers 4
```

(where `dropbox-public` is a rclone remote configured against the folder share URL above; or use the [Dropbox web UI](https://www.dropbox.com/scl/fo/gx1mta4adubja4bsxxmgm/AOni3YRX8TS1E-saUblX-eo?rlkey=0ksw8a5hhkxhy3sljqyycr7bu&dl=0) "Download as ZIP" — ~339 GB).

## VCF cohorts

| ⭐ | Path glob | Files | Size | Recreate |
|---|---|---:|---:|---|
| ⭐ | `projection/A2_lastz/Pv4_cohort_on_GCA_*.vcf.gz` | 7 | 148 GB | `scripts/overnight/07b_a2_redo.sh 6` — needs chains (in git, gzipped) + softmasked FASTAs (Dropbox). ~3 hrs. |
|   | `projection/B_graph/Pv4_cohort_on_GCA_*.vcf.gz`  | 7 | 81 GB  | `scripts/overnight/10c_path_b_final.sh` — needs 8-way PGGB graph + MalariaGEN per-chr. ~6 hrs. |
|   | `projection/A1_wfmash/Pv4_cohort_on_GCA_*.vcf.gz` | 7 | 65 GB  | wfmash → paftools view → CrossMap → bcftools concat (driver: `scripts/run_a1_concat_fix.sh`). ~3 hrs. |

## MalariaGEN per-chromosome VCFs

| ⭐ | Path | Files | Size | Recreate |
|---|---|---:|---:|---|
| ⭐ | `projection/A1_wfmash/mg_renamed/Pv4_*.vcf.gz` | 16 | 24 GB | `bcftools annotate --rename-chrs inputs/annotations/PvP01_plasmodb_to_genbank.tsv` on `/media/anton/scratch/malariagen_pv4/Pv4_PvP01_*.vcf.gz`. |

## Pairwise alignments

| ⭐ | Path | Files | Size | Recreate |
|---|---|---:|---:|---|
| ⭐ | `projection/A2_kegalign/axt/{A}__vs__{B}.axt` | 56 (28 pairs × 2 directions) | 19 GB | `scripts/run_kegalign.sh` — needs 8 softmasked FASTAs + GPU node. ~6 hrs. |

## Multiz multi-way alignments

| ⭐ | Path glob | Files | Size | Recreate |
|---|---|---:|---:|---|
| ⭐ | `work/07_multiz/{hinge}/{hinge}.multiz.maf` | 8 | ~24 GB | `scripts/overnight/01b_multiz_repair.sh` — needs A2 AXTs. ~3 hrs per hinge. |
|   | `work/07_multiz/{hinge}/{hinge}_vs_{other}.maf` | 56 | ~14 GB | `axtToMaf` on A2 AXTs. Pairwise inputs to multiz. |

## Soft-masked assemblies

| ⭐ | Path | Files | Size | Recreate |
|---|---|---:|---:|---|
| ⭐ | `genomes/softmasked/{GCA}.fa` (+ `.fai`, `.mmi`) | 32 | 389 MB | `v2/mask_one.sh` per assembly (longdust ∪ sdust → bedtools maskfasta -soft). ~10 min. |

## Reference graphs

| ⭐ | Path | Size | Recreate |
|---|---|---:|---|
| ⭐ | 8-way PGGB graph (v2 path) | ~975 MB | `pggb -i pv_panSN.fa.gz -n 8 -s 5000 -p 90 -t 32`. ~6 hrs. Lives in v2/pggb_out/; v3 has pointer file in `inputs/pggb/README.md`. |
|   | `cactus_2way/pggb_out/pggb_2way.fa.gz.*.{og,gfa}` | ~300 MB | `pggb -i cactus_2way/pggb_in/pggb_2way.fa.gz -n 2 -s 5000 -p 90 -t 16`. ~30 min. |

## Tooling outputs

| ⭐ | Path | Size | Recreate |
|---|---|---:|---|
|   | `panagram_2way/index/` | 737 MB | See `panagram_2way/NOTES.md` (KMC k=21 + panagram --prepare + snakemake). ~2 min. |
|   | `writeup/browser/` | 917 MB | IGV.js + BAM tracks from chain/PAF/graph blocks. ~30 min. |

## Per-gene workdir archives (work_archives/)

These started as in-tree directories with thousands of small files. Tarred + gzipped to single archives. Archives smaller than 50 MB stay in git as `*_archive.tar.gz` files alongside the would-be directory; the rest staged here.

| ⭐ | Archive | Approx size (gz) | Contents |
|---|---|---:|---|
| ⭐ | `work/06_msa/core_v3_trees_archive.tar.gz` | ~60 MB | 1,584 IQ-TREE per-gene workdirs (treefile, iqtree, log, ckp.gz, bionj). |
| ⭐ | `work/06_msa/core_relaxed_trees_archive.tar.gz` | ~160 MB | 4,215 IQ-TREE per-gene workdirs (relaxed set). |
| ⭐ | `work/06_msa/core_v3_hyphy_archive.tar.gz` | ~95 MB | HyPhy bulk BUSTED 1,584 JSON results. |
| ⭐ | `work/06_msa/core_relaxed_hyphy_archive.tar.gz` | ~400 MB | HyPhy bulk BUSTED 4,215 JSON results. |
|   | `work/08_orthofinder_archive.tar.gz` | ~90 MB | OrthoFinder3 MCL output (orthogroups + sequences + counts). |

**Recreate**: re-run `scripts/build_per_gene_tree.py` / `scripts/run_hyphy_bulk.py` / `scripts/run_orthofinder3.sh` on `core_v3/` / `core_relaxed/` (also tar-gzipped, in git).

## Verification after upload

For ⭐ files, MD5 against `writeup/LARGE_FILES_DROPBOX.tsv` after copying to Dropbox:

```bash
cd /media/anton/scratch/Pv4_dropbox_staging
md5sum -c <(awk -v FS='\t' 'NR>1 && $4=="yes" {print $3"  "$1}' /media/anton/data/sandbox/Pv4/v3/writeup/LARGE_FILES_DROPBOX.tsv | sed 's|^[^/]*/||')
```

## Suggested Dropbox layout

```
Dropbox/Pv4_v3/
├── projection/
│   ├── A1_wfmash/        (cohort VCFs, 65 GB)
│   ├── A2_lastz/         (cohort VCFs, 148 GB)
│   ├── B_graph/          (cohort VCFs, 81 GB)
│   ├── mg_renamed/       (per-chr MalariaGEN, 24 GB)
│   └── A2_kegalign_axt/  (pairwise AXTs, 19 GB)
├── multiz/               (8 hinge MAFs + 56 pairwise, 38 GB)
├── genomes/softmasked/   (385 MB)
├── graphs/               (8-way + 2-way PGGB, 1.3 GB)
├── work_archives/        (per-gene tarballs, ~700 MB)
├── panagram_index/       (737 MB)
└── browser/              (IGV.js + BAMs, 917 MB)
```

Upload command (rclone):
```bash
rclone copy /media/anton/scratch/Pv4_dropbox_staging dropbox:Pv4_v3/ --progress --transfers 4
```
