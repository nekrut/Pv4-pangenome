# Pv4 test panel — end-to-end pipeline test dataset

A curated subset of the full Pv4 v3 inputs that exercises **every pipeline phase (A–K)** end-to-end in **~2–3 h** — between the 30-min smoke test (3 strains, 1 chromosome, skips multiz) and the ~24-h / 339-GB full run. It is real data, regenerable from the script in this directory. See [`../PIPELINE_EXPLANATION.md`](../PIPELINE_EXPLANATION.md) for what each phase does.

## What's in it

**5 strains × 3 chromosomes** (+ optional organelles):

| | strains | role |
|---|---|---|
| anchors | PvP01 (ref), PvW1, PAM | PlasmoDB-curated; projected as annotation sources |
| non-anchors | PvT01, MHC087 | NCBI-annotated; contribute sequence to alignment/graph/orthology |

| PvP01 chr | accession | size | why it's here |
|---|---|---|---|
| chr04 | LT635615.1 | ~1.01 Mb | **Pvs230** (PVP01_0415800) — BUSTED / BUSTED-MH selection example |
| chr05 | LT635616.2 | ~1.52 Mb | **dhfr-ts** (PVP01_0526600, pyrimethamine) + AARP (selection) |
| chr14 | LT635625.2 | ~3.15 Mb | **dhps** (PVP01_1429500, sulfadoxine) — drug-resistance QC + hub default position |

Subtelomeres of all three carry the PIR/PHIST/Pv-fam multigene families that stress Phase C.4 (TOGA rescue), Phase E (orthology disagreement), and Phase F (divergent-paralog MSAs). Phase I (multiz) runs here — the smoke test skips it.

The cohort VCF is the MalariaGEN Pv4 cohort downsampled to **~200 samples** on these three chromosomes (Phase J input + the drug-resistance QC bridge).

## Files

Committed to git (small):

- `contig_map.tsv` — per-strain contig orthologs for chr04/05/14, derived by parsing the in-repo PvP01-as-target cleaned chains; chr14 rows match the hub `defaultPos` work exactly. (PAM merges chr04+chr05 onto one contig — the builder dedupes.)
- `species.conf` — the test config (5 strains, ref PvP01, anchors PvP01/PvW1/PAM, MIN_INTACT 4/3).
- `make_test_data.sh` — regenerates the bundle.
- `samples_200.txt` — the cohort sample subset (written by the VCF stage).
- this `README.md`.

Materialized bundle (heavy, **not** in git — on Dropbox under `Pv4_v3/test_data/`):

- `inputs/assemblies/{S}.fa` (+ `.fai`) — ~30 MB total
- `inputs/annotations/{S}.gff3` + `cohort_chrom_rename.tsv`
- `inputs/proteomes/{S}.proteins.fa`
- `inputs/cohort_vcf/Pv4test_{04,05,14}_v1.vcf.gz` (+ `.tbi`) — ~0.5–1 GB

## How to build

Host needs only `docker` + `bash` (+ `rclone` to upload). All bio-tools run in pinned biocontainers.

```bash
# sequence parts (assemblies, annotations, proteomes, chrom map) — local data only
bash v3/pipeline/make_test_data.sh --stage seq

# cohort VCF — needs the 3 PvP01-coord MalariaGEN per-chr files
bash v3/pipeline/make_test_data.sh --stage vcf \
     --cohort-src /path/to/malariagen_pv4        # holds Pv4_PvP01_{04,05,14}_v1.vcf.gz

# upload the heavy bundle
rclone copy v3/pipeline/test_data dropbox:Pv4_v3/test_data/ --transfers 4
```

## How to run the pipeline on it

```bash
cp v3/pipeline/test_data/species.conf $WORK/pipeline/species.conf   # WORK from species.conf
# stage the built inputs into $WORK/inputs/ (or point WORK at test_data)
bash $WORK/pipeline/run_all.sh
```

Verify: `verify_essentials.sh` passes on the subset; projected MalariaGEN variants land in the dhps (chr14) and dhfr-ts (chr5) CDS of each non-ref strain; Pvs230 (chr4) yields a BUSTED JSON whose p-value shifts under `--multiple-hits Double`; `hubCheck` passes on the 5-strain hub.

## Provenance / public sources

- **Assemblies (GenBank):** PvP01 `GCA_900093555.2`, PvW1 `GCA_914969965.1`, PAM `GCA_949152365.1`, PvT01 `GCA_900093545.1`, MHC087 `GCA_040114635.1` — via NCBI Datasets.
- **Annotations:** PlasmoDB-68 (PvP01, PvW1, PAM); NCBI RefSeq/GenBank (PvT01, MHC087).
- **Cohort VCF:** MalariaGEN Pv4 (1,895 samples vs PvP01) — <https://www.malariagen.net/> Pv4 data release.
- **Contig map:** derived from `work/01_chains/GCA_900093555.2.<strain>.cleaned.chain.gz` (in repo).
