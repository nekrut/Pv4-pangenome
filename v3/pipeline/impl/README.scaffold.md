# Pk v1 — P. knowlesi pangenome + selection-scan pipeline

Seven-strain *Plasmodium knowlesi* panel (taxid 5850). Implements all 11 phases
from `../Pv4/v3/pipeline/LOCAL.md`, parameterized for Pk.

## Quick start

```bash
# 1. Edit species.conf — set WORK, SMOKE_CHROM after downloading assemblies
#    COHORT_VCF_DIR stays TODO until a Pk population VCF exists.

# 2. Download assemblies + annotations
bash pipeline/setup/fetch_assemblies.sh

# 3. Update SMOKE_CHROM in species.conf (chromosome 1 accession of Pk_H)

# 4. Smoke test (~45 min on 32 cores, 1 GPU)
bash pipeline/smoke_test.sh

# 5. Full run (~24 hrs wall)
bash pipeline/run_all.sh
```

## Files

| File | Purpose |
|---|---|
| `pipeline/species.conf` | All tunable parameters: strains, accessions, thresholds |
| `pipeline/setup/fetch_assemblies.sh` | NCBI Datasets download + gffread proteome extraction |
| `pipeline/00_validate_conf.sh` | Pre-flight config sanity checks |
| `pipeline/01_inventory.sh` … `11_project_vcf.sh` | Phase scripts (idempotent) |
| `pipeline/run_all.sh` | Orchestrator (runs 00 → 11 in dependency order) |
| `pipeline/smoke_test.sh` | 3-strain × 300 kb smoke test including TOGA2 |
| `pipeline/lib/run_in_container.sh` | Docker wrapper (pinned images) |
| `pipeline/lib/run_in_apptainer.sh` | Singularity/HPC alternative |
| `pipeline/lib/verify_essentials.sh` | Validate all 27 *-essential outputs |
| `pipeline/lib/fix_gff_chroms.sh` | Rewrite GFF chr names to match FASTA |
| `pipeline/scripts/` | Python helpers (triage, merge, consensus, MSA, PanSN) |

## Strains

| Short name | GenBank accession | Role |
|---|---|---|
| Pk_H | GCA_000006355.3 | **Reference** (H strain, most complete) |
| Pk_ANKA | GCA_009792815.1 | Anchor |
| Pk_YH1 | GCA_014858985.1 | Anchor |
| Pk_TAM | GCA_021201615.1 | Anchor |
| Pk_A1H1 | GCA_014858965.1 | Non-anchor |
| Pk_COF | GCA_019833925.1 | Non-anchor |
| Pk_SIM | GCA_963506685.1 | Non-anchor |

## Known TODOs in species.conf

- `COHORT_VCF_DIR` — no MalariaGEN-style population VCF exists for *P. knowlesi* yet.
  Phase J (CrossMap VCF projection) will be skipped automatically until this is set.
- `SMOKE_CHROM` — update to the GenBank accession of chromosome 1 in Pk_H after download.
- `CHROM_RENAME` — built automatically by `fix_gff_chroms.sh` after assemblies are downloaded.

## Container runtime

All tools run inside pinned Docker images (see `run_in_container.sh`).
On HPC clusters without Docker, set `RUN_IN_CONTAINER=apptainer` to use
Singularity/Apptainer instead.

## Reference

Source recipe: `/media/anton/data/sandbox/Pv4/v3/pipeline/LOCAL.md`
