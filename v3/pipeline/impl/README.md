# Clean LOCAL pipeline — implementation (work in progress)

The runnable phase scripts for the LOCAL (container) pipeline described in [`../LOCAL.md`](../LOCAL.md) and [`../PIPELINE_EXPLANATION.md`](../PIPELINE_EXPLANATION.md). Derived from the P. knowlesi scaffold (`/media/anton/data/sandbox/Pk/v1/pipeline`, generalized + species-agnostic via `species.conf`) and carrying the **22 fixes** from the first end-to-end run against the Pv4 test panel — see [`../SCAFFOLD_FIXES.md`](../SCAFFOLD_FIXES.md) for the full list and rationale.

## Status (first end-to-end run, Pv4 test panel)

| Phase | State |
|---|---|
| A inventory (mash + BUSCO) | ✅ green |
| B soft-mask | ✅ green |
| C.1–3 align + chains | ✅ green |
| C.4 annotation projection | ✅ green (Liftoff; TOGA skipped when rescue queue empty) |
| D PGGB graph | ✅ green |
| E consensus orthology | ✅ green (1,992 orthogroups) — but rbest/graph edge sources contributed 0; see SCAFFOLD_FIXES |
| F codon/protein MSAs | ⚠️ runs (build_msa.py on host, MSA tool tags fixed) but builds 0 MSAs — unresolved, precisely narrowed in SCAFFOLD_FIXES |
| G IQ-TREE / H BUSTED / I multiz / J VCF projection | not yet exercised |

## Layout

- `00_validate_conf.sh`, `01`–`11` phase scripts, `run_all.sh`, `smoke_test.sh`
- `lib/run_in_container.sh` — Docker wrapper (**fixed**: `-i`, valid/individual image tags, python3→pyfaidx, bgzip→samtools, …)
- `lib/run_in_apptainer.sh` — Singularity variant, **NOT yet fixed** (same stale tags as the original; apply the SCAFFOLD_FIXES tag changes before HPC use)
- `lib/verify_essentials.sh`, `lib/fix_gff_chroms.sh`
- `scripts/*.py` — phase helpers (triage, merge, consensus, MSA)
- `setup/{build_anchor_inputs,fetch_assemblies}.sh`
- `species.conf.pv4test.example` — the config used for the test run (also in `../test_data/species.conf`)

## How to run

```bash
WORK=/path/to/workspace          # holds pipeline/ + inputs/
cp species.conf.pv4test.example $WORK/pipeline/species.conf   # then edit for your panel
# stage inputs/{assemblies,annotations,proteomes,cohort_vcf}/ (see ../test_data)
export WORK
bash $WORK/pipeline/00_validate_conf.sh
bash $WORK/pipeline/run_all.sh
```

## Known unresolved (need decisions, not just fixes)

- **F — `build_msa.py` container orchestration.** It's written to run on the *host* and call `run_in_container.sh` for mafft/pal2nal/gffread, but `07_msa.sh` runs it *inside* the pyfaidx container (`cmd python3`) → docker-in-docker, wrapper path mis-resolves. Decide: run it on the host, or build one image with all MSA tools. Until then F builds 0 MSAs.
- **E — degraded consensus.** Orthogroups came from projection edges only; rbest-chain and graph-co-membership sources returned 0 edges (likely missing `.bed` inputs). Phase completes but the 3-way consensus is effectively 1-way.
- **TOGA rescue path untested** — `ghcr.io/hillerlab/toga:latest` is not publicly pullable; the CESAR rescue tail can't run for divergent panels until a working TOGA image exists. (Harmless for close panels: the rescue queue is empty and TOGA is skipped.)
- **Pk-isms remain** — log prefixes, `/tmp/Pk_*` temp names, `pk.{og,gfa}` symlink names, the apptainer wrapper. Cosmetic; doesn't affect function.
