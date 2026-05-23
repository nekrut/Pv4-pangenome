# Pv4 Dataset: Sample downloads and QC

## Context

MalariaGEN *Plasmodium vivax* (Pv4) public dataset from https://www.malariagen.net/

Source metadata: https://www.malariagen.net/wp-content/uploads/2023/11/Pv4_samples.txt

## Summary

- **Total samples in manifest:** 1,895
- **QC pass:** 1,072 (56.5%)
- **QC fail:** 823
- **With ENA accessions:** 1,642 (86.6%)
- **QC pass + ENA accession:** 897

## Metadata quality

All 15 columns fully populated except **ENA accession column (13% missing).** Geographic (lat/long), temporal (year), and provenance (study, site) fields are complete.

| Field | Status |
|-------|--------|
| Sample, Study, Site, Admin division, Country, Lat, Long, Year | Complete |
| Individual grouping, Population, % callable, QC pass, Exclusion reason, Returning traveller | Complete |
| ENA accession | 253/1895 missing (~13%) |

## Download plan: QC-pass samples with ENA accessions

**Scope:** 1,172 runs (2,344 paired-end FASTQ files) from 897 samples

**Data size:** 1.82 TB

**Status:** In progress (12 parallel curl streams)

- **Log:** `.loom/download.log`
- **Output layout:** `fastq/<RUN>/<RUN>_{1,2}.fastq.gz`
- **Started:** 2026-05-01 20:25:40Z
- **ETA:** ~6–12 hours (depends on FTP throughput)

**Known issue:** 1 run has no FASTQ on ENA: `ERR2299660` (will be skipped)

### Monitor progress

```bash
# Current disk usage
df -h /media/anton/data/sandbox/Pv4

# Downloaded file count
find fastq -name '*.fastq.gz' | wc -l

# Total size downloaded
du -sh fastq/

# Last 20 log entries (live)
tail -f .loom/download.log
```

## Plan A: fastp QC on Pv4 paired-end reads [local]

See [`PlanA.md`](PlanA.md). Status: steps 1–3 in progress (running), step 4 (multiqc aggregation) pending.

## NCBI assembly survey (Plan B step 3)

Surveyed 20 *P. vivax* assemblies on NCBI (19 unique after deduping GCA/GCF Sal-I): 9 chromosome-level, 8 scaffold, 2 contig. Plot saved to `assembly_survey.png`, table to `assembly_survey.tsv`.

**Selected Tier A (9 chromosome-level assemblies)** for PGGB input:

GCA_900093555.2 (PvP01, reference), GCA_900178095.1 (PvHMP-013), GCA_900093535.1 (PvC01), GCA_900093545.1 (PvT01), GCA_014843675.1 (NB45), GCA_014843685.1 (LZCH1476), GCA_014843935.1 (LZCH1886), GCA_014843945.1 (LZCH1720), GCA_000002415.2 (Salvador I).

All 9 have N50 in 1.5–1.8 Mb band; clean candidates for graph construction.

## Plan B: Pangenome-based variant calling [local]

See [`PlanB.md`](PlanB.md) for the **final, execution-ready plan**.

Decisions locked: ploidy=1, Tier A (9 chr-level), PGGB `-p 90`, organelles dropped via length filter ≥100 kb, opportunistic watcher, `vg call -a` only, explicit `vg gbwt` index pipeline (avoids known `vg autoindex` GFA bug vg#4302).

Validation: commands cross-checked against vg wiki + GitHub issues via web_search (no inter-agent tool available; web search used as the second-opinion review per request).

Status: steps 1–5 in progress.

- [x] Step 1: conda env `.loom/bio_env/` created (pggb, vg 1.73.0, samtools, bcftools, seqkit, tabix)
- [x] Step 2: NCBI summary done (`.loom/pvivax_assemblies.jsonl`)
- [x] Step 3: Survey + plot done (`assembly_survey.tsv`, `assembly_survey.png`)
- [x] Step 4: Tier A downloaded, organelles dropped (≥100 kb filter), PanSN-renamed → `pggb_input.fa.gz` (191 seqs, 224 Mb)
- [x] Step 5: PGGB complete → `pggb_out/*.smooth.fix.gfa` (125M)
- [x] Step 6: vg index complete → `vg_idx/pv.gbz` (32M, reference_samples=GCA_900093555.2), `pv.dist` (78M), `pv.min` (285M), `pv.snarls` (2.1M)
- [x] Step 7 test: ERR021983 pipeline validated end-to-end (giraffe→surject→pack→call); 224K variant sites, 171K SNPs
- [x] Step 7b: Watcher running (PID in `.loom/watcher.pid`), dispatching samples opportunistically

Helper scripts: `.loom/run_giraffe.sh`, `.loom/giraffe_watcher.sh`

```loom-session
id: 019de532-5eb6-72d2-88e8-20015590bd15
started_at: 2026-05-12T19:15:39.940Z
ended_at: 2026-05-12T19:16:58.324Z
notebook: notebook.md
orphaned_active_steps: 0
```
