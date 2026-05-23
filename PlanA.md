# Plan A: fastp QC on Pv4 paired-end reads [local]

Quality-filter and trim adapters from downloaded *P. vivax* paired-end FASTQ files. Process complete pairs in parallel, preserving originals, with per-sample JSON + HTML reports for QC summaries.

## Steps

- [x] 1. **Install fastp** {#plan-a-step-1} — Set up fastp in conda environment
  - Routing: local
  - Tool: conda (`.loom/env/`, channels: bioconda, conda-forge)
  - Result: fastp 0.22.0 installed
- [x] 2. **Identify complete pairs** {#plan-a-step-2} — Find runs with both R1 and R2 files
  - Routing: local
  - Tool: bash (`find fastq -type d -mindepth 1 -maxdepth 1`)
- [x] 3. **Run fastp in parallel** {#plan-a-step-3} — Process paired files, output to `fastp_qc/`
  - Routing: local
  - Tool: fastp via `.loom/run_fastp.sh`, parallelized with `xargs -P18`
- [ ] 4. **Aggregate QC reports** {#plan-a-step-4} — Summarize fastp JSON outputs across samples
  - Routing: local
  - Tool: multiqc or custom python/jq

## Parameters

| Step | Tool | Parameter | Default | Value | Description |
| --- | --- | --- | --- | --- | --- |
| 1 | conda | channels | `-c bioconda -c conda-forge` | `-c bioconda -c conda-forge` | Channel priority for fastp |
| 1 | conda | env path | `-p .loom/env` | `-p .loom/env` | Per-analysis isolated environment |
| 3 | fastp | --qualified_quality_phred | 15 | 20 | Minimum Phred quality for "qualified" base |
| 3 | fastp | --cut_front | off | on | Sliding-window quality trimming from 5' end |
| 3 | fastp | --cut_tail | off | on | Sliding-window quality trimming from 3' end |
| 3 | fastp | --cut_window_size | 4 | 4 | Sliding window size for quality trimming |
| 3 | fastp | --cut_mean_quality | 20 | 20 | Mean Q threshold within sliding window |
| 3 | fastp | --detect_adapter_for_pe | off | on | Auto-detect Illumina/BGI/MGI PE adapters |
| 3 | fastp | --length_required | 15 | 30 | Discard reads shorter than N after trimming |
| 3 | fastp | --compression | 4 | 6 | gzip compression level (1=fast, 9=small) |
| 3 | fastp | --thread | 2 | 1 | Per-job threads (parallelism handled by xargs) |
| 3 | xargs | -P | - | 18 | Concurrent fastp jobs (system has 20 cores) |

## Command (per-sample, run via `.loom/run_fastp.sh`)

```bash
conda run -p .loom/env fastp \
  -i fastq/$RUN/${RUN}_1.fastq.gz \
  -I fastq/$RUN/${RUN}_2.fastq.gz \
  -o fastp_qc/$RUN/${RUN}_1.fastq.gz \
  -O fastp_qc/$RUN/${RUN}_2.fastq.gz \
  --json fastp_qc/$RUN/${RUN}.fastp.json \
  --html fastp_qc/$RUN/${RUN}.fastp.html \
  --qualified_quality_phred 20 \
  --cut_front --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --detect_adapter_for_pe \
  --length_required 30 \
  --compression 6 \
  --thread 1
```

## Parallel launcher

```bash
find fastq -type d -mindepth 1 -maxdepth 1 \
  | xargs -n1 basename \
  | xargs -P18 -I{} .loom/run_fastp.sh {} \
  > .loom/fastp.log 2>&1 &
```

## Output layout

```
fastp_qc/
├── <RUN>/
│   ├── <RUN>_1.fastq.gz       # trimmed R1
│   ├── <RUN>_2.fastq.gz       # trimmed R2
│   ├── <RUN>.fastp.json       # machine-readable QC report
│   └── <RUN>.fastp.html       # human-readable QC report
```

## Notes

- Source FASTQ in `fastq/<RUN>/` is preserved (not modified).
- Some runs in this dataset contain malformed records (sequence/quality length mismatch); fastp logs the offending reads but continues.
- Plan is incremental — as more downloads complete (`fastq/<RUN>/` populated), re-running the launcher picks up new runs and skips already-processed ones if output exists (add `[ -f fastp_qc/$RUN/$RUN.fastp.json ] && exit 0` to `.loom/run_fastp.sh` for true idempotence).
- First sample (ERR021983) completed in 17 s on this hardware: 191K → 69.6K reads pass filter, 1,378 reads adapter-trimmed.
