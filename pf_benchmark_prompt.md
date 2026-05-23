You're picking up a Parabricks GPU Giraffe vs CPU vg giraffe benchmark on a
P. falciparum pangenome. The workstation has the toolchain prepared from a
prior run on simulated genomes. Hard requirements and gotchas are listed
below — read them all before starting.

## Inputs I'll give you (ask if not provided)

- PGGB graph (GFA, possibly gzipped): <PATH>
- Real Illumina paired-end reads: <R1.fq.gz> <R2.fq.gz>
- Reference sample name as it appears in the PanSN path names of the
  GFA (e.g. "3D7", "Pf3D7", "PfDD2"): <NAME>
- Sample ID label for the BAM read group: <SAMPLE>

If any of these is missing, stop and ask. Don't substitute, don't pull
from public sources without confirmation.

## Hardware

20-core / 40-thread Xeon W-2255 (nproc reports 20), 252 GB RAM,
1× NVIDIA RTX A5000 24 GB, NVMe. Driver 595.x, CUDA 13.2 host.
Only one GPU; force CUDA_VISIBLE_DEVICES=0.

## Toolchain already installed

- vg v1.70.0 at /media/anton/data/sandbox/udt1/.loom/vg_v170/bin/vg
  ** USE THIS EXACT VERSION ** — older vg emits minimizer v10,
  newer vg emits GBZ v2; both are rejected by Parabricks 4.7.0-1.
  vg 1.70 is the only version whose indexes Parabricks 4.7.0 will load.
- samtools at /media/anton/data/sandbox/udt1/.loom/env/bin/samtools
- wgsim available (not needed unless you're simulating)
- Parabricks image: nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1 (pulled,
  ~9 GB). Verify with `docker images | grep parabricks` before starting.
- Docker GPU passthrough: confirm with
    docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi -L
  If that fails, stop — toolkit isn't configured.

## Working directory

Create /media/anton/data/sandbox/udt1/bench_pf/ and do everything there.
Copy (do not symlink) input files into it. Parabricks runs in Docker
with -v $(pwd):/work, so it can only see files in that directory.

## Pipeline

### Step 1 — Stage inputs

Confirm GFA and reads exist and are readable. If GFA is gzipped, gunzip
it inside bench_pf/. Note its size.

### Step 2 — Build vg 1.70 indexes from the GFA

    VG=/media/anton/data/sandbox/udt1/.loom/vg_v170/bin/vg

    $VG autoindex --workflow giraffe --gfa graph.gfa --prefix idx -t 18

This produces idx.giraffe.gbz, idx.shortread.withzip.min, idx.dist,
idx.shortread.zipcodes.

Then mark the reference sample (REQUIRED — Parabricks rejects graphs
with no reference-sense paths):

    $VG gbwt -Z idx.giraffe.gbz --set-reference <NAME> \
             --gbz-format -g idx.r.gbz
    $VG minimizer -d idx.dist -z idx.r.zipcodes \
             -o idx.r.min idx.r.gbz

If --set-reference fails because <NAME> isn't a sample in the graph,
list paths with `$VG paths --list graph.gfa | head` and ask me.

### Step 3 — CPU baseline (single trial, /usr/bin/time -v)

    /usr/bin/time -v $VG giraffe \
        -Z idx.r.gbz -m idx.r.min -d idx.dist \
        --zipcode-name idx.r.zipcodes \
        -f R1.fq.gz -f R2.fq.gz -p -t 20 \
        --output-format BAM > cpu.bam

Capture Elapsed, Maximum resident set size, Percent of CPU from stderr.

### Step 4 — GPU run B (conservative)

In parallel, log nvidia-smi to a file:

    nvidia-smi dmon -s um -d 1 -c 9999 > gpu_B.dmon.log &

Then run:

    /usr/bin/time -v docker run --rm --gpus '"device=0"' \
        -e CUDA_VISIBLE_DEVICES=0 \
        -v $(pwd):/work -w /work \
        nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1 \
        pbrun giraffe \
            --gbz-name idx.r.gbz \
            --minimizer-name idx.r.min \
            --dist-name idx.dist \
            --zipcodes-name idx.r.zipcodes \
            --in-fq R1.fq.gz R2.fq.gz \
            --out-bam gpu_B.bam \
            --num-gpus 1 --nstreams 2 --num-cpu-threads-per-gpu 18

Kill the dmon process after the run finishes.

### Step 5 — GPU run C (aggressive)

Same as step 4 but `--nstreams 3` and log to gpu_C.dmon.log /
gpu_C.bam. Watch for VRAM crossing ~22 GB — A5000 is 24 GB.
If pbrun errors with OOM or driver issues, stop and report; do
not retry with smaller streams.

### Step 6 — Correctness comparison

For each of cpu.bam, gpu_B.bam, gpu_C.bam compute:

  total reads, mapped (-F 4), MAPQ≥30 (-F 4 -q 30),
  median TLEN of properly paired, per-chromosome dispatch (top 14)

Use samtools at /media/anton/data/sandbox/udt1/.loom/env/bin/samtools.
Numbers should match CPU within 1%. Larger drift is a real finding,
flag it; do not assume it's noise.

### Step 7 — GPU utilization summary

For each gpu_*.dmon.log, parse the SM% column (col 2) and FB memory
column (col 8). Report N samples, median SM%, p90 SM%, max SM%,
max VRAM.

### Step 8 — Write bench_pf/REPORT.md

One markdown file with:
- versions (vg, Parabricks, driver), reads source + coverage,
  graph size, index sizes
- exact commands run
- table: wall time, host RAM peak (from /usr/bin/time -v),
  GPU SM% median/max, VRAM peak — for each of the three runs
- table: correctness deltas vs CPU baseline
- any errors or unexpected behavior, with full log excerpts

## Stop conditions (per the original spec)

- Docker can't see GPU → stop, report.
- Parabricks rejects an index (GBZ version, minimizer version,
  missing reference paths) → stop, report. Do NOT silently rebuild
  with different vg versions. The toolchain pin (vg 1.70 + the
  --set-reference rebuild) is what worked once; if it stops working
  on this graph, that's the finding.
- VRAM OOMs during run C → stop after run C, report. Don't drop to
  nstreams=2 and retry; B already covered that config.
- Any BAM correctness metric differs by >5% between CPU and GPU →
  stop after collecting BAMs, report. Don't editorialize, just numbers.

## What to leave out

- No "let me also try X" experiments. Six runs total: CPU once,
  GPU B once, GPU C once. If a run dies, one retry only.
- No silent retries after the first failed run of a config.
- No editing the GFA or reads beyond gunzip.
- No assumptions about which sample is the reference — ask if unsure.
