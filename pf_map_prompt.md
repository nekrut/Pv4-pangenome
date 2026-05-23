# GPU read mapping against a PGGB pangenome with Parabricks Giraffe

Goal: take a PGGB pangenome graph and map a set of paired-end short reads
against it on a single GPU using NVIDIA Parabricks Giraffe. One BAM per
sample. Nothing else — no CPU comparison, no variant calling.

## Inputs I'll provide

- PGGB graph (GFA, may be .gfa or .gfa.gz): <GFA_PATH>
- Reads (paired-end FASTQ.gz). Either:
    (a) a single sample: <R1.fq.gz> <R2.fq.gz>, or
    (b) a directory: <READS_DIR> with files named *_R1.fq.gz / *_R2.fq.gz
        (one pair per sample, same prefix)
- Reference sample name as it appears in the PanSN path names of the
  graph (the part before the first `#`). Example: "3D7", "Pf3D7",
  "PfDD2_v3". This must match exactly.
- Output directory: <OUT_DIR>

If any of these is missing or ambiguous, stop and ask. Don't substitute,
don't pull from public sources without confirmation.

## Required toolchain

You'll need:

- **vg version 1.70.0 exactly.** Older vg emits minimizer v10, newer vg
  emits GBZ v2; Parabricks 4.7.0-1 rejects both. Install with:

      conda create -p ./env_vg170 -c bioconda -c conda-forge -y vg=1.70.0
      VG=./env_vg170/bin/vg
      $VG version  # must say v1.70.0

- **Parabricks 4.7.0 or newer.** Pull:

      docker pull nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1

- **Docker GPU passthrough.** Verify:

      docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi -L

  Should print one GPU line. If it fails, the host needs
  `nvidia-container-toolkit` installed + `nvidia-ctk runtime configure
  --runtime=docker` + `systemctl restart docker`. Stop and report if
  it doesn't work; do not try workarounds.

## Working directory

Create <OUT_DIR>. Copy (not symlink) the GFA into it — Parabricks runs
in a docker container with `-v $(pwd):/work`, so it can only see files
under that directory. Reads can stay where they are; you'll bind-mount
their directory too.

## Pipeline

### Step 1 — Verify inputs

Confirm the GFA exists and the ref sample appears in path names:

    $VG paths --list <GFA_PATH> | head -20

If <NAME> isn't a sample, stop and list the available samples. Don't
guess.

If the GFA is gzipped, gunzip it next to the original (Parabricks and
vg autoindex both want plain .gfa).

### Step 2 — Build vg 1.70 indexes (once per graph)

    cd <OUT_DIR>
    $VG autoindex --workflow giraffe --gfa graph.gfa --prefix idx -t $(nproc)

That produces idx.giraffe.gbz, idx.shortread.withzip.min, idx.dist,
idx.shortread.zipcodes.

Then mark the reference sample (REQUIRED — Parabricks rejects graphs
with no reference-sense paths and you'll get the cryptic error
"No reference or non-alt-allele generic paths available in the graph"):

    $VG gbwt -Z idx.giraffe.gbz --set-reference <NAME> \
             --gbz-format -g idx.r.gbz
    $VG minimizer -d idx.dist -z idx.r.zipcodes \
             -o idx.r.min idx.r.gbz

The final index set you'll feed Parabricks: idx.r.gbz, idx.r.min,
idx.dist, idx.r.zipcodes.

### Step 3 — Map each sample

For each (R1, R2) pair, with output base name <SAMPLE>:

    nvidia-smi dmon -s um -d 1 -c 9999 > <SAMPLE>.dmon.log &
    DMON=$!

    /usr/bin/time -v docker run --rm --gpus '"device=0"' \
        -e CUDA_VISIBLE_DEVICES=0 \
        -v <OUT_DIR>:/work \
        -v <reads parent dir>:/reads:ro \
        -w /work \
        nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1 \
        pbrun giraffe \
            --gbz-name idx.r.gbz \
            --minimizer-name idx.r.min \
            --dist-name idx.dist \
            --zipcodes-name idx.r.zipcodes \
            --in-fq /reads/<R1> /reads/<R2> \
            --out-bam <SAMPLE>.bam \
            --num-gpus 1 --nstreams 2 --num-cpu-threads-per-gpu 18 \
            2> <SAMPLE>.time.log

    kill $DMON 2>/dev/null || true

`--nstreams 2 --num-cpu-threads-per-gpu 18` is the conservative,
known-good config on a 24 GB GPU. Going to `--nstreams 3` uses ~22 GB
VRAM and may OOM on graphs larger than a few hundred Mb of input
sequence; only try it once and only if you want speed and have headroom.

### Step 4 — Per-sample summary

For each output BAM, record from the time log:

- elapsed wall clock
- peak host RSS
- exit status

And from the dmon log, the max VRAM (FB) and max SM% — quick sanity
check that the GPU was actually used.

### Step 5 — Write a single results file

`<OUT_DIR>/MAPPING_SUMMARY.md` with:

- versions: vg, Parabricks (`docker images`), nvidia driver
- graph: file path, size, number of nodes/edges/paths
- one row per sample: wall, peak RAM, max VRAM, max SM%,
  total mapped, MAPQ≥30, output BAM path
- any errors verbatim

## Stop conditions

- Docker can't see GPU → stop, report.
- Parabricks rejects an index file (GBZ version, minimizer version,
  missing reference paths) → stop, report. Do NOT silently rebuild
  with a different vg version. The vg 1.70 pin is intentional.
- VRAM OOM during a run → stop, report. Switch to nstreams=2 if you
  were on 3, but only retry once. If 2 OOMs too, it's a graph-size
  finding, not something to chase.
- A sample errors out → log it, skip to the next sample. Don't abort
  the whole batch on one bad input.

## What to leave out

- No "let me also try the CPU version" comparison. GPU only.
- No editing of the GFA or reads beyond gunzip.
- No assumptions about the reference sample name. Ask if unsure.
- No retries with different parameters after the first run of a sample
  succeeds. Move on.
