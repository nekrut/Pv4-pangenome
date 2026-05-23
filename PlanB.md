# Plan B: Pangenome-based variant calling for Pv4 [local] — FINAL

Build a *P. vivax* pangenome from 9 chromosome-level NCBI assemblies, map fastp-trimmed Pv4 reads opportunistically with `vg giraffe`, call haploid variants with `vg call`, and merge into a cohort VCF.

This document is **execution-ready**: every step has its literal command. Run from `/media/anton/data/sandbox/Pv4/` with `.loom/env/` active. Steps are idempotent (skip-if-output-exists guards in scripts).

---

## Decisions (locked)

| # | Topic | Decision |
|---|-------|----------|
| 1 | Ploidy | **1** (haploid; clonal *P. vivax* asexual stages) |
| 2 | Reference set | **Tier A** — 9 chromosome-level assemblies |
| 3 | PGGB `-p` (mapping identity) | **90** (capture `vir` family divergence) |
| 4 | Organelles | **Drop** — length filter ≥100 kb pre-PGGB removes apicoplast (~30 kb) and mitochondrion (~6 kb) |
| 5 | Mapping cadence | **Opportunistic** — watcher dispatches giraffe per-sample as fastp completes |
| 6 | Variant caller | **`vg call`** only, with `-a` (genotype every snarl, for clean cohort merge) |
| 7 | Index pipeline | Explicit `vg gbwt` / `vg snarls` / `vg minimizer` (skip `vg autoindex` — known broken on PGGB GFAs, vg#4302) |
| 8 | Mapping output | giraffe → **GAM**; `vg pack` consumes GAM; `vg surject` produces inspection BAM in parallel |

## Validation notes

- Index pipeline cross-checked against vg wiki "Mapping short reads with Giraffe" and "SV Genotyping and variant calling" pages.
- `vg call -a` semantics confirmed: "Genotype every snarl, including reference calls (use to compare multiple samples)."
- `vg pack -Q 5` MQ floor is the documented default for short reads.

---

## PGGB input set (Tier A, 9 chromosome-level assemblies)

| Accession | Isolate / name | Submitter | Scaffold N50 |
|---|---|---|---|
| GCA_900093555.2 | **PvP01** (reference path) | WTSI | 1.76 Mb |
| GCA_900178095.1 | PvHMP-013 | SC | 1.62 Mb |
| GCA_900093535.1 | PvC01 | WTSI | 1.59 Mb |
| GCA_900093545.1 | PvT01 | WTSI | 1.56 Mb |
| GCA_014843675.1 | NB45 | USF | 1.53 Mb |
| GCA_014843685.1 | LZCH1476 | USF | 1.63 Mb |
| GCA_014843935.1 | LZCH1886 | USF | 1.62 Mb |
| GCA_014843945.1 | LZCH1720 | USF | 1.62 Mb |
| GCA_000002415.2 | Salvador I | TIGR | 1.68 Mb |

Stored at `.loom/pggb_input_set.txt` (one accession per line).

---

## Output layout

```
refs/                          # raw FASTAs from NCBI (one per accession)
  <accession>.fa
refs_pansn/                    # PanSN-renamed, ≥100kb contigs only
  <accession>.fa
pggb_input.fa.gz               # concatenated PanSN-renamed FASTA (9 haps)
pggb_input.fa.gz.fai
pggb_input.fa.gz.gzi
pggb_out/                      # PGGB outputs
  pggb_input.fa.gz.*.smooth.final.gfa   # the graph
  pggb_input.fa.gz.*.smooth.final.og    # odgi binary
  pggb_input.fa.gz.*.smooth.final.GCA_900093555.2.vcf.gz  # PGGB-emitted VCF vs PvP01
  *.png                                  # 1D/2D viz
vg_idx/
  pv.gbz                       # GBZ graph (giraffe input)
  pv.dist                      # distance index
  pv.min                       # minimizer index
  pv.snarls                    # snarl decomposition (vg call input)
  pv_ref_paths.txt             # PvP01 path names for surjection
gam/
  <RUN>.gam                    # giraffe output (graph alignments)
bam/
  <RUN>.bam                    # surjected to PvP01, sorted, indexed
  <RUN>.bam.bai
pack/
  <RUN>.pack                   # vg pack read-support
vcf/
  <RUN>.vcf.gz                 # vg call per-sample VCF
  <RUN>.vcf.gz.tbi
cohort/
  pv4.merged.vcf.gz            # bcftools merge across samples
  pv4.filtered.vcf.gz          # QUAL ≥ 20, DP band [5, 3×median]
.loom/
  run_giraffe.sh               # per-sample mapper+caller
  giraffe_watcher.sh           # opportunistic dispatcher
  pggb_input_set.txt
  fastp.log, download.log
```

---

## Steps

### Step 1 — Install tools {#plan-b-step-1}

Routing: local · Tool: conda

```bash
cd /media/anton/data/sandbox/Pv4
conda install -p .loom/env -c bioconda -c conda-forge -y \
  pggb vg samtools bcftools seqkit fastix tabix multiqc
.loom/env/bin/pggb --version
.loom/env/bin/vg version | head -1
```

Already installed in earlier turns: `ncbi-datasets-cli`, `matplotlib`, `seaborn`, `pandas`, `jq`, `fastp`. Idempotent.

---

### Step 2 — NCBI summary {#plan-b-step-2} ✅ DONE

Already completed. Output: `.loom/pvivax_assemblies.jsonl` (20 records).

---

### Step 3 — Survey + plot {#plan-b-step-3} ✅ DONE

Already completed. Outputs: `assembly_survey.tsv`, `assembly_survey.png`. Tier A selected.

---

### Step 4 — Download Tier A FASTAs + PanSN-rename + drop organelles {#plan-b-step-4}

Routing: local · Tool: `datasets download` + `seqkit` + `fastix`

```bash
cd /media/anton/data/sandbox/Pv4
mkdir -p refs refs_pansn

# 4a. Download zipped genomes from NCBI
.loom/env/bin/datasets download genome accession \
  --inputfile .loom/pggb_input_set.txt \
  --include genome \
  --filename .loom/tier_a.zip

# 4b. Extract; flatten to refs/<accession>.fa
unzip -o .loom/tier_a.zip -d .loom/tier_a_unzipped
for d in .loom/tier_a_unzipped/ncbi_dataset/data/GC*/; do
  acc=$(basename "$d")
  fa=$(ls "$d"/*.fna 2>/dev/null | head -1)
  [ -n "$fa" ] && cp "$fa" "refs/${acc}.fa"
done
ls -1 refs/   # expect 9 files

# 4c. Drop contigs <100 kb (organelles + small unplaced) and PanSN-rename
#     New name: <accession>#1#<original_contig_name>
for fa in refs/*.fa; do
  acc=$(basename "$fa" .fa)
  .loom/env/bin/seqkit seq -m 100000 "$fa" \
    | .loom/env/bin/fastix -p "${acc}#1#" - \
    > "refs_pansn/${acc}.fa"
done

# 4d. Concatenate, bgzip, index
cat refs_pansn/*.fa | .loom/env/bin/bgzip -@ 18 > pggb_input.fa.gz
.loom/env/bin/samtools faidx pggb_input.fa.gz

# 4e. Sanity check
echo "Haplotypes in input:"
grep -c '^>' refs_pansn/*.fa | sed 's|refs_pansn/||;s/.fa:/\t/'
echo "Total bases:"
.loom/env/bin/seqkit stats pggb_input.fa.gz
```

Expected: 9 haplotypes, ~250–270 Mb total (9 × ~28–30 Mb).

---

### Step 5 — PGGB build {#plan-b-step-5}

Routing: local · Tool: `pggb`

```bash
cd /media/anton/data/sandbox/Pv4
mkdir -p pggb_out

.loom/env/bin/pggb \
  -i pggb_input.fa.gz \
  -o pggb_out \
  -n 9 \
  -p 90 \
  -s 5000 \
  -k 23 \
  -G 700,900,1100 \
  -t 18 \
  -V GCA_900093555.2:# \
  2>&1 | tee .loom/pggb.log
```

**Run time estimate:** ~30–90 min on 18 cores for 9 × 30 Mb haplotypes.

Outputs in `pggb_out/`:
- `*.smooth.final.gfa` — the graph
- `*.smooth.final.og` — odgi binary
- `*.smooth.final.GCA_900093555.2.vcf.gz` — PGGB-emitted VCF vs PvP01 (records SVs/SNPs from the graph itself, *not* per-sample reads)
- `*.png` — 1D/2D visualizations

```bash
# Define for downstream steps
GFA=$(ls pggb_out/*.smooth.final.gfa | head -1)
echo "$GFA" > .loom/pggb_gfa_path.txt
```

---

### Step 6 — Build giraffe indices from the GFA {#plan-b-step-6}

Routing: local · Tool: `vg gbwt`, `vg snarls`, `vg minimizer`

```bash
cd /media/anton/data/sandbox/Pv4
mkdir -p vg_idx
GFA=$(cat .loom/pggb_gfa_path.txt)

# 6a. GFA → GBZ (graph + haplotype index together)
.loom/env/bin/vg gbwt --gbz-format -g vg_idx/pv.gbz -G "$GFA" --progress

# 6b. Distance index (giraffe + vg call need this)
.loom/env/bin/vg gbwt --dist-name vg_idx/pv.dist vg_idx/pv.gbz --progress

# 6c. Minimizer index (giraffe needs this)
.loom/env/bin/vg minimizer -d vg_idx/pv.dist -o vg_idx/pv.min vg_idx/pv.gbz -t 18

# 6d. Snarl decomposition (vg call needs this)
.loom/env/bin/vg snarls -T vg_idx/pv.gbz > vg_idx/pv.snarls

# 6e. List PvP01 paths for giraffe surjection target
.loom/env/bin/vg paths -L -x vg_idx/pv.gbz \
  | grep '^GCA_900093555\.2#' \
  > vg_idx/pv_ref_paths.txt
wc -l vg_idx/pv_ref_paths.txt   # expect 14 (PvP01 nuclear chromosomes)
```

---

### Step 7 — Per-sample giraffe + pack + call (helper script) {#plan-b-step-7}

Routing: local · Tool: `vg giraffe`, `vg pack`, `vg call`, `vg surject`, `samtools`, `bcftools`

This is the per-sample pipeline. The watcher (step 7b) calls this once per sample.

**`.loom/run_giraffe.sh`** (created by this plan, idempotent):

```bash
#!/bin/bash
set -euo pipefail
RUN="$1"

# Idempotent skip
if [ -f "vcf/${RUN}.vcf.gz.tbi" ]; then
  echo "[$RUN] already done"
  exit 0
fi

# Inputs (fastp output)
R1="fastp_qc/${RUN}/${RUN}_1.fastq.gz"
R2="fastp_qc/${RUN}/${RUN}_2.fastq.gz"
[ -f "$R1" ] && [ -f "$R2" ] || { echo "[$RUN] fastp output missing; skip"; exit 0; }

mkdir -p gam bam pack vcf

# 7.1 Map to graph (GAM output)
if [ ! -f "gam/${RUN}.gam" ]; then
  .loom/env/bin/vg giraffe \
    -Z vg_idx/pv.gbz \
    -d vg_idx/pv.dist \
    -m vg_idx/pv.min \
    -f "$R1" -f "$R2" \
    --sample "$RUN" \
    --read-group "ID:${RUN} SM:${RUN} LB:${RUN} PL:ILLUMINA" \
    --output-format GAM \
    --threads 4 \
    --progress \
    > "gam/${RUN}.gam"
fi

# 7.2 Surject to BAM (PvP01 backbone) for inspection
if [ ! -f "bam/${RUN}.bam.bai" ]; then
  .loom/env/bin/vg surject \
    -x vg_idx/pv.gbz \
    --threads 4 \
    --bam-output \
    --sample "$RUN" \
    --read-group "ID:${RUN} SM:${RUN} LB:${RUN} PL:ILLUMINA" \
    --prune-low-cplx \
    --interleaved \
    -F vg_idx/pv_ref_paths.txt \
    "gam/${RUN}.gam" \
    | .loom/env/bin/samtools sort -@ 4 -o "bam/${RUN}.bam" -
  .loom/env/bin/samtools index -@ 4 "bam/${RUN}.bam"
fi

# 7.3 Read-support pack (for vg call)
if [ ! -f "pack/${RUN}.pack" ]; then
  .loom/env/bin/vg pack \
    -x vg_idx/pv.gbz \
    -g "gam/${RUN}.gam" \
    -o "pack/${RUN}.pack" \
    -Q 5 \
    -t 4
fi

# 7.4 Variant call (haploid, all snarls)
if [ ! -f "vcf/${RUN}.vcf.gz" ]; then
  .loom/env/bin/vg call \
    vg_idx/pv.gbz \
    -k "pack/${RUN}.pack" \
    -r vg_idx/pv.snarls \
    -s "$RUN" \
    -a \
    --ploidy 1 \
    -t 4 \
    | .loom/env/bin/bgzip -@ 4 \
    > "vcf/${RUN}.vcf.gz"
  .loom/env/bin/tabix -p vcf "vcf/${RUN}.vcf.gz"
fi

echo "[$RUN] done"
```

Make executable:

```bash
chmod +x .loom/run_giraffe.sh
```

**Test on one sample first (before launching watcher):**

```bash
ls fastp_qc/ | head -1   # pick a completed run
.loom/run_giraffe.sh ERR021983   # or whichever
```

---

### Step 7b — Watcher loop (opportunistic dispatch) {#plan-b-step-7b}

**`.loom/giraffe_watcher.sh`**:

```bash
#!/bin/bash
# Opportunistic dispatcher: poll fastp_qc/, run pipeline on any new completed sample.
# Idempotent: run_giraffe.sh skips already-finished runs.

cd /media/anton/data/sandbox/Pv4
mkdir -p .loom/watcher_log

while true; do
  ls fastp_qc/ 2>/dev/null | while read RUN; do
    if [ -f "fastp_qc/${RUN}/${RUN}.fastp.json" ] \
    && [ ! -f "vcf/${RUN}.vcf.gz.tbi" ]; then
      echo "$RUN"
    fi
  done | xargs -P4 -I{} bash -c '
    .loom/run_giraffe.sh "$1" >> .loom/watcher_log/$(date +%Y%m%d).log 2>&1 \
      || echo "FAIL $1 $(date)" >> .loom/watcher_log/failures.log
  ' _ {}
  sleep 60
done
```

Make executable and launch in background:

```bash
chmod +x .loom/giraffe_watcher.sh
nohup .loom/giraffe_watcher.sh > .loom/watcher_log/main.log 2>&1 &
echo $! > .loom/watcher.pid
```

Stop with: `kill $(cat .loom/watcher.pid)`.

**Concurrency math:** xargs -P4 × giraffe -t 4 = 16 cores reserved for mapping. Leaves 4 cores free for fastp watcher and other work.

---

### Step 8 — Surject already in step 7.2 {#plan-b-step-8} ✅ folded into Step 7

No separate action.

---

### Step 9 — Variant call already in step 7.4 {#plan-b-step-9} ✅ folded into Step 7

No separate action.

---

### Step 10 — Cohort merge + filter {#plan-b-step-10}

Run after the watcher finishes (or periodically — `bcftools merge` is idempotent for the same input set).

```bash
cd /media/anton/data/sandbox/Pv4
mkdir -p cohort

# 10a. Merge per-sample VCFs
ls vcf/*.vcf.gz > .loom/vcf_list.txt
.loom/env/bin/bcftools merge \
  --file-list .loom/vcf_list.txt \
  --threads 18 \
  -Oz -o cohort/pv4.merged.vcf.gz
.loom/env/bin/tabix -p vcf cohort/pv4.merged.vcf.gz

# 10b. Compute median DP for the coverage-band filter
MEDIAN_DP=$(.loom/env/bin/bcftools query -f '%INFO/DP\n' cohort/pv4.merged.vcf.gz \
  | sort -n | awk '{a[NR]=$1} END {print a[int(NR/2)]}')
MAX_DP=$((MEDIAN_DP * 3))
echo "Median DP=$MEDIAN_DP, max DP=$MAX_DP"

# 10c. Filter: QUAL ≥ 20, 5 ≤ DP ≤ 3×median
.loom/env/bin/bcftools filter \
  -e "QUAL<20 || INFO/DP<5 || INFO/DP>${MAX_DP}" \
  -s LowQualOrDP \
  --threads 18 \
  -Oz -o cohort/pv4.filtered.vcf.gz \
  cohort/pv4.merged.vcf.gz
.loom/env/bin/tabix -p vcf cohort/pv4.filtered.vcf.gz

# 10d. Cohort stats
.loom/env/bin/bcftools stats cohort/pv4.filtered.vcf.gz > cohort/pv4.stats.txt
echo "PASS variants: $(.loom/env/bin/bcftools view -f PASS cohort/pv4.filtered.vcf.gz | grep -vc '^#')"
```

---

## Execution order summary

For a Haiku-class executor, run in this exact order, checking exit code after each:

```
1. step-1   conda install
2. step-4   download + PanSN + bgzip
3. step-5   pggb build              [longest single step, 30–90 min]
4. step-6   vg index                [~10 min]
5. step-7   chmod + test on 1 sample
6. step-7b  launch watcher in background
   (let it run alongside fastp + downloads; samples accumulate in vcf/)
7. step-10  cohort merge + filter   [run when no new samples expected]
```

Steps 2 and 3 of Plan B are already complete (NCBI summary + survey plot).

## Failure modes & recovery

| Symptom | Cause | Fix |
|---|---|---|
| `pggb` OOMs | 9 × 30 Mb is small; should not happen on this hardware. If it does → reduce `-G` to `700` only. | Edit step 5 `-G 700` |
| `vg gbwt -G` errors on GFA | PGGB GFA may have non-canonical lines. | `vg convert -g <gfa> -t 18 -p > pv.vg` then `vg gbwt --gbz-format -g pv.gbz pv.vg` |
| `vg surject` exits with code 1 on some BAMs | known issue on long indels in graph | Add `--prune-low-cplx` (already in script) |
| Empty VCF for a sample | low coverage or unmapped library | Inspect `bam/<RUN>.bam` → samtools flagstat |
| Watcher running hot but no progress | giraffe stuck on a malformed FASTQ | check `.loom/watcher_log/failures.log`, manually re-run |
| `fastix` not found | conda install failed silently | `conda install -p .loom/env -c bioconda -y fastix` |

## Things this plan does **not** do (out of scope)

- Annotation (snpEff / VEP)
- Population genetics analyses (PCA, F<sub>ST</sub>, admixture)
- Per-region (e.g. drug resistance loci) targeted summaries
- Copy-number variation calling (graph captures presence/absence, not multi-copy)
- Mitochondrial / apicoplast variant calling (organelles dropped at step 4c by design)

These are natural follow-ups once `cohort/pv4.filtered.vcf.gz` exists.
