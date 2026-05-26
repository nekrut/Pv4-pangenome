# LOCAL.md — Pangenome + Selection-Scan Pipeline

Zero-ambiguity implementation recipe. Authored in the style of `LLM-eval-paper/plan/PLAN.md` for use as an LLM-implementation benchmark and as a complete reference for human operators.

## Goal

Reproduce the 27 `*`-marked essential outputs from `writeup/OUTLINE.md` for any species panel of N (5–15) haploid assemblies, given:

1. N assembly FASTAs (`inputs/assemblies/{S}.fa`)
2. N annotation GFF3s (`inputs/annotations/{S}.gff3`)
3. Per-chromosome cohort VCFs (`inputs/cohort_vcf/{species}_{chrom}.vcf.gz`)
4. A chromosome-rename map (`inputs/annotations/{ref}_to_genbank.tsv`)

Implementer must produce 11 phase scripts plus a top-level orchestrator. Each phase script is self-contained, idempotent, and exits 0 only if all `*`-marked outputs for that phase exist and pass structural validation.

## Boilerplate (top of every phase script)

```bash
#!/usr/bin/env bash
set -euo pipefail

# Source the per-species config. species.conf MUST exist and define:
#   SPECIES, WORK, STRAINS (bash array), REF_STRAIN, ANCHOR_STRAINS,
#   COHORT_VCF_DIR, CHROM_RENAME, N_CORES, GPU, VAR_ANTIGEN_RE
source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

# All tool invocations go through the container wrapper.
# Default to docker; override RUN_IN_CONTAINER=apptainer for HPC clusters.
cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }

log() { echo "$(date +%H:%M:%S) [${BASH_SOURCE[0]##*/}] $1" >&2; }

retry_once() {
  if ! eval "$*"; then log "retry: $*"; eval "$*"; fi
}
```

**Convention used in this recipe**: code blocks show bare tool names (`mash sketch ...`, `bcftools view ...`). The implementer MUST wrap each invocation via the `cmd` alias defined above: `cmd mash sketch ...`, `cmd bcftools view ...`. This keeps the recipe readable while enforcing that every tool runs in its pinned container. The Python helper scripts in `pipeline/scripts/` (e.g. `phase_c2_triage.py`, `build_msa.py`) run via `cmd python3 pipeline/scripts/<script>.py ...` and inherit pysam/pandas/numpy from the bcftools image.

## TOOL_INVENTORY (container-first)

**Every tool runs in a container** — Docker daemon required (≥24), no conda envs, no system installs beyond `bash`, `python3`, `docker`, `samtools` (host-side for index lookups), `gh` (for issue filing only). This eliminates env-drift between machines.

Two container types used:

1. **Pinned bioconda biocontainers** (most tools): `quay.io/biocontainers/{tool}:{version}--{build}`. Pinning to the exact build hash ensures reproducibility.
2. **Tool-author canonical images** (graph/alignment tools without bioconda recipes): `ghcr.io/pangenome/pggb:latest`, `quay.io/biocontainers/kegalign-gpu`, `ghcr.io/hillerlab/toga:latest`.

The pipeline ships a wrapper `pipeline/lib/run_in_container.sh` that picks the right image based on tool name + handles `--gpus`, `-v $WORK:$WORK`, and `-u $(id -u):$(id -g)` for the user-owned outputs.

| Phase | Tool | Container image | Notes |
|---|---|---|---|
| A | mash | `quay.io/biocontainers/mash:2.3--he348c14_4` | |
| A | busco | `ezlabgva/busco:v5.7.1_cv1` | docker hub, lineage downloaded once into `$WORK/work/00_inventory/busco/busco_downloads/` |
| B | longdust | `quay.io/biocontainers/longdust:0.1.0--h5b5514e_0` | bioconda — newly available; if missing on toolshed, fall back to lh3/longdust source build |
| B | sdust | `quay.io/biocontainers/seqkit-sdust:2.8.0--h9ee0642_0` | seqkit's sdust binding; or use `minimap2` image which bundles sdust |
| B | bedtools, samtools | `quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2`, `quay.io/biocontainers/samtools:1.20--h50ea8bc_0` | |
| C | kegalign-gpu | `quay.io/biocontainers/kegalign-gpu:1.0.0--hdfd78af_0` | needs `--gpus device=$GPU`. CPU fallback: `quay.io/biocontainers/lastz:1.04.22--h0c08fa6_1` |
| C | UCSC kentUtils (axtChain, chainNet, ...) | `quay.io/biocontainers/ucsc-chainstitchid:469--h664eb37_0` (and sibling images per tool) | OR build a single multi-tool image; v3 used static binaries dropped into `tools/`, but for new species containerize. |
| C | liftoff | `quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0` | |
| C | agat | `quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0` | |
| C | TOGA2 / CESAR2 | `ghcr.io/hillerlab/toga:latest` (or build from clone — see TOGA repo) | heaviest install; the container is ~3 GB |
| C, F | gffread | `quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4` | |
| D | pggb (wfmash, seqwish, smoothxg, odgi) | `ghcr.io/pangenome/pggb:202412130311080800a17` | pin by digest, not `:latest`. The PGGB image bundles wfmash, seqwish, smoothxg, odgi. |
| F | mafft | `quay.io/biocontainers/mafft:7.526--h031d066_0` | |
| F | pal2nal | `quay.io/biocontainers/pal2nal:14.1--pl5321hdfd78af_5` | |
| F | trimal | `quay.io/biocontainers/trimal:1.4.1--h9948957_8` | |
| G | iqtree3 | `quay.io/biocontainers/iqtree:3.0.0--hdcf5f25_0` (binary named `iqtree3` inside) | If not yet on bioconda, build from https://github.com/iqtree/iqtree3 |
| H | hyphy | `quay.io/biocontainers/hyphy:2.5.62--he91c24d_0` | |
| I | multiz, axtToMaf | `quay.io/biocontainers/multiz:11.2--h470a237_0`, kentUtils image for axtToMaf | |
| J | bcftools, CrossMap | `quay.io/biocontainers/bcftools:1.20--h8b25389_0`, `quay.io/biocontainers/crossmap:0.6.5--pyh7cba7a3_0` | |

Host-side dependencies (NOT containerized):
- `docker` ≥24 with GPU runtime (`nvidia-container-toolkit`) if GPU used
- `bash` ≥4, `python3` ≥3.10 (for the Python helper scripts in `pipeline/scripts/` — `phase_c2_triage.py`, `build_msa.py`, etc. These import `pysam`, `pandas`, `numpy` — but they run inside the bcftools or pggb image which already has them. Run via `run_in_container.sh python3 ...`.)
- `gh` (GitHub CLI) — only for issue filing, not part of the pipeline proper

### `pipeline/lib/run_in_container.sh` (wrapper)

```bash
#!/usr/bin/env bash
set -euo pipefail
# Usage: run_in_container.sh <tool_name> <args...>
TOOL=$1; shift
case $TOOL in
  mash)            IMG="quay.io/biocontainers/mash:2.3--he348c14_4" ;;
  busco)           IMG="ezlabgva/busco:v5.7.1_cv1" ;;
  longdust|sdust)  IMG="quay.io/biocontainers/longdust:0.1.0--h5b5514e_0" ;;
  bedtools)        IMG="quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2" ;;
  samtools)        IMG="quay.io/biocontainers/samtools:1.20--h50ea8bc_0" ;;
  kegalign)        IMG="quay.io/biocontainers/kegalign-gpu:1.0.0--hdfd78af_0" ;;
  axtChain|chainSort|chainPreNet|chainNet|netChainSubset|chainStitchId|chainSwap|axtToMaf) \
                   IMG="quay.io/biocontainers/ucsc-kent-tools:469--h664eb37_0" ;;
  liftoff)         IMG="quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0" ;;
  agat)            IMG="quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0" ;;
  toga|cesar)      IMG="ghcr.io/hillerlab/toga:latest" ;;
  gffread)         IMG="quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4" ;;
  pggb|wfmash|seqwish|smoothxg|odgi) \
                   IMG="ghcr.io/pangenome/pggb:202412130311080800a17" ;;
  mafft)           IMG="quay.io/biocontainers/mafft:7.526--h031d066_0" ;;
  pal2nal)         IMG="quay.io/biocontainers/pal2nal:14.1--pl5321hdfd78af_5" ;;
  trimal)          IMG="quay.io/biocontainers/trimal:1.4.1--h9948957_8" ;;
  iqtree3)         IMG="quay.io/biocontainers/iqtree:3.0.0--hdcf5f25_0" ;;
  hyphy)           IMG="quay.io/biocontainers/hyphy:2.5.62--he91c24d_0" ;;
  multiz)          IMG="quay.io/biocontainers/multiz:11.2--h470a237_0" ;;
  bcftools)        IMG="quay.io/biocontainers/bcftools:1.20--h8b25389_0" ;;
  CrossMap)        IMG="quay.io/biocontainers/crossmap:0.6.5--pyh7cba7a3_0" ;;
  python3)         IMG="quay.io/biocontainers/bcftools:1.20--h8b25389_0" ;;  # bcftools image has python3 + pysam + pandas
  *) echo "unknown tool: $TOOL" >&2; exit 1 ;;
esac
GPU_ARG=""
[[ ${USE_GPU:-0} == 1 ]] && GPU_ARG="--gpus device=${GPU:-0}"
exec docker run --rm $GPU_ARG \
  -v "$WORK:$WORK" -v "${SCRATCH:-/tmp}:${SCRATCH:-/tmp}" \
  -w "$WORK" \
  -u "$(id -u):$(id -g)" \
  --entrypoint "" \
  "$IMG" \
  $TOOL "$@"
```

All phase script invocations below use `$WORK/pipeline/lib/run_in_container.sh <tool>` instead of bare tool names. The script defaults the user to `$(id -u):$(id -g)` so outputs aren't root-owned (this was a real pain point in v3 — gffutils sqlite DBs and OrthoFinder outputs ended up root-owned and required sudo cleanup).

### Singularity / Apptainer alternative

If `docker` is unavailable (HPC clusters typically forbid it), the same images run under Singularity:
```bash
singularity exec --bind $WORK:$WORK --nv docker://quay.io/biocontainers/kegalign-gpu:1.0.0--hdfd78af_0 \
  kegalign $WORK/genomes/softmasked/$A.fa $WORK/genomes/softmasked/$B.fa ...
```
Ship a sibling `pipeline/lib/run_in_apptainer.sh` with the same `TOOL → IMG` mapping but using `singularity exec` instead of `docker run`. Phase scripts source `$RUN_IN_CONTAINER` (env var) to pick which wrapper to use.

## species.conf — template

```bash
# species.conf
SPECIES=Pv4
WORK=/media/anton/data/sandbox/Pv4/v3
STRAINS=(PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087)
REF_STRAIN=PvP01
ANCHOR_STRAINS=(PvP01 Sal-I PvW1 PAM PvSY56)
COHORT_VCF_DIR=/media/anton/scratch/malariagen_pv4
COHORT_CHROM_GLOB='Pv4_PvP01_*_v1.vcf.gz'
CHROM_RENAME=$WORK/inputs/annotations/PvP01_plasmodb_to_genbank.tsv
N_CORES=32
GPU=0
DOCKER_GPU="--gpus device=$GPU"
MIN_INTACT_STRICT=7
MIN_INTACT_RELAXED=5
VAR_ANTIGEN_RE='PIR|PHIST|Pv-fam|MSP|DBP|EBA|RBP|AMA|RAP|SERA|TRAg|STP1|RESA'
# Smoke-test subset (optional)
SMOKE_CHROM=LT635612.2
SMOKE_STRAINS=(PvP01 Sal-I PvW1)
```

Validation step (before any phase): the orchestrator MUST call `bash pipeline/00_validate_conf.sh` which checks every `${STRAINS[@]}` has a readable `inputs/assemblies/{S}.fa` and `inputs/annotations/{S}.gff3`, that `REF_STRAIN` is in `STRAINS`, that `CHROM_RENAME` is a tab-separated 2-column file, and that `N_CORES` ≤ `$(nproc)`.

---

## 1. Phase A — Inventory (mash + BUSCO)

### 1.1 Mash N×N distance matrix

```
mash sketch -p $N_CORES -k 21 -s 10000 -o $WORK/work/00_inventory/mash/sketch \
  ${STRAINS[@]/#/$WORK/inputs/assemblies/}.fa
mash dist -t $WORK/work/00_inventory/mash/sketch.msh $WORK/work/00_inventory/mash/sketch.msh \
  > $WORK/work/00_inventory/mash/dist.tsv
```

- ⭐ Output: `work/00_inventory/mash/dist.tsv` (N+1 rows, tab-separated, first row is column names with strain assemblies, first column is row names)
- Guard: `[[ -s $WORK/work/00_inventory/mash/dist.tsv ]] || { mash sketch ...; mash dist ...; }`
- Validation: `awk -v n=${#STRAINS[@]} 'NR==1 && NF==n+1 {ok=1} END {exit !ok}' $WORK/work/00_inventory/mash/dist.tsv`
- Gotchas:
  - `mash sketch -s 10000` (sketch size) is critical; the default `-s 1000` is too coarse for ~25 Mb apicomplexan genomes. v3 used 10,000.
  - The braced expansion `${STRAINS[@]/#/...}.fa` prepends path to each strain name and appends `.fa`. Verify shell version ≥ bash 4.0.
- Wall time: ~1 min for 8 strains × 25 Mb.

### 1.2 BUSCO completeness per strain

```
for S in "${STRAINS[@]}"; do
  out=$WORK/work/00_inventory/busco/${S}_proteins
  [[ -d $out ]] && continue
  # Run BUSCO on the proteome (not the genome — much faster for gene-dense apicomplexans)
  busco -i $WORK/inputs/proteomes/${S}.proteins.fa \
        -l plasmodium_odb10 \
        -m prot \
        -o ${S}_proteins \
        --out_path $WORK/work/00_inventory/busco/ \
        --download_path $WORK/work/00_inventory/busco/busco_downloads/ \
        --cpu $N_CORES
done
```

- Output (per strain): `work/00_inventory/busco/{S}_proteins/short_summary.specific.plasmodium_odb10.{S}_proteins.txt`
- Guard: per-strain `[[ -d ... ]] && continue` inside loop above.
- Validation: `grep -q "^	C:[0-9]\+\.[0-9]\+%" $out/short_summary.*.txt`
- Gotchas:
  - Lineage `plasmodium_odb10` is species-specific. For *P. falciparum* use the same. For non-Plasmodium, change to the appropriate odb10 (use `busco --list-datasets`).
  - The `--download_path` is required if running offline; without it BUSCO tries to redownload every invocation.
  - Outputs are owned by uid 2000 inside the docker container — beware later cleanup needs `sudo`.

---

## 2. Phase B — Soft-masking

For each strain `S` in `STRAINS`:

### 2.1 Identify low-complexity regions (union of longdust + sdust)

```
mkdir -p $WORK/genomes/mask
longdust -t $N_CORES $WORK/inputs/assemblies/${S}.fa | \
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' > $WORK/genomes/mask/${S}.longdust.bed
sdust $WORK/inputs/assemblies/${S}.fa > $WORK/genomes/mask/${S}.sdust.bed
cat $WORK/genomes/mask/${S}.longdust.bed $WORK/genomes/mask/${S}.sdust.bed | \
  bedtools sort -i - | bedtools merge -i - > $WORK/genomes/mask/${S}.union.bed
```

- Intermediate outputs: `genomes/mask/{S}.{longdust,sdust,union}.bed`
- Guard per-strain: `[[ -s $WORK/genomes/mask/${S}.union.bed ]] || { longdust ...; sdust ...; bedtools merge ...; }`
- Validation: `[[ -s $WORK/genomes/mask/${S}.union.bed ]] && awk '$3 > $2' $WORK/genomes/mask/${S}.union.bed | wc -l | (read n; [[ $n -gt 0 ]])`

### 2.2 Apply soft mask

```
bedtools maskfasta -soft \
  -fi $WORK/inputs/assemblies/${S}.fa \
  -bed $WORK/genomes/mask/${S}.union.bed \
  -fo $WORK/genomes/softmasked/${S}.fa
samtools faidx $WORK/genomes/softmasked/${S}.fa
cut -f1,2 $WORK/genomes/softmasked/${S}.fa.fai > $WORK/genomes/softmasked/${S}.sizes
```

- ⭐ Output (per strain): `genomes/softmasked/{S}.fa` + `.fa.fai` + `.sizes`
- Guard: `[[ -f $WORK/genomes/softmasked/${S}.sizes ]] || { bedtools maskfasta ...; samtools faidx ...; cut -f1,2 ... ; }`
- Validation: `samtools quickcheck -v $WORK/genomes/softmasked/${S}.fa`
- Gotchas:
  - `-soft` lowercases the masked nucleotides; without it `-N` would mask hard (replace with N) and break downstream kmer tools.
  - The `.sizes` file (2-col tab: name, length) is required by `chainNet` and `chainPreNet` later — bake it now.
- Wall time: ~30s per 25 Mb strain, ~5 min total for 8 strains.

---

## 3. Phase C — Pairwise alignment + chain pipeline

### 3.1 KegAlign GPU pairwise lastZ — N×(N-1)/2 unordered pairs

```
mkdir -p $WORK/projection/A2_kegalign/axt
for i in "${!STRAINS[@]}"; do
  for j in "${!STRAINS[@]}"; do
    [[ $i -ge $j ]] && continue
    A=${STRAINS[i]}; B=${STRAINS[j]}
    out=$WORK/projection/A2_kegalign/axt/${A}__vs__${B}.axt
    [[ -s $out ]] && continue
    docker run --rm $DOCKER_GPU \
      -v $WORK:$WORK -w $WORK \
      quay.io/biocontainers/kegalign-gpu:latest \
      kegalign \
        $WORK/genomes/softmasked/${A}.fa \
        $WORK/genomes/softmasked/${B}.fa \
        --strand both \
        --hsp_threshold 5000 \
        --gapped_threshold 6000 \
        --inner 2000 \
        --ydrop 15000 \
        --output_format axt \
        > $out
  done
done
```

- ⭐ Output: `projection/A2_kegalign/axt/{A}__vs__{B}.axt` for every unordered pair (28 pairs for N=8)
- Guard: per-pair `[[ -s $out ]] && continue` in loop.
- Validation: `awk 'NR<5 && /^[0-9]+ /' $out | wc -l | (read n; [[ $n -ge 1 ]])` (first non-comment line must be a numbered alignment header)
- Gotchas:
  - `--strand both` produces both forward and reverse alignments in one AXT — chainSwap is still needed later for reciprocal-best, this just controls which target strand each HSP is reported against.
  - `--hsp_threshold 5000` is the v3 default; lower values produce many spurious short alignments that bloat chain output 10×.
  - GPU memory: A10G/A100 with ≥16 GB handles up to ~100 Mb each side. For larger genomes, split into chr-pair runs and concatenate.
- Wall time: ~3 min per pair on A10G; 28 pairs sequential ≈ 1.5 hrs. Use `xargs -P 1` (NOT parallel — KegAlign saturates the GPU).
- CPU fallback (no GPU): substitute `lastz $A.fa[multiple] $B.fa[multiple] --masking=50 --hspthresh=4500 --gappedthresh=6000 --inner=2000 --ydrop=15000 --format=axt > $out`. ~10× slower.

### 3.2 Chain pipeline — cleaned chains (canonical UCSC chain build)

For each unordered pair `{A, B}` produce both directional chains (`A→B` and `B→A`):

```
for AXT in $WORK/projection/A2_kegalign/axt/*.axt; do
  pair=$(basename $AXT .axt)
  A=${pair%%__vs__*}; B=${pair##*__vs__}
  for dir in "A_to_B" "B_to_A"; do
    if [[ $dir == "A_to_B" ]]; then SRC=$A; TGT=$B; AXT_DIRECTED=$AXT
    else SRC=$B; TGT=$A; AXT_DIRECTED=$WORK/projection/A2_kegalign/axt/${B}__vs__${A}.axt
    # If reverse-directed axt does not exist, use chainSwap below instead.
    fi
    OUT=$WORK/work/01_chains/${SRC}.${TGT}.cleaned.chain
    [[ -s $OUT ]] && continue
    $WORK/tools/axtChain -linearGap=loose $AXT_DIRECTED \
        $WORK/genomes/softmasked/${SRC}.fa \
        $WORK/genomes/softmasked/${TGT}.fa \
        /dev/stdout | \
      $WORK/tools/chainSort stdin /dev/stdout | \
      $WORK/tools/chainPreNet stdin \
        $WORK/genomes/softmasked/${SRC}.sizes \
        $WORK/genomes/softmasked/${TGT}.sizes /dev/stdout | \
      $WORK/tools/chainNet stdin \
        $WORK/genomes/softmasked/${SRC}.sizes \
        $WORK/genomes/softmasked/${TGT}.sizes \
        /tmp/${SRC}_to_${TGT}.net /dev/null | \
      $WORK/tools/netChainSubset /tmp/${SRC}_to_${TGT}.net /dev/stdin /dev/stdout | \
      $WORK/tools/chainStitchId /dev/stdin $OUT
    rm -f /tmp/${SRC}_to_${TGT}.net
  done
done
```

- ⭐ Output: `work/01_chains/{src}.{tgt}.cleaned.chain` (56 files for N=8: 28 pairs × 2 directions)
- Guard: per-direction `[[ -s $OUT ]] && continue`.
- Validation: `head -1 $OUT | awk '$1=="chain" && NF==13 {print}' | wc -l | (read n; [[ $n -eq 1 ]])`
- Gotchas:
  - The 7-stage pipeline `axtChain → chainSort → chainPreNet → chainNet → netChainSubset → chainStitchId` is the **canonical UCSC chain build**, NOT optional. Skipping `chainPreNet` or `netChainSubset` produces ~10× more chains and breaks `CrossMap` projection in Phase J because tiny overlapping chains get picked first.
  - `axtChain -linearGap=loose` is correct for inter-species (≥80% divergence). For intra-species use `-linearGap=medium`. Apicomplexan strains are intra-species; v3 used `loose` because the 8 strains span >5% divergence in some regions.
  - The `.sizes` files MUST come from Phase B; do NOT use the `.fa.fai` directly (different column layout — `.fa.fai` has 5 columns, sizes file has 2).
  - `chainNet`'s second positional output (the target-net) is unused for our needs; redirect to `/dev/null`.
  - `chainStitchId` is the final step — it renumbers chain IDs to be unique. Without it, downstream `CrossMap` can pick the wrong chain for ambiguous regions.

### 3.3 Chain pipeline — rbest chains (reciprocal best)

For each unordered pair, the rbest chain is the intersection of A→B-best and B→A-best:

```
for i in "${!STRAINS[@]}"; do
  for j in "${!STRAINS[@]}"; do
    [[ $i -ge $j ]] && continue
    A=${STRAINS[i]}; B=${STRAINS[j]}
    OUT=$WORK/work/01_chains/${A}.${B}.rbest.chain
    [[ -s $OUT ]] && continue
    # Step 1: swap A→B chain to be B→A (i.e., source/target roles flipped)
    $WORK/tools/chainSwap $WORK/work/01_chains/${A}.${B}.cleaned.chain /dev/stdout | \
      $WORK/tools/chainSort stdin $WORK/work/01_chains/${A}.${B}.swap.chain
    # Step 2: chainNet on the swapped chain → "best B-target alignments"
    $WORK/tools/chainNet $WORK/work/01_chains/${A}.${B}.swap.chain \
        $WORK/genomes/softmasked/${B}.sizes \
        $WORK/genomes/softmasked/${A}.sizes \
        /tmp/B_to_A_net /dev/null
    $WORK/tools/netChainSubset /tmp/B_to_A_net \
        $WORK/work/01_chains/${A}.${B}.swap.chain /dev/stdout | \
      $WORK/tools/chainStitchId /dev/stdin $WORK/work/01_chains/${A}.${B}.swap.cleaned.chain
    # Step 3: swap back, retain only chains that survived both nets
    $WORK/tools/chainSwap $WORK/work/01_chains/${A}.${B}.swap.cleaned.chain /dev/stdout | \
      $WORK/tools/chainSort stdin $OUT
    # Cleanup
    rm -f $WORK/work/01_chains/${A}.${B}.swap*.chain /tmp/B_to_A_net
  done
done
```

- ⭐ Output: `work/01_chains/{A}.{B}.rbest.chain` (28 files for N=8)
- Guard: `[[ -s $OUT ]] && continue`.
- Validation: `[[ -s $OUT ]] && head -1 $OUT | awk '$1=="chain"' | wc -l | (read n; [[ $n -eq 1 ]])`
- Gotcha: rbest chains have ~30-50% fewer alignments than cleaned chains. This is expected — rbest is more conservative.

---

## 4. Phase C (cont.) — Annotation projection (Liftoff + TOGA2/CESAR2 + merge)

For each anchor strain `A` in `ANCHOR_STRAINS`, project the anchor's annotation onto every other strain `Q`:

### 4.1 Liftoff projection

```
mkdir -p $WORK/work/02a_liftoff/${A}-as-ref
cat > /tmp/feature_types.txt <<EOF
protein_coding_gene
ncRNA_gene
pseudogene
EOF
for Q in "${STRAINS[@]}"; do
  [[ $Q == $A ]] && continue
  OUT=$WORK/work/02a_liftoff/${A}-as-ref/${Q}.liftoff.gff3
  [[ -s $OUT ]] && continue
  liftoff \
    -g $WORK/inputs/annotations/${A}.fixed.gff3 \
    -f /tmp/feature_types.txt \
    -o $OUT \
    -dir $WORK/work/02a_liftoff/${A}-as-ref/${Q}_intermediate \
    -p $N_CORES \
    -copies -sc 0.95 \
    $WORK/inputs/assemblies/${Q}.fa \
    $WORK/inputs/assemblies/${A}.fa
done
```

- Intermediate output: `work/02a_liftoff/{A}-as-ref/{Q}.liftoff.gff3`
- Guard: per-pair `[[ -s $OUT ]] && continue`.
- Validation: `[[ $(awk -F'\t' '$3=="gene"' $OUT | wc -l) -gt 100 ]]` (sanity: at least 100 genes lifted)
- Gotchas:
  - `-f feature_types.txt` is REQUIRED. Without it Liftoff lifts everything including `region`, `chromosome`, `assembly_unit` — these are not genes and break Phase E.
  - The first positional FASTA is the TARGET (`${Q}.fa`), the second is the REFERENCE (`${A}.fa`) — easy to invert and produce nonsense.
  - `-copies -sc 0.95` enables gene-copy detection (paralog finder). v3 used this to identify duplicated PIR/PHIST.
  - `${A}.fixed.gff3` is the chromosome-renamed version. If `${A}.gff3` uses PlasmoDB chromnames but `${A}.fa` uses GenBank LT-numbers, run `bash pipeline/lib/fix_gff_chroms.sh ${A}` first. The script reads the FASTA, builds a renaming table, and rewrites the GFF.

### 4.2 Triage (8-rule classification)

```
for Q in "${STRAINS[@]}"; do
  [[ $Q == $A ]] && continue
  OUT_DIR=$WORK/work/02b_triage/${A}-as-ref/${Q}
  [[ -s $OUT_DIR/needs_cesar2.bed ]] && continue
  mkdir -p $OUT_DIR
  python3 $WORK/pipeline/scripts/phase_c2_triage.py \
    --liftoff $WORK/work/02a_liftoff/${A}-as-ref/${Q}.liftoff.gff3 \
    --ref_gff $WORK/inputs/annotations/${A}.fixed.gff3 \
    --output_dir $OUT_DIR \
    --rules clean,low_cov,partial,frameshift,split,short_exon,extra_copy,truncated
done
```

- Output per pair: `work/02b_triage/{A}-as-ref/{Q}/{needs_cesar2.bed, liftoff_clean.gff3, triage.tsv, summary.json}`
- Guard: per-pair `[[ -s $OUT_DIR/needs_cesar2.bed ]]`.
- The 8 rules (from `phase_c2_triage.py` in v3):
  1. `clean` — Liftoff identity ≥ 0.95, full-length, no frameshift → into `liftoff_clean.gff3`
  2. `low_cov` — Liftoff coverage < 0.5 → into `needs_cesar2.bed` (CESAR will retry)
  3. `partial` — Liftoff full but ≤ 0.95 identity → into `needs_cesar2.bed`
  4. `frameshift` — CDS phase inconsistent → into `needs_cesar2.bed`
  5. `split` — gene lifted to multiple contigs → flagged in `triage.tsv`, not in CESAR queue
  6. `short_exon` — internal exon < 20 bp post-projection → into `needs_cesar2.bed`
  7. `extra_copy` — additional Liftoff copy in target with > 0.85 identity → into `triage.tsv` as `extra_copy`
  8. `truncated` — > 30% of CDS missing → into `needs_cesar2.bed`

### 4.3 TOGA2 / CESAR2 rescue (for genes in `needs_cesar2.bed`)

```
for Q in "${STRAINS[@]}"; do
  [[ $Q == $A ]] && continue
  OUT_DIR=$WORK/work/02c_toga/${A}-as-ref/${Q}
  [[ -s $OUT_DIR/annotation.gff3 ]] && continue
  mkdir -p $OUT_DIR
  python3 /opt/TOGA/toga.py \
    $WORK/inputs/annotations/${A}.bed12 \
    $WORK/work/01_chains/${A}.${Q}.cleaned.chain \
    $WORK/genomes/softmasked/${A}.fa $WORK/genomes/softmasked/${Q}.fa \
    --pn $OUT_DIR \
    --cb $WORK/inputs/annotations/${A}.isoforms.tsv \
    --u12 "" \
    --kt \
    --filter_bed $WORK/work/02b_triage/${A}-as-ref/${Q}/needs_cesar2.bed \
    --nc $N_CORES
done
```

- Output per pair: `work/02c_toga/{A}-as-ref/{Q}/{annotation.gff3, classification.tsv}`
- Guard: per-pair `[[ -s $OUT_DIR/annotation.gff3 ]]`.
- Gotchas:
  - TOGA2 requires `$A.bed12` (UCSC BED12 of CDS) + `$A.isoforms.tsv` (gene-id ↔ transcript-id mapping). Build these once per anchor via `scripts/build_anchor_inputs.sh ${A}`.
  - The `.chain` MUST be the cleaned chain (Phase C.2 output), not the AXT — TOGA2 uses chain segments to define CDS-aware projections.
  - `--filter_bed` restricts CESAR to only the genes Liftoff didn't handle cleanly. Without it, CESAR re-projects everything (slow, ~6 hrs per anchor instead of ~2 hrs).
  - TOGA2 writes some files as root inside the docker container. v3 ran TOGA2 directly (not in docker) to avoid this.

### 4.4 Merge Liftoff_clean + TOGA2 outputs

```
mkdir -p $WORK/work/02d_merged/${A}-as-ref
for Q in "${STRAINS[@]}"; do
  [[ $Q == $A ]] && continue
  OUT_GFF=$WORK/work/02d_merged/${A}-as-ref/${Q}.annotation.gff3
  OUT_TSV=$WORK/work/02d_merged/${A}-as-ref/${Q}.classification.tsv
  [[ -s $OUT_GFF && -s $OUT_TSV ]] && continue
  python3 $WORK/pipeline/scripts/phase_c4_merge.py \
    --liftoff_clean $WORK/work/02b_triage/${A}-as-ref/${Q}/liftoff_clean.gff3 \
    --toga_gff $WORK/work/02c_toga/${A}-as-ref/${Q}/annotation.gff3 \
    --toga_class $WORK/work/02c_toga/${A}-as-ref/${Q}/classification.tsv \
    --triage $WORK/work/02b_triage/${A}-as-ref/${Q}/triage.tsv \
    --output_gff $OUT_GFF \
    --output_class $OUT_TSV
done
```

- ⭐ Output: `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` and `.classification.tsv` (one of each per non-anchor strain per anchor — for N=8 with 5 anchors, 5×7 = 35 of each)
- Guard: per-pair both files non-empty.
- Validation: `[[ $(awk -F'\t' '$3=="gene"' $OUT_GFF | wc -l) -ge $(awk -F'\t' '$3=="gene"' $WORK/inputs/annotations/${A}.fixed.gff3 | wc -l | awk '{print int($1*0.7)}') ]]` (sanity: at least 70% of anchor genes have a projected counterpart)
- Gotchas:
  - The merge favors Liftoff for `clean` genes and TOGA2 for everything in `needs_cesar2.bed`. For `extra_copy` rows the merge keeps both calls (one main + one labelled `extra_copy=true` in `classification.tsv`).
  - The merge MUST also tag each gene in classification.tsv with one of: `LIFTOFF_CLEAN, CESAR_RESCUE, CESAR_PARTIAL, MISSING, SPLIT, EXTRA_COPY`. Phase E reads this column to decide consensus.

---

## 5. Phase D — PGGB pangenome graph

### 5.1 PanSN rename + concat

```
mkdir -p $WORK/work/pggb_in
for S in "${STRAINS[@]}"; do
  out=$WORK/work/pggb_in/${S}_pansn.fa
  [[ -s $out ]] && continue
  python3 $WORK/pipeline/scripts/pansn_rename.py \
    $WORK/genomes/softmasked/${S}.fa $S \
    > $out
done
# Concatenate and bgzip
cat $WORK/work/pggb_in/*_pansn.fa | bgzip -c > $WORK/work/pggb_in/all_pansn.fa.gz
samtools faidx $WORK/work/pggb_in/all_pansn.fa.gz
```

The `pansn_rename.py` script (16 lines):
```python
import sys
in_fa, sample = sys.argv[1:]
for line in open(in_fa):
    if line.startswith('>'):
        contig = line[1:].strip().split()[0]
        print(f'>{sample}#1#{contig}')
    else:
        print(line, end='')
```

- Output: `work/pggb_in/all_pansn.fa.gz` (+ `.fai`, `.gzi`)
- Guard: `[[ -s $WORK/work/pggb_in/all_pansn.fa.gz.gzi ]] || { ... }`
- Gotcha: PanSN format is `SAMPLE#HAPLOTYPE#CONTIG` (haploid → haplotype always `1`). Without PanSN names, `odgi paths` later picks arbitrary names and breaks ortho-by-path-co-membership queries in Phase E.

### 5.2 PGGB build

```
docker run --rm -v $WORK:$WORK -w $WORK ghcr.io/pangenome/pggb:latest \
  pggb -i $WORK/work/pggb_in/all_pansn.fa.gz \
       -n ${#STRAINS[@]} \
       -s 5000 -p 90 -k 23 \
       -t $N_CORES \
       -o $WORK/work/pggb_out/
# Symlink canonical names
ln -sf $(ls $WORK/work/pggb_out/*.smooth.fix.gfa | head -1) $WORK/inputs/pggb/pv.gfa
ln -sf $(ls $WORK/work/pggb_out/*.smooth.fix.og  | head -1) $WORK/inputs/pggb/pv.og
```

- ⭐ Output: `inputs/pggb/pv.gfa`, `inputs/pggb/pv.og`
- Guard: `[[ -L $WORK/inputs/pggb/pv.og && -s $WORK/inputs/pggb/pv.og ]] || { docker ... ; ln ... ; }`
- Validation: `docker run --rm -v $WORK:$WORK ghcr.io/pangenome/pggb:latest odgi stats -i $WORK/inputs/pggb/pv.og | grep -q "^#length"`
- Gotchas:
  - `-n` MUST equal `${#STRAINS[@]}`. PGGB uses it to determine the all-vs-all wfmash mapping count.
  - `-p 90` (mapping identity) + `-s 5000` (block size) is suitable for intra-species apicomplexan panels. For more divergent panels reduce to `-p 80`.
  - `-k 23` is the smoothxg-haplotype block size; default is fine.
  - PGGB writes ~5–10 intermediate files into `-o`. The canonical end product is `*.smooth.fix.{og,gfa}`.
- Wall time: ~3–6 hours for 8 × 25 Mb genomes on 32 cores.

---

## 6. Phase E — Consensus orthology

Three independent ortholog calls combined by union-find on (gene, gene) edges with ≥90% reciprocal overlap.

### 6.1 Source 1: Liftoff projections (per anchor → all non-anchors)

Already produced in Phase C.4 (`work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3`). For each anchor gene `g_A` and each query strain `Q`, the projected gene ID in `Q` is the `Name=` attribute in the merged GFF.

### 6.2 Source 2: Chain-based reciprocal best (rbest chains from Phase C.3)

```
mkdir -p $WORK/work/03_consensus
python3 $WORK/pipeline/scripts/phase_e_rbest_overlap.py \
  --chains "$WORK/work/01_chains/*.rbest.chain" \
  --annotations "$WORK/inputs/annotations/*.bed" \
  --strains "${STRAINS[*]}" \
  --min_overlap 0.90 \
  --output $WORK/work/03_consensus/rbest_edges.tsv
```

`rbest_edges.tsv` columns: `strain_a, gene_a, strain_b, gene_b, overlap_a, overlap_b`.

### 6.3 Source 3: Graph path co-membership (PGGB)

```
docker run --rm -v $WORK:$WORK ghcr.io/pangenome/pggb:latest \
  odgi paths -i $WORK/inputs/pggb/pv.og --haplotypes \
  > $WORK/work/03_consensus/graph_paths.tsv
python3 $WORK/pipeline/scripts/phase_e_graph_edges.py \
  --paths $WORK/work/03_consensus/graph_paths.tsv \
  --annotations "$WORK/inputs/annotations/*.bed" \
  --strains "${STRAINS[*]}" \
  --output $WORK/work/03_consensus/graph_edges.tsv
```

### 6.4 Union-find consensus

```
python3 $WORK/pipeline/scripts/phase_e_consensus.py \
  --liftoff_dir $WORK/work/02d_merged \
  --rbest $WORK/work/03_consensus/rbest_edges.tsv \
  --graph $WORK/work/03_consensus/graph_edges.tsv \
  --anchors "${ANCHOR_STRAINS[*]}" \
  --strains "${STRAINS[*]}" \
  --ref $REF_STRAIN \
  --output $WORK/work/03_consensus/ortholog_table.tsv
```

- ⭐ Output: `work/03_consensus/ortholog_table.tsv` (14-column TSV: `orthogroup_id, label, n_strains, max_copies, ${STRAINS[@]}, graph_strains, graph_mean_pav`)
- Guard: `[[ $(wc -l < $WORK/work/03_consensus/ortholog_table.tsv) -gt 1000 ]]` (sanity: > 1000 orthogroups expected for apicomplexan panel)
- Validation: cells use `|` for cross-anchor aliases, `,` for paralog clusters, `-` for absent.
- The merge algorithm (from v3):
  1. Add a node per (strain, gene) seen in any annotation.
  2. Add an edge for every (g_A in anchor A, g_Q in Q) from Liftoff with `classification ∈ {LIFTOFF_CLEAN, CESAR_RESCUE}`.
  3. Add an edge for every rbest chain pair (strain_a.gene_a, strain_b.gene_b) with both overlaps ≥ 0.90.
  4. Add an edge for every graph-path co-membership pair (overlap ≥ 0.90 of CDS over a shared graph path).
  5. Union-find over edges → orthogroup IDs.
  6. For each orthogroup, count strains and max copies; emit one row.

---

## 7. Phase F — codon + protein MSAs (strict + relaxed sets)

For each `min_intact ∈ {MIN_INTACT_STRICT, MIN_INTACT_RELAXED}` and each orthogroup with ≥ `min_intact` 1-to-1 members:

```
for SETNAME in "strict:${MIN_INTACT_STRICT}" "relaxed:${MIN_INTACT_RELAXED}"; do
  NAME=${SETNAME%%:*}
  MIN=${SETNAME##*:}
  if [[ $NAME == "strict" ]]; then OUT_DIR=$WORK/work/06_msa/core_v3
  else OUT_DIR=$WORK/work/06_msa/core_relaxed; fi
  [[ -d $OUT_DIR && $(ls $OUT_DIR/*.codon.aln.fa 2>/dev/null | wc -l) -gt 100 ]] && continue
  mkdir -p $OUT_DIR
  python3 $WORK/pipeline/scripts/build_msa.py \
    --ortho $WORK/work/03_consensus/ortholog_table.tsv \
    --strains "${STRAINS[*]}" \
    --ref $REF_STRAIN \
    --min_intact $MIN \
    --cds_from_gff_fa \
    --aligner mafft-linsi \
    --backtranslate pal2nal \
    --threads $N_CORES \
    --output $OUT_DIR
done
```

`build_msa.py` for each orthogroup:
1. For each strain `S` in the orthogroup, extract CDS via `gffread -x /tmp/${S}.${og}.cds.fa -g $WORK/genomes/softmasked/${S}.fa $WORK/work/02d_merged/${anchor}-as-ref/${S}.annotation.gff3` filtered to gene_id = the strain's entry in the orthogroup.
2. Concat all strain CDSs into `/tmp/${og}.cds.fa`.
3. Translate via `gffread -y /tmp/${og}.aa.fa -g ...` (same source).
4. Protein MSA: `mafft --localpair --maxiterate 1000 --thread 2 /tmp/${og}.aa.fa > $OUT_DIR/${og}.protein.aln.fa` (LINSI mode).
5. Codon backtranslation: `pal2nal.pl $OUT_DIR/${og}.protein.aln.fa /tmp/${og}.cds.fa -output fasta > $OUT_DIR/${og}.codon.aln.fa`.
6. Validate: codon MSA length = 3 × protein MSA length; reject if not.

Cleaned variants in parallel:
- `trimal -in $OUT_DIR/${og}.codon.aln.fa -out $OUT_DIR/${og}.codon.cleaned.fa -automated1` → goes into the `core_v3_clean/` / `core_relaxed_clean/` sibling dir.

- ⭐ Output: `work/06_msa/core_v3/{PVP01_id}.codon.aln.fa` and `.protein.aln.fa` (strict: ~1,584 each); `work/06_msa/core_relaxed/{PVP01_id}.codon.aln.fa` and `.protein.aln.fa` (relaxed: ~4,215 each).
- Guard: per-orthogroup `[[ -s $OUT_DIR/${og}.codon.aln.fa && -s $OUT_DIR/${og}.protein.aln.fa ]] && continue`.
- Validation: `awk 'NR==2 {print length($0)}' $OUT_DIR/${og}.codon.aln.fa | (read n; [[ $((n % 3)) -eq 0 ]])`
- Gotchas:
  - MAFFT-LINSI (`--localpair --maxiterate 1000`) is ~10× slower than default `mafft` but produces alignments suitable for selection analysis. Default `mafft` will silently produce gappy garbage on divergent paralog families.
  - `pal2nal.pl` requires the CDS to be a multiple of 3 with no internal stops. If the CDS has a stop codon, `pal2nal` skips that strain — `build_msa.py` MUST strip internal stops first via `awk` (replace `taa|tag|tga` mid-CDS with `nnn`).
  - The orthogroup label `PVP01_NNNNNNN` is the v3 convention (reference strain's gene ID is used as the orthogroup name when 1-to-1 with anchor). For non-Plasmodium species, swap to whatever the reference annotation uses.
  - Run codon + protein MSAs in parallel: `ls $WORK/work/03_consensus/orthogroup_inputs/* | xargs -P $((N_CORES/2)) -I {} python3 build_msa.py --single_og {} ...`.
- Wall time: ~2 hrs for strict, ~5 hrs for relaxed (sharded across all cores).

---

## 8. Phase G — ML trees (IQ-TREE3)

```
for SET in core_v3 core_relaxed; do
  OUT_BASE=$WORK/work/06_msa/${SET}_trees
  mkdir -p $OUT_BASE
  ls $WORK/work/06_msa/$SET/*.codon.aln.fa | \
    xargs -P $((N_CORES/2)) -I {} bash -c '
      f={}; gene=$(basename $f .codon.aln.fa)
      out_dir='$OUT_BASE'/$gene
      [[ -s $out_dir/${gene}.treefile ]] && exit 0
      mkdir -p $out_dir
      # IQ-TREE3 fails with < 4 unique sequences; fall back to no-bootstrap
      n_unique=$(awk "/^>/{next} {print}" $f | sort -u | wc -l)
      if [[ $n_unique -ge 4 ]]; then
        iqtree3 -s $f -m MFP -B 1000 -T 2 -pre $out_dir/$gene 2>$out_dir/${gene}.log
      else
        iqtree3 -s $f -m MFP -T 2 -pre $out_dir/$gene 2>$out_dir/${gene}.log
      fi
    '
done
```

- ⭐ Output: `work/06_msa/{set}_trees/{gene}/{gene}.treefile` per gene.
- Guard: per-gene `[[ -s $out_dir/${gene}.treefile ]] && exit 0`.
- Validation: `[[ -s $out_dir/${gene}.treefile ]] && awk '/^[(]/' $out_dir/${gene}.treefile | wc -l | (read n; [[ $n -ge 1 ]])`
- Gotchas:
  - IQ-TREE3 fails hard with `< 4` unique sequences. The fallback (no `-B 1000`) is mandatory for variant-antigen genes where most strains have identical paralog copies.
  - `-T 2` (2 threads per tree) × `xargs -P $((N_CORES/2))` saturates without oversubscribing.
  - `-m MFP` triggers ModelFinder Plus (run-time + memory hit). For panels with > 10,000 genes, consider `-m GTR+I+G` to skip model selection (~3× faster).
  - Each tree's working dir contains the `.iqtree` report, `.ckp.gz` checkpoint, `.log`, `.bionj`, `.mldist`, `.contree`, `.splits.nex` — only `.treefile` is `*`-essential. The rest are ~5–10 MB each; for 5,800 genes that's 30+ GB. Tar+gzip the workdirs (Phase F summary section in OUTLINE).
- Wall time: ~2 hrs for strict, ~6 hrs for relaxed.

---

## 9. Phase H — HyPhy bulk BUSTED

```
for SET in core_v3 core_relaxed; do
  OUT_BASE=$WORK/work/06_msa/${SET}_hyphy
  [[ $SET == "core_v3" ]] && OUT_BASE=$WORK/work/06_msa/core_v3_hyphy/bulk
  mkdir -p $OUT_BASE
  ls $WORK/work/06_msa/$SET/*.codon.aln.fa | \
    xargs -P $((N_CORES)) -I {} bash -c '
      f={}; gene=$(basename $f .codon.aln.fa)
      out_dir='$OUT_BASE'/$gene
      [[ -s $out_dir/busted.json ]] && exit 0
      mkdir -p $out_dir
      tree='$WORK'/work/06_msa/'$SET'_trees/$gene/${gene}.treefile
      [[ ! -s $tree ]] && exit 0
      hyphy busted \
        --alignment $f \
        --tree $tree \
        --output $out_dir/busted.json \
        --srv No --branches All \
        > $out_dir/busted.log 2>&1
    '
done
```

- ⭐ Output: `work/06_msa/core_v3_hyphy/bulk/{gene}/busted.json` and `work/06_msa/core_relaxed_hyphy/{gene}/busted.json`.
- Guard: per-gene `[[ -s $out_dir/busted.json ]] && exit 0`.
- Validation: `python3 -c "import json,sys; d=json.load(open(sys.argv[1])); assert 'test results' in d, 'missing test results'" $out_dir/busted.json`
- Gotchas:
  - `--srv No` disables site-to-site rate variation (faster, more interpretable; matches v3 settings).
  - `--branches All` tests selection on the full tree; for branch-specific selection use `--branches Internal` or a labeled tree.
  - HyPhy BUSTED requires the tree from Phase G. If the tree is missing (per-gene IQ-TREE failed), HyPhy silently produces an empty `.json` — that's why the validation step parses the JSON.
  - For top-hit follow-up (MNS check), re-run with `--multiple-hits Double` and compare p-values. v3 found Pvs230 went p=0.056 → 0.500 with multi-hit modeling.
- Wall time: ~1.5 hrs for strict, ~4 hrs for relaxed (on 32 cores).

---

## 10. Phase I — Multiz multi-way MAFs

For each strain `H` (the hinge), build an N-way MAF anchored on `H`'s coordinates:

```
mkdir -p $WORK/work/07_multiz
for H in "${STRAINS[@]}"; do
  OUT_DIR=$WORK/work/07_multiz/${H}
  OUT_FINAL=$OUT_DIR/${H}.multiz.maf
  [[ -s $OUT_FINAL ]] && continue
  mkdir -p $OUT_DIR
  # Step 1: build all pairwise MAFs (H vs everyone else)
  for Q in "${STRAINS[@]}"; do
    [[ $Q == $H ]] && continue
    PAIR_MAF=$OUT_DIR/${H}_vs_${Q}.maf
    [[ -s $PAIR_MAF ]] && continue
    # Use the pair's AXT (either direction; pick the one where H is the source)
    if [[ -s $WORK/projection/A2_kegalign/axt/${H}__vs__${Q}.axt ]]; then
      INPUT_AXT=$WORK/projection/A2_kegalign/axt/${H}__vs__${Q}.axt
    else
      INPUT_AXT=$WORK/projection/A2_kegalign/axt/${Q}__vs__${H}.axt
      # Need to reverse roles via axtSwap — UCSC kentUtils has no axtSwap directly; use chainSwap on the chain
      # For our purposes, just use the existing axt with the pair correctly named
    fi
    $WORK/tools/axtToMaf -tPrefix=${H}. -qPrefix=${Q}. \
      $INPUT_AXT \
      $WORK/genomes/softmasked/${H}.sizes \
      $WORK/genomes/softmasked/${Q}.sizes \
      $PAIR_MAF
  done
  # Step 2: progressively fold all pairwise MAFs into a multi-way MAF rooted at H
  # multiz uses a guide tree (or pairwise distance); use mash distance from Phase A
  python3 $WORK/pipeline/scripts/multiz_progressive.py \
    --hinge $H \
    --pair_dir $OUT_DIR \
    --strains "${STRAINS[*]}" \
    --output $OUT_FINAL
done
```

- ⭐ Output: `work/07_multiz/{H}/{H}.multiz.maf` (one per strain — N total = 8 for v3)
- Guard: per-hinge `[[ -s $OUT_FINAL ]] && continue`.
- Validation: `awk '/^a / {n++} END {if (n<100) exit 1}' $OUT_FINAL` (≥ 100 alignment blocks).
- Gotchas:
  - `multiz` only handles pairwise MAFs as input. The `multiz_progressive.py` driver runs `multiz pair1.maf pair2.maf > tmp1.maf`, then `multiz tmp1.maf pair3.maf > tmp2.maf`, etc., growing the multi-way MAF one strain at a time.
  - The ORDER matters: closest strain first (use mash distances). Worst order can drop 30–40% of alignment blocks.
  - `-tPrefix=H.` / `-qPrefix=Q.` give MAF blocks the right PanSN-like names.
  - The intermediate `H_vs_Q.maf` files are kept (they're not `*`-essential but feed multiz).
- Wall time: ~3 hrs per hinge × N hinges = ~24 hrs total. Run hinges in parallel with `xargs -P 4` if you have 128 GB+ RAM.

---

## 11. Phase J — Cohort VCF projection (A2 chains)

### 11.1 Rename MalariaGEN-style cohort VCF to GenBank chromosome names

```
mkdir -p $WORK/projection/A1_wfmash/mg_renamed
for SRC_VCF in $COHORT_VCF_DIR/$COHORT_CHROM_GLOB; do
  CHR_TAG=$(basename $SRC_VCF | sed 's/.*_\([0-9]\+\|API\|MIT\).*/\1/')
  OUT=$WORK/projection/A1_wfmash/mg_renamed/${SPECIES}_${CHR_TAG}.vcf.gz
  [[ -s $OUT.csi ]] && continue
  bcftools annotate --rename-chrs $CHROM_RENAME $SRC_VCF -Oz -o $OUT
  bcftools index $OUT
done
```

- ⭐ Output: `projection/A1_wfmash/mg_renamed/{species}_{chr}.vcf.gz` (one per chromosome + API + MIT; 16 for Pv4)
- Guard: per-chr `[[ -s $OUT.csi ]]`.
- Validation: `bcftools view -h $OUT 2>/dev/null | grep -q "^##contig=<ID=LT" || bcftools view -h $OUT | grep -q "^##contig=<ID=$(head -1 $CHROM_RENAME | cut -f2)"`
- Gotcha: `--rename-chrs` reads a 2-col TSV (`OLD<TAB>NEW`). The map MUST cover every chrom in the source VCF or `bcftools` errors. Check with `comm -23 <(bcftools view -h $SRC_VCF | awk -F'[=,>]' '/##contig=<ID=/ {print $3}' | sort) <(cut -f1 $CHROM_RENAME | sort)` → empty result means full coverage.

### 11.2 CrossMap project per non-reference target

```
mkdir -p $WORK/projection/A2_lastz
TMP_SORT=/media/anton/scratch/path_a2_sort_tmp
mkdir -p $TMP_SORT

for TGT in "${STRAINS[@]}"; do
  [[ $TGT == $REF_STRAIN ]] && continue
  OUT_COHORT=$WORK/projection/A2_lastz/Pv4_cohort_on_${TGT}.vcf.gz
  [[ -s $OUT_COHORT.csi ]] && continue
  CHAIN=$WORK/work/01_chains/${REF_STRAIN}.${TGT}.cleaned.chain
  TGT_FA=$WORK/genomes/softmasked/${TGT}.fa
  PER_CHR_DIR=$WORK/projection/A2_lastz/per_chr/${TGT}
  mkdir -p $PER_CHR_DIR

  # Project each per-chr VCF
  for SRC_VCF in $WORK/projection/A1_wfmash/mg_renamed/${SPECIES}_*.vcf.gz; do
    CHR_TAG=$(basename $SRC_VCF .vcf.gz | sed "s/^${SPECIES}_//")
    OUT_CHR=$PER_CHR_DIR/${SPECIES}_${CHR_TAG}_on_${TGT}.vcf.gz
    [[ -s $OUT_CHR.csi ]] && continue
    # CrossMap writes unsorted; pipe through bcftools sort
    CrossMap vcf $CHAIN $SRC_VCF $TGT_FA /dev/stdout 2>/dev/null | \
      bcftools sort -T $TMP_SORT -Oz -o $OUT_CHR -
    bcftools index $OUT_CHR
  done
  # Concat per-chr to one cohort VCF
  bcftools concat -a $PER_CHR_DIR/*.vcf.gz -Oz -o $OUT_COHORT
  bcftools index $OUT_COHORT
done
```

- ⭐ Output: `projection/A2_lastz/Pv4_cohort_on_{TGT}.vcf.gz` (+ `.csi`) for every non-ref target (N-1 = 7 for v3)
- Guard: per-target `[[ -s $OUT_COHORT.csi ]] && continue`.
- Validation: `bcftools view -H $OUT_COHORT | wc -l | (read n; [[ $n -ge 100000 ]])` (sanity: ≥ 100k variants survive projection).
- Gotchas:
  - `CrossMap vcf` writes unsorted output (variants may be interleaved across chromosomes due to chain-spanning multi-segment liftover). `bcftools sort` is NOT optional — without it `bcftools index` fails.
  - The third positional arg to `CrossMap vcf` is the **TARGET** FASTA (provides REF for the lifted record). Easy to invert.
  - `bcftools concat -a` (allow duplicates) is required because per-chr VCFs may have overlap from inversion liftover.
  - The `--cleaned.chain` from Phase C.2 is the right input. The `.rbest.chain` is overly conservative and drops ~40% of liftable variants.
  - Sort temp dir must be on fast disk (NVMe or scratch); the default `/tmp` (often tmpfs) runs out of RAM for 25 GB sorts.
- Wall time: ~30 min per target × 7 = ~3.5 hrs (parallel with `xargs -P 4` if I/O permits).

---

## 12. Phase K — UCSC track-hub publishing

Builds a self-contained UCSC track hub from the upstream phase outputs (chains, multiz MAFs, merged annotations, HyPhy BUSTED JSONs, ortholog table). Output lives under `$WORK/ucsc_hub/` in a layout consumable directly by `genome.ucsc.edu` via `hubUrl=`.

### Inputs (all from prior phases)

| Source | Used by |
|---|---|
| `work/07_multiz/{hinge}/{hinge}.multiz.maf` | bigMaf |
| `work/01_chains/{src}.{tgt}.cleaned.chain` | chain track |
| `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` | gene BigBed |
| `work/06_msa/{core_v3,core_relaxed}_hyphy/.../busted.json` | selection BigBed12+5 |
| `work/03_consensus/ortholog_table.tsv` | orthogroup membership BigBed |
| `genomes/softmasked/{S}.fa.fai` | chrom sizes for `bedToBigBed` |
| `inputs/annotations/{anchor}.bed12` | reference gene coords for selection track |

### Additional tools (containerized like the rest)

The kentUtils image used in earlier phases bundles these; pin the same digest:

```bash
# image: quay.io/biocontainers/ucsc-kent-tools:469--h664eb37_0
cmd mafToBigMaf     # MAF → BED3+1 for bigMaf
cmd bedToBigBed     # BED → bigBed (used for bigMaf, BigBed12, BigBed12+5)
cmd mafIndex        # MAF positional index
cmd gff3ToGenePred  # GFF3 → genePred
cmd genePredToBed   # genePred → BED12
cmd hgFindSpec      # OPTIONAL: search-field config (not needed for v1)
```

AutoSql schemas. Save these to `$WORK/ucsc_hub/`:

`bigMaf.as` — fetch once from `https://genome.ucsc.edu/goldenPath/help/examples/bigMaf.as`.

`bigChain.as` + `bigLink.as` — canonical UCSC schemas for the bigChain track type. Save both verbatim from `tools/bigChain.as` and `tools/bigLink.as` (the chain track in trackDb is `type bigChain {TGT_ASSEMBLY}` and needs both `.bb` files: data via `bigDataUrl` and links via `linkDataUrl`).

`bigSelectionPlus5.as`:
```
table bigSelectionPlus5
"BUSTED selection result, one entry per orthogroup mapped to reference gene"
(
string  chrom;                "Reference chromosome"
uint    chromStart;           "0-based start"
uint    chromEnd;              "End"
string  name;                 "Orthogroup ID"
uint    score;                "0-1000, scaled from -log10(qvalue)"
char[1] strand;               "+ or -"
uint    thickStart;           "CDS start"
uint    thickEnd;             "CDS end"
uint    reserved;             "RGB color (use itemRgb)"
int     blockCount;
int[blockCount] blockSizes;
int[blockCount] chromStarts;
string  orthogroup_id;
uint    n_strains;
float   busted_pvalue;
float   busted_qvalue_fdr;
string  gene_family;
)
```

### K.1 bigMaf + mafIndex per hinge

```bash
mkdir -p $WORK/ucsc_hub
for HINGE in "${STRAINS[@]}"; do
  ACC=$(strain_to_accession $HINGE)   # helper: read from species.conf STRAIN_ACCESSIONS
  OUT_DIR=$WORK/ucsc_hub/$ACC
  mkdir -p $OUT_DIR
  [[ -s $OUT_DIR/${HINGE}.multiz.maf.bb ]] && continue
  # Local MAF (may be gzipped if Phase I gzip-piped to staging — handle both)
  MAF=$WORK/work/07_multiz/${HINGE}/${HINGE}.multiz.maf
  [[ -s ${MAF}.gz ]] && MAF=${MAF}.gz
  TMP=/tmp/${HINGE}.maf
  if [[ ${MAF} == *.gz ]]; then gunzip -c $MAF > $TMP; else cp $MAF $TMP; fi
  cmd mafToBigMaf $ACC $TMP /dev/stdout | sort -k1,1 -k2,2n | \
    cmd bedToBigBed -type=bed3+1 -as=$WORK/ucsc_hub/bigMaf.as -tab \
      stdin $WORK/genomes/softmasked/${HINGE}.fa.fai \
      $OUT_DIR/${HINGE}.multiz.maf.bb
  cmd mafIndex $TMP $OUT_DIR/${HINGE}.multiz.maf.bb.bai
  rm -f $TMP
done
```

- Output: `ucsc_hub/{ACC}/{HINGE}.multiz.maf.bb` + `.bai`
- Guard: `[[ -s $OUT_DIR/${HINGE}.multiz.maf.bb ]] && continue`
- Validation: `cmd bigBedSummary $OUT_DIR/${HINGE}.multiz.maf.bb | grep -q "items:"`
- Wall: ~10 min per hinge × 8 = ~80 min serial. Don't parallelize (~3 GB temp MAF per hinge × P would blow scratch).

### K.2 BigBed12 from merged anchor annotations

For each anchor `A` ∈ `ANCHOR_STRAINS` and target `Q` ∈ `STRAINS` (including `A` itself for the self-annotation case):

```bash
for A in "${ANCHOR_STRAINS[@]}"; do
  for Q in "${STRAINS[@]}"; do
    ACC_Q=$(strain_to_accession $Q)
    OUT=$WORK/ucsc_hub/$ACC_Q/annot_from_${A}.bb
    [[ -s $OUT ]] && continue
    # Source GFF: identity case uses inputs/, projection case uses 02d_merged
    if [[ $Q == $A ]]; then
      GFF=$WORK/inputs/annotations/${A}.fixed.gff3
    else
      GFF=$WORK/work/02d_merged/${A}-as-ref/${Q}.annotation.gff3
    fi
    [[ ! -s $GFF ]] && { log "missing GFF: $GFF"; continue; }
    cmd gff3ToGenePred $GFF /tmp/${A}_on_${Q}.gp
    cmd genePredToBed /tmp/${A}_on_${Q}.gp /tmp/${A}_on_${Q}.bed
    sort -k1,1 -k2,2n /tmp/${A}_on_${Q}.bed > /tmp/${A}_on_${Q}.sorted.bed
    cmd bedToBigBed -type=bed12 -tab /tmp/${A}_on_${Q}.sorted.bed \
      $WORK/genomes/softmasked/${Q}.fa.fai $OUT
    rm -f /tmp/${A}_on_${Q}.{gp,bed,sorted.bed}
  done
done
```

- Output (per anchor × target): `ucsc_hub/{ACC_Q}/annot_from_{A}.bb`
- For N=8 strains, ANCHOR_STRAINS=5: 5 × 8 = 40 BigBed files. Each ~few MB.
- Guard: per-pair `[[ -s $OUT ]] && continue`.
- Validation: `cmd bigBedInfo $OUT | grep -q "itemCount:"`

### K.3 BigBed12+5 selection track (BUSTED)

Per reference assembly (build PvP01 first; replicate to others later):

```bash
REF=$REF_STRAIN
ACC_REF=$(strain_to_accession $REF)
for SET in core_v3 core_relaxed; do
  OUT=$WORK/ucsc_hub/$ACC_REF/selection_${SET}.bb
  [[ -s $OUT ]] && continue
  cmd python3 $WORK/pipeline/scripts/build_selection_bigbed.py \
    --ortho       $WORK/work/03_consensus/ortholog_table.tsv \
    --hyphy_dir   $WORK/work/06_msa/${SET}_hyphy \
    --ref_bed12   $WORK/inputs/annotations/${REF}.bed12 \
    --ref_strain  $REF \
    --family_re   "$VAR_ANTIGEN_RE" \
    --output_bed  /tmp/selection_${SET}.bed
  sort -k1,1 -k2,2n /tmp/selection_${SET}.bed > /tmp/selection_${SET}.sorted.bed
  cmd bedToBigBed -type=bed12+5 -tab \
    -as=$WORK/ucsc_hub/bigSelectionPlus5.as \
    /tmp/selection_${SET}.sorted.bed \
    $WORK/genomes/softmasked/${REF}.fa.fai $OUT
  rm -f /tmp/selection_${SET}.{bed,sorted.bed}
done
```

The Python helper `build_selection_bigbed.py` (~150 lines):
1. Reads `ortholog_table.tsv`; for each OG with a 1:1 assignment in the reference strain, get the reference gene ID.
2. Reads `${REF}.bed12.gz`; build gene_id → BED12 row map.
3. For each OG: open `{hyphy_dir}/{og}/busted.json` (or `{hyphy_dir}/bulk/{og}/busted.json` for strict); parse `["test results"]["p-value"]`.
4. Compute BH-FDR q across all OGs.
5. For each OG with a reference gene match: emit BED12+5 with score=scaled(-log10(q)), itemRgb by q-bin (red q<0.01, orange q<0.05, yellow q<0.10, gray else), family from regex on the gene's PlasmoDB description, n_strains from ortholog table row.

- Output: `ucsc_hub/{ACC_REF}/selection_{set}.bb`
- Guard: per-set `[[ -s $OUT ]] && continue`.

### K.4 Orthogroup membership BigBed12

```bash
REF=$REF_STRAIN; ACC_REF=$(strain_to_accession $REF)
OUT=$WORK/ucsc_hub/$ACC_REF/orthogroup_membership.bb
[[ -s $OUT ]] || {
  cmd python3 $WORK/pipeline/scripts/build_orthogroup_bigbed.py \
    --ortho     $WORK/work/03_consensus/ortholog_table.tsv \
    --ref_bed12 $WORK/inputs/annotations/${REF}.bed12 \
    --ref_strain $REF \
    --output    /tmp/og_mem.bed
  sort -k1,1 -k2,2n /tmp/og_mem.bed > /tmp/og_mem.sorted.bed
  cmd bedToBigBed -type=bed12 -tab /tmp/og_mem.sorted.bed \
    $WORK/genomes/softmasked/${REF}.fa.fai $OUT
  rm -f /tmp/og_mem.{bed,sorted.bed}
}
```

Color by `n_strains`: 1 → dark red, 4 → orange, 8 → green.

### K.5 Hub manifests (`hub.txt`, `genomes.txt`, per-assembly `trackDb.txt`)

`genomes.txt` for an assembly hub needs the full 9-field record per assembly (`twoBitPath`, `defaultPos`, `groups`, `description`, `organism`, `scientificName`, `htmlPath`). Empty `defaultPos` or missing `.2bit` makes the hub fail to load.

```bash
cat > $WORK/ucsc_hub/hub.txt <<EOF
hub BRC_Pangenome_${SPECIES}_v1
shortLabel ${SPECIES} pangenome v1
longLabel  BRC ${SPECIES} ${#STRAINS[@]}-strain pangenome (v1) — alignments, annotations, selection
genomesFile genomes.txt
email anton@nekrut.org
descriptionUrl http://nekrut.org/${SPECIES,,}brc.html
EOF

> $WORK/ucsc_hub/genomes.txt
for S in "${STRAINS[@]}"; do
  ACC=$(strain_to_accession $S)
  DEFAULT_POS=$(default_pos_for $ACC)   # e.g. dhps locus for P. vivax: LT635625.2:1264700-1277700
  cat >> $WORK/ucsc_hub/genomes.txt <<EOF
genome $ACC
trackDb $ACC/trackDb.txt
groups $ACC/groups.txt
description ${SPECIES} ${S} (${ACC})
twoBitPath $ACC/${ACC}.2bit
organism ${SPECIES// /_}
defaultPos ${DEFAULT_POS}
scientificName ${SPECIES}
htmlPath $ACC/description.html

EOF
done

# 2bit links + per-assembly groups.txt
for S in "${STRAINS[@]}"; do
  ACC=$(strain_to_accession $S)
  ln -sf $WORK/projection/A2_kegalign/2bit/${ACC}.2bit $WORK/ucsc_hub/$ACC/${ACC}.2bit
  cp $WORK/pipeline/scripts/groups.txt $WORK/ucsc_hub/$ACC/groups.txt
done
```

For each `${ACC}/trackDb.txt`, emit one standalone bigMaf + 2 composites (3 on the reference strain). **bigMaf and bigChain cannot share a composite** — UCSC composites require members of one `type`.

```bash
for S in "${STRAINS[@]}"; do
  ACC=$(strain_to_accession $S)
  cmd python3 $WORK/pipeline/scripts/build_trackdb.py \
    --assembly $ACC \
    --strain   $S \
    --hub_dir  $WORK/ucsc_hub \
    --strains "${STRAINS[*]}" \
    --anchors "${ANCHOR_STRAINS[*]}" \
    --output $WORK/ucsc_hub/$ACC/trackDb.txt
done
```

`build_trackdb.py` per-assembly output (in this order):
- `track {name}_multiz` — standalone, `type bigMaf`
- `track brc_pangenome_chains` — composite, `type bigChain` with 7 sub-tracks (`type bigChain {TGT_ACC}` each, `bigDataUrl` + `linkDataUrl` to the bigChain.bb + bigChain.link.bb pair)
- `track brc_pangenome_annot` — composite, `type bigBed 12` with 4 sub-tracks (one per anchor's projection)
- `track brc_pangenome_select` — reference strain only, `type bigBed 12`, 3 sub-tracks (selection_strict, selection_relaxed, orthogroup_membership)

- Validation: `cmd hubCheck -level=warn http://localhost/.../hub.txt` (HubCheck is in the kentUtils image)

### K.6 Chain files → bigChain conversion

UCSC track hubs need indexed binary `bigChain` files, not gzipped chain text. Convert each `.chain.gz` to a `.bigChain.bb` + `.bigChain.link.bb` pair using the schemas in `tools/bigChain.as` and `tools/bigLink.as`.

```bash
# Schemas (write once)
cp $WORK/pipeline/scripts/bigChain.as $WORK/ucsc_hub/
cp $WORK/pipeline/scripts/bigLink.as  $WORK/ucsc_hub/

# Per-pair conversion (56 pairs for an 8-strain panel)
for S in "${STRAINS[@]}"; do
  ACC=$(strain_to_accession $S)
  mkdir -p $WORK/ucsc_hub/$ACC/chains
  for T in "${STRAINS[@]}"; do
    [[ $T == $S ]] && continue
    ACC_T=$(strain_to_accession $T)
    IN=$WORK/work/01_chains/${S}.${T}.cleaned.chain.gz
    OUT_BB=$WORK/ucsc_hub/$ACC/chains/${ACC}_to_${ACC_T}.bigChain.bb
    OUT_LK=$WORK/ucsc_hub/$ACC/chains/${ACC}_to_${ACC_T}.bigChain.link.bb
    SIZES=$WORK/work/01_chains/${ACC}.sizes
    [[ -s $OUT_BB && -s $OUT_LK ]] && continue
    # 1) parse chain → bigChain.bed + bigLink.bed (one pass each)
    cmd python3 $WORK/pipeline/scripts/chain_to_bigChain.py \
      $IN $WORK/tmp/${ACC}_to_${ACC_T}.bigChain.bed.raw \
          $WORK/tmp/${ACC}_to_${ACC_T}.bigLink.bed.raw
    sort -k1,1 -k2,2n $WORK/tmp/${ACC}_to_${ACC_T}.bigChain.bed.raw \
         > $WORK/tmp/${ACC}_to_${ACC_T}.bigChain.bed
    sort -k1,1 -k2,2n $WORK/tmp/${ACC}_to_${ACC_T}.bigLink.bed.raw  \
         > $WORK/tmp/${ACC}_to_${ACC_T}.bigLink.bed
    # 2) bedToBigBed for each
    cmd bedToBigBed -type=bed6+6 -as=$WORK/ucsc_hub/bigChain.as -tab \
         $WORK/tmp/${ACC}_to_${ACC_T}.bigChain.bed $SIZES $OUT_BB
    cmd bedToBigBed -type=bed4+1 -as=$WORK/ucsc_hub/bigLink.as -tab \
         $WORK/tmp/${ACC}_to_${ACC_T}.bigLink.bed  $SIZES $OUT_LK
  done
done

# Keep raw .chain.gz as download-only artifacts (not a track)
for S in "${STRAINS[@]}"; do
  ACC=$(strain_to_accession $S)
  for T in "${STRAINS[@]}"; do
    [[ $T == $S ]] && continue
    ACC_T=$(strain_to_accession $T)
    ln -sf $WORK/work/01_chains/${S}.${T}.cleaned.chain.gz \
           $WORK/ucsc_hub/$ACC/chains/${ACC}_to_${ACC_T}.chain.gz
  done
done
```

- Run-time: ~5 s per pair, ~5 min total for 56 pairs
- Output: 56 × `.bigChain.bb` (~9 MB total) + 56 × `.bigChain.link.bb` (~60 MB)
- The bigChain track in `trackDb.txt` references both: `bigDataUrl` → `.bigChain.bb`, `linkDataUrl` → `.bigChain.link.bb`

### K.7 hubCheck + push to datacache

```bash
# Local verify
cmd hubCheck -level=warn file://$WORK/ucsc_hub/hub.txt

# Push (rsync to datacache or rclone to S3-style endpoint)
rsync -avz $WORK/ucsc_hub/ user@hgdownload.soe.ucsc.edu:/mirror/hubs/BRC/pangenome_${SPECIES,,}_v1/
# OR
~/.local/bin/rclone copy $WORK/ucsc_hub remote:BRC_hubs/pangenome_${SPECIES,,}_v1/
```

### Wall time + outputs

| Step | Wall | Output size |
|---|---|---|
| K.1 bigMaf × 8 hinges | ~80 min | ~16 GB |
| K.2 BigBed12 × 40 (anchor×target) | ~5 min | ~150 MB |
| K.3 selection BigBed × 2 sets | ~3 min | ~10 MB |
| K.4 orthogroup BigBed | ~1 min | ~5 MB |
| K.5 manifests + 2bit links | ~30 s | text + symlinks |
| K.6 chain → bigChain × 56 | ~5 min | ~70 MB |
| K.7 hubCheck + push | ~3 min | — |
| **Total Phase K** | **~95 min** | **~16.3 GB** |

### Idempotency

Every step has `[[ -s $OUT ]] && continue`. Re-running Phase K resumes where left off. To force rebuild a specific track:
```bash
rm $WORK/ucsc_hub/{ACC}/selection_core_v3.bb
bash pipeline/12_ucsc_hub.sh
```

### Phase K outputs in the data model

Per the issue update, these map to `Pangenome.selection_tracks`, `multiz_alignments`, `merged_annotations`, and `chain_files` slots in `pangenomes.yml`. The hub is the user-facing surface; the catalog entry is the BRC-analytics-side bookkeeping.

---

## 13. Orchestrator (run_all.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname $0)/species.conf"
cd $WORK
mkdir -p logs

bash pipeline/00_validate_conf.sh
log() { echo "$(date +%H:%M:%S) $1" | tee -a logs/orchestrator.log; }

log "Phase A — inventory"
bash pipeline/01_inventory.sh > logs/01_inventory.log 2>&1
log "Phase B — mask"
bash pipeline/02_mask.sh > logs/02_mask.log 2>&1

# Phases C (align+chain) and C (annot) can run in parallel
log "Phase C.1-3 — align + chain"; bash pipeline/03_align_chain.sh > logs/03_align_chain.log 2>&1 &
PID_C=$!
log "Phase C.4 — annotation projection"; bash pipeline/04_annotate_project.sh > logs/04_annotate_project.log 2>&1 &
PID_ANNOT=$!
wait $PID_C $PID_ANNOT

# Phases D and E (E depends on C, can start as soon as C.1-3 done; D depends only on B)
log "Phase D — PGGB"; bash pipeline/05_pggb.sh > logs/05_pggb.log 2>&1 &
PID_D=$!
wait $PID_ANNOT  # E needs annotations
log "Phase E — consensus orthology"; bash pipeline/06_consensus.sh > logs/06_consensus.log 2>&1
wait $PID_D

# Phase F (depends on E + D + C.4)
log "Phase F — MSAs"
bash pipeline/07_msa.sh > logs/07_msa.log 2>&1

# Phase G and H in parallel
log "Phase G — IQ-TREE"; bash pipeline/08_trees.sh > logs/08_trees.log 2>&1 &
PID_G=$!
wait $PID_G   # H needs trees
log "Phase H — HyPhy"; bash pipeline/09_hyphy.sh > logs/09_hyphy.log 2>&1 &
PID_H=$!

# Phase I (depends on C.1; independent of E/F/G/H)
log "Phase I — Multiz"
bash pipeline/10_multiz.sh > logs/10_multiz.log 2>&1

# Phase J (depends on C.2)
log "Phase J — VCF projection"
bash pipeline/11_project_vcf.sh > logs/11_project_vcf.log 2>&1

wait $PID_H
log "ALL PHASES COMPLETE"
bash pipeline/lib/verify_essentials.sh
```

- Guard: each phase script is itself idempotent — re-running `run_all.sh` after a crash resumes where it left off.
- Final validation: `pipeline/lib/verify_essentials.sh` reads `writeup/OUTLINE.md`, extracts every `*`-marked path, and verifies all 27 outputs exist + pass structural validation. Exit 0 only if all 27 pass.

---

## 14. Idempotency summary

Every phase script has the same shape:
1. Source `species.conf`, `cd $WORK`.
2. Iterate over the unit of work (per-strain, per-pair, per-orthogroup).
3. Compute `OUT_PATH` for this unit.
4. `[[ -s $OUT_PATH ]] && continue` — skip if already done.
5. Run the command. Validate the output. If invalid, `rm -f $OUT_PATH; continue` (next attempt will redo).
6. After the loop, `verify_essentials.sh` confirms all `*`-marked outputs exist.

No phase mutates outputs of earlier phases. No phase writes outside `$WORK/{genomes,inputs,projection,work,writeup,logs}/`. Re-running any phase after partial completion either no-ops or finishes the remaining work.

---

## 15. Wall-time budget (32 cores, 1 GPU, 128 GB RAM, 1 TB NVMe scratch)

| Phase | Wall time | CPU-hrs |
|---|---|---|
| A inventory | 5 min | 0.5 |
| B mask | 5 min | 0.5 |
| C.1-3 align + chain | 1.5 hrs (GPU) + 10 min (chains) | 20 |
| C.4 annotation | 12 hrs (TOGA2 dominates) | 350 |
| D PGGB | 5 hrs | 160 |
| E consensus | 30 min | 4 |
| F MSAs | 5 hrs | 130 |
| G IQ-TREE | 7 hrs | 220 |
| H HyPhy | 5 hrs | 160 |
| I multiz | 24 hrs (sequential hinges) or 6 hrs (parallel) | 200 |
| J VCF projection | 3.5 hrs | 110 |
| K UCSC hub | 1.5 hrs | 10 |
| **Total** | **~25.5 hrs wall**, ~1,410 CPU-hrs | |

Critical path: C.4 (annotation) and I (multiz) dominate. Parallelize C.4 across anchors and I across hinges for ~30% wall-time reduction. K is I/O-bound (`mafToBigMaf` on the multi-GB MAFs); leave serial.

---

## 16. Smoke test (smoke_test.sh)

For a 30-minute validation on new species, restrict to `${SMOKE_STRAINS[@]}` (subset of 3) and one chromosome `${SMOKE_CHROM}`. Smoke test runs phases A → B → C.1-3 → C.4 → D → E → F → G → H → J (skips I — multiz needs all strains). Expected outputs: ~50 orthogroup MSAs, ~50 trees, ~50 BUSTED jsons, ~3,000 lifted variants per target.

```bash
#!/usr/bin/env bash
set -euo pipefail
source "$(dirname $0)/species.conf"
cd $WORK
SAVED_STRAINS=("${STRAINS[@]}")
SAVED_ANCHORS=("${ANCHOR_STRAINS[@]}")
STRAINS=("${SMOKE_STRAINS[@]}")
ANCHOR_STRAINS=("${SMOKE_STRAINS[@]:0:1}")  # only first
export STRAINS ANCHOR_STRAINS

# Subset assemblies + annotations to the smoke chromosome
for S in "${STRAINS[@]}"; do
  if [[ ! -s $WORK/inputs/assemblies_smoke/${S}.fa ]]; then
    mkdir -p $WORK/inputs/assemblies_smoke
    samtools faidx $WORK/inputs/assemblies/${S}.fa $SMOKE_CHROM > $WORK/inputs/assemblies_smoke/${S}.fa
    awk -F'\t' -v c=$SMOKE_CHROM '$1==c || /^#/' $WORK/inputs/annotations/${S}.fixed.gff3 > $WORK/inputs/annotations_smoke/${S}.fixed.gff3
  fi
done
# Override paths
WORK_FULL=$WORK
WORK=$WORK/smoke_test
mkdir -p $WORK
ln -sfn $WORK_FULL/inputs/assemblies_smoke $WORK/inputs/assemblies
ln -sfn $WORK_FULL/inputs/annotations_smoke $WORK/inputs/annotations
# Run phases A-H, J
bash pipeline/01_inventory.sh
bash pipeline/02_mask.sh
bash pipeline/03_align_chain.sh
bash pipeline/04_annotate_project.sh
bash pipeline/05_pggb.sh
bash pipeline/06_consensus.sh
bash pipeline/07_msa.sh
bash pipeline/08_trees.sh
bash pipeline/09_hyphy.sh
bash pipeline/11_project_vcf.sh
bash pipeline/lib/verify_essentials.sh --smoke
```

- Wall: ~30 min on the v3 hardware.
- Exit 0 only if all `*`-essential outputs (smoke-subset versions) exist and pass validation.

---

## 17. Failure modes seen in v3 (lessons)

### F.1 — Liftoff fails silently when GFF chroms don't match FASTA
- Symptom: `n matched to gff = 0` in liftoff log; output GFF has zero `gene` features.
- Cause: PlasmoDB GFFs use `PvP01_NN_v1` chrnames; assembly FASTAs use `LT635NNN.M`.
- Fix: run `pipeline/lib/fix_gff_chroms.sh ${S}` before Phase C.1. It rewrites the GFF chrom column based on a 2-column rename table or builds one from FASTA headers.

### F.2 — `pal2nal` drops strains with internal stop codons
- Symptom: codon MSA has fewer rows than protein MSA.
- Cause: pseudogene CDS has TAA/TAG/TGA mid-sequence.
- Fix: strip internal stops in `build_msa.py` before `pal2nal` (replace with NNN).

### F.3 — CrossMap output unsorted, breaks `bcftools index`
- Symptom: `bcftools index: failed to create index`.
- Cause: CrossMap interleaves chromosomes during chain-spanning liftover.
- Fix: pipe through `bcftools sort -T $SCRATCH` before `index`.

### F.4 — Path B header-injection bug
- Symptom: per-chr cohort VCFs are 16 KB stubs (header only, no records).
- Cause: Python rewriter dropped `##contig=` lines that `bcftools sort` needs for chromosome ordering.
- Fix: inject `##contig=<ID=...,length=...,assembly=...>` lines from the target `.fa.fai` into the rewritten VCF header.

### F.5 — `rm -f $vcf` ran unconditionally after sort failure
- Symptom: lost the unsorted CrossMap output without ever producing a sorted version.
- Cause: shell script deleted the intermediate VCF inside the sort branch's `else` arm.
- Fix: only delete the intermediate after the sorted `.vcf.gz` is verified non-empty AND indexed.

### F.6 — IQ-TREE3 hangs on < 4 unique sequences with `-B 1000`
- Symptom: iqtree process at 100% CPU, never produces output.
- Cause: bootstrap support needs ≥ 4 unique seqs; with fewer, IQ-TREE loops trying to resample.
- Fix: count unique seqs before invoking; drop `-B 1000` for genes with < 4.

### F.7 — GENESPACE `parse_annotations` returned 0 matches
- Symptom: GENESPACE pan-gene matrix has 0 rows.
- Cause: GFF attribute column lacks the `Name=` field GENESPACE expects.
- Fix: substitute OrthoFinder3 (the v3 chosen path). Filed `jtlovell/GENESPACE#206`.

### F.8 — Pvs230 MNS event masquerades as positive selection
- Symptom: single-hit BUSTED p=0.056 at codon 720; alleles D720A + D720N at 13.3% in MalariaGEN.
- Cause: the two "parallel diversifying" alleles are actually a single multi-nucleotide substitution (D720T → reaching A or N).
- Fix: re-run BUSTED with `--multiple-hits Double`. p shifts to 0.500.
- Lesson: BUSTED-MH (multiple-hits) should be the default for top hits before publication.

---

## 18. Output verification (`pipeline/lib/verify_essentials.sh`)

```bash
#!/usr/bin/env bash
set -uo pipefail
source "${WORK:-$(pwd)}/pipeline/species.conf"
cd $WORK

MISSING=0
check() {
  local label="$1"; local path="$2"; local validator="$3"
  if [[ ! -e $path ]]; then
    echo "MISSING: $label ($path)"
    MISSING=$((MISSING+1))
    return
  fi
  if ! eval "$validator"; then
    echo "INVALID: $label ($path)"
    MISSING=$((MISSING+1))
    return
  fi
  echo "OK:      $label"
}

# Loop through all *-marked rows in OUTLINE.md, parse paths, and validate
# This is a sketch — full implementation reads each table row.

check "Mash matrix"          "work/00_inventory/mash/dist.tsv" '[[ $(wc -l < work/00_inventory/mash/dist.tsv) -gt 1 ]]'
for S in "${STRAINS[@]}"; do
  check "softmasked $S"      "genomes/softmasked/${S}.fa"      "samtools quickcheck -v genomes/softmasked/${S}.fa"
done
check "PGGB graph .og"        "inputs/pggb/pv.og"               '[[ -s inputs/pggb/pv.og ]]'
check "PGGB graph .gfa"       "inputs/pggb/pv.gfa"              '[[ -s inputs/pggb/pv.gfa ]]'
# ... (continue for all 27 *-essentials) ...
check "Ortho table"           "work/03_consensus/ortholog_table.tsv" '[[ $(wc -l < work/03_consensus/ortholog_table.tsv) -gt 1000 ]]'
# Multiz hinges
for H in "${STRAINS[@]}"; do
  check "Multiz $H"           "work/07_multiz/${H}/${H}.multiz.maf" "awk '/^a / {n++} END {exit !(n>=100)}' work/07_multiz/${H}/${H}.multiz.maf"
done
# A2 cohort VCFs
for TGT in "${STRAINS[@]}"; do
  [[ $TGT == $REF_STRAIN ]] && continue
  check "A2 cohort $TGT"      "projection/A2_lastz/Pv4_cohort_on_${TGT}.vcf.gz" "bcftools view -h projection/A2_lastz/Pv4_cohort_on_${TGT}.vcf.gz | grep -q '^##contig='"
done

if [[ $MISSING -eq 0 ]]; then
  echo "ALL 27 ESSENTIALS VERIFIED"
  exit 0
else
  echo "$MISSING essential outputs missing or invalid"
  exit 1
fi
```

---

## 19. Adapt to new species — Pf3D7 example

| `species.conf` field | v3 (P. vivax) | Pf3D7 (P. falciparum) |
|---|---|---|
| `SPECIES` | `Pv4` | `Pf3D7` |
| `STRAINS` | `(PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087)` | `(3D7 HB3 7G8 IT Dd2 CD01 GA01 GB4)` |
| `REF_STRAIN` | `PvP01` | `3D7` |
| `ANCHOR_STRAINS` | `(PvP01 Sal-I PvW1 PAM PvSY56)` | `(3D7 HB3 IT Dd2)` |
| `COHORT_VCF_DIR` | `/media/anton/scratch/malariagen_pv4` | `/media/anton/scratch/malariagen_pf6` |
| `COHORT_CHROM_GLOB` | `'Pv4_PvP01_*_v1.vcf.gz'` | `'Pf6_3D7_*.vcf.gz'` |
| `CHROM_RENAME` | `inputs/annotations/PvP01_plasmodb_to_genbank.tsv` | `inputs/annotations/3D7_plasmodb_to_genbank.tsv` |
| `VAR_ANTIGEN_RE` | `'PIR\|PHIST\|Pv-fam\|MSP\|DBP\|EBA\|RBP\|AMA\|RAP\|SERA\|TRAg\|STP1\|RESA'` | `'PfEMP1\|var\|rifin\|stevor\|PHIST\|surfin\|SURFIN\|EMP3'` |
| BUSCO lineage (in script 01) | `plasmodium_odb10` | `plasmodium_odb10` |
| KegAlign HSP threshold | 5000 | 4000 (Pf is more AT-rich, more spurious matches; lower to filter) |

Everything else (scripts, conda envs, validation logic) is unchanged.

---

## Open questions / unresolved

1. KegAlign GPU vs CPU lastZ — should the pipeline auto-detect and fall back? Currently the implementer must read `GPU=` from species.conf and branch manually.
2. TOGA2 install — pip / git clone / docker? Each has trade-offs (TOGA2 has many Python deps; docker is ~3 GB image).
3. Smoke test currently skips Phase I (multiz). Should it also include a 2-strain smoke version of multiz, or is that just busywork?
4. `phase_e_consensus.py` is ~500 lines in v3; should the recipe inline it or assume it's pre-shipped with the pipeline?
5. PanSN renaming — should we ship the 16-line `pansn_rename.py` inline in this doc or as a `pipeline/scripts/` file? Inline is faithful to PLAN.md style but bloats the recipe.
