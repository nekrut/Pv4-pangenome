#!/usr/bin/env bash
set -euo pipefail
# Smoke test: run all pipeline phases on SMOKE_STRAINS (3 strains) and
# first 300 kb of SMOKE_CHROM, including TOGA2.  Wall budget ~45 min.
#
# Skips Phase I (multiz needs all strains; impractical in smoke with 3).
# Validates a subset of essentials at the end.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${WORK:-$SCRIPT_DIR/..}/pipeline/species.conf"
cd "$WORK"

log() { echo "$(date +%H:%M:%S) [smoke_test] $1" >&2; }

if [[ "${SMOKE_CHROM}" == TODO* ]]; then
  log "ERROR: SMOKE_CHROM not set in species.conf — update after downloading assemblies" >&2
  log "Set SMOKE_CHROM to the GenBank accession of chromosome 1 in the $REF_STRAIN assembly" >&2
  exit 1
fi

log "=== Smoke test START ==="
log "  SMOKE_STRAINS=${SMOKE_STRAINS[*]}"
log "  SMOKE_CHROM=$SMOKE_CHROM"
log "  Region: first 300 kb"

# ---------------------------------------------------------------------------
# Build smoke inputs (subset assembly + annotation)
# ---------------------------------------------------------------------------
SMOKE_WORK="$WORK/smoke_test"
mkdir -p "$SMOKE_WORK"
mkdir -p "$SMOKE_WORK/inputs/assemblies_smoke"
mkdir -p "$SMOKE_WORK/inputs/annotations_smoke"

for S in "${SMOKE_STRAINS[@]}"; do
  SMOKE_FA="$WORK/inputs/assemblies_smoke/${S}.fa"
  SMOKE_GFF="$WORK/inputs/annotations_smoke/${S}.fixed.gff3"

  if [[ ! -s "$SMOKE_FA" ]]; then
    log "Extracting first 300 kb of $SMOKE_CHROM for $S..."
    # Extract the smoke chromosome, then truncate to 300 kb
    samtools faidx "$WORK/inputs/assemblies/${S}.fa" \
      "${SMOKE_CHROM}:1-300000" > "$SMOKE_FA"
  fi

  if [[ ! -s "$SMOKE_GFF" ]]; then
    log "Subsetting annotation for $S smoke region..."
    awk -F'\t' -v c="$SMOKE_CHROM" \
      '$1==c || /^#/' \
      "$WORK/inputs/annotations/${S}.fixed.gff3" \
      > "$SMOKE_GFF"
  fi
done

# ---------------------------------------------------------------------------
# Override WORK + STRAINS to smoke subset
# ---------------------------------------------------------------------------
WORK_FULL="$WORK"
export WORK="$SMOKE_WORK"
export STRAINS=("${SMOKE_STRAINS[@]}")
export ANCHOR_STRAINS=("${SMOKE_STRAINS[0]}")   # only first strain as anchor
export REF_STRAIN="${SMOKE_STRAINS[0]}"
export N_CORES="${N_CORES}"
export MIN_INTACT_STRICT=2
export MIN_INTACT_RELAXED=2

# Symlink inputs into SMOKE_WORK
ln -sfn "$WORK_FULL/inputs/assemblies_smoke" "$SMOKE_WORK/inputs/assemblies"
ln -sfn "$WORK_FULL/inputs/annotations_smoke" "$SMOKE_WORK/inputs/annotations"

# Write a minimal species.conf for the smoke run inside the smoke_work dir
mkdir -p "$SMOKE_WORK/pipeline"
cat > "$SMOKE_WORK/pipeline/species.conf" <<CONF
SPECIES=${SPECIES}_smoke
WORK=$SMOKE_WORK
STRAINS=(${SMOKE_STRAINS[@]})
REF_STRAIN=${SMOKE_STRAINS[0]}
ANCHOR_STRAINS=(${SMOKE_STRAINS[0]})
COHORT_VCF_DIR=TODO_NO_COHORT_VCF_AVAILABLE
COHORT_CHROM_GLOB='${COHORT_CHROM_GLOB}'
CHROM_RENAME=$WORK_FULL/inputs/annotations/${REF_STRAIN}_plasmodb_to_genbank.tsv
N_CORES=${N_CORES}
GPU=${GPU:-0}
USE_GPU=${USE_GPU:-0}
DOCKER_GPU="--gpus device=${GPU:-0}"
MIN_INTACT_STRICT=2
MIN_INTACT_RELAXED=2
VAR_ANTIGEN_RE='${VAR_ANTIGEN_RE}'
SMOKE_CHROM=${SMOKE_CHROM}
SMOKE_STRAINS=(${SMOKE_STRAINS[@]})
CONF

# Symlink pipeline directory so scripts can find their helpers
ln -sfn "$WORK_FULL/pipeline/lib"     "$SMOKE_WORK/pipeline/lib"
ln -sfn "$WORK_FULL/pipeline/scripts" "$SMOKE_WORK/pipeline/scripts"
ln -sfn "$WORK_FULL/pipeline/setup"   "$SMOKE_WORK/pipeline/setup"

log "Smoke config written; running phases A → H + J (skipping I/multiz)..."

# Run each phase targeting the smoke WORK dir
export RUN_IN_CONTAINER="${RUN_IN_CONTAINER:-container}"

for PHASE in 01 02 03 04 05 06 07 08 09 11; do
  SCRIPT="$WORK_FULL/pipeline/${PHASE}_"*.sh
  SCRIPT_PATH=$(echo $SCRIPT)   # expand glob to path
  [[ -f "$SCRIPT_PATH" ]] || { log "WARNING: script not found: $SCRIPT" >&2; continue; }
  PHASE_NAME=$(basename "$SCRIPT_PATH" .sh)
  log "--- Phase $PHASE_NAME ---"
  # Run with smoke WORK override
  WORK="$SMOKE_WORK" bash "$SCRIPT_PATH" 2>&1 | tee "$SMOKE_WORK/logs/smoke_${PHASE}.log" >&2
  log "--- Phase $PHASE_NAME DONE ---"
done

log "Smoke phases complete — running verify_essentials.sh --smoke"
WORK="$SMOKE_WORK" bash "$WORK_FULL/pipeline/lib/verify_essentials.sh" --smoke

log "=== Smoke test COMPLETE ==="
