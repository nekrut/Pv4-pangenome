#!/usr/bin/env bash
set -euo pipefail
# Phase J — Cohort VCF projection via CrossMap + bcftools.
#
# 11.1 Rename source cohort VCFs to GenBank chromosome names.
# 11.2 CrossMap each per-chr VCF onto each non-reference target strain.
# 11.3 bcftools concat per-chr → per-target cohort VCF.
#
# Skips gracefully if COHORT_VCF_DIR is a TODO placeholder.

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/logs"

# Check if cohort VCF is available
if [[ "${COHORT_VCF_DIR}" == TODO* ]]; then
  log "COHORT_VCF_DIR is a TODO placeholder — Phase J skipped"
  log "When a P. knowlesi population VCF becomes available, set COHORT_VCF_DIR in species.conf"
  exit 0
fi

if [[ ! -d "$COHORT_VCF_DIR" ]]; then
  log "ERROR: COHORT_VCF_DIR does not exist: $COHORT_VCF_DIR" >&2
  exit 1
fi

if [[ "${CHROM_RENAME}" == TODO* ]] || [[ ! -s "$CHROM_RENAME" ]]; then
  log "ERROR: CHROM_RENAME file missing — needed for Phase J" >&2
  exit 1
fi

mkdir -p "$WORK/projection/A1_wfmash/mg_renamed"
TMP_SORT="${SCRATCH:-/tmp}/Pk_vcfsort_$$"
mkdir -p "$TMP_SORT"
trap 'rm -rf "$TMP_SORT"' EXIT

# ---------------------------------------------------------------------------
# 11.1 Rename cohort VCF chroms to GenBank names
# ---------------------------------------------------------------------------
log "Phase J.1 — renaming cohort VCF chromosomes..."

for SRC_VCF in "$COHORT_VCF_DIR/"$COHORT_CHROM_GLOB; do
  [[ -f "$SRC_VCF" ]] || continue
  CHR_TAG=$(basename "$SRC_VCF" .vcf.gz | sed -E "s/^${SPECIES}_//")
  OUT="$WORK/projection/A1_wfmash/mg_renamed/${SPECIES}_${CHR_TAG}.vcf.gz"
  [[ -s "${OUT}.csi" ]] && continue
  cmd bcftools annotate --rename-chrs "$CHROM_RENAME" "$SRC_VCF" -Oz -o "$OUT"
  cmd bcftools index "$OUT"
  log "  Renamed $CHR_TAG → $OUT"
done

# Validate first output
FIRST_RENAMED=$(ls "$WORK/projection/A1_wfmash/mg_renamed/"*.vcf.gz 2>/dev/null | head -1)
if [[ -z "$FIRST_RENAMED" ]]; then
  log "ERROR: No renamed VCFs produced — check COHORT_CHROM_GLOB pattern" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# 11.2 CrossMap projection + 11.3 concat per target
# ---------------------------------------------------------------------------
log "Phase J.2 — CrossMap projection to each non-reference target..."

for TGT in "${STRAINS[@]}"; do
  [[ "$TGT" == "$REF_STRAIN" ]] && continue
  OUT_COHORT="$WORK/projection/A2_lastz/${SPECIES}_cohort_on_${TGT}.vcf.gz"
  [[ -s "${OUT_COHORT}.csi" ]] && continue

  CHAIN="$WORK/work/01_chains/${REF_STRAIN}.${TGT}.cleaned.chain"
  TGT_FA="$WORK/genomes/softmasked/${TGT}.fa"
  PER_CHR_DIR="$WORK/projection/A2_lastz/per_chr/${TGT}"
  mkdir -p "$PER_CHR_DIR"

  if [[ ! -s "$CHAIN" ]]; then
    log "WARNING: chain $CHAIN missing — skipping CrossMap for $TGT" >&2
    continue
  fi

  log "  CrossMap → $TGT..."
  for SRC_VCF in "$WORK/projection/A1_wfmash/mg_renamed/${SPECIES}_"*.vcf.gz; do
    [[ -f "$SRC_VCF" ]] || continue
    CHR_TAG=$(basename "$SRC_VCF" .vcf.gz | sed "s/^${SPECIES}_//")
    OUT_CHR="$PER_CHR_DIR/${SPECIES}_${CHR_TAG}_on_${TGT}.vcf.gz"
    [[ -s "${OUT_CHR}.csi" ]] && continue

    # CrossMap to a temp file (cleaner than a cross-container stdout pipe),
    # then bcftools sort (CrossMap output is unsorted — F.3 lesson).
    RAW_VCF="$PER_CHR_DIR/${SPECIES}_${CHR_TAG}_on_${TGT}.raw.vcf"
    cmd CrossMap vcf "$CHAIN" "$SRC_VCF" "$TGT_FA" "$RAW_VCF"
    cmd bcftools sort -T "$TMP_SORT" -Oz -o "$OUT_CHR" "$RAW_VCF"
    rm -f "$RAW_VCF" "${RAW_VCF}.unmap"
    cmd bcftools index "$OUT_CHR"
  done

  # Concat per-chr → full cohort VCF
  PER_CHR_VCFS=("$PER_CHR_DIR/${SPECIES}_"*"_on_${TGT}.vcf.gz")
  if [[ ${#PER_CHR_VCFS[@]} -eq 0 ]]; then
    log "WARNING: no per-chr VCFs for $TGT" >&2
    continue
  fi
  cmd bcftools concat -a "${PER_CHR_VCFS[@]}" -Oz -o "$OUT_COHORT"
  cmd bcftools index "$OUT_COHORT"

  # Validate
  N_VARS=$(cmd bcftools view -H "$OUT_COHORT" | wc -l)
  if [[ $N_VARS -lt 100000 ]]; then
    log "WARNING: only $N_VARS variants in ${TGT} cohort VCF (expected >= 100k)" >&2
  else
    log "  CrossMap $TGT OK — $N_VARS variants"
  fi
done

log "Phase J complete"
