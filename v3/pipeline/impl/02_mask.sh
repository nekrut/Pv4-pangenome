#!/usr/bin/env bash
set -euo pipefail
# Phase B — Soft-masking: longdust + sdust → union BED → bedtools maskfasta.
# Output per strain: genomes/softmasked/{S}.fa + .fa.fai + .sizes

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/genomes/mask"
mkdir -p "$WORK/genomes/softmasked"
mkdir -p "$WORK/logs"

for S in "${STRAINS[@]}"; do
  FA_IN="$WORK/inputs/assemblies/${S}.fa"
  SIZES="$WORK/genomes/softmasked/${S}.sizes"

  # Idempotency guard on final output
  if [[ -s "$SIZES" ]]; then
    log "$S already softmasked — skipping"
    continue
  fi

  log "Masking $S..."

  # 2.1 Identify low-complexity regions
  LONGDUST_BED="$WORK/genomes/mask/${S}.longdust.bed"
  SDUST_BED="$WORK/genomes/mask/${S}.sdust.bed"
  UNION_BED="$WORK/genomes/mask/${S}.union.bed"

  if [[ ! -s "$LONGDUST_BED" ]]; then
    cmd longdust "$FA_IN" | \
      awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' > "$LONGDUST_BED"
  fi

  if [[ ! -s "$SDUST_BED" ]]; then
    cmd sdust "$FA_IN" > "$SDUST_BED"
  fi

  if [[ ! -s "$UNION_BED" ]]; then
    cat "$LONGDUST_BED" "$SDUST_BED" | \
      cmd bedtools sort -i - | \
      cmd bedtools merge -i - > "$UNION_BED"
  fi

  # Validate union BED
  if ! awk '$3 > $2' "$UNION_BED" | wc -l | (read n; [[ $n -gt 0 ]]); then
    log "WARNING: union BED for $S has no intervals — check input FASTA" >&2
  fi

  # 2.2 Apply soft mask
  FA_OUT="$WORK/genomes/softmasked/${S}.fa"
  if [[ ! -s "$FA_OUT" ]]; then
    cmd bedtools maskfasta -soft \
      -fi "$FA_IN" \
      -bed "$UNION_BED" \
      -fo "$FA_OUT"
  fi

  # Index + sizes
  cmd samtools faidx "$FA_OUT"
  cut -f1,2 "${FA_OUT}.fai" > "$SIZES"

  # Structural validation: faidx index exists (faidx fails on malformed FASTA)
  if [[ ! -s "${FA_OUT}.fai" ]]; then
    log "ERROR: faidx index missing for $S" >&2
    rm -f "$FA_OUT" "${FA_OUT}.fai" "$SIZES"
    exit 1
  fi

  log "$S softmasked OK — $(wc -l < "$SIZES") sequences"
done

log "Phase B complete"
