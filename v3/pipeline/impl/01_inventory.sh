#!/usr/bin/env bash
set -euo pipefail
# Phase A — Inventory: mash N×N distance matrix + BUSCO per-strain completeness.

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

retry_once() { if ! eval "$*"; then log "retry: $*"; eval "$*"; fi; }

mkdir -p "$WORK/work/00_inventory/mash"
mkdir -p "$WORK/work/00_inventory/busco"
mkdir -p "$WORK/logs"

# ---------------------------------------------------------------------------
# 1.1 Mash N×N distance matrix
# ---------------------------------------------------------------------------
MASH_OUT="$WORK/work/00_inventory/mash/dist.tsv"
MASH_SKETCH="$WORK/work/00_inventory/mash/sketch.msh"

if [[ ! -s "$MASH_OUT" ]]; then
  log "Mash sketch (k=21, s=10000, p=$N_CORES)..."
  FA_ARGS=()
  for S in "${STRAINS[@]}"; do
    FA_ARGS+=("$WORK/inputs/assemblies/${S}.fa")
  done
  retry_once cmd mash sketch -p "$N_CORES" -k 21 -s 10000 \
    -o "$WORK/work/00_inventory/mash/sketch" \
    "${FA_ARGS[@]}"

  log "Mash dist..."
  cmd mash dist -t "$MASH_SKETCH" "$MASH_SKETCH" > "$MASH_OUT"
fi

# Validation
N="${#STRAINS[@]}"
if ! awk -v n="$N" 'NR==1 && NF==n+1 {ok=1} END {exit !ok}' "$MASH_OUT"; then
  log "ERROR: mash dist.tsv has wrong number of columns (expected $((N+1)))" >&2
  rm -f "$MASH_OUT"; exit 1
fi
log "Mash matrix OK — $(wc -l < "$MASH_OUT") rows"

# ---------------------------------------------------------------------------
# 1.2 BUSCO per-strain completeness (on proteomes)
# ---------------------------------------------------------------------------
mkdir -p "$WORK/inputs/proteomes"

for S in "${STRAINS[@]}"; do
  BUSCO_OUT="$WORK/work/00_inventory/busco/${S}_proteins"
  [[ -d "$BUSCO_OUT" ]] && continue

  PROTEOME="$WORK/inputs/proteomes/${S}.proteins.fa"
  if [[ ! -s "$PROTEOME" ]]; then
    log "WARNING: proteome not found for $S ($PROTEOME) — extract via gffread first" >&2
    log "Run: cmd gffread $WORK/inputs/annotations/${S}.fixed.gff3 -g $WORK/inputs/assemblies/${S}.fa -y $PROTEOME" >&2
    continue
  fi

  log "BUSCO $S (plasmodium_odb10, prot mode)..."
  cmd busco \
    -i "$PROTEOME" \
    -l plasmodium_odb10 \
    -m prot \
    -o "${S}_proteins" \
    --out_path "$WORK/work/00_inventory/busco/" \
    --download_path "$WORK/work/00_inventory/busco/busco_downloads/" \
    --cpu "$N_CORES"

  # Validate
  SUMMARY_GLOB="$BUSCO_OUT/short_summary.specific.plasmodium_odb10.${S}_proteins.txt"
  if ! grep -q '^	C:[0-9]' $SUMMARY_GLOB 2>/dev/null; then
    log "WARNING: BUSCO summary looks empty for $S" >&2
  else
    log "BUSCO $S OK"
  fi
done

log "Phase A complete"
