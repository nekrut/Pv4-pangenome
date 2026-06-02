#!/usr/bin/env bash
set -euo pipefail
# Phase D — PGGB pangenome graph.
#
# 5.1 PanSN rename each strain's softmasked FASTA + bgzip + index
# 5.2 PGGB build → inputs/pggb/pk.{gfa,og}

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/work/pggb_in"
mkdir -p "$WORK/work/pggb_out"
mkdir -p "$WORK/inputs/pggb"
mkdir -p "$WORK/logs"

# ---------------------------------------------------------------------------
# 5.1 PanSN rename per strain + concat + bgzip
# ---------------------------------------------------------------------------
log "Phase D.1 — PanSN rename..."

for S in "${STRAINS[@]}"; do
  OUT="$WORK/work/pggb_in/${S}_pansn.fa"
  [[ -s "$OUT" ]] && continue
  log "  PanSN rename $S..."
  cmd python3 "$WORK/pipeline/scripts/pansn_rename.py" \
    "$WORK/genomes/softmasked/${S}.fa" "$S" \
    > "$OUT"
done

ALL_GZ="$WORK/work/pggb_in/all_pansn.fa.gz"
if [[ ! -s "${ALL_GZ}.gzi" ]]; then
  log "  Concatenating and bgzipping all PanSN FASTAs..."
  cat "$WORK/work/pggb_in/"*_pansn.fa | \
    cmd bgzip -c > "$ALL_GZ"
  cmd samtools faidx "$ALL_GZ"
fi

# ---------------------------------------------------------------------------
# 5.2 PGGB build
# ---------------------------------------------------------------------------
if [[ ! -L "$WORK/inputs/pggb/pk.og" ]] || [[ ! -s "$WORK/inputs/pggb/pk.og" ]]; then
  log "Phase D.2 — PGGB build (this takes 3-6 hours)..."

  N_STRAINS="${#STRAINS[@]}"
  cmd pggb \
    -i "$ALL_GZ" \
    -n "$N_STRAINS" \
    -s 5000 \
    -p 90 \
    -k 23 \
    -t "$N_CORES" \
    -o "$WORK/work/pggb_out/"

  # Symlink canonical output names
  SMOOTH_GFA=$(ls "$WORK/work/pggb_out/"*.smooth.final.gfa 2>/dev/null | head -1)
  SMOOTH_OG=$(ls "$WORK/work/pggb_out/"*.smooth.final.og 2>/dev/null | head -1)

  if [[ -z "$SMOOTH_GFA" ]] || [[ -z "$SMOOTH_OG" ]]; then
    log "ERROR: PGGB output GFA/OG not found" >&2; exit 1
  fi

  ln -sf "$SMOOTH_GFA" "$WORK/inputs/pggb/pk.gfa"
  ln -sf "$SMOOTH_OG"  "$WORK/inputs/pggb/pk.og"
fi

# Validation: odgi stats must emit a header line
if ! cmd odgi stats -i "$WORK/inputs/pggb/pk.og" 2>/dev/null | grep -q '^#length'; then
  log "ERROR: PGGB .og validation failed (odgi stats)" >&2; exit 1
fi
log "PGGB graph validated OK"

log "Phase D complete"
