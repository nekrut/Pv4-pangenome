#!/usr/bin/env bash
# Build bigMaf + mafIndex for each multiz hinge, 3-concurrent.
# Bypass mafToBigMaf (rejects overlapping blocks); use a custom Python BED3+1 emitter.
set -uo pipefail

TOOLS=/media/anton/data/sandbox/Pv4/v3/tools
HUB=/media/anton/data/sandbox/Pv4/v3/ucsc_hub
STAGING=/media/anton/scratch/Pv4_dropbox_staging
TMP=/media/anton/scratch/maf_tmp
mkdir -p "$TMP"

process_one() {
  local HINGE=$1 ACC=$2
  local OUT_DIR="$HUB/$ACC"
  local OUT_BB="$OUT_DIR/${HINGE}.multiz.maf.bb"
  local OUT_BAI="$OUT_DIR/${HINGE}.multiz.maf.bb.bai"
  local LOG="$OUT_DIR/${HINGE}.maf.log"
  local MAF_GZ="$STAGING/multiz/${HINGE}/${HINGE}.multiz.maf.gz"
  local RAW="$TMP/${HINGE}.raw.maf"
  local FIXED="$TMP/${HINGE}.fixed.maf"
  local BED="$TMP/${HINGE}.bed"
  local SORTED="$TMP/${HINGE}.sorted.bed"
  local SIZES="$TMP/${ACC}.sizes"

  mkdir -p "$OUT_DIR"
  if [[ -s "$OUT_BB" && -s "$OUT_BAI" && $(stat -c '%s' "$OUT_BB") -gt 5000000 ]]; then
    echo "[$HINGE] SKIP (already built)"
    return 0
  fi

  echo "[$HINGE] start"
  rm -f "$OUT_BB" "$OUT_BAI" "$LOG"

  gunzip -c "$MAF_GZ" > "$RAW"
  python3 "$HUB/process_maf.py" "$ACC" "$RAW" "$FIXED" 2>&1 | tee -a "$LOG"
  awk '{print $1"\t"$2}' "$STAGING/softmasked/${ACC}.fa.fai" > "$SIZES"

  echo "[$HINGE] custom BED3+1 emit..."
  python3 "$HUB/maf_to_bigmaf_bed.py" "$ACC" "$FIXED" "$BED" 2>&1 | tee -a "$LOG"

  echo "[$HINGE] sort + bedToBigBed..."
  sort -k1,1 -k2,2n "$BED" > "$SORTED"
  if ! "$TOOLS/bedToBigBed" -type=bed3+1 -as="$HUB/bigMaf.as" -tab "$SORTED" "$SIZES" "$OUT_BB" 2>>"$LOG"; then
    echo "[$HINGE] bedToBigBed FAILED, see $LOG"
    return 1
  fi

  if [[ ! -s "$OUT_BB" ]]; then
    echo "[$HINGE] FAILED at bedToBigBed (no output)"
    return 1
  fi

  echo "[$HINGE] mafIndex..."
  "$TOOLS/mafIndex" "$FIXED" "$OUT_BAI" 2>>"$LOG" || echo "[$HINGE] mafIndex failed (non-fatal)"

  echo "[$HINGE] DONE: $(du -h "$OUT_BB" | cut -f1) bb"
  rm -f "$RAW" "$FIXED" "$BED" "$SORTED" "$SIZES"
}
export -f process_one
export TOOLS HUB STAGING TMP

# Force rebuild PvP01 too
rm -f "$HUB/GCA_900093555.2/PvP01.multiz.maf.bb"

echo "PvP01 GCA_900093555.2
Sal-I GCA_000002415.2
PvW1 GCA_914969965.1
PAM GCA_949152365.1
PvSY56 GCA_003402215.1
PvT01 GCA_900093545.1
PvC01 GCA_900093535.1
MHC087 GCA_040114635.1" | xargs -P 3 -L 1 bash -c 'process_one "$@"' _

echo
echo "=== FINAL bigMaf files ==="
ls -la "$HUB"/*/*.multiz.maf.bb 2>/dev/null
