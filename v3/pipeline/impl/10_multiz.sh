#!/usr/bin/env bash
set -euo pipefail
# Phase I — Multiz multi-way MAFs.
# For each strain H (hinge), build pairwise H-vs-Q MAFs from AXT, then progressively
# fold into N-way MAF via multiz.
# Output: work/07_multiz/{H}/{H}.multiz.maf

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/work/07_multiz"
mkdir -p "$WORK/logs"

for H in "${STRAINS[@]}"; do
  OUT_DIR="$WORK/work/07_multiz/${H}"
  OUT_FINAL="$OUT_DIR/${H}.multiz.maf"
  [[ -s "$OUT_FINAL" ]] && continue

  log "Multiz hinge=$H..."
  mkdir -p "$OUT_DIR"

  # Step 1: build all pairwise MAFs (H vs every Q)
  for Q in "${STRAINS[@]}"; do
    [[ "$Q" == "$H" ]] && continue
    PAIR_MAF="$OUT_DIR/${H}_vs_${Q}.maf"
    [[ -s "$PAIR_MAF" ]] && continue

    # Pick correct AXT (H must be the target/reference in the MAF)
    if [[ -s "$WORK/projection/A2_kegalign/axt/${H}__vs__${Q}.axt" ]]; then
      INPUT_AXT="$WORK/projection/A2_kegalign/axt/${H}__vs__${Q}.axt"
      cmd axtToMaf \
        -tPrefix="${H}." \
        -qPrefix="${Q}." \
        "$INPUT_AXT" \
        "$WORK/genomes/softmasked/${H}.sizes" \
        "$WORK/genomes/softmasked/${Q}.sizes" \
        "$PAIR_MAF"
    elif [[ -s "$WORK/projection/A2_kegalign/axt/${Q}__vs__${H}.axt" ]]; then
      INPUT_AXT="$WORK/projection/A2_kegalign/axt/${Q}__vs__${H}.axt"
      # The AXT has Q as target; swap naming so H is the reference in the MAF
      cmd axtToMaf \
        -tPrefix="${Q}." \
        -qPrefix="${H}." \
        "$INPUT_AXT" \
        "$WORK/genomes/softmasked/${Q}.sizes" \
        "$WORK/genomes/softmasked/${H}.sizes" \
        "$PAIR_MAF"
    else
      log "WARNING: no AXT for $H vs $Q — skipping this pair" >&2
    fi
  done

  # Step 2: order other strains by mash distance from H (closest first)
  # Use a simple approach: sort STRAINS by lexicographic order minus H as fallback.
  # Ideally parse dist.tsv — done inline here.
  MASH_TSV="$WORK/work/00_inventory/mash/dist.tsv"
  ORDERED=()
  if [[ -s "$MASH_TSV" ]]; then
    # Extract column index for H, then sort remaining strains by distance to H
    ORDERED_STR=$(python3 - "$H" "$MASH_TSV" "${STRAINS[@]}" <<'PYEOF'
import sys, csv
hinge = sys.argv[1]
tsv   = sys.argv[2]
all_s = sys.argv[3:]
with open(tsv) as fh:
    rows = list(csv.reader(fh, delimiter='\t'))
header = rows[0]  # row 0 is col names (mash uses full paths as headers)
# Find column for hinge (header entry contains hinge name)
h_col = next((i for i, h in enumerate(header) if hinge in h), -1)
if h_col < 0:
    print(' '.join(s for s in all_s if s != hinge))
    sys.exit(0)
# Find row index for each strain in all_s, read distance to hinge
dists = {}
for row in rows[1:]:
    for s in all_s:
        if s != hinge and s in row[0]:
            try: dists[s] = float(row[h_col])
            except (ValueError, IndexError): pass
ordered = sorted([s for s in all_s if s != hinge], key=lambda s: dists.get(s, 999))
print(' '.join(ordered))
PYEOF
    )
    read -ra ORDERED <<< "$ORDERED_STR"
  else
    for S in "${STRAINS[@]}"; do [[ "$S" != "$H" ]] && ORDERED+=("$S"); done
  fi

  # Step 3: progressive multiz fold
  PREV_MAF="$OUT_DIR/${H}_vs_${ORDERED[0]}.maf"
  if [[ ! -s "$PREV_MAF" ]]; then
    log "WARNING: first pairwise MAF missing for hinge $H vs ${ORDERED[0]}" >&2
    continue
  fi

  TMP_MAF="$PREV_MAF"
  for IDX in "${!ORDERED[@]}"; do
    [[ $IDX -eq 0 ]] && continue
    Q="${ORDERED[$IDX]}"
    NEXT_MAF="$OUT_DIR/${H}_vs_${Q}.maf"
    [[ ! -s "$NEXT_MAF" ]] && continue
    NEW_TMP=$(mktemp "$OUT_DIR/multiz_tmp_XXXXXX.maf")
    cmd multiz "$TMP_MAF" "$NEXT_MAF" 1 > "$NEW_TMP"
    # Remove old tmp (keep pair MAFs)
    [[ "$TMP_MAF" != "$OUT_DIR/${H}_vs_"* ]] && rm -f "$TMP_MAF"
    TMP_MAF="$NEW_TMP"
  done
  mv "$TMP_MAF" "$OUT_FINAL"

  # Validate
  N_BLOCKS=$(awk '/^a /{n++} END{print n}' "$OUT_FINAL")
  if [[ "${N_BLOCKS:-0}" -lt 100 ]]; then
    log "WARNING: multiz MAF for $H has only $N_BLOCKS alignment blocks" >&2
  else
    log "  Multiz $H OK — $N_BLOCKS alignment blocks"
  fi
done

log "Phase I complete"
