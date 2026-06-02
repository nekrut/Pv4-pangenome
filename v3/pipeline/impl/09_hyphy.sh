#!/usr/bin/env bash
set -euo pipefail
# Phase H — HyPhy bulk BUSTED per-gene selection scan.
# Output: work/06_msa/core_v3_hyphy/bulk/{gene}/busted.json
#         work/06_msa/core_relaxed_hyphy/{gene}/busted.json

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/logs"

run_hyphy_set() {
  local SET="$1"
  local MSA_DIR="$WORK/work/06_msa/${SET}"
  local TREE_DIR="$WORK/work/06_msa/${SET}_trees"
  local OUT_BASE

  if [[ "$SET" == "core_v3" ]]; then
    OUT_BASE="$WORK/work/06_msa/core_v3_hyphy/bulk"
  else
    OUT_BASE="$WORK/work/06_msa/${SET}_hyphy"
  fi
  mkdir -p "$OUT_BASE"

  N_SOURCE=$(ls "$MSA_DIR/"*.codon.aln.fa 2>/dev/null | wc -l)
  if [[ $N_SOURCE -eq 0 ]]; then
    log "WARNING: no codon alignments in $MSA_DIR — skipping HyPhy" >&2
    return
  fi
  log "HyPhy BUSTED for $SET ($N_SOURCE genes, $N_CORES parallel)..."

  ls "$MSA_DIR/"*.codon.aln.fa | \
    xargs -P "$N_CORES" -I {} bash -c '
      set -euo pipefail
      f={}
      gene=$(basename "$f" .codon.aln.fa)
      out_dir="'"$OUT_BASE"'/$gene"
      [[ -s "$out_dir/busted.json" ]] && exit 0
      tree="'"$TREE_DIR"'/$gene/${gene}.treefile"
      [[ ! -s "$tree" ]] && exit 0   # no tree → skip
      mkdir -p "$out_dir"
      '"$WORK"'/pipeline/lib/run_in_'"${RUN_IN_CONTAINER:-container}"'.sh hyphy busted \
        --alignment "$f" \
        --tree "$tree" \
        --output "$out_dir/busted.json" \
        --srv No \
        --branches All \
        > "$out_dir/busted.log" 2>&1 || true
      # Validate JSON has test results key
      if [[ -s "$out_dir/busted.json" ]]; then
        python3 -c "
import json, sys
try:
    d = json.load(open(sys.argv[1]))
    assert '\''test results'\'' in d, '\''missing test results key'\''
except Exception as e:
    print(f'\''INVALID: {sys.argv[1]}: {e}'\'', file=sys.stderr)
    sys.exit(1)
" "$out_dir/busted.json" 2>/dev/null || rm -f "$out_dir/busted.json"
      fi
    '

  N_DONE=$(find "$OUT_BASE" -name "busted.json" 2>/dev/null | wc -l)
  log "HyPhy BUSTED $SET done — $N_DONE JSONs"
}

run_hyphy_set "core_v3"
run_hyphy_set "core_relaxed"

log "Phase H complete"
