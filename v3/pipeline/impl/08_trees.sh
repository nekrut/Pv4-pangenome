#!/usr/bin/env bash
set -euo pipefail
# Phase G — ML trees with IQ-TREE3 (model selection + bootstrap).
# Output: work/06_msa/{set}_trees/{gene}/{gene}.treefile

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/logs"

run_trees() {
  local SET="$1"
  local MSA_DIR="$WORK/work/06_msa/${SET}"
  local TREE_DIR="$WORK/work/06_msa/${SET}_trees"
  mkdir -p "$TREE_DIR"

  N_SOURCE=$(ls "$MSA_DIR/"*.codon.aln.fa 2>/dev/null | wc -l)
  if [[ $N_SOURCE -eq 0 ]]; then
    log "WARNING: no codon alignments in $MSA_DIR — skipping tree building" >&2
    return
  fi
  log "IQ-TREE3 for $SET ($N_SOURCE alignments, $(( N_CORES / 2 )) parallel)..."

  ls "$MSA_DIR/"*.codon.aln.fa | \
    xargs -P "$(( N_CORES / 2 ))" -I {} bash -c '
      set -euo pipefail
      f={}
      gene=$(basename "$f" .codon.aln.fa)
      out_dir="'"$TREE_DIR"'/$gene"
      [[ -s "$out_dir/${gene}.treefile" ]] && exit 0
      mkdir -p "$out_dir"
      # Count unique sequences (IQ-TREE hangs on < 4 unique with -B 1000)
      n_unique=$(awk "/^>/{next} {print}" "$f" | sort -u | wc -l)
      if [[ $n_unique -ge 4 ]]; then
        '"$WORK"'/pipeline/lib/run_in_'"${RUN_IN_CONTAINER:-container}"'.sh iqtree3 \
          -s "$f" -m MFP -B 1000 -T 2 -pre "$out_dir/$gene" \
          > "$out_dir/${gene}.log" 2>&1 || true
      else
        '"$WORK"'/pipeline/lib/run_in_'"${RUN_IN_CONTAINER:-container}"'.sh iqtree3 \
          -s "$f" -m MFP -T 2 -pre "$out_dir/$gene" \
          > "$out_dir/${gene}.log" 2>&1 || true
      fi
    '

  N_DONE=$(find "$TREE_DIR" -name "*.treefile" 2>/dev/null | wc -l)
  log "IQ-TREE3 $SET done — $N_DONE treefiles"
}

run_trees "core_v3"
run_trees "core_relaxed"

log "Phase G complete"
