#!/usr/bin/env bash
set -euo pipefail
# Phase F — codon + protein MSAs (strict + relaxed sets).
# For each orthogroup with >= MIN_INTACT_STRICT (or MIN_INTACT_RELAXED) 1-to-1 members:
#   build_msa.py → {gene}.codon.aln.fa + {gene}.protein.aln.fa
# Then trimal → core_v3_clean / core_relaxed_clean sibling dirs.

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/logs"

run_msa_set() {
  local SETNAME="$1"    # "strict" or "relaxed"
  local MIN="$2"
  local OUT_DIR="$3"
  local CLEAN_DIR="${OUT_DIR%/}_clean"

  N_DONE=$(ls "$OUT_DIR/"*.codon.aln.fa 2>/dev/null | wc -l || true)
  if [[ -d "$OUT_DIR" && $N_DONE -gt 100 ]]; then
    log "MSA set $SETNAME: $N_DONE files already exist — skipping build step"
  else
    log "MSA set $SETNAME (min_intact=$MIN)..."
    mkdir -p "$OUT_DIR"

    # Gather per-query merged GFF paths (use REF anchor)
    QUERY_GFFS=()
    QUERY_FASTAS=()
    QUERIES=()
    for Q in "${STRAINS[@]}"; do
      [[ "$Q" == "$REF_STRAIN" ]] && continue
      GFF="$WORK/work/02d_merged/${REF_STRAIN}-as-ref/${Q}.annotation.gff3"
      FA="$WORK/genomes/softmasked/${Q}.fa"
      if [[ -s "$GFF" ]] && [[ -s "$FA" ]]; then
        QUERIES+=("$Q")
        QUERY_GFFS+=("$GFF")
        QUERY_FASTAS+=("$FA")
      fi
    done

    # Run sharded across N_CORES/2 shards for parallelism
    N_SHARDS=$(( N_CORES / 2 ))
    [[ $N_SHARDS -lt 1 ]] && N_SHARDS=1

    for SHARD in $(seq 0 $((N_SHARDS - 1))); do
      python3 "$WORK/pipeline/scripts/build_msa.py" \
        --work "$WORK" \
        --ortho "$WORK/work/03_consensus/ortholog_table.tsv" \
        --strains "${STRAINS[*]}" \
        --ref "$REF_STRAIN" \
        --ref-gff "$WORK/inputs/annotations/${REF_STRAIN}.fixed.gff3" \
        --ref-fasta "$WORK/genomes/softmasked/${REF_STRAIN}.fa" \
        --queries "${QUERIES[@]}" \
        --query-gff "${QUERY_GFFS[@]}" \
        --query-fasta "${QUERY_FASTAS[@]}" \
        --merged-base "$WORK/work/02d_merged/${REF_STRAIN}-as-ref" \
        --min-intact "$MIN" \
        --out-dir "$OUT_DIR" \
        --threads-per-gene 2 \
        --shard "$SHARD" \
        --num-shards "$N_SHARDS" &
    done
    wait
    log "MSA set $SETNAME: $(ls "$OUT_DIR/"*.codon.aln.fa 2>/dev/null | wc -l || true) codon MSAs built"
  fi

  # TrimAl cleaning pass
  N_CLEAN=$(ls "$CLEAN_DIR/"*.codon.cleaned.fa 2>/dev/null | wc -l || true)
  N_SOURCE=$(ls "$OUT_DIR/"*.codon.aln.fa 2>/dev/null | wc -l || true)
  if [[ -d "$CLEAN_DIR" && $N_CLEAN -ge $N_SOURCE && $N_SOURCE -gt 0 ]]; then
    log "TrimAl $SETNAME: $N_CLEAN cleaned files already exist — skipping"
  else
    mkdir -p "$CLEAN_DIR"
    log "TrimAl $SETNAME cleaning ($N_SOURCE alignments)..."
    for ALN in "$OUT_DIR/"*.codon.aln.fa; do
      [[ -f "$ALN" ]] || continue
      GENE=$(basename "$ALN" .codon.aln.fa)
      OUT_CLEAN="$CLEAN_DIR/${GENE}.codon.cleaned.fa"
      [[ -s "$OUT_CLEAN" ]] && continue
      cmd trimal -in "$ALN" -out "$OUT_CLEAN" -automated1
    done
    log "TrimAl $SETNAME done"
  fi
}

run_msa_set "strict"  "$MIN_INTACT_STRICT"  "$WORK/work/06_msa/core_v3"
run_msa_set "relaxed" "$MIN_INTACT_RELAXED" "$WORK/work/06_msa/core_relaxed"

log "Phase F complete"
