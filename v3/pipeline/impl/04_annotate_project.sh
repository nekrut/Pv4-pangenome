#!/usr/bin/env bash
set -euo pipefail
# Phase C.4 — Annotation projection: Liftoff + triage + TOGA2 rescue + merge.
#
# For each ANCHOR_STRAIN A, project A's annotation onto every non-A strain Q.
# Outputs (per anchor–query pair):
#   work/02a_liftoff/{A}-as-ref/{Q}.liftoff.gff3
#   work/02b_triage/{A}-as-ref/{Q}/needs_cesar2.bed + liftoff_clean.gff3 + triage.tsv + summary.json
#   work/02c_toga/{A}-as-ref/{Q}/annotation.gff3 + classification.tsv
#   work/02d_merged/{A}-as-ref/{Q}.annotation.gff3 + .classification.tsv

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/logs"

# Feature-type filter for Liftoff (file written once)
FEATURE_FILE=$(mktemp /tmp/Pk_liftoff_features_XXXXXX.txt)
trap 'rm -f "$FEATURE_FILE"' EXIT
printf 'protein_coding_gene\nncRNA_gene\npseudogene\n' > "$FEATURE_FILE"

for A in "${ANCHOR_STRAINS[@]}"; do
  # Ensure fixed GFF + BED12 + isoforms exist for this anchor
  if [[ ! -s "$WORK/inputs/annotations/${A}.fixed.gff3" ]]; then
    log "Building fixed GFF for anchor $A..."
    bash "$WORK/pipeline/lib/fix_gff_chroms.sh" "$A"
  fi
  if [[ ! -s "$WORK/inputs/annotations/${A}.bed12" ]]; then
    log "Building BED12 + isoforms for anchor $A..."
    bash "$WORK/pipeline/setup/build_anchor_inputs.sh" "$A"
  fi

  mkdir -p "$WORK/work/02a_liftoff/${A}-as-ref"
  mkdir -p "$WORK/work/02b_triage/${A}-as-ref"
  mkdir -p "$WORK/work/02c_toga/${A}-as-ref"
  mkdir -p "$WORK/work/02d_merged/${A}-as-ref"

  for Q in "${STRAINS[@]}"; do
    [[ "$Q" == "$A" ]] && continue

    # ------------------------------------------------------------------
    # 4.1 Liftoff projection
    # ------------------------------------------------------------------
    LIFTOFF_OUT="$WORK/work/02a_liftoff/${A}-as-ref/${Q}.liftoff.gff3"
    if [[ ! -s "$LIFTOFF_OUT" ]]; then
      log "Liftoff ${A}→${Q}..."
      mkdir -p "$WORK/work/02a_liftoff/${A}-as-ref/${Q}_intermediate"
      cmd liftoff \
        -g "$WORK/inputs/annotations/${A}.fixed.gff3" \
        -f "$FEATURE_FILE" \
        -o "$LIFTOFF_OUT" \
        -dir "$WORK/work/02a_liftoff/${A}-as-ref/${Q}_intermediate" \
        -p "$N_CORES" \
        -copies -sc 0.95 \
        "$WORK/inputs/assemblies/${Q}.fa" \
        "$WORK/inputs/assemblies/${A}.fa"
    fi

    # Validate
    GENE_COUNT=$(awk -F'\t' '$3=="gene" || $3=="protein_coding_gene"' "$LIFTOFF_OUT" 2>/dev/null | wc -l)
    if [[ "$GENE_COUNT" -lt 100 ]]; then
      log "ERROR: Liftoff ${A}→${Q} produced only $GENE_COUNT genes — check GFF chrom names" >&2
      exit 1
    fi
    log "  Liftoff ${A}→${Q}: $GENE_COUNT genes"

    # ------------------------------------------------------------------
    # 4.2 Triage
    # ------------------------------------------------------------------
    TRIAGE_DIR="$WORK/work/02b_triage/${A}-as-ref/${Q}"
    if [[ ! -s "$TRIAGE_DIR/needs_cesar2.bed" ]]; then
      log "Triage ${A}→${Q}..."
      mkdir -p "$TRIAGE_DIR"
      cmd python3 "$WORK/pipeline/scripts/phase_c2_triage.py" \
        --liftoff-gff "$LIFTOFF_OUT" \
        --query-fasta "$WORK/genomes/softmasked/${Q}.fa" \
        --reference-bed "$WORK/inputs/annotations/${A}.bed12" \
        --output-dir "$TRIAGE_DIR" \
        --query-name "$Q" \
        --core-identity-min 0.95 \
        --core-coverage-min 0.90 \
        --family-identity-min 0.85 \
        --subtelomere-bp 100000
    fi

    # ------------------------------------------------------------------
    # 4.3 TOGA2 rescue
    # ------------------------------------------------------------------
    TOGA_DIR="$WORK/work/02c_toga/${A}-as-ref/${Q}"
    if [[ ! -s "$TOGA_DIR/annotation.gff3" ]]; then
      log "TOGA2 ${A}→${Q}..."
      mkdir -p "$TOGA_DIR"
      CHAIN="$WORK/work/01_chains/${A}.${Q}.cleaned.chain"
      if [[ ! -s "$CHAIN" ]]; then
        log "WARNING: chain $CHAIN missing — skipping TOGA2 for ${A}→${Q}" >&2
      elif [[ ! -s "$TRIAGE_DIR/needs_cesar2.bed" ]]; then
        log "  no genes flagged for CESAR rescue — skipping TOGA2 for ${A}→${Q}"
      else
        # TOGA2 runs inside its container; note --kt keeps temporary files for debugging
        cmd toga \
          "$WORK/inputs/annotations/${A}.bed12" \
          "$CHAIN" \
          "$WORK/genomes/softmasked/${A}.fa" \
          "$WORK/genomes/softmasked/${Q}.fa" \
          --pn "$TOGA_DIR" \
          --cb "$WORK/inputs/annotations/${A}.isoforms.tsv" \
          --u12 "" \
          --kt \
          --filter_bed "$TRIAGE_DIR/needs_cesar2.bed" \
          --nc "$N_CORES"
      fi
    fi

    # ------------------------------------------------------------------
    # 4.4 Merge Liftoff-clean + TOGA2 outputs
    # ------------------------------------------------------------------
    MERGED_GFF="$WORK/work/02d_merged/${A}-as-ref/${Q}.annotation.gff3"
    MERGED_TSV="$WORK/work/02d_merged/${A}-as-ref/${Q}.classification.tsv"
    if [[ ! -s "$MERGED_GFF" ]] || [[ ! -s "$MERGED_TSV" ]]; then
      log "Merge ${A}→${Q}..."
      cmd python3 "$WORK/pipeline/scripts/phase_c4_merge.py" \
        --query "$Q" \
        --triage-dir "$TRIAGE_DIR" \
        --toga-dir "$TOGA_DIR" \
        --out-dir "$WORK/work/02d_merged/${A}-as-ref" \
        --ref-bed "$WORK/inputs/annotations/${A}.bed12"
    fi

    # Validate merged GFF
    REF_GENES=$(awk -F'\t' '$3=="gene" || $3=="protein_coding_gene"' \
      "$WORK/inputs/annotations/${A}.fixed.gff3" | wc -l)
    MERGED_GENES=$(awk -F'\t' '$3=="gene" || $3=="protein_coding_gene"' \
      "$MERGED_GFF" | wc -l)
    THRESHOLD=$(( REF_GENES * 7 / 10 ))
    if [[ "$MERGED_GENES" -lt "$THRESHOLD" ]]; then
      log "WARNING: merged annotation ${A}→${Q} has only $MERGED_GENES genes" \
          "(expected >= $THRESHOLD from $REF_GENES reference genes)" >&2
    else
      log "  Merged ${A}→${Q}: $MERGED_GENES genes OK"
    fi

  done   # end Q loop
done     # end A loop

log "Phase C.4 complete"
