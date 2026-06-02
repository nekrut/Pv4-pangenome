#!/usr/bin/env bash
set -euo pipefail
# Orchestrator: run all 11 phases in dependency order.
# Phases C.1-3 (align+chain) and C.4 (annotation) run in parallel.
# Phase D (PGGB) starts after B; Phase E starts after C.4 + D.
# Phases G and H are sequential (H needs G's trees).
# Phase I (multiz) depends only on C.1-3 AXTs.
# Phase J (VCF projection) depends on C.2 cleaned chains.
#
# Resume-safe: every phase is idempotent — re-running after crash continues
# where it left off.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${WORK:-$SCRIPT_DIR/..}/pipeline/species.conf"
cd "$WORK"

mkdir -p "$WORK/logs"

log() {
  MSG="$(date +%H:%M:%S) $1"
  echo "$MSG"
  echo "$MSG" >> "$WORK/logs/orchestrator.log"
}

# Validate config first
log "Validating species.conf..."
bash "$WORK/pipeline/00_validate_conf.sh"

log "=== Pk v1 pipeline starting ==="
log "  SPECIES=$SPECIES  STRAINS=${STRAINS[*]}"
log "  REF_STRAIN=$REF_STRAIN  N_CORES=$N_CORES"

# Phase A — inventory
log "Phase A — inventory"
bash "$WORK/pipeline/01_inventory.sh" > "$WORK/logs/01_inventory.log" 2>&1
log "Phase A DONE"

# Phase B — soft masking
log "Phase B — soft masking"
bash "$WORK/pipeline/02_mask.sh" > "$WORK/logs/02_mask.log" 2>&1
log "Phase B DONE"

# Phases C.1-3 (align+chain) and C.4 (annotation) in parallel
log "Phase C.1-3 — align + chain [background]"
bash "$WORK/pipeline/03_align_chain.sh" > "$WORK/logs/03_align_chain.log" 2>&1 &
PID_C=$!

log "Phase C.4 — annotation projection [background]"
bash "$WORK/pipeline/04_annotate_project.sh" > "$WORK/logs/04_annotate_project.log" 2>&1 &
PID_ANNOT=$!

# Phase D — PGGB (depends only on B; can start immediately)
log "Phase D — PGGB pangenome [background]"
bash "$WORK/pipeline/05_pggb.sh" > "$WORK/logs/05_pggb.log" 2>&1 &
PID_D=$!

# Wait for C.1-3 (needed for Phase I and J)
wait "$PID_C"
log "Phase C.1-3 DONE"

# Phase I — Multiz (depends on C.1-3 AXTs; independent of E/F/G/H)
log "Phase I — Multiz multi-way MAFs [background]"
bash "$WORK/pipeline/10_multiz.sh" > "$WORK/logs/10_multiz.log" 2>&1 &
PID_I=$!

# Phase J — VCF projection (depends on C.2 chains)
log "Phase J — cohort VCF projection [background]"
bash "$WORK/pipeline/11_project_vcf.sh" > "$WORK/logs/11_project_vcf.log" 2>&1 &
PID_J=$!

# Wait for C.4 annotation before Phase E
wait "$PID_ANNOT"
log "Phase C.4 DONE"

# Wait for D before Phase E
wait "$PID_D"
log "Phase D DONE"

# Phase E — consensus orthology (depends on C.4 + D)
log "Phase E — consensus orthology"
bash "$WORK/pipeline/06_consensus.sh" > "$WORK/logs/06_consensus.log" 2>&1
log "Phase E DONE"

# Phase F — MSAs (depends on E + D + C.4)
log "Phase F — codon + protein MSAs"
bash "$WORK/pipeline/07_msa.sh" > "$WORK/logs/07_msa.log" 2>&1
log "Phase F DONE"

# Phase G — ML trees (depends on F)
log "Phase G — IQ-TREE3 ML trees"
bash "$WORK/pipeline/08_trees.sh" > "$WORK/logs/08_trees.log" 2>&1
log "Phase G DONE"

# Phase H — HyPhy BUSTED (depends on G)
log "Phase H — HyPhy BUSTED selection scan"
bash "$WORK/pipeline/09_hyphy.sh" > "$WORK/logs/09_hyphy.log" 2>&1
log "Phase H DONE"

# Wait for I and J
wait "$PID_I"
log "Phase I DONE"
wait "$PID_J"
log "Phase J DONE"

log "=== ALL PHASES COMPLETE ==="

# Final verification
log "Running verify_essentials.sh..."
bash "$WORK/pipeline/lib/verify_essentials.sh"
