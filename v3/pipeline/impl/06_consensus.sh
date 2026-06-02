#!/usr/bin/env bash
set -euo pipefail
# Phase E — Consensus orthology.
# Combines three evidence sources via union-find to produce ortholog_table.tsv.

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/work/03_consensus"
mkdir -p "$WORK/logs"

ORTHO_OUT="$WORK/work/03_consensus/ortholog_table.tsv"

if [[ -s "$ORTHO_OUT" ]]; then
  NLINES=$(wc -l < "$ORTHO_OUT")
  if [[ $NLINES -gt 1000 ]]; then
    log "Ortholog table already exists ($NLINES rows) — skipping"
    exit 0
  fi
fi

# ---------------------------------------------------------------------------
# 6.1 Per-strain gene BED (required by E.2 rbest + E.3 graph) — derive from GFF
# ---------------------------------------------------------------------------
for S in "${STRAINS[@]}"; do
  BED="$WORK/inputs/annotations/${S}.bed"
  [[ -s "$BED" ]] && continue
  GFF="$WORK/inputs/annotations/${S}.gff3"
  [[ -s "$GFF" ]] || GFF="$WORK/inputs/annotations/${S}.fixed.gff3"
  [[ -s "$GFF" ]] || { log "WARNING: no GFF for $S — cannot build BED" >&2; continue; }
  awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="gene"||$3=="protein_coding_gene"||$3=="ncRNA_gene"||$3=="pseudogene"{
      id=$9; sub(/.*ID=/,"",id); sub(/;.*/,"",id); sub(/^gene-/,"",id)
      if (id!="") print $1, $4-1, $5, id
    }' "$GFF" > "$BED"
  log "  Generated ${S}.bed ($(wc -l < "$BED") genes)"
done

# ---------------------------------------------------------------------------
# 6.2 Chain-based reciprocal-best edges
# ---------------------------------------------------------------------------
RBEST_EDGES="$WORK/work/03_consensus/rbest_edges.tsv"
if [[ ! -s "$RBEST_EDGES" ]]; then
  log "Phase E.2 — rbest chain overlap edges..."
  CHAIN_GLOB="$WORK/work/01_chains/*.rbest.chain"
  ANNOT_GLOB="$WORK/inputs/annotations/*.bed"
  cmd python3 "$WORK/pipeline/scripts/phase_e_rbest_overlap.py" \
    --chains "$CHAIN_GLOB" \
    --annotations "$ANNOT_GLOB" \
    --strains "${STRAINS[*]}" \
    --min_overlap 0.90 \
    --output "$RBEST_EDGES"
fi

# ---------------------------------------------------------------------------
# 6.3 Graph path co-membership edges
# ---------------------------------------------------------------------------
GRAPH_PATHS="$WORK/work/03_consensus/graph_paths.tsv"
GRAPH_EDGES="$WORK/work/03_consensus/graph_edges.tsv"

if [[ ! -s "$GRAPH_PATHS" ]]; then
  log "Phase E.3 — odgi path haplotypes..."
  cmd odgi paths -i "$WORK/inputs/pggb/pk.og" --haplotypes \
    > "$GRAPH_PATHS"
fi

if [[ ! -s "$GRAPH_EDGES" ]]; then
  log "Phase E.3 — graph co-membership edges..."
  ANNOT_GLOB="$WORK/inputs/annotations/*.bed"
  cmd python3 "$WORK/pipeline/scripts/phase_e_graph_edges.py" \
    --paths "$GRAPH_PATHS" \
    --annotations "$ANNOT_GLOB" \
    --strains "${STRAINS[*]}" \
    --output "$GRAPH_EDGES"
fi

# ---------------------------------------------------------------------------
# 6.4 Union-find consensus
# ---------------------------------------------------------------------------
log "Phase E.4 — union-find consensus orthogroups..."
cmd python3 "$WORK/pipeline/scripts/phase_e_consensus.py" \
  --liftoff_dir "$WORK/work/02d_merged" \
  --rbest "$RBEST_EDGES" \
  --graph "$GRAPH_EDGES" \
  --anchors "${ANCHOR_STRAINS[*]}" \
  --strains "${STRAINS[*]}" \
  --ref "$REF_STRAIN" \
  --output "$ORTHO_OUT"

# Validate
NLINES=$(wc -l < "$ORTHO_OUT")
if [[ $NLINES -le 1000 ]]; then
  log "ERROR: ortholog_table.tsv has only $NLINES lines (expected >1000)" >&2
  exit 1
fi
log "Ortholog table OK — $NLINES orthogroups"

log "Phase E complete"
