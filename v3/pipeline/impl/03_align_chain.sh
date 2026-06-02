#!/usr/bin/env bash
set -euo pipefail
# Phase C.1-3 — KegAlign pairwise AXT + canonical UCSC chain pipeline + rbest chains.
#
# Outputs:
#   projection/A2_kegalign/axt/{A}__vs__{B}.axt   (N*(N-1)/2 unordered pairs)
#   work/01_chains/{src}.{tgt}.cleaned.chain        (N*(N-1) directed chains)
#   work/01_chains/{A}.{B}.rbest.chain              (N*(N-1)/2 unordered)

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

cmd() { "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" "$@"; }
log() { echo "$(date +%H:%M:%S) [$(basename "${BASH_SOURCE[0]}")] $1" >&2; }

mkdir -p "$WORK/projection/A2_kegalign/axt"
mkdir -p "$WORK/work/01_chains"
mkdir -p "$WORK/logs"

# ---------------------------------------------------------------------------
# 3.1 KegAlign GPU pairwise AXT (N*(N-1)/2 unordered pairs)
# ---------------------------------------------------------------------------
log "Phase C.1 — KegAlign pairwise alignment..."

for i in "${!STRAINS[@]}"; do
  for j in "${!STRAINS[@]}"; do
    [[ $i -ge $j ]] && continue
    A="${STRAINS[$i]}"
    B="${STRAINS[$j]}"
    AXT_OUT="$WORK/projection/A2_kegalign/axt/${A}__vs__${B}.axt"
    [[ -s "$AXT_OUT" ]] && continue

    log "  KegAlign $A vs $B..."
    if [[ "${USE_GPU:-0}" == "1" ]]; then
      # GPU mode
      USE_GPU=1 cmd kegalign \
        "$WORK/genomes/softmasked/${A}.fa" \
        "$WORK/genomes/softmasked/${B}.fa" \
        --strand both \
        --hsp_threshold 5000 \
        --gapped_threshold 6000 \
        --inner 2000 \
        --ydrop 15000 \
        --output_format axt \
        > "$AXT_OUT"
    else
      # CPU fallback via lastZ (from same ucsc-kent-tools image for axt output)
      # LastZ is not in the kent-tools image; use a dedicated lastz image.
      # Redefine cmd locally to override the case statement for lastz
      docker run --rm \
        -v "${WORK}:${WORK}" \
        -w "${WORK}" \
        -u "$(id -u):$(id -g)" \
        --entrypoint "" \
        "quay.io/biocontainers/lastz:1.04.52--h7b50bb2_1" \
        lastz \
          "${WORK}/genomes/softmasked/${A}.fa[multiple]" \
          "${WORK}/genomes/softmasked/${B}.fa[multiple]" \
          --masking=50 \
          --hspthresh=4500 \
          --gappedthresh=6000 \
          --inner=2000 \
          --ydrop=15000 \
          --format=axt \
        > "$AXT_OUT"
    fi

    # Validate AXT
    if ! grep -qE '^[0-9]+ ' "$AXT_OUT"; then
      log "ERROR: AXT validation failed for $A vs $B" >&2
      rm -f "$AXT_OUT"; exit 1
    fi
    log "  $A vs $B AXT OK"
  done
done

# ---------------------------------------------------------------------------
# 3.2 Chain pipeline — cleaned chains (N*(N-1) directed)
# ---------------------------------------------------------------------------
log "Phase C.2 — building cleaned chains..."

for AXT in "$WORK/projection/A2_kegalign/axt/"*.axt; do
  [[ -f "$AXT" ]] || continue
  PAIR=$(basename "$AXT" .axt)
  A="${PAIR%%__vs__*}"
  B="${PAIR##*__vs__}"

  # A->B: full chain build from the AXT (axt target = A, query = B)
  CHAIN_AB="$WORK/work/01_chains/${A}.${B}.cleaned.chain"
  if [[ ! -s "$CHAIN_AB" ]]; then
    log "  Chain ${A}->${B}..."
    NET_TMP=$(mktemp /tmp/Pk_net_XXXXXX.net)
    PRENET_TMP=$(mktemp /tmp/Pk_prenet_XXXXXX.chain)
    cmd axtChain -linearGap=loose -faT -faQ "$AXT" \
        "$WORK/genomes/softmasked/${A}.fa" \
        "$WORK/genomes/softmasked/${B}.fa" \
        /dev/stdout | \
      cmd chainSort stdin /dev/stdout | \
      cmd chainPreNet stdin \
        "$WORK/genomes/softmasked/${A}.sizes" \
        "$WORK/genomes/softmasked/${B}.sizes" \
        "$PRENET_TMP"
    cmd chainNet "$PRENET_TMP" \
        "$WORK/genomes/softmasked/${A}.sizes" \
        "$WORK/genomes/softmasked/${B}.sizes" \
        "$NET_TMP" /dev/null
    cmd netChainSubset "$NET_TMP" "$PRENET_TMP" /dev/stdout | \
      cmd chainStitchId /dev/stdin "$CHAIN_AB"
    rm -f "$NET_TMP" "$PRENET_TMP"
    if ! head -1 "$CHAIN_AB" | awk '$1=="chain" && NF==13 {ok=1} END{exit !ok}'; then
      log "ERROR: cleaned chain validation failed for ${A}->${B}" >&2
      rm -f "$CHAIN_AB"; exit 1
    fi
    log "  ${A}->${B} cleaned chain OK"
  fi

  # B->A: swap the A->B chain (AXT is directional; no re-alignment)
  CHAIN_BA="$WORK/work/01_chains/${B}.${A}.cleaned.chain"
  if [[ ! -s "$CHAIN_BA" ]]; then
    log "  Chain ${B}->${A} (chainSwap)..."
    cmd chainSwap "$CHAIN_AB" /dev/stdout | cmd chainSort stdin "$CHAIN_BA"
    if ! head -1 "$CHAIN_BA" | awk '$1=="chain" && NF==13 {ok=1} END{exit !ok}'; then
      log "ERROR: swapped chain validation failed for ${B}->${A}" >&2
      rm -f "$CHAIN_BA"; exit 1
    fi
    log "  ${B}->${A} cleaned chain OK"
  fi
done

# ---------------------------------------------------------------------------
# 3.3 rbest chains (N*(N-1)/2 unordered pairs)
# ---------------------------------------------------------------------------
log "Phase C.3 — building rbest chains..."

for i in "${!STRAINS[@]}"; do
  for j in "${!STRAINS[@]}"; do
    [[ $i -ge $j ]] && continue
    A="${STRAINS[$i]}"
    B="${STRAINS[$j]}"
    RBEST_OUT="$WORK/work/01_chains/${A}.${B}.rbest.chain"
    [[ -s "$RBEST_OUT" ]] && continue

    log "  rbest $A ↔ $B..."
    SWAP_CHAIN="$WORK/work/01_chains/${A}.${B}.swap.chain"
    SWAP_CLEANED="$WORK/work/01_chains/${A}.${B}.swap.cleaned.chain"
    NET_TMP2=$(mktemp /tmp/Pk_rbestnet_XXXXXX.net)
    trap 'rm -f "$NET_TMP2" "$SWAP_CHAIN" "$SWAP_CLEANED"' EXIT

    # Step 1: swap A→B chain to B→A
    cmd chainSwap "$WORK/work/01_chains/${A}.${B}.cleaned.chain" /dev/stdout | \
      cmd chainSort stdin "$SWAP_CHAIN"

    # Step 2: chainNet on swapped → B-target net
    cmd chainNet "$SWAP_CHAIN" \
        "$WORK/genomes/softmasked/${B}.sizes" \
        "$WORK/genomes/softmasked/${A}.sizes" \
        "$NET_TMP2" /dev/null
    cmd netChainSubset "$NET_TMP2" "$SWAP_CHAIN" /dev/stdout | \
      cmd chainStitchId /dev/stdin "$SWAP_CLEANED"

    # Step 3: swap back → rbest chain
    cmd chainSwap "$SWAP_CLEANED" /dev/stdout | \
      cmd chainSort stdin "$RBEST_OUT"

    rm -f "$SWAP_CHAIN" "$SWAP_CLEANED" "$NET_TMP2"

    # Validate
    if ! head -1 "$RBEST_OUT" | awk '$1=="chain" {ok=1} END{exit !ok}'; then
      log "ERROR: rbest chain validation failed for $A ↔ $B" >&2
      rm -f "$RBEST_OUT"; exit 1
    fi
    log "  $A ↔ $B rbest chain OK"
  done
done

log "Phase C.1-3 complete"
