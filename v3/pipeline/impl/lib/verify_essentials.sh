#!/usr/bin/env bash
set -uo pipefail
# Verify all *-essential outputs exist and pass structural validation.
# Exit 0 iff all 27 pass; exit 1 otherwise.
#
# Usage:
#   pipeline/lib/verify_essentials.sh            # full run
#   pipeline/lib/verify_essentials.sh --smoke     # smoke subset (fewer strains)

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

SMOKE=0
[[ "${1:-}" == "--smoke" ]] && SMOKE=1

MISSING=0

check() {
  local label="$1"
  local path="$2"
  local validator="$3"
  if [[ ! -e "$path" ]]; then
    echo "MISSING: $label  ($path)"
    MISSING=$((MISSING + 1))
    return
  fi
  if ! eval "$validator" 2>/dev/null; then
    echo "INVALID: $label  ($path)"
    MISSING=$((MISSING + 1))
    return
  fi
  echo "OK:      $label"
}

# ------------------------------------------------------------------
# Phase A — inventory
# ------------------------------------------------------------------
check "Mash distance matrix" \
  "work/00_inventory/mash/dist.tsv" \
  '[[ $(wc -l < work/00_inventory/mash/dist.tsv) -gt 1 ]]'

for S in "${STRAINS[@]}"; do
  check "BUSCO short_summary ${S}" \
    "work/00_inventory/busco/${S}_proteins" \
    '[[ -d work/00_inventory/busco/${S}_proteins ]]'
done

# ------------------------------------------------------------------
# Phase B — soft-masking
# ------------------------------------------------------------------
for S in "${STRAINS[@]}"; do
  check "Softmasked FASTA ${S}" \
    "genomes/softmasked/${S}.fa" \
    '[[ -s genomes/softmasked/${S}.fa.fai ]]'
  check "Sizes file ${S}" \
    "genomes/softmasked/${S}.sizes" \
    '[[ $(wc -l < genomes/softmasked/${S}.sizes) -ge 1 ]]'
done

# ------------------------------------------------------------------
# Phase C — pairwise AXT alignments
# ------------------------------------------------------------------
for i in "${!STRAINS[@]}"; do
  for j in "${!STRAINS[@]}"; do
    [[ $i -ge $j ]] && continue
    A="${STRAINS[$i]}"
    B="${STRAINS[$j]}"
    check "AXT ${A}__vs__${B}" \
      "projection/A2_kegalign/axt/${A}__vs__${B}.axt" \
      'grep -qE "^[0-9]+ " projection/A2_kegalign/axt/${A}__vs__${B}.axt'
  done
done

# Phase C — cleaned chains
for S1 in "${STRAINS[@]}"; do
  for S2 in "${STRAINS[@]}"; do
    [[ "$S1" == "$S2" ]] && continue
    check "Cleaned chain ${S1}→${S2}" \
      "work/01_chains/${S1}.${S2}.cleaned.chain" \
      'head -1 work/01_chains/${S1}.${S2}.cleaned.chain | awk '\''$1=="chain" && NF==13{ok=1} END{exit !ok}'\'''
  done
done

# Phase C — rbest chains (unordered pairs)
for i in "${!STRAINS[@]}"; do
  for j in "${!STRAINS[@]}"; do
    [[ $i -ge $j ]] && continue
    A="${STRAINS[$i]}"
    B="${STRAINS[$j]}"
    check "Rbest chain ${A}↔${B}" \
      "work/01_chains/${A}.${B}.rbest.chain" \
      'head -1 work/01_chains/${A}.${B}.rbest.chain | awk '\''$1=="chain"{ok=1} END{exit !ok}'\'''
  done
done

# Phase C — merged annotations
for A in "${ANCHOR_STRAINS[@]}"; do
  for Q in "${STRAINS[@]}"; do
    [[ "$Q" == "$A" ]] && continue
    check "Merged annotation ${A}→${Q}" \
      "work/02d_merged/${A}-as-ref/${Q}.annotation.gff3" \
      '[[ $(awk -F'\''\t'\'' '\''$3=="gene" || $3=="protein_coding_gene"'\'' work/02d_merged/${A}-as-ref/${Q}.annotation.gff3 | wc -l) -gt 100 ]]'
    check "Classification TSV ${A}→${Q}" \
      "work/02d_merged/${A}-as-ref/${Q}.classification.tsv" \
      '[[ $(wc -l < work/02d_merged/${A}-as-ref/${Q}.classification.tsv) -gt 100 ]]'
  done
done

# ------------------------------------------------------------------
# Phase D — PGGB
# ------------------------------------------------------------------
check "PGGB graph .og" \
  "inputs/pggb/pk.og" \
  '[[ -s inputs/pggb/pk.og ]]'
check "PGGB graph .gfa" \
  "inputs/pggb/pk.gfa" \
  '[[ -s inputs/pggb/pk.gfa ]]'

# ------------------------------------------------------------------
# Phase E — consensus orthology
# ------------------------------------------------------------------
check "Ortholog table" \
  "work/03_consensus/ortholog_table.tsv" \
  '[[ $(wc -l < work/03_consensus/ortholog_table.tsv) -gt 1000 ]]'

# ------------------------------------------------------------------
# Phase F — MSAs (check at least 50 files exist)
# ------------------------------------------------------------------
check "Core strict codon MSAs" \
  "work/06_msa/core_v3" \
  '[[ $(ls work/06_msa/core_v3/*.codon.aln.fa 2>/dev/null | wc -l) -ge 50 ]]'
check "Core relaxed codon MSAs" \
  "work/06_msa/core_relaxed" \
  '[[ $(ls work/06_msa/core_relaxed/*.codon.aln.fa 2>/dev/null | wc -l) -ge 50 ]]'

# ------------------------------------------------------------------
# Phase G — ML trees
# ------------------------------------------------------------------
check "IQ-TREE strict treefiles" \
  "work/06_msa/core_v3_trees" \
  '[[ $(find work/06_msa/core_v3_trees -name "*.treefile" 2>/dev/null | wc -l) -ge 50 ]]'
check "IQ-TREE relaxed treefiles" \
  "work/06_msa/core_relaxed_trees" \
  '[[ $(find work/06_msa/core_relaxed_trees -name "*.treefile" 2>/dev/null | wc -l) -ge 50 ]]'

# ------------------------------------------------------------------
# Phase H — HyPhy BUSTED
# ------------------------------------------------------------------
check "HyPhy strict BUSTED JSONs" \
  "work/06_msa/core_v3_hyphy/bulk" \
  '[[ $(find work/06_msa/core_v3_hyphy/bulk -name "busted.json" 2>/dev/null | wc -l) -ge 50 ]]'
check "HyPhy relaxed BUSTED JSONs" \
  "work/06_msa/core_relaxed_hyphy" \
  '[[ $(find work/06_msa/core_relaxed_hyphy -name "busted.json" 2>/dev/null | wc -l) -ge 50 ]]'

# ------------------------------------------------------------------
# Phase I — Multiz multi-way MAFs  (skip in smoke mode)
# ------------------------------------------------------------------
if [[ $SMOKE -eq 0 ]]; then
  for H in "${STRAINS[@]}"; do
    check "Multiz MAF hinge=${H}" \
      "work/07_multiz/${H}/${H}.multiz.maf" \
      'awk '\''/^a /{n++} END{exit !(n>=100)}'\'' work/07_multiz/${H}/${H}.multiz.maf'
  done
fi

# ------------------------------------------------------------------
# Phase J — Cohort VCF projection
# ------------------------------------------------------------------
if [[ "${COHORT_VCF_DIR}" == TODO* ]]; then
  echo "SKIP:    Phase J VCF projection (COHORT_VCF_DIR not yet set)"
else
  for TGT in "${STRAINS[@]}"; do
    [[ "$TGT" == "$REF_STRAIN" ]] && continue
    check "A2 cohort VCF ${TGT}" \
      "projection/A2_lastz/${SPECIES}_cohort_on_${TGT}.vcf.gz" \
      '[[ -s projection/A2_lastz/${SPECIES}_cohort_on_${TGT}.vcf.gz.csi ]]'
  done
fi

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
echo ""
if [[ $MISSING -eq 0 ]]; then
  echo "ALL ESSENTIALS VERIFIED (MISSING=0)"
  exit 0
else
  echo "${MISSING} essential outputs missing or invalid"
  exit 1
fi
