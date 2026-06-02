#!/usr/bin/env bash
set -euo pipefail
# Download P. knowlesi assemblies + annotations via NCBI Datasets CLI.
#
# Requirements: `datasets` CLI (ncbi-datasets-cli) must be in PATH.
#   Install: conda install -c conda-forge ncbi-datasets-cli
#   Or: pip install ncbi-datasets-cli
#
# Usage:
#   bash pipeline/setup/fetch_assemblies.sh
#
# Downloads FASTA + GFF3 for each strain defined in species.conf into:
#   inputs/assemblies/{S}.fa
#   inputs/annotations/{S}.gff3
#
# Also extracts proteomes to inputs/proteomes/{S}.proteins.fa via gffread.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${WORK:-$SCRIPT_DIR/../..}/pipeline/species.conf"
cd "$WORK"

log() { echo "$(date +%H:%M:%S) [fetch_assemblies] $1"; }

mkdir -p "$WORK/inputs/assemblies"
mkdir -p "$WORK/inputs/annotations"
mkdir -p "$WORK/inputs/proteomes"
mkdir -p "$WORK/logs"

if ! command -v datasets &>/dev/null; then
  echo "ERROR: 'datasets' CLI not found." >&2
  echo "Install via: conda install -c conda-forge ncbi-datasets-cli" >&2
  echo "Or visit: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/" >&2
  exit 1
fi

# Strain → GenBank accession mapping
# Must match STRAIN_ACCESSION in species.conf
declare -A ACCESSION=(
    [Pk_H]="GCA_000006355.3"
    [Pk_ANKA]="GCA_009792815.1"
    [Pk_A1H1]="GCA_014858965.1"
    [Pk_YH1]="GCA_014858985.1"
    [Pk_COF]="GCA_019833925.1"
    [Pk_TAM]="GCA_021201615.1"
    [Pk_SIM]="GCA_963506685.1"
)

for S in "${STRAINS[@]}"; do
  ACC="${ACCESSION[$S]:-}"
  if [[ -z "$ACC" ]]; then
    log "WARNING: No accession for $S — skipping"
    continue
  fi

  FA_OUT="$WORK/inputs/assemblies/${S}.fa"
  GFF_OUT="$WORK/inputs/annotations/${S}.gff3"

  if [[ -s "$FA_OUT" ]] && [[ -s "$GFF_OUT" ]]; then
    log "$S already downloaded — skipping"
    continue
  fi

  log "Downloading $S ($ACC)..."
  TMP_DIR=$(mktemp -d /tmp/ncbi_Pk_${S}_XXXXXX)
  trap 'rm -rf "$TMP_DIR"' EXIT

  # datasets download genome accession GCA_XXXXXXX.N --include genome,gff3
  datasets download genome accession "$ACC" \
    --include genome,gff3 \
    --filename "$TMP_DIR/ncbi_dataset.zip"

  unzip -q "$TMP_DIR/ncbi_dataset.zip" -d "$TMP_DIR/"

  # Locate FASTA (usually in ncbi_dataset/data/{ACC}/*.fna)
  FASTA_FILE=$(find "$TMP_DIR/ncbi_dataset/data/${ACC}/" -name "*.fna" | head -1)
  GFF_FILE=$(find "$TMP_DIR/ncbi_dataset/data/${ACC}/" -name "*.gff" -o -name "*.gff3" | head -1)

  if [[ -z "$FASTA_FILE" ]]; then
    log "ERROR: No FASTA found for $S ($ACC)" >&2
    continue
  fi
  if [[ -z "$GFF_FILE" ]]; then
    log "WARNING: No GFF found for $S ($ACC) — annotation will be missing" >&2
  fi

  cp "$FASTA_FILE" "$FA_OUT"
  log "  Copied FASTA → $FA_OUT"

  if [[ -n "$GFF_FILE" ]]; then
    cp "$GFF_FILE" "$GFF_OUT"
    log "  Copied GFF3 → $GFF_OUT"
  fi

  rm -rf "$TMP_DIR"
  trap - EXIT

  log "$S download complete"
done

log ""
log "Post-processing: building fixed GFF + BED12 + isoforms..."

for S in "${STRAINS[@]}"; do
  FA="$WORK/inputs/assemblies/${S}.fa"
  GFF="$WORK/inputs/annotations/${S}.gff3"
  [[ ! -s "$FA" || ! -s "$GFF" ]] && continue

  # Fix GFF chromosome names to match FASTA
  if [[ ! -s "$WORK/inputs/annotations/${S}.fixed.gff3" ]]; then
    log "  fix_gff_chroms.sh $S..."
    bash "$WORK/pipeline/lib/fix_gff_chroms.sh" "$S"
  fi

  # Extract proteome via gffread (containerized)
  PROT="$WORK/inputs/proteomes/${S}.proteins.fa"
  if [[ ! -s "$PROT" ]]; then
    log "  Extracting proteome for $S..."
    "$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh" gffread \
      "$WORK/inputs/annotations/${S}.fixed.gff3" \
      -g "$FA" \
      -y "$PROT" || true
  fi
done

# Build BED12 + isoforms for anchor strains
log "Building anchor BED12 + isoforms..."
for A in "${ANCHOR_STRAINS[@]}"; do
  if [[ -s "$WORK/inputs/annotations/${A}.bed12" ]]; then
    log "  $A BED12 already exists — skipping"
    continue
  fi
  if [[ -s "$WORK/inputs/assemblies/${A}.fa" ]]; then
    bash "$WORK/pipeline/setup/build_anchor_inputs.sh" "$A"
  fi
done

log ""
log "=== Download summary ==="
for S in "${STRAINS[@]}"; do
  FA="$WORK/inputs/assemblies/${S}.fa"
  GFF="$WORK/inputs/annotations/${S}.gff3"
  FIXED="$WORK/inputs/annotations/${S}.fixed.gff3"
  FA_STATUS=$([[ -s "$FA"    ]] && echo "OK" || echo "MISSING")
  GF_STATUS=$([[ -s "$GFF"   ]] && echo "OK" || echo "MISSING")
  FX_STATUS=$([[ -s "$FIXED" ]] && echo "OK" || echo "MISSING")
  log "  $S: FASTA=$FA_STATUS  GFF=$GF_STATUS  fixed.gff3=$FX_STATUS"
done

log "fetch_assemblies.sh complete"
