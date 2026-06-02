#!/usr/bin/env bash
set -euo pipefail
# Rewrite the GFF3 chromosome (column 1) to match the assembly FASTA sequence IDs.
#
# Usage:
#   fix_gff_chroms.sh <STRAIN>
#
# Reads:
#   $WORK/inputs/assemblies/${STRAIN}.fa       — assembly FASTA (source of truth)
#   $WORK/inputs/annotations/${STRAIN}.gff3    — original GFF (PlasmoDB or other naming)
#
# Writes:
#   $WORK/inputs/annotations/${STRAIN}.fixed.gff3   — GFF with renamed chroms
#
# Strategy:
#   1. Extract FASTA chrom IDs (first whitespace-delimited token after '>').
#   2. Extract GFF chrom IDs.
#   3. Build a rename map by matching stripped versions (remove '.1', 'chr', etc.).
#   4. Apply via awk in-place rewrite.
#
# If the GFF chroms already match the FASTA, the script copies unchanged.

# -----------------------------------------------------------------------
# Source config if WORK is not set
# -----------------------------------------------------------------------
if [[ -z "${WORK:-}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  source "${SCRIPT_DIR}/../species.conf"
fi

STRAIN="${1:?Usage: fix_gff_chroms.sh <STRAIN>}"

FASTA="$WORK/inputs/assemblies/${STRAIN}.fa"
GFF_ORIG="$WORK/inputs/annotations/${STRAIN}.gff3"
GFF_FIXED="$WORK/inputs/annotations/${STRAIN}.fixed.gff3"

if [[ ! -f "$FASTA" ]]; then
  echo "ERROR: FASTA not found: $FASTA" >&2; exit 1
fi
if [[ ! -f "$GFF_ORIG" ]]; then
  echo "ERROR: GFF not found: $GFF_ORIG" >&2; exit 1
fi

if [[ -s "$GFF_FIXED" ]]; then
  echo "[$STRAIN] fixed GFF already exists — skipping (delete to rerun)"
  exit 0
fi

RENAME_TSV=$(mktemp /tmp/chrmap_XXXXXX.tsv)
trap 'rm -f "$RENAME_TSV"' EXIT

# Build set of FASTA chrom IDs
FA_CHROMS=$(mktemp /tmp/fa_chroms_XXXXXX.txt)
trap 'rm -f "$RENAME_TSV" "$FA_CHROMS"' EXIT
grep '^>' "$FASTA" | awk '{print substr($1,2)}' | sort -u > "$FA_CHROMS"

# Build set of GFF chrom IDs
GFF_CHROMS=$(mktemp /tmp/gff_chroms_XXXXXX.txt)
trap 'rm -f "$RENAME_TSV" "$FA_CHROMS" "$GFF_CHROMS"' EXIT
awk -F'\t' '!/^#/ && NF>=9 {print $1}' "$GFF_ORIG" | sort -u > "$GFF_CHROMS"

# Check if GFF chroms already match FASTA
MATCHES=$(comm -12 "$FA_CHROMS" "$GFF_CHROMS" | wc -l)
GFF_TOTAL=$(wc -l < "$GFF_CHROMS")

if [[ "$MATCHES" -eq "$GFF_TOTAL" && "$GFF_TOTAL" -gt 0 ]]; then
  echo "[$STRAIN] GFF chroms already match FASTA — copying as-is"
  cp "$GFF_ORIG" "$GFF_FIXED"
  exit 0
fi

echo "[$STRAIN] $MATCHES/$GFF_TOTAL GFF chroms match — building rename map"

# Build rename map: for each GFF chrom, try several normalizations to find
# the matching FASTA chrom.
python3 - "$GFF_CHROMS" "$FA_CHROMS" "$RENAME_TSV" <<'PYEOF'
import sys, re

gff_file, fa_file, out_file = sys.argv[1:]

def normalize(s):
    """Strip common prefixes/suffixes for fuzzy matching."""
    s = s.strip()
    s = re.sub(r'\.1$', '', s)   # remove trailing .1
    s = re.sub(r'^chr', '', s, flags=re.IGNORECASE)  # remove chr prefix
    s = re.sub(r'^chromosome_?', '', s, flags=re.IGNORECASE)
    return s.lower()

fa_chroms = [l.strip() for l in open(fa_file) if l.strip()]
gff_chroms = [l.strip() for l in open(gff_file) if l.strip()]

# Build normalized → original dict for FA
fa_norm = {normalize(c): c for c in fa_chroms}

rename = {}
unmatched = []
for gc in gff_chroms:
    key = normalize(gc)
    if key in fa_norm:
        rename[gc] = fa_norm[key]
    else:
        # Try appending .1
        key2 = normalize(gc + '.1')
        if key2 in fa_norm:
            rename[gc] = fa_norm[key2]
        else:
            unmatched.append(gc)

with open(out_file, 'w') as fh:
    for old, new in rename.items():
        fh.write(f"{old}\t{new}\n")

if unmatched:
    print(f"WARNING: {len(unmatched)} GFF chroms unmatched:", file=sys.stderr)
    for u in unmatched[:10]:
        print(f"  {u}", file=sys.stderr)

print(f"Rename map: {len(rename)} entries ({len(unmatched)} unmatched)")
PYEOF

if [[ ! -s "$RENAME_TSV" ]]; then
  echo "[$STRAIN] ERROR: rename map is empty — cannot fix GFF chroms" >&2
  echo "[$STRAIN] Inspect $GFF_CHROMS vs $FA_CHROMS manually" >&2
  exit 1
fi

echo "[$STRAIN] Applying rename map to GFF..."
awk -F'\t' 'BEGIN{OFS="\t"; while((getline < ARGV[1])>0){map[$1]=$2} ARGV[1]=""}
  /^#/ {print; next}
  NF>=9 {if($1 in map) $1=map[$1]; print}
  NF<9  {print}
' "$RENAME_TSV" "$GFF_ORIG" > "$GFF_FIXED"

N_ORIG=$(awk -F'\t' '!/^#/ && $3=="gene"' "$GFF_ORIG" | wc -l)
N_FIXED=$(awk -F'\t' '!/^#/ && $3=="gene"' "$GFF_FIXED" | wc -l)
echo "[$STRAIN] genes before=$N_ORIG after=$N_FIXED"

if [[ "$N_FIXED" -lt "$((N_ORIG / 2))" ]]; then
  echo "[$STRAIN] ERROR: gene count dropped by >50% after rename — check map" >&2
  rm -f "$GFF_FIXED"
  exit 1
fi

echo "[$STRAIN] Fixed GFF written to $GFF_FIXED"
