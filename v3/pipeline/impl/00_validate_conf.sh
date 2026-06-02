#!/usr/bin/env bash
set -euo pipefail
# Sanity-check species.conf before running any pipeline phase.
# Exit 0 only if all checks pass.

source "${WORK:-$(pwd)}/pipeline/species.conf"
cd "$WORK"

ERRORS=0
err() { echo "CONF ERROR: $1" >&2; ERRORS=$((ERRORS + 1)); }
ok()  { echo "OK:         $1"; }

# 1. REF_STRAIN must be in STRAINS
REF_FOUND=0
for S in "${STRAINS[@]}"; do [[ "$S" == "$REF_STRAIN" ]] && REF_FOUND=1 && break; done
[[ $REF_FOUND -eq 1 ]] && ok "REF_STRAIN=$REF_STRAIN is in STRAINS" \
                        || err "REF_STRAIN=$REF_STRAIN is NOT in STRAINS"

# 2. ANCHOR_STRAINS must be a subset of STRAINS
for A in "${ANCHOR_STRAINS[@]}"; do
  FOUND=0
  for S in "${STRAINS[@]}"; do [[ "$S" == "$A" ]] && FOUND=1 && break; done
  [[ $FOUND -eq 1 ]] && ok "ANCHOR_STRAIN=$A in STRAINS" \
                       || err "ANCHOR_STRAIN=$A is NOT in STRAINS"
done

# 3. N_CORES <= nproc
NPROC=$(nproc)
if [[ "$N_CORES" -le "$NPROC" ]]; then
  ok "N_CORES=$N_CORES <= nproc=$NPROC"
else
  err "N_CORES=$N_CORES exceeds nproc=$NPROC (will oversubscribe)"
fi

# 4. Assembly FASTA exists and is readable for each strain
for S in "${STRAINS[@]}"; do
  FA="$WORK/inputs/assemblies/${S}.fa"
  if [[ -s "$FA" ]]; then
    ok "Assembly readable: $S"
  else
    err "Assembly MISSING or empty: $FA"
  fi
done

# 5. Annotation GFF3 exists and is readable for each strain
for S in "${STRAINS[@]}"; do
  GFF="$WORK/inputs/annotations/${S}.gff3"
  if [[ -s "$GFF" ]]; then
    ok "Annotation readable: $S"
  else
    err "Annotation MISSING or empty: $GFF"
  fi
done

# 6. CHROM_RENAME file exists and is a 2-column TSV (if path is not a TODO)
if [[ "${CHROM_RENAME}" == TODO* ]]; then
  echo "WARN:       CHROM_RENAME is a TODO placeholder — Phase J will be skipped"
elif [[ -s "$CHROM_RENAME" ]]; then
  COLS=$(awk -F'\t' 'NR==1{print NF}' "$CHROM_RENAME")
  [[ "$COLS" -eq 2 ]] && ok "CHROM_RENAME is a 2-column TSV" \
                       || err "CHROM_RENAME has $COLS columns (expected 2)"
else
  err "CHROM_RENAME file MISSING or empty: $CHROM_RENAME"
fi

# 7. COHORT_VCF_DIR: warn only (no cohort yet for Pk)
if [[ "${COHORT_VCF_DIR}" == TODO* ]]; then
  echo "WARN:       COHORT_VCF_DIR is a TODO placeholder — Phase J will be skipped"
elif [[ -d "$COHORT_VCF_DIR" ]]; then
  ok "COHORT_VCF_DIR exists: $COHORT_VCF_DIR"
else
  err "COHORT_VCF_DIR MISSING: $COHORT_VCF_DIR"
fi

# 8. WORK directory is writable
if [[ -w "$WORK" ]]; then
  ok "WORK=$WORK is writable"
else
  err "WORK=$WORK is NOT writable"
fi

# 9. Container runtime available
if command -v docker &>/dev/null; then
  ok "docker is available: $(docker --version 2>&1 | head -1)"
elif command -v singularity &>/dev/null; then
  echo "WARN:       docker not found; singularity available — set RUN_IN_CONTAINER=apptainer"
elif command -v apptainer &>/dev/null; then
  echo "WARN:       docker not found; apptainer available — set RUN_IN_CONTAINER=apptainer"
else
  err "No container runtime found (docker, singularity, or apptainer required)"
fi

# 10. Minimum strain count
N="${#STRAINS[@]}"
if [[ $N -ge 5 ]]; then
  ok "Strain count $N (>= 5 required)"
else
  err "Only $N strains defined (need >= 5 for pangenome analysis)"
fi

echo ""
if [[ $ERRORS -eq 0 ]]; then
  echo "All configuration checks PASSED"
  exit 0
else
  echo "$ERRORS configuration check(s) FAILED"
  exit 1
fi
