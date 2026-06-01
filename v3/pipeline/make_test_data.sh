#!/usr/bin/env bash
# make_test_data.sh — materialize the Pv4 TEST PANEL (5 strains × 3 chromosomes).
#
# Builds a small, biology-complete subset of the full Pv4 v3 inputs that
# exercises every pipeline phase (A–K) end-to-end in ~2–3 h. See
# test_data/README.md and pipeline/PIPELINE_EXPLANATION.md.
#
# Host needs only: docker, bash, awk, gunzip, (rclone for the cohort pull).
# All bio-tools run in pinned biocontainers — no local installs.
#
# Usage:
#   bash v3/pipeline/make_test_data.sh [--stage seq|vcf|all] [--cohort-src DIR]
#     --stage seq   assemblies + annotations + proteomes + chrom map (local data; no network)
#     --stage vcf   cohort VCF subset only (needs the MalariaGEN per-chr files)
#     --stage all   both (default)
#     --cohort-src  dir holding Pv4_PvP01_{04,05,14}_v1.vcf.gz (PvP01 coords).
#                   If absent, the VCF stage is skipped with a notice.
set -euo pipefail

# --- paths -------------------------------------------------------------------
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)   # repo root
SRC=$REPO/v3                                                # full v3 inputs
OUT=$REPO/v3/pipeline/test_data                             # bundle lands here
MAP=$OUT/contig_map.tsv                                     # validated, committed
TMP=$(mktemp -d -p "$REPO" .make_test_data.XXXXXX)   # under $REPO so containers can see it
trap 'rm -rf "$TMP"' EXIT

STAGE=all; COHORT_SRC=""
while [[ $# -gt 0 ]]; do case $1 in
  --stage) STAGE=$2; shift 2;;
  --cohort-src) COHORT_SRC=$2; shift 2;;
  *) echo "unknown arg: $1" >&2; exit 1;;
esac; done

# --- panel -------------------------------------------------------------------
STRAINS=(PvP01 PvW1 PAM PvT01 MHC087)
declare -A GCA=( [PvP01]=GCA_900093555.2 [PvW1]=GCA_914969965.1 [PAM]=GCA_949152365.1
                 [PvT01]=GCA_900093545.1 [MHC087]=GCA_040114635.1 )
# strain -> source GFF3 (anchors: fixed/genbank; non-anchors: NCBI)
declare -A GFF=( [PvP01]=$SRC/inputs/annotations/PvP01.genbank.gff3.gz
                 [PvW1]=$SRC/inputs/annotations/PvW1.fixed.gff3.gz
                 [PAM]=$SRC/inputs/annotations/PAM.fixed.gff3.gz
                 [PvT01]=$SRC/annotations/ncbi/GCA_900093545.1.gff.gz
                 [MHC087]=$SRC/annotations/ncbi/GCA_040114635.1.liftoff.gff.gz )
# PvP01 test chromosomes (cohort uses the _v1 PlasmoDB names)
TEST_CHR_ACC=(LT635615.1 LT635616.2 LT635625.2)          # chr04 chr05 chr14
declare -A CHR_V1=( [LT635615.1]=PvP01_04_v1 [LT635616.2]=PvP01_05_v1 [LT635625.2]=PvP01_14_v1 )
N_SAMPLES=200

# --- containers --------------------------------------------------------------
SAMTOOLS=quay.io/biocontainers/samtools:1.20--h50ea8bc_0
BCFTOOLS=quay.io/biocontainers/bcftools:1.20--h8b25389_0
GFFREAD=quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4
dk(){ local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" \
        -v "$REPO:$REPO" -w "$REPO" --entrypoint "" "$img" "$@"; }

log(){ echo "$(date +%H:%M:%S) [make_test_data] $*" >&2; }

mkdir -p "$OUT/inputs/assemblies" "$OUT/inputs/annotations" \
         "$OUT/inputs/proteomes" "$OUT/inputs/cohort_vcf"

# Unique contigs for a strain (PAM merges chr04+chr05 onto one contig → dedupe).
contigs_for(){ awk -v s="$1" '$1==s && $0!~/^#/ {print $4}' "$MAP" | sort -u; }

# ----------------------------------------------------------------------------
seq_stage(){
  log "STAGE seq — assemblies + annotations + proteomes + chrom map"
  for S in "${STRAINS[@]}"; do
    local gca=${GCA[$S]} raw=$SRC/genomes/raw/${GCA[$S]}.fa.gz
    [[ -s $raw ]] || { log "MISSING assembly $raw"; exit 1; }
    mapfile -t CTG < <(contigs_for "$S")
    log "$S: contigs ${CTG[*]}"

    # assembly: gunzip (NCBI gzip != bgzip), index, extract the kept contigs
    gunzip -c "$raw" > "$TMP/$S.full.fa"
    dk "$SAMTOOLS" samtools faidx "$TMP/$S.full.fa"
    dk "$SAMTOOLS" samtools faidx "$TMP/$S.full.fa" "${CTG[@]}" \
        > "$OUT/inputs/assemblies/$S.fa"
    dk "$SAMTOOLS" samtools faidx "$OUT/inputs/assemblies/$S.fa"

    # annotation: keep features on the kept contigs (+ headers)
    printf '%s\n' "${CTG[@]}" > "$TMP/$S.contigs.txt"
    { [[ ${GFF[$S]} == *.gz ]] && gunzip -c "${GFF[$S]}" || cat "${GFF[$S]}"; } | \
      awk -F'\t' 'NR==FNR{keep[$1]=1; next} /^#/{print; next} ($1 in keep)' \
        "$TMP/$S.contigs.txt" - > "$OUT/inputs/annotations/$S.gff3"

    # proteome: regenerate from the subset GFF + subset assembly (only kept genes)
    dk "$GFFREAD" gffread "$OUT/inputs/annotations/$S.gff3" \
        -g "$OUT/inputs/assemblies/$S.fa" -y "$OUT/inputs/proteomes/$S.proteins.fa" \
      || log "$S: gffread proteome step warned (check GFF feature types)"
  done

  # cohort chromosome-rename map (3 rows; cohort is in _v1 PlasmoDB coords)
  : > "$OUT/inputs/annotations/cohort_chrom_rename.tsv"
  for acc in "${TEST_CHR_ACC[@]}"; do
    printf '%s\t%s\n' "${CHR_V1[$acc]}" "$acc" >> "$OUT/inputs/annotations/cohort_chrom_rename.tsv"
  done
  log "wrote cohort_chrom_rename.tsv"
}

# ----------------------------------------------------------------------------
vcf_stage(){
  if [[ -z $COHORT_SRC || ! -d $COHORT_SRC ]]; then
    log "STAGE vcf SKIPPED — no --cohort-src. Provide the 3 PvP01-coord per-chr"
    log "  MalariaGEN files (Pv4_PvP01_{04,05,14}_v1.vcf.gz) and re-run with --stage vcf."
    return 0
  fi
  log "STAGE vcf — downsample to ~$N_SAMPLES samples on chr04/05/14"
  local first=$COHORT_SRC/Pv4_PvP01_04_v1.vcf.gz
  [[ -s $first ]] || { log "MISSING $first"; exit 1; }
  # deterministic ~N-sample spread across the cohort
  local total; total=$(dk "$BCFTOOLS" bcftools query -l "$first" | wc -l)
  local step=$(( total / N_SAMPLES )); (( step < 1 )) && step=1
  dk "$BCFTOOLS" bcftools query -l "$first" | sort | awk -v k="$step" 'NR%k==1' \
      > "$OUT/samples_200.txt"
  log "selected $(wc -l < "$OUT/samples_200.txt") of $total samples (every ${step}th)"

  for acc in "${TEST_CHR_ACC[@]}"; do
    local v1=${CHR_V1[$acc]} tag=${CHR_V1[$acc]#PvP01_}; tag=${tag%_v1}
    local in=$COHORT_SRC/Pv4_PvP01_${tag}_v1.vcf.gz
    local out=$OUT/inputs/cohort_vcf/Pv4test_${tag}_v1.vcf.gz
    [[ -s $in ]] || { log "MISSING $in"; exit 1; }
    dk "$BCFTOOLS" bcftools view -S "$OUT/samples_200.txt" --force-samples \
        -Oz -o "$out" "$in"
    dk "$BCFTOOLS" bcftools index -t "$out"
    log "wrote $(basename "$out")"
  done
  # QC: the drug-resistance / selection landmark sites should remain polymorphic
  log "QC — verify dhps/dhfr/Pvs230 sites still segregate (manual: bcftools view -H ... | check AC>0)"
}

[[ $STAGE == seq || $STAGE == all ]] && seq_stage
[[ $STAGE == vcf || $STAGE == all ]] && vcf_stage

log "done. Bundle: $OUT  (FASTA + GFF + proteomes [+ cohort VCF if built])"
log "Upload heavy parts:  rclone copy $OUT dropbox:Pv4_v3/test_data/ --transfers 4"
