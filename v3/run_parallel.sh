#!/bin/bash
# Parallel pipeline for v3 projection.
# Step 1: rename remaining MalariaGEN VCFs (parallel xargs across chromosomes)
# Step 2 (parallel start): A1 CrossMap liftover 7 targets × 16 chrs; B BED build + odgi position
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
SCRATCH=/media/anton/scratch/malariagen_pv4
RENAME=$V3/projection/A1_wfmash/rename_pvp01_v1_to_genbank.tsv

REMAINING_CHRS="08 09 10 11 12 13 14 API MIT"
ALL_CHRS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 API MIT"
OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"

echo "===== Step 1: parallel rename of remaining $(echo $REMAINING_CHRS | wc -w) chromosomes ====="
rename_one() {
  local chr=$1
  docker run --rm -v /media/anton/data/sandbox/Pv4/v3:/v3 -v /media/anton/scratch/malariagen_pv4:/mg \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bcftools annotate --rename-chrs /v3/projection/A1_wfmash/rename_pvp01_v1_to_genbank.tsv \
                      /mg/Pv4_PvP01_${chr}_v1.vcf.gz -Oz \
                      -o /v3/projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz 2>/dev/null
  docker run --rm -v /media/anton/data/sandbox/Pv4/v3:/v3 \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bcftools index /v3/projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz 2>/dev/null
  echo "renamed $chr"
}
export -f rename_one
echo "$REMAINING_CHRS" | tr ' ' '\n' | xargs -P 9 -I {} bash -c 'rename_one "$@"' _ {}

echo
echo "===== Step 2a: build BED of MG sites from renamed VCFs (parallel) ====="
mkdir -p $V3/projection/B_graph/sites
build_bed_one() {
  local chr=$1
  docker run --rm -v /media/anton/data/sandbox/Pv4/v3:/v3 \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bcftools view /v3/projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz -f PASS 2>/dev/null \
  | awk -v cn="GCA_900093555.2#1#" '!/^#/{print cn $1 "#0\t" $2-1 "\t" $2}' \
  > /media/anton/data/sandbox/Pv4/v3/projection/B_graph/sites/mg_${chr}.bed
}
export -f build_bed_one
echo "$ALL_CHRS" | tr ' ' '\n' | xargs -P 16 -I {} bash -c 'build_bed_one "$@"' _ {}
cat $V3/projection/B_graph/sites/mg_*.bed | sort -k1,1 -k2,2n > $V3/projection/B_graph/sites/all_mg_sites.bed
echo "Total MG PASS sites: $(wc -l < $V3/projection/B_graph/sites/all_mg_sites.bed)"

echo
echo "===== Step 2b: A1 CrossMap liftover (112 jobs, parallel across 8 workers) ====="
mkdir -p $V3/projection/A1_wfmash/lifted
crossmap_one() {
  local OTHER=$1 chr=$2
  mkdir -p /media/anton/data/sandbox/Pv4/v3/projection/A1_wfmash/lifted/${OTHER}
  local out=/media/anton/data/sandbox/Pv4/v3/projection/A1_wfmash/lifted/${OTHER}/Pv4_${chr}_on_${OTHER}.vcf
  [ -s ${out}.gz ] && { echo "$OTHER $chr already done"; return; }
  docker run --rm -v /media/anton/data/sandbox/Pv4/v3:/v3 -w /v3 \
    quay.io/biocontainers/crossmap:0.7.3--pyhdfd78af_0 \
    CrossMap vcf projection/A1_wfmash/PvP01_to_${OTHER}.chain \
                 projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz \
                 genomes/softmasked/${OTHER}.fa \
                 projection/A1_wfmash/lifted/${OTHER}/Pv4_${chr}_on_${OTHER}.vcf 2>/dev/null
  echo "$OTHER $chr done"
}
export -f crossmap_one
# Generate all 112 jobs
JOBS=$(for o in $OTHERS; do for c in $ALL_CHRS; do echo "$o $c"; done; done)
echo "$JOBS" | xargs -P 8 -n 2 bash -c 'crossmap_one "$@"' _

echo
echo "===== Step 3a: A1 concat + index per target ====="
concat_one() {
  local OTHER=$1
  local dir=/media/anton/data/sandbox/Pv4/v3/projection/A1_wfmash/lifted/${OTHER}
  docker run --rm -v /media/anton/data/sandbox/Pv4/v3:/v3 \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "cd /v3/projection/A1_wfmash/lifted/${OTHER} && \
             for f in Pv4_*_on_${OTHER}.vcf; do [ -f \$f ] && bgzip -f \$f && tabix -p vcf \${f}.gz; done && \
             bcftools concat -a Pv4_*_on_${OTHER}.vcf.gz -Oz -o /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz 2>/dev/null && \
             bcftools index /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz"
  echo "concatenated $OTHER"
}
export -f concat_one
echo "$OTHERS" | tr ' ' '\n' | xargs -P 7 -I {} bash -c 'concat_one "$@"' _ {}

echo
echo "===== Step 3b: Path B odgi position (7 targets in parallel) ====="
OG=/media/anton/data/sandbox/Pv4/v2/pggb_out/pggb_input.fa.gz.f705205.c28ecf8.dcad0e6.smooth.fix.og
ODGI=/home/anton/miniconda3/envs/pggb/bin/odgi
odgi_one() {
  local OTHER=$1
  $ODGI position -t 4 \
       -i /media/anton/data/sandbox/Pv4/v2/pggb_out/pggb_input.fa.gz.f705205.c28ecf8.dcad0e6.smooth.fix.og \
       -b /media/anton/data/sandbox/Pv4/v3/projection/B_graph/sites/all_mg_sites.bed \
       -R /media/anton/data/sandbox/Pv4/v3/projection/B_graph/path_lists/${OTHER}.paths.txt \
       > /media/anton/data/sandbox/Pv4/v3/projection/B_graph/sites/mg_on_${OTHER}.tsv \
       2> /media/anton/data/sandbox/Pv4/v3/logs/B_odgi_${OTHER}.log
  n=$(wc -l < /media/anton/data/sandbox/Pv4/v3/projection/B_graph/sites/mg_on_${OTHER}.tsv)
  echo "odgi $OTHER: $n positions"
}
export ODGI
export -f odgi_one
echo "$OTHERS" | tr ' ' '\n' | xargs -P 4 -I {} bash -c 'odgi_one "$@"' _ {}

echo
echo "===== FINAL SUMMARY ====="
for OTHER in $OTHERS; do
  a1_n=$(docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
         bash -c "bcftools view /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz 2>/dev/null | grep -cv '^#'")
  b_n=$(wc -l < $V3/projection/B_graph/sites/mg_on_${OTHER}.tsv 2>/dev/null)
  printf "  %-22s  A1_lifted=%s  B_lifted=%s\n" "$OTHER" "${a1_n:-NA}" "${b_n:-NA}"
done
