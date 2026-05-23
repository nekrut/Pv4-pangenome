#!/bin/bash
# Fix A1 cohort concats: re-tabix any .vcf.gz missing .tbi, then bcftools concat
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
IMG=quay.io/biocontainers/bcftools:1.20--h8b25389_0

fix_one() {
  local OTHER=$1
  local dir=$V3/projection/A1_wfmash/lifted/${OTHER}
  local cohort=$V3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz
  # if cohort already big, skip
  if [ -s $cohort ] && [ $(stat -c%s $cohort) -gt 100000000 ]; then
    echo "[skip] $OTHER cohort already $(du -h $cohort | cut -f1)"; return
  fi
  echo "[fix] $OTHER"
  docker run --rm -v $V3:/v3 $IMG bash -c "
    cd /v3/projection/A1_wfmash/lifted/${OTHER}
    # Sort each per-chr file (CrossMap may emit out-of-order on chain inversions)
    for f in Pv4_*_on_${OTHER}.vcf.gz; do
      bn=\$(basename \$f .vcf.gz)
      if [ ! -s \${bn}.sorted.vcf.gz ]; then
        bcftools sort -Oz -o \${bn}.sorted.vcf.gz \$f 2>/dev/null
        tabix -p vcf -f \${bn}.sorted.vcf.gz
      fi
    done
    rm -f /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz*
    bcftools concat -a Pv4_*_on_${OTHER}.sorted.vcf.gz \
      -Oz -o /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz 2>&1 | tail -5
    bcftools index /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz
    # Clean up sorted intermediates
    rm -f Pv4_*_on_${OTHER}.sorted.vcf.gz Pv4_*_on_${OTHER}.sorted.vcf.gz.tbi
  "
  if [ -s $cohort ]; then
    sz=$(du -h $cohort | cut -f1)
    echo "[done] $OTHER cohort=$sz"
  else
    echo "[FAIL] $OTHER"
  fi
}
export V3 IMG
export -f fix_one

OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_949152365.1 GCA_040114635.1"
echo "$OTHERS" | tr ' ' '\n' | xargs -P 5 -I {} bash -c 'fix_one "$@"' _ {}

echo
echo "===== A1 final cohort sizes ====="
for o in GCA_000002415.2 GCA_900093555.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1; do
  c=$V3/projection/A1_wfmash/Pv4_cohort_on_${o}.vcf.gz
  if [ -s $c ]; then printf "  %-22s  %s\n" "$o" "$(du -h $c | cut -f1)"; fi
done
