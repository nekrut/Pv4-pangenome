#!/bin/bash
set -euo pipefail
V3=/media/anton/data/sandbox/Pv4/v3
SCRATCH=/media/anton/scratch/malariagen_pv4
RENAME=$V3/projection/A1_wfmash/rename_pvp01_v1_to_genbank.tsv
CHROMS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 API MIT"
OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"

# Step 1: rename each MG VCF's CHROM once (16 files)
mkdir -p $V3/projection/A1_wfmash/mg_renamed
for chr in $CHROMS; do
  out=$V3/projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz
  [ -s $out ] && continue
  docker run --rm -v $V3:/v3 -v $SCRATCH:/mg quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bcftools annotate --rename-chrs /v3/projection/A1_wfmash/rename_pvp01_v1_to_genbank.tsv \
                      /mg/Pv4_PvP01_${chr}_v1.vcf.gz -Oz -o /v3/projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz
  docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bcftools index /v3/projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz
  echo "renamed Pv4_${chr}.vcf.gz"
done

# Step 2: per target, liftover all 16 chromosomes then concat
for OTHER in $OTHERS; do
  echo "===== liftover to $OTHER ====="
  outdir=$V3/projection/A1_wfmash/lifted/${OTHER}
  mkdir -p $outdir
  for chr in $CHROMS; do
    [ -s $outdir/Pv4_${chr}_on_${OTHER}.vcf ] && continue
    docker run --rm -v $V3:/v3 -w /v3 quay.io/biocontainers/crossmap:0.7.3--pyhdfd78af_0 \
      CrossMap vcf projection/A1_wfmash/PvP01_to_${OTHER}.chain \
                   projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz \
                   genomes/softmasked/${OTHER}.fa \
                   projection/A1_wfmash/lifted/${OTHER}/Pv4_${chr}_on_${OTHER}.vcf 2>&1 | tail -3
  done
  # Concat
  docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "cd /v3/projection/A1_wfmash/lifted/${OTHER} && \
             for f in Pv4_*_on_${OTHER}.vcf; do bgzip -f \$f; tabix -p vcf \${f}.gz; done && \
             bcftools concat -a Pv4_*_on_${OTHER}.vcf.gz -Oz -o /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz && \
             bcftools index /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz"
  n=$(docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 bash -c "bcftools view /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz | grep -cv '^#'")
  echo "$OTHER: cohort-lifted variants = $n"
done

echo "=== final summary ==="
for OTHER in $OTHERS; do
  n=$(docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 bash -c "bcftools view /v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz | grep -cv '^#'" 2>/dev/null)
  u=$(cat $V3/projection/A1_wfmash/lifted/${OTHER}/*.unmap 2>/dev/null | grep -cv '^#')
  printf "  %-22s  lifted=%d  unmapped=%d\n" "$OTHER" "$n" "$u"
done
