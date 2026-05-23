#!/bin/bash
# Path A2 VCF liftover: lift MalariaGEN Pv4 cohort VCFs from PvP01 coords onto
# each of 7 other reference assemblies using KegAlign-derived chains.
set -eu
V3=/media/anton/data/sandbox/Pv4/v3
PVP01=GCA_900093555.2
OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"
ALL_CHRS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 API MIT"
PARALLEL_J=${1:-8}  # parallel CrossMap jobs

mkdir -p $V3/projection/A2_lastz/chain $V3/projection/A2_lastz/lifted $V3/logs

# Stage KegAlign chains under the expected A2 naming
echo "===== Stage chains ====="
for O in $OTHERS; do
  src=$V3/work/01_chains/${PVP01}.${O}.cleaned.chain
  dst=$V3/projection/A2_lastz/chain/PvP01_to_${O}.chain
  if [ ! -s $dst ]; then ln -sf $src $dst; fi
  echo "  $O: $(wc -l < $src) lines"
done

echo
echo "===== Step 1: CrossMap liftover, $PARALLEL_J parallel ====="
crossmap_one() {
  local O=$1 chr=$2
  local outdir=$V3/projection/A2_lastz/lifted/$O
  mkdir -p $outdir
  local out=$outdir/Pv4_${chr}_on_${O}.vcf
  [ -s ${out}.gz ] && { echo "$O $chr already done"; return 0; }
  docker run --rm -v $V3:/v3 -w /v3 \
    quay.io/biocontainers/crossmap:0.7.3--pyhdfd78af_0 \
    CrossMap vcf projection/A2_lastz/chain/PvP01_to_${O}.chain \
                 projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz \
                 genomes/softmasked/${O}.fa \
                 projection/A2_lastz/lifted/${O}/Pv4_${chr}_on_${O}.vcf \
                 > $V3/logs/a2_${O}_${chr}.log 2>&1
  local rc=$?
  if [ $rc -ne 0 ] || [ ! -s $out ]; then
    echo "$O $chr FAILED (rc=$rc; output empty); see logs/a2_${O}_${chr}.log" >&2
    return 1
  fi
  echo "$O $chr done ($(wc -l < $out) records)"
}
export -f crossmap_one
export V3
JOBS=$(for o in $OTHERS; do for c in $ALL_CHRS; do echo "$o $c"; done; done)
echo "$JOBS" | xargs -P $PARALLEL_J -n 2 bash -c 'crossmap_one "$@"' _

echo
echo "===== Step 2: per-chr sort + bgzip + tabix ====="
sort_one() {
  local O=$1 chr=$2
  local f=$V3/projection/A2_lastz/lifted/$O/Pv4_${chr}_on_${O}.vcf
  [ -s ${f}.gz ] && return 0
  [ ! -s $f ] && return 0
  docker run --rm -v $V3:/v3 \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "cd /v3/projection/A2_lastz/lifted/$O && \
             bcftools sort -Oz -o Pv4_${chr}_on_${O}.sorted.vcf.gz Pv4_${chr}_on_${O}.vcf && \
             mv Pv4_${chr}_on_${O}.sorted.vcf.gz Pv4_${chr}_on_${O}.vcf.gz && \
             tabix -p vcf Pv4_${chr}_on_${O}.vcf.gz"
  rm -f $f
}
export -f sort_one
echo "$JOBS" | xargs -P $PARALLEL_J -n 2 bash -c 'sort_one "$@"' _

echo
echo "===== Step 3: concat per target ====="
concat_one() {
  local O=$1
  local dir=$V3/projection/A2_lastz/lifted/$O
  local out=$V3/projection/A2_lastz/Pv4_cohort_on_${O}.vcf.gz
  [ -s $out ] && { echo "$O already concatenated"; return 0; }
  # Build explicit file list — skip empty files (e.g., chains without API coverage)
  local files=""
  for f in $dir/Pv4_*_on_${O}.vcf.gz; do
    [ ! -s "$f" ] && continue
    local n=$(docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
              bcftools view -H /v3/projection/A2_lastz/lifted/$O/$(basename $f) 2>/dev/null | head -1 | wc -l)
    [ "$n" -eq 0 ] && { echo "  skip empty $(basename $f) for $O"; continue; }
    files="$files /v3/projection/A2_lastz/lifted/$O/$(basename $f)"
  done
  if [ -z "$files" ]; then
    echo "$O concat: no non-empty inputs"; return 1
  fi
  docker run --rm -v $V3:/v3 \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "bcftools concat -a $files -Oz -o /v3/projection/A2_lastz/Pv4_cohort_on_${O}.vcf.gz && \
             bcftools index /v3/projection/A2_lastz/Pv4_cohort_on_${O}.vcf.gz"
  local rc=$?
  if [ $rc -ne 0 ]; then echo "$O concat FAILED rc=$rc"; return 1; fi
  echo "concatenated $O"
}
export -f concat_one
echo "$OTHERS" | tr ' ' '\n' | xargs -P 7 -I {} bash -c 'concat_one "$@"' _ {}

echo
echo "===== Done ====="
ls -lh $V3/projection/A2_lastz/Pv4_cohort_on_*.vcf.gz
