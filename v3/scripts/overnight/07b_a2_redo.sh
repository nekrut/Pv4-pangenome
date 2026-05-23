#!/bin/bash
# Re-run A2 VCF liftover with a single-pipeline CrossMap → bcftools sort -Oz
# to avoid 1 TB of intermediate uncompressed VCFs. Use /media/anton/scratch
# for everything.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
SCRATCH=/media/anton/scratch/A2_lastz_v2
SORT_TMP=/media/anton/scratch/A2_sort_tmp
mkdir -p $SCRATCH/lifted $SORT_TMP logs/overnight
log() { echo "$(date +%H:%M:%S) [a2_redo] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting; output to $SCRATCH"

OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"
ALL_CHRS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 API MIT"
PARALLEL_J=${1:-6}

# Stage chains (already done by previous A2 run — verify)
for O in $OTHERS; do
  src=$V3/work/01_chains/GCA_900093555.2.${O}.cleaned.chain
  dst=$V3/projection/A2_lastz/chain/PvP01_to_${O}.chain
  [ -s $dst ] || cp $src $dst
done

run_one() {
  local O=$1 chr=$2
  local outdir=$SCRATCH/lifted/$O
  mkdir -p $outdir
  local out=$outdir/Pv4_${chr}_on_${O}.vcf.gz
  [ -s ${out} ] && { echo "$O $chr already done"; return 0; }
  # Step 1: CrossMap → tmp .vcf (unsorted)
  local raw=$outdir/Pv4_${chr}_on_${O}.raw.vcf
  docker run --rm -v $V3:/v3 -v $SCRATCH:$SCRATCH -w /v3 \
    quay.io/biocontainers/crossmap:0.7.3--pyhdfd78af_0 \
    CrossMap vcf projection/A2_lastz/chain/PvP01_to_${O}.chain \
                 projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz \
                 genomes/softmasked/${O}.fa \
                 $raw \
                 >> /tmp/a2_redo_${O}_${chr}.log 2>&1
  local rc=$?
  if [ $rc -ne 0 ] || [ ! -s $raw ]; then
    echo "$O $chr FAILED at CrossMap (rc=$rc)"
    return 1
  fi
  # Step 2: bcftools sort -Oz with temp on scratch
  docker run --rm -v $SCRATCH:$SCRATCH -v $SORT_TMP:$SORT_TMP \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bcftools sort -T $SORT_TMP -Oz -o $out $raw \
    >> /tmp/a2_redo_${O}_${chr}.log 2>&1
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s $out ]; then
    echo "$O $chr FAILED at sort (rc=$rc)"
    return 1
  fi
  # Step 3: tabix
  docker run --rm -v $SCRATCH:$SCRATCH \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    tabix -p vcf $out \
    >> /tmp/a2_redo_${O}_${chr}.log 2>&1
  # cleanup raw
  rm -f $raw $out.unmap
  local n=$(docker run --rm -v $SCRATCH:$SCRATCH quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
            bcftools view -H $out 2>/dev/null | wc -l)
  echo "$O $chr done ($n records)"
}
export -f run_one
export V3 SCRATCH SORT_TMP

log "phase 1: per-chr CrossMap+sort, $PARALLEL_J parallel"
JOBS=$(for o in $OTHERS; do for c in $ALL_CHRS; do echo "$o $c"; done; done)
echo "$JOBS" | xargs -P $PARALLEL_J -n 2 bash -c 'run_one "$@"' _

log "phase 2: concat per target"
for O in $OTHERS; do
  out=$V3/projection/A2_lastz/Pv4_cohort_on_${O}.vcf.gz
  [ -s $out ] && { log "$O cohort already concatenated"; continue; }
  # Build file list of non-empty per-chr outputs
  FILES=""
  for chr in $ALL_CHRS; do
    f=$SCRATCH/lifted/$O/Pv4_${chr}_on_${O}.vcf.gz
    [ -s $f ] && FILES="$FILES $f"
  done
  if [ -z "$FILES" ]; then log "$O: no inputs"; continue; fi
  docker run --rm -v $V3:/v3 -v $SCRATCH:$SCRATCH \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "bcftools concat -a $FILES -Oz -o $out && bcftools index $out" \
    >> logs/overnight/a2_redo_concat.log 2>&1
  if [ -s $out ]; then
    log "$O cohort concatenated"
  else
    log "$O cohort CONCAT FAILED"
  fi
done

log "DONE"
