#!/bin/bash
# Phase C.3: TOGA2/CESAR2 fallback for the 6 remaining queries.
# PvW1 already done as smoke test. Run sequentially to avoid CPU contention.
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
IMG=avianalter/toga2:latest

# strain → genbank accession
declare -A ACC
ACC[Sal-I]=GCA_000002415.2
ACC[PvT01]=GCA_900093545.1
ACC[PvC01]=GCA_900093535.1
ACC[PAM]=GCA_949152365.1
ACC[PvSY56]=GCA_003402215.1
ACC[MHC087]=GCA_040114635.1

T_ACC=GCA_900093555.2   # PvP01 anchor

run_toga2() {
  local Q=$1
  local Q_ACC=${ACC[$Q]}
  local outdir=$V3/work/02c_toga/PvP01-as-ref/${Q}
  local chain=$V3/work/01_chains/${T_ACC}.${Q_ACC}.cleaned.chain
  local prep=$V3/work/02c_toga/PvP01_prep_${Q}
  local flagged_bed=$V3/work/02c_toga/PvP01-as-ref-${Q}_needs_cesar2.bed12

  # Check prerequisites
  [ ! -s $chain ]       && { echo "[SKIP $Q] no chain"; return; }
  [ ! -s $flagged_bed ] && { echo "[SKIP $Q] no flagged BED12"; return; }

  # Step 1: prepare-input per query (filters bed12 to only needs_cesar2 genes)
  if [ ! -d $prep ]; then
    echo "[$Q] prepare-input..."
    docker run --rm -v $V3:/v3 -w /v3 $IMG bash -c "
      source /opt/TOGA2/toga2/bin/activate
      python3 /opt/TOGA2/toga2.py prepare-input \
        --ref_2bit /v3/projection/A2_kegalign/2bit/${T_ACC}.2bit \
        --ref_annot ${flagged_bed/$V3/\/v3} \
        --ref_isoforms /v3/inputs/annotations/PvP01.isoforms.tsv \
        --disable_intron_classification \
        --disable_cesar_profiles \
        --output ${prep/$V3/\/v3} 2>&1 | tail -5"
  fi

  # Skip if already done
  if [ -s $outdir/query_annotation.bed ]; then
    echo "[skip $Q] already done ($(wc -l < $outdir/query_annotation.bed) genes)"
    return
  fi

  # Clean and run
  docker run --rm -v $V3:/v3 $IMG rm -rf ${outdir/$V3/\/v3} 2>/dev/null
  mkdir -p $outdir
  echo "[$Q] running TOGA2 ($(wc -l < $prep/toga.transcripts.bed) transcripts)..."
  /usr/bin/time -v docker run --rm -v $V3:/v3 -w /v3 $IMG bash -c "
    source /opt/TOGA2/toga2/bin/activate
    python3 /opt/TOGA2/toga2.py run \
      --ref_2bit /v3/projection/A2_kegalign/2bit/${T_ACC}.2bit \
      --query_2bit /v3/projection/A2_kegalign/2bit/${Q_ACC}.2bit \
      --chain_file ${chain/$V3/\/v3} \
      --ref_annotation ${prep/$V3/\/v3}/toga.transcripts.bed \
      --isoform_file ${prep/$V3/\/v3}/toga.isoforms.tsv \
      --no_u12_file --no_spliceai --no_utr_annotation \
      --ignore_crashed_parallel_batches \
      --output ${outdir/$V3/\/v3} 2>&1 | tail -10" 2> $V3/work/logs/toga2_${Q}.time.log

  if [ -s $outdir/query_annotation.bed ]; then
    ng=$(wc -l < $outdir/query_annotation.bed)
    nint=$(awk -F'\t' '$3=="I"' $outdir/loss_summary.tsv | wc -l)
    nloss=$(awk -F'\t' '$3=="L"' $outdir/loss_summary.tsv | wc -l)
    t=$(grep Elapsed $V3/work/logs/toga2_${Q}.time.log | awk '{print $NF}')
    echo "[done $Q] $ng projections (Intact=$nint, Lost=$nloss) wall=$t"
  else
    echo "[FAIL $Q]"
  fi
}

for Q in Sal-I PvT01 PvC01 PAM PvSY56 MHC087; do
  run_toga2 $Q
done

echo
echo "===== Phase C.3 summary ====="
for Q in PvW1 Sal-I PvT01 PvC01 PAM PvSY56 MHC087; do
  outdir=$V3/work/02c_toga/PvP01-as-ref/${Q}
  if [ -s $outdir/query_annotation.bed ]; then
    ng=$(wc -l < $outdir/query_annotation.bed)
    n_one2one=$(awk -F'\t' '$5=="one2one"' $outdir/orthology_classification.tsv 2>/dev/null | wc -l)
    n_lost=$(awk -F'\t' '$5=="one2zero"' $outdir/orthology_classification.tsv 2>/dev/null | wc -l)
    nint=$(awk -F'\t' '$3=="I"' $outdir/loss_summary.tsv 2>/dev/null | wc -l)
    nloss=$(awk -F'\t' '$3=="L"' $outdir/loss_summary.tsv 2>/dev/null | wc -l)
    printf "  %-8s genes=%5d  intact=%5d  lost=%5d  1to1=%5d  1to0=%5d\n" "$Q" "$ng" "$nint" "$nloss" "$n_one2one" "$n_lost"
  else
    printf "  %-8s NO OUTPUT\n" "$Q"
  fi
done

touch $V3/work/logs/checkpoints/02c_toga.done
