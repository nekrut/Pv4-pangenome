#!/bin/bash
# All-vs-all KegAlign on GPU.
# Workflow per pair:
#  (1) Pre-stage aln_ref.2bit / aln_query.2bit symlinks in workdir
#  (2) kegalign FASTA → emits lastz commands
#  (3) xargs -P 16 runs lastz commands → tmp*.axt
#  (4) cat tmp*.axt → final axt
# 56 directional pairs on single A5000 (sequential).
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
ACCS="GCA_000002415.2 GCA_900093555.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"
IMG=quay.io/biocontainers/kegalign-full:0.1.2.9--hdfd78af_0

mkdir -p $V3/projection/A2_kegalign/{axt,chain,lifted,work}

n_pair=0
for TGT in $ACCS; do
  for QRY in $ACCS; do
    if [ "$TGT" = "$QRY" ]; then continue; fi
    n_pair=$((n_pair+1))
    PAIR=${TGT}__vs__${QRY}
    OUT_AXT=$V3/projection/A2_kegalign/axt/${PAIR}.axt
    WORK=$V3/projection/A2_kegalign/work/${PAIR}
    if [ -s $OUT_AXT ]; then
      echo "[skip $n_pair/56] $PAIR done ($(du -h $OUT_AXT|cut -f1))"
      continue
    fi
    echo "===== [$n_pair/56] KegAlign: $TGT → $QRY ====="
    rm -rf $WORK; mkdir -p $WORK
    # Pre-stage 2bit files with prefix names kegalign's emitted commands expect
    cp $V3/projection/A2_kegalign/2bit/${TGT}.2bit $WORK/aln_ref.2bit
    cp $V3/projection/A2_kegalign/2bit/${QRY}.2bit $WORK/aln_query.2bit
    /usr/bin/time -v docker run --rm --gpus '"device=0"' \
      -v $V3:/v3 -w /v3/projection/A2_kegalign/work/${PAIR} \
      $IMG bash -c "
        kegalign /v3/genomes/softmasked/${TGT}.fa /v3/genomes/softmasked/${QRY}.fa aln_ \
                 --strand both --format axt --num_gpu 1 --num_threads 16 \
                 > commands.txt 2> kegalign.log
        n_cmds=\$(wc -l < commands.txt)
        echo \"emitted \$n_cmds lastz commands\"
        xargs -P 16 -I {} bash -c '{}' < commands.txt
        cat tmp*.axt 2>/dev/null > /v3/projection/A2_kegalign/axt/${PAIR}.axt
      " 2> $V3/logs/A2_kegalign_${PAIR}.time.log
    grep -E 'Elapsed' $V3/logs/A2_kegalign_${PAIR}.time.log | head -1
    if [ -s $OUT_AXT ]; then
      sz=$(du -h $OUT_AXT|cut -f1); n_aln=$(grep -c "^[0-9]" $OUT_AXT 2>/dev/null || echo 0)
      echo "[done] $PAIR  size=$sz  alignments=$n_aln"
    else
      echo "[FAIL] $PAIR"
    fi
    rm -rf $WORK
  done
done

echo
echo "===== axt summary ====="
ls $V3/projection/A2_kegalign/axt/ | wc -l; echo "axt files / 56 expected"
du -sh $V3/projection/A2_kegalign/axt/

echo
echo "===== axt → chain (axtChain + chainSort) ====="
for TGT in $ACCS; do
  for QRY in $ACCS; do
    [ "$TGT" = "$QRY" ] && continue
    AXT=$V3/projection/A2_kegalign/axt/${TGT}__vs__${QRY}.axt
    CHAIN=$V3/projection/A2_kegalign/chain/${TGT}_to_${QRY}.chain
    [ -s $CHAIN ] && continue
    [ ! -s $AXT ] && continue
    $V3/tools/axtChain -linearGap=loose $AXT \
        $V3/projection/A2_kegalign/2bit/${TGT}.2bit \
        $V3/projection/A2_kegalign/2bit/${QRY}.2bit \
        /tmp/$$_${TGT}_${QRY}.unsorted.chain 2>/dev/null
    $V3/tools/chainSort /tmp/$$_${TGT}_${QRY}.unsorted.chain > $CHAIN
    rm -f /tmp/$$_${TGT}_${QRY}.unsorted.chain
    n=$(grep -c '^chain ' $CHAIN)
    echo "$TGT → $QRY  chains=$n"
  done
done

echo
echo "===== A2 CrossMap liftover (PvP01 → 7 others) ====="
TGT=GCA_900093555.2
mkdir -p $V3/projection/A2_kegalign/lifted
for QRY in GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1; do
  CHAIN=$V3/projection/A2_kegalign/chain/${TGT}_to_${QRY}.chain
  [ ! -s $CHAIN ] && { echo "skip $QRY (no chain)"; continue; }
  outdir=$V3/projection/A2_kegalign/lifted/${QRY}; mkdir -p $outdir
  for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 API MIT; do
    out=$outdir/Pv4_${chr}_on_${QRY}.vcf
    [ -s ${out}.gz ] && continue
    docker run --rm -v $V3:/v3 -w /v3 quay.io/biocontainers/crossmap:0.7.3--pyhdfd78af_0 \
      CrossMap vcf projection/A2_kegalign/chain/${TGT}_to_${QRY}.chain \
                   projection/A1_wfmash/mg_renamed/Pv4_${chr}.vcf.gz \
                   genomes/softmasked/${QRY}.fa \
                   projection/A2_kegalign/lifted/${QRY}/Pv4_${chr}_on_${QRY}.vcf 2>/dev/null
  done
  docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "cd /v3/projection/A2_kegalign/lifted/${QRY} && \
             for f in Pv4_*_on_${QRY}.vcf; do [ -f \$f ] && bgzip -f \$f && tabix -p vcf \${f}.gz; done && \
             bcftools concat -a Pv4_*_on_${QRY}.vcf.gz -Oz -o /v3/projection/A2_kegalign/Pv4_cohort_on_${QRY}.vcf.gz 2>/dev/null && \
             bcftools index /v3/projection/A2_kegalign/Pv4_cohort_on_${QRY}.vcf.gz"
  echo "concat $QRY"
done

echo
echo "===== A2 FINAL ====="
for QRY in GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1; do
  n=$(docker run --rm -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
        bash -c "bcftools view /v3/projection/A2_kegalign/Pv4_cohort_on_${QRY}.vcf.gz 2>/dev/null | grep -cv '^#'" 2>/dev/null)
  printf "  PvP01 → %-22s  lifted=%s\n" "$QRY" "${n:-NA}"
done
echo "Total axt preserved: $(ls $V3/projection/A2_kegalign/axt/ | wc -l)/56"
