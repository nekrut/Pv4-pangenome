#!/bin/bash
# Phase D: repeat C.1-C.4 with reciprocal anchors PvW1, PAM, PvSY56.
# Each anchor → 7 queries (including PvP01, the other 6 non-anchor strains).
#
# Per anchor A:
#   1. Build A.bed12 + A.isoforms.tsv + A.family_list.tsv (if PvP01-anchor missing)
#   2. Liftoff A → Q for 7 queries  → work/02a_liftoff/${A}-as-ref/${Q}/
#   3. Triage A → Q (Phase C.2)     → work/02b_triage/${A}-as-ref/${Q}/
#   4. TOGA2 A → Q (Phase C.3)      → work/02c_toga/${A}-as-ref/${Q}/
#   5. Merge A→Q outputs (Phase C.4)→ work/02d_merged/${A}-as-ref/
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
LIFTOFF_IMG=quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0
TOGA_IMG=avianalter/toga2:latest
PYFAIDX_IMG=quay.io/biocontainers/pyfaidx:0.8.1.1--pyhdfd78af_0
GFFREAD_IMG=quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4

# Anchors and their GenBank accessions
declare -A ACC
ACC[Sal-I]=GCA_000002415.2
ACC[PvP01]=GCA_900093555.2
ACC[PvT01]=GCA_900093545.1
ACC[PvC01]=GCA_900093535.1
ACC[PvW1]=GCA_914969965.1
ACC[PAM]=GCA_949152365.1
ACC[PvSY56]=GCA_003402215.1
ACC[MHC087]=GCA_040114635.1

ANCHORS="PvW1 PAM PvSY56"
ALL_STRAINS="PvP01 PvW1 PAM PvSY56 Sal-I PvT01 PvC01 MHC087"

mkdir -p $V3/work/{02a_liftoff,02b_triage,02c_toga,02d_merged}

# Build per-anchor BED12 and isoforms (gffread → BED12)
build_anchor_inputs() {
  local A=$1
  local agff=$V3/inputs/annotations/plasmodb-68/${A}.gff3
  local afa=$V3/inputs/assemblies/${A}.fa
  local bed12=$V3/inputs/annotations/${A}.bed12
  local iso=$V3/inputs/annotations/${A}.isoforms.tsv
  if [ ! -s $bed12 ]; then
    echo "[$A] building bed12 + isoforms"
    docker run --rm -v $V3:/v3 $GFFREAD_IMG \
      gffread --bed $agff -o $V3/inputs/annotations/${A}.bed12_raw 2>/dev/null
    # filter to protein-coding only, strip col 13
    awk -F'\t' '!/^#/ && $3=="protein_coding_gene" {n=split($9,a,";"); for(i=1;i<=n;i++) if(match(a[i],/^ID=/)) {print substr(a[i],4); break}}' $agff | sort -u > /tmp/${A}_pc.txt
    python3 - <<PY
import re
pc=set(open('/tmp/${A}_pc.txt').read().split())
with open('$V3/inputs/annotations/${A}.bed12_raw') as fin, open('$V3/inputs/annotations/${A}.bed12','w') as fout:
    for ln in fin:
        f=ln.rstrip('\n').split('\t')
        if len(f)<13: continue
        m=re.search(r'geneID=([^;]+)', f[12])
        if not m: continue
        if m.group(1) in pc: fout.write('\t'.join(f[:12])+'\n')
PY
    rm -f $V3/inputs/annotations/${A}.bed12_raw
  fi
  if [ ! -s $iso ]; then
    awk -F'\t' '!/^#/ && $3=="mRNA" {
      n=split($9,a,";"); tx=""; gene=""
      for(i=1;i<=n;i++) {
        if (match(a[i],/^ID=/)) tx=substr(a[i],4)
        else if (match(a[i],/^Parent=/)) gene=substr(a[i],8)
      }
      if (tx && gene) print gene "\t" tx
    }' $agff | sort -k1,1 -k2,2 > $iso
  fi
}

# C.1 Liftoff anchor → Q
run_liftoff() {
  local A=$1 Q=$2
  local outdir=$V3/work/02a_liftoff/${A}-as-ref/${Q}
  local outgff=$outdir/${Q}.lifted.gff3
  if [ -s $outgff ]; then echo "[skip liftoff $A→$Q]"; return; fi
  mkdir -p $outdir/intermediate
  echo "[liftoff $A → $Q]"
  /usr/bin/time -v docker run --rm -v $V3:/v3 -w /v3 $LIFTOFF_IMG \
    liftoff -g /v3/inputs/annotations/plasmodb-68/${A}.gff3 \
            -f /v3/work/02a_liftoff/feature_types.txt \
            -o /v3/work/02a_liftoff/${A}-as-ref/${Q}/${Q}.lifted.gff3 \
            -u /v3/work/02a_liftoff/${A}-as-ref/${Q}/${Q}.unmapped.txt \
            -dir /v3/work/02a_liftoff/${A}-as-ref/${Q}/intermediate \
            -copies -sc 0.90 -d 5 -flank 0.1 -polish -p 4 \
            /v3/inputs/assemblies/${Q}.fa /v3/inputs/assemblies/${A}.fa \
            2> $V3/work/logs/liftoff_${A}_${Q}.log
}

# C.2 triage
run_triage() {
  local A=$1 Q=$2
  local outdir=$V3/work/02b_triage/${A}-as-ref/${Q}
  local outsum=$outdir/summary.json
  if [ -s $outsum ]; then echo "[skip triage $A→$Q]"; return; fi
  mkdir -p $outdir
  # Identity thresholds — relax for PvSY56 anchor
  local id_min=0.95
  local fam_id_min=0.85
  if [ "$Q" = "PvSY56" ] || [ "$A" = "PvSY56" ]; then
    id_min=0.85
    fam_id_min=0.75
  fi
  docker run --rm -v $V3:/v3 $PYFAIDX_IMG \
    python3 /v3/scripts/phase_c2_triage.py \
      --liftoff-gff /v3/work/02a_liftoff/${A}-as-ref/${Q}/${Q}.lifted.gff3 \
      --query-fasta /v3/inputs/assemblies/${Q}.fa \
      --reference-bed /v3/inputs/annotations/${A}.bed \
      --family-list /v3/inputs/annotations/plasmodb-68/PvP01.family_list.tsv \
      --core-identity-min $id_min \
      --family-identity-min $fam_id_min \
      --output-dir /v3/work/02b_triage/${A}-as-ref/${Q}/ \
      --query-name ${Q} 2>&1 | tail -2
}

# Build flagged BED12 from triage needs_cesar2.bed (which is BED6)
build_flagged_bed12() {
  local A=$1 Q=$2
  local flagged=$V3/work/02b_triage/${A}-as-ref/${Q}/needs_cesar2.bed
  local out=$V3/work/02c_toga/${A}-as-ref-${Q}_needs_cesar2.bed12
  local bed12=$V3/inputs/annotations/${A}.bed12
  if [ -s $out ] || [ ! -s $flagged ]; then return; fi
  python3 - <<PY
import re, os
flagged=set()
with open('$flagged') as fh:
    for ln in fh:
        f=ln.rstrip('\n').split('\t')
        if len(f)>=4: flagged.add(f[3])
os.makedirs(os.path.dirname('$out'), exist_ok=True)
# Handle two transcript-id conventions:
#   PlasmoDB PvP01:        PVP01_xxx.1  (strip '.N')
#   PlasmoDB PvW1/PAM:     PVW1_xxx_t1  (strip '_tN')
def tx_to_gene(tx):
    m = re.match(r'^(.+)_t\d+\$', tx)
    if m: return m.group(1)
    if re.match(r'^.+\.\d+\$', tx):
        return tx.rsplit('.', 1)[0]
    return tx
n_in = n_kept = 0
with open('$bed12') as fin, open('$out','w') as fout:
    for ln in fin:
        n_in += 1
        f=ln.split('\t'); tx=f[3]
        gene = tx_to_gene(tx)
        if gene in flagged:
            fout.write(ln)
            n_kept += 1
print(f"  [build_flagged_bed12 \$A→\$Q] kept {n_kept}/{n_in}")
PY
}

# C.3 TOGA2 prep-input + run
run_toga2() {
  local A=$1 Q=$2
  local A_ACC=${ACC[$A]}
  local Q_ACC=${ACC[$Q]}
  local outdir=$V3/work/02c_toga/${A}-as-ref/${Q}
  local chain=$V3/work/01_chains/${A_ACC}.${Q_ACC}.cleaned.chain
  local flagged_bed=$V3/work/02c_toga/${A}-as-ref-${Q}_needs_cesar2.bed12
  local prep=$V3/work/02c_toga/${A}_prep_${Q}
  if [ ! -s $chain ]; then echo "[SKIP $A→$Q] no chain"; return; fi
  if [ ! -s $flagged_bed ]; then echo "[SKIP $A→$Q] no flagged BED12"; return; fi
  if [ -s $outdir/query_annotation.bed ]; then echo "[skip toga $A→$Q]"; return; fi
  # prepare-input
  if [ ! -d $prep ]; then
    docker run --rm -v $V3:/v3 -w /v3 $TOGA_IMG bash -c "
      source /opt/TOGA2/toga2/bin/activate
      python3 /opt/TOGA2/toga2.py prepare-input \
        --ref_2bit /v3/projection/A2_kegalign/2bit/${A_ACC}.2bit \
        --ref_annot ${flagged_bed/$V3/\/v3} \
        --ref_isoforms /v3/inputs/annotations/${A}.isoforms.tsv \
        --disable_intron_classification \
        --disable_cesar_profiles \
        --output ${prep/$V3/\/v3} 2>&1 | tail -3"
  fi
  mkdir -p $outdir
  echo "[toga2 $A → $Q] ($(wc -l < $prep/toga.transcripts.bed) transcripts)"
  /usr/bin/time -v docker run --rm -v $V3:/v3 -w /v3 $TOGA_IMG bash -c "
    source /opt/TOGA2/toga2/bin/activate
    python3 /opt/TOGA2/toga2.py run \
      --ref_2bit /v3/projection/A2_kegalign/2bit/${A_ACC}.2bit \
      --query_2bit /v3/projection/A2_kegalign/2bit/${Q_ACC}.2bit \
      --chain_file ${chain/$V3/\/v3} \
      --ref_annotation ${prep/$V3/\/v3}/toga.transcripts.bed \
      --isoform_file ${prep/$V3/\/v3}/toga.isoforms.tsv \
      --no_u12_file --no_spliceai --no_utr_annotation \
      --ignore_crashed_parallel_batches \
      --output ${outdir/$V3/\/v3} 2>&1 | tail -5" 2> $V3/work/logs/toga2_${A}_${Q}.time.log
  if [ -s $outdir/query_annotation.bed ]; then
    ng=$(wc -l < $outdir/query_annotation.bed)
    nint=$(awk -F'\t' '$3=="I"' $outdir/loss_summary.tsv 2>/dev/null | wc -l)
    t=$(grep Elapsed $V3/work/logs/toga2_${A}_${Q}.time.log | awk '{print $NF}')
    echo "[done $A→$Q] genes=$ng intact=$nint wall=$t"
  fi
}

# C.4 merge per anchor
run_merge() {
  local A=$1 Q=$2
  local outdir=$V3/work/02d_merged/${A}-as-ref
  mkdir -p $outdir
  if [ -s $outdir/${Q}.annotation.gff3 ]; then return; fi
  python3 $V3/scripts/phase_c4_merge.py \
    --query $Q \
    --triage-dir $V3/work/02b_triage/${A}-as-ref/${Q} \
    --toga-dir $V3/work/02c_toga/${A}-as-ref/${Q} \
    --out-dir $outdir \
    --ref-bed $V3/inputs/annotations/${A}.bed12 2>&1 | tail -3
}

# === Per-anchor pipeline ===
for ANCHOR in $ANCHORS; do
  echo
  echo "==========================="
  echo "Phase D anchor: $ANCHOR"
  echo "==========================="
  build_anchor_inputs $ANCHOR

  QUERIES=""
  for s in $ALL_STRAINS; do [ "$s" != "$ANCHOR" ] && QUERIES="$QUERIES $s"; done
  echo "Queries for $ANCHOR: $QUERIES"

  # Step C.1: Liftoff (parallel 4 queries at a time)
  echo "--- C.1 Liftoff ($ANCHOR) ---"
  export -f run_liftoff; export V3 LIFTOFF_IMG
  echo $QUERIES | tr ' ' '\n' | xargs -P 3 -I {} bash -c "run_liftoff $ANCHOR {}"

  # Step C.2: Triage (parallel)
  echo "--- C.2 Triage ($ANCHOR) ---"
  export -f run_triage; export PYFAIDX_IMG
  echo $QUERIES | tr ' ' '\n' | xargs -P 4 -I {} bash -c "run_triage $ANCHOR {}"

  # Build flagged BED12 per query
  for Q in $QUERIES; do build_flagged_bed12 $ANCHOR $Q; done

  # Step C.3: TOGA2 (sequential — heavy)
  echo "--- C.3 TOGA2 ($ANCHOR) ---"
  for Q in $QUERIES; do run_toga2 $ANCHOR $Q; done

  # Step C.4: Merge
  echo "--- C.4 Merge ($ANCHOR) ---"
  for Q in $QUERIES; do run_merge $ANCHOR $Q & done; wait
done

touch $V3/work/logs/checkpoints/02d_merged_all_anchors.done
echo "Phase D complete."
