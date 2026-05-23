#!/bin/bash
# Phase C.1: Liftoff fast pass from PvP01 → each of 7 queries.
# All 7 jobs in parallel — each Liftoff job uses internal minimap2 (~10-20 min).
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
IMG=quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0
QUERIES="PvW1 PAM PvSY56 Sal-I PvT01 PvC01 MHC087"

mkdir -p $V3/work/02a_liftoff/PvP01-as-ref

# feature_types.txt — PlasmoDB uses protein_coding_gene etc., not just "gene"
cat > $V3/work/02a_liftoff/feature_types.txt <<'EOF'
protein_coding_gene
ncRNA_gene
pseudogene
EOF

run_liftoff() {
  local Q=$1
  local OUT=$V3/work/02a_liftoff/PvP01-as-ref/${Q}
  mkdir -p $OUT/intermediate
  local outgff=$OUT/${Q}.lifted.gff3
  [ -s $outgff ] && { echo "[skip] $Q already done"; return; }
  echo "[start] Liftoff PvP01 → $Q"
  /usr/bin/time -v docker run --rm -v $V3:/v3 -w /v3 $IMG \
    liftoff \
      -g /v3/inputs/annotations/plasmodb-68/PvP01.gff3 \
      -f /v3/work/02a_liftoff/feature_types.txt \
      -o /v3/work/02a_liftoff/PvP01-as-ref/${Q}/${Q}.lifted.gff3 \
      -u /v3/work/02a_liftoff/PvP01-as-ref/${Q}/${Q}.unmapped.txt \
      -dir /v3/work/02a_liftoff/PvP01-as-ref/${Q}/intermediate \
      -copies \
      -sc 0.90 \
      -d 5 \
      -flank 0.1 \
      -polish \
      -p 4 \
      /v3/inputs/assemblies/${Q}.fa \
      /v3/inputs/assemblies/PvP01.fa \
      2> $V3/work/logs/liftoff_${Q}.log
  ng=$(grep -cP '\tgene\t' $outgff 2>/dev/null || echo 0)
  nu=$(wc -l < $OUT/${Q}.unmapped.txt 2>/dev/null || echo 0)
  echo "[done] $Q  lifted_genes=$ng  unmapped=$nu"
}
export V3 IMG
export -f run_liftoff

echo "$QUERIES" | tr ' ' '\n' | xargs -P 4 -I {} bash -c 'run_liftoff "$@"' _ {}

echo
echo "===== Liftoff summary ====="
for Q in $QUERIES; do
  outgff=$V3/work/02a_liftoff/PvP01-as-ref/${Q}/${Q}.lifted.gff3
  unmap=$V3/work/02a_liftoff/PvP01-as-ref/${Q}/${Q}.unmapped.txt
  if [ -s $outgff ]; then
    ng=$(grep -cP '\tgene\t' $outgff)
    nu=$(wc -l < $unmap)
    printf "  %-8s  lifted=%5d  unmapped=%4d\n" "$Q" "$ng" "$nu"
  else
    printf "  %-8s  FAIL\n" "$Q"
  fi
done
touch $V3/work/logs/checkpoints/02a_liftoff.done
