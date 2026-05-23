#!/bin/bash
# Launch N sharded Phase H jobs in parallel.
set -eu
cd /media/anton/data/sandbox/Pv4/v3
NUM=${1:-8}
mkdir -p logs
for s in $(seq 0 $((NUM-1))); do
  nohup python3 scripts/build_8way_msa_v2.py \
    --queries Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087 \
    --ref-gff inputs/annotations/PvP01.genbank.gff3 \
    --ref-fasta inputs/assemblies/PvP01.fa \
    --query-gff work/02d_merged/PvP01-as-ref/Sal-I.annotation.gff3 \
                work/02d_merged/PvP01-as-ref/PvW1.annotation.gff3 \
                work/02d_merged/PvP01-as-ref/PAM.annotation.gff3 \
                work/02d_merged/PvP01-as-ref/PvSY56.annotation.gff3 \
                work/02d_merged/PvP01-as-ref/PvT01.annotation.gff3 \
                work/02d_merged/PvP01-as-ref/PvC01.annotation.gff3 \
                work/02d_merged/PvP01-as-ref/MHC087.annotation.gff3 \
    --query-fasta inputs/assemblies/Sal-I.fa inputs/assemblies/PvW1.fa \
                  inputs/assemblies/PAM.fa inputs/assemblies/PvSY56.fa \
                  inputs/assemblies/PvT01.fa inputs/assemblies/PvC01.fa \
                  inputs/assemblies/MHC087.fa \
    --merged-base work/02d_merged/PvP01-as-ref \
    --min-intact 7 \
    --out-dir work/06_msa/core_v3 \
    --threads-per-gene 2 \
    --shard $s --num-shards $NUM \
    > logs/phase_h_shard${s}.log 2>&1 &
  echo "launched shard $s/$NUM pid $!"
done
wait
echo "all shards complete"
