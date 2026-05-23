#!/bin/bash
# Relaxed MSA pipeline: min_intact=5 instead of 7. Captures genes that are
# intact in 5+ of 7 query strains (allows up to 2 missing strains). Should
# recover most drug-resistance + vaccine target genes lost from the strict run.
#
# Stages: MSA build (8 shards) → trimAl → IQ-TREE → HyPhy BUSTED
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p logs/overnight work/06_msa/core_relaxed work/06_msa/core_relaxed_clean work/06_msa/core_relaxed_trees work/06_msa/core_relaxed_hyphy
log() { echo "$(date +%H:%M:%S) [relaxed] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting min_intact=5 pipeline"

# Stage 1: 8-shard MSA build
log "stage 1: MSA build (8 parallel shards, min_intact=5)"
for s in $(seq 0 7); do
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
    --min-intact 5 \
    --out-dir work/06_msa/core_relaxed \
    --threads-per-gene 1 \
    --shard $s --num-shards 8 \
    > logs/overnight/relaxed_shard${s}.log 2>&1 &
done
wait
log "stage 1 done: $(ls work/06_msa/core_relaxed/*.codon.aln.fa | wc -l) MSAs"

# Stage 2: trimAl
log "stage 2: trimAl"
run_trim() {
  local f=$1
  local g=$(basename $f .codon.aln.fa)
  local out=work/06_msa/core_relaxed_clean/${g}.codon.cleaned.fa
  [ -s $out ] && return 0
  docker run --rm -v $V3:/v3 -w /v3 \
    quay.io/biocontainers/trimal:1.4.1--h9f5acd7_7 \
    trimal -in $f -out $out -automated1 2>/dev/null
}
export -f run_trim
export V3
ls work/06_msa/core_relaxed/*.codon.aln.fa | xargs -P 12 -I {} bash -c 'run_trim "$@"' _ {}
log "stage 2 done: $(ls work/06_msa/core_relaxed_clean/*.cleaned.fa | wc -l) cleaned"

# Stage 3: IQ-TREE (no -B since some genes may have <4 unique seqs)
log "stage 3: IQ-TREE"
run_tree() {
  local f=$1
  local g=$(basename $f .codon.cleaned.fa)
  local outdir=work/06_msa/core_relaxed_trees/$g
  [ -s $outdir/$g.treefile ] && return 0
  mkdir -p $outdir
  cp $f $outdir/$g.aln.fa
  # Try with bootstrap first; if <4 unique seqs, fall back without -B
  docker run --rm -v $V3:/v3 -w /v3/$outdir \
    quay.io/biocontainers/iqtree:3.1.2--h8471819_0 \
    iqtree -s $g.aln.fa -m MFP -B 1000 -T 2 --quiet --prefix $g \
    2>$V3/logs/overnight/relaxed_iqtree_${g}.log
  if [ ! -s $outdir/$g.treefile ]; then
    docker run --rm -v $V3:/v3 -w /v3/$outdir \
      quay.io/biocontainers/iqtree:3.1.2--h8471819_0 \
      iqtree -s $g.aln.fa -m MFP -T 2 --quiet --prefix $g -redo \
      2>>$V3/logs/overnight/relaxed_iqtree_${g}.log
  fi
}
export -f run_tree
ls work/06_msa/core_relaxed_clean/*.cleaned.fa | xargs -P 8 -I {} bash -c 'run_tree "$@"' _ {}
log "stage 3 done: $(find work/06_msa/core_relaxed_trees -name '*.treefile' | wc -l) trees"

# Stage 4: HyPhy BUSTED (bulk)
log "stage 4: HyPhy BUSTED bulk"
run_busted() {
  local g=$1
  local outdir=work/06_msa/core_relaxed_hyphy/$g
  [ -s $outdir/busted.json ] && return 0
  mkdir -p $outdir
  cp work/06_msa/core_relaxed_clean/${g}.codon.cleaned.fa $outdir/aln.fa
  cp work/06_msa/core_relaxed_trees/${g}/${g}.treefile $outdir/tree.nwk 2>/dev/null
  [ ! -s $outdir/tree.nwk ] && return 1
  docker run --rm -v $V3:/v3 -w /v3/$outdir \
    quay.io/biocontainers/hyphy:2.5.99--h74d3ee0_0 \
    hyphy busted --alignment aln.fa --tree tree.nwk --output busted.json \
    2>>$V3/logs/overnight/relaxed_busted_${g}.log
}
export -f run_busted
find work/06_msa/core_relaxed_trees -name '*.treefile' | xargs -n1 basename | sed 's/.treefile//' | \
  xargs -P 6 -I {} bash -c 'run_busted "$@"' _ {}
log "stage 4 done: $(find work/06_msa/core_relaxed_hyphy -name busted.json | wc -l) BUSTED jsons"

log "ALL RELAXED STAGES COMPLETE"
