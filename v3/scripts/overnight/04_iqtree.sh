#!/bin/bash
# Phase J pt1 — IQ-TREE on cleaned codon MSAs.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/06_msa/core_v3_trees logs/overnight
log() { echo "$(date +%H:%M:%S) [iqtree] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"
TOTAL=$(ls work/06_msa/core_v3_clean/*.cleaned.fa 2>/dev/null | wc -l)
log "$TOTAL cleaned alignments"

run_one() {
  local f=$1
  local g=$(basename $f .codon.cleaned.fa)
  local outdir=work/06_msa/core_v3_trees/$g
  [ -s $outdir/$g.treefile ] && return 0
  mkdir -p $outdir
  cp $f $outdir/$g.aln.fa
  docker run --rm -v $V3:/v3 -w /v3/$outdir \
    quay.io/biocontainers/iqtree:3.1.2--h8471819_0 \
    iqtree -s $g.aln.fa -m MFP -B 1000 -T 2 --quiet --prefix $g 2>>$V3/logs/overnight/iqtree.log
}
export -f run_one
export V3

# Process in parallel: 8 workers * 2 threads each = 16 cores
ls work/06_msa/core_v3_clean/*.cleaned.fa | xargs -P 8 -I {} bash -c 'run_one "$@"' _ {}
DONE=$(find work/06_msa/core_v3_trees -name '*.treefile' 2>/dev/null | wc -l)
log "done; trees: $DONE"
