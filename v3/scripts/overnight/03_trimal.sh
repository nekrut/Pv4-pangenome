#!/bin/bash
# Phase I — trimAl on all codon MSAs from Phase H.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/06_msa/core_v3_clean logs/overnight
log() { echo "$(date +%H:%M:%S) [trimal] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting; input dir: work/06_msa/core_v3"

TOTAL=$(ls work/06_msa/core_v3/*.codon.aln.fa 2>/dev/null | wc -l)
log "$TOTAL input codon alignments"

run_one() {
  local f=$1
  local g=$(basename $f .codon.aln.fa)
  local out=work/06_msa/core_v3_clean/${g}.codon.cleaned.fa
  [ -s $out ] && return 0
  docker run --rm -v $V3:/v3 -w /v3 \
    quay.io/biocontainers/trimal:1.4.1--h9f5acd7_7 \
    trimal -in $f -out $out -automated1 2>/dev/null
}
export -f run_one
export V3

ls work/06_msa/core_v3/*.codon.aln.fa | xargs -P 16 -I {} bash -c 'run_one "$@"' _ {}
log "done; cleaned outputs: $(ls work/06_msa/core_v3_clean/*.cleaned.fa | wc -l)"
