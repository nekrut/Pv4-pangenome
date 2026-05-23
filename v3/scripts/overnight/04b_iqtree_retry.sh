#!/bin/bash
# Re-run IQ-TREE for the 403 genes that failed with "less than 4 sequences" due
# to identical-sequence collapse. Drop -B 1000 (no bootstrap on ultra-conserved
# trees) so iqtree at least produces a treefile.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p logs/overnight
log() { echo "$(date +%H:%M:%S) [iqtree_retry] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

MISSING=$(comm -23 \
  <(ls work/06_msa/core_v3_clean/*.cleaned.fa | xargs -n1 basename | sed 's/.codon.cleaned.fa//' | sort) \
  <(find work/06_msa/core_v3_trees -name '*.treefile' | xargs -n1 basename | sed 's/.treefile//' | sort))
N=$(echo "$MISSING" | wc -l)
log "$N missing trees"

run_one() {
  local g=$1
  local outdir=work/06_msa/core_v3_trees/$g
  [ -s $outdir/$g.treefile ] && return 0
  mkdir -p $outdir
  cp work/06_msa/core_v3_clean/${g}.codon.cleaned.fa $outdir/$g.aln.fa
  # No -B 1000 since these have <4 unique sequences
  docker run --rm -v $V3:/v3 -w /v3/$outdir \
    quay.io/biocontainers/iqtree:3.1.2--h8471819_0 \
    iqtree -s $g.aln.fa -m MFP -T 2 --quiet --prefix $g \
    2>>$V3/logs/overnight/iqtree_retry.log
  if [ ! -s $outdir/$g.treefile ]; then
    echo "$g: iqtree retry FAILED" >> $V3/logs/overnight/STATUS.md
  fi
}
export -f run_one
export V3

echo "$MISSING" | xargs -P 8 -I {} bash -c 'run_one "$@"' _ {}
DONE=$(find work/06_msa/core_v3_trees -name '*.treefile' 2>/dev/null | wc -l)
log "done; total trees: $DONE"
