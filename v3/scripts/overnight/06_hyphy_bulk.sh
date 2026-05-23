#!/bin/bash
# HyPhy BUSTED-only sweep on the non-priority CORE-1:1 alignments.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/06_msa/core_v3_hyphy/bulk logs/overnight
log() { echo "$(date +%H:%M:%S) [hyphy_bulk] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

# All PVP01 IDs with both cleaned aln + tree
ALL_PIDS=$(comm -12 \
  <(ls work/06_msa/core_v3_clean/*.cleaned.fa 2>/dev/null | sed 's|.*/||;s|.codon.cleaned.fa||' | sort) \
  <(find work/06_msa/core_v3_trees -name '*.treefile' 2>/dev/null | sed 's|.*/||;s|.treefile||' | sort))

# Exclude priority list (already covered by hyphy_priority.sh)
if [ -s work/05_priorities/gene_priorities.tsv ]; then
  PRIORITY=$(awk -F'\t' 'NR>1 {for(i=1;i<=NF;i++) if(match($i,/PVP01_[0-9]+/,m)) print m[0]}' work/05_priorities/gene_priorities.tsv | sort -u)
  TO_RUN=$(comm -23 <(echo "$ALL_PIDS") <(echo "$PRIORITY"))
else
  TO_RUN="$ALL_PIDS"
fi
N=$(echo "$TO_RUN" | wc -l)
log "$N genes for bulk BUSTED"

run_one() {
  local pid=$1
  local outdir=work/06_msa/core_v3_hyphy/bulk/${pid}
  mkdir -p $outdir
  [ -s $outdir/busted.json ] && return 0
  cp work/06_msa/core_v3_clean/${pid}.codon.cleaned.fa $outdir/aln.fa
  cp work/06_msa/core_v3_trees/${pid}/${pid}.treefile $outdir/tree.nwk
  docker run --rm -v $V3:/v3 -w /v3/$outdir \
    quay.io/biocontainers/hyphy:2.5.99--h74d3ee0_0 \
    hyphy busted --alignment aln.fa --tree tree.nwk --output busted.json \
    2>>$V3/logs/overnight/hyphy_bulk_${pid}.log \
    || echo "$pid BUSTED FAILED" >> $V3/logs/overnight/STATUS.md
}
export -f run_one
export V3

echo "$TO_RUN" | xargs -P 6 -I {} bash -c 'run_one "$@"' _ {}
DONE=$(find work/06_msa/core_v3_hyphy/bulk -name busted.json | wc -l)
log "done; BUSTED jsons: $DONE"
