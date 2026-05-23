#!/bin/bash
# Phase J pt2 — HyPhy full bundle (BUSTED + aBSREL + MEME) on priority genes.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/06_msa/core_v3_hyphy/priority logs/overnight
log() { echo "$(date +%H:%M:%S) [hyphy_priority] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

# Resolve priority gene list: validator output > researcher draft > fallback
PRIORITY_TSV=""
for p in work/05_priorities/gene_priorities.tsv \
         writeup/gene_priorities.tsv \
         writeup/gene_priorities.draft.tsv \
         work/05_priorities/fallback_priorities.tsv; do
  if [ -s $p ]; then
    PRIORITY_TSV=$p; break
  fi
done
if [ -z "$PRIORITY_TSV" ]; then
  log "no priority list found; aborting"; exit 1
fi
log "using priority list: $PRIORITY_TSV"

# Extract distinct PVP01_xxx IDs (grep — default awk lacks gawk's 3-arg match)
PIDS=$(grep -oE 'PVP01_[0-9]+' $PRIORITY_TSV | sort -u)
N=$(echo "$PIDS" | wc -l)
log "$N priority PVP01 IDs"

run_one() {
  local pid=$1
  local aln=work/06_msa/core_v3_clean/${pid}.codon.cleaned.fa
  local tree=work/06_msa/core_v3_trees/${pid}/${pid}.treefile
  local outdir=work/06_msa/core_v3_hyphy/priority/${pid}
  mkdir -p $outdir
  [ ! -s $aln ] && { echo "$pid: no aln"; return; }
  [ ! -s $tree ] && { echo "$pid: no tree"; return; }
  cp $aln  $outdir/aln.fa
  cp $tree $outdir/tree.nwk
  for METHOD in busted absrel meme fel; do
    [ -s $outdir/${METHOD}.json ] && continue
    docker run --rm -v $V3:/v3 -w /v3/$outdir \
      quay.io/biocontainers/hyphy:2.5.99--h74d3ee0_0 \
      hyphy $METHOD --alignment aln.fa --tree tree.nwk \
            --output ${METHOD}.json 2>>$V3/logs/overnight/hyphy_priority_${pid}.log \
      || echo "$pid $METHOD FAILED" >> $V3/logs/overnight/STATUS.md
  done
}
export -f run_one
export V3

echo "$PIDS" | xargs -P 4 -I {} bash -c 'run_one "$@"' _ {}
DONE=$(find work/06_msa/core_v3_hyphy/priority -name 'busted.json' 2>/dev/null | wc -l)
log "done; BUSTED jsons: $DONE"
