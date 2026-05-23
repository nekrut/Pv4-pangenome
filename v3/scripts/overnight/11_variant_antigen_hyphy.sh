#!/bin/bash
# Step 4 of "telomeric families" plan: full HyPhy bundle on variant-antigen
# genes that survived the relaxed MSA filter (min_intact=5).
# Methods: BUSTED-MH (multi-hit), aBSREL (branch-site), MEME (per-site episodic),
# FEL (per-site pervasive).
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p logs/overnight work/06_msa/core_relaxed_hyphy_va
log() { echo "$(date +%H:%M:%S) [va_hyphy] $1" | tee -a logs/overnight/STATUS.md; }
log "starting"

# Extract variant-antigen PvP01 IDs from family_table that have relaxed MSAs
PIDS=$(python3 - <<'PY'
import csv, os, glob
msa_ids = {os.path.basename(f).replace('.codon.aln.fa', '')
           for f in glob.glob('work/06_msa/core_relaxed/*.codon.aln.fa')}
va_fams = {'PIR','PHIST','Pv-fam','MSP','DBP','EBA','RBP','AMA','RAP','SERA','TRAg','STP1','RESA'}
seen = set()
with open('work/05_families/family_table.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        if row['strain'] == 'PvP01' and row['gene_id'] in msa_ids and row['family'] in va_fams:
            seen.add(row['gene_id'])
print('\n'.join(sorted(seen)))
PY
)
N=$(echo "$PIDS" | wc -l)
log "$N variant-antigen genes with relaxed MSAs"

run_one() {
  local pid=$1
  local outdir=work/06_msa/core_relaxed_hyphy_va/$pid
  mkdir -p $outdir
  local aln=work/06_msa/core_relaxed_clean/${pid}.codon.cleaned.fa
  local tree=work/06_msa/core_relaxed_trees/${pid}/${pid}.treefile
  [ ! -s $aln ] && { echo "$pid: no aln"; return; }
  [ ! -s $tree ] && { echo "$pid: no tree"; return; }
  cp $aln  $outdir/aln.fa
  cp $tree $outdir/tree.nwk
  # 4 methods, each in its own JSON
  for METHOD in busted absrel meme fel; do
    local out=$outdir/${METHOD}.json
    [ -s $out ] && continue
    if [ "$METHOD" = "busted" ]; then
      # BUSTED with multi-hit Double for MNS robustness
      docker run --rm -v $V3:/v3 -w /v3/$outdir \
        quay.io/biocontainers/hyphy:2.5.99--h74d3ee0_0 \
        hyphy busted --alignment aln.fa --tree tree.nwk \
        --multiple-hits Double --output ${METHOD}.json \
        2>>$V3/logs/overnight/va_hyphy_${pid}.log \
        || echo "$pid busted-mh FAILED" >> $V3/logs/overnight/STATUS.md
    else
      docker run --rm -v $V3:/v3 -w /v3/$outdir \
        quay.io/biocontainers/hyphy:2.5.99--h74d3ee0_0 \
        hyphy $METHOD --alignment aln.fa --tree tree.nwk \
        --output ${METHOD}.json \
        2>>$V3/logs/overnight/va_hyphy_${pid}.log \
        || echo "$pid $METHOD FAILED" >> $V3/logs/overnight/STATUS.md
    fi
  done
  echo "$pid done"
}
export -f run_one
export V3

echo "$PIDS" | xargs -P 4 -I {} bash -c 'run_one "$@"' _ {}
DONE=$(find work/06_msa/core_relaxed_hyphy_va -name 'busted.json' 2>/dev/null | wc -l)
log "done; BUSTED-MH jsons: $DONE"
