#!/bin/bash
# OrthoFinder3 — substitute for GENESPACE. Gives orthogroups + species tree.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p logs/overnight
log() { echo "$(date +%H:%M:%S) [orthofinder] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

INPUT=$V3/work/08_orthofinder/input
OUT=$V3/work/08_orthofinder/results
rm -rf $INPUT $OUT 2>/dev/null
mkdir -p $INPUT

# Use the already-cleaned (gene-level-header) proteomes
for S in PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087; do
  src=work/08_genespace/proteomes/$S.fa
  [ -s $src ] || src=work/08_genespace/clean_proteomes/$S.fa
  if [ -s $src ]; then
    cp $src $INPUT/$S.fa
    log "staged $S ($(grep -c '^>' $INPUT/$S.fa) seqs)"
  else
    log "MISSING proteome for $S"
  fi
done

log "running OrthoFinder3"
docker run --rm -v $V3:/v3 -w /v3 \
  quay.io/biocontainers/orthofinder:3.1.4--hdfd78af_0 \
  orthofinder -f /v3/work/08_orthofinder/input -o /v3/work/08_orthofinder/results \
              -t 8 -a 4 -M dendroblast \
  >> logs/overnight/orthofinder.log 2>&1
RC=$?
log "OrthoFinder exit $RC"
if [ -d $OUT ]; then
  RES=$(find $OUT -type d -name "Results_*" | head -1)
  if [ -n "$RES" ]; then
    N=$(wc -l < $RES/Orthogroups/Orthogroups.tsv 2>/dev/null)
    log "orthogroups: $N rows in $RES/Orthogroups/Orthogroups.tsv"
    NSC=$(wc -l < $RES/Orthogroups/Orthogroups_SingleCopyOrthologues.txt 2>/dev/null)
    log "single-copy orthogroups: $NSC"
  fi
fi
log "done"
