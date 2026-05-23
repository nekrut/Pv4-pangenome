#!/bin/bash
# Parallel version: re-use build_one from 10c, but run 6 targets concurrently
# (one xargs slot per target).
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
SCRATCH=/media/anton/scratch
SORT_TMP=$SCRATCH/path_b_sort_tmp
HDR_DIR=/tmp/path_b_headers
mkdir -p $SORT_TMP $HDR_DIR
log() { echo "$(date +%H:%M:%S) [path_b_par] $1" | tee -a logs/overnight/STATUS.md; }

# Targets remaining (skip target 1, which is done)
OTHERS="GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"

# Pre-build header snippets
for O in $OTHERS; do
  [ -s $HDR_DIR/${O}.contigs ] || awk -v acc=$O '{print "##contig=<ID="$1",length="$2",assembly="acc".fa>"}' \
    genomes/softmasked/${O}.fa.fai > $HDR_DIR/${O}.contigs
done

# Source build_one from 10c
. scripts/overnight/10c_path_b_final.sh 2>/dev/null

log "starting parallel build for 6 targets"
echo "$OTHERS" | tr ' ' '\n' | xargs -P 6 -I {} bash -c '
  V3=/media/anton/data/sandbox/Pv4/v3
  cd $V3
  source scripts/overnight/10c_path_b_final.sh > /dev/null 2>&1
  build_one "$@"
' _ {} 2>&1
log "ALL TARGETS DONE"
