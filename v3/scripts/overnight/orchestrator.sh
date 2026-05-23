#!/bin/bash
# Master overnight orchestrator. Runs all remaining plan items with proper dependency
# chains and writes a continuously-updated morning report.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p logs/overnight
STATUS=logs/overnight/STATUS.md
log() { echo "$(date +%Y-%m-%d_%H:%M:%S) | $1" | tee -a $STATUS; }
log "=== overnight orchestrator boot ==="

# === Section 1: launch load/timing logger ===
(
  while true; do
    LA=$(uptime | awk -F'load average: ' '{print $2}')
    NSHARDS=$(pgrep -c -f build_8way_msa_v2 || echo 0)
    NCMAP=$(pgrep -c -f 'CrossMap vcf' || echo 0)
    NDOCK=$(docker ps -q | wc -l)
    NIQTREE=$(pgrep -c -f iqtree || echo 0)
    NHYPHY=$(pgrep -c -f hyphy || echo 0)
    NTRIM=$(pgrep -c -f trimal || echo 0)
    echo "$(date +%H:%M:%S) load=$LA  H_shards=$NSHARDS  cmap=$NCMAP  docker=$NDOCK  iqtree=$NIQTREE  hyphy=$NHYPHY  trimal=$NTRIM" >> logs/overnight/load.log
    sleep 600
  done
) &
LOADGER_PID=$!
log "load logger pid $LOADGER_PID"

# === Section 2: independent jobs (sequential per-script, but in own background processes) ===
log "launching Multiz 8 hinges (sequential)"
nohup bash scripts/overnight/01_multiz_all.sh > logs/overnight/multiz.log 2>&1 &
MULTIZ_PID=$!

log "launching GENESPACE"
nohup bash scripts/overnight/02_genespace.sh > logs/overnight/genespace.log 2>&1 &
GENESPACE_PID=$!

# === Section 3: A2 finish watcher → 3-way comparison ===
(
  log "watcher: A2 finish"
  while pgrep -f run_a2_liftover.sh > /dev/null 2>&1; do sleep 60; done
  log "A2 orchestrator exited; verifying cohort VCFs"
  N=$(ls projection/A2_lastz/Pv4_cohort_on_*.vcf.gz 2>/dev/null | wc -l)
  log "A2: $N / 7 cohort VCFs written"
  if [ $N -ge 1 ]; then
    log "launching A1/A2/B 3-way comparison"
    bash scripts/overnight/07_a2_compare.sh
  fi
) &
A2WATCH_PID=$!

# === Section 4: Phase H → I → J chain ===
(
  log "chain: waiting for Phase H to finish"
  while pgrep -f build_8way_msa_v2 > /dev/null 2>&1; do sleep 60; done
  NMSA=$(ls work/06_msa/core_v3/*.codon.aln.fa 2>/dev/null | wc -l)
  log "Phase H done: $NMSA codon MSAs"

  log "launching Phase I (trimAl)"
  bash scripts/overnight/03_trimal.sh
  NCLEAN=$(ls work/06_msa/core_v3_clean/*.cleaned.fa 2>/dev/null | wc -l)
  log "Phase I done: $NCLEAN cleaned alignments"

  log "launching IQ-TREE"
  bash scripts/overnight/04_iqtree.sh
  NTREE=$(find work/06_msa/core_v3_trees -name '*.treefile' 2>/dev/null | wc -l)
  log "IQ-TREE done: $NTREE trees"

  # Wait briefly for validator agent to finish if it's still running
  for i in $(seq 1 30); do
    [ -s work/05_priorities/gene_priorities.tsv ] && break
    sleep 60
  done
  if [ ! -s work/05_priorities/gene_priorities.tsv ]; then
    log "no validator output after 30 min; using fallback priority list"
    [ -s work/05_priorities/fallback_priorities.tsv ] || python3 scripts/overnight/build_fallback_priorities.py
    cp work/05_priorities/fallback_priorities.tsv work/05_priorities/gene_priorities.tsv
  fi

  log "launching HyPhy priority bundle"
  bash scripts/overnight/05_hyphy_priority.sh

  log "launching HyPhy bulk BUSTED sweep"
  bash scripts/overnight/06_hyphy_bulk.sh
  NJSON=$(find work/06_msa/core_v3_hyphy -name 'busted.json' 2>/dev/null | wc -l)
  log "HyPhy done: $NJSON BUSTED jsons"
) &
HIJ_PID=$!

# === Section 5: wait for all dependent chains ===
log "main chains launched: multiz=$MULTIZ_PID genespace=$GENESPACE_PID a2watch=$A2WATCH_PID hij=$HIJ_PID"
log "orchestrator now sleeping until chains finish"
wait $MULTIZ_PID
log "multiz chain finished"
wait $GENESPACE_PID
log "genespace chain finished"
wait $A2WATCH_PID
log "A2-compare chain finished"
wait $HIJ_PID
log "phase H/I/J chain finished"

# Kill load logger
kill $LOADGER_PID 2>/dev/null || true

# === Section 6: morning report ===
log "writing MORNING_REPORT.md"
python3 scripts/overnight/morning_report.py > logs/overnight/MORNING_REPORT.md 2>&1 || true

log "=== overnight orchestrator complete ==="
