#!/bin/bash
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
log() { echo "$(date +%H:%M:%S) [finalize_main] $1" | tee -a logs/overnight/STATUS.md; }
log "starting; waiting for A2 concats"
while [ "$(ls projection/A2_lastz/Pv4_cohort_on_*.vcf.gz 2>/dev/null | wc -l)" -lt 7 ]; do
  sleep 60
done
log "A2: all 7 cohorts present; running 3-way comparison"
bash scripts/overnight/07_a2_compare.sh
log "running Phase 5 drug-resistance QC"
python3 scripts/overnight/09_phase5_drug_resistance_qc.py >> logs/overnight/phase5.log 2>&1
log "MAIN PLAN COMPLETE"
