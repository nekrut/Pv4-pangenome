#!/bin/bash
# Wait for all repair-phase background jobs to finish, then regenerate the
# morning report.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
log() { echo "$(date +%H:%M:%S) [finalize] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting; waiting for hyphy + multiz_repair to drain"

# Poll until all queued repair jobs have exited
while true; do
  NH=$(pgrep -c -f "hyphy " || echo 0)
  NM=$(pgrep -c -f multiz || echo 0)
  NI=$(pgrep -c -f iqtree || echo 0)
  if [ "$NH" = "0" ] && [ "$NM" = "0" ] && [ "$NI" = "0" ]; then
    break
  fi
  sleep 60
done
log "all repair jobs finished"

# Regenerate morning report
python3 scripts/overnight/morning_report.py > logs/overnight/MORNING_REPORT.md 2>&1 || true
log "MORNING_REPORT.md regenerated"
log "DONE"
