#!/bin/bash
# Rebuild multiz multi-way MAFs for the 3 hinges that failed (PvT01, PvC01, MHC087).
# Iteratively multiz the existing pairwise MAFs. Output to local dir.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p logs/overnight
log() { echo "$(date +%H:%M:%S) [multiz_repair] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

MULTIZ=quay.io/biocontainers/multiz:11.2--h7b50bb2_7
TMP=/media/anton/scratch/multiz_repair_tmp
mkdir -p $TMP

declare -A GCA=(
  [PvP01]=GCA_900093555.2 [Sal-I]=GCA_000002415.2 [PvW1]=GCA_914969965.1
  [PAM]=GCA_949152365.1 [PvSY56]=GCA_003402215.1 [PvT01]=GCA_900093545.1
  [PvC01]=GCA_900093535.1 [MHC087]=GCA_040114635.1
)
STRAINS="PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087"

for HINGE in PvT01 PvC01 MHC087; do
  log "hinge=$HINGE"
  HDIR=work/07_multiz/$HINGE
  PAIRS=""
  for OTHER in $STRAINS; do
    [ "$OTHER" = "$HINGE" ] && continue
    MAF=$HDIR/${HINGE}_vs_${OTHER}.maf
    if [ -s $MAF ]; then
      PAIRS="$PAIRS $MAF"
    fi
  done
  log "  pairs available: $(echo $PAIRS | wc -w)"

  CURRENT=""
  STEP=0
  for P in $PAIRS; do
    if [ -z "$CURRENT" ]; then
      CURRENT=$P
    else
      STEP=$((STEP+1))
      NEXT=$TMP/${HINGE}_cum_${STEP}.maf
      docker run --rm -v $V3:$V3 -v $TMP:$TMP -w $V3 $MULTIZ \
        multiz $CURRENT $P 0 > $NEXT 2>> logs/overnight/multiz_repair_$HINGE.log
      if [ -s $NEXT ]; then
        SIZE=$(du -h $NEXT | awk '{print $1}')
        log "    step $STEP: $SIZE"
        CURRENT=$NEXT
      else
        log "  step $STEP FAILED at $(basename $P)"
        break
      fi
    fi
  done
  if [ -n "$CURRENT" ] && [ -s $CURRENT ]; then
    cp $CURRENT $HDIR/${HINGE}.multiz.maf
    NB=$(grep -c '^a ' $HDIR/${HINGE}.multiz.maf)
    log "  hinge=$HINGE done: $NB multi-way blocks"
    # Clean cum_*.maf intermediates for this hinge
    rm -f $TMP/${HINGE}_cum_*.maf
  else
    log "  hinge=$HINGE FAILED — no final MAF"
  fi
done
log "ALL DONE"
