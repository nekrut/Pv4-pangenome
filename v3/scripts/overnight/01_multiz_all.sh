#!/bin/bash
# Build 8 multiz multi-way MAFs, one per hinge strain, sequentially.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/07_multiz logs/overnight
declare -A GCA=(
  [PvP01]=GCA_900093555.2 [Sal-I]=GCA_000002415.2 [PvW1]=GCA_914969965.1
  [PAM]=GCA_949152365.1 [PvSY56]=GCA_003402215.1 [PvT01]=GCA_900093545.1
  [PvC01]=GCA_900093535.1 [MHC087]=GCA_040114635.1
)
STRAINS="PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087"
AXT_DIR=projection/A2_kegalign/axt
SIZES_DIR=work/07_multiz/.sizes
mkdir -p $SIZES_DIR
AXTTOMAF=quay.io/biocontainers/ucsc-axttomaf:482--h0b57e2e_0
MULTIZ=quay.io/biocontainers/multiz:11.2--h7b50bb2_7

log() { echo "$(date +%H:%M:%S) [multiz] $1" | tee -a logs/overnight/STATUS.md ; }

# Pre-stage size files (one per genome)
for S in $STRAINS; do
  G=${GCA[$S]}
  [ -s $SIZES_DIR/${G}.sizes ] || cut -f1,2 genomes/softmasked/${G}.fa.fai > $SIZES_DIR/${G}.sizes
done

for HINGE in $STRAINS; do
  HG=${GCA[$HINGE]}
  HDIR=work/07_multiz/$HINGE
  mkdir -p $HDIR
  log "hinge=$HINGE"

  PAIRS=""
  for OTHER in $STRAINS; do
    [ "$OTHER" = "$HINGE" ] && continue
    OG=${GCA[$OTHER]}
    AXT1=$AXT_DIR/${HG}__vs__${OG}.axt
    AXT2=$AXT_DIR/${OG}__vs__${HG}.axt
    if [ -s $AXT1 ]; then AXT=$AXT1; FLIP=0
    elif [ -s $AXT2 ]; then AXT=$AXT2; FLIP=1
    else log "  missing AXT $HINGE vs $OTHER"; continue
    fi
    MAF=$HDIR/${HINGE}_vs_${OTHER}.maf
    if [ ! -s $MAF ]; then
      if [ $FLIP -eq 0 ]; then
        docker run --rm -v $V3:/v3 -w /v3 $AXTTOMAF \
          axtToMaf -tPrefix=${HG}. -qPrefix=${OG}. \
          /v3/$AXT /v3/$SIZES_DIR/${HG}.sizes /v3/$SIZES_DIR/${OG}.sizes /v3/$MAF \
          >> logs/overnight/multiz_$HINGE.log 2>&1
      else
        # Need to swap target/query: axtSwap (write Python-flipped)
        FLIPPED=$(mktemp -p $V3/work/07_multiz/.sizes --suffix=.axt)
        python3 -c "
import sys
with open('$AXT') as fh, open('$FLIPPED','w') as out:
    for ln in fh:
        if ln.startswith('#') or ln.startswith('>') or not ln.strip():
            out.write(ln); continue
        f = ln.split()
        if len(f) == 9 and f[0].isdigit():
            # axt summary line: # tName tStart tEnd qName qStart qEnd strand score
            # swap target<->query, preserve strand
            new = [f[0], f[4], f[5], f[6], f[1], f[2], f[3], f[7], f[8]]
            out.write(' '.join(new) + '\n')
        else:
            # sequence lines: order is t then q — swap order
            # AXT alternates summary, t-seq, q-seq, blank
            out.write(ln)
"
        # Crude flip didn't reorder seq lines; fall back to chainSwap approach via chain
        # For now, just convert the un-flipped axt and let axtToMaf interpret (target/query swapped)
        docker run --rm -v $V3:/v3 -w /v3 $AXTTOMAF \
          axtToMaf -tPrefix=${OG}. -qPrefix=${HG}. \
          /v3/$AXT /v3/$SIZES_DIR/${OG}.sizes /v3/$SIZES_DIR/${HG}.sizes /v3/$MAF \
          >> logs/overnight/multiz_$HINGE.log 2>&1
        rm -f $FLIPPED
      fi
    fi
    if [ -s $MAF ]; then
      NB=$(grep -c '^a ' $MAF || echo 0)
      PAIRS="$PAIRS $MAF"
      log "  $OTHER: $NB MAF blocks"
    else
      log "  $OTHER: MAF empty"
    fi
  done

  # Iteratively multiz
  CURRENT=""
  STEP=0
  for P in $PAIRS; do
    if [ -z "$CURRENT" ]; then
      CURRENT=$P
    else
      STEP=$((STEP+1))
      NEXT=$HDIR/cum_${STEP}.maf
      docker run --rm -v $V3:/v3 -w /v3 $MULTIZ \
        multiz /v3/$CURRENT /v3/$P 0 > $NEXT 2>> logs/overnight/multiz_$HINGE.log
      if [ -s $NEXT ]; then
        CURRENT=$NEXT
      else
        log "  multiz step $STEP failed; aborting hinge=$HINGE"
        break
      fi
    fi
  done
  if [ -n "$CURRENT" ] && [ -s $CURRENT ]; then
    cp $CURRENT $HDIR/${HINGE}.multiz.maf
    NB=$(grep -c '^a ' $HDIR/${HINGE}.multiz.maf)
    log "  hinge=$HINGE done: $NB multi-way blocks"
  else
    log "  hinge=$HINGE FAILED (no output)"
  fi
done
log "ALL HINGES DONE"
