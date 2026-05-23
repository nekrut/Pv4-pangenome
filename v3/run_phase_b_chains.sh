#!/bin/bash
# Phase B: build chains from KegAlign axt files.
#
# Per directed pair (T → Q):
#   axt → axtChain → chainSort → chainPreNet → chainNet → netChainSubset → cleaned.chain
#
# Plus per unordered pair: reciprocal-best chain via chainNet from both directions.
#
# Inputs:
#   projection/A2_kegalign/axt/${T}__vs__${Q}.axt    (target=T, query=Q)
#   projection/A2_kegalign/2bit/{T,Q}.2bit
#
# Outputs:
#   work/01_chains/${T}.${Q}.cleaned.chain           directed (T as target)
#   work/01_chains/${T}.${Q}.rbest.chain             reciprocal-best
#
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
ACCS=(GCA_000002415.2 GCA_900093555.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1)

# strain → genbank acc map
declare -A ACC
ACC[Sal-I]=GCA_000002415.2
ACC[PvP01]=GCA_900093555.2
ACC[PvT01]=GCA_900093545.1
ACC[PvC01]=GCA_900093535.1
ACC[PvW1]=GCA_914969965.1
ACC[PAM]=GCA_949152365.1
ACC[PvSY56]=GCA_003402215.1
ACC[MHC087]=GCA_040114635.1

mkdir -p $V3/work/01_chains

# Helper: build a single directed chain
build_directed_chain() {
  local T=$1 Q=$2     # accession-form
  local AXT=$V3/projection/A2_kegalign/axt/${T}__vs__${Q}.axt
  local OUTCHAIN=$V3/work/01_chains/${T}.${Q}.cleaned.chain
  if [ ! -s $AXT ]; then return 1; fi
  if [ -s $OUTCHAIN ]; then return 0; fi
  local T2BIT=$V3/projection/A2_kegalign/2bit/${T}.2bit
  local Q2BIT=$V3/projection/A2_kegalign/2bit/${Q}.2bit
  local TSIZE=$V3/work/01_chains/${T}.sizes
  local QSIZE=$V3/work/01_chains/${Q}.sizes
  [ ! -s $TSIZE ] && twoBitInfo $T2BIT $TSIZE 2>/dev/null || true
  [ ! -s $QSIZE ] && twoBitInfo $Q2BIT $QSIZE 2>/dev/null || true
  local TMP=$V3/work/01_chains/_tmp_${T}_${Q}
  $V3/tools/axtChain -linearGap=loose $AXT $T2BIT $Q2BIT ${TMP}.raw 2>/dev/null
  $V3/tools/chainSort ${TMP}.raw ${TMP}.sort.chain
  # chainPreNet to filter useless chains
  $V3/tools/chainPreNet ${TMP}.sort.chain $TSIZE $QSIZE ${TMP}.prenet.chain 2>/dev/null
  # chainNet — produces target.net + query.net
  $V3/tools/chainNet ${TMP}.prenet.chain $TSIZE $QSIZE ${TMP}.tnet ${TMP}.qnet 2>/dev/null
  # netChainSubset to extract chains in the net (cleaned)
  $V3/tools/netChainSubset ${TMP}.tnet ${TMP}.prenet.chain ${TMP}.subset.chain 2>/dev/null
  # chainStitchId to merge same-ID chain fragments (required by TOGA2)
  $V3/tools/chainStitchId ${TMP}.subset.chain $OUTCHAIN 2>/dev/null
  n_chains=$(grep -c '^chain ' $OUTCHAIN 2>/dev/null || echo 0)
  echo "  [done] $T → $Q  chains=$n_chains"
  rm -f ${TMP}.* 2>/dev/null
  return 0
}

# Helper: build reciprocal-best chain for an unordered pair
build_rbest_chain() {
  local T=$1 Q=$2
  local OUTRB=$V3/work/01_chains/${T}.${Q}.rbest.chain
  if [ -s $OUTRB ]; then return 0; fi
  local CHAIN_TQ=$V3/work/01_chains/${T}.${Q}.cleaned.chain
  local CHAIN_QT=$V3/work/01_chains/${Q}.${T}.cleaned.chain
  if [ ! -s $CHAIN_TQ ] || [ ! -s $CHAIN_QT ]; then return 1; fi
  local TSIZE=$V3/work/01_chains/${T}.sizes
  local QSIZE=$V3/work/01_chains/${Q}.sizes
  local TMP=$V3/work/01_chains/_rbtmp_${T}_${Q}
  # Get target.net from CHAIN_TQ
  $V3/tools/chainPreNet $CHAIN_TQ $TSIZE $QSIZE ${TMP}.tq.prenet 2>/dev/null
  $V3/tools/chainNet ${TMP}.tq.prenet $TSIZE $QSIZE ${TMP}.tq.tnet ${TMP}.tq.qnet 2>/dev/null
  # Swap CHAIN_QT to have T as target, intersect
  $V3/tools/chainSwap $CHAIN_QT ${TMP}.qt.swap
  # Sort by target
  $V3/tools/chainSort ${TMP}.qt.swap ${TMP}.qt.swap.sort
  # Reciprocal best = subset of CHAIN_TQ that's also in the swapped CHAIN_QT's chains
  $V3/tools/netChainSubset ${TMP}.tq.tnet ${TMP}.qt.swap.sort ${TMP}.rb 2>/dev/null
  $V3/tools/chainSort ${TMP}.rb $OUTRB
  n=$(grep -c '^chain ' $OUTRB 2>/dev/null || echo 0)
  echo "  [done] rbest ${T}↔${Q}  chains=$n"
  rm -f ${TMP}.* 2>/dev/null
  return 0
}

export V3
export -f build_directed_chain build_rbest_chain

# Generate all 56 directed pairs
PAIRS=$(python3 -c "
accs = '${ACCS[*]}'.split()
for t in accs:
  for q in accs:
    if t == q: continue
    print(f'{t} {q}')
")

# twoBitInfo helper — use a docker container if local UCSC binaries are missing
twoBitInfo() {
  if [ -x $V3/tools/twoBitInfo ]; then
    $V3/tools/twoBitInfo "$@"
  else
    docker run --rm -v $V3:/v3 quay.io/biocontainers/ucsc-twobitinfo:469--h664eb37_1 \
      twoBitInfo "$@" 2>/dev/null
  fi
}
export -f twoBitInfo

echo "=== Phase B Stage 1: 56 directed chains (parallel -P 8) ==="
echo "$PAIRS" | xargs -P 8 -n 2 bash -c 'build_directed_chain "$@"' _

echo
echo "=== Phase B Stage 2: 28 reciprocal-best chains ==="
RPAIRS=$(python3 -c "
accs = '${ACCS[*]}'.split()
for i,t in enumerate(accs):
  for q in accs[i+1:]:
    print(f'{t} {q}')
")
echo "$RPAIRS" | xargs -P 4 -n 2 bash -c 'build_rbest_chain "$@"' _

echo
echo "===== summary ====="
n_directed=$(ls $V3/work/01_chains/*.cleaned.chain 2>/dev/null | wc -l)
n_rbest=$(ls $V3/work/01_chains/*.rbest.chain 2>/dev/null | wc -l)
echo "directed cleaned: $n_directed / 56"
echo "reciprocal-best: $n_rbest / 28"

touch $V3/work/logs/checkpoints/01_chains.done
