#!/usr/bin/env bash
set -euo pipefail
# Usage: run_in_container.sh <tool_name> <args...>
# Selects the correct pinned Docker image for <tool_name> and runs it with
# the WORK volume mounted, the correct user, and optional GPU passthrough.
#
# Environment variables consumed:
#   WORK         — workspace root (mounted at $WORK:$WORK inside the container)
#   SCRATCH      — scratch dir for large temp files (default /tmp)
#   USE_GPU      — set to 1 to enable --gpus passthrough
#   GPU          — GPU device index/UUID (default 0)

TOOL=$1; shift

case $TOOL in
  mash)
    IMG="staphb/mash:2.3" ;;
  busco)
    IMG="ezlabgva/busco:v5.8.0_cv1" ;;
  longdust)
    IMG="quay.io/biocontainers/longdust:1.4--h577a1d6_0" ;;
  sdust)
    IMG="quay.io/biocontainers/sdust:0.1--h077b44d_2" ;;
  bedtools)
    IMG="quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2" ;;
  samtools)
    IMG="quay.io/biocontainers/samtools:1.20--h50ea8bc_0" ;;
  kegalign)
    IMG="quay.io/biocontainers/kegalign-gpu:1.0.0--hdfd78af_0" ;;
  axtChain)        IMG="quay.io/biocontainers/ucsc-axtchain:482--h0b57e2e_2" ;;
  chainSort)       IMG="quay.io/biocontainers/ucsc-chainsort:482--h0b57e2e_0" ;;
  chainPreNet)     IMG="quay.io/biocontainers/ucsc-chainprenet:482--h0b57e2e_0" ;;
  chainNet)        IMG="quay.io/biocontainers/ucsc-chainnet:482--h0b57e2e_0" ;;
  netChainSubset)  IMG="quay.io/biocontainers/ucsc-netchainsubset:482--h0b57e2e_0" ;;
  chainStitchId)   IMG="quay.io/biocontainers/ucsc-chainstitchid:482--h0b57e2e_0" ;;
  chainSwap)       IMG="quay.io/biocontainers/ucsc-chainswap:482--h0b57e2e_0" ;;
  axtToMaf)        IMG="quay.io/biocontainers/ucsc-axttomaf:482--h0b57e2e_0" ;;
  liftoff)
    IMG="quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0" ;;
  agat)
    IMG="quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0" ;;
  toga|cesar)
    IMG="ghcr.io/hillerlab/toga:latest" ;;
  gffread)
    IMG="quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4" ;;
  pggb|wfmash|seqwish|smoothxg|odgi)
    IMG="ghcr.io/pangenome/pggb:latest" ;;
  mafft)
    IMG="quay.io/biocontainers/mafft:7.221--0" ;;
  pal2nal)
    IMG="quay.io/biocontainers/pal2nal:14--pl526_0" ;;
  trimal)
    IMG="quay.io/biocontainers/trimal:1.5.1--h9948957_0" ;;
  iqtree3)
    IMG="quay.io/biocontainers/iqtree:3.0.0--hdcf5f25_0" ;;
  hyphy)
    IMG="quay.io/biocontainers/hyphy:2.5.62--he91c24d_0" ;;
  multiz)
    IMG="quay.io/biocontainers/multiz:11.2--h470a237_0" ;;
  bcftools)
    IMG="quay.io/biocontainers/bcftools:1.20--h8b25389_0" ;;
  CrossMap)
    IMG="quay.io/biocontainers/crossmap:0.6.5--pyh7cba7a3_0" ;;
  python3)
    IMG="quay.io/biocontainers/pyfaidx:0.9.0.4--pyhdfd78af_0" ;;
  bgzip|tabix)
    IMG="quay.io/biocontainers/samtools:1.20--h50ea8bc_0" ;;
  *)
    echo "run_in_container.sh: unknown tool: $TOOL" >&2; exit 1 ;;
esac

GPU_ARG=""
[[ "${USE_GPU:-0}" == "1" ]] && GPU_ARG="--gpus device=${GPU:-0}"

SCRATCH_DIR="${SCRATCH:-/tmp}"

exec docker run --rm -i ${GPU_ARG} \
  -v "${WORK}:${WORK}" \
  -v "${SCRATCH_DIR}:${SCRATCH_DIR}" \
  -w "${WORK}" \
  -u "$(id -u):$(id -g)" \
  --entrypoint "" \
  "${IMG}" \
  "${TOOL}" "$@"
