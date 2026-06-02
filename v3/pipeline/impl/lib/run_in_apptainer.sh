#!/usr/bin/env bash
set -euo pipefail
# Usage: run_in_apptainer.sh <tool_name> <args...>
# HPC/cluster alternative to run_in_container.sh.  Uses Singularity/Apptainer
# instead of Docker (Docker is typically forbidden on shared HPC systems).
#
# Environment variables consumed:
#   WORK         — workspace root (bound into the container)
#   SCRATCH      — scratch dir (default /tmp)
#   USE_GPU      — set to 1 to pass --nv (NVIDIA GPU) flag
#   APPTAINER_CACHEDIR — where Apptainer caches pulled images (optional)

TOOL=$1; shift

case $TOOL in
  mash)
    IMG="docker://quay.io/biocontainers/mash:2.3--he348c14_4" ;;
  busco)
    IMG="docker://ezlabgva/busco:v5.7.1_cv1" ;;
  longdust|sdust)
    IMG="docker://quay.io/biocontainers/longdust:0.1.0--h5b5514e_0" ;;
  bedtools)
    IMG="docker://quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2" ;;
  samtools)
    IMG="docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0" ;;
  kegalign)
    IMG="docker://quay.io/biocontainers/kegalign-gpu:1.0.0--hdfd78af_0" ;;
  axtChain|chainSort|chainPreNet|chainNet|netChainSubset|chainStitchId|chainSwap|axtToMaf)
    IMG="docker://quay.io/biocontainers/ucsc-kent-tools:469--h664eb37_0" ;;
  liftoff)
    IMG="docker://quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0" ;;
  agat)
    IMG="docker://quay.io/biocontainers/agat:1.4.0--pl5321hdfd78af_0" ;;
  toga|cesar)
    IMG="docker://ghcr.io/hillerlab/toga:latest" ;;
  gffread)
    IMG="docker://quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4" ;;
  pggb|wfmash|seqwish|smoothxg|odgi)
    IMG="docker://ghcr.io/pangenome/pggb:202412130311080800a17" ;;
  mafft)
    IMG="docker://quay.io/biocontainers/mafft:7.526--h031d066_0" ;;
  pal2nal)
    IMG="docker://quay.io/biocontainers/pal2nal:14.1--pl5321hdfd78af_5" ;;
  trimal)
    IMG="docker://quay.io/biocontainers/trimal:1.4.1--h9948957_8" ;;
  iqtree3)
    IMG="docker://quay.io/biocontainers/iqtree:3.0.0--hdcf5f25_0" ;;
  hyphy)
    IMG="docker://quay.io/biocontainers/hyphy:2.5.62--he91c24d_0" ;;
  multiz)
    IMG="docker://quay.io/biocontainers/multiz:11.2--h470a237_0" ;;
  bcftools)
    IMG="docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0" ;;
  CrossMap)
    IMG="docker://quay.io/biocontainers/crossmap:0.6.5--pyh7cba7a3_0" ;;
  python3)
    IMG="docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0" ;;
  bgzip|tabix)
    IMG="docker://quay.io/biocontainers/htslib:1.20--h5efdd21_0" ;;
  *)
    echo "run_in_apptainer.sh: unknown tool: $TOOL" >&2; exit 1 ;;
esac

GPU_ARG=""
[[ "${USE_GPU:-0}" == "1" ]] && GPU_ARG="--nv"

SCRATCH_DIR="${SCRATCH:-/tmp}"

exec singularity exec \
  ${GPU_ARG} \
  --bind "${WORK}:${WORK}" \
  --bind "${SCRATCH_DIR}:${SCRATCH_DIR}" \
  --pwd "${WORK}" \
  "${IMG}" \
  "${TOOL}" "$@"
