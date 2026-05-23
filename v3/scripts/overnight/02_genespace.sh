#!/bin/bash
# GENESPACE pipeline: extract protein FASTAs for all 8 strains, then run GENESPACE.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/08_genespace/{proteomes,gffs,results} logs/overnight
log() { echo "$(date +%H:%M:%S) [genespace] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

# 1. Stage protein FASTAs for all 8 strains
for S in PvP01 Sal-I PvW1 PAM; do
  src=inputs/annotations/plasmodb-68/${S}.proteins.fa
  if [ -s $src ]; then
    cp $src work/08_genespace/proteomes/${S}.fa
    log "staged proteome $S (PlasmoDB)"
  fi
done
# PvSY56 PlasmoDB proteins file
[ -s inputs/annotations/plasmodb-68/PvSY56.proteins.fa ] && \
  cp inputs/annotations/plasmodb-68/PvSY56.proteins.fa work/08_genespace/proteomes/PvSY56.fa && \
  log "staged PvSY56 proteome"

# For PvT01/PvC01/MHC087: extract protein sequences from merged annotation GFFs via gffread
for S in PvT01 PvC01 MHC087; do
  out=work/08_genespace/proteomes/${S}.fa
  [ -s $out ] && { log "$S proteome already exists"; continue; }
  GFF=work/02d_merged/PvP01-as-ref/${S}.annotation.gff3
  FA=inputs/assemblies/${S}.fa
  if [ -s $GFF ] && [ -s $FA ]; then
    docker run --rm -v $V3:/v3 -w /v3 \
      quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4 \
      gffread -y $out -g $FA $GFF 2>> logs/overnight/genespace.log
    log "extracted $S proteome ($(grep -c '^>' $out 2>/dev/null) seqs)"
  fi
done

# 2. Stage GFFs (use PvP01-anchored Liftoff projections for cross-strain consistency)
for S in PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087; do
  if [ "$S" = "PvP01" ]; then
    cp inputs/annotations/plasmodb-68/PvP01.gff3 work/08_genespace/gffs/${S}.gff3
  else
    cp work/02d_merged/PvP01-as-ref/${S}.annotation.gff3 work/08_genespace/gffs/${S}.gff3
  fi
done

# 3. Run GENESPACE
log "running GENESPACE container"
docker run --rm -v $V3:/v3 -w /v3/work/08_genespace \
  doejgi/genespace:latest Rscript -e '
library(GENESPACE)
gpar <- init_genespace(
  wd = "/v3/work/08_genespace/results",
  path2mcscanx = "/MCScanX",
  rawGenomes = list(
    genomes = c("PvP01","Sal-I","PvW1","PAM","PvSY56","PvT01","PvC01","MHC087"),
    paths = list(
      fastaDir = "/v3/work/08_genespace/proteomes",
      gffDir   = "/v3/work/08_genespace/gffs"
    )
  )
)
out <- run_genespace(gpar, overwrite = TRUE)
' >> logs/overnight/genespace.log 2>&1
RC=$?
log "GENESPACE container exit $RC"
if [ -d work/08_genespace/results/pangenome ]; then
  log "pangenome dir exists; files: $(ls work/08_genespace/results/pangenome | wc -l)"
fi
log "done"
