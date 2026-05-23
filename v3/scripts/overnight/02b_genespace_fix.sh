#!/bin/bash
# Fixed GENESPACE pipeline. Uses correct v1.2.3 API.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p logs/overnight
log() { echo "$(date +%H:%M:%S) [genespace_fix] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

WD=$V3/work/08_genespace
RAW=$WD/rawGenomes
rm -rf $WD/results 2>/dev/null
mkdir -p $RAW $WD/results

STRAINS="PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087"

# Lay out the per-strain dirs in GENESPACE's expected structure:
#   rawGenomes/<strain>/<strain>/annotation/{<strain>.gff3, <strain>.fa(proteins)}
log "staging files in rawGenomes layout"
for S in $STRAINS; do
  ADIR=$RAW/$S/$S/annotation
  mkdir -p $ADIR
  # GFF
  if [ "$S" = "PvP01" ] && [ -s inputs/annotations/plasmodb-68/$S.gff3 ]; then
    cp inputs/annotations/plasmodb-68/$S.gff3 $ADIR/$S.gff3
  elif [ -s work/02d_merged/PvP01-as-ref/$S.annotation.gff3 ]; then
    cp work/02d_merged/PvP01-as-ref/$S.annotation.gff3 $ADIR/$S.gff3
  fi
  # Proteome — prefer PlasmoDB, fall back to gffread extraction
  if [ -s inputs/annotations/plasmodb-68/$S.proteins.fa ]; then
    cp inputs/annotations/plasmodb-68/$S.proteins.fa $ADIR/$S.fa
  elif [ -s work/08_genespace/proteomes/$S.fa ]; then
    cp work/08_genespace/proteomes/$S.fa $ADIR/$S.fa
  fi
  log "  $S: gff=$([ -s $ADIR/$S.gff3 ] && echo OK || echo MISS), fa=$([ -s $ADIR/$S.fa ] && echo OK || echo MISS)"
done

# Run GENESPACE
log "running GENESPACE container"
docker run --rm -v $V3:/v3 -w /v3 \
  doejgi/genespace:latest Rscript -e "
library(GENESPACE)
parse_annotations(
  rawGenomeRepo = '/v3/work/08_genespace/rawGenomes',
  genomeDirs = c('PvP01','Sal-I','PvW1','PAM','PvSY56','PvT01','PvC01','MHC087'),
  genomeIDs = c('PvP01','SalI','PvW1','PAM','PvSY56','PvT01','PvC01','MHC087'),
  genespaceWd = '/v3/work/08_genespace/results',
  gffString = 'gff3\$|gff3.gz\$',
  faString = 'fa\$|fasta\$',
  gffIdColumn = 'ID',
  headerEntryIndex = 1,
  overwrite = TRUE
)
gpar <- init_genespace(
  wd = '/v3/work/08_genespace/results',
  path2mcscanx = '/MCScanX',
  path2orthofinder = 'orthofinder',
  path2diamond = 'diamond',
  nCores = 8
)
out <- run_genespace(gpar, overwrite = TRUE)
saveRDS(out, '/v3/work/08_genespace/results/genespace_out.rds')
" >> logs/overnight/genespace_fix.log 2>&1
RC=$?
log "GENESPACE container exit $RC"

# Inspect outputs
if [ -d $WD/results/pangenome ]; then
  log "pangenome dir created; files: $(ls $WD/results/pangenome 2>/dev/null | wc -l)"
fi
if [ -d $WD/results/results ]; then
  log "results subdir: $(ls $WD/results/results 2>/dev/null | wc -l) files"
fi
log "done"
