#!/bin/bash
# A1 vs A2 vs B 3-way comparison: per-site intersection metrics per target reference.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/09_projection_compare logs/overnight
log() { echo "$(date +%H:%M:%S) [a2_compare] $1" | tee -a logs/overnight/STATUS.md ; }
log "starting"

OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"
BCF=/home/anton/miniconda3/envs/bcfmod/bin/bcftools

OUT=work/09_projection_compare/intersection.tsv
echo -e "target\tA1_sites\tA2_sites\tB_sites\tA1_only\tA2_only\tB_only\tA1_A2\tA1_B\tA2_B\tA1_A2_B" > $OUT

for O in $OTHERS; do
  A1=projection/A1_wfmash/Pv4_cohort_on_${O}.vcf.gz
  A2=projection/A2_lastz/Pv4_cohort_on_${O}.vcf.gz
  B=projection/B_graph/Pv4_cohort_on_${O}.vcf.gz
  [ ! -s $A1 ] && { log "$O missing A1"; continue; }
  [ ! -s $A2 ] && { log "$O missing A2"; continue; }
  [ ! -s $B  ] && { log "$O missing B"; continue; }

  WORK=work/09_projection_compare/${O}
  mkdir -p $WORK
  # Extract CHROM POS REF ALT per source
  for SRC in A1 A2 B; do
    case $SRC in A1) VCF=$A1;; A2) VCF=$A2;; B) VCF=$B;; esac
    $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\n' $VCF 2>/dev/null | sort -u > $WORK/${SRC}.sites.tsv
  done

  A1N=$(wc -l < $WORK/A1.sites.tsv); A2N=$(wc -l < $WORK/A2.sites.tsv); BN=$(wc -l < $WORK/B.sites.tsv)
  # Intersections
  A1A2=$(comm -12 $WORK/A1.sites.tsv $WORK/A2.sites.tsv | wc -l)
  A1B=$(comm -12 $WORK/A1.sites.tsv $WORK/B.sites.tsv | wc -l)
  A2B=$(comm -12 $WORK/A2.sites.tsv $WORK/B.sites.tsv | wc -l)
  ALL3=$(comm -12 $WORK/A1.sites.tsv $WORK/A2.sites.tsv | comm -12 - $WORK/B.sites.tsv | wc -l)
  A1only=$((A1N - A1A2 - A1B + ALL3))
  A2only=$((A2N - A1A2 - A2B + ALL3))
  Bonly=$((BN - A1B - A2B + ALL3))
  echo -e "${O}\t${A1N}\t${A2N}\t${BN}\t${A1only}\t${A2only}\t${Bonly}\t${A1A2}\t${A1B}\t${A2B}\t${ALL3}" >> $OUT
  log "$O: A1=$A1N A2=$A2N B=$BN ALL3=$ALL3"
done

log "done; $OUT"
