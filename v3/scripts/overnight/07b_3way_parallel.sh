#!/bin/bash
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p work/09_projection_compare
BCF=/home/anton/miniconda3/envs/bcfmod/bin/bcftools
OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"

query_one() {
  local O=$1 SRC=$2
  local W=work/09_projection_compare/$O
  mkdir -p $W
  local out=$W/${SRC}.sites.tsv
  [ -s $out ] && return 0
  local VCF
  case $SRC in
    A1) VCF=projection/A1_wfmash/Pv4_cohort_on_${O}.vcf.gz ;;
    A2) VCF=projection/A2_lastz/Pv4_cohort_on_${O}.vcf.gz ;;
    B)  VCF=projection/B_graph/Pv4_cohort_on_${O}.vcf.gz ;;
  esac
  $BCF query -f '%CHROM\t%POS\t%REF\t%ALT\n' $VCF 2>/dev/null | LC_ALL=C sort -u -S 2G > $out
  echo "$O.$SRC done"
}
export -f query_one
export BCF V3

JOBS=$(for O in $OTHERS; do for SRC in A1 A2 B; do echo "$O $SRC"; done; done)
echo "21 queries (8 parallel)..."
echo "$JOBS" | xargs -P 8 -n 2 bash -c 'query_one "$@"' _

echo
echo "intersection..."
OUT=work/09_projection_compare/intersection.tsv
echo -e "target\tA1_sites\tA2_sites\tB_sites\tA1_only\tA2_only\tB_only\tA1_A2\tA1_B\tA2_B\tA1_A2_B" > $OUT
for O in $OTHERS; do
  W=work/09_projection_compare/$O
  A1N=$(wc -l < $W/A1.sites.tsv)
  A2N=$(wc -l < $W/A2.sites.tsv)
  BN=$(wc -l < $W/B.sites.tsv)
  A1A2=$(LC_ALL=C comm -12 $W/A1.sites.tsv $W/A2.sites.tsv | wc -l)
  A1B=$(LC_ALL=C comm -12 $W/A1.sites.tsv $W/B.sites.tsv | wc -l)
  A2B=$(LC_ALL=C comm -12 $W/A2.sites.tsv $W/B.sites.tsv | wc -l)
  ALL3=$(LC_ALL=C comm -12 $W/A1.sites.tsv $W/A2.sites.tsv | LC_ALL=C comm -12 - $W/B.sites.tsv | wc -l)
  A1only=$((A1N - A1A2 - A1B + ALL3))
  A2only=$((A2N - A1A2 - A2B + ALL3))
  Bonly=$((BN - A1B - A2B + ALL3))
  echo -e "${O}\t${A1N}\t${A2N}\t${BN}\t${A1only}\t${A2only}\t${Bonly}\t${A1A2}\t${A1B}\t${A2B}\t${ALL3}" >> $OUT
done
echo "DONE"
cat $OUT | column -t
