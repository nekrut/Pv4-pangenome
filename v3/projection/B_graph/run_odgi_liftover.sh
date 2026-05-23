#!/bin/bash
set -euo pipefail
V3=/media/anton/data/sandbox/Pv4/v3
SCRATCH=/media/anton/scratch/malariagen_pv4
OG=/media/anton/data/sandbox/Pv4/v2/pggb_out/pggb_input.fa.gz.f705205.c28ecf8.dcad0e6.smooth.fix.og
ODGI=/home/anton/miniconda3/envs/pggb/bin/odgi
RENAME=$V3/projection/A1_wfmash/rename_pvp01_v1_to_genbank.tsv
CHROMS="01 02 03 04 05 06 07 08 09 10 11 12 13 14 API MIT"
OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"

# Step 1: build cohort site BED with graph-path names
mkdir -p $V3/projection/B_graph/sites
ALL_BED=$V3/projection/B_graph/sites/all_mg_sites.bed
> $ALL_BED
for chr in $CHROMS; do
  # Extract CHROM,POS from MalariaGEN VCF, rename to LT635xxx, prepend graph path prefix
  docker run --rm -v $SCRATCH:/mg -v $V3:/v3 quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bcftools view /mg/Pv4_PvP01_${chr}_v1.vcf.gz -f PASS \
    | awk -v rfile=/v3/projection/A1_wfmash/rename_pvp01_v1_to_genbank.tsv \
      'BEGIN{while((getline ln < rfile)>0){split(ln,a,"\t");map[a[1]]=a[2]}} !/^#/{lt=map[$1]; if(lt) print "GCA_900093555.2#1#" lt "#0\t" $2-1 "\t" $2}' \
    >> $ALL_BED
done
n_sites=$(wc -l < $ALL_BED)
echo "Total MG PASS sites: $n_sites"

# Step 2: per target, run odgi position
for OTHER in $OTHERS; do
  echo "===== odgi position to $OTHER ====="
  /usr/bin/time -v $ODGI position -t 18 \
       -i $OG -b $ALL_BED \
       -R $V3/projection/B_graph/path_lists/${OTHER}.paths.txt \
       > $V3/projection/B_graph/sites/mg_on_${OTHER}.tsv \
       2> $V3/logs/B_odgi_${OTHER}.log
  n_lifted=$(wc -l < $V3/projection/B_graph/sites/mg_on_${OTHER}.tsv)
  printf "%s  lifted=%d/%d (%.1f%%)\n" "$OTHER" "$n_lifted" "$n_sites" "$(awk -v l=$n_lifted -v t=$n_sites 'BEGIN{printf "%.1f", 100*l/t}')"
done

echo
echo "=== Path B step 2 summary ==="
for OTHER in $OTHERS; do
  n=$(wc -l < $V3/projection/B_graph/sites/mg_on_${OTHER}.tsv 2>/dev/null)
  printf "  %-22s  positions_lifted=%d\n" "$OTHER" "$n"
done
