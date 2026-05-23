#!/bin/bash
# BUSCO plasmodium_odb10 baseline for 8 proteomes.
# Light CPU usage allows parallel execution alongside Liftoff/KegAlign.
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
IMG=ezlabgva/busco:v5.8.0_cv1

mkdir -p $V3/work/00_inventory/busco
cd $V3/work/00_inventory/busco

# Download lineage once
if [ ! -d plasmodium_odb10 ]; then
  echo "=== downloading plasmodium_odb10 lineage ==="
  docker run --rm -v $V3/work/00_inventory/busco:/busco_wd -w /busco_wd $IMG \
    busco --download plasmodium_odb10 2>&1 | tail -5
  mv busco_downloads/lineages/plasmodium_odb10 . 2>/dev/null || true
fi

run_busco() {
  local s=$1
  local protfa=""
  case $s in
    Sal-I)   protfa=/v3/inputs/annotations/plasmodb-68/Sal-I.proteins.fa ;;
    PvP01)   protfa=/v3/inputs/annotations/plasmodb-68/PvP01.proteins.fa ;;
    PvW1)    protfa=/v3/inputs/annotations/plasmodb-68/PvW1.proteins.fa ;;
    PAM)     protfa=/v3/inputs/annotations/plasmodb-68/PAM.proteins.fa ;;
    PvSY56)  protfa=/v3/inputs/annotations/plasmodb-68/PvSY56.proteins.fa ;;
    PvT01)   protfa=/v3/inputs/annotations/ncbi-datasets/PvT01.proteins.fa ;;
    PvC01)   protfa=/v3/inputs/annotations/ncbi-datasets/PvC01.proteins.fa ;;
    MHC087)  protfa=/v3/inputs/annotations/ncbi-datasets/MHC087.proteins.fa ;;
  esac
  out=/v3/work/00_inventory/busco/${s}_proteins
  [ -d $V3/work/00_inventory/busco/${s}_proteins ] && [ -s $V3/work/00_inventory/busco/${s}_proteins/short_summary.specific.plasmodium_odb10.${s}_proteins.txt ] && {
    echo "[skip] $s"; return; }
  docker run --rm -v $V3:/v3 -w /v3/work/00_inventory/busco $IMG \
    busco -i $protfa -m proteins -l plasmodium_odb10 --offline \
          --download_path /v3/work/00_inventory/busco/busco_downloads \
          -o ${s}_proteins -c 2 -f 2>$V3/work/logs/busco_${s}.log
  ss=$out/short_summary.specific.plasmodium_odb10.${s}_proteins.txt
  if [ -s $V3/work/00_inventory/busco/${s}_proteins/short_summary.specific.plasmodium_odb10.${s}_proteins.txt ]; then
    line=$(grep 'C:' $V3/work/00_inventory/busco/${s}_proteins/short_summary.specific.plasmodium_odb10.${s}_proteins.txt | head -1)
    printf "  %-8s  %s\n" "$s" "$line"
  fi
}
export V3 IMG
export -f run_busco

echo "$1 $2 $3 $4 $5 $6 $7 $8" | tr ' ' '\n' | xargs -P 4 -I {} bash -c 'run_busco "$@"' _ {}

echo
echo "===== BUSCO summary ====="
for s in PvP01 PvW1 PAM Sal-I PvSY56 PvT01 PvC01 MHC087; do
  ss=$V3/work/00_inventory/busco/${s}_proteins/short_summary.specific.plasmodium_odb10.${s}_proteins.txt
  if [ -s $ss ]; then
    line=$(grep 'C:' $ss | head -1)
    printf "  %-8s  %s\n" "$s" "$line"
  else
    printf "  %-8s  MISSING\n" "$s"
  fi
done
