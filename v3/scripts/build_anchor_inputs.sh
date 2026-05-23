#!/bin/bash
# Build per-anchor BED12 + isoforms + family list for Phase D.
# Normalizes GFF chrom names to match FASTA where needed.
set -uo pipefail

V3=/media/anton/data/sandbox/Pv4/v3
GFFREAD_IMG=quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4

build_one() {
  local A=$1
  local agff_orig=$V3/inputs/annotations/plasmodb-68/${A}.gff3
  local afa=$V3/inputs/assemblies/${A}.fa
  local agff_fixed=$V3/inputs/annotations/${A}.fixed.gff3
  local bed12=$V3/inputs/annotations/${A}.bed12
  local iso=$V3/inputs/annotations/${A}.isoforms.tsv

  # Get FASTA chrom set (first whitespace-token after >)
  local fa_chroms=$(mktemp)
  grep '^>' $afa | awk '{print substr($1,2)}' | sort -u > $fa_chroms

  # Get GFF chrom set
  local gff_chroms=$(mktemp)
  awk -F'\t' '!/^#/ {print $1}' $agff_orig | sort -u > $gff_chroms

  # Check if GFF chroms appear in FASTA (with or without .1 suffix)
  # If GFF chroms need .1, build mapping; else copy as-is
  local needs_suffix=$(awk -F'\t' 'NR==FNR{f[$1]=1; next} {
    if (!($1 in f) && (($1 ".1") in f)) print "yes"; exit
  }' $fa_chroms $gff_chroms)

  if [ "$needs_suffix" = "yes" ]; then
    echo "[$A] adding .1 suffix to GFF chrom names"
    awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} {
      orig=$1; new=$1 ".1"
      # but only if orig is not already in fasta but new is
      $1=new; print
    }' $agff_orig > $agff_fixed
  else
    echo "[$A] GFF chroms already match FASTA — copying as-is"
    cp -f $agff_orig $agff_fixed
  fi

  # gffread --bed
  docker run --rm -v $V3:/v3 -w /v3 $GFFREAD_IMG \
    gffread --bed /v3/inputs/annotations/${A}.fixed.gff3 -o /v3/inputs/annotations/${A}.bed12_raw 2>&1 | tail -2

  # Filter to protein-coding
  awk -F'\t' '!/^#/ && $3=="protein_coding_gene" {n=split($9,a,";"); for(i=1;i<=n;i++) if(match(a[i],/^ID=/)) {print substr(a[i],4); break}}' $agff_fixed | sort -u > /tmp/${A}_pc.txt
  python3 - <<PY
import re
pc=set(open('/tmp/${A}_pc.txt').read().split())
with open('$V3/inputs/annotations/${A}.bed12_raw') as fin, open('$bed12','w') as fout:
    n_in=n_kept=0
    for ln in fin:
        f=ln.rstrip('\n').split('\t')
        n_in+=1
        if len(f)<13: continue
        m=re.search(r'geneID=([^;]+)', f[12])
        if not m: continue
        if m.group(1) in pc:
            fout.write('\t'.join(f[:12])+'\n')
            n_kept+=1
    print(f"  [$A] bed12: kept {n_kept}/{n_in}")
PY
  rm -f $V3/inputs/annotations/${A}.bed12_raw /tmp/${A}_pc.txt

  # Isoforms
  awk -F'\t' '!/^#/ && $3=="mRNA" {
    n=split($9,a,";"); tx=""; gene=""
    for(i=1;i<=n;i++) {
      if (match(a[i],/^ID=/)) tx=substr(a[i],4)
      else if (match(a[i],/^Parent=/)) gene=substr(a[i],8)
    }
    if (tx && gene) print gene "\t" tx
  }' $agff_fixed | sort -k1,1 -k2,2 > $iso

  rm -f $fa_chroms $gff_chroms

  n_bed=$(wc -l < $bed12)
  n_iso=$(wc -l < $iso)
  echo "  [$A] bed12=$n_bed iso=$n_iso"
}

for A in "$@"; do
  build_one $A
done
