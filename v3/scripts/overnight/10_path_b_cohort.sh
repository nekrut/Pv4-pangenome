#!/bin/bash
# Build Path B cohort VCFs by re-labelling MalariaGEN VCF CHROM/POS via
# odgi-position mg_on_*.tsv lookup tables, then concat per target.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p projection/B_graph/cohorts logs/overnight
log() { echo "$(date +%H:%M:%S) [path_b_cohort] $1" | tee -a logs/overnight/STATUS.md; }

OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"
BCF=/home/anton/miniconda3/envs/bcfmod/bin/bcftools
MG_DIR=projection/A1_wfmash/mg_renamed  # already-renamed MG VCFs in LT635xxx CHROMs

build_one() {
  local O=$1
  local TSV=projection/B_graph/sites/mg_on_${O}.tsv
  [ ! -s $TSV ] && { echo "$O: missing $TSV"; return 1; }
  local OUT=projection/B_graph/cohorts/$O
  mkdir -p $OUT
  log "$O: building lookup from $(wc -l < $TSV) odgi positions"

  python3 - <<PY
import csv, gzip, os, sys
from pathlib import Path
import subprocess

# Load lookup: pvp01_chrom_panSN, pvp01_pos → (target_chrom_stripped, target_pos, strand)
lookup = {}
with open('$TSV') as f:
    for ln in f:
        f_ = ln.rstrip().split('\t')
        if len(f_) < 6: continue
        pv_chrom_panSN, pv_start, pv_end, t_start_str, t_end_str, strand = f_
        # PanSN parse
        pv_chrom = pv_chrom_panSN.split('#')[2] if '#' in pv_chrom_panSN else pv_chrom_panSN
        # target positions in "panSN_chrom,pos,strand" format
        t_chrom_panSN, t_pos, t_strand = t_start_str.split(',')
        t_chrom = t_chrom_panSN.split('#')[2] if '#' in t_chrom_panSN else t_chrom_panSN
        lookup[(pv_chrom, int(pv_start))] = (t_chrom, int(t_pos), t_strand)
print(f'$O lookup: {len(lookup)} sites', file=sys.stderr)

# Process each MG per-chr VCF; rewrite CHROM/POS; output per-chr lifted VCF
import glob
mg_dir = 'projection/A1_wfmash/mg_renamed'
out_dir = 'projection/B_graph/cohorts/$O'
os.makedirs(out_dir, exist_ok=True)
RC = lambda s: ''.join({'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}.get(b, 'N') for b in s[::-1])

for vcf_path in sorted(glob.glob(f'{mg_dir}/Pv4_*.vcf.gz')):
    chrom_tag = os.path.basename(vcf_path).replace('Pv4_', '').replace('.vcf.gz', '')
    out_vcf = f'{out_dir}/Pv4_{chrom_tag}_on_$O.vcf'
    if os.path.exists(out_vcf + '.gz'):
        continue
    n_in = 0; n_out = 0
    with gzip.open(vcf_path, 'rt') as fi, open(out_vcf, 'w') as fo:
        for ln in fi:
            if ln.startswith('#'):
                # Don't drop or rewrite contig lines — just copy the file header for now
                if ln.startswith('##contig='):
                    continue
                fo.write(ln)
                continue
            n_in += 1
            f_ = ln.split('\t')
            chrom, pos, _id, ref, alt = f_[0], int(f_[1]), f_[2], f_[3], f_[4]
            key = (chrom, pos - 1)  # BED is 0-based; VCF POS 1-based
            if key not in lookup:
                continue
            t_chrom, t_pos, strand = lookup[key]
            if strand == '-':
                # Reverse complement; assumes biallelic SNV
                ref = RC(ref)
                alt = ','.join(RC(a) for a in alt.split(','))
            f_[0] = t_chrom
            f_[1] = str(t_pos + 1)  # 1-based VCF
            f_[3] = ref
            f_[4] = alt
            fo.write('\t'.join(f_))
            n_out += 1
    print(f'  {chrom_tag}: in={n_in} out={n_out}', file=sys.stderr)
PY

  # bgzip + sort + index per-chr
  log "$O: sorting + indexing per-chr"
  for vcf in $OUT/Pv4_*.vcf; do
    [ ! -s $vcf ] && continue
    [ -s ${vcf}.gz ] && continue
    docker run --rm -v $V3:$V3 -w $V3 \
      quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
      bcftools sort -T /media/anton/scratch -Oz -o ${vcf}.gz $vcf 2>/dev/null && \
      docker run --rm -v $V3:$V3 -w $V3 \
        quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
        tabix -p vcf ${vcf}.gz
    rm -f $vcf
  done

  # Concat
  log "$O: concat"
  FILES=$(ls $OUT/Pv4_*.vcf.gz | sort)
  docker run --rm -v $V3:$V3 -w $V3 \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "bcftools concat -a $FILES -Oz -o $V3/projection/B_graph/Pv4_cohort_on_${O}.vcf.gz && bcftools index $V3/projection/B_graph/Pv4_cohort_on_${O}.vcf.gz" \
    2>>logs/overnight/path_b_concat.log
  if [ -s projection/B_graph/Pv4_cohort_on_${O}.vcf.gz ]; then
    log "$O ✓"
  else
    log "$O FAILED"
  fi
}
export -f build_one
export V3

echo "$OTHERS" | tr ' ' '\n' | xargs -P 4 -I {} bash -c 'build_one {}' 2>&1
log "ALL DONE"
