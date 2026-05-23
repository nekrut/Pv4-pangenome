#!/bin/bash
# Path B cohort VCF build — safe version. Single-target-at-a-time, atomic per-chr
# steps (verify .vcf.gz before deleting .vcf), proper error handling.
set -u
V3=/media/anton/data/sandbox/Pv4/v3
cd $V3
mkdir -p projection/B_graph/cohorts logs/overnight
log() { echo "$(date +%H:%M:%S) [path_b_safe] $1" | tee -a logs/overnight/STATUS.md; }
log "starting"

OTHERS="GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1"
SCRATCH=/media/anton/scratch
SORT_TMP=$SCRATCH/path_b_sort_tmp
mkdir -p $SORT_TMP

build_one() {
  local O=$1
  local OUT=projection/B_graph/cohorts/$O
  mkdir -p $OUT
  local TSV=projection/B_graph/sites/mg_on_${O}.tsv
  if [ ! -s $TSV ]; then log "$O: missing TSV"; return 1; fi

  log "$O START"

  # Python builds all 16 per-chr .vcf files (uncompressed) into $OUT
  python3 - "$O" "$TSV" "$OUT" <<'PYEOF'
import csv, gzip, os, sys, glob
O = sys.argv[1]; TSV = sys.argv[2]; OUT = sys.argv[3]

lookup = {}
with open(TSV) as f:
    for ln in f:
        f_ = ln.rstrip().split('\t')
        if len(f_) < 6: continue
        pv_chrom_panSN, pv_start, pv_end, t_str, t_end_str, strand = f_
        pv_chrom = pv_chrom_panSN.split('#')[2] if '#' in pv_chrom_panSN else pv_chrom_panSN
        t_chrom_panSN, t_pos, t_strand = t_str.split(',')
        t_chrom = t_chrom_panSN.split('#')[2] if '#' in t_chrom_panSN else t_chrom_panSN
        lookup[(pv_chrom, int(pv_start))] = (t_chrom, int(t_pos), t_strand)
print(f'  {O}: lookup loaded ({len(lookup)} sites)', file=sys.stderr, flush=True)

RC = lambda s: ''.join({'A':'T','T':'A','G':'C','C':'G','N':'N'}.get(b.upper(), 'N') for b in s[::-1])

for vcf_path in sorted(glob.glob('projection/A1_wfmash/mg_renamed/Pv4_*.vcf.gz')):
    chrom_tag = os.path.basename(vcf_path).replace('Pv4_', '').replace('.vcf.gz', '')
    out_vcf = f'{OUT}/Pv4_{chrom_tag}_on_{O}.vcf'
    out_gz = out_vcf + '.gz'
    if os.path.exists(out_gz):
        continue
    if os.path.exists(out_vcf):
        os.remove(out_vcf)
    n_in = 0; n_out = 0
    with gzip.open(vcf_path, 'rt') as fi, open(out_vcf, 'w') as fo:
        for ln in fi:
            if ln.startswith('#'):
                if ln.startswith('##contig='):
                    continue
                fo.write(ln); continue
            n_in += 1
            ff = ln.split('\t')
            chrom, pos = ff[0], int(ff[1]); ref = ff[3]; alt = ff[4]
            key = (chrom, pos - 1)
            if key not in lookup: continue
            t_chrom, t_pos, strand = lookup[key]
            if strand == '-':
                ref = RC(ref)
                alt = ','.join(RC(a) for a in alt.split(','))
            ff[0] = t_chrom
            ff[1] = str(t_pos + 1)
            ff[3] = ref
            ff[4] = alt
            fo.write('\t'.join(ff))
            n_out += 1
    print(f'  {O} {chrom_tag}: in={n_in} out={n_out}', file=sys.stderr, flush=True)
PYEOF
  log "$O python done"

  # Sort + bgzip + index per-chr — ATOMIC: only delete .vcf after .vcf.gz exists
  local n_done=0; local n_fail=0
  for vcf in $OUT/Pv4_*.vcf; do
    [ ! -s $vcf ] && continue
    local gz=${vcf}.gz
    [ -s $gz ] && { rm -f $vcf; n_done=$((n_done+1)); continue; }
    local tmp_gz=${vcf}.tmp.gz
    if docker run --rm -v $V3:$V3 -v $SCRATCH:$SCRATCH -w $V3 \
        quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
        bcftools sort -T $SORT_TMP -Oz -o $tmp_gz $vcf 2>>logs/overnight/path_b_${O}.log; then
      mv $tmp_gz $gz
      docker run --rm -v $V3:$V3 -w $V3 \
        quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
        tabix -p vcf $gz 2>>logs/overnight/path_b_${O}.log
      if [ -s $gz ]; then
        rm -f $vcf  # delete uncompressed only after .gz confirmed
        n_done=$((n_done+1))
      else
        n_fail=$((n_fail+1))
        log "$O: sort succeeded but $gz empty"
      fi
    else
      n_fail=$((n_fail+1))
      rm -f $tmp_gz
      log "$O: sort FAILED on $(basename $vcf)"
    fi
  done
  log "$O sort: $n_done done, $n_fail failed"

  # Concat
  local FILES=""
  for f in $OUT/Pv4_*.vcf.gz; do
    [ -s $f ] && FILES="$FILES $f"
  done
  if [ -z "$FILES" ]; then log "$O: no inputs for concat"; return 1; fi
  local COHORT=$V3/projection/B_graph/Pv4_cohort_on_${O}.vcf.gz
  docker run --rm -v $V3:$V3 -w $V3 \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash -c "bcftools concat -a $FILES -Oz -o $COHORT && bcftools index $COHORT" \
    2>>logs/overnight/path_b_${O}.log
  if [ -s $COHORT ]; then
    log "$O ✓ cohort $(du -h $COHORT | cut -f1)"
  else
    log "$O concat FAILED"
    return 1
  fi
}
export -f build_one
export V3 SCRATCH SORT_TMP

# Process sequentially (one at a time) to avoid the earlier docker chaos
for O in $OTHERS; do
  build_one $O
done

log "ALL DONE"
