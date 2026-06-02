#!/usr/bin/env bash
# Build per-anchor BED12 + isoforms TSV for Phase C (TOGA2 input).
# Adapted from /media/anton/data/sandbox/Pv4/v3/scripts/build_anchor_inputs.sh
set -uo pipefail

if [[ -z "${WORK:-}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  source "${SCRIPT_DIR}/../species.conf"
fi

GFFREAD_IMG="quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4"
WRAPPER="$WORK/pipeline/lib/run_in_${RUN_IN_CONTAINER:-container}.sh"

build_one() {
  local A="$1"
  local agff_fixed="$WORK/inputs/annotations/${A}.fixed.gff3"
  local afa="$WORK/inputs/assemblies/${A}.fa"
  local bed12="$WORK/inputs/annotations/${A}.bed12"
  local iso="$WORK/inputs/annotations/${A}.isoforms.tsv"

  if [[ ! -s "$agff_fixed" ]]; then
    echo "[$A] fixed.gff3 not found — run fix_gff_chroms.sh ${A} first" >&2
    return 1
  fi

  # gffread → raw BED12
  local bed12_raw="$WORK/inputs/annotations/${A}.bed12_raw"
  "$WRAPPER" gffread --bed "$agff_fixed" -o "$bed12_raw" 2>&1 | tail -2

  # Filter to protein-coding genes only
  awk -F'\t' '!/^#/ && ($3=="protein_coding_gene" || $3=="gene") {
    n=split($9,a,";")
    for(i=1;i<=n;i++) {
      if(match(a[i],/^ID=/)) {print substr(a[i],4); break}
    }
  }' "$agff_fixed" | sort -u > /tmp/${A}_pc.txt

  python3 - <<PY
import re
pc = set(open('/tmp/${A}_pc.txt').read().split())
with open('${bed12_raw}') as fin, open('${bed12}', 'w') as fout:
    n_in = n_kept = 0
    for ln in fin:
        f = ln.rstrip('\n').split('\t')
        n_in += 1
        if len(f) < 13:
            continue
        m = re.search(r'geneID=([^;]+)', f[12])
        if not m:
            # Try name field (BED4+)
            m2 = re.match(r'^([^\t]+)', f[3] if len(f) > 3 else '')
            if not m2:
                continue
            gid = m2.group(1)
        else:
            gid = m.group(1)
        if gid in pc or not pc:
            fout.write('\t'.join(f[:12]) + '\n')
            n_kept += 1
    print(f"  [${A}] bed12: kept {n_kept}/{n_in}")
PY
  rm -f "$bed12_raw" /tmp/${A}_pc.txt

  # Isoforms TSV: gene_id <tab> transcript_id
  awk -F'\t' '!/^#/ && $3=="mRNA" {
    n=split($9,a,";"); tx=""; gene=""
    for(i=1;i<=n;i++) {
      if (match(a[i],/^ID=/)) tx=substr(a[i],4)
      else if (match(a[i],/^Parent=/)) gene=substr(a[i],8)
    }
    if (tx && gene) print gene "\t" tx
  }' "$agff_fixed" | sort -k1,1 -k2,2 > "$iso"

  local n_bed
  n_bed=$(wc -l < "$bed12")
  local n_iso
  n_iso=$(wc -l < "$iso")
  echo "  [$A] bed12=$n_bed iso=$n_iso"
}

for A in "$@"; do
  build_one "$A"
done
