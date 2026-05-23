#!/usr/bin/env bash
set -euo pipefail

TOOLS=/media/anton/data/sandbox/Pv4/v3/tools
HUB=/media/anton/data/sandbox/Pv4/v3/ucsc_hub
WORK=/media/anton/data/sandbox/Pv4/v3/work
INPUTS=/media/anton/data/sandbox/Pv4/v3/inputs
STAGING=/media/anton/scratch/Pv4_dropbox_staging
TMPDIR=/media/anton/scratch/maf_tmp
mkdir -p "$TMPDIR"

declare -A ACC
ACC[PvP01]=GCA_900093555.2
ACC[PvSY56]=GCA_003402215.1
ACC[Sal-I]=GCA_000002415.2
ACC[PAM]=GCA_949152365.1
ACC[PvT01]=GCA_900093545.1
ACC[PvC01]=GCA_900093535.1
ACC[MHC087]=GCA_040114635.1
ACC[PvW1]=GCA_914969965.1

HINGES="PvP01 Sal-I PvW1 PAM PvSY56 PvT01 PvC01 MHC087"

log() { echo "[$(date +%H:%M:%S)] $*"; }

# ===== K.1 bigMaf processing =====
log "=== K.1: bigMaf processing ==="
process_maf_hinge() {
    local HINGE=$1
    local ACC_H=$2
    local OUT_DIR="$HUB/$ACC_H"
    local OUT_BB="$OUT_DIR/${HINGE}.multiz.maf.bb"
    [[ -s "$OUT_BB" ]] && { echo "  SKIP $HINGE (exists)"; return 0; }
    
    local MAF_GZ="$STAGING/multiz/${HINGE}/${HINGE}.multiz.maf.gz"
    local TMP_RAW="$TMPDIR/${HINGE}.raw.maf"
    local TMP_FIXED="$TMPDIR/${HINGE}.fixed.maf"
    local SIZES="$STAGING/softmasked/${ACC_H}.fa.fai"
    
    echo "  Processing $HINGE ($ACC_H)..."
    gunzip -c "$MAF_GZ" > "$TMP_RAW"
    python3 "$HUB/process_maf.py" "$ACC_H" "$TMP_RAW" "$TMP_FIXED"
    
    # Write 2-col sizes file
    awk '{print $1"\t"$2}' "$SIZES" > "$TMPDIR/${ACC_H}.sizes"
    
    $TOOLS/mafToBigMaf "$ACC_H" "$TMP_FIXED" /dev/stdout 2>"$OUT_DIR/${HINGE}.maf.log" | \
        sort -k1,1 -k2,2n | \
        $TOOLS/bedToBigBed -type=bed3+1 \
            -as="$HUB/bigMaf.as" -tab \
            stdin "$TMPDIR/${ACC_H}.sizes" \
            "$OUT_BB" 2>>"$OUT_DIR/${HINGE}.maf.log"
    
    if [[ -s "$OUT_BB" ]]; then
        # Build mafIndex
        $TOOLS/mafIndex "$TMP_FIXED" "$OUT_BB.bai" 2>>"$OUT_DIR/${HINGE}.maf.log" || true
        echo "  OK: $HINGE -> $OUT_BB"
    else
        echo "  FAILED: $HINGE"
        cat "$OUT_DIR/${HINGE}.maf.log"
    fi
    rm -f "$TMP_RAW" "$TMP_FIXED"
}
export -f process_maf_hinge
export TOOLS HUB WORK STAGING TMPDIR

# Run 3 concurrent
echo "PvP01 GCA_900093555.2
Sal-I GCA_000002415.2
PvW1 GCA_914969965.1
PAM GCA_949152365.1
PvSY56 GCA_003402215.1
PvT01 GCA_900093545.1
PvC01 GCA_900093535.1
MHC087 GCA_040114635.1" | xargs -P 3 -L 1 bash -c 'process_maf_hinge "$@"' _

log "=== K.1 complete ==="

# ===== K.2 Annotation BigBeds =====
log "=== K.2: Annotation BigBeds ==="

# Extract archives to temp
ANNOT_TMP="$TMPDIR/annot"
mkdir -p "$ANNOT_TMP"

for anchor in PvP01 PvSY56 PvW1 PAM; do
    archive="$WORK/02d_merged/${anchor}-as-ref_archive.tar.gz"
    if [[ ! -d "$ANNOT_TMP/${anchor}-as-ref" ]]; then
        log "  Extracting $anchor archive..."
        tar -xzf "$archive" -C "$ANNOT_TMP/"
    fi
done

build_annot_bb() {
    local anchor=$1
    local strain=$2
    local acc_q=$3
    local OUT="$HUB/$acc_q/annot_from_${anchor}.bb"
    [[ -s "$OUT" ]] && { echo "  SKIP annot_from_${anchor} in $acc_q"; return 0; }
    
    local GFF
    local SIZES="$STAGING/softmasked/${acc_q}.fa.fai"
    
    if [[ "$strain" == "$anchor" ]]; then
        # Self-annotation
        if [[ "$anchor" == "PvP01" ]]; then
            GFF="$INPUTS/annotations/PvP01.genbank.gff3.gz"
        else
            GFF="$INPUTS/annotations/${anchor}.fixed.gff3.gz"
        fi
        if [[ ! -f "$GFF" ]]; then
            echo "  SKIP: missing self-annotation GFF for $anchor"
            return 0
        fi
        local tmp_gff="/tmp/${anchor}_self_$$.gff3"
        gunzip -c "$GFF" > "$tmp_gff"
        GFF="$tmp_gff"
    else
        GFF="$ANNOT_TMP/${anchor}-as-ref/${strain}.annotation.gff3"
        if [[ ! -f "$GFF" ]]; then
            echo "  SKIP: missing $GFF"
            return 0
        fi
    fi
    
    local GP="/tmp/${anchor}_on_${strain}_$$.gp"
    local BED="/tmp/${anchor}_on_${strain}_$$.bed"
    local SORTED_BED="/tmp/${anchor}_on_${strain}_$$.sorted.bed"
    local SIZES2="/tmp/${acc_q}_sizes_$$.txt"
    
    awk '{print $1"\t"$2}' "$SIZES" > "$SIZES2"
    
    $TOOLS/gff3ToGenePred "$GFF" "$GP" 2>/dev/null || { echo "  ERROR gff3ToGenePred for $anchor on $strain"; rm -f "$GP" "$BED" "$SORTED_BED" "$SIZES2" "${GFF%%.gff3}*"; return 1; }
    $TOOLS/genePredToBed "$GP" "$BED" 2>/dev/null || { echo "  ERROR genePredToBed for $anchor on $strain"; rm -f "$GP" "$BED" "$SORTED_BED" "$SIZES2" "${GFF%%.gff3}*"; return 1; }
    sort -k1,1 -k2,2n "$BED" > "$SORTED_BED"
    $TOOLS/bedToBigBed -type=bed12 -tab "$SORTED_BED" "$SIZES2" "$OUT" 2>/tmp/${anchor}_on_${strain}_bb.log || { echo "  ERROR bedToBigBed for $anchor on $strain"; cat "/tmp/${anchor}_on_${strain}_bb.log"; rm -f "$GP" "$BED" "$SORTED_BED" "$SIZES2"; return 1; }
    
    rm -f "$GP" "$BED" "$SORTED_BED" "$SIZES2" "/tmp/${anchor}_on_${strain}_bb.log"
    [[ -n "${tmp_gff:-}" ]] && rm -f "$tmp_gff"
    
    if [[ -s "$OUT" ]]; then
        echo "  OK: annot_from_${anchor}.bb in $acc_q"
    else
        echo "  EMPTY: annot_from_${anchor}.bb in $acc_q"
    fi
}
export -f build_annot_bb
export TOOLS HUB WORK INPUTS STAGING TMPDIR ANNOT_TMP

# Build job list
declare -A A_ACC
A_ACC[PvP01]=GCA_900093555.2
A_ACC[PvSY56]=GCA_003402215.1
A_ACC[Sal-I]=GCA_000002415.2
A_ACC[PAM]=GCA_949152365.1
A_ACC[PvT01]=GCA_900093545.1
A_ACC[PvC01]=GCA_900093535.1
A_ACC[MHC087]=GCA_040114635.1
A_ACC[PvW1]=GCA_914969965.1

{
for anchor in PvP01 PvSY56 PvW1 PAM; do
    for strain in PvP01 PvSY56 Sal-I PAM PvT01 PvC01 MHC087 PvW1; do
        acc_q=${A_ACC[$strain]}
        echo "${anchor} ${strain} ${acc_q}"
    done
done
} | xargs -P 4 -L 1 bash -c 'build_annot_bb "$@"' _

log "=== K.2 complete ==="

# ===== K.4 Orthogroup membership BigBed =====
log "=== K.4: Orthogroup membership BigBed ==="
OUT_OG="$HUB/GCA_900093555.2/orthogroup_membership.bb"
if [[ ! -s "$OUT_OG" ]]; then
    python3 "$HUB/build_orthogroup_bb.py" 2>&1 | tail -5
fi

log "=== K.4 complete ==="

log "=== All done ==="
