#!/usr/bin/env python3
"""
Build selection BigBed tracks (strict + relaxed) for PvP01 from HyPhy BUSTED JSONs.
Also builds orthogroup membership BigBed.
"""

import gzip, json, os, subprocess, sys, tarfile, csv
from collections import defaultdict

TOOLS = "/media/anton/data/sandbox/Pv4/v3/tools"
HUB   = "/media/anton/data/sandbox/Pv4/v3/ucsc_hub"
WORK  = "/media/anton/data/sandbox/Pv4/v3/work"
INPUTS = "/media/anton/data/sandbox/Pv4/v3/inputs"
STAGING = "/media/anton/scratch/Pv4_dropbox_staging"

REF_ACC = "GCA_900093555.2"   # PvP01

STRICT_ARCHIVE = f"{WORK}/06_msa/core_v3_hyphy_archive.tar.gz"
RELAXED_ARCHIVE = f"{STAGING}/work_archives/core_relaxed_hyphy_archive.tar.gz"
ORTHOLOG_TABLE  = f"{WORK}/03_consensus/ortholog_table.tsv.gz"
PVP01_BED12     = f"{INPUTS}/annotations/PvP01.bed12.gz"
SIZES_FILE      = f"{STAGING}/softmasked/{REF_ACC}.fa.fai"

AS_FILE = f"{HUB}/bigSelectionPlus5.as"
BEDTOBIGBED = f"{TOOLS}/bedToBigBed"

def load_sizes(fai):
    sizes = {}
    with open(fai) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            sizes[parts[0]] = int(parts[1])
    return sizes

def load_bed12(bedgz):
    """Return dict: gene_id_base -> bed12_fields (list)"""
    bed = {}
    with gzip.open(bedgz, 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 12: continue
            name = f[3]
            # strip .1 isoform suffix
            base = name.rsplit('.', 1)[0] if '.' in name else name
            # keep first occurrence (primary isoform)
            if base not in bed:
                bed[base] = f[:12]
    return bed

def load_ortholog_table(tsv_gz):
    """Return dict: pvp01_gene -> (orthogroup_id, label, n_strains)
    PvP01 column may have multiple genes separated by |; use first one."""
    og_map = {}   # gene_id -> (og_id, label, n_strains)
    og_info = {}  # og_id -> (label, n_strains)
    with gzip.open(tsv_gz, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            og = row['orthogroup_id']
            label = row['label']
            n = int(row['n_strains'])
            og_info[og] = (label, n)
            pvp01_field = row.get('PvP01', '-')
            if pvp01_field == '-' or not pvp01_field:
                continue
            for gene in pvp01_field.split('|'):
                gene = gene.strip()
                # In the table, the gene might include genes from other anchors
                # Keep only PVP01_ genes
                if gene.startswith('PVP01_'):
                    og_map[gene] = (og, label, n)
    return og_map, og_info

def extract_busted_jsons(archive):
    """Return dict: gene_id -> p-value from busted.json"""
    results = {}
    with tarfile.open(archive, 'r:gz') as tf:
        for member in tf.getmembers():
            if not member.name.endswith('busted.json'):
                continue
            # path like core_v3_hyphy/priority/PVP01_XXXXX/busted.json
            parts = member.name.split('/')
            if len(parts) < 2: continue
            gene_id = parts[-2]  # directory name is the gene ID
            fh = tf.extractfile(member)
            if fh is None: continue
            try:
                d = json.load(fh)
                tr = d.get('test results', {})
                pval = tr.get('p-value', None)
                if pval is not None:
                    results[gene_id] = float(pval)
            except Exception as e:
                print(f"  WARN: failed to parse {member.name}: {e}")
    return results

def bh_fdr(pvals_dict):
    """Benjamini-Hochberg FDR. Returns dict: key -> qvalue."""
    items = sorted(pvals_dict.items(), key=lambda x: x[1])
    n = len(items)
    qvals = {}
    prev_q = 1.0
    for i in range(n-1, -1, -1):
        key, p = items[i]
        q = p * n / (i+1)
        q = min(q, prev_q)
        qvals[key] = min(q, 1.0)
        prev_q = q
    return qvals

def qval_to_rgb(q):
    """Color by q-value bin."""
    if q < 0.01:
        return "255,0,0"     # red
    elif q < 0.05:
        return "255,128,0"   # orange
    elif q < 0.10:
        return "200,200,0"   # yellow
    else:
        return "128,128,128" # gray

def qval_to_score(q):
    """0-1000 score, higher is more significant."""
    import math
    if q <= 0:
        return 1000
    try:
        s = int(-math.log10(q) * 100)
    except:
        s = 0
    return max(0, min(1000, s))

def rgb_to_int(rgb_str):
    """Convert '255,0,0' to integer."""
    r, g, b = map(int, rgb_str.split(','))
    return (r << 16) | (g << 8) | b

def is_variant_antigen(gene_id, label):
    """Simple heuristic for variant antigen family."""
    va_labels = {'VIR', 'PIR', 'PHIST', 'DBP', 'RBP', 'SURFIN', 'SERA'}
    for va in va_labels:
        if va.lower() in label.lower():
            return label
    return 'other'

def build_selection_bed(busted_results, qvals, og_map, bed12, sizes, out_bed):
    """Write BED12+5 for selection track."""
    written = 0
    skipped_no_og = 0
    skipped_no_bed = 0
    skipped_no_chrom = 0

    lines = []
    for gene_id, pval in busted_results.items():
        if gene_id not in og_map:
            skipped_no_og += 1
            continue
        og_id, label, n_strains = og_map[gene_id]
        if gene_id not in bed12:
            # Try without the anchor prefix copies
            skipped_no_bed += 1
            continue
        b = bed12[gene_id]
        chrom = b[0]
        if chrom not in sizes:
            skipped_no_chrom += 1
            continue
        qval = qvals.get(gene_id, 1.0)
        rgb = qval_to_rgb(qval)
        score = qval_to_score(qval)
        rgb_int = rgb_to_int(rgb)
        gene_family = is_variant_antigen(gene_id, label)
        # BED12+5: chrom, start, end, name, score, strand, thickStart, thickEnd,
        #           itemRgb, blockCount, blockSizes, blockStarts,
        #           orthogroup_id, n_strains, busted_pvalue, busted_qvalue_fdr, gene_family
        row = '\t'.join([
            b[0], b[1], b[2], og_id, str(score), b[5],
            b[6], b[7], str(rgb_int),
            b[9], b[10], b[11],
            og_id, str(n_strains),
            f"{pval:.6g}", f"{qval:.6g}", gene_family
        ])
        lines.append((b[0], int(b[1]), row))

    # Sort
    lines.sort(key=lambda x: (x[0], x[1]))
    with open(out_bed, 'w') as fh:
        for chrom, start, row in lines:
            fh.write(row + '\n')
            written += 1

    print(f"  Written: {written}, skipped_no_og: {skipped_no_og}, "
          f"skipped_no_bed: {skipped_no_bed}, skipped_no_chrom: {skipped_no_chrom}")
    return written

def bed_to_bigbed(bed_file, sizes_fai, as_file, out_bb, bed_type="bed12+5"):
    """Run bedToBigBed."""
    # Write sizes file in 2-col format
    sizes_2col = out_bb + ".tmp.sizes"
    with open(sizes_fai) as fh, open(sizes_2col, 'w') as out:
        for line in fh:
            parts = line.strip().split('\t')
            out.write(f"{parts[0]}\t{parts[1]}\n")
    cmd = [BEDTOBIGBED, f"-type={bed_type}", f"-as={as_file}", "-tab",
           bed_file, sizes_2col, out_bb]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(sizes_2col)
    if result.returncode != 0:
        print(f"  ERROR bedToBigBed: {result.stderr[:500]}")
        return False
    print(f"  OK: {out_bb}")
    return True

def build_orthogroup_bb(og_map, bed12, sizes, out_dir, sizes_fai):
    """Build orthogroup membership BED12 for PvP01 and convert to BigBed."""
    print("\n=== Building orthogroup membership BigBed ===")
    out_bed = "/tmp/orthogroup_membership.bed"

    # Reverse map: og_id -> list of gene_ids
    og_to_genes = defaultdict(list)
    for gene_id, (og_id, label, n_strains) in og_map.items():
        og_to_genes[og_id].append((gene_id, n_strains))

    lines = []
    for og_id, genes in og_to_genes.items():
        for gene_id, n_strains in genes:
            if gene_id not in bed12:
                continue
            b = bed12[gene_id]
            chrom = b[0]
            if chrom not in sizes:
                continue
            # Color by n_strains: 1=red, 8=green
            r = max(0, int(255 * (1 - (n_strains - 1) / 7)))
            g = max(0, int(255 * ((n_strains - 1) / 7)))
            rgb_int = (r << 16) | g
            score = int(n_strains * 125)  # 125-1000 range
            row = '\t'.join([
                b[0], b[1], b[2], og_id, str(score), b[5],
                b[6], b[7], str(rgb_int),
                b[9], b[10], b[11]
            ])
            lines.append((b[0], int(b[1]), row))

    lines.sort(key=lambda x: (x[0], x[1]))
    with open(out_bed, 'w') as fh:
        for chrom, start, row in lines:
            fh.write(row + '\n')
    print(f"  Written {len(lines)} rows")

    out_bb = f"{out_dir}/orthogroup_membership.bb"
    return bed_to_bigbed(out_bed, sizes_fai, None, out_bb, bed_type="bed12")

def bed_to_bigbed_plain(bed_file, sizes_fai, out_bb):
    """Run bedToBigBed for plain BED12 (no AS file)."""
    sizes_2col = out_bb + ".tmp.sizes"
    with open(sizes_fai) as fh, open(sizes_2col, 'w') as out:
        for line in fh:
            parts = line.strip().split('\t')
            out.write(f"{parts[0]}\t{parts[1]}\n")
    cmd = [BEDTOBIGBED, "-type=bed12", bed_file, sizes_2col, out_bb]
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    os.unlink(sizes_2col)
    if result.returncode != 0:
        print(f"  ERROR bedToBigBed: {result.stderr[:500]}")
        return False
    print(f"  OK: {out_bb}")
    return True


def main():
    out_dir = f"{HUB}/{REF_ACC}"
    os.makedirs(out_dir, exist_ok=True)

    print("Loading genome sizes...")
    sizes = load_sizes(SIZES_FILE)

    print("Loading PvP01 BED12...")
    bed12 = load_bed12(PVP01_BED12)
    print(f"  Loaded {len(bed12)} genes")

    print("Loading ortholog table...")
    og_map, og_info = load_ortholog_table(ORTHOLOG_TABLE)
    print(f"  Loaded {len(og_map)} gene->OG mappings, {len(og_info)} OGs")

    # --- Strict BUSTED ---
    print("\n=== Strict BUSTED (core_v3_hyphy) ===")
    print("Extracting BUSTED JSONs (strict)...")
    strict_results = extract_busted_jsons(STRICT_ARCHIVE)
    print(f"  Found {len(strict_results)} gene results")

    # Compute FDR only on genes that have OG mapping
    strict_for_fdr = {g: p for g, p in strict_results.items() if g in og_map}
    strict_qvals = bh_fdr(strict_for_fdr)
    n_sig_01 = sum(1 for q in strict_qvals.values() if q < 0.01)
    n_sig_05 = sum(1 for q in strict_qvals.values() if q < 0.05)
    print(f"  Significant at q<0.01: {n_sig_01}, q<0.05: {n_sig_05}")

    strict_bed = "/tmp/selection_strict.bed"
    print(f"Building strict selection BED12+5...")
    build_selection_bed(strict_results, strict_qvals, og_map, bed12, sizes, strict_bed)

    print("Converting to BigBed...")
    bed_to_bigbed(strict_bed, SIZES_FILE, AS_FILE,
                  f"{out_dir}/selection_strict.bb")

    # --- Relaxed BUSTED ---
    print("\n=== Relaxed BUSTED (core_relaxed_hyphy) ===")
    print("Extracting BUSTED JSONs (relaxed)...")
    relaxed_results = extract_busted_jsons(RELAXED_ARCHIVE)
    print(f"  Found {len(relaxed_results)} gene results")

    relaxed_for_fdr = {g: p for g, p in relaxed_results.items() if g in og_map}
    relaxed_qvals = bh_fdr(relaxed_for_fdr)
    n_sig_01 = sum(1 for q in relaxed_qvals.values() if q < 0.01)
    n_sig_05 = sum(1 for q in relaxed_qvals.values() if q < 0.05)
    print(f"  Significant at q<0.01: {n_sig_01}, q<0.05: {n_sig_05}")

    relaxed_bed = "/tmp/selection_relaxed.bed"
    print(f"Building relaxed selection BED12+5...")
    build_selection_bed(relaxed_results, relaxed_qvals, og_map, bed12, sizes, relaxed_bed)

    print("Converting to BigBed...")
    bed_to_bigbed(relaxed_bed, SIZES_FILE, AS_FILE,
                  f"{out_dir}/selection_relaxed.bb")

    # --- Orthogroup membership ---
    print("\n=== Orthogroup membership BigBed ===")
    lines = []
    for gene_id, (og_id, label, n_strains) in og_map.items():
        if gene_id not in bed12:
            continue
        b = bed12[gene_id]
        chrom = b[0]
        if chrom not in sizes:
            continue
        r = max(0, int(255 * (1 - (n_strains - 1) / 7)))
        g_c = max(0, int(255 * ((n_strains - 1) / 7)))
        rgb_int = (r << 16) | g_c
        score = int(n_strains * 125)
        row = '\t'.join([
            b[0], b[1], b[2], og_id, str(score), b[5],
            b[6], b[7], str(rgb_int),
            b[9], b[10], b[11]
        ])
        lines.append((b[0], int(b[1]), row))

    lines.sort(key=lambda x: (x[0], x[1]))
    og_bed = "/tmp/orthogroup_membership.bed"
    with open(og_bed, 'w') as fh:
        for chrom, start, row in lines:
            fh.write(row + '\n')
    print(f"  Written {len(lines)} rows to {og_bed}")

    out_bb = f"{out_dir}/orthogroup_membership.bb"
    bed_to_bigbed_plain(og_bed, SIZES_FILE, out_bb)

    print("\n=== Done ===")
    for f in [f"{out_dir}/selection_strict.bb",
              f"{out_dir}/selection_relaxed.bb",
              f"{out_dir}/orthogroup_membership.bb"]:
        if os.path.exists(f):
            sz = os.path.getsize(f)
            print(f"  {f}: {sz:,} bytes")
        else:
            print(f"  MISSING: {f}")

if __name__ == '__main__':
    main()
