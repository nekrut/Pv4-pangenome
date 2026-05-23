#!/usr/bin/env python3
"""Build orthogroup membership BigBed for PvP01."""
import gzip, os, subprocess, sys, csv

TOOLS = "/media/anton/data/sandbox/Pv4/v3/tools"
HUB   = "/media/anton/data/sandbox/Pv4/v3/ucsc_hub"
WORK  = "/media/anton/data/sandbox/Pv4/v3/work"
INPUTS = "/media/anton/data/sandbox/Pv4/v3/inputs"
STAGING = "/media/anton/scratch/Pv4_dropbox_staging"
REF_ACC = "GCA_900093555.2"

def load_bed12(bedgz):
    bed = {}
    with gzip.open(bedgz, 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 12: continue
            name = f[3]
            base = name.rsplit('.', 1)[0] if '.' in name else name
            if base not in bed:
                bed[base] = f[:12]
    return bed

def load_ortholog_table(tsv_gz):
    og_map = {}
    with gzip.open(tsv_gz, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            og = row['orthogroup_id']
            n = int(row['n_strains'])
            pvp01_field = row.get('PvP01', '-')
            if pvp01_field == '-' or not pvp01_field:
                continue
            for gene in pvp01_field.split('|'):
                gene = gene.strip()
                if gene.startswith('PVP01_'):
                    og_map[gene] = (og, n)
    return og_map

out_dir = f"{HUB}/{REF_ACC}"
os.makedirs(out_dir, exist_ok=True)
out_bb = f"{out_dir}/orthogroup_membership.bb"
if os.path.exists(out_bb) and os.path.getsize(out_bb) > 0:
    print(f"Exists: {out_bb}")
    sys.exit(0)

print("Loading BED12...")
bed12 = load_bed12(f"{INPUTS}/annotations/PvP01.bed12.gz")
print(f"  {len(bed12)} genes")

print("Loading ortholog table...")
og_map = load_ortholog_table(f"{WORK}/03_consensus/ortholog_table.tsv.gz")
print(f"  {len(og_map)} gene->OG mappings")

print("Loading sizes...")
sizes = {}
with open(f"{STAGING}/softmasked/{REF_ACC}.fa.fai") as fh:
    for line in fh:
        parts = line.strip().split('\t')
        sizes[parts[0]] = int(parts[1])

lines = []
clipped = 0
for gene_id, (og_id, n_strains) in og_map.items():
    if gene_id not in bed12:
        continue
    b = bed12[gene_id]
    chrom = b[0]
    if chrom not in sizes:
        continue
    # Check if annotation exceeds chrom size (including blocks) - skip if so
    start = int(b[1])
    end = int(b[2])
    thick_start = int(b[6])
    thick_end = int(b[7])
    chrom_size = sizes[chrom]
    block_sizes = [int(x) for x in b[10].rstrip(',').split(',') if x]
    block_starts_rel = [int(x) for x in b[11].rstrip(',').split(',') if x]
    if block_sizes and block_starts_rel:
        true_end = start + block_starts_rel[-1] + block_sizes[-1]
    else:
        true_end = end
    if true_end > chrom_size:
        clipped += 1
        continue  # skip genes that don't fit
    if end > chrom_size or start >= end:
        clipped += 1
        continue
    r = max(0, int(255 * (1 - (n_strains - 1) / 7)))
    g = max(0, int(255 * ((n_strains - 1) / 7)))
    rgb_int = (r << 16) | g
    score = int(n_strains * 125)
    row = '\t'.join([chrom, str(start), str(end), og_id, str(score), b[5],
                     b[6], b[7], str(rgb_int), b[9], b[10], b[11]])
    lines.append((chrom, start, row))

print(f"  Clipped {clipped} entries to chrom size")
lines.sort(key=lambda x: (x[0], x[1]))
out_bed = "/tmp/og_mem.bed"
with open(out_bed, 'w') as fh:
    for c, s, row in lines:
        fh.write(row + '\n')
print(f"  Written {len(lines)} rows to {out_bed}")

sizes_2col = "/tmp/og_sizes.txt"
with open(f"{STAGING}/softmasked/{REF_ACC}.fa.fai") as fh, open(sizes_2col, 'w') as out:
    for line in fh:
        p = line.strip().split('\t')
        out.write(f"{p[0]}\t{p[1]}\n")

cmd = [f"{TOOLS}/bedToBigBed", "-type=bed12", out_bed, sizes_2col, out_bb]
print(f"Running: {' '.join(cmd)}")
r = subprocess.run(cmd, capture_output=True, text=True)
if r.returncode != 0:
    print(f"ERROR: {r.stderr[:500]}")
    sys.exit(1)
os.unlink(out_bed)
os.unlink(sizes_2col)
print(f"OK: {out_bb} ({os.path.getsize(out_bb):,} bytes)")
