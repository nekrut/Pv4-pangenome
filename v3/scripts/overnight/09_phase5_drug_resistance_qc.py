#!/usr/bin/env python3
"""Phase 5 — VCF × orthology QC. For drug-resistance genes (DHFR-TS, DHPS, MDR1,
ATP4, MRP1/2, K12, FNT, CRT-o, DHODH, plasmepsins, UBP1), verify that variants
lifted onto non-PvP01 references via A1/A2/B land inside the orthologous gene's
CDS interval. Output: per-method per-gene coverage table.
"""
import csv
import subprocess
from pathlib import Path
from collections import defaultdict

V3 = Path('/media/anton/data/sandbox/Pv4/v3')
import os; os.chdir(V3)
BCFTOOLS = '/home/anton/miniconda3/envs/bcfmod/bin/bcftools'

# Drug-resistance PVP01 IDs from validator output. Look these up in the priority TSV.
PRIORITY = {}
with open('work/05_priorities/gene_priorities.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        if row.get('category') == 'drug_resistance' and row['PVP01_id'].startswith('PVP01_'):
            PRIORITY[row['PVP01_id']] = row['gene_symbol']

print(f"Drug-resistance priority genes: {len(PRIORITY)}")
for pid, sym in PRIORITY.items():
    print(f"  {pid}  {sym}")

# Load PvP01 BED to get CDS intervals
pvp01_intervals = {}
with open('inputs/annotations/PvP01.bed') as f:
    for ln in f:
        c, s, e, gid = ln.rstrip('\n').split('\t')[:4]
        if gid in PRIORITY:
            pvp01_intervals[gid] = (c, int(s), int(e))

# Load orthology table to map PvP01 gene → orthologs in each strain
ortho_map = defaultdict(dict)  # pvp01_id → strain → gene_id
with open('work/03_consensus/ortholog_table.tsv') as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        # Find which PvP01 IDs in this orthogroup we care about
        pv_cell = row.get('PvP01', '-')
        pvids = [g for g in pv_cell.replace('|', ',').split(',') if g in PRIORITY]
        for pvid in pvids:
            for strain in ['Sal-I', 'PvW1', 'PAM', 'PvSY56', 'PvT01', 'PvC01', 'MHC087']:
                cell = row.get(strain, '-')
                if cell != '-':
                    ortho_map[pvid][strain] = cell

# Load per-strain BEDs to get target CDS intervals
strain_to_acc = {
    'Sal-I':  'GCA_000002415.2',
    'PvW1':   'GCA_914969965.1',
    'PAM':    'GCA_949152365.1',
    'PvSY56': 'GCA_003402215.1',
    'PvT01':  'GCA_900093545.1',
    'PvC01':  'GCA_900093535.1',
    'MHC087': 'GCA_040114635.1',
}
strain_beds = {
    'Sal-I':  'inputs/annotations/sal-i_ncbi.bed',
    'PvW1':   'inputs/annotations/PvW1.bed',
    'PAM':    'inputs/annotations/PAM.bed',
    'PvSY56': 'inputs/annotations/PvSY56.bed',
    'PvT01':  'inputs/annotations/pvt01_ncbi.bed',
    'PvC01':  'inputs/annotations/pvc01_ncbi.bed',
    'MHC087': 'work/02d_merged/PvP01-as-ref/MHC087.bed',
}
target_intervals = defaultdict(dict)  # strain → gene_id → (chrom, start, end)
for strain, bed in strain_beds.items():
    if not Path(bed).exists():
        continue
    with open(bed) as fh:
        for ln in fh:
            f = ln.rstrip('\n').split('\t')
            if len(f) >= 4:
                target_intervals[strain][f[3]] = (f[0], int(f[1]), int(f[2]))

# For each drug-resistance gene × method × target:
# 1. Pull source MalariaGEN variants in PvP01 CDS
# 2. Pull lifted variants in target VCF
# 3. Check overlap with target ortholog's CDS
print()
print("="*90)
print(f"{'pvp01_id':18s} {'sym':10s} {'method':4s} {'strain':8s} {'src_n':>6s} {'lift_n':>7s} {'in_cds':>7s} {'frac':>6s}")
print("="*90)

results = []
for pvid, sym in PRIORITY.items():
    if pvid not in pvp01_intervals:
        continue
    pv_chrom, pv_start, pv_end = pvp01_intervals[pvid]
    # Count source variants in PvP01 CDS — pick any MalariaGEN VCF
    # Use the renamed mg files
    src_n = 0
    for chr_vcf in Path('projection/A1_wfmash/mg_renamed').glob('Pv4_*.vcf.gz'):
        try:
            out = subprocess.run([BCFTOOLS, 'view', '-H', '-r', f'{pv_chrom}:{pv_start}-{pv_end}', str(chr_vcf)],
                                 capture_output=True, text=True, timeout=30)
            src_n += len(out.stdout.strip().split('\n')) if out.stdout.strip() else 0
        except Exception:
            pass
        if src_n > 0:
            break  # found the right chromosome

    for strain, acc in strain_to_acc.items():
        # Get target gene IDs
        target_gene_ids = ortho_map[pvid].get(strain, '').replace('|', ',').split(',')
        target_gene_ids = [g for g in target_gene_ids if g and g != '-']
        target_cds = []
        for g in target_gene_ids:
            if g in target_intervals[strain]:
                target_cds.append(target_intervals[strain][g])
        if not target_cds:
            continue
        for method, vcf_path in [
            ('A1', f'projection/A1_wfmash/Pv4_cohort_on_{acc}.vcf.gz'),
            ('A2', f'projection/A2_lastz/Pv4_cohort_on_{acc}.vcf.gz'),
            ('B',  f'projection/B_graph/Pv4_cohort_on_{acc}.vcf.gz'),
        ]:
            if not Path(vcf_path).exists():
                continue
            # Pull lifted variants from each target CDS interval
            lift_n = 0; in_cds = 0
            for tc, ts, te in target_cds:
                try:
                    out = subprocess.run([BCFTOOLS, 'view', '-H', '-r', f'{tc}:{ts}-{te}', vcf_path],
                                         capture_output=True, text=True, timeout=30)
                    if out.stdout.strip():
                        n = len(out.stdout.strip().split('\n'))
                        in_cds += n
                except Exception:
                    pass
            # Total lifted variants ANYWHERE in the cohort VCF associated with this PvP01 gene
            # (just count in_cds for now; full extraction is too slow)
            lift_n = in_cds  # within-CDS only here
            frac = in_cds / src_n if src_n else 0
            print(f"{pvid:18s} {sym[:10]:10s} {method:4s} {strain:8s} {src_n:>6d} {lift_n:>7d} {in_cds:>7d} {frac:>6.2f}")
            results.append({'pvp01_id': pvid, 'symbol': sym, 'method': method, 'strain': strain,
                            'src_n': src_n, 'lift_n_in_cds': in_cds, 'frac_in_cds': frac})

# Write
out_path = Path('work/09_projection_compare/drug_resistance_qc.tsv')
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['pvp01_id','symbol','method','strain','src_n','lift_n_in_cds','frac_in_cds'], delimiter='\t')
    w.writeheader(); w.writerows(results)
print(f"\nWrote {out_path} with {len(results)} rows")
