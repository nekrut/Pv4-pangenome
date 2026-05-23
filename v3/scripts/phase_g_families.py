#!/usr/bin/env python3
"""Phase G — variant-antigen / gene family annotation.

Skip HMMER. Use PlasmoDB description= field (already family-aware) plus
orthogroup propagation for non-PlasmoDB strains.

Output: work/05_families/family_table.tsv
  gene_id  strain  family  evidence  description
"""
import csv
import re
from pathlib import Path
from collections import defaultdict, Counter

OUT_DIR = Path('work/05_families')
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Strain → PlasmoDB GFF (curated, has descriptions)
PLASMODB_GFFS = {
    'PvP01':  'inputs/annotations/plasmodb-68/PvP01.gff3',
    'Sal-I':  'inputs/annotations/plasmodb-68/Sal-I.gff3',
    'PvW1':   'inputs/annotations/plasmodb-68/PvW1.gff3',
    'PAM':    'inputs/annotations/plasmodb-68/PAM.gff3',
    # PvSY56 was Liftoff-projected — no native PlasmoDB family labels
}

# All strain names
ALL_STRAINS = ['PvP01', 'Sal-I', 'PvW1', 'PAM', 'PvSY56', 'PvT01', 'PvC01', 'MHC087']

# Family classifier — pattern → family label
FAMILY_PATTERNS = [
    # vir/PIR family — many synonyms across PlasmoDB strains
    (re.compile(r'\bPIR protein\b|\bVIR protein\b|\bvir gene\b|\bvir family\b|'
                r'\bvariable surface protein Vir\d*\b|\bvariant interspersed repeat\b', re.I), 'PIR'),
    # PHIST family
    (re.compile(r'\bPHIST\b|\bPhist protein\b|\bPlasmodium helical interspersed subtelomeric\b|'
                r'\bPf-fam-b\b', re.I), 'PHIST'),
    # Pv-fam variants (-a tryptophan-rich, -d, -e, -h)
    (re.compile(r'\bPv-fam-h\b|\bPv-fam-d\b|\bPv-fam-e\b|\bRAD protein\b', re.I), 'Pv-fam'),
    # Other exported/subtelomeric
    (re.compile(r'\bSTP1\b', re.I), 'STP1'),
    (re.compile(r'\bRESA\b|\bring-infected erythrocyte surface antigen\b', re.I), 'RESA'),
    (re.compile(r'\bMSP\b|\bmerozoite surface protein\b', re.I), 'MSP'),
    (re.compile(r'\bRAP\b|\brhoptry-associated protein\b', re.I), 'RAP'),
    (re.compile(r'\bRBP\b|\breticulocyte binding\b', re.I), 'RBP'),
    (re.compile(r'\bAMA[-_]?1\b|\bapical membrane antigen\b', re.I), 'AMA'),
    (re.compile(r'\bEBA\b|\berythrocyte binding antigen\b', re.I), 'EBA'),
    (re.compile(r'\bDBP\b|\bDuffy binding protein\b', re.I), 'DBP'),
    (re.compile(r'\bSERA\b|\bserine repeat antigen\b', re.I), 'SERA'),
    (re.compile(r'\btryptophan[ -]rich antigen\b|\bTRAg\b', re.I), 'TRAg'),
    # tRNA / ncRNA (not antigen families but useful tag)
    (re.compile(r'\btRNA\b', re.I), 'tRNA'),
    (re.compile(r'\brRNA\b|\bribosomal RNA\b', re.I), 'rRNA'),
    (re.compile(r'\bsnoRNA\b|\bsnRNA\b|\bncRNA\b', re.I), 'ncRNA'),
    # Conserved / hypothetical / unknown
    (re.compile(r'\bconserved Plasmodium protein\b|\bconserved unknown function\b', re.I), 'conserved'),
    (re.compile(r'\bunspecified product\b|\bhypothetical protein\b', re.I), 'hypothetical'),
]


def classify(description):
    if not description:
        return 'unannotated'
    for pat, label in FAMILY_PATTERNS:
        if pat.search(description):
            return label
    return 'other'


def parse_gff(path):
    """Yield (gene_id, attributes_dict) for protein_coding_gene + ncRNA_gene + pseudogene."""
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            feat = f[2]
            if feat not in ('protein_coding_gene', 'ncRNA_gene', 'pseudogene', 'gene'):
                continue
            attrs = {}
            for kv in f[8].split(';'):
                kv = kv.strip()
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attrs[k] = v
            gid = attrs.get('ID')
            if gid:
                # Strip 'gene:' / 'transcript:' prefix if any
                gid = gid.split(':', 1)[-1]
                yield gid, attrs


def main():
    rows = []
    family_counts = defaultdict(Counter)

    # Pass 1 — PlasmoDB-annotated strains
    print("Pass 1: PlasmoDB-annotated strains")
    for strain, gff in PLASMODB_GFFS.items():
        if not Path(gff).exists():
            print(f"  skip {strain}: {gff} missing")
            continue
        n = 0
        for gid, attrs in parse_gff(gff):
            desc = attrs.get('description', '')
            from urllib.parse import unquote
            desc = unquote(desc).replace('+', ' ')
            fam = classify(desc)
            rows.append({
                'gene_id': gid,
                'strain': strain,
                'family': fam,
                'evidence': 'plasmodb-description',
                'description': desc[:120],
            })
            family_counts[strain][fam] += 1
            n += 1
        print(f"  {strain:8s}  {n} genes labeled")

    # Pass 2 — orthogroup propagation: if a gene is labeled 'other'/'unannotated'
    # in one strain but the orthogroup has a real family label from another strain,
    # propagate. Also handles non-PlasmoDB strains.
    print("Pass 2: orthogroup propagation across all strains (lift 'other' to family)")
    ortho_path = Path('work/03_consensus/ortholog_table.tsv')
    if not ortho_path.exists():
        print(f"  missing {ortho_path} — skipping pass 2")
    else:
        # Build gene_id → family lookup from pass 1
        gene2fam = {}
        for r in rows:
            gene2fam[r['gene_id']] = r['family']

        # Build the (strain, gene_id) → (family, row_index) lookup
        sg2row = {(r['strain'], r['gene_id']): i for i, r in enumerate(rows)}
        non_plasmodb = ['PvSY56', 'PvT01', 'PvC01', 'MHC087']
        n_propagated = 0
        n_lifted = 0
        with open(ortho_path) as fh:
            r = csv.DictReader(fh, delimiter='\t')
            for og in r:
                # Collect family labels per strain in this OG
                labels = []
                for s in ALL_STRAINS:
                    cell = og.get(s, '-')
                    if cell == '-':
                        continue
                    for gid in re.split(r'[,|]', cell):
                        if (s, gid) in sg2row:
                            f = rows[sg2row[(s, gid)]]['family']
                            if f not in ('other', 'unannotated', 'hypothetical', 'conserved', 'tRNA', 'rRNA', 'ncRNA'):
                                labels.append(f)
                if not labels:
                    continue
                cnt = Counter(labels)
                modal = cnt.most_common(1)[0][0]
                # Apply: lift any 'other/unannotated/hypothetical/conserved' member to modal,
                # add rows for non-PlasmoDB strains that have no entry yet
                for s in ALL_STRAINS:
                    cell = og.get(s, '-')
                    if cell == '-':
                        continue
                    for gid in re.split(r'[,|]', cell):
                        if (s, gid) in sg2row:
                            r_existing = rows[sg2row[(s, gid)]]
                            if r_existing['family'] in ('other', 'unannotated', 'hypothetical', 'conserved'):
                                r_existing['family'] = modal
                                r_existing['evidence'] = r_existing['evidence'] + '+orthogroup-lifted'
                                family_counts[s][modal] += 1
                                family_counts[s][r_existing['family']] -= 1
                                n_lifted += 1
                        elif s in non_plasmodb:
                            rows.append({
                                'gene_id': gid,
                                'strain': s,
                                'family': modal,
                                'evidence': 'orthogroup-propagated',
                                'description': f'propagated from OG={og["orthogroup_id"]}',
                            })
                            sg2row[(s, gid)] = len(rows) - 1
                            family_counts[s][modal] += 1
                            n_propagated += 1
        print(f"  propagated {n_propagated} new labels; lifted {n_lifted} existing 'other'/'hypothetical' labels")

    # Pass 3 — for non-PlasmoDB strains, fill remaining unlabeled genes from their
    # annotation GFF (Liftoff or NCBI). Use whatever description they carry.
    print("Pass 3: residual fill from per-strain annotation GFFs")
    fallback_gffs = {
        'PvSY56': 'work/02d_merged/PvP01-as-ref/PvSY56.annotation.gff3',
        'PvT01':  'work/02d_merged/PvP01-as-ref/PvT01.annotation.gff3',
        'PvC01':  'work/02d_merged/PvP01-as-ref/PvC01.annotation.gff3',
        'MHC087': 'work/02d_merged/PvP01-as-ref/MHC087.annotation.gff3',
    }
    seen = {(r['gene_id'], r['strain']) for r in rows}
    for strain, gff in fallback_gffs.items():
        if not Path(gff).exists():
            continue
        n_added = 0
        for gid, attrs in parse_gff(gff):
            if (gid, strain) in seen:
                continue
            from urllib.parse import unquote
            desc = unquote(attrs.get('description', '')).replace('+', ' ')
            fam = classify(desc)
            rows.append({
                'gene_id': gid,
                'strain': strain,
                'family': fam,
                'evidence': 'annotation-description',
                'description': desc[:120],
            })
            family_counts[strain][fam] += 1
            n_added += 1
        print(f"  {strain:8s}  added {n_added} genes from annotation GFF")

    # Write
    out = OUT_DIR / 'family_table.tsv'
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gene_id', 'strain', 'family', 'evidence', 'description'], delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out}  ({len(rows)} rows)")

    # Per-strain family-count summary
    print("\nPer-strain family counts (top 8):")
    print(f'{"strain":8s} ' + ' '.join(f'{f:>12s}' for f in ['PIR','PHIST','Pv-fam','DBP','MSP','SERA','tRNA','other']))
    for s in ALL_STRAINS:
        c = family_counts[s]
        print(f'{s:8s} ' + ' '.join(f'{c.get(f,0):>12d}' for f in ['PIR','PHIST','Pv-fam','DBP','MSP','SERA','tRNA','other']))

    # Variant-antigen subtable for downstream use
    va_path = OUT_DIR / 'variant_antigens.tsv'
    va_families = {'PIR', 'PHIST', 'Pv-fam', 'STP1', 'RESA', 'TRAg', 'DBP', 'EBA', 'RBP', 'AMA', 'MSP', 'RAP', 'SERA'}
    with open(va_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gene_id', 'strain', 'family', 'evidence', 'description'], delimiter='\t')
        w.writeheader()
        w.writerows(r for r in rows if r['family'] in va_families)
    n_va = sum(1 for r in rows if r['family'] in va_families)
    print(f"\nWrote {va_path}  ({n_va} variant-antigen rows)")


if __name__ == '__main__':
    main()
