#!/usr/bin/env python3
"""Build hardcoded fallback priority gene list by searching PvP01 PlasmoDB GFF
descriptions. Writes work/05_priorities/fallback_priorities.tsv."""
import re
import csv
from pathlib import Path
from urllib.parse import unquote

GFF = 'inputs/annotations/plasmodb-68/PvP01.gff3'
OUT = Path('work/05_priorities/fallback_priorities.tsv')
OUT.parent.mkdir(parents=True, exist_ok=True)

# (name, search regex against description, category, importance)
WANTED = [
    # Drug resistance
    ('DHFR-TS', r'dihydrofolate reductase', 'drug_resistance', 'critical'),
    ('DHPS', r'dihydropteroate synth', 'drug_resistance', 'critical'),
    ('MDR1', r'multidrug resistance protein 1\b|\bMDR1\b', 'drug_resistance', 'critical'),
    ('MRP1', r'multidrug resistance(-| )associated protein 1\b|\bMRP1\b', 'drug_resistance', 'high'),
    ('MRP2', r'multidrug resistance(-| )associated protein 2\b|\bMRP2\b', 'drug_resistance', 'high'),
    ('CRT-O', r'chloroquine resistance transporter|\bCRT\b', 'drug_resistance', 'critical'),
    ('DHODH', r'dihydroorotate dehydrogenase', 'drug_resistance', 'critical'),
    ('CYTB', r'\bcytochrome b\b', 'mitochondrial_marker', 'critical'),
    ('ATP4', r'sodium-? ?pumping ATPase\b|\bATP4\b|\bP-type ATPase 4\b', 'drug_resistance', 'critical'),
    ('K13', r'kelch protein.*13|kelch.*K13', 'drug_resistance', 'high'),
    ('UBP1', r'ubiquitin (carboxy.terminal )?hydrolase\b.*\b1\b|\bUBP1\b', 'drug_resistance', 'high'),
    ('FNT', r'formate.nitrite transporter', 'drug_resistance', 'high'),
    # Vaccine targets
    ('CSP', r'circumsporozoite protein\b', 'vaccine_target', 'critical'),
    ('MSP1', r'merozoite surface protein 1\b', 'vaccine_target', 'critical'),
    ('MSP3', r'merozoite surface protein 3\b', 'vaccine_target', 'high'),
    ('MSP9', r'merozoite surface protein 9\b', 'vaccine_target', 'moderate'),
    ('AMA1', r'apical membrane antigen 1\b|\bAMA-?1\b', 'vaccine_target', 'critical'),
    ('TRAP', r'thrombospondin.related anonymous protein\b|\bTRAP\b', 'vaccine_target', 'high'),
    ('PVS25', r'\bs25\b.*protein|\bPvs25\b|25 kDa ookinete', 'sexual_stage_transmission', 'critical'),
    ('PVS28', r'\bs28\b.*protein|\bPvs28\b', 'sexual_stage_transmission', 'high'),
    ('PVS230', r'\bs230\b|\bPvs230\b', 'sexual_stage_transmission', 'high'),
    ('PVS48/45', r'\bs48[\/_]?45\b|\bPvs48\b', 'sexual_stage_transmission', 'high'),
    # Invasion / erythrocyte binding
    ('DBP', r'Duffy.binding(-| )protein|\bDBP\b', 'erythrocyte_binding', 'critical'),
    ('DBP2', r'Duffy.binding-?like protein 2|\bDBP-?2\b|\bEBP-?2\b', 'erythrocyte_binding', 'high'),
    ('RBP1', r'reticulocyte.binding protein 1\b|\bRBP-?1\b', 'erythrocyte_binding', 'high'),
    ('RBP2a', r'reticulocyte.binding protein 2a\b|\bRBP-?2a\b', 'erythrocyte_binding', 'critical'),
    ('RBP2b', r'reticulocyte.binding protein 2b\b|\bRBP-?2b\b', 'erythrocyte_binding', 'critical'),
    ('RBP2c', r'reticulocyte.binding protein 2c\b|\bRBP-?2c\b', 'erythrocyte_binding', 'high'),
    ('RBP3', r'reticulocyte.binding protein 3\b|\bRBP-?3\b', 'erythrocyte_binding', 'moderate'),
    ('MAEBL', r'\bMAEBL\b|merozoite adhesive erythrocyte', 'invasion', 'high'),
    # Liver stage
    ('LSA1', r'liver.stage antigen 1\b|\bLSA-?1\b', 'liver_stage', 'high'),
    ('LISP1', r'liver.specific protein 1\b|\bLISP-?1\b', 'liver_stage', 'moderate'),
    ('UIS3', r'\bUIS-?3\b|upregulated in infect.*sporozo', 'liver_stage', 'moderate'),
    # Housekeeping controls
    ('TUB-A', r'alpha-?tubulin', 'translation_housekeeping', 'reference'),
    ('TUB-B', r'beta-?tubulin', 'translation_housekeeping', 'reference'),
    ('ACT', r'actin-?I\b|actin protein\b', 'translation_housekeeping', 'reference'),
    ('GAPDH', r'glyceraldehyde.3.phosphate dehydrogenase', 'metabolism', 'reference'),
    ('EF1A', r'elongation factor 1.alpha|\bEF-?1.?alpha\b', 'translation_housekeeping', 'reference'),
    ('HSP70', r'heat shock protein 70|\bHSP-?70\b', 'translation_housekeeping', 'moderate'),
    # AP2 transcription / chromatin
    ('AP2-G', r'transcription factor with AP2 domain.*g|AP2-?G', 'chromatin_regulator', 'high'),
    ('HP1', r'heterochromatin protein 1\b|\bHP-?1\b', 'chromatin_regulator', 'high'),
]


def parse_gene_records(path):
    """Yield (gene_id, description) for each gene record."""
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('protein_coding_gene',):
                continue
            attrs = {}
            for kv in f[8].split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attrs[k] = v
            gid = attrs.get('ID', '').split(':', 1)[-1]
            desc = unquote(attrs.get('description', '')).replace('+', ' ')
            yield gid, desc


def main():
    records = list(parse_gene_records(GFF))
    print(f'parsed {len(records)} PvP01 protein_coding_gene records')
    rows = []
    for name, pat, cat, imp in WANTED:
        rgx = re.compile(pat, re.I)
        hits = [(gid, desc) for gid, desc in records if rgx.search(desc)]
        if not hits:
            print(f'  no match: {name} -> regex {pat!r}')
            rows.append({'gene_symbol': name, 'PVP01_id': 'TBD',
                         'category': cat, 'importance': imp,
                         'description': '(no match in PvP01 GFF)'})
            continue
        # Take all matching hits
        for gid, desc in hits:
            rows.append({'gene_symbol': name, 'PVP01_id': gid,
                         'category': cat, 'importance': imp,
                         'description': desc[:120]})
    with open(OUT, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['gene_symbol','PVP01_id','category','importance','description'], delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print(f'wrote {OUT} with {len(rows)} rows ({len({r["PVP01_id"] for r in rows if r["PVP01_id"]!="TBD"})} distinct PVP01 IDs)')


if __name__ == '__main__':
    main()
