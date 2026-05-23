#!/usr/bin/env python3
"""Morning status report — what completed, what failed, key numbers, suggested next steps.
Run after the orchestrator finishes (or anytime to get a snapshot)."""
import os
import re
import time
from pathlib import Path
from datetime import datetime

V3 = Path('/media/anton/data/sandbox/Pv4/v3')
os.chdir(V3)


def count_files(pattern):
    from glob import glob
    return len(glob(pattern, recursive=True))


def file_size(path):
    p = Path(path)
    return p.stat().st_size if p.exists() else 0


def head(path, n=15):
    p = Path(path)
    if not p.exists():
        return '(missing)'
    with open(p) as f:
        return ''.join(f.readlines()[:n])


def print_sec(title):
    print(f"\n## {title}\n")


# -----------------------------------------------------------------------------
print(f"# Overnight orchestrator — morning report")
print(f"_Generated {datetime.now().isoformat(timespec='seconds')}_")

# Live job inventory (anything still running?)
print_sec("Background processes still alive")
import subprocess
for label, pat in [
    ('Phase H MSA shards', 'build_8way_msa_v2'),
    ('Phase I trimAl',     'trimal'),
    ('IQ-TREE',            'iqtree'),
    ('HyPhy',              'hyphy'),
    ('CrossMap (A2)',      'CrossMap vcf'),
    ('multiz',             'multiz'),
    ('Docker containers',  None),  # special
]:
    if pat is None:
        n = int(subprocess.run('docker ps -q | wc -l', shell=True, capture_output=True, text=True).stdout.strip() or 0)
    else:
        n = int(subprocess.run(f'pgrep -c -f "{pat}" || echo 0', shell=True, capture_output=True, text=True).stdout.strip())
    print(f"- {label}: **{n}** procs")

# Tasks
print_sec("Per-phase status & key numbers")

# Phase H
nmsa = count_files('work/06_msa/core_v3/*.codon.aln.fa')
print(f"### Phase H (codon MSA)\n- codon MSAs built: **{nmsa}** / 2066 target")
print(f"- summary shards: {count_files('work/06_msa/core_v3/summary_shard*.tsv')}")

# Phase I
nclean = count_files('work/06_msa/core_v3_clean/*.cleaned.fa')
print(f"\n### Phase I (trimAl)\n- cleaned alignments: **{nclean}** / {nmsa}")

# Phase J — trees
ntrees = count_files('work/06_msa/core_v3_trees/*/*.treefile')
print(f"\n### Phase J (IQ-TREE)\n- gene trees: **{ntrees}** / {nclean}")
# quick check: do trees have 8 leaves?
qc_8leaf = 0; qc_lt = 0
import glob as g
for tf in g.glob('work/06_msa/core_v3_trees/*/*.treefile')[:50]:
    txt = Path(tf).read_text()
    n_leaves = len([t for t in re.split(r'[(),]', txt) if t.strip().split(':')[0].strip()])
    if n_leaves == 8: qc_8leaf += 1
    else: qc_lt += 1
print(f"- quality (first 50 trees): {qc_8leaf} with 8 leaves, {qc_lt} other")

# Phase J — HyPhy
nbusted_pri = count_files('work/06_msa/core_v3_hyphy/priority/*/busted.json')
nbusted_bulk = count_files('work/06_msa/core_v3_hyphy/bulk/*/busted.json')
nabsrel = count_files('work/06_msa/core_v3_hyphy/priority/*/absrel.json')
nmeme = count_files('work/06_msa/core_v3_hyphy/priority/*/meme.json')
nfel  = count_files('work/06_msa/core_v3_hyphy/priority/*/fel.json')
print(f"\n### Phase J (HyPhy)\n- priority BUSTED: **{nbusted_pri}**, aBSREL: {nabsrel}, MEME: {nmeme}, FEL: {nfel}")
print(f"- bulk BUSTED: **{nbusted_bulk}**")

# Multiz
print(f"\n### Multiz 8-hinge MAFs")
for hinge in ['PvP01','Sal-I','PvW1','PAM','PvSY56','PvT01','PvC01','MHC087']:
    final = Path(f'work/07_multiz/{hinge}/{hinge}.multiz.maf')
    if final.exists():
        nblock = subprocess.run(f'grep -c "^a " {final}', shell=True, capture_output=True, text=True).stdout.strip()
        print(f"- {hinge}: ✓ {nblock} multi-way blocks ({file_size(final)//1024//1024} MB)")
    else:
        print(f"- {hinge}: ✗ no MAF")

# A2 / Path projection
print(f"\n### Path A2 VCF liftover")
n_a2 = count_files('projection/A2_lastz/Pv4_cohort_on_*.vcf.gz')
print(f"- cohort VCFs: **{n_a2}** / 7")
intersection = Path('work/09_projection_compare/intersection.tsv')
print(f"- intersection table: {'exists' if intersection.exists() else 'missing'}  {intersection if intersection.exists() else ''}")
if intersection.exists():
    print('```\n' + intersection.read_text()[:1200] + '\n```')

# GENESPACE
print(f"\n### GENESPACE")
gres = Path('work/08_genespace/results')
if gres.exists():
    nfiles = sum(1 for _ in gres.rglob('*'))
    print(f"- result dir exists; total files {nfiles}")
    pgm = Path('work/08_genespace/results/pangenome')
    if pgm.exists():
        for f in pgm.rglob('*'):
            if f.is_file():
                print(f"  - {f.relative_to(gres)} ({file_size(f)//1024} KB)")
else:
    print(f"- result dir not created (check logs/overnight/genespace.log)")

# Gene priorities
print(f"\n### Gene priorities")
for label, path in [
    ('researcher draft', 'writeup/gene_priorities.draft.tsv'),
    ('validator output', 'work/05_priorities/gene_priorities.tsv'),
    ('fallback', 'work/05_priorities/fallback_priorities.tsv'),
]:
    n = subprocess.run(f'wc -l < {path} 2>/dev/null || echo 0', shell=True, capture_output=True, text=True).stdout.strip()
    print(f"- {label}: {n} lines  ({path})")

# Quality flags
print_sec("Quality flags")
issues = []
# (a) corrupt JSONs
for d in ['work/06_msa/core_v3_hyphy/priority', 'work/06_msa/core_v3_hyphy/bulk']:
    for j in Path(d).rglob('*.json') if Path(d).exists() else []:
        if file_size(j) < 100:
            issues.append(f"small/empty JSON: {j}")
# (b) empty cleaned alns
for f in g.glob('work/06_msa/core_v3_clean/*.cleaned.fa'):
    if file_size(f) < 100:
        issues.append(f"empty cleaned aln: {f}")
# (c) STATUS.md FAILED entries
status_md = Path('logs/overnight/STATUS.md')
if status_md.exists():
    for ln in status_md.read_text().splitlines():
        if 'FAILED' in ln or 'ERROR' in ln:
            issues.append(ln)
if issues:
    print(f"Found **{len(issues)}** issues (showing first 30):")
    for i in issues[:30]:
        print(f"- `{i}`")
else:
    print("- (none detected)")

# Load history
print_sec("Load / timing graph (samples every 10 min)")
load_log = Path('logs/overnight/load.log')
if load_log.exists():
    print("```")
    for ln in load_log.read_text().splitlines()[-60:]:
        print(ln)
    print("```")
else:
    print("- (no samples; orchestrator just started?)")

# Suggested next steps
print_sec("Suggested next steps for tomorrow")
suggestions = []
if nmsa < 2066:
    suggestions.append(f"Phase H incomplete ({nmsa}/2066). Restart any dead shards via `scripts/run_phase_h_sharded.sh`.")
if nclean < nmsa:
    suggestions.append(f"Phase I trimAl incomplete. Rerun `scripts/overnight/03_trimal.sh`.")
if ntrees < nclean:
    suggestions.append(f"IQ-TREE incomplete ({ntrees}/{nclean}). Rerun `scripts/overnight/04_iqtree.sh`.")
if nbusted_pri < 100:
    suggestions.append(f"HyPhy priority bundle thin ({nbusted_pri} BUSTED). Inspect `logs/overnight/hyphy_priority_*.log` for failures.")
if n_a2 < 7:
    suggestions.append(f"A2 incomplete ({n_a2}/7 cohort VCFs).")
if not intersection.exists():
    suggestions.append("A1/A2/B 3-way comparison not generated.")
if not gres.exists():
    suggestions.append("GENESPACE results dir missing — likely container failure; check logs/overnight/genespace.log")
if not Path('work/05_priorities/gene_priorities.tsv').exists():
    suggestions.append("Validator did not produce final priorities; used fallback.")
# Decisions for user
suggestions.append("Decide: PvP01 still says 'Madagascar' in OUTLINE? (already corrected to Peruvian Amazon — verify against De Meulenaere 2023 BMC Genomics).")
suggestions.append("Decide: strain-naming convention PAM vs PvPAM (currently PAM throughout filesystem).")
if not suggestions:
    print("- All pipeline stages completed cleanly. Next: interpret HyPhy results, write up.")
else:
    for s in suggestions:
        print(f"- {s}")

# Files of interest
print_sec("Key files of interest")
key_files = [
    'work/03_consensus/ortholog_table.tsv',
    'work/05_families/family_table.tsv',
    'work/05_families/variant_antigens.tsv',
    'work/05_priorities/gene_priorities.tsv',
    'work/05_priorities/fallback_priorities.tsv',
    'work/09_projection_compare/intersection.tsv',
    'writeup/BLOG1.md',
    'logs/overnight/STATUS.md',
    'logs/overnight/load.log',
]
for f in key_files:
    p = Path(f)
    if p.exists():
        print(f"- ✓ `{f}` ({file_size(f)//1024} KB)")
    else:
        print(f"- ✗ `{f}` (missing)")
