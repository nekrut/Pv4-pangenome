# Panagram 2-way (PvP01 vs PAM) - feasibility notes

## TL;DR

Panagram renders. Server lives at http://127.0.0.1:8050/ after running
`./launch_viewer.sh`. Pangenome tab + Anchor genome tab work; Chromosome tab
needs an explicit chromosome pick from the dropdown.

The PGGB graph at `cactus_2way/pggb_out/*.smooth.fix.gfa` was not used.
Panagram is alignment-free - it builds its own KMC k-mer index from the
input FASTAs.

## Install method that worked

Pip in a Python 3.13 venv. No Docker / Singularity image is published by the
authors; two unrelated community images (`jmcparland/panagram`,
`magwood/panagram`) exist on DockerHub but are not from the Schatz lab and
were not used. The pip install was fast (no compile failures from the KMC
submodule) on this box.

```bash
cd /media/anton/data/sandbox/Pv4/v3/panagram_2way
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip setuptools
git clone --recursive https://github.com/kjenike/panagram.git
cd panagram && pip install .
```

samtools / tabix / bgzip are required at runtime but not pip-installable;
sourced from the existing `bcfmod` conda env (`/home/anton/miniconda3/envs/bcfmod/bin`).
The launcher prepends that to PATH.

## Exact commands run

```bash
# inputs
mkdir -p index/FASTAS index/GFFS
cp inputs/assemblies/PvP01.fa  panagram_2way/index/FASTAS/PvP01.fa
cp inputs/assemblies/PAM.fa    panagram_2way/index/FASTAS/PAM.fa
cp inputs/annotations/plasmodb-68/PvP01.gff3 panagram_2way/index/GFFS/PvP01.gff3
cp inputs/annotations/plasmodb-68/PAM.gff3   panagram_2way/index/GFFS/PAM.gff3
# samples.tsv (tab-separated; columns name, fasta, gff, id, anchor)
cat > panagram_2way/index/samples.tsv <<'EOF'
name	fasta	gff	id	anchor
PvP01	FASTAS/PvP01.fa	GFFS/PvP01.gff3	0	True
PAM	FASTAS/PAM.fa	GFFS/PAM.gff3	1	True
EOF

# anchoring (writes config.yaml + Snakefile)
source venv/bin/activate
export PATH=/home/anton/miniconda3/envs/bcfmod/bin:$PATH
cd panagram_2way/index
panagram index samples.tsv -k 21 -c 8 \
    --kmc.threads 8 --kmc.memory 12 --kmc.use_existing \
    --gff_gene_types protein_coding_gene ncRNA_gene pseudogene \
    --prepare
snakemake --cores 8 -F all

# view
panagram view . --ndebug --port 8050 --host 127.0.0.1
```

The whole index step took ~2 min wall on 8 cores. Peak RAM is bounded by
KMC's `--kmc.memory 12` flag (12 GB).

## Caveats / issues encountered

- **`gff_gene_types` default is wrong for PlasmoDB GFFs.** Default is `gene`;
  PlasmoDB / VEuPathDB uses `protein_coding_gene` (plus `ncRNA_gene`,
  `pseudogene`). With the default, panagram builds the index but the gene
  count column in `anchor/*/chrs.tsv` is zero and the Anchor-genome tab
  shows no gene conservation track. Override with `--gff_gene_types
  protein_coding_gene ncRNA_gene pseudogene` at the `index --prepare` step.

- **PvP01 PlasmoDB GFF matches the FASTA cleanly** (172 unique seqids in the
  GFF body, all present in the FASTA). The companion
  `inputs/annotations/PvP01.genbank.gff3` only covers FLZR* contigs and
  drops the main chromosomes - do not use it.

- **PAM PlasmoDB GFF matches PAM.fa contigs** (`CASCJQ*`) one-to-one.

- **PanSN-named FASTAs from `cactus_2way/pggb_in/` won't work** - KMC
  rejects `#` in genome names. Use the per-genome FASTAs at
  `inputs/assemblies/{PvP01,PAM}.fa` which carry plain
  contig IDs (LT635*, FLZR*, CASCJQ*).

- **Server dies if launched without `setsid`/`nohup` from a backgrounded
  shell.** The `launch_viewer.sh` runs it in the foreground; for detached
  use `nohup ./launch_viewer.sh > view.log 2>&1 & disown` or
  `setsid ./launch_viewer.sh < /dev/null > view.log 2>&1 &`.

- **`panagram view --ndebug` is required** for stable serving; without it
  Flask debug-mode auto-reloader spawns two processes and intermittently
  segfaults on this Python 3.13 / dash 4.1 combination.

- **WebGL warning in headless Chrome** ("WebGL is not supported by your
  browser") is cosmetic - the plots are SVG fallbacks in Plotly and still
  render. In a real browser with WebGL the heatmaps are interactive.

- **Chromosome tab requires explicit dropdown selection.** Deep-link
  arguments (`panagram view . PvP01 LT635612.2 1 1021644`) update the
  coordinate input box but the chromosome dropdown stays empty until the
  user picks one, so the Chromosome tab axes are blank in the initial load.

- **No Docker image published by upstream.** Two community-built images
  exist (`jmcparland/panagram`, `magwood/panagram`) but pip install was
  clean here so containerization gave no advantage.

## Comparison vs odgi viz / BandageNG for the 2-way case

For two haploid genomes panagram and odgi/Bandage answer different
questions and the 2-way overlap is thin. Odgi `viz` (already produced at
`cactus_2way/pggb_out/*.viz_multiqc.png`) shows the path coverage of the
graph linearised along the smoothxg order - it is the right tool for
"where are the bubbles, where do paths diverge", and BandageNG draws the
underlying topology. Panagram does neither: it never aligns sequences. It
counts shared 21-mers in 100-bp windows and renders the resulting
presence/absence matrix as a heatmap per chromosome, plus a pairwise
Jaccard distance, plus annotation-aware conservation tracks. For two
genomes the Jaccard heatmap is a 2x2 cell - not informative - and the
per-chromosome k-mer presence track tells you the same story as the odgi
inv/depth views but at lower resolution and without graph coordinates.
Panagram earns its keep at 5+ genomes where the bitmap shows clade
structure; at N=2 it is a downgrade from odgi for graph questions and
from a dotplot for synteny questions. The Anchor-genome conservation
plot (panagram bins the 21-mer presence across each anchor chromosome
and overlays gene density) is the one panel that has no direct odgi
equivalent and is the main reason to keep panagram around for this pair.

## Files of interest

- `launch_viewer.sh` - start the Dash server.
- `index/` - the panagram index directory (~737 MB; KMC databases dominate at 658 MB).
- `index/samples.tsv`, `index/config.yaml`, `index/Snakefile` - generated by `panagram index --prepare`.
- `index/anchor/{PvP01,PAM}/chrs.tsv` - per-chromosome size + gene_count (sanity-check that the GFF wired up).
- `screenshot_pangenome_tall.png` - Pangenome tab, full layout.
- `cdp_01_pangenome_tab.png`, `cdp_02_anchor_tab.png`, `cdp_03_chromosome_tab.png` - tab-by-tab headless captures.
- `view.log` - last server stdout.

## Repro on a fresh machine

The KMC databases under `index/kmc/` are deterministic given the input
FASTAs and k. To rebuild from scratch:

```bash
rm -rf index/anchor index/kmc index/.snakemake
panagram index samples.tsv -k 21 -c 8 --kmc.threads 8 --kmc.memory 12 \
    --gff_gene_types protein_coding_gene ncRNA_gene pseudogene --prepare
snakemake --cores 8 all
```
