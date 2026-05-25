# Mockup — Pangenome section on the *P. vivax* organism page

After PRs #1244, #1273, #1274, #1277, #1278 the organism detail page renders three sections under the hero: **Workflows → Assemblies**. The Pangenome bundle adds **one more section below Assemblies**. No new tab, no new top-level route — just a fourth scrolling section that mirrors the existing visual pattern.

```
╔══════════════════════════════════════════════════════════════════════════════╗
║  ┌─────────┐  Plasmodium vivax                                               ║  ← existing
║  │ avatar  │  Taxon 5855 · 8 reference assemblies · Last updated 2026-05    ║    HERO
║  │  ◆◆◆◆   │  [Configure Inputs ▾]                                          ║
║  └─────────┘                                                                 ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Workflows                                                                   ║  ← existing
║  ┌────────────────────────────────────────────────────────────────────────┐ ║    SECTION
║  │ ▸ Run BUSCO (lineage assessment)                  [Configure Inputs]  │ ║
║  │ ▸ Run NCBI Datasets dehydrated bundle build       [Configure Inputs]  │ ║
║  │ ▸ Run organism-scoped Hyphy workflow              [Configure Inputs]  │ ║
║  └────────────────────────────────────────────────────────────────────────┘ ║
╠══════════════════════════════════════════════════════════════════════════════╣
║  Assemblies                                                                  ║  ← existing
║  ⓘ  Perform an analysis in the context of an assembly.                       ║    SECTION
║  ┌──────────────────────────────────────┬───────────────┬─────────┬───────┐ ║
║  │ Strain          │ Accession          │  Length (bp)  │ Contigs │ Year  │ ║
║  ├──────────────────────────────────────┼───────────────┼─────────┼───────┤ ║
║  │ PvP01 (ref)     │ GCA_900093555.2    │  29.0 Mb      │   242   │ 2019  │ ║
║  │ Sal-I           │ GCA_000002415.2    │  27.0 Mb      │  2,747  │ 2009  │ ║
║  │ PvW1            │ GCA_914969965.1    │  29.0 Mb      │    19   │ 2022  │ ║
║  │ PAM             │ GCA_949152365.1    │  29.4 Mb      │    28   │ 2023  │ ║
║  │ PvSY56          │ GCA_003402215.1    │  23.8 Mb      │    14   │  —    │ ║
║  │ PvT01           │ GCA_900093545.1    │  29.0 Mb      │   374   │ 2016  │ ║
║  │ PvC01           │ GCA_900093535.1    │  30.2 Mb      │   425   │ 2016  │ ║
║  │ MHC087          │ GCA_040114635.1    │  29.2 Mb      │   126   │ 2024  │ ║
║  └──────────────────────────────────────┴───────────────┴─────────┴───────┘ ║
║                                                                              ║
║  Genome-distance matrix (sourmash)                                           ║  ← NEW
║  ⓘ  Pairwise MinHash distances among the 8 assemblies (k=31, scaled=1000).  ║   sub-card
║  ┌────────────────────────────────────────────────────────────────────────┐ ║
║  │           PvP01 Sal-I PvW1  PAM  PvSY56 PvT01 PvC01 MHC087            │ ║
║  │  PvP01    ▓▓▓▓▓ ░░░   ░░    ░░    ░░     ░░    ░░    ░░               │ ║
║  │  Sal-I    ░░░   ▓▓▓▓▓ ░░░   ░░░   ░░     ░░░   ░░░   ░░░              │ ║
║  │  PvW1     ░░    ░░░   ▓▓▓▓▓ ░░    ░░░    ░░    ░░    ░░               │ ║
║  │  PAM      ░░    ░░░   ░░    ▓▓▓▓▓ ░░     ░░    ░░    ░░               │ ║
║  │  PvSY56   ░░    ░░    ░░░   ░░    ▓▓▓▓▓  ░░    ░░    ░░               │ ║
║  │  PvT01    ░░    ░░░   ░░    ░░    ░░     ▓▓▓▓▓ ░░    ░░               │ ║
║  │  PvC01    ░░    ░░░   ░░    ░░    ░░     ░░    ▓▓▓▓▓ ░░               │ ║
║  │  MHC087   ░░    ░░░   ░░    ░░    ░░     ░░    ░░    ▓▓▓▓▓            │ ║
║  │                                                                        │ ║
║  │  white = 0.00 (identical)  ░░ = 0.01–0.03  ░░░ = 0.03–0.06            │ ║
║  │  [Download sourmash_dist.tsv]   [Download dendrogram.nwk]              │ ║
║  └────────────────────────────────────────────────────────────────────────┘ ║
╠══════════════════════════════════════════════════════════════════════════════╣
║                                                                              ║
║  Pangenome              [bundle: plasmodium-vivax-v1 · 2026-05 ·  ⓘ More]    ║  ← NEW
║                                                                              ║    SECTION
║  ⓘ  Cross-strain alignments, orthologs, and selection scans across the      ║
║      8 reference assemblies above. Browse on UCSC, query by gene name,       ║
║      or download the bundle.                                                 ║
║                                                                              ║
║  ┌────────────────────────────────────────────────────────────────────────┐ ║
║  │  ▾  Browse on UCSC                                                     │ ║
║  │                                                                        │ ║
║  │     Each assembly's UCSC browser carries 4 pangenome tracks:           │ ║
║  │       • Multi-z alignment (8-way)                                      │ ║
║  │       • Pairwise chains to each other strain                           │ ║
║  │       • Cross-strain gene projections (from 4 PlasmoDB anchors)        │ ║
║  │       • MalariaGEN cohort variants (1,895 samples)                     │ ║
║  │                                                                        │ ║
║  │     Click an assembly to open its UCSC browser:                        │ ║
║  │       [ PvP01 ]  [ Sal-I ]  [ PvW1 ]  [ PAM ]                          │ ║
║  │       [ PvSY56 ] [ PvT01 ] [ PvC01 ] [ MHC087 ]                        │ ║
║  │                                                                        │ ║
║  │     Selection tracks (HyPhy BUSTED) live on the PvP01 browser only.    │ ║
║  └────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
║  ┌────────────────────────────────────────────────────────────────────────┐ ║
║  │  ▾  Find a gene                                                        │ ║
║  │                                                                        │ ║
║  │     Look up an orthogroup by gene name, gene symbol, PlasmoDB          │ ║
║  │     description, or sequence:                                          │ ║
║  │                                                                        │ ║
║  │     ┌──────────────────────────────────────────────────────────┐ [Go] │ ║
║  │     │ PVP01_0415800,  dhfr,  "transmission-blocking", or FASTA │       │ ║
║  │     └──────────────────────────────────────────────────────────┘       │ ║
║  │                                                                        │ ║
║  │     7,504 orthogroups · 5,778 single-copy in all 8 strains             │ ║
║  │     · 1,584 with BUSTED scan (strict) · 4,215 (relaxed)                │ ║
║  │                                                                        │ ║
║  │     Examples:  Pvs230 (PVP01_0415800)  ·  dhfr (PVP01_0526600)         │ ║
║  │                ·  mdr1 (PVP01_1010900)  ·  PIR genes (browse list)     │ ║
║  └────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
║  ┌────────────────────────────────────────────────────────────────────────┐ ║
║  │  ▾  Download the bundle                                                │ ║
║  │                                                                        │ ║
║  │     • Pangenome graph (PGGB, 8-way)              [pv.gfa.gz, 75 MB]   │ ║
║  │                                                  [pv.og, 654 MB]      │ ║
║  │     • Pairwise chains (56 cleaned + 28 rbest)    [chains.tar.gz, 75 MB]│ ║
║  │     • Multi-z MAFs (8 hinges)                    [multiz/, 24 GB]     │ ║
║  │     • Cross-strain annotation projections        [annotations.tar.gz, │ ║
║  │       (4 anchors × 7 targets)                     20 MB]              │ ║
║  │     • Consensus ortholog table                   [ortholog_table.tsv, │ ║
║  │                                                   2.2 MB]             │ ║
║  │     • Per-OG codon + protein MSAs                [core_v3.tar.gz,     │ ║
║  │       (strict: 1,584 · relaxed: 4,215)            29 MB + 100 MB]     │ ║
║  │     • Per-OG ML trees + HyPhy BUSTED JSONs       [trees + hyphy, ~2 GB]│ ║
║  │     • MalariaGEN cohort VCFs (lifted to each)    [VCFs, 148 GB]       │ ║
║  │                                                                        │ ║
║  │     All files mirror the GitHub repo + Dropbox folder. MD5 manifest:  │ ║
║  │     [LARGE_FILES_DROPBOX.tsv]                                          │ ║
║  └────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
║  Provenance: Pv4 v3 (May 2026) · GitHub: nekrut/Pv4-pangenome · Docs:       ║
║  PANGENOME · MULTIZ · ORTHOLOGY · MSA_HYPHY · MALARIAGEN_VCF_PROJECTION      ║
║                                                                              ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

## What the gene-lookup → orthogroup detail page looks like

Reached either by typing in the "Find a gene" box above, or by clicking a gene in any UCSC hub's gene track. Route: `/organisms/5855/pangenome/orthogroup/{gene_id}`.

```
╔══════════════════════════════════════════════════════════════════════════════╗
║  ◀ Back to P. vivax pangenome                                                ║
║                                                                              ║
║  Pvs230 (PVP01_0415800)                          [orthogroup OG_0415800]    ║
║  6-cysteine protein · vaccine candidate · 8 strains · single-copy            ║
║  HyPhy BUSTED q = 0.500 (MNS-corrected · originally 0.056)                   ║
║                                                                              ║
║  ┌────────────────────────────────────────────────────────────────────────┐ ║
║  │  ▾  Codon MSA  (MSAViewer.js)                                          │ ║
║  │     [Protein toggle ▾]   [Download FASTA]   [Color: family · clustal]  │ ║
║  │     ────────────────────────────────────────────────────────────────── │ ║
║  │     PvP01    ATGAAGAAACTCTTTTCTGTTACTGTACTATTGCTCATTTGCATT─────GCT…    │ ║
║  │     Sal-I    ATGAAGAAACTCTTTTCTGTTACTGTACTATTGCTCATTTGCATT─────GCT…    │ ║
║  │     PvW1     ATGAAGAAACTCTTTTCTGTTACTGTACTATTGCTCATTTGCATT─────GCT…    │ ║
║  │     PAM      ATGAAGAAACTCTTTTCTGTTACTGTACTATTGCTCATTTGCATT─────GCT…    │ ║
║  │     …                                                                  │ ║
║  └────────────────────────────────────────────────────────────────────────┘ ║
║                                                                              ║
║  ┌────────────────────────────────┐  ┌──────────────────────────────────┐  ║
║  │  ▾  ML tree (phylotree.js)     │  │  ▾  HyPhy BUSTED summary         │  ║
║  │                                │  │                                  │  ║
║  │      PvP01                     │  │  p-value (single-hit)     0.056  │  ║
║  │       │                        │  │  p-value (MNS-corrected)  0.500  │  ║
║  │       ├── Sal-I                │  │  ω-distribution           …      │  ║
║  │       │                        │  │  positively-selected      none   │  ║
║  │       ├──┬─ PvW1                │  │      sites (q<0.05)              │  ║
║  │       │  └─ PAM                │  │                                  │  ║
║  │       …                        │  │  [Open in HyPhy-Vision ↗]        │  ║
║  └────────────────────────────────┘  └──────────────────────────────────┘  ║
║                                                                              ║
║  ┌────────────────────────────────────────────────────────────────────────┐ ║
║  │  ▾  Cross-references                                                   │ ║
║  │                                                                        │ ║
║  │     PlasmoDB:  PVP01_0415800, PVX_003905, PvP01_04_v1:1,234,567-…      │ ║
║  │     UniProt:   A0A1G4HD37                                              │ ║
║  │     UCSC:      Open PvP01 browser at this gene ↗                       │ ║
║  │     MalariaGEN: 47 variants in the 1,895-sample cohort (see VCF link) │ ║
║  └────────────────────────────────────────────────────────────────────────┘ ║
╚══════════════════════════════════════════════════════════════════════════════╝
```

## Implementation hints (for the PR-3 developer)

- The new section reuses the same `Section` / `KeyValueSection` primitives from `app/views/EntityView/` that the Workflows + Assemblies sections already use.
- "Browse on UCSC" group: 8 buttons rendered from the existing `assembly.ucscBrowserUrl` field (from `catalog/build/ts/build-assemblies.ts:92`). No new code path — just a different button label.
- "Find a gene" search box: front-end takes the query string, calls a new `/api/pangenome/orthogroup-lookup?q={name_or_seq}` endpoint that resolves the orthogroup either via the ortholog table (name lookup) or MMseqs2 (sequence lookup) and redirects to the orthogroup detail page.
- "Download the bundle" list: pulls URLs from the new `Pangenome` catalog entry's manifest (the `manifest.json` BRC will rsync from Dropbox).
- The orthogroup detail page is its own route — `pages/data/organisms/[entityId]/pangenome/orthogroup/[ogId]/index.tsx` — with three panels (MSAViewer.js, phylotree.js, BUSTED summary) reading from the BRC TACC storage URL.
