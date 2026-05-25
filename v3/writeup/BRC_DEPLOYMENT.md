# Deploying Pv4 v3 to BRC-analytics — developer sketch

A 1-page map from analysis outputs to BRC deployment slots. Walks left-to-right: the 8 input genomes → the 5 analysis blocks → which files come out → where they go (Git / Dropbox / UCSC hub) → which BRC-analytics data-model field they fill.

## Where this lands in the BRC UI

The natural host for the bundle is the existing **organism detail page** for *P. vivax* (taxon 5855). PRs #1261, #1263, #1274, #1277, #1278 (all merged May 2026) gave organism detail pages tabbed sections and inline workflow configurators. We add one new tab — **"Pangenome"** — to that page. No new top-level route, no separate `PangenomeView`. The tab surfaces the bundle's downloads + per-OG gene browser + cross-links to the UCSC hubs of each member assembly.

The orthogroup detail page (reached by clicking a gene in any UCSC hub track, see Block 4 below) is a sibling route:
`https://brc-analytics.org/organisms/5855/pangenome/orthogroup/{gene_id}`.

## Starting point: 8 input artifacts

All 8 input genomes are already in the BRC-analytics catalog (assembly accessions verified via `catalog/output/assemblies.json`). Sources for each input type below; everything in this table eventually surfaces on the BRC UCSC track hub at `hgdownload.soe.ucsc.edu/hubs/BRC/` keyed by the GCA accession.

| Input                     | Count | Source                                                                                                                                                                                  | What                                                                                                          |
| ------------------------- | ----: | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| Assembly FASTAs           |     8 | NCBI Datasets (`datasets download genome accession <GCA>`) for all 8; PlasmoDB release 68 for the 5 curated strains as alternate source                                                 | PvP01, Sal-I, PvW1, PAM, PvSY56, PvT01, PvC01, MHC087 (full strain ↔ GCA map in [PANGENOME.md](PANGENOME.md)) |
| Annotation GFFs           |     8 | PlasmoDB-68 (`https://plasmodb.org/common/downloads/release-68/{Strain}/gff/data/PlasmoDB-68_{Strain}.gff`) for PvP01, Sal-I, PvW1, PAM, PvSY56; NCBI Datasets for PvT01, PvC01, MHC087 | 5 PlasmoDB-curated + 3 NCBI-Datasets-only (lower confidence)                                                  |
| Per-strain protein FASTAs |     8 | PlasmoDB-68 `AnnotatedProteins.fasta` for the 5 curated strains; `gffread -y` against the genome + GFF for the 3 NCBI-only strains                                                      | One protein record per gene, used by OrthoFinder3 backup and family-label HMMer pass                          |
| Cohort VCF (per chr)      |    16 | MalariaGEN Pv4 release at `https://www.malariagen.net/data_package/open-dataset-plasmodium-vivax-v4` — 14 nuclear chr + API + MIT, one VCF per chr                                      | 1,895 samples in PvP01 coords; PlasmoDB-style chromosome naming (`PvP01_NN_v1`)                               |

Total starting compute: ~24 hours wall on a 32-core + 1 GPU box.

### Inventory artifact — sourmash genome-distance matrix

A small sub-output produced before the 5 analysis blocks proper, surfaced on the **Assemblies** section of the organism page (not a block of its own). Pairwise MinHash distances among the N member assemblies — for v1 we compute it with sourmash (`sourmash sketch dna -p k=31,scaled=1000` + `sourmash compare`); it replaces the mash matrix from the v3 inventory step.

| File                |   Size | Where                                   | BRC data-model slot                               |
| ------------------- | -----: | --------------------------------------- | ------------------------------------------------- |
| `sourmash_dist.tsv` | <10 KB | **Git** — `work/00_inventory/sourmash/` | `Pangenome.distance_matrix` (new slot)            |
| `sourmash_dist.png` |  20 KB | **Git** — `work/00_inventory/sourmash/` | (rendered server-side or in-browser from the TSV) |
| `dendrogram.nwk`    |  <1 KB | **Git** — `work/00_inventory/sourmash/` | `Pangenome.distance_matrix_dendrogram` (new slot) |

Surface as: a heatmap embedded under the Assemblies table on the organism page, with download links for the TSV and the Newick dendrogram. ~30 sec compute on 8 × 25 Mb assemblies, runs as part of the inventory step in `pipeline/01_inventory.sh`.

## Five analysis blocks → output files → BRC deployment

### Block 1 — Pangenome graph

**v1 plan**: just store the 3 files at BRC and expose them as direct downloads. No visualization yet — graph viz is deferred to a later release.

Doc: [PANGENOME.md](PANGENOME.md)

| File                                |   Size | Where                             | BRC data-model slot                        |
| ----------------------------------- | -----: | --------------------------------- | ------------------------------------------ |
| `pv.gfa.gz` (canonical graph, text) |  75 MB | **Git** — `v2/pggb_out/`          | `Pangenome.graph` (GFA1)                   |
| `pv.og` (odgi binary)               | 654 MB | **Dropbox** — `Pv4_v3/pggb_8way/` | `Pangenome.graph_og` (new slot, optional)  |
| `pv.og.lay` (visualization layout)  | 109 MB | **Dropbox** — `Pv4_v3/pggb_8way/` | `Pangenome.graph_lay` (new slot, optional) |

Surface as: 3 download links in the "Pangenome" tab of the *P. vivax* organism page (or a single "Download graph bundle" button). No UCSC track; no in-browser rendering. Graph visualization (odgi viz thumbnail, interactive bubble plot) is a v2 feature.

### Block 2 — Pairwise chains + multi-way MAFs

**v1 plan**: surface the chains + multi-way alignment on every assembly's UCSC browser page. Each assembly's hub gets one composite (`brc_pangenome_align`) with the assembly-as-hinge bigMaf as the parent + 7 chain sub-tracks (one per other strain).

Doc: [MULTIZ.md](MULTIZ.md)

| File                             | Count | Size   | Where                                  | BRC data-model slot                        |
| -------------------------------- | ----: | ------ | -------------------------------------- | ------------------------------------------ |
| `{src}.{tgt}.cleaned.chain.gz`   |    56 | 60 MB  | **Git** — `work/01_chains/`            | `Pangenome.chain_files[]`                  |
| `{src}.{tgt}.rbest.chain.gz`     |    28 | 15 MB  | **Git** — `work/01_chains/`            | `Pangenome.chain_files[]`                  |
| `{hinge}.multiz.maf.gz`          |     8 | 24 GB  | **Dropbox** — `Pv4_v3/multiz/{hinge}/` | `Pangenome.multiz_alignments[]` (new slot) |
| `{hinge}.multiz.maf.bb` (bigMaf) |     8 | 4.6 GB | **Dropbox** — `Pv4_v3/ucsc_hub/{ACC}/` | `SelectionTrack` of type `bigMaf`          |

Surface as: UCSC composite `brc_pangenome_align` on every assembly's hub (8 hubs total). Default visibility: bigMaf `pack`, chains `hide` (user can expand). Each assembly's hub is reachable from the BRC assembly page via the existing `ucscBrowserUrl` mechanism in `catalog/build/ts/build-assemblies.ts`.

### Block 3 — Cross-strain annotations + ortholog table

Doc: [ORTHOLOGY.md](ORTHOLOGY.md)

| File                                                          | Count | Size         | Where                                                      | BRC data-model slot                         |
| ------------------------------------------------------------- | ----: | ------------ | ---------------------------------------------------------- | ------------------------------------------- |
| `{anchor}-as-ref/{Q}.annotation.gff3` + `.classification.tsv` |    56 | 20 MB tar.gz | **Git** — `work/02d_merged/{anchor}-as-ref_archive.tar.gz` | `Pangenome.merged_annotations[]` (new slot) |
| `ortholog_table.tsv` (7,504 orthogroups × 14 cols)            |     1 | 2.2 MB       | **Git** — `work/03_consensus/`                             | `Pangenome.ortholog_table`                  |
| `family_table.tsv` (55,153 rows)                              |     1 | 5.4 MB       | **Git** — `work/05_families/`                              | `Pangenome.family_table` (new slot)         |
| `annot_from_{anchor}.bb` (BigBed12, per assembly)             |    32 | 6 MB total   | **Dropbox** — `Pv4_v3/ucsc_hub/{ACC}/`                     | `SelectionTrack` of type `bigBed`           |

UCSC hub composite `brc_pangenome_annot` (one per assembly): 4 BigBed sub-tracks (one per anchor).

### Block 4 — Per-OG MSAs + ML trees + HyPhy BUSTED

**v1 plan**: build a gene browser with an AI-agent interface. The user reaches the MSA + tree + BUSTED for a gene in two ways:

1. **From the UCSC browser** — clicking any gene feature on the `annot_from_{anchor}` BigBed track (Block 3) navigates to an orthogroup detail page hosted under the *P. vivax* organism page. Wired via the `url` attribute in the hub's `trackDb.txt`: each gene's details-page URL is `https://brc-analytics.org/organisms/5855/pangenome/orthogroup/$$` (UCSC substitutes the feature name for `$$`). The BRC route looks up the orthogroup from the gene ID via the ortholog table from Block 3 and renders the per-OG panel.
2. **From the AI-agent gene browser** — the user types a gene name (any strain's gene ID, gene symbol, or PlasmoDB description fragment) or pastes a sequence. Name lookup uses the ortholog table + family table from Block 3; sequence lookup uses MMseqs2 against the per-strain protein FASTAs from the starting-point table. Both paths return the same per-OG panel as the UCSC click-through.

Per-OG artifacts (alignments, trees, BUSTED JSONs) will live on **TACC storage allocated to BRC** once handed off — the Dropbox / Git locations below are the source for the BRC-side ingest, not the long-term home. The agent fetches them from the TACC URL on demand; no eager rendering of 5,800 panels.

Doc: [MSA_HYPHY.md](MSA_HYPHY.md)

| File                                             | Count | Size       | Where                                                           | BRC data-model slot               |
| ------------------------------------------------ | ----: | ---------- | --------------------------------------------------------------- | --------------------------------- |
| Codon + protein MSAs (strict, `min_intact ≥ 7`)  | 3,168 | 29 MB raw  | **Git** — `work/06_msa/core_v3_archive.tar.gz`                  | `Pangenome.cds_alignments[]`      |
| Codon + protein MSAs (relaxed, `min_intact ≥ 5`) | 8,430 | 100 MB raw | **Git** — `work/06_msa/core_relaxed_archive.tar.gz`             | `Pangenome.cds_alignments[]`      |
| IQ-TREE per-OG workdirs (strict)                 | 1,584 | 177 MB raw | **Git** — `work/06_msa/core_v3_trees_archive.tar.gz`            | `Pangenome.gene_trees` (new slot) |
| IQ-TREE per-OG workdirs (relaxed)                | 4,215 | 474 MB raw | **Dropbox** — `work_archives/core_relaxed_trees_archive.tar.gz` | `Pangenome.gene_trees`            |
| BUSTED JSON (strict)                             | 1,584 | 380 MB raw | **Git** — `work/06_msa/core_v3_hyphy_archive.tar.gz`            | `SelectionTrack.method = BUSTED`  |
| BUSTED JSON (relaxed)                            | 4,215 | 1.5 GB raw | **Dropbox** — `work_archives/core_relaxed_hyphy_archive.tar.gz` | `SelectionTrack.method = BUSTED`  |
| `selection_strict.bb` (BigBed12+5, PvP01)        |     1 | 120 KB     | **Dropbox** — `Pv4_v3/ucsc_hub/GCA_900093555.2/`                | `SelectionTrack` of type `bigBed` |
| `selection_relaxed.bb` (BigBed12+5, PvP01)       |     1 | 284 KB     | **Dropbox** — `Pv4_v3/ucsc_hub/GCA_900093555.2/`                | `SelectionTrack` of type `bigBed` |
| `orthogroup_membership.bb` (BigBed12, PvP01)     |     1 | 60 KB      | **Dropbox** — `Pv4_v3/ucsc_hub/GCA_900093555.2/`                | `SelectionTrack` of type `bigBed` |

Surface as: (a) AI-agent gene browser — query by gene name / symbol / sequence; agent returns the MSA + tree + BUSTED for the matching orthogroup, with deep-links to HyPhy-Vision for richer rendering. (b) UCSC composite `brc_pangenome_select` on the PvP01 assembly hub (3 BigBed sub-tracks: strict + relaxed BUSTED q-values, orthogroup membership). The UCSC tracks let users see selection on the genome browser; the AI agent lets them ask gene-level questions.

### Block 5 — MalariaGEN cohort VCF projection

**v1 plan**: cohort VCFs become first-class hub tracks on every assembly's UCSC browser. PvP01 hub shows the source cohort (16 per-chr VCFs concatenated, in PvP01 coords); each non-PvP01 assembly hub shows that assembly's lifted cohort VCF. UCSC `vcfTabix` track type — needs the `.vcf.gz` + `.vcf.gz.tbi` pair on the hub.

Doc: [MALARIAGEN_VCF_PROJECTION.md](MALARIAGEN_VCF_PROJECTION.md)

| File                                                 | Count | Size   | Where                              | BRC data-model slot                             |
| ---------------------------------------------------- | ----: | ------ | ---------------------------------- | ----------------------------------------------- |
| `Pv4_{chr}.vcf.gz` (renamed source, PvP01 coords)    |    16 | 24 GB  | **Dropbox** — `Pv4_v3/mg_renamed/` | `Pangenome.cohort_vcf_source` (new slot)        |
| `Pv4_cohort_on_{GCA}.vcf.gz` (lifted to each target) |     7 | 148 GB | **Dropbox** — `Pv4_v3/A2_lastz/`   | `Pangenome.cohort_vcf_projections[]` (new slot) |
| `.vcf.gz.tbi` (tabix index per cohort VCF)           |  7+16 | small  | **Dropbox** — alongside the VCFs   | (sibling index for `vcfTabix` track)            |

Pre-deployment work needed: concat the 16 PvP01-coord per-chr VCFs into one `Pv4_cohort.vcf.gz` for the PvP01 hub (PvP01 currently has only the per-chr split). Tabix-index every cohort VCF (`bcftools index -t`). After that the hub `trackDb.txt` adds a `vcfTabix` track per assembly with sample subsetting controls.

Surface as: a `brc_pangenome_cohort` composite on every assembly hub (8 hubs), one `vcfTabix` sub-track per assembly. Population-genomic signal tracks (iHS, FST, π, Tajima's D, AF-by-country BigWigs/BigBeds) derived from these VCFs are a separate v2 workflow.

## BRC-analytics catalog entry (`catalog/source/pangenomes.yml`)

```yaml
- id: plasmodium-vivax-v1
  species_taxonomy_id: 5855
  version: "2026-05"
  reference_anchor: GCA_900093555.2   # PvP01
  member_assemblies:
    - GCA_900093555.2   # PvP01 (anchor)
    - GCA_000002415.2   # Sal-I  (note: catalogued as GCF_000002415.2 — translate in build step)
    - GCA_914969965.1   # PvW1
    - GCA_949152365.1   # PAM
    - GCA_003402215.1   # PvSY56
    - GCA_900093545.1   # PvT01
    - GCA_900093535.1   # PvC01
    - GCA_040114635.1   # MHC087
  # Track lists below are auto-populated by build_pangenomes from the datacache manifest.
```

Track inventory auto-populated from `https://datacache.galaxyproject.org/brc/data/pangenomes/plasmodium-vivax-v1/manifest.json` (which BRC will produce from the Dropbox staging once we hand it over).

## UCSC track-hub layout

```
hgdownload.soe.ucsc.edu/hubs/BRC/pangenome_plasmodium_vivax_v1/
├── hub.txt
├── genomes.txt
└── {ACC}/                                  # one dir per assembly (8 total)
    ├── trackDb.txt                         # 4 composites: align / annot / select / cohort
    ├── {hinge}.multiz.maf.bb               # Block 2 — bigMaf
    ├── {hinge}.multiz.maf.bb.bai           # (optional; we skipped in v1)
    ├── chains/{ACC}_to_{ACC_T}.chain.gz    # Block 2 — 7 per assembly
    ├── annot_from_{anchor}.bb              # Block 3 — 4 per assembly
    ├── selection_{strict,relaxed}.bb       # Block 4 — PvP01 only
    ├── orthogroup_membership.bb            # Block 4 — PvP01 only
    └── Pv4_cohort_on_{ACC}.vcf.gz + .tbi   # Block 5 — 1 cohort VCF per assembly
```

4 hub composites per assembly:

| Composite name         | Sub-tracks                                             | Visibility default           |
| ---------------------- | ------------------------------------------------------ | ---------------------------- |
| `brc_pangenome_align`  | 1 bigMaf + 7 chains                                    | bigMaf `pack`, chains `hide` |
| `brc_pangenome_annot`  | 4 BigBeds (one per anchor's projection)                | 1 `pack`, 3 `hide`           |
| `brc_pangenome_select` | 2 selection BigBeds + 1 orthogroup BigBed (PvP01 only) | `pack`                       |
| `brc_pangenome_cohort` | 1 `vcfTabix` track (MalariaGEN 1,895 samples)          | `dense`                      |

## Deployment PR sequence (recap from brc-analytics#1279)

| PR | Adds                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | Risk                                |
| --- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------- |
| 1  | LinkML schema for `Pangenome` + new slots (graph_og, multiz_alignments, gene_trees, cohort_vcf_*, family_table, merged_annotations). Empty `pangenomes.yml`. `build_pangenomes` step emitting `pangenomes.json` + `assembly-artifacts.json`. GCA↔GCF translation for Sal-I.                                                                                                                                                                                                                                          | Low — purely additive               |
| 2  | Data ingest of Pv-v1 artifacts (above 5 blocks) from Dropbox to TACC + datacache. `publish_ucsc_hub` step generating per-assembly `trackDb.txt` and `hub.txt` / `genomes.txt`. Rsync to `hgdownload.soe.ucsc.edu/hubs/BRC/`. Concat + tabix-index the cohort VCFs for hub vcfTabix tracks.                                                                                                                                                                                                                           | Low — text manifests + static files |
| 3  | **New "Pangenome" tab on the *P. vivax* organism detail page** (`/organisms/5855`). Tab sections: Bundle downloads (Block 1 graph), Member assemblies (links to each UCSC hub from Block 2), Orthogroup browser (Block 3 + 4 — gene name / sequence lookup → per-OG panel with MSAViewer.js + phylotree.js + HyPhy-Vision deep-link), Cohort VCF (Block 5 download + per-assembly hub links). Plus the orthogroup detail route at `/organisms/5855/pangenome/orthogroup/{id}` that the UCSC click-through points at. | Medium — UI                         |
| 4  | `POPULATION_GENOMICS` + `SELECTION_ANALYSIS` workflow categories. IWC workflow registration for `pangenome-build-and-project` + `hyphy-selection-screen-from-ortholog-table`.                                                                                                                                                                                                                                                                                                                                        | Medium — Galaxy integration         |

Each PR is mergeable independently. Suggested cadence: PR 1 lands first (schema only), then PR 2 in parallel with PR 3, PR 4 last.

## What developers need from me

1. **The Dropbox folder URL** (already shared in [LARGE_FILES_DROPBOX.md](LARGE_FILES_DROPBOX.md)): `https://www.dropbox.com/scl/fo/gx1mta4adubja4bsxxmgm/AOni3YRX8TS1E-saUblX-eo?rlkey=0ksw8a5hhkxhy3sljqyycr7bu&dl=0`
2. **The GitHub repo**: `https://github.com/nekrut/Pv4-pangenome`
3. **The 6 doc files in `v3/writeup/`** that describe each block in detail — the per-block "Outputs" tables in those docs are the source of truth for what BRC catalog entries should list.
4. **MD5 manifest** at `v3/writeup/LARGE_FILES_DROPBOX.tsv` (132 entries with path + size + md5 + essential flag).

What I need from developers:

- A staging area on `datacache.galaxyproject.org/brc/data/pangenomes/plasmodium-vivax-v1/` to rsync the Dropbox contents into.
- A `manifest.json` generator (probably an `aws s3 ls` equivalent on datacache) so the `build_pangenomes` step can auto-populate track lists.
- The 4-PR sequence above implemented.

That's the whole picture.
