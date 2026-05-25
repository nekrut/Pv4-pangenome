# Orthologous gene groups

The consensus ortholog table is the cross-strain backbone of v3 — every per-gene downstream analysis (MSAs, ML trees, HyPhy selection scans, UCSC selection tracks, family-stratified plots) keys off it. We built it by triangulating three independent orthology signals and merging with union-find.

**Related documents.** [PANGENOME.md → How the graph fed downstream analyses](PANGENOME.md#how-the-graph-fed-downstream-analyses) covers evidence stream 3 (graph-path co-membership) and the graph PAV cross-validation that adds `graph_strains` / `graph_mean_pav` columns to the consensus table; [MULTIZ.md → How these alignments fed downstream analyses](MULTIZ.md#how-these-alignments-fed-downstream-analyses) covers the cleaned chains (input to Liftoff/TOGA2 in stream 1) and the rbest chains (input to stream 2); [MALARIAGEN_VCF_PROJECTION.md → How these projected VCFs fed downstream analyses](MALARIAGEN_VCF_PROJECTION.md#how-these-projected-vcfs-fed-downstream-analyses) covers the cross-strain drug-resistance QC that uses gene IDs from this table; [MSA_HYPHY.md](MSA_HYPHY.md) covers the per-orthogroup MSAs, ML trees, and HyPhy BUSTED scans that consume rows from this table; [MICROSYNTENY.md](MICROSYNTENY.md) covers the subtelomeric ribbon plots that visualize orthogroup membership across strains.

## How we built it

Three orthology evidence streams, each producing (strain_a.gene_a, strain_b.gene_b) edges. Union-find over all edges from all three streams → orthogroup IDs.

### Stream 1 — Liftoff + TOGA2/CESAR2 projections

**Inputs** (per anchor strain `A` and non-anchor target `Q`):

| Input                 | Path                                              | What it is                                                                                                                                                                  |
| --------------------- | ------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Anchor annotation     | `inputs/annotations/{A}.fixed.gff3`               | GFF3 from PlasmoDB-68 (PvP01, PvW1, PAM, PvSY56), with chromosome names fixed to match the assembly FASTA                                                                   |
| Anchor assembly       | `genomes/softmasked/{A}.fa`                       | Soft-masked FASTA from the masking step                                                                                                                                     |
| Target assembly       | `genomes/softmasked/{Q}.fa`                       | Soft-masked FASTA of the non-anchor strain                                                                                                                                  |
| Cleaned chain         | `work/01_chains/{A}.{Q}.cleaned.chain.gz`         | Directional UCSC chain from anchor→target; built by the KegAlign chain pipeline ([MULTIZ.md → How we built them](MULTIZ.md#how-we-built-them), "Chain pipeline" subsection) |
| Anchor BED12          | `inputs/annotations/{A}.bed12.gz`                 | Required by TOGA2/CESAR2 to define CDS coordinates of every anchor transcript                                                                                               |
| Anchor isoforms table | `inputs/annotations/{A}.isoforms.tsv.gz`          | Required by TOGA2; 2-col gene_id ↔ transcript_id map                                                                                                                        |
| Feature filter        | `protein_coding_gene`, `ncRNA_gene`, `pseudogene` | Liftoff's `-f` whitelist; without it Liftoff lifts `region` / `chromosome` records which break downstream tools                                                             |

The projection runs in two passes:

1. **Liftoff** does the bulk of the projection. Clean transfers (identity ≥ 0.95, full-length, no frameshift) are kept as-is.
2. **TOGA2/CESAR2** rescues the genes Liftoff marked as low-coverage, partial, frameshifted, or with internal short exons — CESAR's CDS-aware alignment recovers many of these.

The merged outputs per (anchor, target) pair land at `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` along with a per-gene `classification.tsv` tagging each gene as `LIFTOFF_CLEAN`, `CESAR_RESCUE`, `CESAR_PARTIAL`, `MISSING`, `SPLIT`, or `EXTRA_COPY`.

Four anchors (PvP01, PvW1, PAM, PvSY56 — the strains with curated PlasmoDB annotations) × 7 non-anchor targets = 28 projections. Each anchor gene's `Name=` attribute in the projected GFF gives the orthologous gene ID in the target. **Edge emitted**: `(anchor.gene_anchor, target.projected_id)`.

### Stream 2 — Reciprocal-best chain hits

**Inputs** (per unordered pair of strains `{A, B}`):

| Input               | Path                                                                 | What it is                                                                                                                                                                                                    |
| ------------------- | -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| rbest chain         | `work/01_chains/{A}.{B}.rbest.chain.gz`                              | Reciprocal-best chain from the KegAlign chain pipeline; the chainSwap + re-net intersection ([MULTIZ.md → How we built them](MULTIZ.md#how-we-built-them), "Chain pipeline — rbest" subsection)               |
| Per-strain gene BED | `inputs/annotations/{A}.bed12.gz`, `inputs/annotations/{B}.bed12.gz` | BED12 of every gene's CDS coordinates — derived from the strain's annotation GFF via `gff3ToGenePred` + `genePredToBed`. For non-PlasmoDB strains, derived from the Liftoff projection in `work/02d_merged/`. |
| Overlap threshold   | `0.90` (hard-coded)                                                  | Both `g_A` and `g_B` must have ≥ 90 % of their CDS length covered by the same rbest chain segment                                                                                                             |

For each unordered pair of strains, the script walks every rbest chain segment, intersects it with the BED12 coordinates of both strains' genes, and emits a pair `(g_A, g_B)` whenever both reciprocal-overlap thresholds are met. **Edge emitted**: `(strain_a.gene_a, strain_b.gene_b)`.

This catches orthologs that the annotation projection missed — typically genes that Liftoff transferred to the wrong locus or skipped because of repeat-region collisions. It's also the only stream that doesn't depend on a chosen anchor.

### Stream 3 — PGGB graph path co-membership

**Inputs**:

| Input                  | Path                                                                                                         | What it is                                                                                                                                                                                                                  |
| ---------------------- | ------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| PGGB graph (.og)       | `inputs/pggb/pv.og`                                                                                          | 8-way PGGB graph in odgi binary format — symlinked to the v2 build ([PANGENOME.md → How we built the graph](PANGENOME.md#how-we-built-the-graph))                                                                           |
| Per-strain PanSN BED   | per-strain gene BED12 with chrom names prefixed `{SAMPLE}#1#{contig}` to match the graph's PanSN-named paths | Built by rewriting the BED12 chrom column to match `odgi paths` output. The PvP01 version is `work/04_graph_validation/core_1to1.PvP01.pansn.bed`; other strains' BEDs are built on the fly inside `phase_e_graph_edges.py` |
| Node-overlap threshold | `0.90` (hard-coded)                                                                                          | Two genes' CDS regions must traverse ≥ 90 % of the same graph nodes (Jaccard-like overlap of node sets)                                                                                                                     |

For each gene, the script runs `odgi paths` to enumerate which graph nodes the gene's CDS region traverses on that strain's path. For pairs of genes from different strains whose traversal node sets share ≥ 90 %, an edge is added. **Edge emitted**: `(strain_a.gene_a, strain_b.gene_b)`.

This stream is decisive for variant-antigen subtelomeric clusters (PIR, PHIST, vir paralog expansions) — chain rbest is too conservative for tandemly-duplicated loci, but the graph's bubble structure connects paralogs by shared traversal.

### Union-find merge

```
for (strain, gene) in all_annotated_genes:
    add_node((strain, gene))

for edge in stream1 + stream2 + stream3:
    union(edge.a, edge.b)

for connected_component in disjoint_set:
    emit orthogroup_id = OG_NNNNNN
```

For each orthogroup, count strains represented and max-copies per strain. Label as:

| Label                | Definition                                      |     Count |
| -------------------- | ----------------------------------------------- | --------: |
| **CORE-1:1**         | Single copy in all 8 strains                    |     5,778 |
| **CORE-VAR**         | Present in all 8 strains but copy-number varies |       163 |
| **PARTIAL**          | Missing in 1+ strains                           |     1,408 |
| **LINEAGE-SPECIFIC** | Present in ≤ 2 strains                          |       137 |
| **FAMILY**           | Variant-antigen tandem-paralog cluster          |        18 |
| **Total**            |                                                 | **7,504** |

Driver script: `scripts/phase_e_consensus.py` (~500 lines). Wall time: ~30 min on the 8-strain panel.

### Family labels

After the consensus merge we tag every gene with its variant-antigen / functional family using PlasmoDB-description regex on the anchor strains. Families: `PIR`, `PHIST`, `Pv-fam`, `MSP`, `DBP`, `EBA`, `RBP`, `AMA`, `RAP`, `SERA`, `TRAg`, `STP1`, `RESA`. Non-anchor strains inherit the family label through the orthogroup. Skipped HMMER as too slow given PlasmoDB already encodes family-membership in the gene description.

13,135 of the 55,153 gene rows are tagged as variant-antigen family members.

### OrthoFinder3 backup

GENESPACE was the originally-planned synteny-constrained orthology refinement, but `parse_annotations` returned 0 matches on our merged GFFs (filed `jtlovell/GENESPACE#206`). We substituted OrthoFinder3 (MCL clustering on MMseqs2 all-vs-all proteomes) as an independent comparator — the OrthoFinder orthogroups overlap our union-find consensus on 5,329 OGs.

## Outputs

| File                                                                                                                                                                                                                                       |     Rows | Size                      | Where                                                                   |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------: | ------------------------- | ----------------------------------------------------------------------- |
| **`work/03_consensus/ortholog_table.tsv`** — the main 14-column consensus table                                                                                                                                                            |    7,504 | 2.2 MB                    | **Git** (gzipped)                                                       |
| `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` + `.classification.tsv` — per-(anchor,target) merged annotations                                                                                                                     | 56 files | ~5 MB / anchor archive    | **Git** — `work/02d_merged/{anchor}-as-ref_archive.tar.gz` (4 archives) |
| `work/04_graph_validation/core_1to1.PvP01.pansn.bed` — input BED for odgi pav                                                                                                                                                              |    4,186 | <1 MB                     | Git                                                                     |
| `work/04_graph_validation/core_1to1.pav.tsv` — graph PAV per CORE-1:1 region per sample (built by the cross-validation step in [PANGENOME.md → How the graph fed downstream analyses](PANGENOME.md#how-the-graph-fed-downstream-analyses)) |   33,489 | 3.7 MB                    | **Git** (gzipped)                                                       |
| `work/05_families/family_table.tsv` — gene → family tag                                                                                                                                                                                    |   55,153 | 5.4 MB                    | **Git** (gzipped)                                                       |
| `work/05_families/variant_antigens.tsv` — variant-antigen subset                                                                                                                                                                           |   13,135 | 1.2 MB                    | **Git** (gzipped)                                                       |
| `work/05_priorities/gene_priorities.tsv` — drug-resistance + vaccine + research priority gene list                                                                                                                                         |      141 | <100 KB                   | Git                                                                     |
| `work/08_orthofinder/` — OrthoFinder3 substitute output (MCL orthogroups, gene trees, single-copy alignments)                                                                                                                              |        — | 280 MB raw → 88 MB tar.gz | **Dropbox** — `Pv4_v3/work_archives/08_orthofinder_archive.tar.gz`      |
| `work/05_families/manhattan/` — family-stratified Manhattan plots                                                                                                                                                                          |  14 PNGs | <2 MB                     | Git                                                                     |

The 14 columns of `ortholog_table.tsv`: `orthogroup_id, label, n_strains, max_copies`, then one column per strain (`PvP01, Sal-I, PvW1, PAM, PvSY56, PvT01, PvC01, MHC087`), then `graph_strains, graph_mean_pav`. Cell separators: `|` for cross-anchor aliases, `,` for paralog clusters, `-` for absent.

The graph PAV columns are added by the graph cross-validation step described in [PANGENOME.md → How the graph fed downstream analyses](PANGENOME.md#how-the-graph-fed-downstream-analyses) — `graph_strains` is how many strains the graph confirms traverse the region, `graph_mean_pav` is the mean per-sample PAV across strains.

## How orthology fed downstream analyses

Five v3 analyses key off the ortholog table.

**Codon and protein MSAs, ML trees, and HyPhy BUSTED bulk run.** Three sequential per-orthogroup stages — Stage F MSAs filtered by `min_intact` (≥ 7 / ≥ 5 yielding 1,584 / 4,215 orthogroups), Stage G IQ-TREE ML trees per OG, Stage H HyPhy BUSTED tests for positive selection. Full recipe + Pvs230 MNS follow-up worked example: [MSA_HYPHY.md → How we built them](MSA_HYPHY.md#how-we-built-them).

**UCSC hub selection tracks.** The orthogroup_id → PvP01 gene_id map from the ortholog table is the join key for two UCSC tracks: `selection_strict.bb` and `selection_relaxed.bb` (BUSTED q-values as BED12+5), plus `orthogroup_membership.bb` (color by `n_strains`). All three live in the hub at `ucsc_hub/GCA_900093555.2/`. The publishing step (BUSTED-to-BigBed12+5 join via the ortholog table) is documented in [MSA_HYPHY.md → How these outputs fed downstream analyses](MSA_HYPHY.md#how-these-outputs-fed-downstream-analyses).

**Subtelomeric microsynteny.** The chr-end ribbon plots at `writeup/microsynteny/chr{1-14}_{L,R}.png` (28 plots) use the ortholog table to draw connectors between orthologous genes across strains, colored by family label. Drives the variant-antigen narrative (PIR-dense clusters in PvW1, PAM; truncated subtelomeres in Sal-I, PvSY56). Full recipe: [MICROSYNTENY.md](MICROSYNTENY.md).

## Re-running on a different species

The recipe at `v3/pipeline/LOCAL.md` Section 6 walks through this. Prerequisites are documented in their own docs: cleaned + rbest chains and the cross-strain annotation projections come from [MULTIZ.md → How we built them](MULTIZ.md#how-we-built-them) (chains) and Stream 1 of [How we built it](#how-we-built-it) above (annotations); the PGGB graph comes from [PANGENOME.md → How we built the graph](PANGENOME.md#how-we-built-the-graph).

Walk-through:

1. After the soft-masking + KegAlign + chain pipeline ([MULTIZ.md → Re-running on a different species](MULTIZ.md#re-running-on-a-different-species)) and the PGGB graph build ([PANGENOME.md → Re-running on a different species](PANGENOME.md#re-running-on-a-different-species)) have completed, edit `pipeline/species.conf` to set `STRAINS`, `REF_STRAIN`, and `ANCHOR_STRAINS` (the subset with curated annotations — typically reference + 3 others).
2. `bash pipeline/06_consensus.sh` — runs the three-stream consensus and emits `work/03_consensus/ortholog_table.tsv`. Wall: ~30 min.
3. `bash pipeline/07_msa.sh` — MSAs (strict + relaxed). Wall: ~5 hours.
4. `bash pipeline/08_trees.sh` — IQ-TREE per orthogroup. Wall: ~7 hours.
5. `bash pipeline/09_hyphy.sh` — bulk BUSTED. Wall: ~5 hours.

Three scale-dependent parameters:

- For more divergent panels (inter-species), drop the rbest overlap threshold from 0.90 → 0.80 — more orthogroups but slightly more false positives.
- For variant-antigen-heavy genomes, the graph-path stream (stream 3) becomes proportionally more important. If the PGGB graph quality is poor (e.g. inadequate divergence for wfmash to find anchors), the consensus will be biased toward stream 1.
- For panels with very different annotation provenance (e.g. mixing PlasmoDB and NCBI Datasets), Liftoff's `-copies -sc 0.95` should be adjusted down for the NCBI-only strains because NCBI annotations skip many variant-antigen paralogs that PlasmoDB curates.

For the *P. knowlesi* scaffold at `Pk/v1/pipeline/06_consensus.sh`, the `MIN_INTACT_STRICT=6` and `MIN_INTACT_RELAXED=4` thresholds are scaled for 7 strains (vs Pv's 7/5 for 8 strains).

## Galaxy tool-wrap list

The orthology pipeline is mostly composed of small custom Python helpers and one or two existing tools. Custom helpers need wrapping; the existing tools just need workflow assembly.

| Tool                                             | Galaxy state | Note                                                                                                                              |
| ------------------------------------------------ | ------------ | --------------------------------------------------------------------------------------------------------------------------------- |
| **Liftoff** (Stream 1 annotation projection)     | ✅ exists     | IUC; also covered in [MULTIZ.md → Galaxy tool-wrap list](MULTIZ.md#galaxy-tool-wrap-list)                                         |
| **TOGA2 / CESAR2** (Stream 1 rescue pass)        | ⚠ needs wrap | Heaviest container (~3 GB); biocontainer doesn't yet exist; would need to build from upstream                                     |
| **agat_sp_merge_annotations**                    | ✅ exists     | IUC alternative to our custom `phase_c4_merge.py` for the Liftoff + CESAR merge step. Less feature-aware than ours but acceptable |
| **`phase_c2_triage.py`** — 8-rule classification | ⚠ needs wrap | Custom Python — small wrapper                                                                                                     |
| **`phase_c4_merge.py`** — Liftoff + CESAR merge  | ⚠ needs wrap | Custom Python; can substitute agat_sp_merge_annotations if needed                                                                 |
| **`phase_e_rbest_overlap.py`** — stream 2 edges  | ⚠ needs wrap | Custom Python; reads rbest chain + BED                                                                                            |
| **`phase_e_graph_edges.py`** — stream 3 edges    | ⚠ needs wrap | Custom Python; reads odgi paths output                                                                                            |
| **`phase_e_consensus.py`** — union-find merge    | ⚠ needs wrap | Custom Python; emits the 14-column table                                                                                          |
| **odgi paths**                                   | ✅ exists     | IUC; one of the odgi subcommands                                                                                                  |
| **odgi pav** (Phase F)                           | ✅ exists     | IUC                                                                                                                               |
| **gffread** (CDS extraction for Phase F)         | ✅ exists     | IUC                                                                                                                               |
| **MAFFT** (LINSI mode)                           | ✅ exists     | IUC                                                                                                                               |
| **pal2nal**                                      | ✅ exists     | IUC                                                                                                                               |
| **trimAl** (cleaned MSA variant)                 | ✅ exists     | IUC                                                                                                                               |
| **IQ-TREE3**                                     | ⚠ partial    | IQ-TREE2 wrapper exists; needs version bump to v3 (binary is `iqtree3`)                                                           |
| **HyPhy BUSTED**                                 | ✅ exists     | IUC + datamonkey wrappers                                                                                                         |
| **OrthoFinder3**                                 | ✅ exists     | IUC; substitute for the GENESPACE branch we couldn't use                                                                          |

The five custom helpers (`phase_c2_triage.py`, `phase_c4_merge.py`, `phase_e_rbest_overlap.py`, `phase_e_graph_edges.py`, `phase_e_consensus.py`) are all short Python scripts — < 500 lines each — so wrapping is straightforward. Bundle them as a single `pv4-orthology-helpers` tool suite with subcommands and they fit cleanly under one Galaxy package.

Suggested workflow shape: a single `consensus-orthology` workflow that takes:

- Per-strain GFF collection (from Liftoff/TOGA2 workflow output)
- Pairwise rbest chain collection (from chain workflow output)
- One PGGB `.og` graph

…and produces:

- `ortholog_table.tsv` (this doc)
- `family_table.tsv` (the family-label tagging step in "How we built it" above)
- `core_1to1.pav.tsv` (the graph cross-validation step in [PANGENOME.md → How the graph fed downstream analyses](PANGENOME.md#how-the-graph-fed-downstream-analyses))

…all as collection outputs. Category: `COMPARATIVE_GENOMICS` (new).

The downstream MSAs, ML trees, and HyPhy BUSTED runs (see "How orthology fed downstream analyses" above) are separately wrappable as a `selection-screen-from-ortholog-table` workflow — these don't depend on the orthology build, only on the ortholog table + per-strain CDS, so they can be re-run on alternative orthology backends.
