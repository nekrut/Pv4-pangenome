# v3 outline — current status (2026-05-22)

*P. vivax* 8-strain pangenome analysis in `/media/anton/data/sandbox/Pv4/v3/`.

**Two tracks, both complete:**
1. **VCF projection** — lift MalariaGEN Pv4 (1,895-sample VCF, PvP01 coords) onto 7 other strains via 3 independent methods
2. **Orthology + per-gene evolution** — gene-level projections, intactness classification, 8-way codon MSAs, HyPhy bundle

## Inputs

### 8 *P. vivax* assemblies

| Strain    | GenBank accession   | Origin / isolate                                 |   Total length |   Contigs | NCBI annot date   | Reference                    |
| :-------- | :------------------ | :----------------------------------------------- | -------------: | --------: | :---------------- | :--------------------------- |
| Sal-I     | GCA_000002415.2     | Salvador I, El Salvador (monkey-adapted)         |        27.0 Mb |     2,747 | 2009-05-06        | Carlton et al. 2008          |
| PvP01     | GCA_900093555.2     | Papua, Indonesia (modern primary reference)      |        29.0 Mb |       242 | 2019-08-01        | Auburn et al. 2016           |
| PvW1      | GCA_914969965.1     | Mauritania, West African (PvW1 isolate)          |        29.0 Mb |        19 | 2022-07-01        | Benavente et al. 2021        |
| PAM       | GCA_949152365.1     | Peruvian Amazon (Pv01-19; PAM = Peruvian AMazon) |        29.4 Mb |        28 | 2023-03-11        | De Meulenaere et al. 2023    |
| PvSY56    | GCA_003402215.1     | Asian field isolate (Auburn PCR-free)            |        23.8 Mb |        14 | none              | Auburn et al. 2018           |
| PvT01     | GCA_900093545.1     | Thai field isolate (Auburn series)               |        29.0 Mb |       374 | 2016-10-07        | Auburn et al. 2016           |
| PvC01     | GCA_900093535.1     | Cambodian field isolate (Auburn series)          |        30.2 Mb |       425 | 2016-12-23        | Auburn et al. 2016           |
| MHC087    | GCA_040114635.1     | Recent field isolate (MHC087)                    |        29.2 Mb |       126 | none              | unpublished (2024)           |

### PlasmoDB ↔ GenBank correspondence

All bp-identical between PlasmoDB and GenBank releases except Sal-I (PlasmoDB carries +6.3 kb of subtelomeric patch over the 2009 GenBank release). For consistency with the PGGB graph, GenBank versions used throughout. Chromosome naming harmonized: PlasmoDB `PvP01_NN_v2` → GenBank `LT635xxx`; `.1` version suffix re-added where PlasmoDB stripped it.

### Annotation tiers
1. **Tier 1 PlasmoDB-68** (4 strains: Sal-I, PvP01, PvW1, PAM) — community-curated, has variant-antigen family labels
2. **Tier 2 NCBI** (PvT01, PvC01) — automated annotation, no family labels
3. **Tier 3 Liftoff from PvP01** (PvSY56, MHC087) — PlasmoDB PvSY56 GFF was 404; MHC087 had no NCBI annotation

### MalariaGEN Pv4 cohort VCFs
16 per-chr multisample VCFs (1,895 samples, ~25 GB) at `/media/anton/scratch/malariagen_pv4/`. CHROMs renamed `PvP01_NN_v1` → `LT635xxx` (GenBank).

### PGGB graph (8-way, from v2)
Hard-masked (longdust+sdust union, 11–19% per strain) before PGGB. Built with `-s 5000 -p 90 -n 8`. ODGI/GFA at `v3/inputs/pggb/pv.{og,gfa}`, 1,318 paths across 8 PanSN-named samples.

A separate **2-way graph (PvP01 + PAM)** also built (`v3/cactus_2way/pggb_out/`) — apples-to-apples comparison with the 8-way; 256 chr1 blocks vs 324 in 8-way (transitive evidence adds blocks).

---

## Completed phases — outputs

### Phase A — inventory

Build the reusable input layer: pull and normalize per-strain FASTAs + annotations, harmonize chromosome naming across PlasmoDB / NCBI / GenBank conventions, soft-mask low-complexity regions (longdust + sdust union), check assembly + annotation quality (mash distance, BUSCO, AGAT GFF validity), and persist the 8-way PGGB graph from v2.

| Item                        | Path                                                      |
| :-------------------------- | :-------------------------------------------------------- |
| Mash 8×8 distance matrix    | `v3/work/00_inventory/mash/`                              | *
| BUSCO `plasmodium_odb10`    | `v3/work/00_inventory/busco/`                             |
| AGAT GFF sanity check       | `v3/work/00_inventory/agat_check/`                        |
| Soft-masked assemblies      | `v3/genomes/softmasked/{GCA}.fa`                          | *
| Mask BEDs (longdust+sdust)  | `v2/genomes/mask_bed/{GCA}.union.bed`                     |
| PlasmoDB-68 annotations     | `v3/inputs/annotations/plasmodb-68/`                      |
| NCBI annotations            | `v3/inputs/annotations/ncbi-datasets/`                    |
| PvP01 variant-antigen list  | `v3/inputs/annotations/plasmodb-68/PvP01.family_list.tsv` |
| PvP01 gene BED (6,508 rows) | `v3/inputs/annotations/PvP01.bed`                         |
| PvP01 GenBank GFF           | `v3/inputs/annotations/PvP01.genbank.gff3`                |
| PlasmoDB → GenBank rename   | `v3/inputs/annotations/PvP01_plasmodb_to_genbank.tsv`     |
| 8 raw assemblies            | `v3/inputs/assemblies/{strain}.fa`                        |
| 8-way PGGB graph (v2 build) | `v3/inputs/pggb/pv.{og,gfa}` → symlink to v2              | *

### Phase B — pairwise alignment chains

Generate all the pairwise alignments needed for VCF projection and multi-way MAF construction. Two aligner families used: wfmash (segment-based, mash-prefilter, fast and PGGB-consistent) and KegAlign (GPU lastZ drop-in, seed-and-extend, denser anchors). KegAlign runs all-vs-all on the 8 assemblies — 28 directed pairs through axt → axtChain → chainSort → chainPreNet → chainNet → netChainSubset → chainStitchId. Two KegAlign parameter sets retained: default (HoxD70, hspthresh 3000, seed 12of19) and TUNED (matrix +100/-100, hspthresh 4500, seed 14of22) to compare subtelomeric noise control.

| Item                             | Count | Path                                                              |
| :------------------------------- | :---: | :---------------------------------------------------------------- |
| A1 wfmash PAFs                   | 7     | `v3/projection/A1_wfmash/PvP01_vs_{GCA}.paf`                      |
| A1 chains                        | 7     | `v3/projection/A1_wfmash/PvP01_to_{GCA}.chain`                    |
| A2 KegAlign AXTs                 | 28    | `v3/projection/A2_kegalign/axt/{GCA}__vs__{GCA}.axt`              | *
| A2 cleaned chains (directed)     | 56    | `v3/work/01_chains/{src}.{tgt}.cleaned.chain`                     | *
| A2 rbest chains (per pair)       | 28    | `v3/work/01_chains/{src}.{tgt}.rbest.chain`                       | *
| Per-genome sizes                 | 8     | `v3/work/01_chains/{GCA}.sizes`                                   |
| A2 staged chains for liftover    | 7     | `v3/projection/A2_lastz/chain/PvP01_to_{GCA}.chain`               |
| 2bit per assembly                | 8     | `v3/projection/A2_kegalign/2bit/{GCA}.2bit`                       |
| KegAlign TUNED vs default panels | 2     | `v3/writeup/synteny_canonical/PvP01_to_PAM__{1,2}_kegalign_*.png` |

### Phases C+D — annotation projection (PvP01 + 3 anchors = 4 × 7 queries done)

Project gene models from each anchor strain onto every other strain, then merge sources. Each anchor → 7 queries pipeline: (C.1) Liftoff fast pass with PlasmoDB feature types and tandem-copy preservation, (C.2) 8-rule triage flagging genes Liftoff handled poorly, (C.3) TOGA2/CESAR2 fallback for CDS-aware re-projection on triaged genes, (C.4) merge Liftoff + CESAR2 outputs into one annotation per query with source + intactness tags. Done from 4 anchors (PvP01, PvW1, PAM, PvSY56) → 21 reciprocal projections that feed Phase E's union-find consensus.

| Stage                    | Path                                                        |
| :----------------------- | :---------------------------------------------------------- |
| C.1 Liftoff raw GFF      | `v3/work/02a_liftoff/{anchor}-as-ref/{Q}/{Q}.lifted.gff3`   |
| C.2 Triage (8-rule)      | `v3/work/02b_triage/{anchor}-as-ref/{Q}/`                   |
| C.3 TOGA2/CESAR2 rescue  | `v3/work/02c_toga/{anchor}-as-ref/{Q}/`                     |
| C.3 TOGA2 prep dirs      | `v3/work/02c_toga/{anchor}_prep_{Q}/`                       | 
| C.4 Merged annotation    | `v3/work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3`    | *
| C.4 Classification table | `v3/work/02d_merged/{anchor}-as-ref/{Q}.classification.tsv` | *

Anchors covered: `PvP01-as-ref/`, `PvW1-as-ref/`, `PAM-as-ref/`, `PvSY56-as-ref/`. 21/21 reciprocal anchor projections done. Intactness column tagged per gene (`I`/`PI`/`UL`/`L`/`PG`/`M`).

### Phase E — consensus orthology

Reconcile the 21 per-anchor classifications into one cross-strain orthogroup table. Union-find on (anchor#gene_id) nodes, with cross-anchor aliases merged via ≥90% reciprocal position overlap on the QUERY strain's coordinates. Native-BED seeding for non-PvP01 strains absorbs the (PVPAM_xxx ↔ PVW1_xxx ↔ PVP01_xxx) ID-convention collisions that otherwise inflate orthogroup counts ~2×. Output labels each orthogroup: CORE-1:1, CORE-VAR, PARTIAL, LINEAGE-SPECIFIC, FAMILY.

| Item                          | Path                                      |
| :---------------------------- | :---------------------------------------- |
| Orthogroup table (7,504 rows) | `v3/work/03_consensus/ortholog_table.tsv` | *
| Builder script                | `v3/scripts/phase_e_consensus.py`         |

Table columns include the 8 strains + `graph_strains` + `graph_mean_pav` added by Phase F.

Composition: 5,778 CORE-1:1 / 163 CORE-VAR / 1,408 PARTIAL / 137 LINEAGE-SPECIFIC / 18 FAMILY. Union-find on (anchor#gene_id) + position aliasing (≥90% reciprocal overlap).

### Phase F — PGGB cross-validation

Independent quality check using the v2 PGGB graph: for each CORE-1:1 orthogroup, extract the PvP01 gene's BED region and call `odgi pav -S` (sample-grouped presence/absence). Confirms that the orthology call is structurally supported by the graph — every strain's path actually traverses the gene's nodes. Promotes orthology evidence from "annotation-derived" to "annotation + graph-confirmed" and feeds graph_strains + graph_mean_pav columns into the consensus table.

| Item                      | Path                                                    |
| :------------------------ | :------------------------------------------------------ |
| CORE-1:1 BED for odgi pav | `v3/work/04_graph_validation/core_1to1.PvP01.pansn.bed` |
| Graph PAV table           | `v3/work/04_graph_validation/core_1to1.pav.tsv`         |

PAV table = 1,460 regions × 8 samples. Phase F appends `graph_strains` + `graph_mean_pav` columns **in-place** to the Phase E `ortholog_table.tsv`.

4,048 / 4,186 CORE-1:1 genes (96.7%) confirmed graph-traversed in all 8 strains at PAV ≥ 0.5.

### Phase G — gene family labels

Tag every gene with its variant-antigen / functional family using PlasmoDB description-field regex (PIR, PHIST, Pv-fam-h, MSP, DBP, EBA, RBP, AMA, RAP, SERA, TRAg, STP1, RESA) — skipped HMMER as too slow given PlasmoDB already encodes family-membership in the description. For non-PlasmoDB strains, cross-strain family labels are propagated through Phase E orthogroups. Output drives the family-stratified analyses (Manhattan plots colored by family, subtelomeric microsynteny, variant-antigen HyPhy bundle).

| Item                     |   Rows | Path                                       |
| :----------------------- | -----: | :----------------------------------------- |
| Full family table        | 55,153 | `v3/work/05_families/family_table.tsv`     |
| Variant-antigen subtable | 13,135 | `v3/work/05_families/variant_antigens.tsv` |
| Generator script         |      — | `v3/scripts/phase_g_families.py`           |

Per-PvP01 counts: 1,205 PIR + 76 PHIST + 23 MSP + 13 SERA + 2 DBP. Cross-strain via orthogroup propagation.

### Phase H — codon-aware MSAs (two stringency cuts)

Build 8-way codon-aware multiple sequence alignments per orthogroup. Pipeline: pyfaidx + GFF parsing extracts CDS per strain → translation to protein → MAFFT-LINSI protein MSA → pal2nal codon back-translation → output FASTA. Two stringency cuts: **strict** requires the gene to be intact (`I` or `PI`) in all 7 queries (min_intact=7), **relaxed** allows up to 2 missing strains (min_intact=5). Strict gives 1,584 MSAs dominated by housekeeping core; relaxed gives 4,215 MSAs that recover most drug-resistance + vaccine target genes plus a small fraction of variant antigens. Sharded build across 8 parallel workers.

| Item                                | Count | Path                                               |
| :---------------------------------- | ----: | :------------------------------------------------- |
| Strict codon MSAs (min_intact=7)    | 1,584 | `v3/work/06_msa/core_v3/{PVP01}.codon.aln.fa`      | *
| Strict protein MSAs                 | 1,584 | `v3/work/06_msa/core_v3/{PVP01}.protein.aln.fa`    | *
| Strict per-shard summaries          |     8 | `v3/work/06_msa/core_v3/summary_shard{0..7}.tsv`   |
| Relaxed codon MSAs (min_intact=5)   | 4,215 | `v3/work/06_msa/core_relaxed/{PVP01}.codon.aln.fa`   | *
| Relaxed protein MSAs                | 4,215 | `v3/work/06_msa/core_relaxed/{PVP01}.protein.aln.fa` | *
| Earlier prototype runs (superseded) | 3,316 | `v3/work/06_msa/8way_{strict,relaxed}_v2/`           |

Builder: `v3/scripts/build_8way_msa_v2.py` (sharded with `--shard N --num-shards M`).

### Phase I — trimAl QC

Column-level alignment trimming via trimAl `-automated1` to remove poorly-aligned segments before downstream tree inference and selection analysis. Applied to both strict and relaxed MSAs; preserves alignment length when possible, drops gappy / inconsistent columns aggressively where needed. Output is the canonical input to IQ-TREE and HyPhy.

| Item                                | Count | Path                                                           |
| :---------------------------------- | ----: | :------------------------------------------------------------- |
| Cleaned strict MSAs (`-automated1`) | 1,584 | `v3/work/06_msa/core_v3_clean/{PVP01_*}.codon.cleaned.fa`      |
| Cleaned relaxed MSAs                | 4,215 | `v3/work/06_msa/core_relaxed_clean/{PVP01_*}.codon.cleaned.fa` |

### Phase J — IQ-TREE + HyPhy bundle

Per-gene maximum-likelihood trees + selection-analysis bundle. IQ-TREE3 with `-m MFP -B 1000` (ModelFinder Plus + ultrafast bootstrap) on every cleaned MSA; alignments with <4 unique sequences (~25% of CORE-1:1 due to ultra-conservation) retried without bootstrap. HyPhy applied at three coverage levels: (1) **bulk BUSTED** on all 1,584 strict CORE-1:1 trees → genome-wide episodic selection screen; (2) **priority bundle** (BUSTED + aBSREL + MEME + FEL) on 41 curated priority genes that have CORE-1:1 alignments; (3) **BUSTED-MH** (`--multiple-hits Double`) on Pvs230 as MNS-robustness check; (4) **variant-antigen bundle** (BUSTED-MH + aBSREL + MEME + FEL) on the 30 variant-antigen genes that survived relaxed orthology filter.

| Item                                          | Path                                                                            |
| :-------------------------------------------- | :------------------------------------------------------------------------------ |
| **Strict ML trees** (`-m MFP -B 1000`), 1,584 | `v3/work/06_msa/core_v3_trees/{PVP01_*}/{PVP01_*}.treefile`                     | *
| Strict IQ-TREE per-gene workdirs              | `v3/work/06_msa/core_v3_trees/{PVP01_*}/` (full IQ-TREE output set)             | 
| **Relaxed ML trees**, 4,215                   | `v3/work/06_msa/core_relaxed_trees/{PVP01_*}/{PVP01_*}.treefile`                | *
| **HyPhy bulk BUSTED** (1,584 jsons)           | `v3/work/06_msa/core_v3_hyphy/bulk/{PVP01_*}/busted.json`                       | *
| **HyPhy priority bundle** (41 × 4 methods)    | `v3/work/06_msa/core_v3_hyphy/priority/{PVP01_*}/{busted,absrel,meme,fel}.json` |
| **HyPhy BUSTED-MH** on Pvs230 (MNS-robust)    | `v3/work/06_msa/core_v3_hyphy/priority_mh/PVP01_0415800/busted_mh.json`         |
| **HyPhy relaxed bulk BUSTED** (4,215 jsons)   | `v3/work/06_msa/core_relaxed_hyphy/{PVP01_*}/busted.json`                       | *
| **HyPhy variant-antigen bundle** (30 × 4)     | `v3/work/06_msa/core_relaxed_hyphy_va/{PVP01_*}/{busted,absrel,meme,fel}.json`  |

Bulk BUSTED hits: 107 at p<0.05, 49 at p<0.01, 20 at p<0.001, 12 at p<10⁻⁴ (10× null enrichment at p<0.001).

### OrthoFinder3 (GENESPACE substitute)

Independent orthology call as triangulation against the consensus Phase E table. GENESPACE was the original choice but failed 4× on parse_annotations FASTA-header / GFF matching despite repeated header rewrites — switched to OrthoFinder3 (3.1.4) which only needs cleaned proteomes. Required pre-processing: dedupe protein IDs (PAM had 50 ENA→mRNA accession collisions), strip non-canonical AA characters (`.`, `*`, `J`, `U`, `O` → `X`). Provides phylogenetic hierarchical orthogroups + a species tree as a side benefit. Loses GENESPACE's riparian synteny plots, but those exist in other forms in `v3/writeup/synteny*/`.

| Item                                             | Path                                                                                             |
| :----------------------------------------------- | :----------------------------------------------------------------------------------------------- |
| Cleaned proteomes (deduplicated, AAs normalized) | `v3/work/08_orthofinder/input/{strain}.fa`                                                       |
| **Orthogroups table** (5,329)                    | `v3/work/08_orthofinder/results/Results_May21/Orthogroups/Orthogroups.tsv`                       |
| Single-copy orthogroups                          | `v3/work/08_orthofinder/results/Results_May21/Orthogroups/Orthogroups_SingleCopyOrthologues.txt` |
| Per-species summary                              | `v3/work/08_orthofinder/results/Results_May21/Comparative_Genomics_Statistics/`                  |
| Phylogenetic hierarchical orthogroups            | `v3/work/08_orthofinder/results/Results_May21/Phylogenetic_Hierarchical_Orthogroups/`            |
| GENESPACE failed attempts (input dirs preserved) | `v3/work/08_genespace/{rawGenomes,clean_proteomes,proteomes,gffs}/`                              |

### Multiz 8-hinge multi-way MAFs

For each strain H ∈ {8 strains}, build a multi-way MAF anchored on H by iteratively running `multiz` over the 7 pairwise MAFs that include H. Each pairwise MAF derived from KegAlign AXT via `axtToMaf` with target-prefix declared (handles +/- strand and AXT direction). Produces 8 different views of the same alignment data, each with a different reference frame — useful for downstream conservation scoring (phyloP/phastCons), comparative annotation, ancestral state reconstruction, or any analysis that benefits from a non-PvP01-centric coordinate system.

| Hinge  | Path                                         | Multi-way blocks |
| :----- | :------------------------------------------- | ---------------: |
| PvP01  | `v3/work/07_multiz/PvP01/PvP01.multiz.maf`   |          278,823 | *
| Sal-I  | `v3/work/07_multiz/Sal-I/Sal-I.multiz.maf`   |          255,439 | *
| PvW1   | `v3/work/07_multiz/PvW1/PvW1.multiz.maf`     |          230,551 | *
| PAM    | `v3/work/07_multiz/PAM/PAM.multiz.maf`       |          261,828 | *
| PvSY56 | `v3/work/07_multiz/PvSY56/PvSY56.multiz.maf` |          128,220 | *
| PvT01  | `v3/work/07_multiz/PvT01/PvT01.multiz.maf`   |          302,162 | *
| PvC01  | `v3/work/07_multiz/PvC01/PvC01.multiz.maf`   |            ~254k | *
| MHC087 | `v3/work/07_multiz/MHC087/MHC087.multiz.maf` |          267,047 | *

Per-pair intermediate MAFs at `v3/work/07_multiz/{hinge}/{hinge}_vs_{other}.maf`.

### VCF projection — all 3 paths complete

Lift the MalariaGEN Pv4 cohort (1,895 samples in PvP01 coords) onto each of the 7 other reference assemblies via three independent methods, so we can intersect and flag method-of-origin per variant. **A1 wfmash** uses the wfmash PAF → chain → CrossMap pipeline (consistent with the alignment evidence that built the graph). **A2 KegAlign** uses KegAlign chains → CrossMap (denser anchors, captures sub-segment homologies wfmash misses). **B graph-native** uses `odgi position` to translate each MalariaGEN position through the v2 PGGB graph node space to each target reference path, then joins with the source MalariaGEN VCF to produce a cohort VCF in target coords (handles complex rearrangements that linear chains lose). All three produce a cohort VCF per target with the original 1,895 sample genotypes and the new coords.

| Path                                                   | Cohort VCF location                                                 | Size range    |
| :----------------------------------------------------- | :------------------------------------------------------------------ | :------------ |
| **A1 wfmash**                                          | `v3/projection/A1_wfmash/Pv4_cohort_on_{GCA_*}.vcf.gz` (+ `.csi`)   | 3.6–17 GB × 7 |
| **A2 KegAlign chains** (canonical for drug-resistance) | `v3/projection/A2_lastz/Pv4_cohort_on_{GCA_*}.vcf.gz` (+ `.csi`)    | 4.5–23 GB × 7 | *
| **B graph-native**                                     | `v3/projection/B_graph/Pv4_cohort_on_{GCA_*}.vcf.gz` (+ `.csi`)     | ~12 GB × 7    |
| MalariaGEN per-chr (LT635xxx-renamed)                  | `v3/projection/A1_wfmash/mg_renamed/Pv4_{01..14,API,MIT}.vcf.gz`    | source        | *
| A2 per-chr intermediates                               | `/media/anton/scratch/A2_lastz_v2/lifted/{GCA_*}/Pv4_*_on_*.vcf.gz` | scratch       |
| B per-chr intermediates                                | `v3/projection/B_graph/cohorts/{GCA_*}/Pv4_*.vcf.gz`                | per-chr       |
| B per-site position TSVs (1.3M sites/target)           | `v3/projection/B_graph/sites/mg_on_{GCA_*}.tsv`                     | lookup        |
| B vg-deconstruct per-target VCFs                       | `v3/projection/B_graph/pangenome_on_{GCA_*}.vcf`                    | source        |
| B path lists (PanSN per target)                        | `v3/projection/B_graph/path_lists/{GCA_*}.paths.txt`                | source        |

### 3-way intersection + Phase 5 QC

Quantitative comparison of the three projection methods. For each target reference: extract (CHROM, POS, REF, ALT) tuples from each method's cohort VCF, then compute pairwise + triple intersections via `comm`. Produces an intersection table per strain. **Phase 5** adds a focused QC: for 19 drug-resistance genes (DHFR-TS, DHPS, MDR1, ATP4, MRP1/2, CRT-o, DHODH, plasmepsins, UBP1, FNT, K12, ...) verify lifted variants land inside the orthologous gene's CDS interval in the target assembly via each method. Quantifies which method actually covers the drug-resistance landscape vs. drops it. Result: A2 wins all 19 (median A2 = 0.72, B = 0.45, A1 = 0.13); A1 fails entirely on 9/19. Triple-method consensus retains 31–73% per strain.

| Item                            | Count | Path                                                        |
| :------------------------------ | ----: | :---------------------------------------------------------- |
| 3-way intersection table        |     7 | `v3/work/09_projection_compare/intersection.tsv`            |
| Per-method per-target site TSVs |    21 | `v3/work/09_projection_compare/{GCA_*}/{A1,A2,B}.sites.tsv` |
| Drug-resistance QC table        |   138 | `v3/work/09_projection_compare/drug_resistance_qc.tsv`      |
| 3-way runner                    |     1 | `v3/scripts/overnight/07b_3way_parallel.sh`                 |
| Phase 5 runner                  |     1 | `v3/scripts/overnight/09_phase5_drug_resistance_qc.py`      |

### Priority gene list (curated + validated)

Two-agent pipeline produced a focused list of biologically meaningful genes for downstream selection analysis (alternative to the bulk 1,584-gene screen). **Researcher agent** surveyed the *P. vivax* literature via web search and drafted 141 candidates across 14 categories (drug resistance, vaccine targets, invasion, erythrocyte binding, variant antigen markers, liver stage, sexual stage / transmission, cytoadherence, metabolism, translation housekeeping, apicoplast, mitochondrial, surface antigens, chromatin regulators). **Validator agent** cross-checked every PVP01 ID against `v3/inputs/annotations/plasmodb-68/PvP01.gff3`, corrected 8 fully-wrong IDs, and resolved 81 TBD entries to real locus tags. Output drives the HyPhy priority bundle.

| Item                                                              | Path                                            |
| :---------------------------------------------------------------- | :---------------------------------------------- |
| **Validated priority list** (141 rows, 134 with mapped PVP01 IDs) | `v3/work/05_priorities/gene_priorities.tsv`     |
| Validator report                                                  | `v3/work/05_priorities/validation_report.md`    |
| Researcher draft                                                  | `v3/writeup/gene_priorities.draft.tsv`          |
| Fallback (built from PlasmoDB GFF regex)                          | `v3/work/05_priorities/fallback_priorities.tsv` |

### 3-way A1 ∩ A2 ∩ B intersection

| strain    |        A1 |        A2 |         B |     A1∩A2∩B | consensus % |
| :-------- | --------: | --------: | --------: | ----------: | ----------: |
| Sal-I     |   968,878 | 2,183,149 | 1,298,850 |     296,755 |         31% |
| **PvT01** | 1,960,586 | 2,709,545 | 1,298,884 | **951,334** |     **73%** |
| PvC01     | 1,426,193 | 2,233,609 | 1,298,754 |     516,554 |         40% |
| PvW1      |   457,394 | 2,715,144 | 1,299,075 |     232,831 |         51% |
| PAM       |   422,442 | 2,713,808 | 1,298,852 |     215,527 |         51% |
| PvSY56    |   669,905 | 2,360,586 | 1,244,151 |     309,820 |         46% |
| MHC087    |   770,694 | 2,673,319 | 1,285,260 |     368,514 |         48% |

Output `v3/work/09_projection_compare/intersection.tsv`. A2 ≫ A1 in raw counts (2–6×). B is stable at ~1.3M sites across targets. Triple-consensus retains 31–73%.

### Phase 5 — drug-resistance gene CDS coverage QC

A2 **wins all 19 drug-resistance genes** evaluated. Median A2 coverage = 0.72, B = 0.45, A1 = 0.13. A1 fails entirely (0%) on 9 of 19 drug-resistance genes. Output `v3/work/09_projection_compare/drug_resistance_qc.tsv` (138 rows).

**Recommendation:** prefer A2 cohort VCFs for downstream drug-resistance / vaccine-target site-level analysis. Use A1∩A2∩B intersection as the conservative high-confidence set. A1 weakest — segment-length floor in wfmash skips short / fragmented gene regions.

---

## Per-gene evolutionary analyses

### Priority gene list (141 rows, validated)
- Curated by research agent (web survey of *P. vivax* drug-resistance, vaccine, invasion, liver-stage, sexual-stage literature)
- Validator agent corrected 8 fully-wrong IDs and resolved 81 TBD → real PVP01 locus tags
- 157 distinct PVP01_xxx IDs across the priority categories
- Output `v3/work/05_priorities/gene_priorities.tsv` (validated) + `v3/work/05_priorities/validation_report.md`

### HyPhy results summary

**Bulk BUSTED (1,584 genes, single-hit):**
- Top hits: PVP01_1255400 (p=0; needs GARD review), PVP01_0815500 (p=1.18e-9), PVP01_0929400 (p=3.04e-7), PVP01_0104700 (p=1.86e-6), PVP01_1027200 (p=3.89e-6)

**Priority bundle (41 genes, full BUSTED+aBSREL+MEME+FEL):**
- Strongest positive selection: CSS (BUSTED p=3.9e-4, invasion), RON2 (p=1.1e-3, invasion), AARP (p=1.7e-2, invasion), Pvs230 (p=0.056, vaccine target), RON4 (19 FEL− sites — heavily purified)
- Vaccine target Pvs230 picked up by 3 methods (BUSTED + aBSREL + MEME)

**Variant antigen HyPhy (30 genes from relaxed set, BUSTED-MH):**
- 0 whole-gene hits at p<0.05 after multi-hit modeling
- 20 RAP + 3 TRAg + 2 PIR + 2 RESA + 1 each MSP, AMA, RBP
- Per-site MEME hits scattered across RAP family; strong purifying signal (high FEL−)

### Pvs230 deep-dive (Appendices A + B of `hyphy_report.md`)

**Appendix A: domain mapping** — 6 MEME hits mapped onto Feng et al. 2022 IVP/ICP/C-term boundaries:
- 4 of 6 MEME hits land in IVP (codons 65, 102, 112, 252) — directly reproduces clinical-cohort findings of Doi 2011 + Feng 2022 with 8 reference strains
- 1 MEME+FEL+ double hit at codon 720 in ICP
- 1 MEME hit at codon 2,080 — outside Doi/Feng surveyed range, novel candidate

**Appendix B: MalariaGEN VCF validation** — 1,895 field-sample SNPs at Pvs230 CDS:
- N/S ratio by domain: IVP 6.18, ICP 1.04, C-term 0.89 — textbook positive-selection-on-immune-domain signature
- 5 of 6 MEME hits have direct AF support at the exact codon: C102Y 25.5%, R112H 6.3% + R112C 5.1%, V252M 6.5%, D720A 13.3% + D720N 13.3% (linked), L2080I 0.3%
- 7 of 9 FEL− codons carry MalariaGEN SNPs — **all synonymous**, frequency 2–73%. Textbook purifying selection.

**BUSTED-MH catch:** The D720A + D720N pair (both at 13.3% AF) turned out to be a **single multi-nucleotide substitution (MNS)** linked across 234 samples encoding D720T. Re-running BUSTED with `--multiple-hits Double` shifts Pvs230 whole-gene p from 0.056 → 0.500 (LRT 4.38 → 0.00). **The whole-gene Pvs230 BUSTED signal is fully driven by one MNS event**, not two parallel diversifying events. IVP per-site hits remain valid.

### Subtelomeric microsynteny (28 plots)

Per-chromosome 5' and 3' end ribbon plots across all 8 strains, colored by gene family. Output `v3/writeup/microsynteny/chr{1-14}_{L,R}.png`. chr1 5' shows PIR-dense clusters across strains; PvW1 and PAM preserve the largest PIR cluster expansions; Sal-I and PvSY56 truncated due to fragmented subtelomeres.

---

## Cross-method comparison plots

- `v3/writeup/synteny_3way/PvP01_to_PAM_5way.png` — 4-panel composite (KegAlign default, KegAlign tuned, wfmash -n 8, PGGB 2-way graph blocks)
- `v3/writeup/synteny_3way/chr1_PvP01_PAM_graph.png` — 7-panel chr1 focus (ribbons + bp coverage + subway map)
- `v3/writeup/synteny_3way/dotplot_PvP01_chr1_chr2_vs_PAM.png` — dot plot
- `v3/writeup/synteny_canonical/PvP01_to_PAM__{1..6}.png` — per-method whole-genome ribbon panels
- `v3/writeup/synteny/`, `v3/writeup/synteny_paf/`, `v3/writeup/synteny_graph/` — per-pair PvP01-to-each-strain plots
- `v3/writeup/microsynteny/` — 28 chromosome-end ribbon plots colored by family

---

## Documents produced this session

- `v3/writeup/BLOG1.md` — 1,927-word BRC-analytics blog post #1 (PvP01+PAM 2-genome pangenome walkthrough); reviewed; PAM geography corrected (Peruvian Amazon, not Madagascar)
- `v3/writeup/hyphy_report.md` + `.pdf` — 18-page HyPhy selection-analysis report with 2 Pvs230 appendices; reviewed
- `v3/writeup/tree_analyses.md` — 33 candidate evolutionary analyses across 8 categories (recommended for next round)
- `v3/writeup/WORKFLOW_PLAN.md` — reproducible workflow plan for porting this to Galaxy + rerunning on other species; 13 priority tools needing wrappers
- `v3/writeup/gene_priorities.draft.tsv` + `v3/work/05_priorities/gene_priorities.tsv` — curated 141-row priority list

---

## Open scientific questions / next steps

1. **Variant-antigen family analysis** — PIR/PHIST/Pv-fam analysis needs per-family clustering, not 8-way single-copy ortholog filter. Original plan called for Pfam HMM + microsynteny anchors; we did neither (used description regex + the strict 8-way pipeline).
2. **GARD recombination pre-screen** — bulk BUSTED top hits not GARD-filtered. PVP01_1255400 (p=0) likely a numerical artifact.
3. **BUSTED-MH on bulk** — the Pvs230 MNS lesson generalizes. Other top hits should be re-tested with `--multiple-hits Double` before publication.
4. **Outgroup** — no rooted selection tests possible without *P. cynomolgi* or *P. knowlesi* ortholog injection.
5. **McDonald-Kreitman** — bridge within-species cohort (MalariaGEN) + between-species 8-strain divergence.
6. **HyPhy RELAX** — Sal-I (monkey-adapted lab clade) vs field-isolate clade; drug-exposed vs drug-naive within MalariaGEN.
7. **Geographic stratification of AF analyses** — Pv4 is ~70% Asia, ~25% Americas, ~5% Africa.

---

## Known infrastructure gaps (carry forward into Galaxy migration)

- `odgi untangle` SDSL int_vector assertion crash on the 8-way graph (filed issue #632)
- GENESPACE finicky on parse_annotations headers (4 failures); substituted with OrthoFinder3
- Datamonkey MCP server OAuth flow incomplete in Claude Code (filed `veg/datamonkey-js-server#383`)
- CrossMap output is unsorted — every chain-based liftover pipeline needs explicit `bcftools sort` step
- IQ-TREE rejects bootstrap on alignments with <4 unique sequences — needs fallback to no-bootstrap mode
- 13 tools need Galaxy wrappers for a full Galaxy port (see `WORKFLOW_PLAN.md` for the priority list)
