# Evolutionary analyses we can do with 1,584 *P. vivax* gene trees

_Author: Anton's research agent / Date: 2026-05-20_

## Overview

We have an 8-strain *Plasmodium vivax* pangenome in a PGGB graph, 5,778 CORE-1:1 orthogroups, and—by tomorrow morning—1,584 codon-aware MSAs with matching IQ-TREE2 ML gene trees (bootstrap 1000). Annotation layers include 13,135 variant-antigen family labels, 141 curated drug-resistance/vaccine/immune-evasion priorities, and MalariaGEN Pv4 cohort VCFs lifted onto all seven non-PvP01 references via three independent projection pipelines. HyPhy per-gene dN/dS (BUSTED + aBSREL + MEME + FEL on priorities, BUSTED sweep on the rest) is already queued, as is a GENESPACE refinement.

Here, we catalogue analyses that exploit the gene trees, codon MSAs, projected VCFs, and pangenome graph we already have—organised by what they tell us biologically and ranked at the end.

## 1. Phylogenetic signal & gene tree discordance

### 1.1 ASTRAL-IV coalescent species tree
- **What it tells us**: a quartet-based summary species tree under the multispecies coalescent, with local posterior support per branch—immediately exposes ILS-driven discordance against any concatenation tree we build. Branch lengths come back in coalescent units (and in substitutions-per-site via the integrated CASTLES-II).
- **Tool**: ASTER v1.23+ (ships ASTRAL-IV, ASTRAL-Pro3, wASTRAL, CASTER), invoked as `astral4`. Conda: `bioconda::aster`.
- **Inputs**: concatenated Newick of the 1,584 gene trees (we have them) plus optional bootstrap collections for wASTRAL. No extra data needed.
- **Difficulty**: low—one command, runs in minutes on 8 taxa.
- **Citation**: Zhang, Nielsen, Mirarab 2025, *MBE* msaf172 (doi:10.1093/molbev/msaf172).

### 1.2 Robinson–Foulds distance distribution
- **What it tells us**: empirical distribution of pairwise RF distances between gene trees and against the ASTRAL/concatenation backbone. Bimodality or a long tail of high-RF genes flags loci with discordant histories—candidates for hidden paralogy, lateral transfer of cassettes between vir/PHIST blocks, or assembly errors.
- **Tool**: `ete3 compare` (ETE v3.1.3) or `dendropy.calculate.treecompare.symmetric_difference`.
- **Inputs**: 1,584 gene trees plus a reference tree (ASTRAL output from §1.1).
- **Difficulty**: low—Python loop, ~minutes.
- **Citation**: Huerta-Cepas, Serra, Bork 2016, *MBE* 33:1635–1638 (doi:10.1093/molbev/msw046); Robinson & Foulds 1981, *Math. Biosci.* 53:131–147.

### 1.3 IQ-TREE2 site- and gene-concordance factors
- **What it tells us**: for every branch of the species tree, the fraction of decisive gene trees (gCF) and informative sites (sCF) that support the bipartition. With 8 taxa the species tree has 5 internal branches; we can map gCF/sCF onto each and identify branches where signal is genuinely weak versus branches where strong but conflicting signal exists—the classic ILS-versus-introgression diagnostic. The 2022 sCFL likelihood-based variant in IQ-TREE 2.2.2+ is preferred over the original site-counting sCF.
- **Tool**: IQ-TREE v2.3+ with `-te species.tree -s concat.fa --scfl 100 --gcf gene_trees.nwk`.
- **Inputs**: ASTRAL species tree (§1.1), concatenated codon supermatrix, individual gene-tree file. We have all three.
- **Difficulty**: low–medium—main cost is building the supermatrix.
- **Citation**: Minh et al. 2020, *MBE* 37:2727–2733 (doi:10.1093/molbev/msaa106); Mo et al. 2022, *MBE* 40:msac215 for sCFL.

### 1.4 Quartet Sampling
- **What it tells us**: per-branch quartet concordance (QC), differential (QD), and informativeness (QI) scores—decomposes branch support into three orthogonal axes: signal strength, asymmetry of discordance (introgression versus ILS), and information content. Particularly powerful with small taxon counts because every internal branch has only one quartet.
- **Tool**: Quartet Sampling v1.3.1 (`quartet_sampling.py`).
- **Inputs**: species tree + concatenated alignment + partition file. We have or can build all three.
- **Difficulty**: medium—needs RAxML/IQ-TREE invocations under the hood; runtime scales with replicate count.
- **Citation**: Pease et al. 2018, *Am. J. Bot.* 105:385–403 (doi:10.1002/ajb2.1016).

### 1.5 TreeShrink outlier branch detection
- **What it tells us**: identifies abnormally long terminal or internal branches across the 1,584 gene trees—markers of misalignment, hidden paralogy, or contamination. We can flag suspect orthogroups before they bias downstream dN/dS or ASTRAL.
- **Tool**: TreeShrink v1.3.9.
- **Inputs**: gene-tree collection + matching MSAs (optional).
- **Difficulty**: low.
- **Citation**: Mai & Mirarab 2018, *BMC Genomics* 19(Suppl 5):272 (doi:10.1186/s12864-018-4620-2).

### 1.6 DiscoVista discordance visualisation
- **What it tells us**: heatmap and ridge-plot summaries of gene-tree support for each species-tree branch—turns §1.3 numbers into a paper-ready figure.
- **Tool**: DiscoVista v1.0 (Docker image).
- **Inputs**: gene trees + species tree + annotation file.
- **Difficulty**: low.
- **Citation**: Sayyari, Whitfield, Mirarab 2018, *Mol. Phylogenet. Evol.* 122:110–115 (doi:10.1016/j.ympev.2018.01.019).

## 2. Recombination detection

Recombination violates the single-tree assumption of every codon-based selection test—running GARD or 3SEQ as a pre-screen before HyPhy is the conservative move. With 1,584 alignments we expect a non-trivial fraction to harbour breakpoints; vir/PHIST families almost certainly do.

### 2.1 HyPhy GARD
- **What it tells us**: maximum-likelihood breakpoint positions and the number of distinct topological partitions per alignment, with KH-test support. Genes with GARD-significant breakpoints can either be excluded or split into recombination-free segments before BUSTED/MEME/FEL.
- **Tool**: HyPhy v2.5.65 (`hyphy gard --alignment <fa>`); Datamonkey container available.
- **Inputs**: 1,584 codon MSAs (we have).
- **Difficulty**: medium—GARD is expensive (~minutes to hours per gene on the long vir loci). Easily parallelised.
- **Citation**: Kosakovsky Pond et al. 2006, *MBE* 23:1891–1901 (doi:10.1093/molbev/msl051); Kosakovsky Pond et al. 2006, *Bioinformatics* 22:3096–3098 (doi:10.1093/bioinformatics/btl474).

### 2.2 3SEQ
- **What it tells us**: triplet-based scan for mosaic structure; orders of magnitude faster than GARD, and the 2018 algorithmic update makes 1,584 alignments trivial. Use as a fast triage filter before launching GARD on suspect cases.
- **Tool**: 3SEQ v1.8.
- **Inputs**: codon MSAs.
- **Difficulty**: low—single binary, seconds per gene.
- **Citation**: Lam, Ratmann, Boni 2018, *MBE* 35:247–251 (doi:10.1093/molbev/msx263).

### 2.3 PhiPack PHI / NSS / Max-χ²
- **What it tells us**: three complementary recombination tests in one package—PHI is robust to rate heterogeneity, NSS to phylogenetic structure, Max-χ² to small samples. Triangulating all three reduces false positives.
- **Tool**: PhiPack (Bruen 2006).
- **Inputs**: codon MSAs.
- **Difficulty**: low.
- **Citation**: Bruen, Philippe, Bryant 2006, *Genetics* 172:2665–2681 (doi:10.1534/genetics.105.048975).

### 2.4 RDP5 ensemble scan
- **What it tells us**: GUI-or-CLI driven ensemble of seven recombination detectors (RDP, GENECONV, BOOTSCAN, MAXCHI, CHIMAERA, SISCAN, 3SEQ) with consensus calls and refined breakpoints. Best for the variant-antigen loci where we expect mosaic structure.
- **Tool**: RDP5 v5.51.
- **Inputs**: codon MSAs (Phylip/FASTA).
- **Difficulty**: medium—RDP5 is interactive by default; batch mode requires a Windows VM or Wine on Linux.
- **Citation**: Martin et al. 2021, *Virus Evol.* 7:veaa087 (doi:10.1093/ve/veaa087).

## 3. Topology / tree-distance methods for population history

### 3.1 Per-gene site-concordance landscape
- **What it tells us**: rather than averaging sCF across the genome (§1.3), bin sCF values along chromosomes and visualise—reveals concordance valleys around centromeres, sub-telomeric vir/PHIST blocks, and recombination hotspots. The Pv genome is small (~30 Mb) so the landscape is publication-grade in one figure.
- **Tool**: IQ-TREE2 `--scfl` per gene + custom plot (pandas + matplotlib).
- **Inputs**: per-gene MSA + gene tree + species tree.
- **Difficulty**: medium—loop of 1,584 IQ-TREE runs plus tidying.
- **Citation**: same as §1.3.

### 3.2 Twisst topology weighting
- **What it tells us**: with 8 taxa partitioned into 3–4 groups (PlasmoDB-curated versus NCBI-only; geographic origin; lab-adapted Sal-I versus field), Twisst quantifies the genome-wide weight of each possible topology in sliding windows or per-gene-tree—exposes introgression signal that single-locus sCF cannot resolve.
- **Tool**: Twisst (`twisst.py`).
- **Inputs**: gene trees, taxon partition file.
- **Difficulty**: low–medium—straightforward once partitions are defined.
- **Citation**: Martin & Van Belleghem 2017, *Genetics* 206:429–438 (doi:10.1534/genetics.116.194720).

### 3.3 Chromosomal RF sketch
- **What it tells us**: RF distance from every gene tree to the ASTRAL backbone, plotted along chromosomes—cheap sanity check on §3.1.
- **Tool**: `dendropy` or `ete3` + rolling-median script.
- **Inputs**: gene trees + backbone tree.
- **Difficulty**: low.
- **Citation**: Robinson & Foulds 1981, *Math. Biosci.* 53:131–147.

### 3.4 TreeMix / SpaceMix population graph (skip for now)
- **What it tells us**: admixture-graph inference. With 8 single-strain references and no within-population replicates, TreeMix is underpowered—belongs in the Pv4 cohort layer, not the strain layer.
- **Citation**: Pickrell & Pritchard 2012, *PLoS Genet.* 8:e1002967 (doi:10.1371/journal.pgen.1002967).

## 4. Selection beyond per-gene dN/dS

### 4.1 HyPhy RELAX (test versus reference branches)
- **What it tells us**: a single parameter *k* contrasts the strength of selection on a labelled "test" branch set versus the rest of the tree—*k* > 1 means intensified selection, *k* < 1 relaxed. Two obvious comparisons: (a) Sal-I lab-adapted lineage versus field isolates, asking whether decades of culture relaxed selection on invasion genes; (b) Monkey-derived strains (MHC087) versus human isolates.
- **Tool**: HyPhy v2.5.65 `hyphy relax`.
- **Inputs**: codon MSAs + gene trees with branch labels.
- **Difficulty**: medium—needs branch labelling pipeline; runtime modest.
- **Citation**: Wertheim et al. 2015, *MBE* 32:820–832 (doi:10.1093/molbev/msu400).

### 4.2 HyPhy CONTRAST-FEL
- **What it tells us**: site-by-site detection of dN/dS differences between two or more pre-specified branch sets—complementary to RELAX, which gives a single *k*. Where RELAX answers "is selection different on average?", CONTRAST-FEL answers "which residues drive the difference?".
- **Tool**: HyPhy v2.5.65 `hyphy contrast-fel`.
- **Inputs**: codon MSAs + labelled trees.
- **Difficulty**: medium.
- **Citation**: Kosakovsky Pond et al. 2021, *MBE* 38:1184–1198 (doi:10.1093/molbev/msaa263).

### 4.3 HyPhy SLAC
- **What it tells us**: a counting-based per-site dN/dS estimator that is fast, deterministic, and complementary to the parametric FEL/MEME. Including SLAC in the standard panel doubles confidence in site-level calls when methods agree.
- **Tool**: HyPhy v2.5.65 `hyphy slac`.
- **Inputs**: codon MSAs + gene trees.
- **Difficulty**: low—cheap addition to the existing HyPhy queue.
- **Citation**: Kosakovsky Pond & Frost 2005, *MBE* 22:1208–1222 (doi:10.1093/molbev/msi105).

### 4.4 GARD-aware BUSTED post-hoc
- **What it tells us**: re-runs BUSTED on the recombination-free partitions identified in §2.1—removes the most common cause of inflated false-positive dN/dS calls. Required before publishing any BUSTED hit on a vir or PHIST gene.
- **Tool**: HyPhy `hyphy busted` re-invoked on `--alignment.gard.partitions`.
- **Inputs**: GARD output (§2.1) + codon MSAs.
- **Difficulty**: low—wrapper around the existing queue.
- **Citation**: Murrell et al. 2015, *MBE* 32:1365–1371 (doi:10.1093/molbev/msv035).

### 4.5 Tajima's *D* and Fay–Wu *H* on MalariaGEN VCFs
- **What it tells us**: within-species summary statistics on the projected Pv4 cohort VCFs—negative Tajima's *D* signals population expansion or purifying selection; positive *D* and negative Fay–Wu *H* point to balancing selection (the classic vaccine-antigen signature). Run per gene across the 141 priorities and cross-reference with the BUSTED/MEME hits.
- **Tool**: scikit-allel v1.3.5 (`allel.tajima_d`, `allel.windowed_tajima_d`); VCFtools v0.1.16 for sanity-check.
- **Inputs**: projected Pv4 VCFs on PvP01 (we have all three projections—A1 wfmash, A2 KegAlign, B graph), gene coordinates.
- **Difficulty**: low–medium—Python notebook, hour or two.
- **Citation**: Tajima 1989, *Genetics* 123:585–595; Fay & Wu 2000, *Genetics* 155:1405–1413 (doi:10.1093/genetics/155.3.1405).

### 4.6 McDonald–Kreitman test (within Pv4 polymorphism versus 8-strain divergence)
- **What it tells us**: 2×2 contingency table of nonsynonymous and synonymous changes—polymorphic in Pv4 versus fixed between strains—classifies each gene as under positive, neutral, or purifying selection. The MalariaGEN VCFs supply polymorphism, the 8-strain MSA supplies divergence; the asymptotic-MK extension corrects for weakly deleterious mutations.
- **Tool**: `asymptoticMK` web service or `iMKT` R package; custom counting script around VCFs + MSAs.
- **Inputs**: codon MSAs + gene-coordinate-filtered Pv4 VCFs. We have both.
- **Difficulty**: medium—needs careful synonymous-site annotation and ancestral allele calls (a single outgroup would tighten this).
- **Citation**: McDonald & Kreitman 1991, *Nature* 351:652–654; Murga-Moreno et al. 2019, *Nucleic Acids Res.* 47:W283–W288 (doi:10.1093/nar/gkz372).

## 5. Functional / pathway-level analyses

### 5.1 GO enrichment of selection hits
- **What it tells us**: which GO biological processes, molecular functions, and cellular components are over-represented among BUSTED-significant genes versus the 1,584-gene background—turns per-gene calls into biological narrative.
- **Tool**: topGO v2.62 (R/Bioconductor) with `weight01` and `elim` algorithms.
- **Inputs**: PlasmoDB-68 GO annotations for PvP01 (we have), BUSTED p-values from the queued sweep.
- **Difficulty**: low.
- **Citation**: Alexa, Rahnenführer 2025, *Bioconductor* (doi:10.18129/B9.bioc.topGO); Alexa et al. 2006, *Bioinformatics* 22:1600–1607.

### 5.2 KEGG / MPMP pathway mapping
- **What it tells us**: drug-resistance and metabolic pathway enrichment among positively selected genes; particularly useful for catching co-evolution of antifolate, artemisinin, and 8-aminoquinoline resistance loci.
- **Tool**: `clusterProfiler` v4.16 (R/Bioconductor) with the KEGG `pvx` organism database; cross-check against the Malaria Parasite Metabolic Pathways (MPMP) resource.
- **Inputs**: BUSTED hits with KEGG IDs.
- **Difficulty**: low.
- **Citation**: Wu et al. 2021, *The Innovation* 2:100141 (doi:10.1016/j.xinn.2021.100141).

### 5.3 GSEA-style ranking against Pv-specific functional categories
- **What it tells us**: instead of arbitrary p-value cutoffs, rank all 1,584 genes by BUSTED LRT statistic and ask whether the variant-antigen families (vir, PHIST, Pv-fam-a/b/c), invasion genes (DBP, MSP, RBP), housekeeping ribosomes, and antifolate pathway are enriched at the top. Uses our family-table annotation directly.
- **Tool**: `fgsea` v1.32 (R/Bioconductor).
- **Inputs**: family-table annotations + BUSTED LRT scores.
- **Difficulty**: low.
- **Citation**: Korotkevich et al. 2021, *bioRxiv* 060012 (doi:10.1101/060012).

### 5.4 Co-evolution / coupling analysis on individual gene families
- **What it tells us**: positions whose substitutions co-vary across strains—candidate compensatory pairs, contact residues, or epistatic loci. With 8 sequences direct-coupling analysis is statistically underpowered; we can rescue power by augmenting each MSA with Pv4 cohort variants (treating each polymorphism as an additional row) before running DCA.
- **Tool**: `plmDCA` / `EVcouplings` v0.2; CCMpred for fast pseudolikelihood DCA.
- **Inputs**: codon MSAs augmented with Pv4 polymorphism (we have the components).
- **Difficulty**: medium-high—augmentation pipeline is non-trivial; DCA needs N_eff > 100 per family.
- **Citation**: Morcos et al. 2011, *PNAS* 108:E1293–E1301 (doi:10.1073/pnas.1111471108); Hopf et al. 2019, *Bioinformatics* 35:1582–1584.

## 6. Pangenome-specific analyses

### 6.1 Pan/Core/Cloud frequency versus dN/dS
- **What it tells us**: whether accessory genes (cloud, present in 2–6 of 8 strains) experience different selection regimes than core 1:1 orthologs—the classic pangenome-evolution question. Plot BUSTED ω against `graph_strains` count from our consensus orthology table.
- **Tool**: custom R/Python; `ggplot2` for the figure.
- **Inputs**: `work/03_consensus/ortholog_table.tsv` + BUSTED output. We have the table; BUSTED is queued.
- **Difficulty**: low.
- **Citation**: McInerney, McNally, O'Connell 2017, *Nat. Microbiol.* 2:17040 (doi:10.1038/nmicrobiol.2017.40).

### 6.2 Variant-antigen network structure
- **What it tells us**: build a graph where nodes are vir/PHIST/Pv-fam orthogroups and edges weight pairwise sequence similarity or shared sub-telomeric position—community detection reveals whether vir subfamilies (A through E) form coherent evolutionary modules or freely recombine.
- **Tool**: `igraph` v2.0 + `leidenalg` for Leiden community detection; BLAST or MMseqs2 v15 for the similarity graph.
- **Inputs**: variant-antigen MSAs + family-table labels. We have both.
- **Difficulty**: medium.
- **Citation**: Traag, Waltman, van Eck 2019, *Sci. Rep.* 9:5233 (doi:10.1038/s41598-019-41695-z).

### 6.3 Microsynteny breakpoints versus selection peaks
- **What it tells us**: do genes flanking PGGB-graph breakpoints—where the 8 paths diverge in topology—show elevated dN/dS? A positive correlation would link structural rearrangement to adaptive evolution; absence would argue that breakpoints sit in selectively unconstrained intervals.
- **Tool**: `odgi extract` on `inputs/pggb/pv.og` for breakpoint coordinates; `bedtools` v2.31 for the gene overlap; permutation test in Python.
- **Inputs**: PGGB graph (we have), gene coords, BUSTED output.
- **Difficulty**: medium.
- **Citation**: Garrison et al. 2024 (PGGB), *Nat. Methods* 21:1430–1438 (doi:10.1038/s41592-024-02430-3).

### 6.4 PAV phylogeny versus sequence phylogeny
- **What it tells us**: construct a binary presence-absence matrix across all 13,135 orthogroups, build a Dollo-parsimony or maximum-likelihood PAV tree, and compare against the ASTRAL species tree—mismatches flag massive lineage-specific gene gain/loss, candidate horizontal transfer hotspots, or assembly artefacts.
- **Tool**: `IQ-TREE2 -st BIN` with `+ASC` ascertainment correction; `Count` v10.04 for Dollo parsimony.
- **Inputs**: orthology table (we have).
- **Difficulty**: low–medium.
- **Citation**: Csurös 2010, *Bioinformatics* 26:1910–1912 (doi:10.1093/bioinformatics/btq315).

## 7. Drug-resistance evolution

### 7.1 Phylogenetic mapping of resistance mutations
- **What it tells us**: project known *P. vivax* resistance variants—*dhfr* F57L/S58R/T61M/S117N/I173L, *dhps* A383G/A553G, *mdr1* Y976F, *crt-o* K10 insertions—onto the IQ-TREE gene tree for each locus and infer along which species-tree branch each mutation arose. With 8 reference strains plus the projected Pv4 cohort, we can date relative emergence.
- **Tool**: custom Python on the gene trees + `mapeq` / `treesnatcher` for visualisation; `ggtree` v3.16 for the figure.
- **Inputs**: per-gene trees and MSAs for the resistance loci + Pv4 VCFs.
- **Difficulty**: low–medium.
- **Citation**: Auburn et al. 2019, *eLife* 8:e44073 (doi:10.7554/eLife.44073).

### 7.2 Ancestral state reconstruction on resistance loci
- **What it tells us**: maximum-likelihood (or empirical Bayes) ancestral codons at each internal node of the gene tree—reveals whether a resistance allele is ancestral to the species, multiply-emerged, or lost in specific lineages. Adding a *P. cynomolgi* outgroup (see §8) is critical for rooting.
- **Tool**: PAML v4.10 (`codeml` with `RateAncestor=1`) or HyPhy `AncestralSequences.bf`.
- **Inputs**: rooted gene tree + codon MSA. Rooting needs §8.1.
- **Difficulty**: medium—main lift is the outgroup.
- **Citation**: Yang 2007, *MBE* 24:1586–1591 (doi:10.1093/molbev/msm088).

### 7.3 Selection scan around resistance loci in Pv4
- **What it tells us**: iHS, nSL, and XP-EHH within the projected Pv4 cohort flanking each resistance locus—detects recent hard sweeps. Complements the inter-strain dN/dS by capturing intra-population recent selection that is invisible to 8-strain analyses.
- **Tool**: `selscan` v2.0; `hapbin` for very large cohorts.
- **Inputs**: phased Pv4 VCFs (we have unphased—phasing with `shapeit5` is an extra step).
- **Difficulty**: medium—phasing is the bottleneck.
- **Citation**: Szpiech & Hernandez 2014, *MBE* 31:2824–2827 (doi:10.1093/molbev/msu211).

## 8. Outgroup-dependent extensions

The 8-strain reference set has no outgroup. Most analyses above run on unrooted trees; ancestral state reconstruction (§7.2), McDonald–Kreitman polarisation (§4.6), and rooted dN/dS branch tests all benefit from one. We recommend adding *P. cynomolgi* (sister to *P. vivax*, full PacBio assembly available—GenBank GCA_000956335) as the primary outgroup, with *P. knowlesi* and *P. inui* as backups.

### 8.1 Add *P. cynomolgi* as outgroup
- **What it tells us**: roots all 1,584 gene trees, polarises every substitution, and allows the full Yang/Nielsen branch-site dN/dS suite plus rigorous ancestral reconstruction.
- **Tool**: re-run the orthology + MSA + IQ-TREE pipeline with the *cynomolgi* proteome added; or use ProteinOrtho/OrthoFinder lightweight extension.
- **Inputs**: *P. cynomolgi* genome + GFF (publicly available); existing pipeline.
- **Difficulty**: medium—~half a day of pipeline re-run; downstream analyses become rooted.
- **Citation**: Tachibana et al. 2012, *Nat. Genet.* 44:1051–1055 (doi:10.1038/ng.2375).

### 8.2 Branch-site test for episodic positive selection on the *P. vivax* stem
- **What it tells us**: once rooted with *cynomolgi*, the *P. vivax* MRCA branch becomes a foreground for codeml's branch-site model A test—identifies genes that experienced positive selection on the lineage leading to all extant *P. vivax*. Classic test for host-switch adaptation.
- **Tool**: PAML v4.10 `codeml` model A versus null; or HyPhy aBSREL on the labelled stem.
- **Inputs**: rooted gene trees from §8.1 + codon MSAs.
- **Difficulty**: medium.
- **Citation**: Zhang, Nielsen, Yang 2005, *MBE* 22:2472–2479 (doi:10.1093/molbev/msi237).

## Priorities — top 5 analyses we'd run first

1. **ASTRAL-IV species tree (§1.1)** — single command, gives us the backbone every other discordance analysis needs, and is the natural figure-1 panel of the paper.
2. **HyPhy GARD pre-screen + 3SEQ triage (§2.1, §2.2)** — recombination is the most common source of false positives in the queued BUSTED sweep; running 3SEQ as a fast first pass and GARD on flagged loci will save us from publishing inflated dN/dS calls on vir/PHIST genes.
3. **IQ-TREE2 gCF/sCFL + Twisst landscape (§1.3, §3.1, §3.2)** — converts 1,584 trees into a per-chromosome concordance landscape, immediately visible alongside the BUSTED results.
4. **HyPhy RELAX on Sal-I versus field isolates (§4.1)** — directly tests the long-standing question of whether lab adaptation has relaxed selection on invasion genes; with only 8 taxa this is one of the few branch-comparison tests that has adequate power.
5. **Tajima's *D* + McDonald–Kreitman on the 141 priorities, using projected Pv4 VCFs (§4.5, §4.6)** — bridges the strain-level (8) and cohort-level (thousands) data we have invested in projecting; cross-references with BUSTED hits to separate ancient adaptation from ongoing selection.

## Considered but rejected

- **TreeMix / SpaceMix population graph (§3.4)** — 8 single-strain references give no within-population replicates; the appropriate scale is the Pv4 cohort, not the pangenome panel.
- **GARD-aware BUSTED on every gene without GARD running first** — circular; we must do recombination detection first (§2.1).
- **Phylogenetic Independent Contrasts** — requires continuous phenotypic traits across the 8 strains; we have none beyond geographic origin and a few categorical labels.
- **Standalone DCA on individual 8-row MSAs (§5.4 without augmentation)** — N_eff is far too low for reliable contact prediction. We list the augmented variant only.
- **Full rooted branch-site tests (§8.2) before adding an outgroup** — listed under §8 with the explicit outgroup-addition prerequisite.

## References note

DOIs verified against Oxford Academic, PubMed, and journal records as of May 2026. Tool versions reflect the most recent stable releases—ASTER v1.23, HyPhy v2.5.65, IQ-TREE v2.3.6, scikit-allel v1.3.5, topGO v2.62, RDP5 v5.51.
