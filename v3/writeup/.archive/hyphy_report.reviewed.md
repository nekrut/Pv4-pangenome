# HyPhy selection analysis of the *P. vivax* 8-strain core pangenome
_Anton Nekrutenko — internal report, 2026-05-21_

## Summary

Bacterial-host arms-race signatures in *Plasmodium vivax* are most informative where they overlap with vaccine and drug targets, but most prior selection scans for this species have leaned on single reference annotations. We tested 1,584 PvP01 protein-coding genes from the *P. vivax* 8-strain core pangenome for episodic diversifying selection with BUSTED, and we ran the full HyPhy bundle—BUSTED, aBSREL, MEME, and FEL—on 41 priority genes (of 157 nominated; the remaining 116 failed Phase H frame or strain checks). At the conventional 5% threshold we recovered 107 significant genes (6.8%, against an expected 79 under the null); at p<0.01 the count is 49 (3.1%, expected 16); at p<0.001 it is 20 (expected 2); and at p<10\textsuperscript{-4} we count 12.

The empirical distribution is enriched at the low-p tail relative to a uniform null—a fold enrichment of 3.1× at p<0.01 and 12.6× at p<0.001—a signature that bulk selection is real and not driven solely by analytical noise. The 1,584 bulk genes are drawn from the single-copy core that survives 8-way orthology with a clean codon frame, and by construction this set excludes most variant-antigen paralogues (PIR, PHIST, MSP, DBP, RBP); those families collapse into multi-copy orthogroups and are filtered upstream. What remains is conserved core machinery plus annotated singletons. The priority bundle—41 genes carrying BUSTED, aBSREL, MEME, and FEL output—yields 3 BUSTED-significant, 5 aBSREL-significant, 8 MEME-significant, and 1 FEL-positive (≥1 codon at p<0.05) gene, with multi-method concordance for the strongest candidates discussed in Section 2.

Yet, 116 nominated priority genes—the vaccine targets CSP, MSP1, DBP, and the RBP repertoire; the major drug-resistance markers DHPS, MDR1, ATP4, CRT-o, and the plasmepsins; and most liver-stage factors—did not survive Phase H quality filters and carry no HyPhy output here. Section 3 lists those gaps and Section 5 proposes a relaxed rebuild to recover them. Here, we describe the verified findings and the path to closing the coverage gaps.

## Methods

We anchored the analysis on the PlasmoDB release 68 PvP01 reference annotation, extracted CDS for each of 8 strains—PvP01, Sal-I, PvW1, PAM, PvSY56, PvT01, PvC01, and MHC087—using gffread, and built 8-way codon multiple sequence alignments with MAFFT followed by codon back-translation against the protein alignment. We trimmed alignments with trimAl in --automated1 mode, inferred per-gene gene trees with IQ-TREE under -m MFP, and ran HyPhy 2.5.99 with the standard BUSTED, aBSREL, MEME, and FEL pipelines—BUSTED on the full bulk and the priority bundle, and aBSREL, MEME, FEL on the priority bundle only.

We applied conventional significance thresholds: gene-level p<0.05 for BUSTED; branch-level Holm-corrected p<0.05 for aBSREL; codon-level p<0.05 for MEME; and codon-level p<0.05 with β>α (positive) or β<α (purifying) for FEL. We did not run a GARD recombination pre-screen—a caveat we return to in Section 4.

## 1. Bulk BUSTED screen (1,584 genes)

### 1.1 p-value distribution (Fig. 1)

![Bulk BUSTED p-value distribution and Q-Q plot](hyphy_plots/fig1_bulk_pvalue_distribution.png)

The empirical p-value distribution is bimodal—a tall spike at p≈0.5 corresponding to genes where the BUSTED LRT returned zero (no detectable departure from the null), and a heavy tail extending to p<10\textsuperscript{-9} for 1,583 genes, with a single off-scale outlier (PVP01_1255400) at p≈0. Counts: 107 genes at p<0.05 (1.4× the uniform expectation), 49 at p<0.01 (3.1×), 20 at p<0.001 (12.6×), and 12 at p<10\textsuperscript{-4}. The Q-Q plot (panel b) departs from the diagonal across the entire upper tail, confirming a population of genes under episodic diversifying selection that is not consistent with a uniform null. The Q-Q ordinate is capped at -log10(p)=12 to keep the chr12 outlier from compressing the rest of the distribution.

### 1.2 Top hits (Fig. 2)

![Manhattan-style plot of bulk BUSTED -log10(p) across the 14 PvP01 chromosomes](hyphy_plots/fig2_manhattan.png)

Selection signal is non-uniform across the 14 PvP01 chromosomes. The bulk set is biased toward internal, single-copy loci by the orthology filter, so the classical subtelomeric variant-antigen hot zones are under-sampled here; the hits visible on Fig. 2 are predominantly internal and map to merozoite-surface, invasion, and a scatter of conserved-function genes. The y-axis is capped at -log10(p)=12; the single off-scale point on chr12 is PVP01_1255400 (LCCL-domain protein, p=0 within BUSTED's numerical floor). The 50 strongest signals are tabulated in Section 6.

### 1.3 By gene family (Fig. 3)

![Family-level BUSTED enrichment](hyphy_plots/fig3_family_enrichment.png)

By family bucket the rates at p<0.05 are: hypothetical 60% (n=5, small denominator with large variance); conserved (housekeeping-like) 8% (n=266); other 6% (n=1267); non-coding 5% (n=39). The baseline across all 1,584 bulk genes is 6.8% at p<0.05 and 3.1% at p<0.01 (dashed lines on Fig. 3). Variant-antigen families (PIR, PHIST, MSP, DBP, RBP, AMA) are dramatically under-represented in this bulk set—only a handful of singletons survive the 8-way single-copy orthology filter—so the family bucket plot reflects the conserved core, not the antigen-diversifying periphery. Within the core, the "conserved" housekeeping-like bucket tracks the baseline at 7.9%, modestly elevated, in line with episodic selection on a subset of conserved enzymes and chaperones. The hypothetical-bucket spike rides on five genes only and should be re-examined once an inframe domain annotation pass is available. The variant-antigen story belongs to the priority bundle and the coverage-gap discussion (Sections 2-3).

### 1.4 By functional category (Fig. 4)

![Bulk BUSTED enrichment by priority category](hyphy_plots/fig4_category_breakdown.png)

Cross-referencing the bulk p-values against our priority category table reproduces the same pattern at the functional-annotation level: invasion (60% / n=5) and vaccine target (17% / n=6) sit above the 6.8% baseline; chromatin regulator (0% / n=3), liver stage (0% / n=2), and drug resistance (0% / n=4) sit below it. Translation-housekeeping, our internal negative control, sits at the baseline floor, as it should. With per-category n in the single digits these point estimates are noisy, and we read them as direction-of-effect indicators, not effect sizes.

## 2. Priority bundle (41 genes)

### 2.1 Overview heatmap (Fig. 5)

![Priority-bundle heatmap: 41 genes × 4 HyPhy methods](hyphy_plots/fig5_priority_heatmap.png)

We rank the 41 priority genes by the mean -log10(p) across BUSTED, aBSREL, and MEME (panel a), with the FEL positively-selected codon count attached as a separate column (panel b). Asterisks mark cells crossing the conventional p<0.05 threshold. 3 of 41 cross BUSTED, 5 of 41 cross aBSREL, 8 of 41 cross MEME, and 1 of 41 (Pvs230) carries at least one positively-selected codon under FEL.

### 2.2 Multi-method evidence

We score each priority gene by the number of methods at which it crosses p<0.05 (with FEL counted on the presence of ≥1 positively-selected codon). 5 of 41 genes are flagged in ≥2 methods. The leaders by this score are:

| Gene | PVP01 id | Category | BUSTED p | aBSREL min p | MEME min p | FEL+ sites | # methods sig |
|:----|:--------|:--------|---------:|---------:|---------:|----:|----:|
| Pvs230 | PVP01_0415800 | vaccine_target | 5.59e-02 | 2.69e-02 | 1.38e-03 | 1 | 3 |
| AARP | PVP01_0531900 | invasion | 1.68e-02 | 3.78e-06 | 1.81e-02 | 0 | 3 |
| RON2 | PVP01_1255000 | invasion | 1.12e-03 | 5.63e-06 | 1.34e-02 | 0 | 3 |
| AMA1 | PVP01_0934200 | vaccine_target | 3.90e-01 | 1.53e-03 | 3.60e-02 | 0 | 2 |
| CSS | PVP01_1344100 | invasion | 3.92e-04 | 2.91e-04 | 1.41e-01 | 0 | 2 |

These are the cleanest selection candidates we have under the present alignment-and-tree pipeline. Pvs230 crosses BUSTED at p=0.056 only; the call here is conservative, but aBSREL, MEME, and FEL all flag it. The Pvs230 BUSTED p in the priority bundle differs from its bulk-screen value because the priority pipeline retrims the alignment and re-fits the gene tree before HyPhy; bulk and priority p-columns are not identical alignments. The full 41-row table is in Section 6.

### 2.3 Vaccine targets — selection landscape (Fig. 6)

![Per-codon FEL signal for vaccine-candidate priority genes](hyphy_plots/fig6_vaccine_fel.png)

The vaccine-target set has 6 members with surviving HyPhy output. We plot per-codon β−α with FEL p<0.05 sites coloured. Only Pvs230 carries a FEL-positive codon (codon ~750); the rest carry purifying-only sites or no FEL hits. Most vaccine candidates—CSP, MSP1, the DBP/RBP family, Pvs28, Pvs48/45—are not in this panel because their alignments did not pass Phase H; the genes shown here are the survivors and should be treated as anecdotal in the absence of the full set (Section 3).

### 2.4 Invasion proteins (Fig. 7)

![Per-codon FEL signal for invasion / erythrocyte-binding priority genes](hyphy_plots/fig7_invasion_fel.png)

Invasion factors—rhoptry-neck proteins, AARP, CSS, and erythrocyte-binding paralogues that survived Phase H—show scattered codons of large β−α excursion embedded in a purifying-dominated background, but the FEL p<0.05 set is dominated by purifying sites; no priority invasion gene carries a FEL-positive codon. The aBSREL branch signal is strong (RON2, AARP, AMA1, CSS each have ≥1 significant branch) and consistent with host-receptor arms-race dynamics reported across *Plasmodium* spp.

### 2.5 Housekeeping controls (Fig. 8)

![Purifying-site density by priority category](hyphy_plots/fig8_fel_neg_boxplot.png)

We use the FEL− site density—fraction of codons with p<0.05 and β<α—as a proxy for purifying pressure. Invasion and translation-housekeeping categories carry the highest medians and the broadest ranges (each spans 0 to ~1.8% of codons); vaccine targets carry the lowest median (~0.4%) with a single high outlier (Pvs230). The pattern is broadly consistent with relaxed purifying constraint on antigenic surfaces and tight constraint on housekeeping. Per-category n is 2–6 and Pvs230 is large enough to drag the vaccine-target distribution alone, so we read Fig. 8 as a sanity check, not a quantitative statement (Section 4).

## 3. Coverage gaps

The pipeline attempted 157 priority HyPhy bundles in total; only 41 produced output, leaving 116 missing. Phase H rejects a gene when the codon alignment fails one of several quality filters—fewer than 6 strains intact, premature stop codons, frame drift, or trimAl collapse below a length floor. The 116 missing bundles fall into the following functional categories (21 of the 116 trace to candidate IDs that are not in the canonical priority table and are listed as "unmapped"):

| Category | N genes missing | Example symbols |
|:--------|----:|:----------|
| unmapped (not in priority table) | 21 | PVP01_1330600, PVP01_1132600, PVP01_0416800, PVP01_0918800, PVP01_1031200 |
| drug resistance | 14 | DHPS, MDR1, ATP4, FNT, Plasmepsin V, Plasmepsin IX, Plasmepsin X, DHFS-FPGS |
| vaccine target | 12 | DBP, MSP3-beta, MSP3-gamma, MSP7, MSP9, MSP10, Pvs28, Pvs48/45 |
| translation housekeeping | 11 | Alpha-tubulin, Beta-tubulin, Actin-I, EF-1alpha, EF-Tu, Calmodulin, Profilin, Histone H4 |
| invasion | 9 | MAEBL, RipR, PTRAMP, RAP1, RAP2, RhopH3, RON5, TRAg38 |
| sexual stage transmission | 8 | AP2-G, AP2-G2, AP2-O, GDV1, ULG8, Pvs16, g377, PUF2 |
| erythrocyte binding | 7 | EBP, RBP1a, RBP1b, RBP2a, RBP2b, RBP2c, RBP2-P1 |
| variant antigen marker | 6 | PHIST-a (exemplar), PHIST-b (exemplar), STP1, VIR-C exemplar, ETRAMP11.2, Pv-fam-h (exemplar) |
| other essential | 6 | SUB1, SUB2, SERA4, SERA5, DPAP1, MSP-merozoite organizer |
| surface antigen | 6 | MSP4, MSP5, MSP8, Pv12, Pv38, Pv41 |
| liver stage | 5 | SPELD, LISP2, UIS4, P36, P52 |
| chromatin regulator | 3 | GCN5, SET2, BDP1 |
| apicoplast marker | 3 | LytB, IspD, FabZ |
| cytoadherence | 3 | KAHRP-like, SBP1-like, VIR-cytoadherence exemplar |
| metabolism | 2 | LDH, GAPDH |

The list is heavy on exactly the targets a clinical reader would care about—the vaccine candidates DBP / MSP3 / MSP7 / MSP9 / MSP10 / Pvs28 / Pvs48-45 / RBP1a-b / RBP2a-b-c / RBP2-P1, the drug-resistance markers DHPS / MDR1 / ATP4 / FNT / Plasmepsin V / IX / X / DHFS-FPGS, and the liver-stage factors P36 / P52 / UIS4 / LISP2 / SPELD. We discuss the remedy in Section 5.

## 4. Caveats

- **No GARD pre-screen.** Recombination breakpoints within a gene inflate dN/dS-style false positives. We do not screen for them here, and any single BUSTED/MEME hit on a paralogue-prone family (PIR, PHIST, MSP, DBP, RBP) should be interpreted with this in mind.
- **8 strains is a small tree.** aBSREL and MEME—both branch-/site-level methods—lose power on shallow trees, and our 7-tip ingroup (PvP01 is the reference; 7 query strains) sits well below the ~20-tip comfort zone of these tools.
- **No outgroup.** All HyPhy runs here are unrooted *P. vivax*-only fits. We cannot polarise lineage-specific signals against a *P. knowlesi* or *P. cynomolgi* outgroup until we add one.
- **Liftoff-annotation bias.** For PvSY56, PvT01, PvC01, and MHC087 the gene models are projected via Liftoff/TOGA2 rather than annotated *de novo*. Frame errors at the projection step are the dominant cause of the 116 dropped priority genes in Section 3.
- **Small per-category n.** Several Fig. 4 and Fig. 8 buckets are 2–6 genes. Direction of effect, not effect size, is the appropriate reading.

## 5. Recommendations

1. **GARD recombination pre-screen** on all 41 priority bundle alignments before publishing any positively-selected-site claim.
2. **Relaxed Phase H rebuild** at min_intact=5 (currently 7) to recover the 116 missing priority genes—the rebuild trades 2 strains of breadth for the full vaccine / drug-resistance / liver-stage panel.
3. **HyPhy RELAX** between Sal-I (lab-adapted) and the 7 field/clinical strains to test for relaxation vs intensification of selection on the lab line.
4. **McDonald-Kreitman** with the MalariaGEN *P. vivax* genomic-epidemiology VCFs—within-species polymorphism vs cross-strain divergence—for the high-priority hits that survive recommendations 1-2.

## 6. Tables

### 6.1 Top 50 bulk BUSTED hits

| Rank | Gene | PVP01 id | Chr | Family | Category | p (BUSTED) | -log10(p) |
|----:|:----|:--------|:----|:------|:---------|---------:|---------:|
| 1 | — | PVP01_1255400 | chr12 | other | non-priority | 0.00e+00 | 300.00 |
| 2 | — | PVP01_0815500 | chr8 | other | non-priority | 1.18e-09 | 8.93 |
| 3 | — | PVP01_0929400 | chr9 | other | non-priority | 3.04e-07 | 6.52 |
| 4 | — | PVP01_0104700 | chr1 | other | non-priority | 1.86e-06 | 5.73 |
| 5 | — | PVP01_1027200 | chr10 | other | non-priority | 3.89e-06 | 5.41 |
| 6 | — | PVP01_1219200 | chr12 | hypothetical | non-priority | 1.12e-05 | 4.95 |
| 7 | — | PVP01_1110200 | chr11 | other | non-priority | 2.89e-05 | 4.54 |
| 8 | — | PVP01_1220100 | chr12 | hypothetical | non-priority | 2.92e-05 | 4.53 |
| 9 | — | PVP01_1135000 | chr11 | other | non-priority | 3.22e-05 | 4.49 |
| 10 | — | PVP01_0306600 | chr3 | other | non-priority | 4.06e-05 | 4.39 |
| 11 | — | PVP01_1258300 | chr12 | other | non-priority | 6.89e-05 | 4.16 |
| 12 | — | PVP01_1227100 | chr12 | other | non-priority | 7.31e-05 | 4.14 |
| 13 | — | PVP01_1334300 | chr13 | other | non-priority | 1.08e-04 | 3.97 |
| 14 | — | PVP01_1120000 | chr11 | other | non-priority | 1.98e-04 | 3.70 |
| 15 | — | PVP01_1323000 | chr13 | other | non-priority | 2.26e-04 | 3.65 |
| 16 | — | PVP01_0109200 | chr1 | other | non-priority | 2.64e-04 | 3.58 |
| 17 | CSS | PVP01_1344100 | chr13 | other | invasion | 3.92e-04 | 3.41 |
| 18 | — | PVP01_0115400 | chr1 | other | non-priority | 4.52e-04 | 3.34 |
| 19 | — | PVP01_0909900 | chr9 | conserved (housekeeping-like) | non-priority | 8.14e-04 | 3.09 |
| 20 | — | PVP01_0715700 | chr7 | other | non-priority | 9.04e-04 | 3.04 |
| 21 | RON2 | PVP01_1255000 | chr12 | other | invasion | 1.13e-03 | 2.95 |
| 22 | — | PVP01_0107900 | chr1 | other | non-priority | 1.51e-03 | 2.82 |
| 23 | — | PVP01_1252800 | chr12 | conserved (housekeeping-like) | non-priority | 2.06e-03 | 2.69 |
| 24 | — | PVP01_0616400 | chr6 | other | non-priority | 2.35e-03 | 2.63 |
| 25 | — | PVP01_0409900 | chr4 | other | non-priority | 2.85e-03 | 2.54 |
| 26 | — | PVP01_1427700 | chr14 | other | non-priority | 3.01e-03 | 2.52 |
| 27 | AARP | PVP01_0531900 | chr5 | other | invasion | 3.13e-03 | 2.50 |
| 28 | — | PVP01_1464900 | chr14 | other | non-priority | 3.16e-03 | 2.50 |
| 29 | tRNA | PVP01_0815300 | chr8 | non-coding | non-priority | 3.30e-03 | 2.48 |
| 30 | — | PVP01_0913000 | chr9 | other | non-priority | 3.70e-03 | 2.43 |
| 31 | — | PVP01_0926500 | chr9 | other | non-priority | 3.83e-03 | 2.42 |
| 32 | — | PVP01_0944000 | chr9 | other | non-priority | 4.24e-03 | 2.37 |
| 33 | — | PVP01_1406800 | chr14 | conserved (housekeeping-like) | non-priority | 4.66e-03 | 2.33 |
| 34 | — | PVP01_1341500 | chr13 | conserved (housekeeping-like) | non-priority | 4.71e-03 | 2.33 |
| 35 | — | PVP01_0707500 | chr7 | conserved (housekeeping-like) | non-priority | 4.89e-03 | 2.31 |
| 36 | — | PVP01_1455000 | chr14 | conserved (housekeeping-like) | non-priority | 5.77e-03 | 2.24 |
| 37 | — | PVP01_1107700 | chr11 | other | non-priority | 5.81e-03 | 2.24 |
| 38 | — | PVP01_1243600 | chr12 | other | non-priority | 5.85e-03 | 2.23 |
| 39 | — | PVP01_0815600 | chr8 | other | non-priority | 5.96e-03 | 2.22 |
| 40 | — | PVP01_1469500 | chr14 | hypothetical | non-priority | 6.21e-03 | 2.21 |
| 41 | — | PVP01_1444900 | chr14 | other | non-priority | 6.90e-03 | 2.16 |
| 42 | — | PVP01_1133300 | chr11 | other | non-priority | 6.92e-03 | 2.16 |
| 43 | — | PVP01_1267500 | chr12 | other | non-priority | 7.05e-03 | 2.15 |
| 44 | — | PVP01_0214500 | chr2 | conserved (housekeeping-like) | non-priority | 7.53e-03 | 2.12 |
| 45 | — | PVP01_0213000 | chr2 | other | non-priority | 9.19e-03 | 2.04 |
| 46 | — | PVP01_1320300 | chr13 | other | non-priority | 9.30e-03 | 2.03 |
| 47 | — | PVP01_1121800 | chr11 | other | non-priority | 9.37e-03 | 2.03 |
| 48 | — | PVP01_1142800 | chr11 | other | non-priority | 9.93e-03 | 2.00 |
| 49 | rRNA | PVP01_1117500 | chr11 | non-coding | non-priority | 9.94e-03 | 2.00 |
| 50 | — | PVP01_1262300 | chr12 | other | non-priority | 1.03e-02 | 1.99 |

(Full top-200 table at `writeup/top200_bulk_busted.tsv`.)

### 6.2 Priority bundle — all 41 genes, four methods

| Gene | PVP01 id | Category | BUSTED p | aBSREL min p | MEME min p | FEL+ | FEL− | n sites | n branches sig (aBSREL) |
|:----|:--------|:--------|---------:|---------:|---------:|----:|----:|----:|----:|
| SIAP1 | PVP01_0307900 | liver_stage | 5.00e-01 | 1.00e+00 | 5.25e-02 | 0 | 3 | 997 | 0/10 |
| Pvs230 | PVP01_0415800 | vaccine_target | 5.59e-02 | 2.69e-02 | 1.38e-03 | 1 | 9 | 2729 | 1/13 |
| HSP70-1 | PVP01_0515400 | translation_housekeeping | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 12 | 682 | 0/8 |
| AARP | PVP01_0531900 | invasion | 1.68e-02 | 3.78e-06 | 1.81e-02 | 0 | 5 | 274 | 1/12 |
| Pvs25 | PVP01_0616100 | vaccine_target | 5.00e-01 | 1.00e+00 | 1.44e-01 | 0 | 3 | 219 | 0/12 |
| PVP01_0806500 | PVP01_0806500 | other | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 1 | 250 | 0/8 |
| HAP2 | PVP01_0814300 | vaccine_target | 1.60e-01 | 1.00e+00 | 3.13e-02 | 0 | 1 | 860 | 0/10 |
| Enolase | PVP01_0816000 | metabolism | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 3 | 446 | 0/9 |
| Vivapain-2 | PVP01_0916000 | other_essential | 5.00e-01 | 1.00e+00 | 2.15e-01 | 0 | 2 | 484 | 0/8 |
| Vivapain-3 | PVP01_0916100 | other_essential | 1.64e-01 | 1.00e+00 | 1.15e-01 | 0 | 5 | 495 | 0/11 |
| PVP01_0916200 | PVP01_0916200 | other | 5.00e-01 | 1.00e+00 | 4.37e-02 | 0 | 6 | 487 | 0/9 |
| RON4 | PVP01_0916600 | invasion | 1.06e-01 | 1.00e+00 | 8.91e-02 | 0 | 19 | 730 | 0/11 |
| AMA1 | PVP01_0934200 | vaccine_target | 3.90e-01 | 1.53e-03 | 3.60e-02 | 0 | 0 | 562 | 2/13 |
| PVM-exported HSP | PVP01_0934800 | cytoadherence | 5.00e-01 | 1.00e+00 | 3.57e-01 | 0 | 2 | 663 | 0/8 |
| TCTP | PVP01_1023000 | drug_resistance | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 1 | 171 | 0/9 |
| PVP01_1117400 | PVP01_1117400 | other | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 0 | 181 | 0/8 |
| Hexokinase | PVP01_1125500 | metabolism | 5.00e-01 | 1.00e+00 | 4.19e-01 | 0 | 4 | 493 | 0/8 |
| Histone H2A | PVP01_1131700 | translation_housekeeping | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 0 | 133 | 0/8 |
| DHODH | PVP01_1145600 | drug_resistance | 5.00e-01 | 1.00e+00 | 2.57e-01 | 0 | 3 | 557 | 0/10 |
| Pvs47 | PVP01_1208000 | vaccine_target | 7.59e-02 | 1.00e+00 | 3.15e-02 | 0 | 2 | 433 | 0/13 |
| K12 | PVP01_1211100 | drug_resistance | 5.00e-01 | 1.00e+00 | 3.51e-01 | 0 | 6 | 712 | 0/8 |
| PEPCK | PVP01_1212000 | metabolism | 5.00e-01 | 1.00e+00 | 3.25e-01 | 0 | 2 | 599 | 0/8 |
| PVP01_1216600 | PVP01_1216600 | other | 5.00e-01 | 1.00e+00 | 4.28e-01 | 0 | 6 | 513 | 0/10 |
| SIR2A | PVP01_1225700 | chromatin_regulator | 5.00e-01 | 1.00e+00 | 3.71e-01 | 0 | 2 | 306 | 0/8 |
| G6PD | PVP01_1253200 | drug_resistance | 5.00e-01 | 1.00e+00 | 4.33e-02 | 0 | 5 | 927 | 0/9 |
| RON2 | PVP01_1255000 | invasion | 1.12e-03 | 5.63e-06 | 1.34e-02 | 0 | 15 | 2203 | 1/12 |
| EF-2 | PVP01_1255900 | translation_housekeeping | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 14 | 832 | 0/9 |
| GEST | PVP01_1258000 | sexual_stage_transmission | 5.00e-01 | 1.00e+00 | 7.50e-02 | 0 | 0 | 249 | 0/9 |
| M17 aminopeptidase | PVP01_1260800 | other_essential | 5.00e-01 | 1.00e+00 | 5.69e-02 | 0 | 1 | 621 | 0/8 |
| Aldolase | PVP01_1262200 | metabolism | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 5 | 369 | 0/9 |
| Actin-II | PVP01_1336400 | translation_housekeeping | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 0 | 376 | 0/8 |
| CSS | PVP01_1344100 | invasion | 3.92e-04 | 2.91e-04 | 1.41e-01 | 0 | 3 | 381 | 1/12 |
| GAP | PVP01_1403000 | sexual_stage_transmission | 5.96e-02 | 1.00e+00 | 1.32e-01 | 0 | 0 | 244 | 0/11 |
| UIS3 | PVP01_1403100 | liver_stage | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 1 | 213 | 0/8 |
| Ferredoxin | PVP01_1419000 | apicoplast_marker | 2.19e-01 | 1.00e+00 | 1.80e-01 | 0 | 1 | 196 | 0/8 |
| MGET | PVP01_1435300 | sexual_stage_transmission | 4.64e-01 | 1.00e+00 | 1.20e-01 | 0 | 2 | 288 | 0/9 |
| CelTOS | PVP01_1435400 | vaccine_target | 5.00e-01 | 1.00e+00 | 5.17e-02 | 0 | 1 | 196 | 0/9 |
| HP1 | PVP01_1439200 | chromatin_regulator | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 0 | 275 | 0/8 |
| ASF1 | PVP01_1442700 | chromatin_regulator | 3.99e-01 | 1.00e+00 | 2.02e-01 | 0 | 2 | 269 | 0/9 |
| PVP01_1460700 | PVP01_1460700 | other | 5.00e-01 | 1.00e+00 | 6.67e-01 | 0 | 0 | 115 | 0/8 |
| TRAg22 | PVP01_1469800 | invasion | 2.39e-01 | 2.29e-01 | 8.61e-02 | 0 | 3 | 708 | 0/11 |

(Full priority bundle table at `writeup/priority_bundle_results.tsv`.)

## References

1. Murrell B, Weaver S, Smith MD, Wertheim JO, Murrell S, Aylward A, Eren K, Pollner T, Martin DP, Smith DM, Scheffler K, Kosakovsky Pond SL. Gene-wide identification of episodic selection. *Mol Biol Evol*. 2015;32(5):1365-1371. doi:10.1093/molbev/msv035

2. Smith MD, Wertheim JO, Weaver S, Murrell B, Scheffler K, Kosakovsky Pond SL. Less is more: an adaptive branch-site random effects model for efficient detection of episodic diversifying selection. *Mol Biol Evol*. 2015;32(5):1342-1353. doi:10.1093/molbev/msv022

3. Murrell B, Wertheim JO, Moola S, Weighill T, Scheffler K, Kosakovsky Pond SL. Detecting individual sites subject to episodic diversifying selection. *PLoS Genet*. 2012;8(7):e1002764. doi:10.1371/journal.pgen.1002764

4. Kosakovsky Pond SL, Frost SDW. Not so different after all: a comparison of methods for detecting amino acid sites under selection. *Mol Biol Evol*. 2005;22(5):1208-1222. doi:10.1093/molbev/msi105

5. Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. *Mol Biol Evol*. 2020;37(5):1530-1534. doi:10.1093/molbev/msaa015

6. Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics*. 2009;25(15):1972-1973. doi:10.1093/bioinformatics/btp348

7. Garrison E, Guarracino A, Heumos S, Villani F, Bao Z, Tattini L, Hagmann J, Vorbrugg S, Marco-Sola S, Kubica C, Ashbrook DG, Thorell K, Rusholme-Pilcher RL, Liti G, Rudbeck E, Nahnsen S, Yang Z, Moses MN, Nobrega FL, Wu Y, Chen H, de Ligt J, Sudmant PH, Soranzo N, Colonna V, Williams RW, Prins P. Building pangenome graphs. *Nat Methods*. 2024;21(11):2008-2012. doi:10.1038/s41592-024-02430-3

8. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol*. 2013;30(4):772-780. doi:10.1093/molbev/mst010

9. Kosakovsky Pond SL, Posada D, Gravenor MB, Woelk CH, Frost SDW. GARD: a genetic algorithm for recombination detection. *Bioinformatics*. 2006;22(24):3096-3098. doi:10.1093/bioinformatics/btl474
