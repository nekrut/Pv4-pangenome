# BLOG1.md fact-check report

Generated 2026-05-25. Claims are organized by the categories in the verification request.

---

## 1. Verdict table

| # | Claim (paraphrased) | Verdict | Note |
|---|---|---|---|
| **Category A — Algorithm / tool descriptions** | | | |
| A1 | PGGB pipeline order: wfmash → seqwish → smoothxg | IMPRECISE | Correct order but incomplete: PGGB also runs gfaffix after smoothxg. The PGGB README states "wfmash … seqwish … smoothxg, gfaffix." |
| A2 | wfmash uses MashMap to pre-filter candidate orthologous segments by k-mer-based mash distance | CORRECT | Confirmed: wfmash docs describe MashMap3 minmer-sketch approximate mapping as the first stage before base-level alignment. |
| A3 | Minigraph-Cactus "starts from a designated reference and adds non-reference structural variation on top" | CORRECT | Cactus pangenome.md: "Unlike Progressive Cactus, Minigraph-Cactus does depend on a predetermined reference genome." |
| A4 | Minigraph (Li 2020) is "the structural-variation engine that Minigraph-Cactus is built around" | CORRECT | Minigraph-Cactus uses minigraph to build the initial SV graph, then maps assemblies back to it. The description is accurate in spirit. |
| A5 | KegAlign is "GPU-accelerated lastZ" | CORRECT | KegAlign GitHub confirms it is a GPU-accelerated wrapper that generates and optimizes lastZ commands. |
| A6 | `axtChain` + `chainStitchId` produce chains usable by CrossMap / liftOver | CORRECT | The MULTIZ.md pipeline (6-step UCSC chain build) confirms both tools are used; CrossMap accepts the resulting chain. Command names are correct UCSC kentUtils tools. |
| A7 | `cactus-update-prepare` is a real Minigraph-Cactus subcommand | CORRECT | Confirmed in the Cactus doc/ directory: cactus-update-prepare.md exists and documents an incremental-update workflow. |
| A8 | `odgi position` "translates a coordinate on any path to a node in the graph and back out to a coordinate on any other path" | CORRECT | Consistent with odgi documentation and project internal use (OUTLINE.md describes the same function). |
| A9 | `vg deconstruct` "walks the graph and emits a VCF in which every bubble is a variant record" | CORRECT | Standard vg deconstruct behavior; no contradicting evidence found. |
| A10 | PGGB `-n 8` flag means "each query segment permitted up to eight competing mappings" | WRONG | In PGGB, `-n` / `--n-haplotypes` specifies the number of haplotypes. The competing-mappings count is controlled by a separate flag (`-c` / `--n-mappings`). The description in the blog conflates two distinct parameters. See correction below. |
| **Category B — KegAlign / lastZ parameter defaults** | | | |
| B1 | KegAlign default scoring: "HoxD70 matrix, hspthresh 3000, seed 12-of-19" | CORRECT | All three confirmed against the lastZ documentation: default matrix is HOXD70, default hspthresh is 3000, default seed is 12-of-19. KegAlign inherits these lastZ defaults (no overrides found in the KegAlign README). |
| B2 | Tuned parameter set: "matrix +100/-100, hspthresh 4500, seed 14-of-22" | CORRECT | These are non-default values that tighten specificity; confirmed consistent with the OUTLINE.md tuned-run description. Values are internally documented and plausible. |
| **Category C — Assembly statistics** | | | |
| C1 | PvP01 (GCA_900093555.2): "29.0 Mb across 242 contigs" | CORRECT | Computed from PvP01.fa.fai: sum of column 2 = 29.0 Mb; row count = 242. |
| C2 | PvP01 "assembled by Auburn et al. (2016) from a Papua, Indonesia isolate" | CORRECT | PMC5172418 confirms Papua, Indonesia provenance for PvP01. |
| C3 | PvP01 assembled "using long-read sequencing" | WRONG | Auburn et al. 2016 used Illumina short-read sequencing (75 bp, 100 bp, 250 bp paired-end reads on GAII, HiSeq 2000, MiSeq, and HiSeq 2500 mate-pair libraries). No long-read data were used. See correction below. |
| C4 | Auburn et al. 2016 reference: Wellcome Open Res 2016;1:4, DOI 10.12688/wellcomeopenres.9876.1 | CORRECT | Semantic Scholar and PMC confirm: journal = Wellcome Open Research, volume 1, article 4, DOI confirmed. |
| C5 | PAM (GCA_949152365.1): "28-contig, 29.4 Mb chromosome-scale assembly" | CORRECT | PAM.fa.fai: 28 rows; total size = 29.437 Mb (~29.4 Mb). |
| C6 | PAM is "the Peruvian Amazon Pv01-19 isolate" | CORRECT | PMC10568799 confirms isolate name Pv01-19 from Iquitos, Peruvian Amazon. |
| C7 | De Meulenaere et al. (2023), BMC Genomics 2023;24(1):606, DOI 10.1186/s12864-023-09707-5 | CORRECT | Semantic Scholar confirms: title "A new Plasmodium vivax reference genome for South American isolates", BMC Genomics, 2023, volume 24, DOI correct. Article number needs spot-check (see Uncertainties). |
| **Category D — Gene IDs + chromosome assignments** | | | |
| D1 | *dhps* = PVP01_1429500 on chromosome 14 | CORRECT | GFF3 (PvP01.genbank.gff3.gz) shows PVP01_1429500 on contig LT635625.2 = chromosome 14 (16th entry in fai, 14th nuclear chromosome). |
| D2 | *dhps* encodes dihydropteroate synthase, carries sulfadoxine resistance markers | CORRECT | GFF3 name = "PPPK-DHPS" (hydroxymethyldihydropterin pyrophosphokinase-dihydropteroate synthase). DHPS is the canonical sulfadoxine target. |
| D3 | *dhfr-ts* = PVP01_0526600 on chromosome 5 | CORRECT | GFF3 shows PVP01_0526600 on contig LT635616.2 = chromosome 5 (7th entry in fai, 5th nuclear chromosome). |
| **Category E — Numerical facts about the v3 analysis** | | | |
| E1 | "256 alignment blocks across the full genome" | WRONG | The file `synteny_2way/PvP01_to_PAM_2way.blocks.tsv` has 256 lines total: 1 header row + 255 data rows. The correct count is 255 alignment blocks. |
| E2 | "1,895-sample MalariaGEN Pv4 cohort" | CORRECT | MalariaGEN Pv4 data-package page states "1,895 worldwide samples of Plasmodium vivax." |
| E3 | The v3 analysis lifts the Pv4 cohort from PvP01 onto "seven other strains" (Path A2) | CORRECT | MALARIAGEN_VCF_PROJECTION.md: "For each of the 7 non-PvP01 strain assemblies… and each of the 16 chromosomes." |
| **Category F — References** | | | |
| F1 | Paten et al. 2017 — Genome Res 27(5):665-676, DOI 10.1101/gr.214155.116 | CORRECT | Semantic Scholar: title "Genome graphs and the evolution of genome inference", Genome Research, 2017, volume 27, pages 665-676, DOI confirmed. |
| F2 | Garrison et al. 2024 — Nat Methods 21(11):2008-2012, DOI 10.1038/s41592-024-02430-3 | CORRECT | PubMed 39433878 confirms: Nature Methods, 2024, vol 21, issue 11, pages 2008-2012, DOI confirmed. |
| F3 | Hickey/Paten et al. 2024 — Nat Biotechnol 42(4):663-673, DOI 10.1038/s41587-023-01793-w | CORRECT | Semantic Scholar confirms: "Pangenome graph construction from genome alignments with Minigraph-Cactus", Nat Biotechnol, 2023→2024, vol 42, pages 663-673, DOI confirmed. |
| F4 | Auburn et al. 2016 — Wellcome Open Res 2016;1:4, DOI 10.12688/wellcomeopenres.9876.1 | CORRECT | Confirmed (see C4 above). |
| F5 | De Meulenaere et al. 2023 — BMC Genomics 2023;24(1):606, DOI 10.1186/s12864-023-09707-5 | CORRECT (DOI + journal), UNCERTAIN (article number 606 not independently confirmed from the paper itself — see Uncertainties). |
| F6 | Li et al. 2020 — Genome Biol 21(1):265, DOI 10.1186/s13059-020-02168-z | CORRECT | Semantic Scholar confirms: "The design and construction of reference pangenome graphs with minigraph", Genome Biology, 2020, vol 21, DOI confirmed. Article number 265 not confirmed from Semantic Scholar but is consistent with the BMC numbering scheme; DOI resolves. |

---

## 2. Errors that need fixing

**A10 — PGGB `-n 8` description (WRONG)**

The blog says: "including `-n 8` so that each query segment is permitted up to eight competing mappings"

This is incorrect. In PGGB, `-n` / `--n-haplotypes` specifies the number of haplotypes (genome copies), not competing mappings per segment. PGGB passes this value to smoothxg for graph normalization. The parameter that controls competing mappings is `-c` / `--n-mappings` (wfmash flag). For a two-genome warm-up, `-n 2` would be the appropriate haplotype count; `-n 8` in context means the flag was set to match the eight-way build's haplotype count, not to allow eight competing alignments per segment.

Proposed replacement: "including `-n 8` to match the haplotype count of the eight-way build—PGGB passes this figure to smoothxg for graph normalization"

**C3 — PvP01 sequencing technology (WRONG)**

The blog says: "assembled by Auburn et al. (2016)… using long-read sequencing aimed specifically at the subtelomeric pir family"

Auburn et al. 2016 used Illumina short-read sequencing (75 bp, 100 bp, and 250 bp paired-end reads on GAII/HiSeq 2000/MiSeq, plus HiSeq 2500 mate-pair libraries). No long reads were used. The assembly approach was iterative scaffolding of Illumina data, not long-read assembly.

Proposed replacement: "assembled by Auburn et al. (2016) from a Papua, Indonesia isolate using Illumina short-read sequencing with targeted effort to resolve the subtelomeric pir family"

**E1 — Alignment block count (WRONG)**

The blog says: "The two-way graph yields 256 alignment blocks across the full genome"

`PvP01_to_PAM_2way.blocks.tsv` contains 256 lines = 1 header + 255 data rows. The correct count is **255 alignment blocks**.

Proposed replacement: "The two-way graph yields 255 alignment blocks across the full genome"

---

## 3. Uncertainties

- **De Meulenaere 2023 article number**: the blog cites "BMC Genomics 2023;24(1):606." The DOI and journal are confirmed, but article number 606 could not be independently confirmed from an open full-text page (Springer redirected behind a paywall). The PMC accession (PMC10568799) is confirmed open-access; the author can verify article number 606 against the PMC or journal page directly.

- **PAM isolate name in De Meulenaere**: the blog says "Pv01-19"; PMC10568799 confirms this is correct, so this is now CORRECT, not uncertain.

- **Garrison et al. 2024 author list**: the blog cites "Garrison E, Guarracino A, Heumos S, Villani F, Bao Z, Tattini L, et al." — confirmed correct from Semantic Scholar author list order.

- **KegAlign default hspthresh and seed**: KegAlign's own documentation does not specify `hspthresh` or `seed` defaults independently; the claim that it inherits lastZ defaults is confirmed by the KegAlign GitHub description ("generates and optimizes lastZ commands"). lastZ defaults are confirmed (hspthresh 3000, seed 12-of-19). CORRECT, but technically KegAlign's documentation does not override these, so they are inherited.

- **Auburn 2016 article number / DOI version**: DOI resolves to v1 of the article (wellcomeopenres.9876.1); Wellcome Open Research articles can have versioned DOIs. The v1 DOI is what the blog cites and is the published version.

---

## 4. Suggestions

- **A1 (PGGB pipeline)**: Add gfaffix to the pipeline description for completeness: "wfmash → seqwish → smoothxg → gfaffix." The blog's omission of gfaffix is minor for a blog post audience but gfaffix is part of the canonical PGGB output.

- **A3 (Minigraph-Cactus)**: The description is accurate but could note that the primary top-level command is `cactus-pangenome`, not just "Minigraph-Cactus," to help readers navigate the software.

- **C6/C5 (PAM description)**: "28-contig" is technically correct but slightly undersells the assembly quality. The paper (PMC10568799) describes it as "14 chromosomes plus 12 unassigned contigs" — mentioning "14 chromosomes" in the blog would clarify what "chromosome-scale" means here.

- **D2 (dhps function)**: The blog's characterization as "the sulfadoxine resistance markers" is accurate for surveillance purposes; however, sulfadoxine resistance in *P. vivax* is typically described as tracking *dhps* (sulfadoxine) together with *dhfr* (pyrimethamine) for SP resistance. The blog does mention both *dhps* and *dhfr-ts* in the next sentence, so no correction is needed, but consider adding "SP (sulfadoxine-pyrimethamine) resistance" to give readers a drug name they can look up.

- **Ref [3] year**: Hickey et al. is correctly cited as Nat Biotechnol 2024, but Semantic Scholar records the original submission year as 2023 (online May 2023, print vol 42 2024). The blog's citation (Nat Biotechnol. 2024;42(4):663-673) reflects the final print year, which is the correct citation form.
