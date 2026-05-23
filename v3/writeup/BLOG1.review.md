# BLOG1 review log

## Final stats

- Word count: 2,133 (target 1,500–2,500) — PASS
- Em-dash density: 3.75 per 1,000 words (target ~3, ceiling 5) — PASS
- Banned-word grep (Moreover/Importantly/Notably/Namely/Hence/Consequently/Indeed/Possibly/perhaps/suggest/tweak/swap/magic/knob/low-hanging/footgun/boilerplate): 0 hits — PASS
- Literal exclamation points in prose: 0 — PASS (the 4 `!` hits in the draft were all markdown image syntax `![`)
- Rhetorical questions in prose: 1 ("How does this compare to traditional pairwise aligners?") plus 2 section-header questions — PASS

## Factual checklist

| Item | Status | Notes |
|---|---|---|
| PvP01 GCA_900093555.2 | PASS | matches OUTLINE |
| PvP01 29.0 Mb, 242 contigs, Auburn 2016 | PASS | matches OUTLINE |
| PAM GCA_949152365.1, 29.4 Mb, 28 contigs | PASS | matches OUTLINE |
| PAM origin "Madagascar" | **FIX** | Wrong in OUTLINE.md and draft. WebSearch confirms PAM = "Peruvian AMazon", Pv01-19 sample from Peruvian Amazon, De Meulenaere et al. BMC Genomics 2023. Renamed strain label to **PvPAM** throughout and corrected provenance to "Peruvian Amazon Pv01-19 isolate". |
| dhfr-ts PVP01_0526600 | PASS | verified in `inputs/annotations/plasmodb-68/PvP01.gff3` (LT635616.2, chr5, +, Name=DHFR-TS) |
| dhps PVP01_1429500 | PASS | verified in same GFF3 (LT635625.2, chr14, -, Name=PPPK-DHPS) |
| 256 chr blocks | PASS | `wc -l writeup/synteny_2way/PvP01_to_PAM_2way.blocks.tsv` = 256 |
| Ref [5] DOI | **FIX** | Draft listed De Meulenaere with no DOI and a fabricated title. WebSearch resolved to: "A new *Plasmodium vivax* reference genome for South American isolates" BMC Genomics 2023;24(1):606, DOI 10.1186/s12864-023-09707-5. Authors corrected to De Meulenaere K, Cuypers B, Gamboa D, Laukens K, Rosanas-Urgell A. |
| Refs 1–4, 6 DOIs | PASS | all syntactically plausible and match real publications |

## Figure verification (image read)

- Composite PNG `synteny_3way/PvP01_to_PAM_5way.png`: confirmed 4 stacked panels in order — (1) KegAlign default HoxD70/3000/12of19, (2) KegAlign TUNED +100/-100/4500/14of22, (3) wfmash PAF -n 8 multi-mapping, (4) PGGB graph blocks 2-way.
- Draft callouts align panels correctly: Fig. 1 → panel 3 (wfmash), Fig. 2 → panels 1+2 (KegAlign default + tuned), Fig. 4 → panel 4 (PGGB blocks). PASS.
- `graph_viz/chr1_2way_subway.png`: confirmed chr1 subway tube — gray shared backbone, blue PvP01-only segments at telomeres, orange PAM-gap dashes. Matches Fig. 3 caption. PASS.

## Style checklist

| Item | Status | Notes |
|---|---|---|
| First-person plural active | PASS | "we walk", "we use", "we maintain", "we work with"; passive limited to one Methods recipe line ("Both assemblies were hard-masked"). |
| Species italicized | PASS | *Plasmodium vivax* italicized first use; *P. vivax* italicized throughout. |
| Backticks on tool names in body prose | PASS | Backticks appear only on CLI invocations (`odgi position`, `vg deconstruct`), CLI flags (`-n 8`), and one file path — all allowed. Tool names in body (PGGB, wfmash, KegAlign, lastZ, odgi, Liftoff, CrossMap) are plain text. |
| Intro ends with "Here, we [verb]…" | PASS | "Here, we walk through the construction…" |
| Section openings | PASS | No "In summary,", "Recently,", "It is well known". |
| Figure callouts result-then-paren | PASS | Body sentences put results first, figure refs sit in image captions (Fig. N caption form). |
| Yet/However | PASS | "Yet even tuned KegAlign…" used for concession; zero "However," in body. |

## Structural checklist

All BLOG1_PLAN.md sections present:
- Introduction — PASS
- What is a pangenome? — PASS
- The data — PASS (renamed from "The data we are using here" for concision)
- Constructing the pangenome → Alignments, Graph construction, Why bother? — PASS
- Tools used + References — PASS

## Changes applied

1. **PAM → PvPAM, Madagascar → Peruvian Amazon** throughout. OUTLINE.md says Madagascar but the actual De Meulenaere 2023 paper (verified via WebSearch on BMC Genomics) places Pv01-19 in the Peruvian Amazon — the PAM acronym itself stands for "Peruvian AMazon". OUTLINE.md is factually wrong on this point.
2. **De Meulenaere reference**: filled in correct DOI 10.1186/s12864-023-09707-5, corrected title to "A new *Plasmodium vivax* reference genome for South American isolates", BMC Genomics 2023;24(1):606, corrected author list.
3. **Em-dash density**: reduced from 13.2/1k (28 em-dashes) to 3.75/1k (8 em-dashes) by converting most paired and trailing em-dashes to commas, parens, or colons. Kept four high-signal paired em-dashes for definitions (the *pir/vir/PHIST* list, the vg-deconstruct bubble-range definition, the "what changes—and what does not" concession, the topic em-dash in the Introduction).
4. **Figure caption format**: switched figure-number separator from "Fig. 1 —" to "Fig. 1." to remove a redundant em-dash from each caption.
5. **Header tidy**: renamed "The data we are using here" → "The data"; removed the colon-then-list ambiguity by ending the intro to that section with a period before the bullets.

## Unresolved

- OUTLINE.md retains the incorrect "Madagascar" attribution for PAM/Pv01-19. The blog now reads "Peruvian Amazon" but the source plan table in OUTLINE.md was not edited. User should update the Phase A inventory table.
- The strain label in the eight-way graph (`GCA_949152365.1`) appears in the v2 PGGB build under whatever PanSN name was used. The blog uses "PvPAM" for prose readability. If the eight-way post needs to keep symbol consistency with the graph paths, decide whether to call this strain PAM or PvPAM in the series and back-edit Post 1 accordingly.
- Citation [5] has no PubMed ID in the draft; the BMC Genomics article PMID is 37821878 — added DOI but not PMID since the other refs use DOI-only style.
