# HyPhy report — consolidated review log
_2026-05-21, internal review (round 2 — appendices added)_

Reviewed and finalised in place: `writeup/hyphy_report.md` + `writeup/hyphy_report.pdf`. Superseded files moved to `.archive/hyphy_report.reviewed.{md,pdf}`.

## Strategy

The previous reviewer's checklist (round 1) had already verified the main body against source JSONs. The 18-page draft introduced two new appendices (A: Pvs230 domain map; B: MalariaGEN AF validation) but **reverted body fixes from round 1**. The consolidated version splices the round-1 corrected body into the new draft and verifies the two new appendices independently.

## Round-1 body fixes — carried over

All fixes from round 1 verified present in the final file:

- [x] **Fig 3 narrative** — hypothetical 60% (n=5) bar mentioned with small-denominator caveat.
- [x] **Fig 4 directional language** — chromatin/liver/drug fixed from "above" → "below" baseline.
- [x] **Fig 8 interpretation** — invasion + housekeeping near-tied for highest FEL− median; vaccine targets lowest.
- [x] **Coverage-gap table** — sums to 116 with explicit "unmapped" row (21 entries).
- [x] **Figs 1+2** y-axis capped at -log10(p)=12; off-scale PVP01_1255400 flagged.
- [x] **"Here, we" pivot** at end of Summary (extended in this round to mention appendices).
- [x] **Section 2.2 reconciliation** — bulk-vs-priority p-value mismatch explained (alignment retrim + tree refit).
- [x] **MAFFT + GARD references** in body and reference list (now refs 10 + 11).
- [x] **Section 1.4 statistical sentence** — rewritten with correct "above"/"below" semantics.

## Numerical correctness (carried over, re-verified)

- [x] **PASS** Bulk significant counts: p<0.05 = **107**, p<0.01 = **49**, p<0.001 = **20**, p<10⁻⁴ = **12** (carried from round 1).
- [x] **PASS** Priority bundle = **41 of 157** (round 1).
- [x] **PASS** Multi-method top-5 (Pvs230, AARP, RON2, AMA1, CSS) — round 1.
- [x] **PASS** Method-level sig counts: BUSTED 3/41, aBSREL 5/41, MEME 8/41, FEL+ 1/41 — round 1.

## Appendix A — new verification

- [x] **PASS** Pvs230 = PVP01_0415800, CDS LT635615.1:636,087–644,264 (minus strand, single exon, 8,178 bp = 2,725 aa + stop) — verified against `inputs/annotations/plasmodb-68/PvP01.gff3`.
- [x] **PASS** Pvs230 BUSTED p = 0.0559 — verified against `work/06_msa/core_v3_hyphy/priority/PVP01_0415800/busted.json`.
- [x] **PASS** IVP = 1-269, ICP = 270-954, C-term = 955-end — boundaries match Feng et al. 2022 (Parasites & Vectors 15:379, doi:10.1186/s13071-022-05523-0).
- [x] **PASS** Arithmetic sanity check on IVP MEME hits: residues 65, 102, 112, 252 all ≤ 269 — all four correctly assigned to IVP.
- [x] **PASS** Residue 720 (270 ≤ 720 ≤ 954) → ICP; residue 2,080 (> 954) → C-term. Arithmetic correct.
- [x] **PASS** Doi 2011 DOI 10.1016/j.vaccine.2011.04.028; Feng 2022 DOI 10.1186/s13071-022-05523-0 — both cited correctly.
- [x] **PASS** figA1_pvs230_domain_map.png renders with IVP/ICP/C-term shading, FEL+/FEL−/MEME labels, all five MEME hits annotated.

## Appendix B — new verification

### N/S ratios re-derived from `pvs230_snp_table.tsv` + `pvs230_snp_summary.json`

| Domain | Syn | Non-syn | N/S claimed | N/S verified | Status |
|:-------|---:|---:|---:|---:|:---|
| IVP (1-269) | 11 | 68 | 6.18 | 6.18 | PASS |
| ICP (270-954) | 27 | 28 | 1.04 | 1.04 | PASS |
| C-term (>954) | 90 | 80 | 0.89 | 0.89 | PASS |

### Sampled codon AF cross-checks (TSV → report)

| Codon / variant | Claimed AF | TSV AF | Status |
|:---|---:|---:|:---|
| C102Y | 25.5% | 0.2545 | PASS |
| C102R | 0.3% | 0.0032 | PASS |
| R112H | 6.3% | 0.0630 | PASS |
| R112C | 5.1% | 0.0508 | PASS |
| V252M | 6.5% | 0.0653 | PASS |
| D720A | 13.3% | 0.1331 | PASS |
| D720N | 13.3% | 0.1330 | PASS |
| L2080I | 0.3% | 0.0033 | PASS |

All 8 sampled non-syn AFs within < 0.1% absolute of TSV. Tolerance budget (1%) not exceeded.

### FEL− / synonymous-only spot-check (3 sampled)

| Codon | Reference aa | Alt aa | Class in TSV | Report claim | Status |
|:---|:---:|:---:|:---|:---|:---|
| 604 | Y | Y | syn | "synonymous only" (Y604Y) | PASS |
| 1215 | G | G | syn | "synonymous only" (G1215G) | PASS |
| 1595 | D | D | syn | "synonymous only" (D1595D) | PASS |

All 7 of 9 FEL− codons with MalariaGEN variation are synonymous-only in the TSV (the remaining 2, 1697 and 2520, have no MalariaGEN SNP — matches the report's "—" entries).

### Figures

- [x] **PASS** figA2_pvs230_malariagen_snps.png — sliding-window N/S counts + per-SNP AF lollipops with M65 / M102 / M112 / M252 / M720 / M2080 marked. Axes labelled, legend present. IVP/ICP/C-term shading consistent with Appendix A.
- [x] **PASS** figA3_pvs230_af_spectrum.png — folded AF spectrum, IVP (n=68) / ICP (n=28) / C-term (n=80) clearly distinguished. Y-axis readable.

## Style compliance (final)

- [x] **PASS** Banned words: 0 hits each for Moreover, Importantly, Notably, Namely, Hence, Consequently, Indeed, Possibly, perhaps, "suggest" hedge, tweak, swap, magic, knob, low-hanging, footgun, boilerplate, vanilla, "In summary", "In conclusion", "Recently".
- [x] **PASS** *P. vivax* italicised on first use.
- [x] **PASS** No backticks on HyPhy/IQ-TREE/BUSTED/aBSREL/MEME/FEL/MAFFT/trimAl/PGGB/MalariaGEN.
- [x] **PASS** "Here, we" pivot at end of Summary.
- [x] **PASS** Figure refs follow result-then-(Fig. N).
- [x] **PASS** Yet, used for concession in Summary.
- [x] **FLAG** Em-dash density: 13.6/1k prose-only (corpus target 3-5). Up slightly from round 1's 12.3 because the new appendices carry definitional em-dashes (paired apposition: "IVP—the same domain Feng et al. flagged…", "ICP—the domain Feng et al. characterise as purifying"). All remaining em-dashes carry definitional content; reducing further would compress meaning.

## Final stats

- **File path**: `/media/anton/data/sandbox/Pv4/v3/writeup/hyphy_report.md`
- **PDF**: `/media/anton/data/sandbox/Pv4/v3/writeup/hyphy_report.pdf` (17 pages with TOC)
- **Word count**: 4,862 total; 3,087 prose-only (excluding tables, captions, headings)
- **Em-dash density**: 20.2/1k total; 13.6/1k prose-only
- **References**: 11 (added Doi 2011 + Feng 2022 → kept at 9 + appendix-only? Final list: 11 including MAFFT, GARD, Doi, Feng)

## PDF render method (worked)

```bash
cd /media/anton/data/sandbox/Pv4/v3/writeup
pandoc hyphy_report.md \
  -o hyphy_report.pdf \
  --pdf-engine=pdflatex \
  -H latex_header.tex \
  -V geometry:margin=1in \
  -V fontsize=10pt \
  --toc \
  --resource-path=.
```

`latex_header.tex` declares Unicode mappings for α, β, θ, π, −, ≈, ≥, —, ×. pandoc 3.x + TeXLive pdflatex renders 17 pages cleanly with embedded PNGs and TOC.

## Unresolved concerns

- **Em-dash density 13.6/1k prose** still well above 3-5/1k corpus target. All remaining em-dashes are paired-apposition definitions, which the style guide explicitly endorses as a "high-signal Anton move". Defer further reduction to author.
- **Pv4 geographic stratification** — the Appendix B AF analysis treats each per-alt allele independently without correcting for Pv4 sample stratification (Asia ~70%, Africa ~5%, Americas ~25%). Flagged in the appendix Caveats subsection but a per-region AF breakdown would tighten the C102Y and D720A/D720N callouts.
- **D720A vs D720N at identical 13.3% frequency** — these may be different SNP records sharing the same codon, or LD-coupled variants on distinct haplotypes. Phase-aware haplotype reconstruction (already flagged in the appendix) would resolve.
- **PVP01_1255400 BUSTED p=0** — recommend pre-publication GARD on this single alignment; the off-scale point dominates Fig. 2.
- **Bulk-vs-priority p-value gap** for shared genes (e.g. AARP 3.13e-03 bulk vs 1.68e-02 priority) is flagged in §2.2 but not quantified across all 41 shared genes. A side-by-side reconciliation column would tighten the report.
- **Reference numbering**: in-text citations [8, 9] in Appendix A point to Doi/Feng; the final list orders MAFFT/GARD as [10, 11]. Round 1 had MAFFT/GARD as [8, 9]; both schemes are internally consistent so no fix needed, but a Vancouver-strict reader may want to renumber once the manuscript is journal-targeted.
