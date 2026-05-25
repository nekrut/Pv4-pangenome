# Subtelomeric microsynteny

Per-chromosome ribbon plots of gene order across all 8 strains, colored by variant-antigen family. The plots show that subtelomeric regions in *P. vivax* are dense in PIR/PHIST tandem expansions, and that the expansions vary substantially between strains — Sal-I and PvSY56 carry truncated subtelomeres relative to the curated reference set.

**Related documents.** [ORTHOLOGY.md](ORTHOLOGY.md) defines the orthogroup table that links genes across strains — every ribbon connector in these plots is a (PvP01-gene, target-gene) pair from `ortholog_table.tsv`. The Phase G family-label tagging step described in [ORTHOLOGY.md → How we built it](ORTHOLOGY.md#how-we-built-it) provides the color key. [PANGENOME.md](PANGENOME.md) and [MULTIZ.md](MULTIZ.md) are the upstream alignment layers but are not directly read by this analysis.

## How we built them

For each PvP01 chromosome's 5' and 3' end (first / last 300 kb), one ribbon plot is rendered. 14 chromosomes × 2 ends = 28 plots.

### Inputs

| Input                    | Path                                                                                                                          | What it is                                                                                                                                                                                     |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Ortholog table           | `work/03_consensus/ortholog_table.tsv`                                                                                        | The per-strain gene IDs in each orthogroup. From [ORTHOLOGY.md → Outputs](ORTHOLOGY.md#outputs)                                                                                                |
| Family table             | `work/05_families/family_table.tsv`                                                                                           | Gene → family-tag map (`PIR`, `PHIST`, `MSP`, `DBP`, `EBA`, `RBP`, `AMA`, `RAP`, `SERA`, `TRAg`, `STP1`, `RESA`, or `other`). 55,153 rows. From [ORTHOLOGY.md → Outputs](ORTHOLOGY.md#outputs) |
| Per-strain BEDs          | `inputs/annotations/{strain}.bed` (anchors), derived from `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` (non-anchors) | Gene coordinates per strain                                                                                                                                                                    |
| PvP01 chromosome lengths | hardcoded in `scripts/overnight/12_subtelomeric_microsynteny.py`                                                              | 14 nuclear chromosomes, sizes from the PvP01 reference assembly. Used to define the 5' / 3' 300 kb windows                                                                                     |
| Window size              | `SUBTEL_BP = 300_000` (hardcoded)                                                                                             | The 5' window is `[0, 300_000]` of the PvP01 chromosome; the 3' window is `[chr_len - 300_000, chr_len]`                                                                                       |
| Family colors            | Per-family palette in the script (e.g., PIR = red, PHIST = purple)                                                            | One RGB per family; genes outside any family are gray                                                                                                                                          |

### How each plot is drawn

For each (PvP01 chromosome, end ∈ {L, R}):

1. **Pick the PvP01 window.** `[0, 300 kb]` for the 5' end, `[chr_len - 300 kb, chr_len]` for the 3' end.
2. **Collect PvP01 genes in the window.** From the PvP01 BED, filter to genes whose entire CDS is inside the window.
3. **Find each PvP01 gene's orthogroup.** Look up the gene ID in `ortholog_table.tsv`.
4. **For each non-PvP01 strain, gather every ortholog gene in any of the PvP01-window orthogroups.** Use the strain's column of `ortholog_table.tsv` and the strain's BED to get the gene's coordinates. These ortholog genes may live anywhere on the strain's chromosomes — not necessarily in a 300 kb subtelomere — because subtelomeric content varies between strains.
5. **Plot.** Eight horizontal tracks, one per strain, top-to-bottom. Each gene drawn as a colored rectangle; color = family tag. PvP01 (top track) is at its real coordinates; the other 7 strains' genes are **packed sequentially** in gene-hit order with a fixed visual gene-width (~3 kb) and a 500-bp visual gap between genes. Each strain's hits are grouped by contig; contig boundaries get a 5 kb visual gap.
6. **Draw ribbons.** For every PvP01 gene in the window, a thin diagonal line connects its center to the center of its ortholog in every other strain's track, colored by family.

### What the x-axis means (two coordinate systems on one axis)

This is the most-asked question about the plots, so it deserves an explicit answer:

**PvP01 row (top track)**: real genomic position in bp, with the window's start subtracted to begin at 0. For a 5' plot it spans 0–300,000 bp; for a 3' plot it spans `chr_len - 300_000` to `chr_len`, also shifted to start at 0.

**All other 7 strain rows**: NOT real coordinates. Each gene is drawn at a fixed width (`max(1500 bp, true gene length)`), with 500 bp gap between adjacent genes and 5 kb gap between contigs. The x position is just the cumulative width of preceding orthology-hit genes — purely a gene-order layout.

So inter-gene distances on non-PvP01 tracks are meaningless — gaps don't reflect real spacing in those genomes. The ribbons are doing the actual work: matching colors at matching ortholog endpoints. Real-coordinate alignment would have made non-PvP01 tracks invisible because subtelomere lengths vary by 10× between strains.

Driver script: `scripts/overnight/12_subtelomeric_microsynteny.py` (~250 lines). Wall: ~1 min for all 28 plots.

## Outputs

| File                                                | Count | Size      | Where                                |
| --------------------------------------------------- | ----: | --------- | ------------------------------------ |
| `writeup/microsynteny/chr{1-14}_L.png`              |    14 | ~50 KB ea | **Git** — `v3/writeup/microsynteny/` |
| `writeup/microsynteny/chr{1-14}_R.png`              |    14 | ~50 KB ea | **Git** — `v3/writeup/microsynteny/` |
| `scripts/overnight/12_subtelomeric_microsynteny.py` |     1 | 8 KB      | **Git** — the driver script          |

All 28 PNGs total ~1.5 MB. The script is parameterized via per-strain BED paths and PvP01 chromosome name → length map; swap species by editing both.

## How the plots are used

This is a presentation analysis, not a downstream input — the PNGs feed into the BLOG post, the HyPhy report, and the brc-analytics PangenomeView's variant-antigen panel. There is no downstream computation that consumes them.

Two findings the plots make concrete:

- **PvW1 and PAM preserve the largest PIR cluster expansions** (50+ PIR genes in some subtelomeres). Most other strains have ≤ 20 PIR per subtelomere.
- **Sal-I and PvSY56 truncated subtelomeres.** Both isolates have fewer PIR/PHIST genes at chromosome ends, consistent with fragmented subtelomere assembly (Sal-I is the monkey-adapted lab strain; PvSY56 is the Auburn PCR-free assembly).

## Re-running on a different species

The script is closely tied to *P. vivax* variant-antigen family naming (PIR / PHIST / Pv-fam-h) and to PvP01 as the reference for the ribbon layout. To adapt:

1. Edit the `PVP01_CHRS_FULL` dict at the top of `scripts/overnight/12_subtelomeric_microsynteny.py` — replace with the new reference's chromosome accession → (chr_name, length) map.
2. Edit the `FAMILY_COLORS` dict — add / remove family tags as appropriate for the species. For *P. knowlesi*: drop `PIR`, add `SICAvar`; keep the shared families (`DBP`, `RBP`, `MSP`, `AMA`, `SERA`, `EBA`).
3. Provide a populated `work/05_families/family_table.tsv` (the Phase G family-label step described in [ORTHOLOGY.md → How we built it](ORTHOLOGY.md#how-we-built-it)).
4. Provide a populated `work/03_consensus/ortholog_table.tsv`.
5. Adjust `SUBTEL_BP` if the species has very different subtelomere extents. *P. vivax* has ~300 kb subtelomeres; *P. falciparum* has ~50-100 kb; *Toxoplasma* and other apicomplexans vary.

For the *P. knowlesi* scaffold, set `VAR_ANTIGEN_RE='SICAvar|DBP|RBP|MSP|AMA|SERA|EBA'` in `pipeline/species.conf` and run the script after the ortholog + family steps finish.

Wall: ~1 min for the 28 plots once the inputs are in place.

## Galaxy tool-wrap list

This is a single-script analysis — a thin matplotlib driver that reads three TSVs and emits PNGs. The Galaxy wrap is straightforward.

| Tool                                        | Galaxy state    | Note                                                                                                                                        |
| ------------------------------------------- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `subtelomeric_microsynteny.py` (the driver) | ⚠ needs wrap    | Single Python script with three inputs (ortholog table, family table, per-strain BEDs) and a directory output (28 PNGs). 8 KB; tiny wrapper |
| `matplotlib`, `numpy`, `pandas`             | ✅ in base image | Standard dependencies; available in any Python image                                                                                        |

Suggested workflow shape: a single tool, no subworkflow. Inputs: `ortholog_table.tsv` + `family_table.tsv` + a collection of per-strain BEDs. Output: a collection of 28 PNGs (one per chromosome × end). Category: `COMPARATIVE_GENOMICS`.

For an interactive BRC-analytics frontend, the plots can be re-implemented in D3 or another JS plotting library — the layout is regular enough (8 tracks, fixed gene widths, family-colored rectangles, family-colored ribbons) that an in-browser version reading the same three TSVs would render in real time. Not necessary for v1 — the static PNGs are sufficient. A v2 add-on could let users zoom into any single orthogroup and pull up its MSA + tree + BUSTED result via the cross-link to [MSA_HYPHY.md](MSA_HYPHY.md).
