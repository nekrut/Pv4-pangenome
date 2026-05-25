# Per-orthogroup MSAs, ML trees, and HyPhy selection scans

Three sequential per-gene analyses keyed off the consensus ortholog table: codon + protein MSAs, IQ-TREE ML trees, and HyPhy BUSTED tests for positive selection. Together they produce ~5,800 per-orthogroup result bundles that drive the BRC selection tracks and the HyPhy report.

**Related documents.** [ORTHOLOGY.md](ORTHOLOGY.md) defines the orthogroup table that selects which loci enter this pipeline — every MSA / tree / BUSTED run keys off a row in `ortholog_table.tsv`; [MULTIZ.md → How we built them](MULTIZ.md#how-we-built-them) covers the chains that the upstream Stream 1 annotation projection uses to derive per-strain CDS coordinates for non-PlasmoDB strains; [PANGENOME.md](PANGENOME.md) is unrelated to this pipeline but is the sibling structural-variation layer.

## How we built them

Three sequential stages, all per-orthogroup. Each stage's output is the next stage's input.

### Stage F — codon and protein MSAs

**Inputs**:

| Input             | Path                                                                                                                    | What it is                                                                                                                                                                  |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Ortholog table    | `work/03_consensus/ortholog_table.tsv`                                                                                  | 7,504-row consensus table from [ORTHOLOGY.md → How we built it](ORTHOLOGY.md#how-we-built-it); column `n_strains` is the filter key                                         |
| Per-strain GFF    | `inputs/annotations/{strain}.fixed.gff3` (anchors), `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` (non-anchors) | Source of CDS coordinates for each strain                                                                                                                                   |
| Per-strain FASTA  | `genomes/softmasked/{strain}.fa`                                                                                        | Soft-masked assembly FASTA (the sequence the CDS coordinates point into)                                                                                                    |
| `min_intact` gate | `7` for the **strict** set, `5` for the **relaxed** set                                                                 | The minimum number of strains that must have the orthogroup intact (single-copy, no frameshift) for the OG to qualify for MSA. Filters 7,504 → 1,584 strict / 4,215 relaxed |

For each qualifying orthogroup:

1. **Extract CDS per strain.** `gffread -x /tmp/{strain}.{og}.cds.fa -g genomes/softmasked/{strain}.fa <gff filtered to this orthogroup's gene>`.
2. **Strip internal stop codons.** `pal2nal` drops strains with mid-CDS `TAA`/`TAG`/`TGA`; `build_msa.py` replaces those mid-sequence with `NNN` before passing the CDS to pal2nal. Without this fix, pseudogenes drop out of every alignment they should be in.
3. **Translate.** `gffread -y /tmp/{og}.aa.fa -g {fasta} <gff>` produces the protein sequences.
4. **Protein MSA.** `mafft --localpair --maxiterate 1000 --thread 2 /tmp/{og}.aa.fa > {OG}.protein.aln.fa`. LINSI mode is ~10× slower than default `mafft` but produces alignments suitable for selection analysis — default MAFFT produces gappy garbage on divergent paralog families.
5. **Backtranslate to codons.** `pal2nal.pl {OG}.protein.aln.fa /tmp/{og}.cds.fa -output fasta > {OG}.codon.aln.fa`.
6. **Verify.** Codon MSA length must be exactly 3 × protein MSA length. If not, the orthogroup is rejected (logged in the summary TSV).
7. **Clean variant (sibling output).** `trimAl -in {OG}.codon.aln.fa -out {OG}.codon.cleaned.fa -automated1` and the protein analog, written to the `_clean/` sibling directory.

Two sets emitted in parallel: `core_v3/` (`min_intact ≥ 7`) and `core_relaxed/` (`min_intact ≥ 5`).

Driver script: `scripts/build_8way_msa_v2.py`. Wall: ~2 hours strict, ~5 hours relaxed (sharded across all cores via the per-OG worker pool).

### Stage G — IQ-TREE ML trees

**Inputs**:

| Input              | Path                                     | What it is                                                                                                                                                                   |
| ------------------ | ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Codon MSA          | `work/06_msa/{set}/{PVP01}.codon.aln.fa` | From Stage F. The codon alignment is preferred over protein because it carries more sites and matches the BUSTED input format (Stage H requires concordant tree + codon MSA) |
| Substitution model | `MFP` (ModelFinder Plus auto-selection)  | IQ-TREE picks per-OG. For large panels (> 10k OGs), `GTR+I+G` is the recommended manual override to skip ModelFinder (~3× faster)                                            |
| Bootstrap reps     | `1000` (ultrafast bootstrap, `-B 1000`)  | Hard fallback to no-bootstrap when the alignment has < 4 unique sequences; IQ-TREE hangs at 100 % CPU with `-B 1000` on monomorphic alignments                               |
| Threads per tree   | `2` (`-T 2`)                             | Combined with `xargs -P $((N_CORES/2))` to saturate without oversubscribing                                                                                                  |

Per-orthogroup invocation:

```
n_unique=$(awk '/^>/{next} {print}' $msa | sort -u | wc -l)
if [[ $n_unique -ge 4 ]]; then
  iqtree3 -s $msa -m MFP -B 1000 -T 2 -pre $outdir/$gene
else
  iqtree3 -s $msa -m MFP -T 2 -pre $outdir/$gene
fi
```

Each per-OG workdir contains the `.treefile` (Newick — the only `*`-essential output), `.iqtree` report, `.log`, `.ckp.gz` checkpoint, `.bionj` and `.mldist` (initial trees), `.contree` (consensus), and `.splits.nex`. For 5,800 OGs the per-OG workdirs sum to ~30 GB — we tar+gz the workdirs into single archives per set after the runs finish.

Driver scripts: `scripts/overnight/04_iqtree.sh` + `04b_iqtree_retry.sh`. Wall: ~2 hours strict, ~6 hours relaxed.

### Stage H — HyPhy BUSTED selection scan

**Inputs**:

| Input               | Path                                               | What it is                                                                                                                                                     |
| ------------------- | -------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Codon MSA           | `work/06_msa/{set}/{PVP01}.codon.aln.fa`           | From Stage F                                                                                                                                                   |
| Tree                | `work/06_msa/{set}_trees/{PVP01}/{PVP01}.treefile` | From Stage G. If the tree is missing (per-OG IQ-TREE failed), HyPhy silently produces an empty JSON — we validate post-run by parsing the `"test results"` key |
| Site rate variation | `--srv No`                                         | Disable site-to-site rate variation (faster + more interpretable; the v3 setting)                                                                              |
| Branch set          | `--branches All`                                   | Test selection on the whole tree. For branch-specific selection, use `--branches Internal` or a labeled tree                                                   |

Per-orthogroup invocation:

```
hyphy busted \
  --alignment $msa \
  --tree $tree \
  --output $outdir/busted.json \
  --srv No --branches All \
  > $outdir/busted.log 2>&1
```

Driver script: `scripts/overnight/06_hyphy_bulk.sh`. Wall: ~1.5 hours strict, ~4 hours relaxed (well-parallelized — BUSTED on ~600 codon sites takes ~30 s; 5,800 jobs on 32 cores ≈ 90 min).

### Multi-nucleotide-substitution (MNS) follow-up

Single-hit BUSTED is biased by codon-spanning multi-nucleotide substitutions — what looks like "two parallel diversifying mutations" can actually be one MNS event. For top hits with a tight cluster of close-by alleles in the MalariaGEN cohort, re-run BUSTED with `--multiple-hits Double` and compare p-values.

**Pvs230 worked example.** Bulk BUSTED on `PVP01_0415800` (Pvs230, transmission-blocking vaccine candidate) reported p = 0.056 with the D720A and D720N alleles at 13.3 % combined in MalariaGEN. Re-running with `--multiple-hits Double` shifted p to 0.500 — the "two parallel diversifying alleles" was actually a single MNS event (D720T → reaching A or N depending on which codon position flipped). Output: `work/06_msa/core_v3_hyphy/priority_mh/PVP01_0415800/busted_mh.json`.

The MNS-aware re-test should be the default for any top-hit interpretation before publication.

## Outputs

Per-OG outputs live under `work/06_msa/`. The bulk-run sets are tar+gzipped after completion (5,800 per-OG workdirs in a single git operation is infeasible).

| Set                                        | Stage      | What                                                                          | Files | Raw size | Where                                                                                                       |
| ------------------------------------------ | ---------- | ----------------------------------------------------------------------------- | ----: | -------: | ----------------------------------------------------------------------------------------------------------- |
| `core_v3/`                                 | F          | Strict codon + protein MSAs (`min_intact ≥ 7`)                                | 3,168 |    29 MB | **Git** — `work/06_msa/core_v3_archive.tar.gz`                                                              |
| `core_v3_clean/`                           | F          | trimAl-cleaned strict MSAs                                                    | 3,168 |    20 MB | **Git** — `work/06_msa/core_v3_clean_archive.tar.gz`                                                        |
| `core_relaxed/`                            | F          | Relaxed codon + protein MSAs (`min_intact ≥ 5`)                               | 8,430 |   100 MB | **Git** — `work/06_msa/core_relaxed_archive.tar.gz`                                                         |
| `core_relaxed_clean/`                      | F          | trimAl-cleaned relaxed MSAs                                                   | 8,430 |    70 MB | **Git** — `work/06_msa/core_relaxed_clean_archive.tar.gz`                                                   |
| `core_v3_trees/`                           | G          | Strict IQ-TREE per-OG workdirs (treefile + .iqtree + .log + .ckp.gz + .bionj) | 1,584 |   177 MB | **Git** — `work/06_msa/core_v3_trees_archive.tar.gz` (44 MB compressed)                                     |
| `core_relaxed_trees/`                      | G          | Relaxed IQ-TREE per-OG workdirs                                               | 4,215 |   474 MB | **Dropbox** — `Pv4_v3/work_archives/core_relaxed_trees_archive.tar.gz` (119 MB compressed; too big for Git) |
| `core_v3_hyphy/bulk/`                      | H          | Strict BUSTED JSON per orthogroup                                             | 1,584 |   380 MB | **Git** — `work/06_msa/core_v3_hyphy_archive.tar.gz` (43 MB compressed)                                     |
| `core_relaxed_hyphy/`                      | H          | Relaxed BUSTED JSON per orthogroup                                            | 4,215 |   1.5 GB | **Dropbox** — `Pv4_v3/work_archives/core_relaxed_hyphy_archive.tar.gz` (162 MB compressed; too big for Git) |
| `core_relaxed_hyphy_va/`                   | H          | Variant-antigen-only BUSTED subset (PIR / PHIST / MSP / DBP / EBA …)          |   245 |    26 MB | **Git** — `work/06_msa/core_relaxed_hyphy_va_archive.tar.gz`                                                |
| `core_v3_hyphy/priority_mh/PVP01_0415800/` | H          | Pvs230 BUSTED-MH follow-up                                                    |     1 |    <1 MB | **Git** — `work/06_msa/core_v3_hyphy_archive.tar.gz` (inside the strict bundle)                             |
| `writeup/hyphy_plots/`                     | downstream | Manhattan + family-stratified plots, top-hits panels                          |    12 |     8 MB | **Git** — `writeup/hyphy_plots/`                                                                            |
| `writeup/hyphy_report.{md,pdf}`            | downstream | 18-page report with top-hits table + Pvs230 Appendix A/B                      |     2 |   1.5 MB | **Git** — `writeup/`                                                                                        |

Total: ~1,150 archive files in Git compressed; ~1,160 in Dropbox compressed.

## How these outputs fed downstream analyses

Four downstream views key off the per-OG MSA / tree / BUSTED bundles.

**BUSTED-as-BigBed UCSC selection track.** The orthogroup_id → PvP01 gene_id map from the ortholog table (see [ORTHOLOGY.md → Outputs](ORTHOLOGY.md#outputs)) is the join key. Each BUSTED JSON's `"test results"."p-value"` becomes a score; BH-FDR q-value across all OGs becomes the colour bin. The result is `ucsc_hub/GCA_900093555.2/selection_strict.bb` and `selection_relaxed.bb` (BigBed12+5, color by q-bin). PvP01-only for v1; replicating to the other 7 anchors deferred to v2.

**HyPhy report.** The 18-page `writeup/hyphy_report.pdf` ranks orthogroups by BUSTED q-value with family annotations, includes the Pvs230 BUSTED-MH appendix, and provides the MalariaGEN-frequency overlay for the top 200 hits. Driver: `writeup/build_hyphy_report.py` + `build_pvs230_appendix.py`.

**Family-stratified Manhattan plots.** Per-chromosome BUSTED -log10(q) plots colored by Phase G family label, separated by gene family (PIR, PHIST, MSP, DBP, EBA, RBP, …). Plots at `writeup/hyphy_plots/manhattan_by_family_*.png`. Highlights that variant-antigen genes carry the strongest positive-selection signal, as expected.

**Variant-antigen BUSTED bundle.** The `core_relaxed_hyphy_va/` subset re-runs BUSTED on just the 245 variant-antigen orthogroups (PIR, PHIST, MSP, DBP, EBA, RBP, AMA, SERA, TRAg, STP1, RESA). Same tool, narrower input — separating these from the housekeeping bulk makes per-family selection inferences cleaner. Driver: `scripts/overnight/11_variant_antigen_hyphy.sh`.

## Re-running on a different species

Prereqs: an ortholog table from the [ORTHOLOGY.md → Re-running on a different species](ORTHOLOGY.md#re-running-on-a-different-species) recipe, and the upstream annotations + soft-masked FASTAs it consumed.

Walk-through:

1. Edit `pipeline/species.conf` to set `MIN_INTACT_STRICT` and `MIN_INTACT_RELAXED` (scaled for the panel size — Pv4 used 7/5 for 8 strains; Pk v1 uses 6/4 for 7 strains).
2. `bash pipeline/07_msa.sh` — Stage F MSAs (strict + relaxed in parallel). Wall: ~5 hours total.
3. `bash pipeline/08_trees.sh` — Stage G IQ-TREE per orthogroup, both sets. Wall: ~7 hours.
4. `bash pipeline/09_hyphy.sh` — Stage H bulk BUSTED, both sets. Wall: ~5 hours.

Three scale-dependent parameters worth tuning:

- For very large orthogroup counts (> 10,000 OGs in the relaxed set), substitute `iqtree3 -m GTR+I+G` for `-m MFP` to skip ModelFinder — ~3× wall speedup, with the trade-off that you lose per-OG model selection.
- For tight intra-species data, the BUSTED test statistic has low power. Use the `--multiple-hits Double` flag for top hits before drawing conclusions — see the Pvs230 example above.
- For deeply paralogous gene families (variant antigens), separate them into their own subdirectory and re-run BUSTED with `--branches Internal` or a labelled tree. Mixed-orthology blocks confuse BUSTED's positive-selection inference.

The *P. knowlesi* scaffold at `Pk/v1/pipeline/07_msa.sh`, `08_trees.sh`, `09_hyphy.sh` is a parameterized copy for that species.

### Top-hit interpretation checklist

Before reporting any BUSTED top hit, run through:

1. **Codon MSA quality**. Open the MSA in a viewer; reject if > 30 % gaps in any column, or if any strain is a translation pseudogene fragment that pal2nal kept in by accident.
2. **Tree quality**. Inspect the per-OG tree topology; reject if monophyly violates known species relationships (suggests a hidden paralog).
3. **MNS check**. Re-run BUSTED with `--multiple-hits Double`. If p shifts by > 0.3, the original signal is suspect — flag it.
4. **MalariaGEN frequency overlay**. For each codon flagged by BUSTED, look up the allele frequency in MalariaGEN's 1,895-sample cohort (see [MALARIAGEN_VCF_PROJECTION.md](MALARIAGEN_VCF_PROJECTION.md)). Singleton alleles can drive BUSTED p-values; their selection inference is weak.
5. **Family stratification**. If the gene is in a variant-antigen family (PIR, PHIST, MSP, etc.), check the family-stratified Manhattan plot first — most family members will show positive selection by construction.

## Galaxy tool-wrap list

| Tool                                                     | Galaxy state                  | Note                                                                                                                                                                                  |
| -------------------------------------------------------- | ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **gffread** (CDS extraction)                             | ✅ exists                      | IUC; covered also in [ORTHOLOGY.md → Galaxy tool-wrap list](ORTHOLOGY.md#galaxy-tool-wrap-list)                                                                                       |
| **MAFFT** (LINSI mode)                                   | ✅ exists                      | IUC; verify the `--localpair --maxiterate 1000` flags are exposed                                                                                                                     |
| **pal2nal**                                              | ✅ exists                      | IUC; the wrapper handles `-output fasta`                                                                                                                                              |
| **trimAl**                                               | ✅ exists                      | IUC; verify `-automated1` is selectable                                                                                                                                               |
| **`build_msa.py`** (custom orchestrator with stop-strip) | ⚠ needs wrap                  | Custom Python that chains gffread + MAFFT + pal2nal + the internal-stop fix. Wrap as one tool that takes (ortholog table row, per-strain GFF collection, per-strain FASTA collection) |
| **IQ-TREE3**                                             | ⚠ partial                     | IQ-TREE2 wrapper exists; needs version bump to v3 (binary `iqtree3`). Most flags are unchanged                                                                                        |
| **HyPhy BUSTED**                                         | ✅ exists                      | IUC and datamonkey wrappers. Verify `--srv No`, `--branches All`, and `--multiple-hits Double` are all exposed                                                                        |
| **`build_selection_bigbed.py`** (BUSTED-to-BigBed12+5)   | ⚠ needs wrap                  | Custom Python that joins BUSTED JSONs to per-OG BED12 via the ortholog table, computes BH-FDR, emits BigBed12+5 with q-bin colour. Used by the UCSC publishing step                   |
| **HyPhy-Vision** (per-OG result browser)                 | ✅ external (vision.hyphy.org) | Deep-link via `?json=<datacache URL>` once HyPhy-Vision CORS is confirmed (tracking: `veg/hyphy-vision#892`)                                                                          |
| **MSAViewer.js** (in-browser MSA view)                   | not Galaxy — frontend         | BRC-analytics UI component; reads codon/protein MSAs from datacache                                                                                                                   |
| **phylotree.js** (in-browser tree view)                  | not Galaxy — frontend         | BRC-analytics UI component; reads `.treefile` from datacache                                                                                                                          |

The natural workflow shape is one `per-orthogroup-msa-tree-hyphy` workflow that iterates over a collection of ortholog table rows. Inputs: ortholog table + per-strain GFF collection + per-strain FASTA collection. Outputs: collections of `.codon.aln.fa`, `.protein.aln.fa`, `.treefile`, `busted.json`, plus a summary TSV.

The downstream `BUSTED-to-BigBed` + `HyPhy-Vision deep-link` views are wrappable as a separate `publish-selection-results` workflow that takes the per-OG BUSTED JSONs and produces the UCSC track + the Vision deep-link map. This separation keeps the heavy compute (MSAs + trees + BUSTED) decoupled from the publishing step.

Category in `workflow_categories.yml`: `SELECTION_ANALYSIS` (new — see brc-analytics issue #1279 PR plan).
