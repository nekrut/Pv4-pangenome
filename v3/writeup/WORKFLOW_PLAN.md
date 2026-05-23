# Comparative-genomics workflow plan — local → Galaxy

A reproducible, species-agnostic pipeline that takes N assemblies of a single species through pangenome construction, pairwise alignments, orthology, variation projection, and per-gene evolutionary analysis. Built from lessons learned on the 8-strain *Plasmodium vivax* analysis; designed to rerun on any eukaryotic pathogen the BRCs cover.

## Inputs (parametric per species)

| Input | Per-species value | Optional? |
|---|---|---|
| Assemblies (FASTA) | N ≥ 2, ideally 4–10 | required |
| GFF3 annotations | as many as available; rest projected via Liftoff | required for ≥1 |
| Reference proteomes | matching annotated assemblies | required for ≥1 |
| Cohort VCFs (population) | e.g. MalariaGEN Pv4 | optional (gates projection phase) |
| Priority gene list | curated locus tags + categories | optional (gates targeted HyPhy bundle) |
| Outgroup assembly | sister species | optional (enables rooted selection tests) |

All other inputs (chromosome lists, masking BEDs, k-mer settings) are derived.

## Phase A — Input preparation

**Tasks:** harmonize FASTA contig naming (PanSN format `sample#hap#ctg#frag`), mask low-complexity, build size/index files.

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| seqkit | header rewrite, FASTA stats | ✓ |
| samtools (faidx) | index | ✓ |
| sdust | short low-complexity | ✓ (wrapped) |
| **longdust** | long low-complexity windows | **needs wrapper** |
| bedtools (merge, maskfasta) | union mask + hard-mask | ✓ |

**Output:** softmasked + hardmasked FASTA per assembly; per-genome contig sizes; PanSN-named FASTA bundle.

## Phase B — Pangenome construction

**Tasks:** build a sequence-based PGGB graph from all N assemblies.

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| **PGGB** (orchestrator) | wfmash → seqwish → smoothxg | **needs wrapper** (or invoke as Docker step) |
| wfmash | all-pairs alignment | ✓ (wrapped) |
| seqwish | PAF → graph induction | ✓ |
| smoothxg | window-based normalization | ✓ |
| odgi | graph manipulation, stats, viz | ✓ (partially wrapped — need `odgi pav`, `odgi position`, `odgi viz`) |
| vg | deconstruct, paths | ✓ (basic), need `vg deconstruct` exposed |

**Key parameters to expose in Galaxy:** `-s` segment length (5000 for *P. vivax*), `-p` identity (90%), `-n` haplotypes (= N), block selection. Document defaults that work for compact (<50 Mb), AT-rich, polyploid-aware single-species sets.

**Output:** `.og`, `.gfa`, `.smooth.fix.{og,gfa}`; multi-strain MAF blocks via custom Python parsing of GFA shared-node runs (the `odgi untangle` SDSL bug we hit forced this — file upstream).

## Phase C — Pairwise alignments (chains for projection)

Three alignment paths produce three coordinate-projection methods. From this session: **KegAlign (A2) is the clear winner** for drug-resistance gene coverage; keep all three for QC triangulation.

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| wfmash → paftools.js → chain | A1 chain | wfmash ✓, paftools.js ✓ |
| **KegAlign** (GPU lastZ) | A2 chain via axt → axtChain → chainSort → chainPreNet → chainNet → netChainSubset → chainStitchId | **partially wrapped — needs Galaxy GPU profile + axtChain pipeline wrapped as a single tool** |
| odgi position + custom join | B graph-native projection | needs wrapper for the per-site position → VCF rewrite step |
| UCSC kentutils (axtToMaf, axtChain, chainSort, chainNet, netChainSubset, chainStitchId) | chain post-processing | ✓ (already in Galaxy under `ucsc-` prefix; consolidate into a chain-pipeline subworkflow) |
| paftools.js (minimap2 distribution) | PAF → chain | ✓ |

**Output per pair:** `.cleaned.chain` + `.rbest.chain` + `.axt` + sizes files.

**Multi-way alignments (multiz):**

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| axtToMaf | per-pair AXT → MAF | ✓ |
| **multiz** | iterative multi-way MAF builder | **needs wrapper** (biocontainer exists at `quay.io/biocontainers/multiz:11.2`) |

For each candidate hinge (reference strain), iteratively multiz the N−1 pairwise MAFs → one multi-way MAF per hinge. 8-hinge sweep on *P. vivax* produced 128k–302k multi-way blocks per hinge.

## Phase D — Orthology

**Multi-source consensus** (the strict-7 cut we used was too aggressive; default to **min-intact = N − 2** to retain variant antigens).

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| Liftoff | annotation projection from reference | ✓ (recently wrapped) |
| **TOGA2 / CESAR2** | CDS-aware annotation projection | **needs wrapper** (Docker image at `avianalter/toga2:latest`, 6.8 GB) |
| OrthoFinder3 | orthogroup discovery | ✓ |
| MMseqs2 (easy-rbh) | fast pairwise RBH for triangulation | ✓ |
| gffread | extract CDS + protein from GFF+FASTA | ✓ |
| custom Phase-E union-find merger | reciprocal anchor reconciliation + position-aliasing dedup | **needs wrapper** (Python script with classification.tsv input) |
| custom Phase-G family tagger | description-based PIR/PHIST/MSP/etc. annotation | **needs wrapper** (regex-driven) |

**Pitfalls we hit (document in wrapper):**
- Liftoff IDs (`PVP01_xxx`) collide across strains → use position-based 50% reciprocal overlap to merge cross-anchor aliases
- NCBI vs PlasmoDB gene-ID conventions need a strain-aware mapping table
- Native-BED seeding is required for non-PvP01 strains, or aliases miss the union-find
- `proteomes` need cleaning before OrthoFinder: strip headers to single ID, deduplicate (PAM had 50 ENA→mRNA collisions), replace non-canonical AAs (`.`, `*`, `J`, `U`, `O` → `X`)

## Phase E — Variant projection (population VCF → all references)

Compare three paths; declare the consensus set as high-confidence.

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| CrossMap (HKUST) | chain-based VCF liftover | ✓ |
| bcftools (sort, concat, view, index, liftover plugin) | post-processing | ✓ (sort + concat); `+liftover` plugin needs the SCORE bundle (Genovese) — **wrapper missing** |
| odgi position | per-site coordinate translation through graph | needs wrapper for the population-VCF rewrite step (custom Python join: position TSV + source VCF → target VCF) |

**Critical bug from session — document in wrapper:**
- CrossMap output is unsorted → ALWAYS pipe through `bcftools sort -Oz -T <scratch>` before `tabix`
- VCF headers must include target-genome `##contig=` lines or `bcftools sort` rejects records
- Strand-flipped chain blocks need REF/ALT reverse-complement during graph-position rewrite

## Phase F — Per-gene evolutionary analyses

### F.1 Codon-aware MSA (per CORE-1:1 orthogroup)

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| pyfaidx + Python | per-strain CDS extraction from merged GFF | ✓ |
| MAFFT | protein alignment | ✓ |
| MACSE | codon-aware MSA (frameshift-aware) | needs wrapper (biocontainer exists) |
| pal2nal | protein → codon backtranslation | ✓ |
| trimAl `-automated1` | column-level trimming | ✓ |

**Default: min_intact = N − 2.** Strict-(N) excludes most variant antigens.

### F.2 Tree inference

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| IQ-TREE3 (`-m MFP`) | ML tree | ✓ |

**Document fallback:** drop `-B 1000` for alignments with <4 unique sequences (IQ-TREE rejects bootstrap otherwise). ~25% of CORE-1:1 *P. vivax* genes hit this.

### F.3 Selection analysis

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| HyPhy (BUSTED, aBSREL, MEME, FEL, RELAX, SLAC, Contrast-FEL) | per-gene + per-site selection | ✓ (BUSTED/aBSREL/MEME/FEL); **RELAX + Contrast-FEL need wrappers** |
| **BUSTED --multiple-hits Double** (MH variant) | MNS-robust whole-gene test | **needs wrapper** — important! MNS events generalize across species; we showed Pvs230's "two parallel selection events" was a single MNS, p shifted 0.056 → 0.500 |
| GARD | recombination pre-screen | ✓ |
| 3SEQ | breakpoint detection | needs wrapper |
| ASTRAL-IV | species tree from gene trees | needs wrapper |
| Twisst | topology weighting | needs wrapper |
| HyPhy SLAC, Contrast-FEL | additional site tests | partial — Contrast-FEL needs wrapper |

### F.4 Within-species selection from cohort VCF

| Tool | Role | Galaxy-ready? |
|---|---|:-:|
| Custom Python (codon table + MalariaGEN VCF) | per-codon N/S, domain stratification, MNS detection via phasing | **needs wrapper** — generalizable as a "per-gene codon-level cohort QC" |
| McDonald-Kreitman test scripts | polymorphism × divergence | needs wrapper |

## Phase G — Reporting

Generate per-species writeup containing: assembly inventory, masking stats, orthology counts (per-strain + cross-strain), 3-way projection intersection, HyPhy results, family heatmaps, microsynteny plots.

Local Python tooling we built: subtelomeric microsynteny renderer, family heatmap, FEL/MEME per-site lollipops, Manhattan p-value plot, p-value Q-Q, family-enrichment bars. All matplotlib — straightforward to wrap as report-generation tools.

## Galaxy migration — tools that need wrapping (priority-ordered)

These are the gaps preventing a full Galaxy rerun on a new species:

| # | Tool | Why critical |
|---|---|---|
| 1 | **PGGB orchestrator** | central pangenome step; currently a Docker shell script |
| 2 | **multiz** | multi-way MAF builder; no Galaxy wrapper exists |
| 3 | **TOGA2/CESAR2** | CDS-aware projection where Liftoff fails |
| 4 | **longdust** | masking step; sdust+longdust union is what PGGB needs |
| 5 | **odgi pav** + **odgi position** | graph-native variant projection; required for Path B |
| 6 | **odgi viz** (subway / 1D) | reporting figure |
| 7 | **HyPhy BUSTED-MH** (multi-hit) | MNS-robust whole-gene selection; mandatory caveat tool |
| 8 | **HyPhy RELAX + Contrast-FEL** | branch-set selection tests |
| 9 | **MACSE** | codon-aware MSA where MAFFT+pal2nal struggles on frameshift |
| 10 | **bcftools `+liftover` plugin (SCORE bundle)** | VCF liftover via chain |
| 11 | **Phase E reconciler** (Python) | union-find orthology merger with position aliasing |
| 12 | **Path B population-VCF rewriter** | join odgi position output + MalariaGEN VCF → target-coord cohort VCF |
| 13 | **GARD breakpoint scanner** + **ASTRAL-IV** | recombination pre-screen + coalescent species tree |

Other steps (CrossMap, bcftools sort/concat, MAFFT, pal2nal, trimAl, IQ-TREE, OrthoFinder3, gffread, Liftoff, BUSTED/aBSREL/MEME/FEL, samtools, bedtools, seqkit) are already in Galaxy and just need to be assembled into a parameterized workflow.

## Local-vs-Galaxy execution plan

### On this workstation (first pass)

Single driver script: `scripts/run_species.sh <species_name> <input_dir>` that chains the phases. Output to `work/<species_name>/{01_..09_}`. Use Docker per-tool with bind-mounted scratch for sort tempfiles. Parallelize at the per-gene level (Phase F.1, F.2, F.3) with `xargs -P 8` and at the per-strain level for VCF projection.

**Lessons baked into the driver:**
- Pre-build per-genome contig-header snippets so bcftools sort accepts CrossMap output
- Always atomic .vcf → .vcf.gz (only `rm -f .vcf` after `.vcf.gz` validated)
- IQ-TREE retry without `-B 1000` on alignments with <4 unique sequences
- HyPhy BUSTED `--multiple-hits Double` for any whole-gene call before publishing
- Liftoff feature_types must include `protein_coding_gene ncRNA_gene pseudogene` (default `gene` matches nothing in PlasmoDB)
- OrthoFinder proteomes must be deduplicated + non-canonical AAs replaced

### Galaxy workflow (second pass)

One **invocation workflow** per species that takes a "species manifest" YAML and produces all phase outputs. Subworkflows per phase. Per-tool wrappers from the priority list above unlock the unwrapped steps.

The manifest is the only thing that changes per species — chromosome lists, masking parameters, alignment thresholds all derive from a sensible default + manifest overrides.

## Estimated effort

- **Local driver script** (collating what exists): ~1 week
- **Wrapping priority tools 1-13**: ~3-4 weeks if done in parallel with focused IUC/contributors (multiz, PGGB, odgi pav/position/viz, HyPhy variants are the largest)
- **Galaxy workflow assembly**: ~1 week once wrappers exist
- **Per-species reruns**: ~24-48 hours wall on this hardware tier; faster with GPU node for KegAlign A2 chains

## Outstanding science gaps (carry forward)

1. **Variant antigen family analysis** (telomeric PIR/PHIST/Pv-fam) — needs per-family clustering pipeline, not the 8-way single-copy ortholog filter
2. **GARD recombination pre-screen** on bulk BUSTED hits before publishing
3. **Outgroup-rooted selection tests** — needs a *P. cynomolgi* or *P. knowlesi* ortholog injection
4. **McDonald-Kreitman population × divergence** test using projected cohort VCFs
5. **Per-region geographic stratification** of within-species AF calls
