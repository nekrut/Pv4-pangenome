# LOCAL.md, in plain language

A readable walkthrough of what the pangenome + selection-scan pipeline actually computes. `LOCAL.md` is the machine recipe — every command, guard, and pinned container. This document covers the same analyses at the level of "what is done and why," without the bash.

The pipeline is species-agnostic: it takes N haploid assemblies (5–15) with their annotations and a cohort VCF, and everything species-specific lives in one config file. The worked example throughout is **Pv4** — the eight-strain *Plasmodium vivax* panel the pipeline was built on — shown in parentheses where a concrete number helps.

## What it does

Given N assemblies with their annotations, the pipeline: aligns all pairs and builds cleaned UCSC chains; projects each anchor strain's annotation onto every other strain; builds a PGGB pangenome graph; reconciles three independent ortholog signals into one orthogroup table; builds codon-aware MSAs per orthogroup; infers per-gene trees and runs BUSTED selection tests; folds the pairwise alignments into per-strain multiz MAFs; re-maps a cohort VCF from the reference's coordinates onto every other assembly; and packages everything as a UCSC track hub.

## Inputs

Set once in `species.conf`:

- **Assemblies** — one FASTA per strain (Pv4: 8).
- **Annotations** — one GFF3 per strain.
- **Cohort VCF** — population variant calls, one file per chromosome, all in the reference's coordinates (Pv4: the 1,895-sample MalariaGEN cohort, in PvP01 coordinates).
- **Chromosome-name map** — needed only when the annotation's sequence names differ from the assembly's.

The config also names the strain set, the reference, and the **anchor subset** — the strains whose own annotation is trusted enough to use as a projection source. Anchors lend their gene sets to every other strain (Phase C.4); non-anchors contribute sequence but not their lower-confidence annotation. The reference is always one of the anchors (Pv4: 5 PlasmoDB-curated anchors, 3 NCBI-only non-anchors).

Swap the config and the same pipeline runs on another species — see the last section.

### Building the chromosome-name map

Whether you need this depends on where the annotation came from. An NCBI annotation usually already names sequences by their GenBank accession, so no map is needed. A community database (Pv4: PlasmoDB; also Ensembl, the VEuPathDB family, an in-house annotation) typically names sequences its own way (`PvP01_14_v1`) while the assembly FASTA uses the accession (`LT635625.2`) — then you need a 2-column map to reconcile them, or the lift tools match nothing silently.

The map only has to cover the reference assembly — both the reference annotation and the cohort VCF live in reference coordinates, while the target assemblies are already accession-named. It's a 2-column TSV, annotation name `<TAB>` accession.

The canonical source is the NCBI assembly report (`*_assembly_report.txt`), which carries both names per sequence — the submitter's `Sequence-Name` and the `GenBank-Accn`:

```bash
awk -F'\t' '!/^#/ {print $1 "\t" $5}' <accession>_assembly_report.txt > ref_to_genbank.tsv
```

When there's no report, match the two FASTAs by sequence length — the same assembly under two naming schemes is byte-identical, so length is a unique key (this is what `fix_gff_chroms.sh` falls back to). Watch for version drift: some databases bump a sequence-version suffix between releases (PlasmoDB: `PvP01_14_v1` → `_v2`) while the accession stays fixed, so build the map from the same release you're projecting.

## Essential outputs

27 in total. The ones that carry the science:

- **`ortholog_table.tsv`** — per orthogroup: which strains, copy number, where the three signals agreed.
- **Codon + protein MSAs** per orthogroup, strict and relaxed sets.
- **BUSTED result** per gene.
- **Cohort VCF re-mapped** onto each non-reference assembly.
- **UCSC track hub** tying it together.

## Dependency shape

```
 inputs: assemblies + annotations + cohort VCF
    │
    ▼
 A inventory ─ B mask ─┬─ C.1-3 align + chains
                       │      └─ C.4 annotation projection
                       ├─ D PGGB graph
                       └─ (chains + graph + projections feed:)
                          E consensus orthology
                             └─ F codon/protein MSAs
                                  ├─ G IQ-TREE
                                  │     └─ H HyPhy BUSTED
                                  └─ (separately)
                          I multiz MAFs           (needs C.1-3 axt)
                          J cohort VCF projection (needs C.1-3 chains)
                          K UCSC hub              (needs C.4, E, G, H, I, J)
    │
    ▼
 outputs: ortholog table · MSAs · BUSTED · re-mapped VCFs · UCSC hub
```

Independent phases run concurrently; the orchestrator blocks only where one phase consumes another's output.

## Phase A — Inventory

```
 assemblies                         proteomes
     │                                  │
     ▼                                  ▼
 sourmash sketch                     BUSCO
 (k=31, scaled=1000)              (clade odb10)
     │                                  │
     ▼                                  ▼
 {S}.sig.gz ──> sourmash compare    completeness %
     │               │                  │
     ▼               ▼                  ▼
 [BRC catalog]   compare.csv        QC: trust the
                 (N×N similarity)   annotation?
                     │
                     ▼
              [Phase I fold order]
```

`sourmash sketch dna -p k=31,scaled=1000` per strain, then `sourmash compare` gives the N×N similarity matrix. The same tool runs in the BRC deployment, so the per-strain `.sig.gz` sketches computed here are the exact artifacts the catalog ingests — compute once, reuse there. The matrix orders the multiz fold in Phase I (closest = highest similarity). One caveat downstream: `sourmash compare` reports similarity (1.0 = identical), not distance, so the fold-order step sorts descending.

`BUSCO` runs on each proteome (not the genome — faster, and these annotations are gene-dense) against the lineage closest to the panel (Pv4: `plasmodium_odb10`). The completeness scores are the sanity check on annotation quality before anything downstream trusts it.

## Phase B — Soft-masking

```
 assembly {S}.fa
     │
     ├──> longdust ──> {S}.longdust.bed ┐
     │                                  ├──> cat | sort | merge ──> {S}.union.bed
     └──> sdust    ──> {S}.sdust.bed    ┘                              │
                                                                       ▼
 assembly + union.bed ──> bedtools maskfasta -soft ──> softmasked/{S}.fa
                                                            ├──> .fa.fai
                                                            └──> .sizes ──> [Phases C, I]
```

Per strain, take the union of `longdust` and `sdust` intervals and soft-mask with `bedtools maskfasta -soft`. Soft, not hard — lowercase, not N, because the k-mer tools downstream need the bases present. The step also writes the 2-column `.sizes` file that `chainNet` / `chainPreNet` require later; cheaper to bake it now than regenerate it.

## Phase C — Pairwise alignment and chains

### C.1–3, alignment → chains

```
 softmasked A.fa, B.fa
     │
     ▼
 KegAlign (GPU)  ──or──  lastz (CPU)
     │
     ▼
 A__vs__B.axt ─────────────────────────────────> [Phase I: pairwise MAFs]
     │
     ▼
 axtChain -linearGap=loose
   → chainSort → chainPreNet → chainNet → netChainSubset → chainStitchId
     │
     ▼
 {A}.{B}.cleaned.chain ──> [C.4 TOGA · Phase J VCF · Phase K bigChain]
     │
     ▼  swap → net → subset → swap back
 {A}.{B}.rbest.chain   ──> [Phase E: orthology edges]
```

All N(N−1)/2 unordered pairs (Pv4: 28) through KegAlign on GPU (`--strand both --hsp_threshold 5000 --gapped_threshold 6000 --inner 2000 --ydrop 15000`, AXT out), or `lastz` on CPU at ~10× the cost. `hsp_threshold 5000` is deliberate — lower values bloat chain output ~10× with spurious short HSPs.

Then the canonical UCSC chain build, both directions per pair: `axtChain -linearGap=loose → chainSort → chainPreNet → chainNet → netChainSubset → chainStitchId`. `loose` suits panels that diverge several percent in places; use `medium` for tighter intra-species sets. The pipeline is not optional — dropping `chainPreNet` or `netChainSubset` yields ~10× more chains and breaks CrossMap in Phase J, which then picks tiny overlapping chains first.

Also build reciprocal-best chains per pair (swap → net → subset → swap back). These keep only mutual best alignments — ~30–50% fewer, more conservative. The pipeline uses cleaned chains for variant projection (rbest drops ~40% of liftable variants) and rbest for orthology evidence (where the conservatism is the point).

### C.4, annotation projection

```
 anchor A annotation ──> Liftoff (A → Q) ──> triage (phase_c2_triage.py)
                                                  │
                ┌─────────────────────────────────┴──────┐
                ▼                                        ▼
              clean                          low_cov / partial / frameshift /
                │                            short_exon / truncated
                │                                        │
                │                                        ▼
                │                                 needs_cesar2.bed
                │                                        │
                │                                        ▼
                │                          TOGA2 / CESAR2 (cleaned chain)
                │                          exon realign + classify
                ▼                                        ▼
        liftoff_clean.gff3                       annotation.gff3
                │                                        │
                └────────────────┬───────────────────────┘
                                 ▼
                    merge (phase_c4_merge.py)        (split, extra_copy:
                                 │                    flagged in triage.tsv)
                                 ▼
        {Q}.annotation.gff3 + classification.tsv
        tags: LIFTOFF_CLEAN / CESAR_RESCUE / ...  ──> [Phase E]
```

The native annotations are heterogeneous — different groups, different methods, not comparable gene-for-gene. So C.4 discards them as a *source* and re-derives each anchor's gene set in every other genome. The result is one self-consistent gene set per anchor, re-expressed across the panel. It runs over every (anchor `A`, target `Q`) pair (Pv4: 5 anchors × 7 others = 35 projections), each independent; the reference's own curated GFF is the input for its row. Four sub-steps per pair.

**Liftoff — the fast path.** Alignment-based feature transfer: align each reference gene's sequence to `Q` with minimap2, then place it where coverage × identity is maximized *while preserving exon/intron structure* (the gene lifts as a unit, and competing mappings are resolved against each other). `-copies -sc 0.95` adds paralog detection — extra copies in `Q` scoring ≥0.95, which matters for the expanded multigene families. Restricted to `protein_coding_gene / ncRNA_gene / pseudogene`; fed the chromosome-renamed GFF so seqids match. Liftoff is right when sequence is conserved and structure intact, and degrades on divergent, frameshifted, partial, or split genes.

**Triage — route by how well Liftoff did.** `phase_c2_triage.py` scores each lifted gene (coverage, identity, frameshift, internal exon sizes, copy number, contig spread) and bins it into eight classes. **clean** (≥0.95 identity, full-length, in-frame) goes straight to `liftoff_clean.gff3`. **low_cov / partial / frameshift / short_exon / truncated** go to `needs_cesar2.bed` — the rescue worklist. **split** (multiple contigs) is flagged only, not rescued. **extra_copy** is recorded as a paralog.

**TOGA2 / CESAR2 — rescue the hard tail.** A different approach: chain-based ortholog inference, not sequence liftover. It uses the cleaned chain from C.1–3 to locate each reference gene's syntenic region in `Q`, then CESAR2 realigns the coding exons with explicit splice-site and reading-frame awareness — reconstructing CDS structure across the divergence Liftoff choked on, and classifying each projection (intact / partial / missing). `--filter_bed needs_cesar2.bed` restricts it to the triaged genes (~2 h/anchor instead of ~6 h re-projecting everything).

**Merge.** `phase_c4_merge.py` takes the Liftoff call for clean genes, the CESAR call for rescued ones, keeps both for extra copies, and tags every gene with its provenance — `LIFTOFF_CLEAN / CESAR_RESCUE / CESAR_PARTIAL / MISSING / SPLIT / EXTRA_COPY`.

The division of labor is the point: Liftoff does coordinate transfer on intact genes (the cheap majority), TOGA/CESAR does ortholog reconstruction on broken or divergent ones (the hard tail), and triage routes between them. The provenance tag is the hand-off to Phase E, which trusts only `LIFTOFF_CLEAN` and `CESAR_RESCUE` genes as orthology edges. Output: a projected, traceable annotation in every strain, per anchor.

## Phase D — PGGB graph

```
 softmasked {S}.fa  (all N)
     │
     ▼
 PanSN rename (SAMPLE#1#CONTIG)
     │
     ▼
 concat + bgzip ──> all_pansn.fa.gz
     │
     ▼
 pggb -n N -s 5000 -p 90 -k 23
   (wfmash → seqwish → smoothxg → gfaffix → odgi)
     │
     ▼
 pv.{og,gfa} ──> [Phase E: graph co-membership edges]
```

PanSN-rename (`SAMPLE#1#CONTIG`), concat, bgzip, then `pggb -n N -s 5000 -p 90 -k 23`. `-n` must equal the strain count; `-p 90 -s 5000` suits a compact intra-species panel (drop to `-p 80` for more divergent sets, raise `-s` for larger genomes). The canonical product is `*.smooth.fix.{og,gfa}`, symlinked to `inputs/pggb/pv.{og,gfa}`. The graph is the third orthology signal and a standalone artifact. ~5 h on 32 cores for Pv4.

The input is the soft-masked FASTAs from Phase B — but PGGB neither needs nor uses masking. wfmash is case-insensitive at the minimizer level, so the soft mask is inert here: PGGB effectively builds from full sequence, and repeat tangling is controlled by `-p` / `-s` / segment length, not by masking the input. Hard masking would be actively wrong — runs of N fragment the graph (wfmash can't seed across them, seqwish breaks paths), dropping exactly the low-complexity regions and corrupting the surrounding topology. So the same soft-masked set serves both phases: it shapes the Phase C pairwise alignments (lastz/KegAlign do skip lowercase seeds) and is harmlessly ignored here.

## Phase E — Consensus orthology

```
 C.4 projection ──────> projection edges  (tag ∈ {LIFTOFF_CLEAN, CESAR_RESCUE})
 C.1-3 rbest chains ──> rbest edges        (≥0.90 overlap both ways)
 D graph paths ───────> graph edges        (≥0.90 CDS overlap)
        │
        ▼   all edges over nodes = (strain, gene)
   union-find  ──>  connected components = orthogroups
        │
        ▼
   ortholog_table.tsv  (one row per orthogroup)  ──> [Phase F]
```

No single method calls orthologs reliably across a panel with expanded multigene families, so E combines three independent signals. The model: every (strain, gene) is a **node**; each signal contributes **edges** between nodes it believes are orthologous; union-find collapses the connected components into orthogroups. It's not a majority vote — any one high-confidence edge merges two genes — so the quality gate lives on each edge, not on the count.

The three edge sets, each with a different failure mode:

**Projection edges** — from C.4. An anchor gene links to the target gene it projected onto, kept only where the merge tagged it `LIFTOFF_CLEAN` or `CESAR_RESCUE`. Directional (anchor → target), so two targets connect transitively through their shared anchor. Strong where the annotation transfer was clean; blind to genes that exist in no anchor.

**Reciprocal-best edges** — `phase_e_rbest_overlap.py` over the rbest chains from C.1–3, keeping pairs whose CDS overlap is ≥0.90 *both ways*. Sequence-level and symmetric, covering all pairs rather than only anchor-involving ones, and independent of annotation quality. Conservative by construction (rbest already dropped the ambiguous alignments); blind in regions that align poorly.

**Graph co-membership edges** — `odgi paths --haplotypes` exports each strain's path through the PGGB graph, and `phase_e_graph_edges.py` links genes whose CDS traverse the same graph path with ≥0.90 overlap. Catches orthology through the pangenome structure, independent of pairwise alignment; blind where the graph tangles (collapsed repeats, paralog families).

`phase_e_consensus.py` builds the node set from every annotation, adds all three edge sets, runs union-find, and emits one row per component to `ortholog_table.tsv` (14 columns: id, label, n_strains, max_copies, per-strain members, graph stats). Cell notation: `|` separates cross-anchor aliases, `,` separates paralog copies, `-` marks absent. A few thousand orthogroups for a panel this size. The three sources are complementary, not redundant — each covers the others' blind spots, and a pair linked by all three is the most robust kind of call.

## Phase F — Codon and protein MSAs

```
 ortholog_table.tsv  +  softmasked FASTAs  +  merged GFFs
     │
     ▼
 gffread -x : pull CDS per member ──> strip internal stops (→ NNN)
     │
     ▼
 translate ──> MAFFT-LINSI (protein MSA)
     │
     ▼
 pal2nal back-translate ──> {og}.codon.aln.fa (+ .protein.aln.fa)
     │                              │
     │                              ▼  trimal -automated1
 validate: codon len = 3×protein   {og}.codon.cleaned.fa
     │
     ▼
 strict set (~1,584)  +  relaxed set (~4,215) ──> [Phases G, H]
```

Per orthogroup, per intactness threshold: pull CDS (`gffread -x`), translate, align proteins with MAFFT-LINSI (`--localpair --maxiterate 1000`), back-translate with `pal2nal`. LINSI is ~10× slower than default MAFFT but default MAFFT produces gappy garbage on divergent paralog families. `pal2nal` drops any strain with an internal stop, so `build_msa.py` strips internal stops to NNN first. Validation: codon length = 3 × protein length.

Two sets by `min_intact`: **strict** (Pv4: ~1,584 orthogroups) and **relaxed** (~4,215). `trimal -automated1` produces cleaned siblings. Orthogroups are named by the reference gene ID when 1:1 with the reference (Pv4: `PVP01_...`).

## Phase G — IQ-TREE

```
 {og}.codon.aln.fa
     │
     ▼
 count unique sequences
     │
   ┌─┴──────────────┐
   ▼                ▼
 ≥4 unique       <4 unique
 iqtree3 -m MFP  iqtree3 -m MFP
 -B 1000         (no bootstrap)
   │                │
   └───────┬────────┘
           ▼
     {gene}.treefile ──> [Phase H]
```

Per gene, `iqtree3 -m MFP -B 1000 -T 2`. Genes with <4 unique sequences drop `-B 1000` — bootstrap resampling hangs below four, and large near-identical gene families routinely fall under that. Only the `.treefile` is essential; the rest of each workdir is tar+gzipped (thousands of genes × ~5–10 MB otherwise).

## Phase H — HyPhy BUSTED

```
 {og}.codon.aln.fa  +  {gene}.treefile
     │
     ▼
 hyphy busted --srv No --branches All
     │
     ▼
 busted.json ──> validate (parse "test results")
     │
     ▼  top hits only
 re-run --multiple-hits Double
     │
     ▼
 selection call ──> [Phase K: selection track]
```

Per gene, `hyphy busted --srv No --branches All` over the Phase F alignment and Phase G tree. If the tree is missing, BUSTED writes an empty JSON — hence the parse-the-JSON validation. For top hits, re-run with `--multiple-hits Double` before trusting them: a single multi-nucleotide substitution can imitate diversifying selection (in Pv4, Pvs230 went p=0.056 → 0.500 under MH). BUSTED-MH should be the pre-publication default for top hits.

## Phase I — Multiz MAFs

```
 C.1-3 axt (H as source)               Phase A compare.csv
     │                                      │
     ▼                                      ▼
 axtToMaf -tPrefix=H. -qPrefix=Q.       fold order (closest first)
     │                                      │
     ▼                                      │
 H_vs_Q.maf  (one per other strain)         │
     │                                      │
     └─────────────────┬────────────────────┘
                       ▼
            multiz_progressive.py
            (merge closest, fold next, ... on H's scaffold)
                       │
                       ▼
             {H}.multiz.maf  ──> [Phase K: bigMaf]
             (one per hinge; Pv4: 8)
```

A multiz MAF is one strain's coordinate frame with the rest of the panel aligned into it, block by block. I builds N of them — one per strain as the **hinge** (the reference row) — so the alignment can be browsed from any strain's point of view. Two sub-steps per hinge.

**Pairwise MAFs.** Convert the C.1–3 pairwise alignments where the hinge `H` is the source into MAF with `axtToMaf -tPrefix=H. -qPrefix=Q.` (the prefixes give blocks PanSN-like names). One `H`-vs-`Q` MAF per other strain; these are kept as the intermediates that feed the fold.

**Progressive fold.** multiz only merges *pairwise* MAFs, so `multiz_progressive.py` grows the multi-way alignment one strain at a time: merge the two closest pairwise MAFs, then fold the next strain's pairwise MAF into that running alignment, and so on — every block threaded onto `H`'s coordinate scaffold. Order is closest-first by the Phase A sourmash similarities, because multiz threads each new strain transitively through what's already aligned: fold a distant strain early and its blocks fail to thread and get dropped. The wrong order loses 30–40% of blocks.

One MAF per strain (Pv4: 8), feeding the bigMaf tracks in Phase K. Slowest phase — ~24 h serial for Pv4, ~6 h with hinges run in parallel (each hinge is independent, but each holds a multi-GB MAF in flight, so parallelism is RAM-bound).

## Phase J — Cohort VCF projection

```
 cohort VCF (per chr, reference coords)
     │
     ▼
 bcftools annotate --rename-chrs ──> GenBank-named VCF
     │
     ▼
 CrossMap vcf  (cleaned chain ref→target, target FASTA)
     │
     ▼
 bcftools sort (fast scratch) ──> bcftools index ──> per-chr VCF
     │
     ▼
 bcftools concat -a ──> cohort_on_{TGT}.vcf.gz
     │
     ▼
 [one per non-ref target; Pv4: 7]
```

Rename the cohort VCF chromosomes to the assembly accessions (`bcftools annotate --rename-chrs`), then `CrossMap vcf` per non-reference target using the **cleaned** chain (rbest is too conservative here). CrossMap writes unsorted output — chain-spanning liftover interleaves chromosomes — so pipe through `bcftools sort` (on fast scratch, not tmpfs; a large cohort sort exhausts RAM otherwise) before indexing, or indexing fails. Concat per-chr with `-a` to allow the duplicates inversion liftover produces. Output: the full cohort in each non-reference coordinate system (Pv4: 7 targets, ≥100k variants surviving each).

## Phase K — UCSC track hub

```
 multiz MAFs    ─> mafToBigMaf → bedToBigBed ────────────────> bigMaf .bb (+ mafIndex)
 cleaned chains ─> chain_to_bigChain.py → bedToBigBed ───────> bigChain .bb + link .bb
 merged GFFs    ─> gff3ToGenePred → genePredToBed → bedToBigBed ─> BigBed12
 BUSTED jsons   ─> build_selection_bigbed.py ───────────────> BigBed12+5
 ortholog table ─> build_orthogroup_bigbed.py ──────────────> BigBed12
        │
        ▼   all .bb  +  hub.txt / genomes.txt / trackDb.txt
   ucsc_hub/{ACC}/
        │   constraints: composites single-type (bigMaf ≠ bigChain);
        │                genomes.txt needs real defaultPos + working twoBitPath
        ▼
   hubCheck ──> push
```

Convert each upstream output to the hub's binary formats and write the manifests (`hub.txt`, `genomes.txt`, per-assembly `trackDb.txt`).

- multiz MAFs → **bigMaf** (`mafToBigMaf` → `bedToBigBed -type=bed3+1`), plus a `mafIndex`.
- cleaned chains → **bigChain**: `chain_to_bigChain.py` emits the bigChain.bed (6+6) and bigLink.bed (4+1), each `bedToBigBed`'d. The track references both — `bigDataUrl` to the chain `.bb`, `linkDataUrl` to the link `.bb`. The browser will not load gzipped chain text.
- merged annotations → BigBed12 (`gff3ToGenePred → genePredToBed → bedToBigBed`).
- BUSTED results → BigBed12+5, scored by −log10(q), colored by q-bin.
- ortholog table → BigBed12, colored by n_strains.

Two hub constraints that cost real time: composites must be single-type, so bigMaf and bigChain live in separate tracks (not one composite); and `genomes.txt` needs a real `defaultPos` (any landmark locus — Pv4 uses *dhps*, a drug-resistance gene) and a working `twoBitPath` per assembly, or the hub fails to load. `hubCheck` validates before push. ~95 min, ~16 GB for Pv4.

## How it runs

`run_all.sh` runs the phases in dependency order, parallelizes the independent ones, and `wait`s at the joins. Every step is idempotent — `[[ -s $OUT ]] && continue` — so a killed run resumes cleanly and nothing recomputes. Every tool runs in a pinned container (Docker, or Singularity on clusters), so there are no local installs and no env drift between machines. Every step validates its own output, and `verify_essentials.sh` confirms all 27 essentials exist and pass before the run is called complete.

## Cost

32 cores, 1 GPU, 128 GB RAM, fast NVMe scratch (Pv4): ~25.5 h wall, ~1,410 CPU-hrs. The critical path is C.4 (TOGA dominates, ~12 h) and I (multiz, ~24 h serial); both shrink with parallelism. A 30-minute smoke test (3 strains, 1 chromosome, skipping multiz) catches broken inputs before a full run.

## Failure modes from v3

- **Liftoff lifted 0 genes** — annotation seqids didn't match the FASTA's accessions. Build the chromosome-name map first.
- **Codon MSAs short a few strains** — pseudogene CDS with internal stops break `pal2nal`. Strip internal stops to NNN first.
- **CrossMap output wouldn't index** — written unsorted. `bcftools sort` before index; sort on fast disk.
- **Lost an intermediate** — `rm -f` ran in the sort-failure branch. Only delete after the sorted output is verified.
- **IQ-TREE hung** — `<4` unique sequences with `-B 1000`. Count uniques; drop bootstrap below four.
- **False positive selection** — two "parallel" alleles were one MNS event. Re-test top hits with `--multiple-hits Double`.

## New species

Edit `species.conf`: species name, strain set, reference, anchors, cohort VCF path/glob, chromosome-rename map, BUSCO lineage, and the multigene-family regex (the lineage-specific families that confuse orthology — Pv4: `PIR|PHIST|Pv-fam|MSP|...`; *P. falciparum*: `var|rifin|stevor|...`). Tuning that may move with the panel: `hsp_threshold` (lower for AT-rich genomes — Pf uses 4000 vs Pv's 5000), the PGGB `-p`/`-s` pair for divergence and genome size, and `axtChain -linearGap`. Scripts, containers, and validation logic are unchanged.
