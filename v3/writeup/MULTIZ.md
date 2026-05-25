# Pairwise and multi-way alignments

The 8-strain pairwise alignments (KegAlign) and the 8 multi-way MAFs (multiz) are the second alignment layer in v3, sitting alongside the PGGB graph. They drive the canonical drug-resistance variant projection, every annotation projection step, the consensus orthology second source, and the UCSC multi-z hub track.

**Related documents.** [PANGENOME.md](PANGENOME.md) covers the parallel PGGB graph layer (built from the same soft-masked FASTAs); [ORTHOLOGY.md](ORTHOLOGY.md) covers the consensus ortholog table that consumes our cleaned chains for annotation projection and our rbest chains as one of the orthology evidence streams; [MALARIAGEN_VCF_PROJECTION.md](MALARIAGEN_VCF_PROJECTION.md) covers the population-cohort variant projection — the canonical (Path A2) workflow that lifts MalariaGEN samples through our cleaned chains.

## How we built them

### Pairwise alignments — KegAlign

KegAlign is a GPU-accelerated lastZ implementation. For each of the 28 unordered strain pairs we ran one KegAlign job and recovered the AXT output (one block per HSP):

```
docker run --gpus device=0 --rm -v $WORK:$WORK \
  quay.io/biocontainers/kegalign-gpu:1.0.0--hdfd78af_0 \
  kegalign \
    genomes/softmasked/{A}.fa \
    genomes/softmasked/{B}.fa \
    --strand both \
    --hsp_threshold 5000 \
    --gapped_threshold 6000 \
    --inner 2000 \
    --ydrop 15000 \
    --output_format axt \
  > projection/A2_kegalign/axt/{A}__vs__{B}.axt
```

Parameters chosen to match UCSC lastZ conventions:

- `--strand both` — emit both forward and reverse alignments
- `--hsp_threshold 5000` — initial HSP score gate (raised from default 3000 to keep spurious short hits out of the 28 GB downstream output)
- `--gapped_threshold 6000` — gapped-extension score gate
- `--inner 2000`, `--ydrop 15000` — standard UCSC tuning for pairwise mammalian-scale runs, works fine for ~25 Mb apicomplexan genomes

GPU wall time on an NVIDIA RTX A5000 (24 GB): ~3 min per pair, ~1.5 hours sequential for 28 pairs. CPU fallback via plain lastZ: ~30 min per pair.

Soft-masked FASTAs went in (longdust ∪ sdust → bedtools maskfasta -soft, Phase B). KegAlign honors soft-masking by skipping seeds where ≥50 % of seed bases are lowercase.

### Chain pipeline (the bridge from AXT to everything downstream)

Each AXT becomes one cleaned chain in each direction via the canonical UCSC 7-step pipeline:

```
axtChain -linearGap=loose  $axt  $src.fa  $tgt.fa  /dev/stdout |
  chainSort stdin /dev/stdout |
  chainPreNet stdin $src.sizes $tgt.sizes /dev/stdout |
  chainNet stdin $src.sizes $tgt.sizes /tmp/net /dev/null |
  netChainSubset /tmp/net /dev/stdin /dev/stdout |
  chainStitchId /dev/stdin work/01_chains/{src}.{tgt}.cleaned.chain
```

Reciprocal-best chains are then derived by swapping and re-netting:

```
chainSwap $cleaned /dev/stdout |
  chainSort stdin /dev/stdout |
  chainNet stdin $tgt.sizes $src.sizes /tmp/B_to_A.net /dev/null
netChainSubset /tmp/B_to_A.net $swapped /dev/stdout |
  chainStitchId /dev/stdin /dev/stdout |
  chainSwap /dev/stdin /dev/stdout |
  chainSort stdin work/01_chains/{A}.{B}.rbest.chain
```

`-linearGap=loose` is the right choice for divergence in the 5–15 % range that we see across the Pv strains.

### Multi-way alignments — multiz

One MAF per strain-as-hinge. For each hinge `H`, the 7 pairwise AXTs against the other strains get converted to MAF and then progressively folded into one multi-way MAF:

```
# Step 1: AXT → MAF per pair
for Q in (other 7 strains):
  axtToMaf -tPrefix=$H. -qPrefix=$Q. \
    projection/A2_kegalign/axt/${H}__vs__${Q}.axt \
    genomes/softmasked/${H}.sizes \
    genomes/softmasked/${Q}.sizes \
    work/07_multiz/${H}/${H}_vs_${Q}.maf

# Step 2: progressive multiz fold
# Order pairs by mash distance (closest to most distant)
multiz pair1.maf pair2.maf > cum_1.maf
multiz cum_1.maf pair3.maf > cum_2.maf
...
multiz cum_6.maf pair7.maf > ${H}.multiz.maf
```

The fold order matters — closest first, most distant last. We pulled the order from the mash distance matrix (`work/00_inventory/mash/dist.tsv` — built in the inventory step described in `pipeline/LOCAL.md` Section 1). Worst-case fold order drops 30–40 % of alignment blocks.

Multi-z version was 11.2 (visible in MAF block headers: `# multiz.v11.2 ...`). Wall time: ~3 hours per hinge × 8 hinges = ~24 hours serial, ~6 hours at 4 hinges concurrent.

## Outputs

| File                                            | Count | Size each      | Total  | Where                                                                   |
| ----------------------------------------------- | ----: | -------------- | ------ | ----------------------------------------------------------------------- |
| `projection/A2_kegalign/axt/{A}__vs__{B}.axt`   |    28 | ~680 MB        | 19 GB  | **Dropbox** — `Pv4_v3/A2_kegalign_axt/*.gz`                             |
| `work/01_chains/{src}.{tgt}.cleaned.chain.gz`   |    56 | ~1 MB          | 60 MB  | **Git**                                                                 |
| `work/01_chains/{src}.{tgt}.rbest.chain.gz`     |    28 | ~500 KB        | 15 MB  | **Git**                                                                 |
| `work/07_multiz/{hinge}/{hinge}_vs_{other}.maf` |    56 | ~600 MB–1.5 GB | 38 GB  | **Dropbox** — `Pv4_v3/multiz/{hinge}/*.gz` (gzip-piped during transit)  |
| `work/07_multiz/{hinge}/{hinge}.multiz.maf`     |     8 | ~3 GB          | 24 GB  | **Dropbox** — `Pv4_v3/multiz/{hinge}/*.gz`                              |
| `ucsc_hub/{ACC}/{hinge}.multiz.maf.bb`          |     8 | 14 MB–851 MB   | 4.6 GB | **Dropbox** — `Pv4_v3/ucsc_hub/{ACC}/` (also in Dropbox staging mirror) |
| Per-hinge multiz logs                           |     8 | small          | <1 MB  | **Git** — alongside the .bb files in `ucsc_hub/{ACC}/`                  |

The cleaned and rbest chains are small enough to live in the GitHub repo gzipped. Everything else — raw AXTs, pairwise MAFs, multi-way MAFs, bigMafs — sits on Dropbox.

## How these alignments fed downstream analyses

Five v3 analyses consume KegAlign and/or multiz outputs.

**Annotation projection across strains (consumes cleaned chains).** Liftoff and TOGA2/CESAR2 both consume the cleaned chains. Liftoff uses the chain to define the reference → target coordinate mapping for clean gene transfers; TOGA2/CESAR2 uses the same chain to do CDS-aware projection on the genes Liftoff marks as `needs_cesar2`. The merged outputs land in `work/02d_merged/{anchor}-as-ref/{Q}.annotation.gff3` for 4 anchor strains × 7 targets = 28 projections. Full recipe: [ORTHOLOGY.md → How we built it](ORTHOLOGY.md#how-we-built-it), Stream 1.

**Consensus orthology evidence stream 2 (consumes rbest chains).** Reciprocal-best chains feed one of the three orthology evidence streams. For each pair of strains, we ask which gene pairs share an rbest chain segment covering ≥ 90 % of both CDS. Edges from these reciprocal hits go into the union-find. The other two evidence streams are Liftoff projections (above) and graph-path co-membership (see [PANGENOME.md → How the graph fed downstream analyses](PANGENOME.md#how-the-graph-fed-downstream-analyses)). Full recipe: [ORTHOLOGY.md → How we built it](ORTHOLOGY.md#how-we-built-it), Stream 2.

**Path A2 cohort VCF projection (the canonical one — consumes PvP01→target cleaned chains).** MalariaGEN's 1,895-sample cohort VCF gets lifted from PvP01 coordinates onto each non-PvP01 reference via CrossMap, using the PvP01→{target} cleaned chains as the liftover map. **This is the canonical projection** — Path A2 won the 19-drug-resistance-gene QC over Path A1 wfmash and Path B graph-native, with a median per-gene score of 0.72 vs 0.13 (A1) and 0.45 (B). Outputs at `projection/A2_lastz/Pv4_cohort_on_{GCA}.vcf.gz`, 4.5–23 GB each, all on Dropbox. Full recipe: [MALARIAGEN_VCF_PROJECTION.md → How we did it](MALARIAGEN_VCF_PROJECTION.md#how-we-did-it).

**UCSC multi-z hub track (consumes multi-way MAFs).** Each of the 8 multi-way MAFs becomes one `bigMaf` track in the UCSC hub, accessible at `dropbox:Pv4_v3/ucsc_hub/{ACC}/{hinge}.multiz.maf.bb`. The chains also appear as chain sub-tracks under the same composite — one chain-to-each-target track per assembly. `mafToBigMaf` was bypassed (its overlap check rejects multiz's legitimate overlaps); we wrote a 50-line Python emitter that ships the raw block text as BED3+1. Hub-build recipe (no separate doc yet): `pipeline/LOCAL.md` Section 12.

**Future — phastCons / phyloP conservation tracks (would consume multi-way MAFs).** The multi-way MAFs are also the natural input for PHAST conservation scoring (phyloFit on 4dSites → phastCons HMM → bigWig). With 8 closely-related strains the signal is weak; adding *P. cynomolgi* or *P. knowlesi* as outgroups would help. Not done in v3 — flagged in the pipeline plan.

## Re-running on a different species

The recipe at `v3/pipeline/LOCAL.md` Sections 3 (KegAlign + chain pipeline) and 10 (multiz) walks through this. Quick walk-through:

1. Drop N soft-masked FASTAs into `genomes/softmasked/{strain}.fa` (the soft-masking step is described in `pipeline/LOCAL.md` Section 2 — longdust ∪ sdust → bedtools maskfasta -soft). One contig per chromosome works best.
2. Edit `pipeline/species.conf`: set `STRAINS`, `REF_STRAIN`, `N_CORES`, and `USE_GPU=1` if you have a GPU.
3. `bash pipeline/03_align_chain.sh` — runs all 28 KegAlign pairs serially (the GPU is the bottleneck, parallelism doesn't help), then the 7-stage chain pipeline for cleaned + rbest. Wall: ~2 hours.
4. `bash pipeline/10_multiz.sh` — converts AXTs to pairwise MAFs and folds them progressively per hinge. Wall: ~6 hours at 4-concurrent hinges, ~24 hours serial.

Two scale-dependent parameters:

- For closely related panels (intra-species, ≤ 5 % divergence), `--hsp_threshold 5000` is right. For inter-species (≥ 10 % divergence), drop to 4000 to keep enough seeds.
- For genomes much larger than 25 Mb, multi-z's progressive fold gets memory-heavy on the last few stages. Plan ≥ 128 GB RAM per concurrent hinge.

No GPU available — KegAlign falls back to the lastZ CPU implementation via `quay.io/biocontainers/lastz`. Same parameters, ~10× slower per pair.

The chain pipeline depends on UCSC kentUtils binaries (`axtChain`, `chainSort`, `chainPreNet`, `chainNet`, `netChainSubset`, `chainStitchId`, `chainSwap`) — these all run inside `quay.io/biocontainers/ucsc-kent-tools` so the host doesn't need them installed.

The *P. knowlesi* scaffold at `Pk/v1/pipeline/` is a working copy of this recipe for that species (7 strains, H-strain as reference).

## Galaxy tool-wrap list

For surfacing the alignment + chain + multiz build as a Galaxy workflow, these tools need wrappers (or version-bumped wrappers):

| Tool                                                            | Galaxy state             | Note                                                                                                                         |
| --------------------------------------------------------------- | ------------------------ | ---------------------------------------------------------------------------------------------------------------------------- |
| **KegAlign-GPU**                                                | ⚠ needs wrap             | The biocontainer exists, but no Galaxy wrapper — the heaviest step. Needs GPU job runner config (`gpus="auto"` requirement). |
| **lastZ** (CPU fallback)                                        | ✅ exists                 |                                                                                                                              |
| **axtChain**                                                    | ✅ exists                 | IUC; verify -linearGap flag works                                                                                            |
| **chainSort, chainPreNet, chainNet**                            | ⚠ partial                | Some subcommands are wrapped in `ucsc-tools-meta-package`; others not                                                        |
| **netChainSubset, chainStitchId, chainSwap**                    | ⚠ needs wrap             | All from kentUtils; tiny wrappers each                                                                                       |
| **axtToMaf**                                                    | ✅ exists                 | IUC                                                                                                                          |
| **multiz** (multi-way fold)                                     | ⚠ exists, refresh needed | An older Galaxy wrapper exists; verify on current Galaxy and confirm n-way input handling                                    |
| **multiz_progressive driver** (the per-hinge fold-order Python) | ⚠ needs wrap             | Custom Python helper that reads mash distances and runs multiz in the right order; small wrapper                             |
| **mafToBigMaf** (for UCSC hub publishing)                       | ⚠ needs wrap             | The kentUtils image has it; no Galaxy wrapper. We bypassed it anyway in v3 because of the overlap-block issue.               |
| **bedToBigBed, mafIndex**                                       | ✅ exists / ⚠ partial     | For hub publishing                                                                                                           |

The KegAlign-GPU wrap is the rate-limiter. Without a GPU runner in the Galaxy instance, the CPU lastZ fallback works but multiplies wall time by ~10× — fine for the 25 Mb apicomplexan scale, painful for anything larger.

The chain pipeline is best wrapped as a single subworkflow (the 7-stage piped chain build is fragile if split across separate Galaxy steps because each kentUtils binary expects `stdin` / `/dev/stdout` semantics that Galaxy collection iteration breaks). Subworkflow with internal hidden steps is the natural shape.

Multi-z, similarly, is best wrapped as a per-hinge subworkflow: collection of 7 pairwise MAFs in → one multi-way MAF out, with the fold order driven by the mash distance matrix passed as an auxiliary input.
