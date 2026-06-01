# Galaxy execution

Galaxy companion to **[PIPELINE_EXPLANATION.md](PIPELINE_EXPLANATION.md)** — the same phases (A–K), executed as a Galaxy workflow for usegalaxy.org / usegalaxy.eu. For what each phase does and why, read that doc; this one gives the Galaxy tool IDs, workflow steps, and wrappers to build for each phase. The local-container companion is **[LOCAL.md](LOCAL.md)**.

Status legend: ✅ = tool exists in the toolshed (mostly `iuc`); ⚠️ = wrapper needs writing or updating. Thirteen wrappers are new — they're collected with a priority order under [Galaxy port](#galaxy-port-wrappers-deliverables-open-issues) at the end.

End product: a reusable `.ga` workflow that produces the 27 OUTLINE essentials in a published Galaxy history.

## Inputs

```yaml
class: GalaxyWorkflow
name: pangenome_selection_pipeline
inputs:
  assemblies:     { type: collection, collection_type: list, label: "N genome FASTAs (one per strain)" }
  annotations:    { type: collection, collection_type: list, label: "N annotation GFF3s, parallel to assemblies" }
  cohort_vcfs:    { type: collection, collection_type: list, label: "Per-chromosome cohort VCFs (multi-sample)" }
  ref_strain:     { type: text,  label: "Reference strain name (one of the strain names)" }
  anchor_strains: { type: text,  label: "Comma-separated anchor strains for annotation projection" }
  chrom_rename:   { type: data,  format: tabular, label: "Chromosome rename map (OLD<TAB>NEW); see PIPELINE_EXPLANATION" }
```

Collection naming is the binding key: assemblies and annotations are parallel lists indexed by strain name. The species-specific knobs (`ref_strain`, `anchor_strains`, BUSCO lineage, the multigene-family regex) are workflow-form parameters — the same `species.conf` fields LOCAL.md sets, surfaced in the run form.

## Phase A — Inventory

| Tool | Status | Notes |
|---|---|---|
| `sourmash_sketch` | ✅ iuc | k=31, scaled=1000; one `.sig` per assembly. Same sketches the BRC catalog ingests. |
| `sourmash_compare` | ✅ iuc | N×N similarity matrix (CSV) |
| `busco` | ✅ iuc | lineage is species-configurable |

```yaml
- id: sourmash_sketch
  tool_id: toolshed.g2.bx.psu.edu/repos/iuc/sourmash_sketch/sourmash_sketch
  in: { sequences: assemblies }
- id: sourmash_compare
  tool_id: toolshed.g2.bx.psu.edu/repos/iuc/sourmash_compare/sourmash_compare
  in: { signatures: sourmash_sketch/signatures }
- id: busco
  tool_id: busco
  in: { sequences: assemblies, lineage: <clade>_odb10 }
```

No wrappers needed.

## Phase B — Soft-masking

| Tool | Status | Notes |
|---|---|---|
| `longdust` | ⚠️ wrap | single static binary; new wrapper |
| `sdust` | ✅ iuc | within the minimap2 suite |
| `bedtools_merge` | ✅ iuc | union of the two interval sets |
| `bedtools_maskfastabed` | ✅ iuc | `soft_mask: true` |

```yaml
# parallel per strain via map-over
- id: longdust   { tool_id: longdust,   map_over: assemblies }
- id: sdust      { tool_id: sdust,      map_over: assemblies }
- id: union_bed  { tool_id: bedtools_merge,        map_over: [longdust, sdust] }
- id: mask       { tool_id: bedtools_maskfastabed, map_over: [assemblies, union_bed], params: { soft_mask: true } }
```

Wrappers: **longdust**.

## Phase C — Pairwise alignment and chains

### C.1–3, alignment → chains

| Tool | Status | Notes |
|---|---|---|
| `kegalign_gpu` | ⚠️ wrap | biggest compute step; needs a GPU job runner |
| `lastz` | ✅ iuc | CPU fallback (~10× slower) |
| `axtChain`, `chainSort`, `chainPreNet`, `chainNet`, `netChainSubset`, `chainStitchId`, `chainSwap` | ⚠️ wrap | UCSC kentUtils; some in toolshed at inconsistent versions |

```yaml
- id: pair_collections        # tiny helper: build NxN pair-collections from the strain list
  tool_id: __pair_strains__
  in: { masked: mask }
  params: { include_self: false, both_directions: true }
- id: kegalign
  tool_id: kegalign_gpu
  map_over: pair_collections
  requires: gpu
- id: axt_to_chain  { tool_id: axtChain, map_over: kegalign }
- id: chain_net_clean { tool_id: __chain_pipeline_clean__, map_over: axt_to_chain }  # ⭐ cleaned.chain
- id: chain_rbest     { tool_id: __chain_pipeline_rbest__, map_over: axt_to_chain }  # ⭐ rbest.chain
```

Wrappers: **kegalign_gpu**, the **chain pipeline tools** (six kentUtils), and the **`__pair_strains__`** helper. `__chain_pipeline_clean__` / `__chain_pipeline_rbest__` bundle the canonical UCSC chain stages (see PIPELINE_EXPLANATION C.1–3) as one logical tool each.

### C.4, annotation projection

| Tool | Status | Notes |
|---|---|---|
| `liftoff` | ✅ iuc | the fast-path projector |
| `phase_c2_triage` | ⚠️ wrap | routes lifted genes into clean / rescue / flag |
| `toga2` / `cesar2` | ⚠️ wrap | rescue engine; CLI exists, no Galaxy tool yet |
| `agat_sp_merge_annotations` | ✅ iuc | merge step (substitutes our `phase_c4_merge.py`) |

```yaml
- id: project_per_anchor
  subworkflow: project_one_anchor      # liftoff → triage → toga2 → merge
  map_over: anchor_strains_collection
  in: { anchor_strain: ..., anchor_annot: ..., query_strains: ..., chains: chain_net_clean }
```

Wrappers: **toga2/cesar2**, **phase_c2_triage**, **merge** (or use `agat_sp_merge_annotations`). For closely-related strains Liftoff covers most genes, so TOGA2 is the lowest-priority wrapper — the subworkflow can ship Liftoff-only first and add the rescue branch later.

## Phase D — PGGB graph

| Tool | Status | Notes |
|---|---|---|
| `pansn_rename` | ⚠️ wrap | five-line helper (`SAMPLE#1#CONTIG`) |
| `__concat_seqs__` | ⚠️ wrap | tiny helper: concat masked FASTAs into one multifasta |
| `pggb` | ✅ iuc | one big multifasta + `n` |

```yaml
- id: pansn_rename { tool_id: pansn_rename, map_over: mask }
- id: concat_fasta { tool_id: __concat_seqs__, in: { fastas: pansn_rename } }
- id: pggb
  tool_id: pggb
  in: { sequences: concat_fasta, n: <N>, s: 5000, p: 90, k: 23 }
```

Wrappers: **pansn_rename**, **`__concat_seqs__`**.

## Phase E — Consensus orthology

| Tool | Status | Notes |
|---|---|---|
| `phase_e_consensus` | ⚠️ wrap | port of `scripts/phase_e_consensus.py`; one Python script |

```yaml
- id: consensus
  tool_id: phase_e_consensus
  in:
    annotations: project_per_anchor/annotation_gff3   # projection edges
    chains:      chain_rbest                           # reciprocal-best edges
    graph:       pggb/og                               # graph co-membership edges (via odgi paths)
```

Wrappers: **phase_e_consensus**. `odgi paths` is bundled in the pggb container, so the graph-edge extraction can live inside this tool.

## Phase F — Codon and protein MSAs

| Tool | Status | Notes |
|---|---|---|
| `gffread` | ✅ iuc | pull CDS per member |
| `mafft` | ✅ iuc | configure LINSI (`--localpair --maxiterate 1000`) |
| `pal2nal` | ✅ iuc | codon back-translation |
| `trimal` | ✅ iuc | `-automated1` cleaned siblings |

```yaml
- id: msa_strict  { subworkflow: build_msa_pipeline, in: { min_intact: 7, ortho: consensus, mode: codon } }
- id: msa_relaxed { subworkflow: build_msa_pipeline, in: { min_intact: 5, ortho: consensus, mode: codon } }
```

Wrappers: **build_msa_pipeline** subworkflow — wraps gffread + mafft + pal2nal + trimal as one logical unit, including the strip-internal-stops step (see PIPELINE_EXPLANATION F).

## Phase G — IQ-TREE

| Tool | Status | Notes |
|---|---|---|
| `iqtree` | ✅ iuc | iqtree2 exists; bump to iqtree3 |

```yaml
- id: trees_strict  { tool_id: iqtree, map_over: msa_strict,  params: { model: MFP, bootstrap: 1000 } }
- id: trees_relaxed { tool_id: iqtree, map_over: msa_relaxed, params: { model: MFP, bootstrap: 1000 } }
```

Wrappers: **iqtree3** version bump. The <4-unique-sequence fallback (drop bootstrap) should be handled inside the wrapper or a conditional, since the bootstrap hang is otherwise silent.

## Phase H — HyPhy BUSTED

| Tool | Status | Notes |
|---|---|---|
| `hyphy_busted` | ✅ iuc / datamonkey | `--srv No --branches All` |

```yaml
- id: hyphy_strict  { tool_id: hyphy_busted, map_over: msa_strict }
- id: hyphy_relaxed { tool_id: hyphy_busted, map_over: msa_relaxed }
```

No new wrapper, but **bundle** the per-gene JSONs into a single collapsed collection element — Galaxy histories choke past ~10k datasets, and the relaxed set alone is ~4,200 genes.

## Phase I — Multiz MAFs

| Tool | Status | Notes |
|---|---|---|
| `axtToMaf` | ✅ iuc | part of UCSC tools |
| `multiz` | ⚠️ verify | an old toolshed wrapper exists; check it runs on current Galaxy |

```yaml
- id: multiz_per_hinge
  subworkflow: multiz_hinge_subworkflow   # axtToMaf → progressive fold, order from sourmash similarity
  map_over: strain_list
  in: { hinge_strain: ..., pairwise_axts: kegalign }
```

Wrappers: **multiz-modern** — verify/refresh the existing wrapper. Fold order comes from the Phase A `compare.csv` (closest = highest similarity).

## Phase J — Cohort VCF projection

| Tool | Status | Notes |
|---|---|---|
| `bcftools_annotate` | ✅ iuc | `--rename-chrs` |
| `crossmap_vcf` | ✅ iuc | cleaned chain, target FASTA |
| `bcftools_sort`, `bcftools_concat` | ✅ iuc | sort before index; concat `-a` |

```yaml
- id: rename_cohort
  tool_id: bcftools_annotate
  map_over: cohort_vcfs
  params: { rename_chrs: chrom_rename }
- id: crossmap_per_target
  subworkflow: crossmap_target           # crossmap → bcftools sort → index, per chr
  map_over: target_strain_collection
  in: { chains: chain_net_clean, vcfs: rename_cohort }
- id: concat_per_target
  tool_id: bcftools_concat
  map_over: crossmap_per_target/per_chr
```

No new wrappers. Use the **cleaned** chain (not rbest), and keep the sort-before-index ordering — same constraints as LOCAL.md.

## Phase K — UCSC track hub

A separate `pangenome-publish-ucsc-hub` workflow that consumes the outputs of the build/project/selection workflows + the chain bundle and emits a `hub.txt` collection ready to push to a range-capable host (e.g. `hgdownload.soe.ucsc.edu/hubs/BRC/`). **Dropbox links do not work as a hub host** — UCSC needs HTTP byte-range support.

| Tool | Status | Notes |
|---|---|---|
| `mafToBigMaf` | ⚠️ wrap | kentUtils; verify subcommand is wrapped |
| `mafIndex` | ⚠️ wrap | likely missing |
| `bedToBigBed` | ✅ iuc | |
| `gff3ToGenePred`, `genePredToBed` | ✅ iuc | |
| `faToTwoBit` | ✅ iuc | `.2bit` referenced from `genomes.txt` |
| `chain_to_bigChain` | ⚠️ wrap | custom: chain → bigChain.bed (6+6) + bigLink.bed (4+1). Hubs need `type bigChain`, not `.chain.gz`. |
| `build_selection_bigbed` | ⚠️ wrap | custom: BUSTED JSONs + ortholog table + BED12 → BED12+5 |
| `build_orthogroup_bigbed` | ⚠️ wrap | custom: ortholog table → BED12 |
| `build_trackdb` | ⚠️ wrap | custom: per-assembly `trackDb.txt`; bigMaf and bigChain **cannot** share a composite |
| `build_genomes_txt` | ⚠️ wrap | custom: 9-field records; `defaultPos` must be a real `chrN:start-end` |
| `hubCheck` | ⚠️ wrap | validate before publish |

```yaml
- id: maf_to_bigmaf  { tool_id: mafToBigMaf, map_over: multiz_alignments }   # one per hinge
- id: maf_index      { tool_id: mafIndex,    map_over: multiz_alignments }
- id: fa_to_2bit     { tool_id: faToTwoBit,  map_over: assemblies }
- id: chain_to_bigchain { tool_id: chain_to_bigChain, map_over: chain_files }     # 56 pairs (8×7)
- id: bigchain_bed_to_bb { tool_id: bedToBigBed, map_over: chain_to_bigchain/bigChain_bed }  # -type=bed6+6 -as=bigChain.as
- id: biglink_bed_to_bb  { tool_id: bedToBigBed, map_over: chain_to_bigchain/bigLink_bed }   # -type=bed4+1 -as=bigLink.as
- id: gff_to_bed12   { tool_id: gff3ToGenePred + genePredToBed (subworkflow), map_over: merged_annotations }
- id: selection_bigbed   { tool_id: build_selection_bigbed, in: { hyphy_jsons: ..., ortholog_table: consensus, ref_bed12: ... } }
- id: orthogroup_bigbed  { tool_id: build_orthogroup_bigbed }
- id: trackdb_per_assembly { tool_id: build_trackdb,     map_over: assemblies }
- id: genomes_txt        { tool_id: build_genomes_txt }
- id: hubcheck           { tool_id: hubCheck }
- id: publish            { tool_id: rsync_or_rclone_to_datacache }   # or emit a tarball for manual rsync
```

The hub layout needs per-assembly directories: encode this with a `list:list` collection (outer = assembly accession, inner = track type), and have `build_trackdb` take the outer collection name as the accession.

## How it runs

- **Collection iteration** is what makes the NxN-pairwise and per-anchor loops natural; the cost is a few tiny helper tools (`__pair_strains__`, `__concat_seqs__`) to build pair/concat collections from the input list.
- **GPU**: the `kegalign_gpu` step needs `gpus="auto"` and a GPU job runner; config differs between usegalaxy.org and usegalaxy.eu.
- **Dataset count**: bundle per-gene HyPhy and tree outputs into collapsed collections — histories degrade past ~10k datasets.
- **Skip vs v3**: don't port the Path B graph-native VCF projection (Path A2 won the drug-resistance QC and is the production path); skip GENESPACE (v3 fell back to OrthoFinder3, whose wrapper exists — or rely on Phase E consensus alone). Both cuts reduce wrapper work for v1.
- Wall time on a GPU + 32-core Galaxy instance: ~24 h, same as local.

## New species

1. Upload N FASTAs + N GFF3s into a new history (Upload tool with collection-naming rules).
2. Upload the per-chromosome cohort VCFs (collection) and `chrom_rename.tsv`.
3. Import the published `pangenome_selection_pipeline.ga`.
4. Run workflow — Galaxy maps inputs by collection.
5. Set species params in the form (`ref_strain`, `anchor_strains`, BUSCO lineage, family regex; for AT-rich genomes lower the KegAlign HSP threshold).
6. Submit.

Everything else (tool versions, the workflow graph) is unchanged — the same parameterization LOCAL.md's `species.conf` carries, surfaced on the run form.

## Galaxy port: wrappers, deliverables, open issues

**Wrappers to write (13)**, highest impact first:

1. **kegalign_gpu** — GPU-required, biggest compute step; without it Galaxy falls back to lastZ-CPU.
2. **chain pipeline tools** — axtChain, chainNet, chainSwap, netChainSubset, chainStitchId, chainSort, chainPreNet (some in toolshed at inconsistent versions).
3. **phase_e_consensus** — direct port of the Python script.
4. **toga2/cesar2** — most complex, needs a container; lowest priority since Liftoff covers most genes for closely-related strains.
5. **longdust** — single binary.
6. **pansn_rename** + **`__concat_seqs__`** + **`__pair_strains__`** — five-line helpers.
7. **multiz-modern** — verify/refresh the old wrapper.
8. **build_msa_pipeline** subworkflow — gffread + mafft + pal2nal + trimal.
9. **merge_annotation** (or use `agat_sp_merge_annotations`) and **phase_c2_triage**.
10. Phase K: **chain_to_bigChain, build_selection_bigbed, build_orthogroup_bigbed, build_trackdb, build_genomes_txt, mafToBigMaf, mafIndex, hubCheck**.

**Deliverables**

1. `pipeline/galaxy/workflow.ga` — exported workflow JSON
2. `pipeline/galaxy/tool_wrappers/*.xml` — the new wrappers
3. `pipeline/galaxy/test_history.tar` — small smoke-test history (3 strains, 1 chr) for CI
4. `pipeline/galaxy/SHED.md` — toolshed submission checklist
5. `pipeline/galaxy/README.md` — usegalaxy.org user-facing instructions

**Open issues**

1. KegAlign wrapper — GPU runner config differs per server; needs `gpus="auto"` and tests on both usegalaxy.org and usegalaxy.eu.
2. TOGA2 — wrap, or ship Liftoff-only for v1?
3. multiz — last toolshed update is old; may need a maintenance pass.
4. Collection-naming conventions vs our `{strain}` glob — needs a helper to convert.
5. Where does the wrapped pipeline live — main toolshed, or a satellite project (like `usegalaxy.org/duplex`)?
