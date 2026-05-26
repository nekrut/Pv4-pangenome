# Galaxy pipeline — workflow + tool wrappers

Same 27 OUTLINE essentials, executed as a Galaxy workflow. Designed for usegalaxy.org / usegalaxy.eu deployment.

## Goal

A reusable `.ga` workflow file that takes:
- N genome FASTAs (collection)
- N annotation GFF3s (collection, parallel to FASTAs)
- per-chromosome cohort VCFs (collection)
- A small species-config YAML

and produces the 27 essentials in a published Galaxy history.

## Tool inventory

| Phase | Galaxy tool ID | Status | Notes |
|---|---|---|---|
| A | `toolshed.g2.bx.psu.edu/.../mash_sketch` + `mash_screen` | ✅ exists | iuc, IUC-maintained |
| A | `quay.io/.../busco` | ✅ exists | iuc |
| B | `bedtools_maskfastabed` | ✅ exists | iuc |
| B | `longdust` | ⚠️ wrap | new wrapper needed; binary is single static exe |
| B | `sdust` (from minimap2 tooling) | ✅ exists (within minimap2 suite) | iuc |
| C | `kegalign-gpu` | ⚠️ wrap | new wrapper needed; existing biocontainer; needs GPU job runner config |
| C | `lastz` | ✅ exists | iuc; CPU fallback if no GPU |
| C | `axtChain`, `chainNet`, `chainSwap`, `netChainSubset`, `chainStitchId`, `chainSort`, `chainPreNet` | ⚠️ wrap | UCSC tools; some wrappers exist (e.g. `axtChain` in toolshed), others missing |
| C | `liftoff` | ✅ exists | iuc |
| C | `toga2` / `cesar2` | ⚠️ wrap | new wrapper; TOGA2 has a CLI but no Galaxy tool yet |
| C | `agat_sp_merge_annotations` | ✅ exists | for the merge step (substitutes our merge_annotation.py) |
| D | `pggb` | ✅ exists | iuc; takes one big multifasta + n |
| D | `seqkit_pansn_rename` or custom | ⚠️ wrap | small helper |
| E | (custom: phase E consensus) | ⚠️ wrap | port `phase_e_consensus.py` to Galaxy tool |
| F | `mafft` | ✅ exists | iuc; configure LINSI mode |
| F | `pal2nal` | ✅ exists | iuc |
| F | `trimal` | ✅ exists | iuc |
| F | `gffread` | ✅ exists | iuc |
| G | `iqtree3` | ✅ exists (iqtree2 exists; iqtree3 needs version bump) | iuc |
| H | `hyphy_busted` | ✅ exists | datamonkey wraps; also iuc |
| I | `multiz` (chained MAF inputs) | ⚠️ exists ? | older Galaxy wrapper exists; verify compatibility |
| I | `axtToMaf` | ✅ exists | iuc (part of UCSC tools) |
| J | `bcftools_annotate` | ✅ exists | iuc, supports --rename-chrs |
| J | `crossmap_vcf` | ✅ exists | iuc |
| J | `bcftools_sort`, `bcftools_concat` | ✅ exists | iuc |

**Wrappers to write (13)**: longdust, kegalign-gpu, axtChain (re-wrap), chainNet, chainSwap, netChainSubset, chainStitchId, toga2/cesar2, pansn_rename, phase_e_consensus, multiz-modern, build_msa_pipeline (gffread+mafft+pal2nal chained), merge_annotation.

## Workflow structure (`.ga` skeleton)

```yaml
class: GalaxyWorkflow
name: pangenome_selection_pipeline
inputs:
  assemblies:
    type: collection
    collection_type: list:paired   # OR list of fastas
    label: "N genome FASTAs (one per strain)"
  annotations:
    type: collection
    collection_type: list
    label: "N annotation GFF3s, parallel to assemblies"
  cohort_vcfs:
    type: collection
    collection_type: list
    label: "Per-chromosome cohort VCFs (multi-sample)"
  ref_strain:
    type: text
    label: "Reference strain name (one of the strain names)"
  chrom_rename:
    type: data
    format: tabular
    label: "Chromosome rename map (OLD<TAB>NEW)"
  anchor_strains:
    type: text
    label: "Comma-separated anchor strains for annotation projection"

steps:
  # Phase A
  - id: mash_sketch
    tool_id: toolshed.g2.bx.psu.edu/repos/iuc/mash/mash_sketch
    in:
      sequences: { source: assemblies }
  - id: mash_dist
    tool_id: mash_dist
    in: { sketch: mash_sketch/sketch }
  - id: busco
    tool_id: busco
    in:
      sequences: { source: assemblies }
      lineage: plasmodium_odb10           # species-configurable

  # Phase B (parallel per strain via map-over)
  - id: longdust
    tool_id: longdust
    map_over: assemblies
  - id: sdust
    tool_id: sdust
    map_over: assemblies
  - id: union_bed
    tool_id: bedtools_merge
    map_over: [longdust, sdust]
  - id: mask
    tool_id: bedtools_maskfastabed
    map_over: [assemblies, union_bed]
    params: { soft_mask: true }

  # Phase C — KegAlign + chain pipeline (pairwise)
  - id: pair_collections
    tool_id: __pair_strains__
    in: { masked: mask }
    params: { include_self: false, both_directions: true }
  - id: kegalign
    tool_id: kegalign_gpu
    map_over: pair_collections
    requires: gpu
  - id: axt_to_chain
    tool_id: axtChain
    map_over: kegalign
  - id: chain_net_clean
    tool_id: __chain_pipeline_clean__
    map_over: axt_to_chain
  - id: chain_rbest
    tool_id: __chain_pipeline_rbest__
    map_over: axt_to_chain
  # ⭐ outputs: cleaned.chain, rbest.chain

  # Phase C — annotation projection (per anchor)
  - id: project_per_anchor
    subworkflow: project_one_anchor
    map_over: anchor_strains_collection
    in:
      anchor_strain: ...
      anchor_annot: ...
      query_strains: ...

  # Phase D — PGGB
  - id: pansn_rename
    tool_id: pansn_rename
    map_over: mask
  - id: concat_fasta
    tool_id: __concat_seqs__
    in: { fastas: pansn_rename }
  - id: pggb
    tool_id: pggb
    in:
      sequences: concat_fasta
      n: 8                              # ${#STRAINS[@]}
      s: 5000
      p: 90

  # Phase E
  - id: consensus
    tool_id: phase_e_consensus
    in:
      annotations: project_per_anchor/annotation_gff3
      chains: chain_net_clean
      graph: pggb/og

  # Phase F (sub-workflow per orthogroup, with min_intact filter)
  - id: msa_strict
    subworkflow: build_msa_pipeline
    in: { min_intact: 7, ortho: consensus, mode: codon }
  - id: msa_relaxed
    subworkflow: build_msa_pipeline
    in: { min_intact: 5, ortho: consensus, mode: codon }

  # Phase G — IQ-TREE per gene
  - id: trees_strict
    tool_id: iqtree
    map_over: msa_strict
    params: { model: MFP, bootstrap: 1000 }
  - id: trees_relaxed
    tool_id: iqtree
    map_over: msa_relaxed

  # Phase H — HyPhy bulk
  - id: hyphy_strict
    tool_id: hyphy_busted
    map_over: msa_strict
  - id: hyphy_relaxed
    tool_id: hyphy_busted
    map_over: msa_relaxed

  # Phase I — Multiz
  - id: multiz_per_hinge
    subworkflow: multiz_hinge_subworkflow
    map_over: strain_list
    in: { hinge_strain: ..., pairwise_axts: kegalign }

  # Phase J — VCF projection
  - id: rename_cohort
    tool_id: bcftools_annotate
    map_over: cohort_vcfs
    params: { rename_chrs: chrom_rename }
  - id: crossmap_per_target
    subworkflow: crossmap_target
    map_over: target_strain_collection
    in: { chains: chain_net_clean, vcfs: rename_cohort }
  - id: concat_per_target
    tool_id: bcftools_concat
    map_over: crossmap_per_target/per_chr
```

Galaxy collection patterns make pairwise (NxN) and per-anchor (per-strain) iteration natural, but require some custom helper tools (`__pair_strains__`, `__concat_seqs__`) to build the pair-collections from the input list. These are tiny wrappers.

## Wrapper-development priority

Order of effort (highest impact first):

1. **kegalign_gpu** — GPU-required, biggest compute step; without it Galaxy falls back to lastZ-CPU (10× slower)
2. **chain pipeline tools** (axtChain, chainNet, chainSwap, netChainSubset, chainStitchId) — six small wrappers, all UCSC kentUtils; some already in toolshed but inconsistent versions
3. **phase_e_consensus** — direct port of `scripts/phase_e_consensus.py`; one Python script in a wrapper
4. **toga2** — most complex, requires container; lower priority since Liftoff alone covers most genes for closely-related strains
5. **longdust** — single binary; one weekend's work
6. **pansn_rename** — five-line Python script
7. **multiz-modern** — verify existing wrapper still works on newer Galaxy; refresh if broken
8. **build_msa_pipeline subworkflow** — wraps gffread + mafft + pal2nal + trimal as one logical unit

## UCSC track-hub publishing (Phase K — new section)

Mirrors LOCAL.md section 12. As a Galaxy workflow, this is a separate `pangenome-publish-ucsc-hub` workflow that takes the outputs of `pangenome-build-and-project` + `hyphy-selection-screen` + the chain bundle, and produces a `hub.txt` collection ready for rsync to `hgdownload.soe.ucsc.edu/hubs/BRC/`.

### Galaxy tool wrappers needed (additions)

| Tool | Status | Notes |
|---|---|---|
| `mafToBigMaf` | ⚠️ wrap | Part of kentUtils; existing `ucsc-kent-tools` wrappers in toolshed but not all subcommands. Verify `mafToBigMaf` is wrapped or extend. |
| `bedToBigBed` | ✅ exists | iuc |
| `mafIndex` | ⚠️ wrap | Likely missing |
| `gff3ToGenePred`, `genePredToBed` | ✅ exists (part of UCSC tools) | iuc |
| `chain_to_bigChain` | ⚠️ wrap | Custom Python tool: parses chain → bigChain.bed (6+6) + bigLink.bed (4+1) for bedToBigBed. UCSC track hubs need `type bigChain`, not text `.chain.gz`. |
| `faToTwoBit` | ✅ exists | iuc — needed for the `.2bit` files referenced from `genomes.txt`. |
| `build_selection_bigbed` | ⚠️ wrap | Custom Python tool: BUSTED JSONs + ortholog table + BED12 → BED12+5 |
| `build_orthogroup_bigbed` | ⚠️ wrap | Similar custom Python tool |
| `build_trackdb` | ⚠️ wrap | Custom Python tool: emits per-assembly `trackDb.txt`. Output layout: 1 standalone bigMaf + bigChain composite + bigBed12 composite (+ selection composite on reference strain). bigMaf and bigChain cannot share a composite. |
| `build_genomes_txt` | ⚠️ wrap | Custom Python tool: emits 9-field `genomes.txt` entries (genome, trackDb, groups, description, twoBitPath, organism, defaultPos, scientificName, htmlPath). `defaultPos` must be a real `chrN:start-end` per assembly. |
| `hubCheck` | ⚠️ wrap | Validates hub structure before publishing |

### Workflow structure

```yaml
steps:
  - id: maf_to_bigmaf
    tool_id: mafToBigMaf
    map_over: multiz_alignments    # one per hinge
  - id: maf_index
    tool_id: mafIndex
    map_over: multiz_alignments

  - id: fa_to_2bit
    tool_id: faToTwoBit
    map_over: assemblies

  - id: chain_to_bigchain
    tool_id: chain_to_bigChain        # emits .bigChain.bed + .bigLink.bed
    map_over: chain_files             # 56 pairs (8 × 7)
  - id: bigchain_bed_to_bb
    tool_id: bedToBigBed              # -type=bed6+6 -as=bigChain.as
    map_over: chain_to_bigchain/bigChain_bed
  - id: biglink_bed_to_bb
    tool_id: bedToBigBed              # -type=bed4+1 -as=bigLink.as
    map_over: chain_to_bigchain/bigLink_bed

  - id: gff_to_bed12
    tool_id: gff3ToGenePred + genePredToBed (subworkflow)
    map_over: merged_annotations

  - id: selection_bigbed
    tool_id: build_selection_bigbed
    in:
      hyphy_jsons: collection from hyphy-selection-screen workflow
      ortholog_table: phase E output
      ref_bed12: from inputs

  - id: orthogroup_bigbed
    tool_id: build_orthogroup_bigbed

  - id: trackdb_per_assembly
    tool_id: build_trackdb            # 1 standalone bigMaf + bigChain composite + bigBed12 composite
    map_over: assemblies

  - id: genomes_txt
    tool_id: build_genomes_txt        # 9-field per-assembly record + defaultPos

  - id: hubcheck
    tool_id: hubCheck

  - id: publish
    tool_id: rsync_or_rclone_to_datacache  # or just emit a tarball for manual rsync
```

### Galaxy collection trick

The hub layout requires per-assembly directories. Galaxy collections can encode this via `list:list` (outer = assembly, inner = track type). The wrapper for `build_trackdb` should accept the outer collection name as the assembly accession.

## Test on new species via Galaxy

1. Upload N FASTAs + N GFF3s into a new history (use the "Upload" tool with rules for collection naming)
2. Upload cohort VCFs (collection of per-chr .vcf.gz)
3. Upload chrom_rename.tsv
4. Import the published `pangenome_selection_pipeline.ga` workflow
5. Click "Run workflow" — Galaxy maps inputs by collection
6. Configure species params via the workflow form (ref_strain, anchor_strains, lineage for BUSCO)
7. Submit

Wall time on a Galaxy instance with GPU + 32-core: ~24 hrs (same as local).

## What to do differently from v3

- **Don't replicate the Path B graph-native VCF projection** in the Galaxy version yet — v3 found A2 KegAlign chains (Path A2) won the drug-resistance QC. Path A2 is the production path. Skip Path B in Galaxy v1 to reduce wrapper complexity; add later if needed.
- **Skip GENESPACE** — v3 substituted OrthoFinder3 after GENESPACE failures. OrthoFinder3 wrapper exists in toolshed. Use that, or use Phase E consensus alone.
- **Bundle the per-gene HyPhy results** as a single collection element rather than thousands of individual datasets — Galaxy histories choke on >10k datasets. Use the "collapse collection" pattern.

## Deliverables for Galaxy port

1. `pipeline/galaxy/workflow.ga` — exported workflow JSON
2. `pipeline/galaxy/tool_wrappers/*.xml` — 13 new wrappers
3. `pipeline/galaxy/test_history.tar` — small smoke-test history (3 strains, 1 chr) for CI
4. `pipeline/galaxy/SHED.md` — toolshed submission checklist for the new wrappers
5. `pipeline/galaxy/README.md` — usegalaxy.org user-facing instructions

## Open issues / unresolved

1. KegAlign Galaxy wrapper — GPU job runner config differs per server; needs `gpus="auto"` requirement and tests on both usegalaxy.org and usegalaxy.eu
2. TOGA2 wrapper — should we wrap it or skip TOGA2 in the Galaxy version and rely on Liftoff only?
3. multiz — last toolshed update is old; may need a maintenance pass before reuse
4. Workflow collection-naming conventions — Galaxy collection rules vs our `{strain}` glob — needs a helper tool to convert
5. Where does the wrapped pipeline live? — usegalaxy.org main toolshed, or a satellite project (similar to `usegalaxy.org/duplex`)?
