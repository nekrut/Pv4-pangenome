# Pangenome

The 8-strain *Plasmodium vivax* pangenome graph was built in v2 with PGGB and is the structural foundation for several v3 analyses.

**Related documents.** [MULTIZ.md](MULTIZ.md) covers the parallel KegAlign pairwise alignment + multiz multi-way layer; [ORTHOLOGY.md](ORTHOLOGY.md) covers the consensus ortholog table that consumes graph paths as one of its evidence streams; [MALARIAGEN_VCF_PROJECTION.md](MALARIAGEN_VCF_PROJECTION.md) covers the population-cohort variant projection (the graph-native branch is mentioned there as a non-canonical comparator).

## How we built the graph

Input — 8 haploid assemblies, PanSN-renamed (`SAMPLE#1#CONTIG`) and concatenated into one bgzipped multifasta:

| Strain | Accession       | Source                                   |
| ------ | --------------- | ---------------------------------------- |
| PvP01  | GCA_900093555.2 | Papua, Indonesia — modern reference      |
| Sal-I  | GCA_000002415.2 | Salvador I, El Salvador (monkey-adapted) |
| PvW1   | GCA_914969965.1 | Mauritania                               |
| PAM    | GCA_949152365.1 | Peruvian Amazon                          |
| PvSY56 | GCA_003402215.1 | Asian field isolate (Auburn PCR-free)    |
| PvT01  | GCA_900093545.1 | Thai field isolate (Auburn)              |
| PvC01  | GCA_900093535.1 | Cambodian field isolate (Auburn)         |
| MHC087 | GCA_040114635.1 | Recent field isolate, 2024               |

Build command (from `v2/pggb_out/*.params.yml`):

```
pggb -i pggb_in/pggb_input.fa.gz -o pggb_out -s 5000 -p 90 -n 8 -k 23 -t 18 -V GCA_900093555.2 -Y '#'
```

The flags:

- `-s 5000` — wfmash segment length, 5 kb
- `-p 90` — wfmash mapping identity, 90 %
- `-n 8` — number of haplotypes (one per strain, all haploid)
- `-k 23` — seqwish min-match length, 23 bp
- `-V GCA_900093555.2` — call vg deconstruct against PvP01 to emit a per-strain VCF
- `-Y '#'` — PanSN delimiter
- `-t 18` — 18 threads

Internally: wfmash → seqwish → smoothxg (POA-aware smoothing, block-id-min 0.9) → gfaffix → odgi viz/layout. Wall time on the build was ~50 minutes.

PGGB version was the `ghcr.io/pangenome/pggb` image as of mid-May 2026 (the params.yml header carries the full digest).

## Outputs

All produced under `v2/pggb_out/`. The 75-MB gzipped GFA stays in the GitHub repo; the binary `.og` and `.og.lay` go to Dropbox.

| File                                                         |        Size | What it is                                                                                                                                                | Where                             |
| ------------------------------------------------------------ | ----------: | --------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------- |
| `pggb_input.fa.gz.f705205.c28ecf8.dcad0e6.smooth.fix.gfa.gz` |       75 MB | The canonical pangenome graph in GFA1 (gzipped).                                                                                                          | **Git** — `v2/pggb_out/`          |
| `*.smooth.fix.og`                                            |      654 MB | odgi binary of the same graph (for fast queries).                                                                                                         | **Dropbox** — `Pv4_v3/pggb_8way/` |
| `*.smooth.fix.og.lay`                                        |      109 MB | odgi layout for visualization.                                                                                                                            | **Dropbox** — `Pv4_v3/pggb_8way/` |
| `*.smooth.fix.affixes.tsv.gz`                                |        20 B | Empty-ish affix-collapse log (no shared affixes called).                                                                                                  | Git                               |
| `*.smooth.fix.GCA_900093555.2.vcf`                           |         0 B | Empty VCF from the `-V` flag (deconstruct issued by pggb itself returned no records — we ran `vg deconstruct` per-target directly in Phase J, see below). | Git                               |
| `*.smooth.fix.og.lay.draw.png`, `*_multiqc.png`, viz panels  | ~2 MB total | Layout PNGs and odgi viz panels (multiQC-ready).                                                                                                          | Git                               |
| `*.params.yml`, `*.log`                                      |       ~1 MB | Run parameters and build log.                                                                                                                             | Git                               |
| `pggb_input.fa.gz` (the PanSN-renamed multifasta input)      |       74 MB | Source FASTA fed into pggb.                                                                                                                               | Git (in `pggb_in/`)               |

## How the graph fed downstream analyses

Three v3 analyses consume the graph directly.

**Consensus orthology — graph-path co-membership as evidence stream 3.** One of three orthology evidence sources is graph-path co-membership. For each pair of strains, we ask via `odgi paths` which gene CDS regions traverse the same set of nodes, with at least 90 % of CDS overlap on a shared path. The other two sources are Liftoff/TOGA2 annotation projections and reciprocal-best chain hits. Union-find over the three edge sets produced the consensus table at `work/03_consensus/ortholog_table.tsv` — 7,504 orthogroups, with `graph_strains` and `graph_mean_pav` columns added by the next step. See [ORTHOLOGY.md → How we built it](ORTHOLOGY.md#how-we-built-it) for the full three-stream merge.

**Graph PAV cross-validation.** For every CORE-1:1 orthogroup we extract the PvP01 gene's BED region and call `odgi pav -S` (sample-grouped presence/absence). The output `work/04_graph_validation/core_1to1.pav.tsv` (33,489 rows × 6 columns) tells us, per orthogroup, what fraction of each strain's path actually traverses the gene's nodes. 4,048 of 4,186 CORE-1:1 genes (96.7 %) are graph-traversed in all 8 strains at PAV ≥ 0.5. This is the QC that promotes orthology calls from "annotation-derived" to "annotation + graph-confirmed". See [ORTHOLOGY.md → Outputs](ORTHOLOGY.md#outputs) for the `graph_strains` and `graph_mean_pav` columns in the consensus table.

**Graph-native cohort VCF projection (the non-canonical branch).** For each non-PvP01 reference, `vg deconstruct -P <ref>` on the GBZ-converted graph produces a pangenome VCF (8-sample) in that reference's coordinates. We then translate MalariaGEN's 1,895-sample variant sites with `odgi position` against the same graph and emit a cohort VCF per target. This is the "graph-native" branch of the 3-way (A1 wfmash / A2 KegAlign / B graph) projection comparison. **Path B lost to Path A2 on the drug-resistance QC** (median score 0.45 vs A2's 0.72), so A2 is canonical for genotype work — but Path B is preserved as a cross-validation comparator. The canonical (A2) path is documented in [MALARIAGEN_VCF_PROJECTION.md](MALARIAGEN_VCF_PROJECTION.md).

The graph is also the underlying structure for the chr1 subway-tube and chr1-zoom visualizations (`scripts/chr1_subway_tube.py`, `scripts/chr1_subway_zoom.py`) and for the panagram k-mer browser of the 2-way PvP01+PAM graph (`panagram_2way/`).

## Re-running on a different species

Same parameters work for any apicomplexan-scale haploid panel (5–15 assemblies, ~25 Mb each). The Pv4-pangenome repo has the recipe at `v3/pipeline/LOCAL.md` Section 5 (Phase D). Quick walk-through:

1. Drop N assembly FASTAs into `inputs/assemblies/{strain}.fa`. They must be one contig per chromosome, ideally with hard-masking dropped (soft-masking is fine — PGGB handles it).
2. Edit one config file `pipeline/species.conf`: set `STRAINS`, `REF_STRAIN`, and `N_CORES`. That's it.
3. Run `bash pipeline/02_mask.sh` (longdust ∪ sdust → soft-mask the FASTAs), then `bash pipeline/05_pggb.sh`. The latter does PanSN renaming, concatenation, and the pggb call. Expect ~3–6 hours on a 32-core box for 8 × 25 Mb genomes.

For different scales, two parameters merit attention:

- For closely related strains (intra-species like Pv), `-p 90 -s 5000` is right.
- For inter-species panels (e.g. *Plasmodium* clade with *P. cynomolgi*, *P. knowlesi* in the mix), drop to `-p 80 -s 5000` to keep mappings across the 10–15 % divergence boundary.
- For larger genomes (>100 Mb), bump segment-length to `-s 10000` to keep wfmash from exploding.

The pipeline is container-first — every tool runs via the wrapper `pipeline/lib/run_in_container.sh`, pinned to specific bioconda or `ghcr.io/pangenome/pggb` digests. No conda envs needed on the host.

The *P. knowlesi* scaffold at `Pk/v1/pipeline/` is a working copy of the recipe parameterized for that species (7 strains, H as reference). The user fills in `species.conf` accession list, runs `bash pipeline/setup/fetch_assemblies.sh`, then `bash pipeline/run_all.sh`.

## Galaxy tool-wrap list

For surfacing the build as a Galaxy workflow, these tools need wrappers (or version-bumped wrappers):

| Tool                                                | Galaxy state             | Note                                                           |
| --------------------------------------------------- | ------------------------ | -------------------------------------------------------------- |
| **pggb**                                            | ⚠ exists, refresh needed | The IUC wrapper lags the upstream image; pin a specific digest |
| **wfmash**                                          | ✅ exists                 |                                                                |
| **seqwish**                                         | ✅ exists                 |                                                                |
| **smoothxg**                                        | ✅ exists                 |                                                                |
| **odgi** (paths, pav, position, viz, layout, stats) | ✅ exists                 | Verify subcommand coverage matches the recipe                  |
| **gfaffix**                                         | ⚠ needs wrap             | Single Rust binary; trivial wrapper                            |
| **vg deconstruct**                                  | ⚠ needs wrap             | vg is wrapped for some subcommands; deconstruct may be missing |
| **vg convert** (GFA → GBZ)                          | ⚠ needs wrap             | Same situation                                                 |
| **PanSN rename**                                    | ⚠ needs wrap             | 16-line helper, custom tool                                    |
| **multi-FASTA concat + bgzip**                      | ✅ exists                 | Use existing fasta-merge + bgzip wrappers chained              |

The pggb step is the only heavy compute (~3–6 hours per build); everything else is sub-hour. Galaxy collection iteration over `{sample}` → per-sample PanSN-renamed FASTA → concat → pggb is the natural pipeline pattern.

Once the build wrapper is in place, the downstream consumers — consensus orthology ([ORTHOLOGY.md](ORTHOLOGY.md)), graph PAV cross-validation (described above in "How the graph fed downstream analyses"), and graph-native cohort VCF projection (the non-canonical comparator to [MALARIAGEN_VCF_PROJECTION.md](MALARIAGEN_VCF_PROJECTION.md)) — all become separate workflows that take the pggb `.og` as input. These three are independently wrap-able.
