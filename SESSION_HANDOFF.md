# Session handoff — 2026-05-26

Pickup notes for resuming the Pv4 v3 *P. vivax* pangenome publishing work on a different machine. Everything that matters is in the GitHub repo or referenced by URL. Local file paths in this doc are the pre-handoff machine's paths (`/media/anton/data/sandbox/Pv4/…`); on the new machine, the equivalent will be the path you clone the repo to.

## Where we are

The Pv4 v3 8-strain *P. vivax* pangenome analysis is complete and being prepared for BRC-analytics deployment. The data, docs, blog, and UCSC hub are all built. Active threads are around external review (BRC-analytics issue, Galaxy skills PR) and waiting on BRC-side implementation.

## What's live online

### GitHub repositories

| Repo | Branch | Purpose |
|---|---|---|
| [`nekrut/Pv4-pangenome`](https://github.com/nekrut/Pv4-pangenome) | `main` | Canonical analysis repo — code, small data, docs, pipeline plans, blog, UCSC hub config |
| [`galaxyproject/brc-analytics`](https://github.com/galaxyproject/brc-analytics) | — | BRC-analytics catalog (target for deployment PRs) |
| [`galaxyproject/galaxy-skills`](https://github.com/galaxyproject/galaxy-skills) | `main` | Galaxy-related Claude/agent skills (trackhubs PR open here) |

### Open issues + PRs

| Item | URL | State |
|---|---|---|
| BRC-analytics issue #1279 — "Add *P. vivax* pangenome bundle" | <https://github.com/galaxyproject/brc-analytics/issues/1279> | open; rewritten; carries the full plan + UCSC file listing + corrections |
| galaxy-skills PR #18 — trackhubs skill | <https://github.com/galaxyproject/galaxy-skills/pull/18> | open, awaiting review |

### Dropbox folder

Large outputs (multiz MAFs, cohort VCFs, bigMaf, bigChain bundle, the full UCSC hub) live on Dropbox in folder `Pv4_v3/`. The folder-level share URL is in `v3/writeup/LARGE_FILES_DROPBOX.md` and `v3/README.md`. Per-file share URLs (~270 of them) live in `/tmp/dropbox_links.tsv` and `/tmp/dropbox_links_v2.tsv` on the source machine — **these need to be regenerated on the new machine** with `rclone link` against the same `dropbox:Pv4_v3/` remote (see "Resuming on a new machine" below).

## What was done this session (recent → older)

1. **UCSC hub configuration fix** (commit `bf02523f`)
   - Converted all 56 chain files from `.chain.gz` to `bigChain.bb` + `bigChain.link.bb` (UCSC track hubs need the indexed binary, not gzipped text)
   - Split the previous mixed-type `brc_pangenome_align` composite into a standalone bigMaf + separate `brc_pangenome_chains` composite of `type bigChain` (UCSC requires composites to be single-type)
   - Set `defaultPos` per hub to the dhps locus (`PVP01_1429500` chr14 ortholog) — every hub now jumps to the same biology-relevant region on first load
   - Added `.2bit` symlinks into `projection/A2_kegalign/2bit/` for every hub
   - Wrote `v3/tools/chain_to_bigChain.py` + `bigChain.as` + `bigLink.as` (kentUtils-free conversion path)
   - Updated LOCAL.md + GALAXY.md Phase K sections; added `BRC_DEPLOYMENT.md` corrections
   - Posted [correction comment on brc-analytics#1279](https://github.com/galaxyproject/brc-analytics/issues/1279#issuecomment-4548150201)
   - `v3/writeup/UCSC_CORRECTIONS.md` mirrors the comment

2. **Galaxy-skills `trackhubs` PR** (PR #18 above)
   - New skill family `trackhubs/` with SKILL.md + 4 references (composite rules, chain→bigChain, genomes.txt fields, hubCheck debugging) distilled from this exact review cycle

3. **BLOG1 fact-check + image fixes + PDF**
   - Per-panel image embeds (was repeating composite)
   - Updatability paragraph: Cactus uses `cactus-update-prepare` for incremental adds; PGGB has no incremental path (wfmash is decomposable in principle but seqwish/smoothxg always rebuild from full PAF)
   - Fixed: -n 8 flag (haplotype count, not multi-mapping), PvP01 sequencing (Illumina short-read, not long-read), 256→255 block count, gfaffix added to pipeline list
   - PDF built and reviewed

4. **BRC-analytics deployment doc set**
   - `v3/writeup/BRC_DEPLOYMENT.md` — 5-block sketch: data-model entries, 8 UCSC hubs, organism-page Pangenome section, Galaxy-backed services, 4-PR sequence
   - `v3/writeup/PANGENOME_MOCKUP.md` — ASCII mockup of the organism-page Pangenome section + orthogroup detail page
   - **Pattern named**: "Galaxy as query backend for data portals" — sourmash distance lookups and MMseqs2 protein homology searches are Galaxy workflows BRC runs in the background, not BRC-hosted indexes

5. **Per-output documentation set** — the project "memory" lives here:
   - `v3/writeup/PANGENOME.md` — PGGB graph build (wfmash → seqwish → smoothxg → gfaffix → odgi)
   - `v3/writeup/MULTIZ.md` — KegAlign GPU lastZ + UCSC chain pipeline → multiz progressive fold
   - `v3/writeup/MALARIAGEN_VCF_PROJECTION.md` — 1,895-sample MalariaGEN cohort lift to 7 non-PvP01 assemblies
   - `v3/writeup/ORTHOLOGY.md` — 3-stream union-find consensus orthology
   - `v3/writeup/MSA_HYPHY.md` — per-orthogroup codon MSAs + BUSTED + BUSTED-MH
   - `v3/writeup/MICROSYNTENY.md` — neighborhood conservation analysis
   - All cross-linked with section anchors; tables monospace-aligned; written in Anton's voice (not AI jargon)

## Where to pick up

No explicit pending task. Likely next moves, in rough order of urgency:

1. **Wait for response on brc-analytics#1279.** The issue + UCSC files correction comment is the current handoff to the BRC team. They need to review and then a BRC developer opens the 4-PR sequence sketched in `v3/writeup/BRC_DEPLOYMENT.md`.

2. **Wait for galaxy-skills PR #18 review.** Self-contained skill; no follow-up work expected unless reviewers request changes.

3. **Optional: build the `.bb.bai` mafIndex companions** for each of the 8 bigMaf files. Skipped in v1 — hubs work without them but are slower at whole-genome zoom. ~5 min per assembly. Tool: `mafIndex` from kentUtils (in `v3/tools/mafIndex`).

4. **Optional: lay out the P. knowlesi pipeline scaffold.** `Pk/` directory has a stub started; LOCAL.md / GALAXY.md were generalized to make a second-species build mostly parameter changes. Not yet exercised end-to-end.

## File locations cheat-sheet

```
Pv4-pangenome/
├── README.md                     ← top-level orientation, Dropbox folder link
├── SESSION_HANDOFF.md            ← this file
├── v2/pggb_out/                  ← v2 8-way PGGB graph (gzipped GFA + small index files only)
└── v3/
    ├── README.md                 ← analysis-level overview
    ├── pipeline/                 ← LOCAL.md (containers + bash), GALAXY.md (workflow port plan)
    ├── tools/                    ← kent utilities + custom scripts (chain_to_bigChain.py, bigChain.as, bigLink.as)
    ├── projection/               ← VCF projection paths (A1 wfmash, A2 KegAlign, B graph-native)
    │   └── A2_kegalign/2bit/     ← .2bit files for all 8 assemblies; bigChain build inputs
    ├── work/01_chains/           ← canonical .chain.gz (Git) + .sizes per assembly
    ├── work/07_multiz/           ← multiz output (intermediates gitignored; bigMaf on Dropbox)
    ├── ucsc_hub/                 ← deployable hub (small files in Git; bigMafs gitignored, on Dropbox)
    │   ├── hub.txt
    │   ├── genomes.txt
    │   └── GCA_*/                ← 8 hubs, each with trackDb.txt + 2bit + .bb files + chains/
    └── writeup/                  ← all docs (analysis memory, deployment plans, blog, correction notes)
```

## Resuming on a new machine

Minimum environment to be able to do anything useful:

```bash
# 1. Clone the repo
git clone git@github.com:nekrut/Pv4-pangenome.git
cd Pv4-pangenome

# 2. gh CLI auth (for issue/PR work)
gh auth login

# 3. (optional) rclone with Dropbox remote
# Only needed for re-uploading or regenerating per-file share URLs.
# The folder-level share URL in v3/README.md already gives you read access without rclone.
~/.local/bin/rclone config              # set up `dropbox:` remote
~/.local/bin/rclone ls dropbox:Pv4_v3/  # smoke test

# 4. (optional, only if rebuilding) container engine for the kent tools
# v3/tools/ has prebuilt linux binaries committed.
# For other platforms, use the quay.io/biocontainers/ucsc-kent-tools image referenced in
# v3/pipeline/LOCAL.md.

# 5. (optional, only if re-running phases) conda envs
# v3/pipeline/LOCAL.md lists each env's contents.
# Recreate with: conda create -n bcfnew -c bioconda -c conda-forge bcftools htslib lastz
#                conda create -n orth   -c bioconda mmseqs2 macse pal2nal orthofinder sonicparanoid
#                etc.
```

The Dropbox folder share URL (in `v3/README.md` and `v3/writeup/LARGE_FILES_DROPBOX.md`) gives anyone direct download access to all heavy outputs — no rclone needed for read-only consumption.

## What's NOT in this repo

- Multiz `.maf.bb` files (~4.6 GB total) — gitignored, on Dropbox
- Multiz `.maf.gz` raw alignments (~24 GB) — Dropbox only
- Cohort VCFs (`Pv4_cohort_on_*.vcf.gz`, ~115 GB) — Dropbox only
- Per-orthogroup HyPhy JSONs (~2 GB) — Dropbox only
- Galaxy histories — never staged; LOCAL.md is the canonical recipe
- The full UCSC hub including the 4.6 GB of bigMaf files — Dropbox `ucsc_hub/`

All of these can be regenerated from the recipe in `v3/pipeline/LOCAL.md` against the inputs in `v3/genomes/`.

## Useful URLs (one-stop)

- Repo: <https://github.com/nekrut/Pv4-pangenome>
- BRC issue: <https://github.com/galaxyproject/brc-analytics/issues/1279>
- BRC correction comment: <https://github.com/galaxyproject/brc-analytics/issues/1279#issuecomment-4548150201>
- Galaxy-skills PR: <https://github.com/galaxyproject/galaxy-skills/pull/18>
- Dropbox folder: see `v3/README.md` (folder share URL is the rclone link to `Pv4_v3/`)
