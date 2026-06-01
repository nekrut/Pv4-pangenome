# Clean-pipeline scaffold — bugs found by the first end-to-end test run

The clean LOCAL pipeline (11 phase scripts + `run_all.sh` + `lib/run_in_container.sh` + helpers) exists as the **P. knowlesi scaffold** (`/media/anton/data/sandbox/Pk/v1/pipeline/`), generalized to be species-agnostic via `species.conf`. It had **never been run end-to-end**. Running it against the Pv4 test panel (5 strains × 3 chromosomes; see `test_data/`) surfaced the bugs below — this is exactly what the test data is for. Fixes were applied in the ephemeral test workspace (`/media/anton/data/sandbox/Pv4test`); they must be folded into the canonical scaffold before it ships as the clean Pv4 LOCAL pipeline.

Host context: docker present, RTX A5000 GPU present, bio-tools available only via containers.

## Phase A — Inventory (mash + BUSCO) — GREEN

1. **`lib/*.sh` not executable** → `chmod +x pipeline/*.sh pipeline/lib/*.sh`. (git/cp can drop the bit.)
2. **Stale `mash` image tag** `quay.io/biocontainers/mash:2.3--he348c14_4` (404) → `staphb/mash:2.3` (or `quay…/mash:2.3--hb105d93_10`).
3. **Stale `busco` image tag** `ezlabgva/busco:v5.7.1_cv1` → `ezlabgva/busco:v5.8.0_cv1`.

## Phase B — Soft-masking — GREEN

4. **`WORK` not exported** — phase scripts source `species.conf` (sets `WORK`) but don't export it, so the `run_in_container.sh` child dies on `set -u` with `WORK: unbound`. Fix: `export WORK SCRATCH USE_GPU GPU SPECIES N_CORES` in `species.conf`. (Would have broken every phase.)
5. **`longdust`/`sdust` collapsed to one image + stale tag** — they're separate biocontainers, and longdust is now v1.4: `longdust → quay…/longdust:1.4--h577a1d6_0`, `sdust → quay…/sdust:0.1--h077b44d_2`.
6. **`longdust -t $N_CORES` misuse** — longdust's `-t` is a *score threshold* (FLOAT, default 0.6), not threads; `-t 16` suppressed all output (empty masks). longdust has no threads flag. Fix: drop `-t`.
7. **`docker run` missing `-i`** — every containerized tool that reads piped stdin (`bedtools sort -i -`, the whole chain pipeline) got no input and the upstream producer died with SIGPIPE (141). Fix: `docker run --rm -i …` in `run_in_container.sh`. **Most impactful fix.**
8. **`samtools quickcheck` used to validate FASTA** — quickcheck is for BAM/CRAM/VCF and rejects plain FASTA. Fix: validate via the faidx index instead (`[[ -s ${FA_OUT}.fai ]]`; faidx already fails on malformed FASTA).

## Phase C.1–3 — Pairwise align + chains

9. **No `kegalign`/`kegalign-gpu` biocontainer exists** on quay → use the script's lastz CPU fallback (`USE_GPU=0`). Fine for small test genomes; GPU KegAlign untested for lack of an image.
10. **Stale lastz tag** in the fallback branch `lastz:1.04.22--h0c08fa6_1` → `lastz:1.04.52--h7b50bb2_1`.
11. **AXT validator too strict** — `awk 'NR<10 && /^[0-9]+ /'` only checks the first 9 lines, but lastz writes a 13-line `#` comment header, so the first block is line 14. Valid AXTs failed. Fix: `grep -qE '^[0-9]+ ' "$AXT_OUT"`.
12. **`ucsc-kent-tools` meta-image won't pull** — every tag 401s on quay. Map each kent tool to its **individual** biocontainer instead: `axtChain→ucsc-axtchain:482--h0b57e2e_2`, `chainSort→ucsc-chainsort:482--h0b57e2e_0`, `chainPreNet→ucsc-chainprenet`, `chainNet→ucsc-chainnet`, `netChainSubset→ucsc-netchainsubset`, `chainStitchId→ucsc-chainstitchid`, `chainSwap→ucsc-chainswap`, `axtToMaf→ucsc-axttomaf` (all `482--h0b57e2e_0`).
13. **`axtChain` given FASTA without `-faT -faQ`** — defaults to expecting 2bit/nib for the tNibDir/qNibDir args (`... is not a 2bit file or a directory`). Fix: `axtChain -linearGap=loose -faT -faQ …`.
15. **Directional AXT misused for B→A** — the scaffold reran `axtChain` on the *same* (A→B) AXT with the FASTA args swapped to build the B→A chain. An AXT is directional (target=A, query=B), so axtChain can't find B's contigs in A.fa (`can not find sequence name … from fasta file`). Fix: build A→B once, then `chainSwap` it to get B→A (no re-alignment). (Its own comment "orientation is implicit" was the wrong assumption.)
14. **Chain pipeline (3.2) mis-wired** — `netChainSubset` was fed `chainNet`'s stdout (empty; chainNet writes to *files*), so the net referenced chain IDs absent from the empty chain input (`hashMustFindVal: '8' not found`). `netChainSubset` must take the **pre-net** chain (the same chains given to `chainNet`). Fix: save the chainPreNet output to a temp file, run chainNet on it, then `netChainSubset net.tmp prenet.tmp | chainStitchId`. (The rbest block 3.3 was already wired correctly.)

## Phases D (PGGB) + C.4 (annotation)

16. **`python3` mapped to an image with no python3** — the wrapper routed `python3` to the bcftools biocontainer "ships python3 + pysam + pandas + numpy" (it does not; `python3` isn't on its PATH), so every helper-script step died with `exec: "python3": executable file not found`. This blocked **four phases** (D PanSN-rename, C.4 triage/merge, E consensus). The only third-party dep across all helpers is `pyfaidx` (lazy-imported in `phase_c2_triage.py`; everything else is stdlib) — so map `python3 → quay.io/biocontainers/pyfaidx:0.9.0.4--pyhdfd78af_0` (full python + pyfaidx). (Liftoff itself worked: PvP01→PvW1 lifted 963 genes.)
17. **Stale htslib (bgzip/tabix) tag** `quay.io/biocontainers/htslib:1.20--h5efdd21_0` (404). The samtools image already ships `bgzip` + `tabix` — map `bgzip|tabix → quay.io/biocontainers/samtools:1.20--h50ea8bc_0`.

18. **Stale PGGB image tag** `ghcr.io/pangenome/pggb:202412130311080800a17` (404) → `ghcr.io/pangenome/pggb:latest` (pin a real digest before shipping).
19. **TOGA run unconditionally even with an empty rescue queue** — `04_annotate_project.sh` invoked `toga` whenever its output didn't exist yet, without checking `needs_cesar2.bed`. For closely-related strains the triage queue is empty (all genes lift clean), so TOGA had nothing to do — yet it ran and hit the unavailable `ghcr.io/hillerlab/toga:latest` image (`error from registry: denied`). Fix: guard with `elif [[ ! -s "$TRIAGE_DIR/needs_cesar2.bed" ]]; then skip`. **Separate open issue:** the TOGA image is not publicly pullable, so the CESAR rescue path is untested — needs a working TOGA container (or a documented build) before the rescue tail can run for divergent panels.

## Recurring pattern

The scaffold's ~18 pinned container tags were written speculatively and many are stale/removed/unpullable. Before shipping: **verify every tag pulls** (`docker pull`), prefer individual biocontainers over meta-packages, and pin to digests. The Pv4 docs already switched Phase A mash→sourmash; the scaffold still uses mash (works, but inconsistent with the Pv4 LOCAL.md).

## Status at time of writing

A ✅, B ✅, C.1–3 in progress (chain pipeline). Phases C.4 (TOGA — heaviest, anchor BED12/isoforms auto-built), D (PGGB), E, F, G, H, I, J not yet exercised.
