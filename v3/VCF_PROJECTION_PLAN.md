# VCF projection plan — MalariaGEN PvP01 cohort onto the other 7 assemblies

## Context

MalariaGEN's Pv4 callset (1,895 samples, ~25 GB of per-chromosome VCFs on
`/media/anton/scratch/malariagen_pv4/`) is in PvP01 coordinates. To use these
calls for any analysis anchored on Sal-I, PvW1, PAM, MHC087, PvC01, PvT01, or
PvSY56, we need to project them into each target reference's coordinate system.

This plan compares **three projection methods** so we can pick the most reliable
one (or annotate sites with method-of-origin):

- **A1**: pairwise wfmash → chain → `bcftools +liftover`
- **A2**: pairwise lastz  → chain → `bcftools +liftover`
- **B**:  graph-native via `odgi position` + `vg deconstruct` against the v2 PGGB graph

Soft-masking (rather than the hard-masking v2 used) is requested so the chains
include indel-rich repeat boundaries where short reads still produce useful
SNP calls.

## Inputs

| input | path |
|---|---|
| MalariaGEN per-chromosome VCFs (16 files) | `/media/anton/scratch/malariagen_pv4/Pv4_PvP01_{01..14,API,MIT}_v1.vcf.gz` |
| Raw 8-assembly FASTAs | `/media/anton/data/sandbox/Pv4/v2/genomes/raw/*.fa` |
| Existing union mask BEDs (longdust + sdust) | `/media/anton/data/sandbox/Pv4/v2/genomes/mask_bed/{ACC}.union.bed` |
| v2 PGGB graph (OG + GFA + GBZ) | `/media/anton/data/sandbox/Pv4/v2/pggb_out/...smooth.fix.{gfa,og}`, `/media/anton/data/sandbox/Pv4/v2/vg_idx/pv.gbz` |

## Setup

1. **Soft-mask all 8 assemblies** (reuse v2 union BEDs):
   ```
   for acc in GCA_000002415.2 GCA_900093555.2 GCA_900093545.1 GCA_900093535.1 \
              GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1; do
     bedtools maskfasta \
       -fi /media/anton/data/sandbox/Pv4/v2/genomes/raw/${acc}.fa \
       -bed /media/anton/data/sandbox/Pv4/v2/genomes/mask_bed/${acc}.union.bed \
       -fo /media/anton/data/sandbox/Pv4/v3/genomes/softmasked/${acc}.fa \
       -soft
     samtools faidx /media/anton/data/sandbox/Pv4/v3/genomes/softmasked/${acc}.fa
   done
   ```

2. **New conda env with bcftools ≥ 1.20 + lastz**:
   ```
   conda create -n bcfnew -c bioconda -c conda-forge "bcftools>=1.20" htslib lastz -y
   ```

3. **Rename MalariaGEN VCFs' CHROM to match `PvP01_NN_v1` if needed**
   (already correct; only sanity-check).

## Path A1 — wfmash pairwise → chain → bcftools +liftover

For each of the 7 non-PvP01 references:

```
TARGET=GCA_900093555.2  # PvP01
OTHER=GCA_000002415.2   # e.g., Sal-I; loop over the other 6 too

# 1. Pairwise alignment with masking + best-hit-per-segment to suppress runaway
wfmash -s 5000 -p 90 -n 1 -X -Y'#' -t 18 \
       v3/genomes/softmasked/${TARGET}.fa \
       v3/genomes/softmasked/${OTHER}.fa \
       > v3/projection/A1_wfmash/PvP01_vs_${OTHER}.paf 2> v3/logs/A1_${OTHER}.wfmash.log

# 2. Filter PAF: drop ultra-short blocks and runaway aspect ratios
awk 'BEGIN{FS=OFS="\t"} $11>=1000 && ($11/$2)<=5 && ($2/$11)<=5' \
    v3/projection/A1_wfmash/PvP01_vs_${OTHER}.paf \
    > v3/projection/A1_wfmash/PvP01_vs_${OTHER}.filtered.paf

# 3. PAF -> chain
k8 paftools.js view -f chain \
   v3/projection/A1_wfmash/PvP01_vs_${OTHER}.filtered.paf \
   > v3/projection/A1_wfmash/PvP01_to_${OTHER}.chain

# 4. Liftover each chromosome VCF then concat
for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 API MIT; do
  bcftools +liftover --no-version -Oz -- \
    -s v3/genomes/softmasked/${TARGET}.fa \
    -f v3/genomes/softmasked/${OTHER}.fa \
    -c v3/projection/A1_wfmash/PvP01_to_${OTHER}.chain \
    /media/anton/scratch/malariagen_pv4/Pv4_PvP01_${chr}_v1.vcf.gz \
    > v3/projection/A1_wfmash/Pv4_${chr}_on_${OTHER}.vcf.gz
done
bcftools concat v3/projection/A1_wfmash/Pv4_*_on_${OTHER}.vcf.gz -Oz \
    -o v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz
bcftools index v3/projection/A1_wfmash/Pv4_cohort_on_${OTHER}.vcf.gz
```

**Why these wfmash flags**:
- `-s 5000 -p 90`: same segment-length and identity threshold as the v2 PGGB build → chains are consistent with what the graph captured
- `-n 1`: only the best mapping per segment, so a single locus on PvP01 maps to at most one locus per other assembly → suppresses runaway many-to-many alignments in tandem-repeat regions
- `-Y'#'`: skip self-mapping when query and target share a PanSN-style prefix (defensive)
- Soft-masked input is silently handled by wfmash (lowercase ≠ N; wfmash treats both as DNA but downstream MashMap-style sketching naturally avoids masked regions in seeding)

## Path A2 — lastz pairwise → chain → bcftools +liftover

```
TARGET=v3/genomes/softmasked/GCA_900093555.2.fa
OTHER=v3/genomes/softmasked/GCA_000002415.2.fa

# 1. Align with masking awareness + runaway-protection flags
lastz ${TARGET}[multiple] ${OTHER}[multiple] \
      --masking=50 \
      --hspthresh=4500 --gappedthresh=6000 \
      --inner=2000 --ydrop=15000 \
      --format=axt --chain \
      --output=v3/projection/A2_lastz/PvP01_vs_${OTHER}.axt 2> v3/logs/A2_${OTHER}.lastz.log

# 2. axt -> chain via UCSC tools (axtChain or axtToChain)
samtools faidx ${TARGET}; cut -f1,2 ${TARGET}.fai > PvP01.sizes
samtools faidx ${OTHER}; cut -f1,2 ${OTHER}.fai > ${OTHER}.sizes
axtToChain v3/projection/A2_lastz/PvP01_vs_${OTHER}.axt \
           PvP01.sizes ${OTHER}.sizes \
           v3/projection/A2_lastz/PvP01_to_${OTHER}.unsorted.chain
chainSort v3/projection/A2_lastz/PvP01_to_${OTHER}.unsorted.chain \
          > v3/projection/A2_lastz/PvP01_to_${OTHER}.chain

# 3. bcftools +liftover (identical to A1)
```

**Why these lastz flags**:
- `--masking=50`: skip seeds where ≥50% of the bases are soft-masked → respects our lowercase masking
- `--hspthresh=4500 --gappedthresh=6000`: standard UCSC pairwise-alignment thresholds; **caps runaway alignment growth in low-complexity regions** by raising the score floor
- `--inner=2000 --ydrop=15000`: bounded gap extension; another layer of runaway-protection
- The user explicitly asked to watch for runaway alignments; the combination of `--masking` + the threshold flags is the standard mechanism in the UCSC pairwise pipeline.

`axtToChain` + `chainSort` come from the UCSC kentUtils package (install via
`conda install -n bcfnew -c bioconda ucsc-axttochain ucsc-chainsort`).

## Path B — Graph-native via odgi + vg deconstruct

Reuses the existing v2 PGGB graph (no realignment needed):

```
# 1. For each non-PvP01 reference path in the v2 graph, produce a per-ref deconstruct VCF
for ref in GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 \
           GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1; do
  vg deconstruct -P ${ref} -t 18 \
     /media/anton/data/sandbox/Pv4/v2/vg_idx/pv.gbz \
     > v3/projection/B_graph/pangenome_on_${ref}.vcf 2> v3/logs/B_deconstruct_${ref}.log
done

# 2. For each MalariaGEN PASS variant, look up its position on each other ref via odgi position
bcftools query -f '%CHROM\t%POS0\t%POS\t%CHROM:%POS\n' \
    /media/anton/scratch/malariagen_pv4/Pv4_PvP01_*_v1.vcf.gz \
    | sort -k1,1 -k2,2n > v3/projection/B_graph/mg_sites.bed

for ref in GCA_000002415.2 GCA_900093545.1 GCA_900093535.1 \
           GCA_914969965.1 GCA_949152365.1 GCA_003402215.1 GCA_040114635.1; do
  odgi position \
    -i /media/anton/data/sandbox/Pv4/v2/pggb_out/pggb_input.fa.gz.f705205.c28ecf8.dcad0e6.smooth.fix.og \
    -b v3/projection/B_graph/mg_sites.bed \
    -r ${ref} -t 18 \
    > v3/projection/B_graph/mg_sites_on_${ref}.tsv 2> v3/logs/B_odgi_${ref}.log
done
```

Then join the odgi-position output (PvP01-pos → other-ref-pos) against
the original MalariaGEN VCF to produce a Path-B liftover VCF.

**Strengths/weaknesses**:
- Handles complex rearrangements that chain files lose
- Limited to regions where the target path exists in the graph (we did mask
  ~17% of each genome in v2; those bases are absent from the v3 soft-masked
  graph paths). For variants in those regions, Path B silently drops them —
  expect lower yield than A1/A2 in subtelomeres.

## Comparison

Per target reference:

| metric | A1 | A2 | B |
|---|---|---|---|
| input variants (PvP01) | N | N | N |
| variants successfully lifted | n_A1 | n_A2 | n_B |
| dropped — chain gap / no graph path | ... | ... | ... |
| dropped — REF mismatch | ... | ... | ... |
| dropped — multi-allelic ambiguity | ... | ... | ... |

Then pairwise concordance:
```
bcftools isec -p A1_vs_A2 A1.vcf.gz A2.vcf.gz       # 11 = in both, 10 = A1-only, 01 = A2-only
bcftools isec -p A1_vs_B  A1.vcf.gz B.vcf.gz
bcftools isec -p A2_vs_B  A2.vcf.gz B.vcf.gz
bcftools isec -p all3 -n=3 A1.vcf.gz A2.vcf.gz B.vcf.gz   # consensus
```

## Verification

The known drug-resistance codons (from `Pv4_drug_resistance_marker_genotypes.txt`)
must land on the orthologous gene's CDS in every target reference, irrespective
of method:

- *dhfr-ts* codons 57, 58, 61, 111, 117, 173 (PvP01 chr 5: ~1075584–1079966)
- *dhps* codon 383 (PvP01 chr 14: 1269756–1272657)
- *mdr1* CN locus (PvP01 chr 10: 476677–483934)

Cross-validated against `ORTHOLOGY_PLAN.md` Phase 5 (the projection-vs-ortholog
QC bridge).

## Deliverables

- 7 cohort-level VCFs per method × 3 methods = 21 VCFs in `v3/projection/{A1,A2,B}/`
- One `comparison_summary.tsv` reporting per-pair, per-method success rates and
  inter-method concordance
- One markdown report in `v3/projection/REPORT.md` summarizing findings
