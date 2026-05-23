#!/usr/bin/env python3
"""Generate the Pvs230 domain-mapping figure for the appendix."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = Path("/media/anton/data/sandbox/Pv4/v3")
PRIO = ROOT / "work/06_msa/core_v3_hyphy/priority/PVP01_0415800"
OUT = ROOT / "writeup/hyphy_plots/figA1_pvs230_domain_map.png"

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

# Parse alignment, build ref-codon -> protein-residue map
seqs = {}
cur = None
for line in open(PRIO / "aln.fa"):
    line = line.rstrip()
    if line.startswith(">"):
        cur = line[1:].split()[0]
        seqs[cur] = []
    else:
        seqs[cur].append(line)
for k in seqs:
    seqs[k] = "".join(seqs[k])

ref = "PvP01_REF"
ref_seq = seqs[ref]
codon_to_refpos = {}
ref_pos = 0
for i in range(len(ref_seq) // 3):
    codon = ref_seq[3 * i:3 * i + 3]
    if codon == "---":
        codon_to_refpos[i + 1] = None
    else:
        ref_pos += 1
        codon_to_refpos[i + 1] = ref_pos
prot_len = ref_pos

# Parse FEL and MEME
def parse_fel():
    d = json.load(open(PRIO / "fel.json"))
    h = [x[0] for x in d["MLE"]["headers"]]
    rows = d["MLE"]["content"]["0"]
    ia, ib, ip = h.index("alpha"), h.index("beta"), h.index("p-value")
    out = []
    for ci, r in enumerate(rows, 1):
        a, b, pv = float(r[ia]), float(r[ib]), float(r[ip])
        out.append((ci, a, b, pv))
    return out


def parse_meme():
    d = json.load(open(PRIO / "meme.json"))
    h = [x[0] for x in d["MLE"]["headers"]]
    rows = d["MLE"]["content"]["0"]
    ip = h.index("p-value")
    return [(ci, float(r[ip])) for ci, r in enumerate(rows, 1)]


fel = parse_fel()
meme = parse_meme()

# Build x/y for plot — use PvP01 protein residue on x; codons that are gap-in-ref get dropped
def to_resseries(records, value_getter):
    xs, ys = [], []
    for rec in records:
        ci = rec[0]
        rp = codon_to_refpos[ci]
        if rp is None:
            continue
        xs.append(rp)
        ys.append(value_getter(rec))
    return np.array(xs), np.array(ys)


# FEL: -log10(p) signed by sign(beta-alpha)
def fel_signed_logp(rec):
    _, a, b, p = rec
    p = max(p, 1e-10)
    s = -np.log10(p)
    return s if b > a else -s


fel_x, fel_y = to_resseries(fel, fel_signed_logp)
meme_x, meme_y = to_resseries(meme, lambda r: -np.log10(max(r[1], 1e-10)))

# Figure: two-panel — top FEL, bottom MEME, with IVP/ICP/C-term shading
fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=True)

IVP_END = 269
ICP_END = 954

def shade_domains(ax):
    ax.axvspan(1, IVP_END, color="#ffd6a8", alpha=0.55, lw=0, label="IVP (1-269)")
    ax.axvspan(IVP_END + 1, ICP_END, color="#cfe2f3", alpha=0.55, lw=0, label="ICP (270-954)")
    ax.axvspan(ICP_END + 1, prot_len, color="#e8e8e8", alpha=0.55, lw=0, label="C-terminal (>954)")

# Top — FEL
ax = axes[0]
shade_domains(ax)
ax.scatter(fel_x, fel_y, s=3, c="grey", alpha=0.35, lw=0)
# Highlight significant
sig_pos_mask = np.array([
    (codon_to_refpos[r[0]] is not None) and (r[3] < 0.05) and (r[2] > r[1])
    for r in fel
])
sig_neg_mask = np.array([
    (codon_to_refpos[r[0]] is not None) and (r[3] < 0.05) and (r[2] < r[1])
    for r in fel
])
fel_x_full = np.array([codon_to_refpos[r[0]] for r in fel if codon_to_refpos[r[0]] is not None])
fel_y_full = np.array([fel_signed_logp(r) for r in fel if codon_to_refpos[r[0]] is not None])
# Recompute masks aligned to filtered arrays
keep_idx = [i for i, r in enumerate(fel) if codon_to_refpos[r[0]] is not None]
fel_kept = [fel[i] for i in keep_idx]
mask_pos = np.array([(r[3] < 0.05) and (r[2] > r[1]) for r in fel_kept])
mask_neg = np.array([(r[3] < 0.05) and (r[2] < r[1]) for r in fel_kept])
ax.scatter(fel_x_full[mask_pos], fel_y_full[mask_pos], s=24, c="firebrick", label=f"FEL+ p<0.05 (n={mask_pos.sum()})", zorder=5)
ax.scatter(fel_x_full[mask_neg], fel_y_full[mask_neg], s=18, c="#1f77b4", label=f"FEL- p<0.05 (n={mask_neg.sum()})", zorder=4)
ax.axhline(0, color="k", lw=0.6)
ax.axhline(-np.log10(0.05), color="firebrick", ls=":", lw=0.6)
ax.axhline(np.log10(0.05), color="#1f77b4", ls=":", lw=0.6)
ax.set_ylabel("FEL signed -log10(p)\n(positive = beta>alpha)")
ax.set_title("Pvs230 (PVP01_0415800) per-codon selection signal mapped to PvP01 protein coordinates")
# annotate positive hits
for x, y in zip(fel_x_full[mask_pos], fel_y_full[mask_pos]):
    ax.annotate(f"res {int(x)}", (x, y), fontsize=7, xytext=(4, 4), textcoords="offset points")
ax.legend(loc="upper right", frameon=False, fontsize=7, ncol=2)

# Bottom — MEME
ax = axes[1]
shade_domains(ax)
ax.scatter(meme_x, meme_y, s=3, c="grey", alpha=0.35, lw=0)
meme_kept = [(codon_to_refpos[r[0]], -np.log10(max(r[1], 1e-10)), r[1]) for r in meme if codon_to_refpos[r[0]] is not None]
mx = np.array([t[0] for t in meme_kept])
my = np.array([t[1] for t in meme_kept])
mp = np.array([t[2] for t in meme_kept])
m_sig = mp < 0.05
ax.scatter(mx[m_sig], my[m_sig], s=24, c="darkgreen", label=f"MEME p<0.05 (n={m_sig.sum()})", zorder=5)
ax.axhline(-np.log10(0.05), color="darkgreen", ls=":", lw=0.6)
ax.set_ylabel("MEME -log10(p)")
ax.set_xlabel("PvP01 protein residue (1-based)")
# annotate hits
for x, y in zip(mx[m_sig], my[m_sig]):
    ax.annotate(f"res {int(x)}", (x, y), fontsize=7, xytext=(4, 4), textcoords="offset points")
ax.legend(loc="upper right", frameon=False, fontsize=7)

# Add a domain legend below
patches = [
    mpatches.Patch(color="#ffd6a8", alpha=0.55, label="IVP — positively selected in Feng 2022"),
    mpatches.Patch(color="#cfe2f3", alpha=0.55, label="ICP — purifying in Feng 2022"),
    mpatches.Patch(color="#e8e8e8", alpha=0.55, label="C-terminal (outside Feng surveyed region)"),
]
axes[1].legend(handles=patches + [mpatches.Patch(color="darkgreen", label=f"MEME p<0.05 (n={m_sig.sum()})")],
               loc="upper right", frameon=False, fontsize=7, ncol=2)

fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print(f"Wrote {OUT}")
print(f"Protein length: {prot_len}")
print(f"FEL+ hits: {[(int(fel_x_full[i]), float(fel_y_full[i])) for i in range(len(fel_x_full)) if mask_pos[i]]}")
print(f"FEL- hits: {[int(fel_x_full[i]) for i in range(len(fel_x_full)) if mask_neg[i]]}")
print(f"MEME hits (residue, p): {[(int(mx[i]), float(mp[i])) for i in range(len(mx)) if m_sig[i]]}")
