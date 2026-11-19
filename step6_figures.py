"""
Step 6: manuscript figures (PNG 300 dpi + vector PDF) from the pipeline outputs.
Run from the data folder after steps 1-5:  python step6_figures.py
Reads : out/coq_with_alphamissense.csv, out3/structure_features.csv, out5/coq4_functional_variants.csv
Writes: figures/Figure1 .. Figure5 (.png and .pdf)
Cross-validation AUROC values in Figure 4A are copied from out4/added_value_log.txt (step 4).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HI, LO = 0.564, 0.34
C_P, C_B, C_V, C_G = "#C0392B", "#2E6DA4", "#9AA0A6", "#E69F00"   # colorblind-safe-ish
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": "DejaVu Sans", "pdf.fonttype": 42})
os.makedirs("figures", exist_ok=True)
GENES = ["PDSS1","PDSS2","COQ2","COQ3","COQ4","COQ5","COQ6","COQ7","COQ8A","COQ8B","COQ9","COQ10A","COQ10B"]

def save(fig, name):
    fig.savefig(f"figures/{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"figures/{name}.pdf", bbox_inches="tight")
    plt.close(fig)

allv = pd.read_csv("out/coq_with_alphamissense.csv")
sf = pd.read_csv("out3/structure_features.csv")
lab = sf[sf.label.isin(["P", "B"])].copy()
lab["y"] = (lab.label == "P").astype(int)

# ---------- Figure 1: dataset ----------
funnel = [("ClinVar rows (all genes)", 9355006), ("Target genes", 8126), ("GRCh38", 4037),
          ("Single-nucleotide variants", 3667), ("True missense", 1718),
          ("P/LP, B/LB or VUS", 1563), ("Review stars >= 1 (P, B)", 1544)]
fig, ax = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw={"width_ratios": [1, 1.25]})
labels = [f for f, _ in funnel][1:]; vals = [v for _, v in funnel][1:]
ax[0].barh(range(len(vals))[::-1], vals, color="#5B7C99")
for i, v in zip(range(len(vals))[::-1], vals):
    ax[0].text(v + 80, i, f"{v:,}", va="center", fontsize=8)
ax[0].set_yticks(range(len(vals))[::-1]); ax[0].set_yticklabels(labels)
ax[0].set_xlabel("Variants (after each filter)"); ax[0].set_xlim(0, 9500)
ax[0].set_title("A  Filtering of ClinVar rows in the 13 genes\n(9,355,006 rows read)", loc="left", fontsize=9)
cnt = allv.groupby(["gene", "label"]).size().unstack(fill_value=0).reindex(GENES).fillna(0)
x = np.arange(len(GENES)); w = 0.4
ax[1].bar(x - w/2, cnt.get("P", 0), w, color=C_P, label="Pathogenic / likely pathogenic")
ax[1].bar(x + w/2, cnt.get("B", 0), w, color=C_B, label="Benign / likely benign")
ax[1].axhline(10, color="k", ls=":", lw=0.8); ax[1].text(len(GENES) - 0.5, 10.6, "10 variants", ha="right", fontsize=7)
ax[1].set_xticks(x); ax[1].set_xticklabels(GENES, rotation=45, ha="right")
ax[1].set_ylabel("Labeled variants"); ax[1].legend(frameon=False, fontsize=8)
ax[1].set_title("B  Labeled variants per gene (VUS not shown; 1,395 in total)", loc="left", fontsize=9)
fig.tight_layout(); save(fig, "Figure1")

# ---------- Figure 2: discrimination ----------
def roc_points(y, s):
    thr = np.unique(np.r_[s, 0, 1])[::-1]
    P, N = (y == 1).sum(), (y == 0).sum()
    return [((s >= t) & (y == 0)).sum() / N for t in thr], [((s >= t) & (y == 1)).sum() / P for t in thr]
y, s = lab.y.values, lab.am_pathogenicity.values
fpr, tpr = roc_points(y, s)
fig, ax = plt.subplots(1, 2, figsize=(9, 4))
ax[0].plot(fpr, tpr, color="k", lw=1.8); ax[0].plot([0, 1], [0, 1], ls=":", color="grey")
ax[0].text(0.55, 0.15, "AUROC 0.976\n(95% CI 0.953-0.992)\nn = 144 (67 P, 77 B)", fontsize=8)
ax[0].set_xlabel("False positive rate"); ax[0].set_ylabel("True positive rate")
ax[0].set_title("A  AlphaMissense, ClinVar-labeled variants", loc="left", fontsize=9)
rng = np.random.default_rng(1)
for lbl, col, pos in [(0, C_B, 0), (1, C_P, 1)]:
    v = lab[lab.y == lbl].am_pathogenicity.values
    ax[1].scatter(pos + rng.uniform(-0.18, 0.18, len(v)), v, s=14, color=col, alpha=0.75, linewidths=0)
    ax[1].hlines(np.median(v), pos - 0.3, pos + 0.3, color="k", lw=1.5)
ax[1].axhline(HI, color="grey", ls="--", lw=0.8); ax[1].axhline(LO, color="grey", ls="--", lw=0.8)
ax[1].text(1.45, HI + 0.01, "0.564", fontsize=7, ha="right"); ax[1].text(1.45, LO + 0.01, "0.34", fontsize=7, ha="right")
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["Benign / likely benign", "Pathogenic / likely pathogenic"])
ax[1].set_xlim(-0.5, 1.5); ax[1].set_ylabel("AlphaMissense score")
ax[1].set_title("B  Score by ClinVar class (bars: medians)", loc="left", fontsize=9)
fig.tight_layout(); save(fig, "Figure2")

# ---------- Figure 3: by gene, discordant labeled ----------
fig, ax = plt.subplots(figsize=(10, 4.8))
disc_idx = lab[((lab.y == 0) & (lab.am_pathogenicity > HI)) | ((lab.y == 1) & (lab.am_pathogenicity < LO))].index
gi = {g: i for i, g in enumerate(GENES)}
rng = np.random.default_rng(2)
lab["xpos"] = [gi[g] + (0.17 if yy == 1 else -0.17) + rng.uniform(-0.1, 0.1) for g, yy in zip(lab.gene, lab.y)]
for lbl, col in [(0, C_B), (1, C_P)]:
    d = lab[lab.y == lbl]
    ax.scatter(d.xpos, d.am_pathogenicity, s=16, color=col, alpha=0.8, linewidths=0,
               label="Benign / likely benign" if lbl == 0 else "Pathogenic / likely pathogenic")
ax.axhline(HI, color="grey", ls="--", lw=0.8); ax.axhline(LO, color="grey", ls="--", lw=0.8)
ax.text(len(GENES) - 0.4, HI + 0.012, "0.564", fontsize=7, ha="right", color="grey")
ax.text(len(GENES) - 0.4, LO + 0.012, "0.34", fontsize=7, ha="right", color="grey")
placed = {}
for i in disc_idx:
    r = lab.loc[i]
    k = (r.gene, r.y); n = placed.get(k, 0); placed[k] = n + 1
    ax.scatter([r.xpos], [r.am_pathogenicity], s=70, facecolors="none", edgecolors="k", linewidths=1.1)
    ax.annotate(r.protein_variant, (r.xpos, r.am_pathogenicity), xytext=(8, 4 - 11 * n if r.y == 1 else 4 + 11 * n),
                textcoords="offset points", fontsize=7)
ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=45, ha="right")
ax.set_ylabel("AlphaMissense score")
ax.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
fig.text(0.01, 0.98, "Labeled variants by gene; circled and named: the 9 discordant variants", fontsize=9, va="top")
fig.tight_layout(rect=(0, 0, 1, 0.93)); save(fig, "Figure3")

# ---------- Figure 4: structure adds nothing ----------
fig, ax = plt.subplots(figsize=(5.6, 4))
sets = ["AlphaMissense\nalone", "Structure\nalone", "AlphaMissense\n+ structure"]
cv = [0.972, 0.868, 0.967]; logo = [0.966, 0.833, 0.962]
xx = np.arange(3); w = 0.38
b1 = ax.bar(xx - w/2, cv, w, color="#5B7C99", label="5-fold CV (10 repeats)")
b2 = ax.bar(xx + w/2, logo, w, color="#A9C0D6", label="Leave-one-gene-out")
for b in list(b1) + list(b2):
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.005, f"{b.get_height():.3f}", ha="center", fontsize=7)
ax.set_ylim(0.6, 1.1); ax.set_xticks(xx); ax.set_xticklabels(sets); ax.set_ylabel("AUROC (y-axis starts at 0.6)")
ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=2)
fig.tight_layout(); save(fig, "Figure4")

# ---------- Figure 5: COQ4 yeast-tested variants ----------
fig, ax = plt.subplots(figsize=(5.6, 3.6))
f5 = pd.read_csv("out5/coq4_functional_variants.csv")
colmap = {"lost": C_P, "partial": C_G}
f5 = f5.sort_values("am_pathogenicity", ascending=False).reset_index(drop=True)
ax.barh(range(len(f5))[::-1], f5.am_pathogenicity, color=[colmap[c] for c in f5.yeast_call])
for i, r in zip(range(len(f5))[::-1], f5.itertuples()):
    ax.text(r.am_pathogenicity + 0.015, i, f"{r.am_pathogenicity:.3f}", va="center", fontsize=8)
ax.set_yticks(range(len(f5))[::-1]); ax.set_yticklabels([f"COQ4 {v}" for v in f5.protein_variant])
ax.axvline(HI, color="grey", ls="--", lw=0.8); ax.axvline(LO, color="grey", ls="--", lw=0.8)
ax.set_xlim(0, 1.1); ax.set_xlabel("AlphaMissense score")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=C_P, label="Yeast growth lost"), Patch(color=C_G, label="Partial function")],
          frameon=False, fontsize=8, loc="lower right")
fig.tight_layout(); save(fig, "Figure5")
print("Wrote figures/Figure1-5 (.png, .pdf)")
