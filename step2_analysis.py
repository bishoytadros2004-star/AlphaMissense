"""
Step 2: benchmark AlphaMissense on ClinVar-labeled CoQ-pathway missense variants.

Input : out/coq_with_alphamissense.csv  (from step1_labels.py)
Output: out2/metrics.csv, out2/vus_ranked.csv, out2/analysis_log.txt, out2/roc_pr.png
Run   : python step2_analysis.py
All randomness uses a fixed seed (SEED) so results are exactly reproducible.
"""
import os
import numpy as np
import pandas as pd

SEED = 20261118
N_BOOT = 2000
# published AlphaMissense class cut-offs: >0.564 likely pathogenic, <0.34 likely benign
HI, LO = 0.564, 0.34
IN = os.path.join("out", "coq_with_alphamissense.csv")
OUT = "out2"
os.makedirs(OUT, exist_ok=True)
log = []
def note(s=""):
    print(s); log.append(str(s))

def auroc(y, s):
    """AUROC via rank statistic (ties = 0.5)."""
    y = np.asarray(y); s = np.asarray(s)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return (gt + 0.5 * eq) / (len(pos) * len(neg))

def auprc(y, s):
    """Average precision (step-wise), higher score = more pathogenic."""
    y = np.asarray(y); s = np.asarray(s)
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y); k = np.arange(1, len(y) + 1)
    prec = tp / k
    return (prec * y).sum() / max(y.sum(), 1)

def within_gene_auroc(df):
    """Pooled within-gene concordance: only P-vs-B pairs from the SAME gene.
    Removes any advantage from simply scoring a whole gene higher."""
    num = den = 0.0
    for _, g in df.groupby("gene"):
        pos = g.loc[g.y == 1, "s"].values; neg = g.loc[g.y == 0, "s"].values
        if len(pos) == 0 or len(neg) == 0:
            continue
        num += (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
        den += len(pos) * len(neg)
    return num / den if den else np.nan, den

def boot_ci(df, fn, rng):
    """Stratified bootstrap: resample within each gene x label cell."""
    cells = [g for _, g in df.groupby(["gene", "y"])]
    vals = []
    for _ in range(N_BOOT):
        b = pd.concat([c.iloc[rng.integers(0, len(c), len(c))] for c in cells])
        v = fn(b)
        vals.append(v[0] if isinstance(v, tuple) else v)
    vals = np.array(vals, dtype=float)
    return np.nanpercentile(vals, [2.5, 97.5])

d = pd.read_csv(IN)
note(f"Input rows: {len(d)}")
lab = d[d.label.isin(["P", "B"])].copy()
vus = d[d.label == "VUS"].copy()
lab["y"] = (lab.label == "P").astype(int)
n_miss = lab.am_pathogenicity.isna().sum()
note(f"Labeled variants: {len(lab)} (P={int(lab.y.sum())}, B={int((1-lab.y).sum())}); "
     f"without AlphaMissense score (excluded): {n_miss}")
note(f"  missing by class: P={int(lab[lab.am_pathogenicity.isna()].y.sum())}, "
     f"B={int((1-lab[lab.am_pathogenicity.isna()].y).sum())}")
lab = lab.dropna(subset=["am_pathogenicity"]).rename(columns={"am_pathogenicity": "s"})
rng = np.random.default_rng(SEED)

rows = []
def report(name, df, fn_auc=auroc, within=False):
    if within:
        a, npairs = within_gene_auroc(df)
        ci = boot_ci(df, within_gene_auroc, rng)
        extra = f"{int(npairs)} same-gene P-B pairs"
    else:
        a = fn_auc(df.y, df.s)
        ci = boot_ci(df, lambda b: auroc(b.y, b.s), rng)
        extra = f"n={len(df)} (P={int(df.y.sum())}, B={int((1-df.y).sum())})"
    note(f"{name:42s} AUROC {a:.3f} (95% CI {ci[0]:.3f}-{ci[1]:.3f})  [{extra}]")
    rows.append(dict(analysis=name, metric="AUROC", value=a, ci_lo=ci[0], ci_hi=ci[1], detail=extra))

note("\n== Discrimination ==")
report("Pooled, all genes", lab)
report("Within-gene only (gene-adjusted)", lab, within=True)
report("Pooled, excluding COQ8A", lab[lab.gene != "COQ8A"])
report("COQ8A alone", lab[lab.gene == "COQ8A"])

ap = auprc(lab.y, lab.s)
prev = lab.y.mean()
ci_ap = boot_ci(lab, lambda b: auprc(b.y, b.s), rng)
note(f"{'Pooled AUPRC':42s} {ap:.3f} (95% CI {ci_ap[0]:.3f}-{ci_ap[1]:.3f}); "
     f"no-skill baseline = prevalence {prev:.3f}")
rows.append(dict(analysis="Pooled AUPRC", metric="AUPRC", value=ap, ci_lo=ci_ap[0], ci_hi=ci_ap[1],
                 detail=f"baseline {prev:.3f}"))

note("\n== Published AlphaMissense cut-offs on ClinVar-labeled variants ==")
called_p = lab.s > HI; called_b = lab.s < LO
amb = ~(called_p | called_b)
note(f"Ambiguous band ({LO}-{HI}): {amb.sum()} of {len(lab)} labeled variants "
     f"({100*amb.mean():.1f}%); P in band={int(lab[amb].y.sum())}, B in band={int((1-lab[amb].y).sum())}")
tp = int(((called_p) & (lab.y == 1)).sum()); fn = int(((~called_p) & (lab.y == 1)).sum())
tn = int(((called_b) & (lab.y == 0)).sum()); fp = int(((~called_b) & (lab.y == 0)).sum())
note(f"Pathogenic called (>{HI}): sensitivity {tp}/{tp+fn} = {tp/(tp+fn):.3f}")
note(f"Benign called (<{LO}): specificity {tn}/{tn+fp} = {tn/(tn+fp):.3f}  (ambiguous counted as not-benign)")

note("\n== Benign variants scored as pathogenic (possible mislabels or model errors) ==")
bb = lab[(lab.y == 0) & (lab.s > HI)][["gene", "Name", "s", "stars"]]
note(bb.to_string(index=False) if len(bb) else "none")
note("\n== Pathogenic variants scored below 0.34 (model misses) ==")
pm = lab[(lab.y == 1) & (lab.s < LO)][["gene", "Name", "s", "stars"]]
note(pm.to_string(index=False) if len(pm) else "none")

note("\n== VUS (AlphaMissense only; NOT a clinical reclassification) ==")
v = vus.dropna(subset=["am_pathogenicity"]).copy()
note(f"VUS scored: {len(v)} of {len(vus)}; unscored: {len(vus)-len(v)}")
v["am_band"] = np.where(v.am_pathogenicity > HI, "likely_pathogenic",
                np.where(v.am_pathogenicity < LO, "likely_benign", "ambiguous"))
tab = v.groupby(["gene", "am_band"]).size().unstack(fill_value=0)
note(tab.to_string())
note(f"VUS in likely_pathogenic band: {(v.am_band=='likely_pathogenic').sum()} of {len(v)}")
v.sort_values("am_pathogenicity", ascending=False).to_csv(os.path.join(OUT, "vus_ranked.csv"), index=False)

pd.DataFrame(rows).to_csv(os.path.join(OUT, "metrics.csv"), index=False)

try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
    thr = np.unique(np.r_[lab.s.values, 0, 1])[::-1]
    P = lab.y.sum(); N = len(lab) - P
    tpr = [((lab.s >= t) & (lab.y == 1)).sum() / P for t in thr]
    fpr = [((lab.s >= t) & (lab.y == 0)).sum() / N for t in thr]
    ax[0].plot(fpr, tpr); ax[0].plot([0, 1], [0, 1], "k:")
    ax[0].set_xlabel("False positive rate"); ax[0].set_ylabel("True positive rate"); ax[0].set_title("ROC")
    for lbl, c in [(0, "tab:blue"), (1, "tab:red")]:
        ax[1].hist(lab[lab.y == lbl].s, bins=20, alpha=.6, color=c, label="B/LB" if lbl == 0 else "P/LP")
    ax[1].set_xlabel("AlphaMissense score"); ax[1].legend(); ax[1].set_title("Score distribution")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "roc_pr.png"), dpi=200)
except Exception as e:
    note(f"(plot skipped: {e})")

note(f"\nSeed={SEED}, bootstrap replicates={N_BOOT}")
open(os.path.join(OUT, "analysis_log.txt"), "w").write("\n".join(log))
