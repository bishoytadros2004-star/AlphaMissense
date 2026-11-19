"""
Step 4: does AlphaFold structure add anything beyond AlphaMissense, and how does the
VUS 'Tier 1' rule behave on the labeled variants?  Cross-validated, fixed seed.
Input: out3/structure_features.csv   Output: out4/added_value_log.txt
"""
import os
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

SEED, HI = 20261118, 0.564
os.makedirs("out4", exist_ok=True); log = []
def note(s=""): print(s); log.append(str(s))

f = pd.read_csv(os.path.join("out3", "structure_features.csv"))
lab = f[f.label.isin(["P", "B"])].reset_index(drop=True)
y = (lab.label == "P").astype(int).values
note(f"Labeled variants with AlphaMissense + structure: {len(lab)} (P={y.sum()}, B={(1-y).sum()})")

sets = {"AlphaMissense only": ["am_pathogenicity"],
        "Structure only (pLDDT+burial)": ["plddt", "neighbors"],
        "AlphaMissense + structure": ["am_pathogenicity", "plddt", "neighbors"]}
def model(): return make_pipeline(StandardScaler(), LogisticRegression(max_iter=1000))

note("\n== Test 1a: repeated 5-fold CV (10 repeats), AUROC per repeat ==")
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
res = {k: [[] for _ in range(10)] for k in sets}
# pooled out-of-fold predictions per repeat
oof = {k: np.zeros((10, len(lab))) for k in sets}
for i, (tr, te) in enumerate(cv.split(lab, y)):
    rep = i // 5
    for k, cols in sets.items():
        m = model().fit(lab.loc[tr, cols], y[tr])
        oof[k][rep, te] = m.predict_proba(lab.loc[te, cols])[:, 1]
auc = {k: np.array([roc_auc_score(y, oof[k][r]) for r in range(10)]) for k in sets}
for k in sets:
    note(f"{k:34s} mean AUROC {auc[k].mean():.3f}  (range over repeats {auc[k].min():.3f}-{auc[k].max():.3f})")
dd = auc["AlphaMissense + structure"] - auc["AlphaMissense only"]
note(f"Gain from adding structure to AlphaMissense: mean {dd.mean():+.3f} (range {dd.min():+.3f} to {dd.max():+.3f})")

note("\n== Test 1b: leave-one-gene-out (train on other genes, predict held-out gene) ==")
for k, cols in sets.items():
    p = np.full(len(lab), np.nan)
    for g in lab.gene.unique():
        te = (lab.gene == g).values
        if (y[~te] == 1).sum() < 5 or (y[~te] == 0).sum() < 5: continue
        m = model().fit(lab.loc[~te, cols], y[~te]); p[te] = m.predict_proba(lab.loc[te, cols])[:, 1]
    ok = ~np.isnan(p)
    note(f"{k:34s} AUROC {roc_auc_score(y[ok], p[ok]):.3f}  (n={ok.sum()})")

note("\n== Test 2: Tier-1 rule on labeled variants (burial cut-off learned inside each training fold) ==")
rng_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
t1 = np.zeros((10, len(lab)), bool)
for i, (tr, te) in enumerate(rng_cv.split(lab, y)):
    rep = i // 5
    cut = lab.loc[tr][y[tr] == 1].neighbors.quantile(0.25)
    d = lab.loc[te]
    t1[rep, te] = ((d.am_pathogenicity > HI) & (d.plddt >= 70) & (d.neighbors >= cut)).values
am_only = (lab.am_pathogenicity > HI).values
def sens_spec(pred):
    return pred[y == 1].mean(), 1 - pred[y == 0].mean(), (pred & (y == 1)).sum() / max(pred.sum(), 1)
s = [sens_spec(t1[r]) for r in range(10)]
note(f"AM > {HI} only         : sensitivity {sens_spec(am_only)[0]:.3f}, specificity {sens_spec(am_only)[1]:.3f}, precision {sens_spec(am_only)[2]:.3f}")
note(f"Tier-1 rule (mean/10 reps): sensitivity {np.mean([a[0] for a in s]):.3f}, specificity {np.mean([a[1] for a in s]):.3f}, precision {np.mean([a[2] for a in s]):.3f}")
note("(precision here uses the 46.5% pathogenic prevalence of the labeled set; real VUS prevalence is unknown and likely much lower)")
open(os.path.join("out4", "added_value_log.txt"), "w").write("\n".join(log))
