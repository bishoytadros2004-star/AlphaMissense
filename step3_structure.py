"""
Step 3b: add AlphaFold structural context to every scored variant.

Needs: out/coq_with_alphamissense.csv (step 1) and af_models/*.pdb (step 3a).
Run  : python step3_structure.py
Out  : out3/structure_features.csv, out3/structure_log.txt

Features per variant (all from the AlphaFold model; numpy only, no extra installs):
  plddt     per-residue AlphaFold confidence (stored in the PDB B-factor column)
  neighbors number of other residues whose CB (CA for Gly) lies within 10 A = burial proxy
Numbering safety: residue numbers come from AlphaMissense's protein_variant (UniProt
numbering, same as the AlphaFold models). The script checks the reference amino acid
against the model and DROPS any variant that does not match, so a numbering mismatch
can never be used silently.
"""
import glob, os, re
import numpy as np
import pandas as pd

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H",
       "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
       "TYR":"Y","VAL":"V"}
HI, LO, CUT = 0.564, 0.34, 10.0
OUT = "out3"; os.makedirs(OUT, exist_ok=True)
log = []
def note(s=""):
    print(s); log.append(str(s))

def auroc(y, s):
    y = np.asarray(y); s = np.asarray(s, dtype=float)
    p, n = s[y == 1], s[y == 0]
    return ((p[:, None] > n[None, :]).sum() + 0.5 * (p[:, None] == n[None, :]).sum()) / (len(p) * len(n))

def load_model(path):
    res = {}
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        name = line[12:16].strip(); rn = line[17:20]; num = int(line[22:26])
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        b = float(line[60:66])
        r = res.setdefault(num, {"aa": AA3.get(rn, "X"), "plddt": b})
        if name == "CA": r["ca"] = xyz
        if name == "CB": r["cb"] = xyz
    nums = sorted(res)
    coords = np.array([res[n].get("cb", res[n].get("ca")) for n in nums])
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    nb = (d < CUT).sum(axis=1) - 1
    return {n: dict(aa=res[n]["aa"], plddt=res[n]["plddt"], neighbors=int(nb[i])) for i, n in enumerate(nums)}

d = pd.read_csv(os.path.join("out", "coq_with_alphamissense.csv"))
d = d.dropna(subset=["am_pathogenicity", "protein_variant"]).copy()
pv = d.protein_variant.str.extract(r"^([A-Z])(\d+)([A-Z])$")
d["ref_aa"], d["resnum"], d["alt_aa"] = pv[0], pd.to_numeric(pv[1]), pv[2]
d = d.dropna(subset=["resnum"]); d["resnum"] = d["resnum"].astype(int)

models = {}
for p in glob.glob(os.path.join("af_models", "*.pdb")):
    models[os.path.basename(p).split("_")[0]] = load_model(p)
note(f"Models loaded: {sorted(models)}")

feat, bad = [], 0
for _, r in d.iterrows():
    m = models.get(r.gene, {}).get(r.resnum)
    if m is None or m["aa"] != r.ref_aa:
        bad += 1; continue
    feat.append(dict(r, plddt=m["plddt"], neighbors=m["neighbors"]))
f = pd.DataFrame(feat)
note(f"Variants with AlphaMissense score: {len(d)}; matched to model with identical reference residue: {len(f)}; "
     f"dropped for numbering/residue mismatch: {bad}")
f.to_csv(os.path.join(OUT, "structure_features.csv"), index=False)

lab = f[f.label.isin(["P", "B"])].copy(); lab["y"] = (lab.label == "P").astype(int)
note("\n== Do structure features alone separate P from B? (AUROC, unadjusted; descriptive) ==")
for c in ["plddt", "neighbors", "am_pathogenicity"]:
    note(f"{c:18s} AUROC {auroc(lab.y, lab[c]):.3f}   median P={lab[lab.y==1][c].median():.1f}  median B={lab[lab.y==0][c].median():.1f}")

note("\n== Discordant labeled variants with structure ==")
disc = lab[((lab.y == 0) & (lab.am_pathogenicity > HI)) | ((lab.y == 1) & (lab.am_pathogenicity < LO))]
note(disc[["gene", "protein_variant", "label", "am_pathogenicity", "plddt", "neighbors", "stars"]].to_string(index=False))

# exploratory VUS tiers; burial cut-off = 25th percentile of pathogenic variants (heuristic, data-derived)
nb_cut = lab[lab.y == 1].neighbors.quantile(0.25)
v = f[f.label == "VUS"].copy()
conf = v.plddt >= 70
v["tier"] = np.select(
    [(v.am_pathogenicity > HI) & conf & (v.neighbors >= nb_cut),
     (v.am_pathogenicity > HI)],
    ["1_AM-pathogenic+confident+buried", "2_AM-pathogenic_only"], default="3_other")
note(f"\nVUS tiers (exploratory heuristic; burial cut-off neighbors >= {nb_cut:.0f}, pLDDT >= 70):")
note(v.tier.value_counts().sort_index().to_string())
v.sort_values(["tier", "am_pathogenicity"], ascending=[True, False]).to_csv(os.path.join(OUT, "vus_tiered.csv"), index=False)
open(os.path.join(OUT, "structure_log.txt"), "w").write("\n".join(log))
