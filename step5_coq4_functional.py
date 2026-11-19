"""
Step 5: AlphaMissense scores for the COQ4 variants that have yeast complementation evidence,
including G95D and P193S, which are not in our ClinVar-labeled set.
Run from the data folder:  python step5_coq4_functional.py
Output: out5/coq4_functional_variants.csv
Evidence calls below are CURATED, not measured here: confirm each against the cited source.
"""
import gzip, os
import pandas as pd

UNIPROT = "Q9Y3A0"   # COQ4 (verified by step3_download_alphafold.py)
EVIDENCE = {
    "R240C": ("lost", "UniProt variant page; Brea-Calvo 2015 AJHG; Mero 2021 via OMIM (unverified)"),
    "P193S": ("lost", "Mero 2021 via OMIM summary (unverified; primary paper not read)"),
    "G95D":  ("lost", "Mero 2021 via OMIM summary (unverified; primary paper not read)"),
    "R102H": ("partial", "Sarto thesis Padova 2023/24 abstract; Mero 2021 via OMIM (unverified)"),
    "C231Y": ("partial", "Sarto thesis Padova 2023/24 abstract"),
}
os.makedirs("out5", exist_ok=True)
with gzip.open("AlphaMissense_hg38.tsv.gz", "rt") as fh:
    skip = 0
    for line in fh:
        if line.startswith("#CHROM"):
            break
        skip += 1
hits = []
for ch in pd.read_csv("AlphaMissense_hg38.tsv.gz", sep="\t", skiprows=skip, chunksize=2_000_000, dtype=str):
    ch.columns = [c.lstrip("#") for c in ch.columns]
    ucol = next((c for c in ch.columns if "uniprot" in c.lower()), None)
    if ucol is None:
        raise SystemExit(f"No UniProt column found; columns are {list(ch.columns)}")
    sel = ch[(ch[ucol] == UNIPROT) & (ch["protein_variant"].isin(EVIDENCE))]
    if len(sel):
        hits.append(sel)
if not hits:
    raise SystemExit("No matches found; check UniProt ID / column names")
am = pd.concat(hits)
am["am_pathogenicity"] = am["am_pathogenicity"].astype(float)
am = am.groupby("protein_variant", as_index=False).agg(am_pathogenicity=("am_pathogenicity", "mean"),
                                                        am_class=("am_class", "first"), n_rows=("am_pathogenicity", "size"))
am["yeast_call"] = am.protein_variant.map(lambda v: EVIDENCE[v][0])
am["evidence_source"] = am.protein_variant.map(lambda v: EVIDENCE[v][1])
am.to_csv("out5/coq4_functional_variants.csv", index=False)
print(am.to_string(index=False))
missing = set(EVIDENCE) - set(am.protein_variant)
if missing:
    print("Not found in AlphaMissense:", sorted(missing))
