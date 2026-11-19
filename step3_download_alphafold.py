"""
Step 3a: download AlphaFold models for the 13 CoQ-pathway proteins and verify each
UniProt ID really belongs to the expected gene (IDs below are from memory; the
check will flag any that are wrong instead of silently using the wrong protein).

Run on your own computer (needs internet):  python step3_download_alphafold.py
Output: af_models/<GENE>_<UNIPROT>.pdb  and  af_models/download_log.csv
"""
import csv, json, os, urllib.request

TARGETS = {
    "PDSS1": "Q5T2R2", "PDSS2": "Q86YH6", "COQ2": "Q96H96", "COQ3": "Q9NZJ6",
    "COQ4": "Q9Y3A0", "COQ5": "Q5HYK3", "COQ6": "Q9Y2Z9", "COQ7": "Q99807",
    "COQ8A": "Q8NI60", "COQ8B": "Q96D53", "COQ9": "O75208",
    "COQ10A": "Q96MF6", "COQ10B": "Q9H8M1",
}
ALIASES = {"COQ8A": {"COQ8A", "ADCK3"}, "COQ8B": {"COQ8B", "ADCK4"}}
os.makedirs("af_models", exist_ok=True)
rows = []
for gene, acc in TARGETS.items():
    try:
        with urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}", timeout=60) as r:
            meta = json.load(r)[0]
        api_gene = str(meta.get("gene", "")).upper()
        ok = api_gene in ALIASES.get(gene, {gene})
        status = "OK" if ok else f"GENE MISMATCH (API says {api_gene or 'unknown'})"
        seq_len = len(meta.get("uniprotSequence", ""))
        if ok:
            pdb_url = meta["pdbUrl"]
            with urllib.request.urlopen(pdb_url, timeout=120) as r2:
                open(f"af_models/{gene}_{acc}.pdb", "wb").write(r2.read())
        rows.append([gene, acc, api_gene, seq_len, meta.get("latestVersion", ""), status])
    except Exception as e:
        rows.append([gene, acc, "", "", "", f"FAILED: {e}"])
    print(rows[-1])
with open("af_models/download_log.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["gene", "uniprot", "api_gene", "length", "version", "status"]); w.writerows(rows)
print("\nDone. Anything not 'OK' must be fixed before we use it.")
