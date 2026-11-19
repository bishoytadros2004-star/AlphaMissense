"""
Step 1: build the labeled missense variant table for CoQ10 pathway genes.

Run from the folder that holds the two downloaded files:
    python step1_labels.py
    python step1_labels.py --min-stars 2          # stricter label quality
    python step1_labels.py --skip-alphamissense   # counts only (fast)

Outputs (in ./out):
    clinvar_missense_labels.csv   one row per variant (label P / B / VUS)
    label_counts.csv              per-gene counts, missense only
    coq_with_alphamissense.csv    labels joined to AlphaMissense (unless skipped)
    run_log.txt                   parameters and filter counts, for the Methods section
"""
import argparse
import gzip
import os
import re
import sys

import pandas as pd

GENES = ["PDSS1", "PDSS2", "COQ2", "COQ3", "COQ4", "COQ5", "COQ6", "COQ7",
         "COQ8A", "COQ8B", "COQ9", "COQ10A", "COQ10B"]

CLINVAR = "variant_summary.txt.gz"
ALPHAMIS = "AlphaMissense_hg38.tsv.gz"
OUT = "out"

STARS = {
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classification for the individual variant": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}

P_TERMS = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
B_TERMS = {"Benign", "Likely benign", "Benign/Likely benign"}
VUS_TERMS = {"Uncertain significance"}

# true amino-acid substitution, e.g. p.Arg123Trp (excludes Ter, fs, del, dup, =, ext)
MISSENSE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")


def is_missense(name):
    m = MISSENSE_RE.search(str(name))
    if not m:
        return False
    return m.group(1) != "Ter" and m.group(3) not in ("Ter",) and m.group(1) != m.group(3)


def label_of(sig):
    if sig in P_TERMS:
        return "P"
    if sig in B_TERMS:
        return "B"
    if sig in VUS_TERMS:
        return "VUS"
    return None  # conflicting, drug response, risk factor, etc. -> excluded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-stars", type=int, default=1,
                    help="minimum ClinVar review stars for P and B labels (default 1)")
    ap.add_argument("--skip-alphamissense", action="store_true")
    args = ap.parse_args()

    for f in [CLINVAR] + ([] if args.skip_alphamissense else [ALPHAMIS]):
        if not os.path.exists(f):
            sys.exit(f"ERROR: missing input file {f} (run from the data folder)")
    os.makedirs(OUT, exist_ok=True)
    log = []

    def note(msg):
        print(msg)
        log.append(msg)

    note(f"min review stars for P/B: {args.min_stars}")

    cols = ["GeneSymbol", "ClinicalSignificance", "ReviewStatus", "Assembly", "Type",
            "Name", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    # read in chunks and keep only target-gene rows, so the full 9M-row table
    # never sits in memory (a plain read_csv ran out of memory on this machine)
    gene_pat = "|".join(GENES)
    parts = []
    n_read = 0
    for chunk in pd.read_csv(CLINVAR, sep="\t", compression="gzip", usecols=cols,
                             dtype=str, chunksize=250_000):
        n_read += len(chunk)
        hit = chunk["GeneSymbol"].fillna("").str.contains(gene_pat, regex=True)
        if hit.any():
            parts.append(chunk[hit])
    df = pd.concat(parts, ignore_index=True)
    note(f"ClinVar rows read: {n_read:,}")

    # exact gene match; some rows list several genes separated by ';'
    gene_set = set(GENES)
    df["gene"] = df["GeneSymbol"].fillna("").apply(
        lambda s: next((g for g in s.split(";") if g in gene_set), None))
    df = df[df["gene"].notna()]
    note(f"rows in target genes: {len(df):,}")

    df = df[df["Assembly"] == "GRCh38"]
    note(f"after GRCh38 filter: {len(df):,}")

    df = df[df["Type"] == "single nucleotide variant"]
    note(f"after SNV filter: {len(df):,}")

    df = df[df["Name"].apply(is_missense)]
    note(f"after true-missense filter: {len(df):,}")

    df["label"] = df["ClinicalSignificance"].apply(label_of)
    df = df[df["label"].notna()]
    note(f"after removing conflicting/other classes: {len(df):,}")

    df["stars"] = df["ReviewStatus"].map(STARS)
    unknown = df[df["stars"].isna()]["ReviewStatus"].unique()
    if len(unknown):
        note(f"WARNING unrecognized ReviewStatus values (treated as 0 stars): {list(unknown)}")
    df["stars"] = df["stars"].fillna(0).astype(int)
    # stars gate applies to P and B; VUS kept for the later ranking step
    keep = (df["label"] == "VUS") | (df["stars"] >= args.min_stars)
    df = df[keep]
    note(f"after review-star filter: {len(df):,}")

    df["chrom"] = df["Chromosome"].astype(str)
    df["pos"] = pd.to_numeric(df["PositionVCF"], errors="coerce")
    df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)
    df = df.rename(columns={"ReferenceAlleleVCF": "ref", "AlternateAlleleVCF": "alt"})
    df["chrom"] = df["chrom"].apply(lambda c: c if c.startswith("chr") else "chr" + c)
    df = df.drop_duplicates(subset=["gene", "chrom", "pos", "ref", "alt"])
    note(f"unique variants: {len(df):,}")

    out_cols = ["gene", "chrom", "pos", "ref", "alt", "Name", "label", "stars"]
    df[out_cols].to_csv(f"{OUT}/clinvar_missense_labels.csv", index=False)

    counts = (df.pivot_table(index="gene", columns="label", values="pos",
                             aggfunc="count", fill_value=0)
                .reindex(GENES).fillna(0).astype(int))
    for c in ("P", "B", "VUS"):
        if c not in counts:
            counts[c] = 0
    counts = counts[["P", "B", "VUS"]]
    counts.to_csv(f"{OUT}/label_counts.csv")
    note("\nMissense-only counts per gene (P = P/LP, B = B/LB):")
    note(counts.to_string())
    note(f"\nPOOLED  P={counts['P'].sum()}  B={counts['B'].sum()}  VUS={counts['VUS'].sum()}")

    if not args.skip_alphamissense:
        keys = df[["chrom", "pos", "ref", "alt"]].drop_duplicates()
        hits = []
        total = 0
        # find header row (lines starting '#' are comments; the column header is '#CHROM...')
        with gzip.open(ALPHAMIS, "rt") as fh:
            skip = 0
            for line in fh:
                if line.startswith("#CHROM"):
                    break
                skip += 1
        reader = pd.read_csv(ALPHAMIS, sep="\t", compression="gzip", skiprows=skip,
                             chunksize=2_000_000, dtype=str)
        for chunk in reader:
            chunk.columns = [c.lstrip("#") for c in chunk.columns]
            total += len(chunk)
            chunk = chunk.rename(columns={"CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt"})
            chunk["pos"] = chunk["pos"].astype(int)
            # vectorized join on genomic key: fast, and keeps only our variants
            sel = chunk.merge(keys, on=["chrom", "pos", "ref", "alt"], how="inner")
            if len(sel):
                hits.append(sel)
            if total % 20_000_000 < 2_000_000:
                print(f"  scanned {total:,} AlphaMissense rows, matched {sum(len(h) for h in hits):,}")
        note(f"\nAlphaMissense rows scanned: {total:,}")
        if not hits:
            sys.exit("ERROR: no AlphaMissense matches; check chromosome naming / genome build")
        am = pd.concat(hits)
        am["am_pathogenicity"] = am["am_pathogenicity"].astype(float)
        dup = am.duplicated(subset=["chrom", "pos", "ref", "alt"], keep=False).sum()
        note(f"AlphaMissense rows with multiple transcripts for same variant: {dup}")
        am = (am.groupby(["chrom", "pos", "ref", "alt"], as_index=False)
                .agg(am_pathogenicity=("am_pathogenicity", "mean"),
                     am_class=("am_class", "first"),
                     protein_variant=("protein_variant", "first")))
        merged = df[out_cols].merge(am, on=["chrom", "pos", "ref", "alt"], how="left")
        missing = merged["am_pathogenicity"].isna().sum()
        note(f"variants with no AlphaMissense score: {missing} of {len(merged)}")
        merged.to_csv(f"{OUT}/coq_with_alphamissense.csv", index=False)

    with open(f"{OUT}/run_log.txt", "w") as fh:
        fh.write("\n".join(log))
    print(f"\nWrote outputs to ./{OUT}/")


if __name__ == "__main__":
    main()
