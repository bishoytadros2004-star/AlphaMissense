# AlphaMissense on CoQ10 biosynthesis genes

Code, logs and figures for a short report evaluating AlphaMissense on ClinVar missense variants in 13 CoQ10
biosynthesis genes (PDSS1, PDSS2, COQ2-COQ7, COQ8A, COQ8B, COQ9, COQ10A, COQ10B).

**Status:** manuscript in preparation. Archive DOI: TODO (Zenodo).

## Reproduce
1. Download the two large inputs listed in `DATA_SOURCES.md` into the repository root (they are git-ignored).
2. Run `run_all.bat` (Windows) or run the steps below in order. Steps 1 and 5 scan large files and take minutes.
3. Outputs appear in `out/`, `out2/` ... `out5/` and `figures/`. Fixed seed: 20261118.

| Step | Script | What it does |
|---|---|---|
| 1 | `step1_labels.py` | Missense-only, GRCh38, review-star filtered ClinVar labels; joins AlphaMissense |
| 2 | `step2_analysis.py` | AUROC/AUPRC with stratified bootstrap CIs, within-gene analysis, VUS bands |
| 3a | `step3_download_alphafold.py` | Downloads AlphaFold models and checks UniProt IDs against gene names |
| 3b | `step3_structure.py` | pLDDT and burial features with a reference-residue check |
| 4 | `step4_added_value.py` | Cross-validated test of whether structure adds to AlphaMissense |
| 5 | `step5_coq4_functional.py` | AlphaMissense scores for COQ4 variants with yeast evidence (curated calls) |
| 6 | `step6_figures.py` | Manuscript figures (PNG and PDF) |

## What is in `results/`
Final logs and tables from the analysis run: `results_summary.md` (all numbers used in the paper), `run_log.txt`,
`label_counts.csv`, `metrics.csv`, `analysis_log.txt`, `added_value_log.txt`, `discordant_variants_table.csv`,
`coq4_functional_variants.csv`.

## Not included, on purpose
Raw ClinVar/AlphaMissense files, AlphaFold models, and per-variant AlphaMissense score tables for all variants,
pending a check of the AlphaMissense license.

## Limitations
- Small labeled set (149 variants); per-gene estimates are not possible.
- Functional (yeast) calls are curated from secondary sources and must be verified before use.
- Only AlphaMissense is evaluated.

## AI assistance
An AI assistant (Claude, Anthropic) helped write and test these scripts. The author reviewed the analyses and
takes responsibility for the content.

## License
Code: MIT (see `LICENSE`). Third-party data keep their own licenses.
