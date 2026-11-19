# Data sources and versions


## Inputs (not included in this repository)
| File | Source | Retrieved | Size (bytes) | SHA-256 |
|---|---|---|---|---|
| variant_summary.txt.gz | ClinVar (NCBI) | 2026-11-18 | 458,480,521 | TODO |
| AlphaMissense_hg38.tsv.gz | AlphaMissense predictions | 2026-11-18 | 642,961,469 
| af_models/*.pdb | AlphaFold Protein Structure Database (model version 6 at download) | downloaded by step3_download_alphafold.py | - | not needed |

Compute checksums on Windows with: `certutil -hashfile variant_summary.txt.gz SHA256`

## Genes and UniProt identifiers (checked against AlphaFold database gene names)
PDSS1 Q5T2R2; PDSS2 Q86YH6; COQ2 Q96H96; COQ3 Q9NZJ6; COQ4 Q9Y3A0; COQ5 Q5HYK3; COQ6 Q9Y2Z9; COQ7 Q99807;
COQ8A Q8NI60; COQ8B Q96D53; COQ9 O75208; COQ10A Q96MF6; COQ10B Q9H8M1.

## Software
Python 3.11 (Windows). Record exact versions with `pip freeze > requirements-frozen.txt` and commit that file.
Random seed for all resampling: 20261118.

## Curated (not computed) inputs
`step5_coq4_functional.py` contains curated yeast-complementation calls for COQ4 variants. They come from
secondary sources (UniProt variant page, OMIM summary of Mero et al. 2021, an abstract of a Padova thesis)
and must be verified against the original papers before they are relied on.
