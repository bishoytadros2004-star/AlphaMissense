# Results summary (single source of truth for the Results section)
Data: ClinVar variant_summary (downloaded 2026-11-18), AlphaMissense hg38 scores, AlphaFold DB models (13 proteins).
All numbers below were produced by the scripts in this folder; seed 20261118.

## Dataset (step1_labels.py)
- 9,355,006 ClinVar rows -> 8,126 in 13 CoQ-pathway genes -> GRCh38 4,037 -> SNV 3,667 -> true missense 1,718
  -> P/LP, B/LB or VUS only 1,563 -> review stars >=1 for P/B 1,544 unique variants.
- Labeled: 69 P/LP, 80 B/LB; VUS 1,395. AlphaMissense rows scanned: 71,697,556.
- 31 variants unscored by AlphaMissense (5 labeled, 26 VUS). Labeled unscored: COQ2 L25V, COQ2 G29C, COQ6 S21F (N-terminal),
  COQ4 M1T, COQ7 M1L (start-loss). 23 COQ2 VUS lie in the N-terminal ~50 residues (genomic chr4:83284763-83284899).
  AlphaFold/UniProt COQ2 model is 371 aa; COQ2 has multiple start codons (hypothesis: scored isoform starts later; NOT confirmed).
- Genes with >=10 P and >=10 B: only COQ8A (30/12) near; per-gene benchmark not feasible.

## AlphaMissense discrimination (step2_analysis.py; stratified bootstrap, 2000 reps)
- Pooled (n=144: 67 P, 77 B): AUROC 0.976 (95% CI 0.953-0.992); AUPRC 0.968 (0.941-0.990), prevalence 0.465.
- Within-gene only (630 same-gene P-B pairs): AUROC 0.981 (0.951-0.998).
- Excluding COQ8A (n=102): 0.958 (0.923-0.987). COQ8A alone: 1.000 (degenerate CI).
- Review stars >=2 only (n=64: 25 P, 39 B): AUROC 0.963 (0.917-0.994).
- Published cut-offs: ambiguous band 0.34-0.564 holds 5/144 (3.5%); sensitivity 59/67 = 0.881; specificity (<0.34) 71/77 = 0.922.
- 9 discordant labeled variants: see discordant_variants_table.csv.

## VUS (AlphaMissense only; not clinical reclassification)
- 1,369 scored; 314 in likely-pathogenic band; 26 unscored. Ranked list: vus_ranked.csv.

## Structure (step3_structure.py, step4_added_value.py)
- 1,513/1,513 scored variants matched reference residue in AlphaFold models (0 dropped).
- Alone: burial (neighbors within 10 A) AUROC 0.879; pLDDT 0.688.
- 10x repeated 5-fold CV: AlphaMissense 0.972; structure only 0.868; AlphaMissense+structure 0.967 (change -0.005).
- Leave-one-gene-out: 0.966 / 0.833 / 0.962.
- Exploratory 'Tier 1' rule (AM>0.564, pLDDT>=70, buried): sensitivity 0.690 vs 0.881 for AM alone; specificity 0.958 vs 0.935. Not recommended.

## Limitations to state
Small labeled set (149); ClinVar labels may be noisy/circular; one predictor only; binary labels hide partial-function variants;
functional evidence sparse (R240C loss of function; R102H partial); N-terminal COQ2 / start-loss variants unscored.
