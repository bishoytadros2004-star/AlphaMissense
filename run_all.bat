@echo off
REM Reproduces every result. Put variant_summary.txt.gz and AlphaMissense_hg38.tsv.gz in this folder first.
REM Steps 1 and 5 scan the large files and take several minutes each.
pip install -r requirements.txt
python step1_labels.py
python step2_analysis.py
python step3_download_alphafold.py
python step3_structure.py
python step4_added_value.py
python step5_coq4_functional.py
python step6_figures.py
