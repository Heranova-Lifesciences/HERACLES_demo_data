# S18

What the mapping flag changes. (A) Fragments detected under -v 0 and -v 1, split into shared and run-specific blocks. (B) Overlap of the two nominally significant sets. The 14,888 fragments unique to -v 1 carry 1.0% of the reads mapped to tRNA references.

## code
- `make_figR3_3p3_mapping_v2.py`

## data
- `ablation_v2_report.json` — 1 KB

## figure
- `S18_mapping_ablation.pdf`
- `S18_mapping_ablation.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
