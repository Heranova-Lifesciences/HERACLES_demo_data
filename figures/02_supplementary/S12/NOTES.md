# S12

FDR-positive features against the number of tests, for all three datasets at both aggregation levels.

## code
- `make_figR3_2p2_fdr_scale.py`

## data
- `TableS25_permuted_label_controls.csv` — 13 KB

## figure
- `S12_fdr_vs_tests.pdf`
- `S12_fdr_vs_tests.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
