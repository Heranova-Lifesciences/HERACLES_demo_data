# S04

Permuted-label controls: nominal calls under the observed assignment against permuted labels, for the clustering and no-clustering workflows, in the in-house cohort and the two public datasets. Blue bars are observed, grey bars the permuted median with the full range as an error bar.

## code
- `make_figR1_permutation_v2.py`

## data
- `TableS16_permuted_label_controls.csv` — 13 KB
- `TableS25b_permutation_two_workflows.csv` — 0 KB

## figure
- `S4_permutation.pdf`
- `S4_permutation.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
