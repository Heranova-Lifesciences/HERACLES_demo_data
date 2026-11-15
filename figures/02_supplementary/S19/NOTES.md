# S19

Label-aware re-clustering. Condition-specific fragments and clusters in the three datasets (A), with the permuted-label control and the union-of-representatives sensitivity analysis (B).

## code
- `make_figR3_3p4_labelaware_v4.py`

## data
- (none: values are hard-coded in the script)

## figure
- `S19_labelaware.pdf`
- `S19_labelaware.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
