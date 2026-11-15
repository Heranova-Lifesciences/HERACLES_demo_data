# S11

The 102 features tested in both the fragment- and cluster-level matrices, stratified by whether the representative inherits any members: paired -log10 p (A) and the gain per feature (B).

## code
- `make_figR3_2C_LOO.py`

## data
- `TableS20_fig2C_pairs_stratified.csv` — 11 KB

## figure
- `S11_stratified.pdf`
- `S11_stratified.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
