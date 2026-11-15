# S08

Effect of the minimum length ratio (0.90, 0.94, 0.95, 0.98) on cluster topology (A) and on the number of nominally significant DE clusters (B).

## code
- `make_figS15.py`

## data
- `threshold_summary.json` — 0 KB

## figure
- `S8_min_length_ratio.pdf`
- `S8_min_length_ratio.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
