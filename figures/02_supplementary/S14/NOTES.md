# S14

Nominal p of each published marker under the full fit and under every leave-one-library-out refit, cluster level. Cell colour is -log10 p, open circles mark BH-adjusted p < 0.05, and NA marks a refit in which the marker was filtered out. Row labels follow the revised Figure 4 marker names.

## code
- `make_figR3_2C_LOO.py`

## data
- `TableS22_LOO_marker_stability.csv` — 18 KB

## figure
- `S14_LOO_stability.pdf`
- `S14_LOO_stability.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
