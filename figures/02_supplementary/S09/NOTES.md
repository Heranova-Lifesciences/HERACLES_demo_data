# S09

Marker-level re-analysis across the four tools, redrawn from Figure 4A/B, with degenerate fits and all-zero groups marked. Rows are published markers, columns are tools, and each cell gives the log2 fold change with its nominal p.

## code
- `make_fig4ab_v4_panels.py`

## data
- `TableS19_marker_benchmark.csv` — 4 KB
- `TableS26_comparator_marker_rerun.csv` — 8 KB

## figure
- `S9_comparator_rerun.pdf`
- `S9_comparator_rerun.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
