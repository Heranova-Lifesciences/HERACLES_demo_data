# Figure2

Panels A and B are taken from the original figure (`source_fig2_original.pdf`); panel C is the newly drawn paired p-value plot built from `TableS20_fig2C_pairs_stratified.csv`.

## code
- `compose_main_figures.py`
- `make_mainfig_panels.py`

## data
- `TableS19_marker_benchmark.csv` — 4 KB
- `TableS20_fig2C_pairs_stratified.csv` — 11 KB
- `source_fig2_original.pdf` — 555 KB

## figure
- `Figure2_revised.pdf`
- `Figure2_revised.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
