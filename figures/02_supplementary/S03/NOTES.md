# S03

Cluster-size structure. (A) Number of clusters by member count, with the 10,089 singleton clusters (59.6%) marked. (B) Total cluster count distribution.

## code
- `make_figures.py`

## data
- `cluster_summary_detailed.csv` — 19.1 MB, in `00_shared_data/`

## figure
- `S3_singleton_structure.pdf`
- `S3_singleton_structure.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
