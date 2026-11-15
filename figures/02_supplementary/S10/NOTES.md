# S10

Parental tRNAs per cluster. (A) 16,252 clusters derive from a single parental tRNA and 679 (4.0%) contain fragments from more than one. (B) Effect of blocking cross-parent joins on cluster topology.

## code
- `make_figS21_cross_parent.py`

## data
- `TableS14_cluster_members.tsv` — 8.6 MB, in `00_shared_data/`
- `cluster_summary_blocked.csv` — 19.3 MB, in `00_shared_data/`

## figure
- `S10_cross_parent.pdf`
- `S10_cross_parent.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
