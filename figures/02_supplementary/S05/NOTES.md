# S05

Representative versus union-of-members target sets. (A) Jaccard index between the representative's target set and the union of its members' target sets. (B) Per-cluster distribution of that Jaccard index.

## code
- `make_figures.py`

## data
- `per_cluster_overlap.csv` — 21 KB

## figure
- `S5_rep_vs_members.pdf`
- `S5_rep_vs_members.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
