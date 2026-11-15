# S07

Sequence similarity between cluster representatives. (A) Distribution of similarity between representatives of different clusters on the same parent tRNA, computed from interval overlap with no alignment step. (B) Whether each cluster is closer to its own representative than to the nearest representative of another cluster.

## code
- `make_figS14_similarity_hist.py`

## data
- `cluster_summary_detailed_46.csv` — 19.1 MB, in `00_shared_data/`
- `cluster_separation.csv` — 3.1 MB, in `00_shared_data/`

## figure
- `S7_representative_identity.pdf`
- `S7_representative_identity.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
