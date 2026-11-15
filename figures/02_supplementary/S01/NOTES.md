# S01

Composition of the 26,518 member-to-representative relationships (A) and the terminal offset between each member and its representative (B).

## code
- `make_figS2_v2.py`

## data
- `member_vs_rep.json` — 698 KB

## figure
- `S1_cluster_relationships.pdf`
- `S1_cluster_relationships.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
