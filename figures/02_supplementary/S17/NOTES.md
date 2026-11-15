# S17

Parameter grid. (A) Cluster count per setting against the default of 16,931, with the change annotated. (B) Singleton share and multi-member cluster count per setting.

## code
- `make_figR3_grid.py`

## data
- `TableS31_param_grid_topology.csv` — 1 KB

## figure
- `S17_grid_sensitivity.pdf`
- `S17_grid_sensitivity.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
