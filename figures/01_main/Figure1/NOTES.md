# Figure1

Workflow schematic of HERACLES, drawn as vector graphics with pymupdf. No external data input.

## code
- `make_fig1_vector.py`

## data
- (none: values are hard-coded in the script)

## figure
- `Figure1_vector.pdf`
- `Figure1_vector.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
