# S13

Figure 4C/4D redrawn: log2 fold change for the five Xiong and three Lu markers (C), and -log10 nominal P with BH as a separate series and its own legend (D).

## code
- `make_mainfig_fig4cd_v4.py`

## data
- `TableS19_marker_benchmark.csv` — 4 KB

## figure
- `S13_fig4cd_redrawn.pdf`
- `S13_fig4cd_redrawn.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
