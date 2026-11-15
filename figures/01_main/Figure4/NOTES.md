# Figure4

2 x 2 composite. A and B are the marker-by-tool heatmaps for Xiong 2024 and Lu 2021; C and D are the fragment- and cluster-level log2 fold change and the -log10 nominal P / BH panel with its own legend.

## code
- `compose_fig4_2x2.py`
- `make_fig4ab_v4_panels.py`
- `make_mainfig_fig4cd_v4.py`

## data
- `Fig4_panelA.png` — 480 KB
- `Fig4_panelB.png` — 299 KB
- `Fig4_panelC.png` — 364 KB
- `Fig4_panelD.png` — 269 KB
- `TableS19_marker_benchmark.csv` — 4 KB
- `TableS26_comparator_marker_rerun.csv` — 8 KB

## figure
- `Figure4_revised_v3.pdf`
- `Figure4_revised_v3.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
