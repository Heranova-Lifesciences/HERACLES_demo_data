# S06

Origin of the released CDS blacklist: three-way Venn of the top-10% predicted tsRNA-target connectivity sets from the in-house, Xiong 2024 and Lu 2021 cohorts. The 1,377-gene three-way intersection is the released blacklist. Region counts are hard-coded in the script; the underlying three-way summary and DE table are in `data/`.

## code
- `make_figS20_blacklist_origin_venn.py`

## data
- `three_way_summary.json` — 1 KB
- `three_way_de_cluster.csv` — 2.2 MB, in `00_shared_data/`

## figure
- `S6_blacklist_origin_venn.pdf`
- `S6_blacklist_origin_venn.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
