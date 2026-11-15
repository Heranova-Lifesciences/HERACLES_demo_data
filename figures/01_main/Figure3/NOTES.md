# Figure3

Panel A is the workflow diagram. Panels B and C are KEGG over-representation bubble plots for CDS and 3'UTR target genes: the x axis is the Enrichr combined score, the colour bar is -log10 nominal P, and bubble size is the gene count.

## code
- `compose_figure3_refstyle.py`

## data
- `KEGG_enrichment_3UTR.csv` — 51 KB
- `e29_KEGG_2021_Human.csv` — 11 KB

## figure
- `Figure3_revised.pdf`
- `Figure3_revised.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
