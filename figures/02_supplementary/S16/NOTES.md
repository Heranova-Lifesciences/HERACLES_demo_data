# S16

Pathway-level view of the same shuffled control: pairwise Jaccard of the KEGG pathway sets of the real and 20 shuffled replicates at e = -29, using the top-20 terms per library. A nominal-count variant is available from the same script with `--metric nom`.

## code
- `make_figR3_KEGG_shuffle_heatmap.py`

## data
- `enrichr/` — 63 files

## figure
- `S16_pathway_jaccard.pdf`
- `S16_pathway_jaccard.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
