# S15

Dinucleotide-shuffled sequence controls. (A) Jaccard of the top-5% target gene list of the real sequences against each of 20 shuffled replicates at three energy thresholds, with the shuffled-versus-shuffled baseline and the difference annotated. (B) Proportion of the original target genes recovered across the shuffled replicates.

## code
- `make_figR3_shuffle_rep242.py`

## data
- `genes/` — 63 files

## figure
- `S15_shuffle_thresholds.pdf`
- `S15_shuffle_thresholds.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
