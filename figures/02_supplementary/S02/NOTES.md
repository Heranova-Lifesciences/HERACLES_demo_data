# S02

The read-length boundary. (A) Unique fragments by length, with the 18-nt annotation floor and the 50-nt cycle limit marked. (B) Fraction of fragments ending exactly at the annotated tRNA 3' terminus, falling to 0 of 1,235 at 50 nt. (C) Distance from that terminus for 43-49-nt and 50-nt fragments.

## code
- `make_figures.py`

## data
- `cluster_summary_detailed.csv` — 19.1 MB, in `00_shared_data/`

## figure
- `S2_read_length_boundary.pdf`
- `S2_read_length_boundary.png`

## How to run
The scripts read their inputs through paths relative to the project root (`cluster_summary_detailed.csv`, `server_results/...`, `output/...`). Copy the files out of `00_shared_data/` into the paths the script expects, then run it with Python 3 and the packages listed in `requirements.txt` at the top level of this package.
