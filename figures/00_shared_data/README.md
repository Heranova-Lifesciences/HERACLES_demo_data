# 00_shared_data

Input files larger than 1 MB are stored here once, rather than being duplicated
inside the figure folders that use them.

| File | Size | Used by |
|---|---|---|
| `cluster_summary_detailed.csv` | 19 MB | S02, S03 (`make_figures.py`) |
| `cluster_summary_detailed_46.csv` | 19 MB | S07 (representative similarity) |
| `cluster_summary_blocked.csv` | 19 MB | S10 (cross-parent control) |
| `TableS14_cluster_members.tsv` | 8.2 MB | S10 |
| `cluster_separation.csv` | 3.0 MB | S07 |
| `three_way_de_cluster.csv` | 2.2 MB | S06 |

To regenerate a figure, copy the file back to the relative path the script
expects (for example `cluster_summary_detailed.csv` in the project root) and run
the script from there.
