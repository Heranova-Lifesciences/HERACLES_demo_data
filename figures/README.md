# Figure generation package — HERACLES revised manuscript

One folder per figure. Each contains three parts:

- `code/` — the script(s) that generate the figure
- `data/` — the input files the script reads
- `figure/` — the delivered PDF and PNG

Input files larger than 1 MB are kept once in `00_shared_data/` instead of being
duplicated across figure folders; the `NOTES.md` in each folder states which
shared file it needs.

## Contents

| Figure | Script(s) | Main inputs |
|---|---|---|
| Figure1 | `make_fig1_vector.py` | none (vector drawing) |
| Figure2 | `compose_main_figures.py`, `make_mainfig_panels.py` | `TableS19_marker_benchmark.csv`, `TableS20_fig2C_pairs_stratified.csv`, `source_fig2_original.pdf` |
| Figure3 | `compose_figure3_refstyle.py` | `e29_KEGG_2021_Human.csv`, `KEGG_enrichment_3UTR.csv` |
| Figure4 | `compose_fig4_2x2.py`, `make_fig4ab_v4_panels.py`, `make_mainfig_fig4cd_v4.py` | `TableS19_marker_benchmark.csv`, `TableS26_comparator_marker_rerun.csv`, `Fig4_panelA-D.png` |
| S01 | `make_figS2_v2.py` | `member_vs_rep.json` |
| S02 | `make_figures.py` | `cluster_summary_detailed.csv` (shared) |
| S03 | `make_figures.py` | `cluster_summary_detailed.csv` (shared) |
| S04 | `make_figR1_permutation_v2.py` | `TableS16_permuted_label_controls.csv`, `TableS25b_permutation_two_workflows.csv` |
| S05 | `make_figures.py` | `per_cluster_overlap.csv` |
| S06 | `make_figS20_blacklist_origin_venn.py` | `three_way_summary.json`, `three_way_de_cluster.csv` (shared) |
| S07 | `make_figS14_similarity_hist.py` | `cluster_summary_detailed_46.csv`, `cluster_separation.csv` (shared) |
| S08 | `make_figS15.py` | `threshold_summary.json` |
| S09 | `make_fig4ab_v4_panels.py` | `TableS19_marker_benchmark.csv`, `TableS26_comparator_marker_rerun.csv` |
| S10 | `make_figS21_cross_parent.py` | `TableS14_cluster_members.tsv`, `cluster_summary_blocked.csv` (shared) |
| S11 | `make_figR3_2C_LOO.py` | `TableS20_fig2C_pairs_stratified.csv` |
| S12 | `make_figR3_2p2_fdr_scale.py` | `TableS25_permuted_label_controls.csv` |
| S13 | `make_mainfig_fig4cd_v4.py` | `TableS19_marker_benchmark.csv` |
| S14 | `make_figR3_2C_LOO.py` | `TableS22_LOO_marker_stability.csv` |
| S15 | `make_figR3_shuffle_rep242.py` | `genes/` (60 text files) |
| S16 | `make_figR3_KEGG_shuffle_heatmap.py` | `enrichr/` (63 CSV files) |
| S17 | `make_figR3_grid.py` | `TableS31_param_grid_topology.csv` |
| S18 | `make_figR3_3p3_mapping_v2.py` | `ablation_v2_report.json` |
| S19 | `make_figR3_3p4_labelaware_v4.py` | none (values hard-coded in the script) |

## How to run

The scripts read their inputs through paths relative to the project root, for
example `cluster_summary_detailed.csv`, `server_results/results/...` and
`output/...`. To regenerate a figure, place the files from that figure's `data/`
(and any shared file listed in its `NOTES.md`) at the paths the script expects,
then run it with Python 3 and the packages in `requirements.txt`.

Note that several scripts perform a compositing step: `compose_fig4_2x2.py` and
`compose_main_figures.py` paste existing panel images rather than drawing them,
so the panel PNGs/PDFs in `data/` are required as well as the plotting scripts.

## Environment

Python 3 with numpy, pandas, matplotlib, matplotlib-venn, pymupdf and Pillow
(see `requirements.txt`).

All data files were produced by the HERACLES pipeline and the reviewer-response
analysis scripts. Raw FASTQ files for the in-house cohort are available through
controlled access under HRA019031 and are not redistributed here.
