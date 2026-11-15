# In-house endometriosis cohort — count matrices

Processed tsRNA count matrices for the 46 peripheral-blood libraries of the
in-house endometriosis cohort (20 cases, 26 controls) analysed in the
HERACLES manuscript.

| File | Rows | Columns | Size |
|---|---|---|---|
| `fragment_counts.tsv` | 43,554 tsRNA fragments | 1 ID + 46 samples | 6.9 MB |
| `cluster_counts.tsv` | 16,931 clusters | 1 ID + 46 samples | 2.7 MB |
| `metadata.tsv` | 46 samples | `sample`, `group` | — |

## Format

Both matrices are tab-separated with integer read counts, one row per feature
and one column per library.

**`fragment_counts.tsv`** — the non-clustered (fragment-level) matrix. The first
column, `tsRNA_id`, is the HERACLES fragment identifier and encodes the parent
tRNA, the coordinates on that tRNA, the fragment class and the sequence:

```
tRNA-Ala-AGC-1-1:18-35:tRF-5:ggggtatagctcagtggt
 └ parent tRNA ┘ └ coords ┘ └class┘ └── sequence ──┘
```

**`cluster_counts.tsv`** — the clustered matrix produced by the default
directional clustering (abundance ratio 2.0, k-mer size 10, minimum length
ratio 0.94). Each row is a cluster, named after its representative fragment;
the value is the summed count of every fragment assigned to that cluster. The
two matrices therefore carry the same total: 40,135,349 reads.

`metadata.tsv` gives the group of each library (`Disease` or `Control`).

## How these were produced

Both matrices are direct outputs of the HERACLES pipeline (version 1.0.0) run
on the 46 collapsed FASTA files: `-v 0` Bowtie mapping against GtRNAdb hg38,
fragments of 18-50 nt, followed by the directional clustering module. No
filtering beyond the pipeline defaults was applied, so the fragment matrix
contains all 43,554 annotated fragments and the cluster matrix all 16,931
clusters, including the 10,089 singletons.

## Data availability

Raw FASTQ files are available through controlled access under **HRA019031**
and are not redistributed here. The matrices above are processed counts and
contain no personal identifiers; libraries are labelled with internal sample
codes only.
