# Reproducibility materials for the OmicNetR manuscript

This directory contains the exact analysis scripts used for the manuscript:

**OmicNetR: An R Package for Interpretable Network-Based Integration of Multi-Omics Data Using Sparse Canonical Correlation Analysis**

## Software environment

- R 4.5.1 (2025-06-13)
- OmicNetR 0.1.1
- mixOmics 6.34.0
- ggplot2 4.0.3
- igraph 2.3.3
- patchwork 1.3.2
- dplyr 1.2.1
- tidyr 1.3.2
- readr 2.2.0

The benchmark script writes `sessionInfo.txt` for the complete computational environment.

## Scripts

### 1. Simulation benchmark and Figure 2

Run:

```r
source("scripts/01_simulation_benchmark_and_figure2.R")
```

This script performs the OmicNetR-versus-mixOmics benchmark using 1,200 independently simulated datasets:
- n = 30, 60, 120
- signal strengths = 0.5, 1.0, 1.5, 2.0
- 100 replicates per condition
- 800 genes and 150 metabolites
- 20 true signal genes and 20 true signal metabolites
- 70% training / 30% held-out testing
- 2 canonical components
- master seed = 20260824

It generates the files underlying Figure 2, Supplementary Figure S1, and Supplementary Tables S1-S4, including benchmark summaries, paired F1 comparisons, held-out correlation comparisons, bootstrap stability summaries, runtime plots, and `sessionInfo.txt`.

### 2. MC903 real-data analysis and Figure 4

The script expects the following processed input files in the working directory:

```text
RNA_seq_normalised.csv
Metabolomics_raw.csv
```

Run:

```r
source("scripts/02_MC903_real_data_and_figure4.R")
```

The script aligns the matched transcriptomic/metabolomic samples, preprocesses the metabolomics matrix, runs OmicNetR sCCA, generates the individual outputs used for Figure 3, performs post-hoc OA/DIL orientation of Component 1, computes the group-comparison statistics, creates the directional networks, and compiles publication-quality Figure 4.

Primary output directory:

```text
OmicNetR_MC903_FINAL/
```

Important outputs include:
- `Figure3A_Integrated_VIP_circle.png`
- `Figure3B_Heatmap_integrated.png`
- `Integrated_top20_network_edges.csv`
- `Component1_sample_scores_oriented.csv`
- `Component1_group_score_means.csv`
- `Component1_group_association_statistics.csv`
- `Component1_orientation_information.csv`
- `Figure_4_FINAL_PUBLICATION_ALIGNED.png`
- `Figure_4_FINAL_PUBLICATION_ALIGNED.pdf`

### 3. Compile publication-quality Figure 3

Run the MC903 real-data script first.

Then run:

```r
source("scripts/03_compile_figure3.R")
```

This uses the Figure 3 component outputs from `OmicNetR_MC903_FINAL/` and generates:
- `Figure_3_FINAL_PUBLICATION.png`
- `Figure_3_FINAL_PUBLICATION.pdf`

## Recommended run order

```text
1. scripts/01_simulation_benchmark_and_figure2.R
2. scripts/02_MC903_real_data_and_figure4.R
3. scripts/03_compile_figure3.R
```

## Real-data availability

Only place the processed MC903 input matrices in the public repository if they are approved for public release. If they cannot be made public, do not upload them; keep the `data/README.md` file and provide the appropriate data-access statement in the manuscript.

## Interpretation note

The Figure 3 network uses signed loading-product scores for association direction and absolute loading-product magnitude for edge strength. Figure 4 instead visualizes opposite directions of a post-hoc oriented latent component and should not be interpreted as direct treatment effects or causal molecular interactions.
