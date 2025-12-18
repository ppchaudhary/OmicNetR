
knitr::opts_chunk\$set( collapse = TRUE, comment = “\#\>”, fig.path =
“man/figures/README-”, out.width = “100%” )

# Load the LOCAL package during knitting

devtools::load_all()

\##OmicNetR

\#OmicNetR is an R package for the integrated analysis of multi-omics
datasets using Sparse Canonical Correlation Analysis (sCCA).

\#Installation

\#You can install the development version of OmicNetR from GitHub:

# install.packages(“devtools”)

# devtools::install_github(“ppchaudhary/OmicNetR”)

\#Quick Start Example

\#This example demonstrates how to generate integrated networks and
importance plots.

library(OmicNetR)

# 1. Generate synthetic data

set.seed(123) data \<- generate_dummy_omics( n_samples = 60, n_genes =
800, n_metabolites = 150 )

# 2. Run sCCA

scca_model \<- omic_scca( data$X, data$Y, penalty_X = 0.7, penalty_Y =
0.7 )

# 3. Plot bipartite network

net_data \<- scca_to_network( scca_model, weight_threshold = 0.05 )
plot_bipartite_network(net_data)

# 4. Plot circular importance

plot_pathway_circle( scca_model, top_features = 30 )

Biological Interpretation

Nodes

Blue circles represent genes

Orange circles represent metabolites

Edges

Green lines indicate positive correlations

Red lines indicate negative correlations

Contact

Developed by Prem Prashant Chaudhary GitHub ID: ppchaudhary
