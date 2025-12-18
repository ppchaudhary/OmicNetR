#!/usr/bin/env Rscript
# tools/render_readme.R
#
# Purpose:
#   1) Generate README figures into man/figures/
#   2) Render README.Rmd -> README.md
#
# Run from the package root:
#   Rscript tools/render_readme.R

stopifnot(file.exists("DESCRIPTION"))
dir.create(file.path("man", "figures"), recursive = TRUE, showWarnings = FALSE)

# Install/load deps needed for plotting in the README figures.
# (Keep this script "developer-only" — it's not run during R CMD check.)
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("OmicNetR", quietly = TRUE)) {
  message("Installing OmicNetR from local source...")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install(upgrade = "never")
}

library(OmicNetR)
library(ggplot2)

set.seed(123)

omics_data <- generate_dummy_omics(
  n_samples = 60,
  n_genes = 800,
  n_metabolites = 150,
  n_linked = 20
)

X <- omics_data$X
Y <- omics_data$Y

aligned <- align_omics(X, Y)

scca_model <- omic_scca(
  X = aligned$X,
  Y = aligned$Y,
  n_components = 2,
  penalty_X = 0.70,
  penalty_Y = 0.70
)

net_data_comp1 <- scca_to_network(scca_model, comp_select = 1, weight_threshold = 0.01)

# Keep top 50 strongest edges for a clean plot
if (nrow(net_data_comp1) > 50) {
  net_data_comp1 <- net_data_comp1[order(abs(net_data_comp1$Weight_Product), decreasing = TRUE), ]
  net_data_comp1 <- net_data_comp1[1:50, ]
}

# A) Bipartite network plot (base R)
png(file.path("man", "figures", "network_plot.png"), width = 1000, height = 1000, res = 150)
plot_bipartite_network(
  net_data_comp1,
  gene_color = "#1F77B4",
  metabolite_color = "#FF7F0E",
  layout_type = "fr"
)
dev.off()

# B) Circle plot (ggplot2)
circle_plot <- plot_pathway_circle(scca_model, top_features = 30, pathway_db = "KEGG")
ggplot2::ggsave(
  filename = file.path("man", "figures", "circle_plot.png"),
  plot = circle_plot,
  width = 7,
  height = 7,
  dpi = 300
)

# C) Heatmap (base R)
png(file.path("man", "figures", "heatmap_plot.png"), width = 1200, height = 1200, res = 150)
OmicNetR:::plot_correlation_heatmap(scca_model = scca_model, X = aligned$X, Y = aligned$Y, top_n = 25)
dev.off()

# Render README
if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown")
rmarkdown::render("README.Rmd", output_format = "github_document")

message("Done! Updated README.md and wrote images to man/figures/.")
