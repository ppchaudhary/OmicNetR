library(OmicNetR)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# 1. Setup and Data Generation
message("Step 1: Generating dummy omics data...")
omics_data <- generate_dummy_omics(
  n_samples = 60, 
  n_genes = 800, 
  n_metabolites = 150, 
  n_linked = 20  
)

# Crucial: Keep the original matrices for the heatmap later
X <- omics_data$X 
Y <- omics_data$Y 

# 2. Sparse CCA Model Fitting
message("Step 2: Fitting sCCA model...")
scca_model <- omic_scca(X = X, Y = Y, n_components = 2, penalty_X = 0.70, penalty_Y = 0.70)

# 3. Network Generation (The "Handshake")
message("Step 3: Generating network edge list...")
net_data_comp1 <- scca_to_network(scca_model, comp_select = 1, weight_threshold = 0.01)

# Keep top 50 strongest edges for a clean plot using Base R
if(nrow(net_data_comp1) > 50) {
  net_data_comp1 <- net_data_comp1[order(abs(net_data_comp1$Weight_Product), decreasing = TRUE), ]
  net_data_comp1 <- net_data_comp1[1:50, ]
}

# ----------------------------------------------------------------------
# 4. Saving and Rendering Visuals to man/figures
# ----------------------------------------------------------------------
out_dir <- file.path(tempdir(), "OmicNetR_example_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message(paste0("Step 4: Saving plots to: ", out_dir))

# Create the directory if it doesn't exist (root of your package)

# A. Bipartite Network Plot
# We use png() because this function uses Base R/igraph plotting
png(file.path(out_dir, "network_plot.png"), width = 1000, height = 1000, res = 150)
plot_bipartite_network(
  net_data_comp1,
  gene_color = "#1F77B4", 
  metabolite_color = "#FF7F0E",
  layout_type = "fr"
)
dev.off()

# B. Feature Importance Circle Plot
# This returns a ggplot2 object, so we use ggsave
pathway_circle_plot <- plot_pathway_circle(
  scca_model,
  top_features = 30, 
  pathway_db = "KEGG"
)
ggplot2::ggsave(file.path(out_dir, "circle_plot.png"), plot = pathway_circle_plot, width = 7, height = 7, dpi = 300)

# C. Global Correlation Heatmap
# This uses the original matrices (X and Y) to calculate correlations 
# for the features selected by the sCCA model.
png(file.path(out_dir, "heatmap_plot.png"), width = 1000, height = 1000, res = 150)
plot_correlation_heatmap(
  scca_model = scca_model, 
  X = X, 
  Y = Y, 
  top_n = 25
)
dev.off()

message("\nAnalysis workflow successfully completed!")
message("\nFiles generated in:")
message(out_dir)
message("- network_plot.png")
message("- circle_plot.png")
message("- heatmap_plot.png")
