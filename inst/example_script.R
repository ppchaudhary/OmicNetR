# ----------------------------------------------------------------------
# 1. Setup and Data Generation
# ----------------------------------------------------------------------

library(OmicNetR)
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

message("Step 1: Generating dummy omics data...")
omics_data <- generate_dummy_omics(
  n_samples = 60, 
  n_genes = 800, 
  n_metabolites = 150, 
  n_linked = 20  # Increased linked features for better visualization
)

X <- omics_data$X 
Y <- omics_data$Y 

# ----------------------------------------------------------------------
# 2. Sparse CCA Model Fitting (Relaxed Penalties)
# ----------------------------------------------------------------------

message("Step 2: Fitting sCCA model...")

# We lowered penalties from 0.95 to 0.70. 
# This ensures the model "picks" enough features to form a network.
scca_model <- omic_scca(
  X = X, 
  Y = Y, 
  n_components = 2,  
  penalty_X = 0.70,  
  penalty_Y = 0.70   
)

cat("\nCanonical Correlations (Explained Variance):\n")
print(scca_model$canonical_correlations)

# ----------------------------------------------------------------------
# 3. Network Generation (The "Handshake")
# ----------------------------------------------------------------------

message("Step 3: Generating network edge list...")

# We lowered the threshold to 0.01 to ensure we see the initial connections.
net_data_comp1 <- scca_to_network(
  scca_model, 
  comp_select = 1,      
  weight_threshold = 0.01 
)

# If the network is still empty, stop and warn
if(nrow(net_data_comp1) == 0) {
  stop("Network is empty. Try lowering the penalty_X/Y even further.")
}

# OPTIONAL: Keep only the Top 50 strongest edges for a very clean plot
net_data_comp1 <- net_data_comp1 %>%
  arrange(desc(abs(Weight_Product))) %>%
  head(50)

cat("\nPlotting the top", nrow(net_data_comp1), "strongest interactions.\n")

# ----------------------------------------------------------------------
# 4. Visualization
# ----------------------------------------------------------------------

message("Step 4: Rendering plots...")

# A. Bipartite Network Plot
bipartite_plot <- plot_bipartite_network(
  net_data_comp1,
  gene_color = "#1F77B4", 
  metabolite_color = "#FF7F0E",
  layout_type = "fr"
)

# Display the network
print(bipartite_plot)

# B. Feature Importance Circle Plot
pathway_circle_plot <- plot_pathway_circle(
  scca_model,
  top_features = 25, 
  pathway_db = "KEGG"
)

# Display the circle plot
print(pathway_circle_plot)

message("\nAnalysis workflow successfully completed!")