#' @importFrom stats rnorm
#' @title Generate Dummy Multi-Omics Data
#' @description Generates synthetic, linked RNA-seq and Metabolomics datasets.
#' @param n_samples Number of samples (rows).
#' @param n_genes Number of genes (columns in X).
#' @param n_metabolites Number of metabolites (columns in Y).
#' @param n_linked Number of features linked by a hidden variable.
#' @return A list containing X (RNA-seq matrix), Y (Metabolomics matrix), and metadata.
#' @export
generate_dummy_omics <- function(n_samples = 50, n_genes = 1000, n_metabolites = 200, n_linked = 10) {
  
  # 1. Create a hidden biological signal (Z)
  Z <- stats::rnorm(n_samples) # Represents a latent factor (e.g., disease severity)
  
  # 2. Create RNA-seq matrix (X)
  X <- matrix(stats::rnorm(n_samples * n_genes), nrow = n_samples, ncol = n_genes)
  colnames(X) <- paste0("Gene_", 1:n_genes)
  
  # Inject signal into the first 'n_linked' genes
  for (i in 1:n_linked) {
    sign <- ifelse(i %% 2 == 0, 1, -1)
    X[, i] <- X[, i] + sign * 1.5 * Z + stats::rnorm(n_samples, 0, 0.5)
  }
  
  # 3. Create Metabolomics matrix (Y)
  Y <- matrix(stats::rnorm(n_samples * n_metabolites), nrow = n_samples, ncol = n_metabolites)
  colnames(Y) <- paste0("Met_", 1:n_metabolites)
  
  # Inject signal into the first 'n_linked' metabolites
  for (i in 1:n_linked) {
    sign <- ifelse(i %% 2 == 0, -1, 1)
    Y[, i] <- Y[, i] + sign * 1.8 * Z + stats::rnorm(n_samples, 0, 0.5)
  }
  
  # 4. Create Metadata
  metadata <- data.frame(
    SampleID = paste0("S", 1:n_samples),
    Group = factor(ifelse(Z > 0.5, "High_Z", "Low_Z")),
    Z_Score = Z
  )
  rownames(X) <- rownames(Y) <- metadata$SampleID
  
  # Standardize the data before output
  X <- scale(X)
  Y <- scale(Y)
  
  return(list(X = X, Y = Y, metadata = metadata))
}
