#' @title Align Multi-Omic Datasets
#' @description Ensures that X and Y matrices have matching samples in the exact same order.
#' @param X Matrix or data frame (Samples x Features).
#' @param Y Matrix or data frame (Samples x Features).
#' @return A list containing the aligned X and Y matrices.
#' @export
align_omics <- function(X, Y) {
  common_samples <- intersect(rownames(X), rownames(Y))
  
  if (length(common_samples) == 0) {
    stop("Error: No matching sample names found between X and Y. Check your rownames.")
  }
  
  message(paste("Successfully aligned", length(common_samples), "matching samples."))
  
  # Reorder and subset
  X_aligned <- as.matrix(X[common_samples, , drop = FALSE])
  Y_aligned <- as.matrix(Y[common_samples, , drop = FALSE])
  
  return(list(X = X_aligned, Y = Y_aligned))
}

#' @title Perform Sparse Canonical Correlation Analysis (sCCA)
#' @description Fits a sparse PLS model in canonical mode to identify shared variation.
#' @param X Normalized RNA-seq matrix (Samples x Features).
#' @param Y Normalized Metabolomics matrix (Samples x Features).
#' @param n_components Number of components.
#' @param penalty_X Sparsity for X (0 to 1, where 1 is most sparse).
#' @param penalty_Y Sparsity for Y (0 to 1, where 1 is most sparse).
#' @export
omic_scca <- function(X, Y, n_components = 2, penalty_X = 0.9, penalty_Y = 0.9) {
  
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package 'mixOmics' is required. Please install it.")
  }
  
  # Calculate number of variables to keep
  keep_x <- max(round(ncol(X) * (1 - penalty_X)), 1)
  keep_y <- max(round(ncol(Y) * (1 - penalty_Y)), 1)
  
  message(paste("Model Optimization: Keeping", keep_x, "genes and", keep_y, "metabolites."))
  
  scca_result <- mixOmics::spls(
    X, Y, 
    ncomp = n_components, 
    mode = "canonical", 
    keepX = rep(keep_x, n_components), 
    keepY = rep(keep_y, n_components)
  )
  
  result_list <- list(
    canonical_correlations = scca_result$expl_var, 
    loadings = list(X = scca_result$loadings$X, Y = scca_result$loadings$Y),
    variates = list(X = scca_result$variates$X, Y = scca_result$variates$Y),
    penalties = list(penalty_X = penalty_X, penalty_Y = penalty_Y)
  )
  
  class(result_list) <- "OmicNetR_sCCA"
  return(result_list)
}

#' @title Convert sCCA Loadings to Network Edges
#' @description Generates an edge list for network plotting from sCCA loadings.
#' @param scca_model The result object from omic_scca().
#' @param comp_select Which canonical component to use.
#' @param weight_threshold Minimum absolute product of weights to include an edge.
#' @export
scca_to_network <- function(scca_model, comp_select = 1, weight_threshold = 0.05) {
  
  
  u <- scca_model$loadings$X[, comp_select] 
  v <- scca_model$loadings$Y[, comp_select] 
  
  u_nonzero <- u[u != 0]
  v_nonzero <- v[v != 0]
  
  if(length(u_nonzero) == 0 | length(v_nonzero) == 0) {
    stop("No features selected. Try decreasing the penalty_X or penalty_Y.")
  }
  
  edges <- expand.grid(
    Gene = names(u_nonzero), 
    Metabolite = names(v_nonzero), 
    stringsAsFactors = FALSE
  )
  
  edges$Weight_Product <- u_nonzero[edges$Gene] * v_nonzero[edges$Metabolite]
  
  edges_filt <- edges[abs(edges$Weight_Product) >= weight_threshold, , drop = FALSE]
  
  
  edges_filt$Interaction_Type <- ifelse(edges_filt$Weight_Product > 0, "Positive", "Negative")
  
  return(edges_filt)
}
