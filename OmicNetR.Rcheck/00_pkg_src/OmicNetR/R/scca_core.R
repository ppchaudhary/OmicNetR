#' @title Align Multi-Omic Datasets
#' @description Ensures that X and Y matrices have matching samples in the exact same order.
#' @param X Matrix or data frame (Samples x Features).
#' @param Y Matrix or data frame (Samples x Features).
#' @return A list containing the aligned X and Y matrices.
#' @export
align_omics <- function(X, Y) {
  common_samples <- intersect(rownames(X), rownames(Y))
  if (length(common_samples) == 0) {
    stop("Error: No matching sample names found between X and Y.")
  }
  X_aligned <- as.matrix(X[common_samples, , drop = FALSE])
  Y_aligned <- as.matrix(Y[common_samples, , drop = FALSE])
  return(list(X = X_aligned, Y = Y_aligned))
}

#' @title Perform Sparse Canonical Correlation Analysis (sCCA)
#' @description Fits a sparse PLS model in canonical mode.
#' @param X Normalized RNA-seq matrix.
#' @param Y Normalized Metabolomics matrix.
#' @param n_components Number of components.
#' @param penalty_X Sparsity for X (0 to 1).
#' @param penalty_Y Sparsity for Y (0 to 1).
#' @importFrom mixOmics spls
#' @export
omic_scca <- function(X, Y, n_components = 2, penalty_X = 0.9, penalty_Y = 0.9) {
  keep_x <- max(round(ncol(X) * (1 - penalty_X)), 1)
  keep_y <- max(round(ncol(Y) * (1 - penalty_Y)), 1)
  
  scca_result <- mixOmics::spls(
    X, Y, ncomp = n_components, mode = "canonical",
    keepX = rep(keep_x, n_components), keepY = rep(keep_y, n_components)
  )
  
  res <- list(
    canonical_correlations = scca_result$expl_var,
    loadings = list(X = scca_result$loadings$X, Y = scca_result$loadings$Y),
    variates = list(X = scca_result$variates$X, Y = scca_result$variates$Y)
  )
  class(res) <- "OmicNetR_sCCA"
  return(res)
}

#' @title Convert sCCA Loadings to Network Edges
#' @description Generates an edge list using Base R to avoid global variable notes.
#' @param scca_model The result object from omic_scca().
#' @param comp_select Which canonical component to use.
#' @param weight_threshold Minimum absolute product of weights to include.
#' @return A data frame of edges.
#' @export
scca_to_network <- function(scca_model, comp_select = 1, weight_threshold = 0.05) {
  u <- scca_model$loadings$X[, comp_select]
  v <- scca_model$loadings$Y[, comp_select]
  
  u_nonzero <- u[u != 0]
  v_nonzero <- v[v != 0]
  
  if (length(u_nonzero) == 0 || length(v_nonzero) == 0) {
    stop("No features selected. Try decreasing the penalty.")
  }
  
  edges <- expand.grid(
    Gene = names(u_nonzero),
    Metabolite = names(v_nonzero),
    stringsAsFactors = FALSE
  )
  
  # FIX: Use [[]] to prevent "no visible binding for global variable" NOTE
  edges[["Weight_Product"]] <- u_nonzero[edges[["Gene"]]] * v_nonzero[edges[["Metabolite"]]]
  
  # FIX: Base R subsetting to remove 'dplyr' and pipe dependency
  edges_filt <- edges[abs(edges[["Weight_Product"]]) >= weight_threshold, ]
  
  if (nrow(edges_filt) > 0) {
    edges_filt[["Interaction_Type"]] <- ifelse(edges_filt[["Weight_Product"]] > 0, "Positive", "Negative")
  }
  
  return(edges_filt)
}