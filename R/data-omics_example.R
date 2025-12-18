#' Example Multi-Omics Dataset
#'
#' A simulated multi-omics dataset included in the OmicNetR package.
#' This dataset is intended for demonstrating data alignment,
#' sparse CCA analysis, and network visualization functions.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{X}{A numeric matrix of gene expression values
#'   (samples in rows, genes in columns).}
#'   \item{Y}{A numeric matrix of metabolite abundances
#'   (samples in rows, metabolites in columns).}
#' }
#'
#' @details
#' The dataset is small by design and should not be used for
#' biological inference. It is provided solely for examples,
#' vignettes, and unit testing.
#'
#' @source
#' Simulated data generated within the OmicNetR package.
#'
#' @keywords datasets
"omics_example"
