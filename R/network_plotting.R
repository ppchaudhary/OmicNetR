# ----------------------------------------------------------------------
# Plot Bi-partite sCCA Weight Network
# ----------------------------------------------------------------------

#' @title Plot Bi-partite sCCA Weight Network
#' @description Optimized version using Base-R igraph engine to prevent memory exhaustion.
#' @param net_data The edge list data frame from scca_to_network().
#' @param gene_color Color for gene nodes.
#' @param metabolite_color Color for metabolite nodes.
#' @param layout_type igraph layout to use (default "fr").
#' @return A graph object (invisibly).
#' @importFrom igraph graph_from_data_frame V E degree layout_with_fr layout_as_bipartite
#' @importFrom graphics legend plot
#' @export
plot_bipartite_network <- function(net_data, gene_color = "#1F77B4", metabolite_color = "#FF7F0E", layout_type = "fr") {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network plotting.")
  }
  
  used_nodes <- unique(c(net_data$Gene, net_data$Metabolite))
  nodes <- data.frame(
    name = used_nodes,
    type = ifelse(used_nodes %in% net_data$Gene, "Gene", "Metabolite"),
    stringsAsFactors = FALSE
  )
  
  g <- igraph::graph_from_data_frame(
    d = net_data[, c("Gene", "Metabolite")], 
    vertices = nodes, 
    directed = FALSE
  )
  
  igraph::V(g)$color <- ifelse(igraph::V(g)$type == "Gene", gene_color, metabolite_color)
  igraph::V(g)$size <- igraph::degree(g) * 1.5 + 4
  igraph::V(g)$label.cex <- 0.7
  igraph::V(g)$label.font <- 2
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$frame.color <- "white"
  
  igraph::E(g)$color <- ifelse(net_data$Interaction_Type == "Positive", "#2ca02c", "#d62728")
  igraph::E(g)$width <- abs(net_data$Weight_Product) * 8
  
  message("Rendering network using Base-R Graphics Engine...")
  
  l <- if(layout_type == "fr") {
    igraph::layout_with_fr(g, weights = abs(net_data$Weight_Product))
  } else {
    igraph::layout_as_bipartite(g)
  }
  
  graphics::plot(g, 
                 layout = l,
                 vertex.label = igraph::V(g)$name,
                 vertex.label.dist = 1.2,
                 main = "OmicNetR: Bipartite Network",
                 sub = paste("Top", nrow(net_data), "sCCA interactions"))
  
  graphics::legend("bottomleft", legend=c("Gene", "Metabolite"), 
                   pt.bg=c(gene_color, metabolite_color), pch=21, pt.cex=2, bty="n", title="Layers")
  
  return(invisible(g))
}

# ----------------------------------------------------------------------
# Canonical Loading Pathway Circle Plot
# ----------------------------------------------------------------------

#' @title Canonical Loading Pathway Circle Plot
#' @description Visualizes top feature importance in a radial layout.
#' @param scca_model The result object from omic_scca().
#' @param top_features Number of most weighted features to map.
#' @param pathway_db Conceptual database name for labeling.
#' @return A ggplot2 object.
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar theme_minimal scale_fill_manual labs theme element_text
#' @importFrom stats reorder
#' @export
plot_pathway_circle <- function(scca_model, top_features = 40, pathway_db = "KEGG") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  gene_loadings <- abs(scca_model$loadings$X[, 1])
  met_loadings <- abs(scca_model$loadings$Y[, 1])
  all_loadings <- c(gene_loadings, met_loadings)
  names(all_loadings) <- c(rownames(scca_model$loadings$X), rownames(scca_model$loadings$Y))
  
  top_list <- sort(all_loadings, decreasing = TRUE)[1:min(top_features, length(all_loadings))]
  
  df <- data.frame(
    Feature = names(top_list),
    Loading = as.numeric(top_list),
    Type = ifelse(names(top_list) %in% rownames(scca_model$loadings$X), "Gene", "Metabolite"),
    stringsAsFactors = FALSE
  )
  
  circle_p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(Feature, Loading), y = Loading, fill = Type)) +
    ggplot2::geom_bar(stat = "identity", width = 0.8) +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = c("Gene" = "#1F77B4", "Metabolite" = "#FF7F0E")) +
    ggplot2::labs(
      title = "OmicNetR: Feature Importance",
      subtitle = paste("Top", top_features, "features (Abs. sCCA Loadings)"),
      x = "", y = ""
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7))
  
  return(circle_p)
}

# ----------------------------------------------------------------------
# Global Gene-Metabolite Correlation Heatmap
# ----------------------------------------------------------------------

#' @title Global Gene-Metabolite Correlation Heatmap
#' @description Visualizes the correlation matrix with a gradient color scale.
#' @param scca_model The result object from omic_scca().
#' @param X Aligned RNA-seq matrix.
#' @param Y Aligned Metabolomics matrix.
#' @param top_n Number of top features from each omic to include.
#' @importFrom graphics image axis layout par plot.new plot.window rect text rasterImage
#' @importFrom grDevices colorRampPalette as.raster
#' @importFrom stats cor
#' @export
plot_correlation_heatmap <- function(scca_model, X, Y, top_n = 20) {
  
  # 1. Identify top features
  u_top <- names(sort(abs(scca_model$loadings$X[,1]), decreasing = TRUE)[1:min(top_n, nrow(scca_model$loadings$X))])
  v_top <- names(sort(abs(scca_model$loadings$Y[,1]), decreasing = TRUE)[1:min(top_n, nrow(scca_model$loadings$Y))])
  
  cor_mat <- stats::cor(X[, u_top, drop = FALSE], Y[, v_top, drop = FALSE])
  
  # 2. Setup colors (Blue to Red)
  col_pal <- grDevices::colorRampPalette(c("#1F77B4", "white", "#D62728"))(100)
  
  # 3. Layout: Split screen for Heatmap (left) and Color Scale (right)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1.2))
  graphics::par(mar = c(8, 8, 4, 1))
  
  # 4. Draw Heatmap (zlim ensures white = 0 correlation)
  graphics::image(1:nrow(cor_mat), 1:ncol(cor_mat), cor_mat, 
                  col = col_pal, axes = FALSE, xlab = "", ylab = "",
                  main = "OmicNetR: Feature Correlation Landscape",
                  zlim = c(-1, 1))
  
  graphics::axis(1, at = 1:nrow(cor_mat), labels = u_top, las = 2, cex.axis = 0.7)
  graphics::axis(2, at = 1:ncol(cor_mat), labels = v_top, las = 2, cex.axis = 0.7)
  
  # 5. Draw Color Gradient Legend
  graphics::par(mar = c(8, 1, 4, 3))
  legend_image <- grDevices::as.raster(matrix(rev(col_pal), ncol = 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(-1, 1))
  graphics::rasterImage(legend_image, 0, -1, 0.4, 1)
  
  # Add labels to the scale
  graphics::axis(4, at = c(-1, 0, 1), labels = c("-1", "0", "1"), las = 1, cex.axis = 0.8)
  graphics::text(1.3, 0, labels = "Correlation (r)", srt = 270, xpd = TRUE, font = 2)
  
  # Reset layout
  graphics::layout(1)
  return(invisible(cor_mat))
}