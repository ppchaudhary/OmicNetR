#' @title Plot Bi-partite sCCA Weight Network
#' @importFrom igraph graph_from_data_frame V E degree layout_with_fr layout_as_bipartite
#' @importFrom grDevices png dev.off
#' @importFrom graphics legend plot
#' @export
plot_bipartite_network <- function(net_data, gene_color = "#1F77B4", metabolite_color = "#FF7F0E", layout_type = "fr") {
  
  used_nodes <- unique(c(net_data$Gene, net_data$Metabolite))
  nodes <- data.frame(
    name = used_nodes,
    type = ifelse(used_nodes %in% net_data$Gene, "Gene", "Metabolite"),
    stringsAsFactors = FALSE
  )
  
  g <- igraph::graph_from_data_frame(d = net_data, vertices = nodes, directed = FALSE)
  
  igraph::V(g)$color <- ifelse(igraph::V(g)$type == "Gene", gene_color, metabolite_color)
  igraph::V(g)$size <- igraph::degree(g) * 1.2 + 5
  igraph::E(g)$color <- ifelse(net_data$Interaction_Type == "Positive", "#2ca02c", "#d62728")
  igraph::E(g)$width <- abs(net_data$Weight_Product) * 10
  
  l <- if(layout_type == "fr") igraph::layout_with_fr(g) else igraph::layout_as_bipartite(g)
  
  # Internal plotting
  plot(g, layout = l, main = "OmicNetR Bipartite Network")
  graphics::legend("bottomleft", legend=c("Gene", "Metabolite"), 
                   fill=c(gene_color, metabolite_color), bty="n")
  
  return(invisible(g))
}

#' @title Canonical Loading Pathway Circle Plot
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot segments text
#' @importFrom stats reorder
#' @export
plot_pathway_circle <- function(scca_model, top_features = 40, pathway_db = "KEGG") {
  
  gene_loadings <- abs(scca_model$loadings$X[, 1])
  met_loadings <- abs(scca_model$loadings$Y[, 1])
  all_loadings <- c(gene_loadings, met_loadings)
  
  top_list <- sort(all_loadings, decreasing = TRUE)[1:min(top_features, length(all_loadings))]
  
  n <- length(top_list)
  theta <- seq(0, 2 * pi, length.out = n + 1)[-1]
  radius <- top_list / max(top_list)
  
  # Base R Plotting
  graphics::plot(0, 0, type = "n", xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), 
                 axes = FALSE, xlab = "", ylab = "", main = paste(pathway_db, "Feature Importance"))
  
  for (i in 1:n) {
    x <- radius[i] * cos(theta[i])
    y <- radius[i] * sin(theta[i])
    graphics::segments(0, 0, x, y, col = "gray80")
    graphics::text(1.1 * cos(theta[i]), 1.1 * sin(theta[i]), names(top_list)[i], cex = 0.6)
  }
  return(invisible(NULL))
}