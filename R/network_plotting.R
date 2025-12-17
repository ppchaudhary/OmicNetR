#' @title Plot Bi-partite sCCA Weight Network
#' @description Optimized version using Base-R igraph engine to prevent memory exhaustion.
#' @param net_data The edge list data frame from scca_to_network().
#' @param gene_color Color for gene nodes.
#' @param metabolite_color Color for metabolite nodes.
#' @param layout_type igraph layout to use (default "fr").
#' @return A graph object (invisibly).
#' @export
plot_bipartite_network <- function(net_data, gene_color = "#1F77B4", metabolite_color = "#FF7F0E", layout_type = "fr") {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network plotting.")
  }
  
  # 1. Prune nodes: Only include nodes that are actually in the filtered edge list
  used_nodes <- unique(c(net_data$Gene, net_data$Metabolite))
  nodes <- data.frame(
    name = used_nodes,
    type = ifelse(used_nodes %in% net_data$Gene, "Gene", "Metabolite"),
    stringsAsFactors = FALSE
  )
  
  # 2. Create the igraph object
  g <- igraph::graph_from_data_frame(
    d = net_data[, c("Gene", "Metabolite")], 
    vertices = nodes, 
    directed = FALSE
  )
  
  # 3. Assign visual attributes directly to the igraph object
  # Vertex (Node) attributes
  igraph::V(g)$color <- ifelse(igraph::V(g)$type == "Gene", gene_color, metabolite_color)
  igraph::V(g)$size <- igraph::degree(g) * 1.5 + 4
  igraph::V(g)$label.cex <- 0.7
  igraph::V(g)$label.font <- 2
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$frame.color <- "white"
  
  # Edge attributes (Positive = Green, Negative = Red)
  igraph::E(g)$color <- ifelse(net_data$Interaction_Type == "Positive", "#2ca02c", "#d62728")
  # Scale widths for visibility
  igraph::E(g)$width <- abs(net_data$Weight_Product) * 8
  
  # 4. Use the Base R plotting engine (Extremely memory efficient)
  message("Rendering network using Base-R Graphics Engine...")
  
  # Prepare the layout
  l <- if(layout_type == "fr") {
    igraph::layout_with_fr(g, weights = abs(net_data$Weight_Product))
  } else {
    igraph::layout_as_bipartite(g)
  }
  
  # Execute the plot
  plot(g, 
       layout = l,
       vertex.label = igraph::V(g)$name,
       vertex.label.dist = 1.2,
       main = "OmicNetR: Bipartite Network",
       sub = paste("Top", nrow(net_data), "sCCA interactions"))
  
  # Add a small legend for the nodes
  legend("bottomleft", legend=c("Gene", "Metabolite"), 
         pt.bg=c(gene_color, metabolite_color), pch=21, pt.cex=2, bty="n", title="Layers")
  
  # Return the graph object invisibly to prevent double-printing
  return(invisible(g))
}

#' @title Canonical Loading Pathway Circle Plot
#' @description Visualizes top feature importance in a radial layout.
#' @param scca_model The result object from omic_scca().
#' @param top_features Number of most weighted features to map.
#' @param pathway_db Conceptual database name for labeling.
#' @return A ggplot2 object.
#' @export
plot_pathway_circle <- function(scca_model, top_features = 40, pathway_db = "KEGG") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  # Extract absolute loadings
  gene_loadings <- abs(scca_model$loadings$X[, 1])
  met_loadings <- abs(scca_model$loadings$Y[, 1])
  
  all_loadings <- c(gene_loadings, met_loadings)
  names(all_loadings) <- c(rownames(scca_model$loadings$X), rownames(scca_model$loadings$Y))
  
  # Select top features
  top_list <- sort(all_loadings, decreasing = TRUE)[1:min(top_features, length(all_loadings))]
  
  df <- data.frame(
    Feature = names(top_list),
    Loading = as.numeric(top_list),
    Type = ifelse(names(top_list) %in% rownames(scca_model$loadings$X), "Gene", "Metabolite"),
    stringsAsFactors = FALSE
  )
  
  # Circle Plot remains in ggplot2 as it is much less memory-intensive than networks
  circle_p <- ggplot2::ggplot(df, aes(x = reorder(Feature, Loading), y = Loading, fill = Type)) +
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