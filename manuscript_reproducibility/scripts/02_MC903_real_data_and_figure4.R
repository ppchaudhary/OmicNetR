
## ============================================================
## OmicNetR real-data workflow: MC903 atopic dermatitis
## Corrected directional interpretation + publication Figure 4
## ============================================================

## ----------------------------
## 0. INSTALL / LOAD PACKAGES
## ----------------------------
cran_pkgs <- c(
  "OmicNetR",
  "ggplot2",
  "patchwork",
  "igraph",
  "png",
  "grid"
)

missing_pkgs <- cran_pkgs[
  !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(OmicNetR)
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(png)
  library(grid)
})

## ----------------------------
## 1. USER SETTINGS
## ----------------------------
RNA_FILE <- "RNA_seq_normalised.csv"
MET_FILE <- "Metabolomics_raw.csv"

OUTDIR <- "OmicNetR_MC903_FINAL"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

N_COMPONENTS <- 2
PENALTY_X <- 0.70
PENALTY_Y <- 0.70

TOP_NETWORK_EDGES <- 20
TOP_HEATMAP_FEATURES <- 25
TOP_VIP_FEATURES <- 30

set.seed(20260825)

## ============================================================
## 2. RNA-SEQ INPUT
## ============================================================
rna <- read.csv(
  RNA_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_ids <- as.character(rna$SampleID)
rna$SampleID <- NULL

X0 <- as.matrix(rna)
storage.mode(X0) <- "numeric"

# X = samples x genes
X <- t(X0)

# Sample names are the RNA-seq column names (OA1..OA5, DIL1..DIL5)
rownames(X) <- colnames(rna)

# Clean gene names
gene_ids <- trimws(gene_ids)
gene_ids[
  is.na(gene_ids) | gene_ids == ""
] <- paste0(
  "gene_",
  which(is.na(gene_ids) | gene_ids == "")
)

colnames(X) <- make.unique(gene_ids)

# Scale gene features
X <- scale(X)
X <- as.matrix(X)

stopifnot(!is.null(colnames(X)))
stopifnot(sum(duplicated(colnames(X))) == 0)

## ============================================================
## 3. METABOLOMICS INPUT
## ============================================================
met <- read.csv(
  MET_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_cols <- grep(
  "^(OA|DIL)[0-9]+$",
  colnames(met),
  value = TRUE
)

stopifnot(length(sample_cols) >= 10)

feature_names <- as.character(met[[1]])

Yraw <- t(
  as.matrix(
    met[, sample_cols, drop = FALSE]
  )
)

Yraw <- apply(
  Yraw,
  2,
  function(z) {
    as.numeric(
      gsub(",", "", z)
    )
  }
)

storage.mode(Yraw) <- "numeric"

rownames(Yraw) <- sample_cols

feature_names <- trimws(feature_names)
feature_names[
  is.na(feature_names) | feature_names == ""
] <- paste0(
  "met_",
  which(is.na(feature_names) | feature_names == "")
)

colnames(Yraw) <- make.unique(feature_names)

## ============================================================
## 4. METABOLOMICS PREPROCESSING
## ============================================================

# Treat zeros as missing
Yraw[Yraw == 0] <- NA

# Remove completely missing metabolite features
Yraw <- Yraw[
  ,
  colSums(!is.na(Yraw)) > 0,
  drop = FALSE
]

# Half-minimum imputation within each feature
for (j in seq_len(ncol(Yraw))) {

  if (anyNA(Yraw[, j])) {

    mn <- suppressWarnings(
      min(
        Yraw[, j],
        na.rm = TRUE
      )
    )

    if (is.finite(mn)) {
      Yraw[
        is.na(Yraw[, j]),
        j
      ] <- 0.5 * mn
    } else {
      Yraw[
        is.na(Yraw[, j]),
        j
      ] <- 0
    }
  }
}

# TIC-like sample normalization
sample_totals <- rowSums(Yraw)

Ynorm <- Yraw /
  sample_totals *
  median(sample_totals)

# Log transform + autoscale
Y <- scale(
  log2(Ynorm + 1)
)

Y <- as.matrix(Y)

rownames(Y) <- rownames(Yraw)
colnames(Y) <- make.unique(
  trimws(colnames(Yraw))
)

stopifnot(!is.null(colnames(Y)))
stopifnot(sum(duplicated(colnames(Y))) == 0)

## ============================================================
## 5. ALIGN MATCHED SAMPLES
## ============================================================
common <- intersect(
  rownames(X),
  rownames(Y)
)

X <- X[
  common,
  ,
  drop = FALSE
]

Y <- Y[
  common,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(X),
    rownames(Y)
  )
)

X_used <- X
Y_used <- Y

group <- ifelse(
  grepl("^OA", rownames(X_used)),
  "OA",
  "DIL"
)

metadata <- data.frame(
  Sample = rownames(X_used),
  Group = group,
  stringsAsFactors = FALSE
)

write.csv(
  metadata,
  file.path(OUTDIR, "MC903_sample_metadata.csv"),
  row.names = FALSE
)

cat("\nMatched samples:\n")
print(metadata)

## ============================================================
## 6. RUN OmicNetR sCCA
## ============================================================
scca_model <- omic_scca(
  X = X_used,
  Y = Y_used,
  n_components = N_COMPONENTS,
  penalty_X = PENALTY_X,
  penalty_Y = PENALTY_Y
)

saveRDS(
  scca_model,
  file.path(OUTDIR, "MC903_OmicNetR_scca_model.rds")
)

## ============================================================
## 7. INTEGRATED NETWORK
## ============================================================
net0 <- scca_to_network(
  scca_model,
  comp_select = 1,
  weight_threshold = 0
)

# Preserve raw loading-product score
net0$Weight_Product_raw <- net0$Weight_Product

# Top integrated edges for Figure 3 / exploratory output
net_top <- net0[
  order(
    abs(net0$Weight_Product_raw),
    decreasing = TRUE
  ),
]

net_top <- head(
  net_top,
  TOP_NETWORK_EDGES
)

if (nrow(net_top) > 1) {
  net_top$Weight_Product <- as.numeric(
    scale(
      abs(net_top$Weight_Product_raw)
    )
  )
}

png(
  file.path(
    OUTDIR,
    "Figure3C_network_integrated_top20.png"
  ),
  width = 1800,
  height = 1800,
  res = 250
)

plot_bipartite_network(
  net_top,
  gene_color = "#1F77B4",
  metabolite_color = "#FF7F0E",
  layout_type = "fr"
)

dev.off()

write.csv(
  net_top,
  file.path(
    OUTDIR,
    "Integrated_top20_network_edges.csv"
  ),
  row.names = FALSE
)

## ============================================================
## 8. INTEGRATED VIP CIRCLE
## ============================================================
save_vip_circle <- function(
    model,
    prefix,
    top_features = 30,
    pathway_db = "KEGG") {

  ok <- try({

    p <- plot_pathway_circle(
      model,
      top_features = top_features,
      pathway_db = pathway_db
    )

    ggsave(
      file.path(
        OUTDIR,
        paste0(prefix, "_VIP_circle.png")
      ),
      plot = p,
      width = 7,
      height = 7,
      dpi = 600,
      bg = "white"
    )

    TRUE

  }, silent = TRUE)

  if (!isTRUE(ok)) {

    message(
      prefix,
      ": KEGG-based VIP circle failed; retrying with pathway_db='NONE'."
    )

    ok2 <- try({

      p <- plot_pathway_circle(
        model,
        top_features = top_features,
        pathway_db = "NONE"
      )

      ggsave(
        file.path(
          OUTDIR,
          paste0(prefix, "_VIP_circle.png")
        ),
        plot = p,
        width = 7,
        height = 7,
        dpi = 600,
        bg = "white"
      )

      TRUE

    }, silent = TRUE)

    if (!isTRUE(ok2)) {
      message(
        prefix,
        ": VIP circle failed."
      )
    }
  }
}

save_vip_circle(
  scca_model,
  "Figure3A_Integrated",
  top_features = TOP_VIP_FEATURES,
  pathway_db = "KEGG"
)

## ============================================================
## 9. INTEGRATED HEATMAP
## ============================================================
save_heatmap <- function(
    model,
    Xmat,
    Ymat,
    out_png,
    top_n = 25) {

  png(
    file.path(
      OUTDIR,
      out_png
    ),
    width = 1800,
    height = 1800,
    res = 250
  )

  plot_correlation_heatmap(
    scca_model = model,
    X = Xmat,
    Y = Ymat,
    top_n = top_n
  )

  dev.off()
}

save_heatmap(
  scca_model,
  X_used,
  Y_used,
  "Figure3B_Heatmap_integrated.png",
  top_n = TOP_HEATMAP_FEATURES
)

## ============================================================
## 10. CORRECTED TREATMENT-ASSOCIATED COMPONENT ORIENTATION
## ============================================================

# Component 1 canonical sample scores
scores_X <- as.numeric(
  scca_model$variates$X[, 1]
)

scores_Y <- as.numeric(
  scca_model$variates$Y[, 1]
)

scores_df <- data.frame(
  Sample = rownames(X_used),
  Score_X = scores_X,
  Score_Y = scores_Y,
  Group = group,
  stringsAsFactors = FALSE
)

cat("\nRaw Component 1 score means before orientation:\n")

score_summary_raw <- aggregate(
  cbind(
    Score_X,
    Score_Y
  ) ~ Group,
  data = scores_df,
  FUN = mean
)

print(score_summary_raw)

# Original loadings
gene_load <- scca_model$loadings$X[, 1]
met_load <- scca_model$loadings$Y[, 1]

# Determine direction from the transcriptomic canonical score.
# Canonical signs are mathematically arbitrary.
mean_OA_X <- mean(
  scores_df$Score_X[
    scores_df$Group == "OA"
  ]
)

mean_DIL_X <- mean(
  scores_df$Score_X[
    scores_df$Group == "DIL"
  ]
)

ORIENTATION_FLIPPED <- FALSE

if (mean_OA_X < mean_DIL_X) {

  ORIENTATION_FLIPPED <- TRUE

  gene_load <- -gene_load
  met_load <- -met_load

  scores_X <- -scores_X
  scores_Y <- -scores_Y

  message(
    "Component 1 sign flipped so that the positive direction corresponds to OA."
  )

} else {

  message(
    "Component 1 orientation retained; positive direction already corresponds to OA."
  )
}

# Store oriented scores
scores_df$Score_X_oriented <- scores_X
scores_df$Score_Y_oriented <- scores_Y

score_summary_oriented <- aggregate(
  cbind(
    Score_X_oriented,
    Score_Y_oriented
  ) ~ Group,
  data = scores_df,
  FUN = mean
)

cat("\nOriented Component 1 score means:\n")
print(score_summary_oriented)

write.csv(
  scores_df,
  file.path(
    OUTDIR,
    "Component1_sample_scores_oriented.csv"
  ),
  row.names = FALSE
)

write.csv(
  score_summary_oriented,
  file.path(
    OUTDIR,
    "Component1_group_score_means.csv"
  ),
  row.names = FALSE
)

## ============================================================
## 11. GROUP ASSOCIATION TESTS FOR COMPONENT 1
## ============================================================

# Welch t-test
t_test_X <- t.test(
  Score_X_oriented ~ Group,
  data = scores_df
)

t_test_Y <- t.test(
  Score_Y_oriented ~ Group,
  data = scores_df
)

# Wilcoxon rank-sum test as a small-sample sensitivity analysis
wilcox_X <- wilcox.test(
  Score_X_oriented ~ Group,
  data = scores_df,
  exact = FALSE
)

wilcox_Y <- wilcox.test(
  Score_Y_oriented ~ Group,
  data = scores_df,
  exact = FALSE
)

# Cohen's d helper
cohens_d_two_groups <- function(values, groups, group1 = "OA", group2 = "DIL") {

  x1 <- values[groups == group1]
  x2 <- values[groups == group2]

  n1 <- length(x1)
  n2 <- length(x2)

  s1 <- sd(x1)
  s2 <- sd(x2)

  pooled_sd <- sqrt(
    (
      (n1 - 1) * s1^2 +
      (n2 - 1) * s2^2
    ) /
      (n1 + n2 - 2)
  )

  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }

  (
    mean(x1) - mean(x2)
  ) / pooled_sd
}

d_X <- cohens_d_two_groups(
  scores_df$Score_X_oriented,
  scores_df$Group
)

d_Y <- cohens_d_two_groups(
  scores_df$Score_Y_oriented,
  scores_df$Group
)

score_stats <- data.frame(
  Layer = c("Transcriptomic canonical score", "Metabolomic canonical score"),
  OA_mean = c(
    mean(scores_df$Score_X_oriented[scores_df$Group == "OA"]),
    mean(scores_df$Score_Y_oriented[scores_df$Group == "OA"])
  ),
  DIL_mean = c(
    mean(scores_df$Score_X_oriented[scores_df$Group == "DIL"]),
    mean(scores_df$Score_Y_oriented[scores_df$Group == "DIL"])
  ),
  Welch_t_p = c(
    t_test_X$p.value,
    t_test_Y$p.value
  ),
  Wilcoxon_p = c(
    wilcox_X$p.value,
    wilcox_Y$p.value
  ),
  Cohens_d_OA_minus_DIL = c(
    d_X,
    d_Y
  )
)

print(score_stats)

write.csv(
  score_stats,
  file.path(
    OUTDIR,
    "Component1_group_association_statistics.csv"
  ),
  row.names = FALSE
)

## ============================================================
## 12. FIGURE 4A: COMPONENT 1 SAMPLE SCORE PLOT
## ============================================================

# Use transcriptomic canonical scores as the orientation anchor
p_scores <- ggplot(
  scores_df,
  aes(
    x = Group,
    y = Score_X_oriented,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 3.2,
    alpha = 0.9
  ) +
  theme_classic(
    base_size = 15
  ) +
  guides(
    fill = "none"
  ) +
  labs(
    title = "Treatment association with cross-omics Component 1",
    x = NULL,
    y = "Oriented Component 1 canonical score"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    )
  )

ggsave(
  file.path(
    OUTDIR,
    "Figure4A_Component1_OA_DIL_scores.png"
  ),
  p_scores,
  width = 5.5,
  height = 5.2,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    OUTDIR,
    "Figure4A_Component1_OA_DIL_scores.pdf"
  ),
  p_scores,
  width = 5.5,
  height = 5.2,
  device = cairo_pdf,
  bg = "white"
)

## ============================================================
## 13. ASSIGN LOADING-DERIVED DIRECTIONS
## ============================================================

# Work on a copy so integrated network is preserved
net_dir <- net0

net_dir$Gene_sign <- sign(
  gene_load[
    net_dir$Gene
  ]
)

net_dir$Met_sign <- sign(
  met_load[
    net_dir$Metabolite
  ]
)

net_dir$Direction <- "Mixed"

net_dir$Direction[
  net_dir$Gene_sign > 0 &
    net_dir$Met_sign > 0
] <- "OA_associated"

net_dir$Direction[
  net_dir$Gene_sign < 0 &
    net_dir$Met_sign < 0
] <- "DIL_associated"

cat("\nDirectional edge counts:\n")
print(
  table(
    net_dir$Direction
  )
)

write.csv(
  net_dir,
  file.path(
    OUTDIR,
    "Component1_directional_network_all_edges.csv"
  ),
  row.names = FALSE
)

## ============================================================
## 14. FIGURE 4B / 4C: DIL- AND OA-ASSOCIATED NETWORKS
##     Publication version with FIXED label sizes
## ============================================================

net_OA <- subset(
  net_dir,
  Direction == "OA_associated"
)

net_DIL <- subset(
  net_dir,
  Direction == "DIL_associated"
)

net_OA <- net_OA[
  order(
    abs(net_OA$Weight_Product_raw),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

net_DIL <- net_DIL[
  order(
    abs(net_DIL$Weight_Product_raw),
    decreasing = TRUE
  ),
  ,
  drop = FALSE
]

net_OA <- head(
  net_OA,
  TOP_NETWORK_EDGES
)

net_DIL <- head(
  net_DIL,
  TOP_NETWORK_EDGES
)

write.csv(
  net_DIL,
  file.path(
    OUTDIR,
    "Figure4B_DIL_associated_top20_edges.csv"
  ),
  row.names = FALSE
)

write.csv(
  net_OA,
  file.path(
    OUTDIR,
    "Figure4C_OA_associated_top20_edges.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# Custom fixed-style bipartite network plotting function
# This ensures identical internal label size in both panels.
# ------------------------------------------------------------

plot_fixed_bipartite <- function(
    net_df,
    out_png,
    panel_title,
    label_cex = 0.85,
    node_cex_gene = 1.55,
    node_cex_met = 1.55,
    edge_width_min = 1.0,
    edge_width_max = 6.0,
    width_px = 1800,
    height_px = 1800,
    res = 250) {

  if (nrow(net_df) == 0) {
    stop(
      paste(
        "No edges available for:",
        panel_title
      )
    )
  }

  # Build unique node table
  genes <- unique(
    as.character(net_df$Gene)
  )

  mets <- unique(
    as.character(net_df$Metabolite)
  )

  vertices <- data.frame(
    name = c(genes, mets),
    type = c(
      rep("Gene", length(genes)),
      rep("Metabolite", length(mets))
    ),
    stringsAsFactors = FALSE
  )

  # Build graph
  edges <- data.frame(
    from = as.character(net_df$Gene),
    to = as.character(net_df$Metabolite),
    stringsAsFactors = FALSE
  )

  g <- igraph::graph_from_data_frame(
    d = edges,
    directed = FALSE,
    vertices = vertices
  )

  # Node appearance
  V(g)$color <- ifelse(
    V(g)$type == "Gene",
    "#1F77B4",
    "#FF7F0E"
  )

  V(g)$size <- ifelse(
    V(g)$type == "Gene",
    18 * node_cex_gene,
    18 * node_cex_met
  )

  V(g)$frame.color <- "black"
  V(g)$frame.width <- 0.8

  # Fixed node-label styling
  V(g)$label.color <- "black"
  V(g)$label.cex <- label_cex
  V(g)$label.family <- "sans"

  # Edge widths based on absolute loading-product score
  ww <- abs(
    net_df$Weight_Product_raw
  )

  if (all(!is.finite(ww))) {
    ww <- rep(1, length(ww))
  }

  ww[
    !is.finite(ww)
  ] <- min(
    ww[is.finite(ww)],
    na.rm = TRUE
  )

  if (
    length(unique(ww)) <= 1
  ) {
    ew <- rep(
      (edge_width_min + edge_width_max) / 2,
      length(ww)
    )
  } else {
    ew <- edge_width_min +
      (
        (ww - min(ww)) /
          (max(ww) - min(ww))
      ) *
      (
        edge_width_max - edge_width_min
      )
  }

  E(g)$width <- ew
  E(g)$color <- "#2A9D2F"

  # Deterministic Fruchterman-Reingold layout
  set.seed(20260825)
  lay <- igraph::layout_with_fr(
    g,
    niter = 2000
  )

  png(
    file.path(
      OUTDIR,
      out_png
    ),
    width = width_px,
    height = height_px,
    res = res,
    bg = "white"
  )

  par(
    mar = c(5.5, 2.5, 5.5, 2.5),
    xpd = NA
  )

  plot(
    g,
    layout = lay,
    main = panel_title,
    vertex.label.cex = label_cex,
    vertex.label.color = "black",
    vertex.label.family = "sans",
    vertex.size = V(g)$size,
    vertex.color = V(g)$color,
    vertex.frame.color = "black",
    edge.width = E(g)$width,
    edge.color = E(g)$color,
    edge.curved = 0,
    asp = 1,
    rescale = TRUE
  )

  # Legend with identical size across panels
  legend(
    "bottomleft",
    legend = c(
      "Gene",
      "Metabolite"
    ),
    pch = 21,
    pt.bg = c(
      "#1F77B4",
      "#FF7F0E"
    ),
    pt.cex = 1.6,
    cex = 1.0,
    bty = "n",
    title = "Layers"
  )

  mtext(
    "Top 20 loading-derived associations",
    side = 1,
    line = 2.8,
    cex = 1.0
  )

  dev.off()
}

# DIL must be left in final compiled Figure 4
plot_fixed_bipartite(
  net_df = net_DIL,
  out_png = "Figure4B_DIL_associated_network_PUBLICATION.png",
  panel_title = "DIL-associated loading direction",
  label_cex = 0.85
)

# OA must be right in final compiled Figure 4
plot_fixed_bipartite(
  net_df = net_OA,
  out_png = "Figure4C_OA_associated_network_PUBLICATION.png",
  panel_title = "OA-associated loading direction",
  label_cex = 0.85
)

## ============================================================
## 15. EXPLORATORY DIRECTIONAL HEATMAPS
##     (saved separately; recommended for supplementary use)
## ============================================================

top_names <- function(
    load_vec,
    n = 25) {

  load_vec <- load_vec[
    is.finite(load_vec)
  ]

  if (length(load_vec) == 0) {
    return(character(0))
  }

  names(
    sort(
      abs(load_vec),
      decreasing = TRUE
    )
  )[
    seq_len(
      min(
        n,
        length(load_vec)
      )
    )
  ]
}

make_subset_model <- function(
    model,
    keep_genes,
    keep_mets) {

  m <- model

  m$loadings$X <- m$loadings$X[
    intersect(
      rownames(m$loadings$X),
      keep_genes
    ),
    ,
    drop = FALSE
  ]

  m$loadings$Y <- m$loadings$Y[
    intersect(
      rownames(m$loadings$Y),
      keep_mets
    ),
    ,
    drop = FALSE
  ]

  if (
    nrow(m$loadings$X) < 2 ||
      nrow(m$loadings$Y) < 2
  ) {
    stop(
      "Subset model has too few features for directional heatmap."
    )
  }

  m
}

genes_OA_assoc <- top_names(
  gene_load[
    gene_load > 0
  ],
  n = TOP_HEATMAP_FEATURES
)

mets_OA_assoc <- top_names(
  met_load[
    met_load > 0
  ],
  n = TOP_HEATMAP_FEATURES
)

genes_DIL_assoc <- top_names(
  gene_load[
    gene_load < 0
  ],
  n = TOP_HEATMAP_FEATURES
)

mets_DIL_assoc <- top_names(
  met_load[
    met_load < 0
  ],
  n = TOP_HEATMAP_FEATURES
)

scca_model_OA_assoc <- make_subset_model(
  scca_model,
  genes_OA_assoc,
  mets_OA_assoc
)

scca_model_DIL_assoc <- make_subset_model(
  scca_model,
  genes_DIL_assoc,
  mets_DIL_assoc
)

png(
  file.path(
    OUTDIR,
    "Supplementary_Heatmap_OA_associated_direction.png"
  ),
  width = 1800,
  height = 1800,
  res = 250
)

plot_correlation_heatmap(
  scca_model = scca_model_OA_assoc,
  X = X_used,
  Y = Y_used,
  top_n = TOP_HEATMAP_FEATURES
)

dev.off()

png(
  file.path(
    OUTDIR,
    "Supplementary_Heatmap_DIL_associated_direction.png"
  ),
  width = 1800,
  height = 1800,
  res = 250
)

plot_correlation_heatmap(
  scca_model = scca_model_DIL_assoc,
  X = X_used,
  Y = Y_used,
  top_n = TOP_HEATMAP_FEATURES
)

dev.off()

## ============================================================
## 16. OPTIONAL GROUP-SUBSET ANALYSES
##     Exploratory only; n = 5 per group.
##     Do not use these as primary treatment-specific evidence.
## ============================================================

oa_samples <- rownames(X_used)[
  grepl(
    "^OA",
    rownames(X_used)
  )
]

dil_samples <- rownames(X_used)[
  grepl(
    "^DIL",
    rownames(X_used)
  )
]

X_OA_only <- X_used[
  oa_samples,
  ,
  drop = FALSE
]

Y_OA_only <- Y_used[
  oa_samples,
  ,
  drop = FALSE
]

X_DIL_only <- X_used[
  dil_samples,
  ,
  drop = FALSE
]

Y_DIL_only <- Y_used[
  dil_samples,
  ,
  drop = FALSE
]

# These fits are intentionally wrapped in try()
# because the sample size is very small.
scca_OA_only <- try(
  omic_scca(
    X = X_OA_only,
    Y = Y_OA_only,
    n_components = N_COMPONENTS,
    penalty_X = PENALTY_X,
    penalty_Y = PENALTY_Y
  ),
  silent = TRUE
)

scca_DIL_only <- try(
  omic_scca(
    X = X_DIL_only,
    Y = Y_DIL_only,
    n_components = N_COMPONENTS,
    penalty_X = PENALTY_X,
    penalty_Y = PENALTY_Y
  ),
  silent = TRUE
)

if (!inherits(scca_OA_only, "try-error")) {

  save_vip_circle(
    scca_OA_only,
    "Exploratory_OA_only",
    top_features = TOP_VIP_FEATURES,
    pathway_db = "KEGG"
  )

  save_heatmap(
    scca_OA_only,
    X_OA_only,
    Y_OA_only,
    "Exploratory_OA_only_heatmap.png",
    top_n = TOP_HEATMAP_FEATURES
  )
}

if (!inherits(scca_DIL_only, "try-error")) {

  save_vip_circle(
    scca_DIL_only,
    "Exploratory_DIL_only",
    top_features = TOP_VIP_FEATURES,
    pathway_db = "KEGG"
  )

  save_heatmap(
    scca_DIL_only,
    X_DIL_only,
    Y_DIL_only,
    "Exploratory_DIL_only_heatmap.png",
    top_n = TOP_HEATMAP_FEATURES
  )
}

## ============================================================
## 17. CREATE FINAL COMBINED PUBLICATION FIGURE 4
##
## FINAL LAYOUT:
##
##                 A. Component 1 canonical scores
##
##                    DIL             OA
##                     ↓               ↓
##
##        B. DIL-associated   C. OA-associated
##              network            network
##
## ============================================================

# ------------------------------------------------------------
# Force Panel A order: DIL LEFT, OA RIGHT
# ------------------------------------------------------------

scores_df$Group <- factor(
  scores_df$Group,
  levels = c(
    "DIL",
    "OA"
  )
)

p_scores_final <- ggplot(
  scores_df,
  aes(
    x = Group,
    y = Score_X_oriented,
    fill = Group
  )
) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_jitter(
    width = 0.08,
    size = 3.2,
    alpha = 0.90
  ) +
  theme_classic(
    base_size = 15
  ) +
  guides(
    fill = "none"
  ) +
  labs(
    title = "Treatment association with cross-omics Component 1",
    x = NULL,
    y = "Oriented Component 1 canonical score"
  ) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.y = element_text(
      size = 14
    ),
    axis.text.x = element_text(
      size = 14,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 12
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )

# Save standalone corrected Panel A too
ggsave(
  file.path(
    OUTDIR,
    "Figure4A_Component1_OA_DIL_scores_PUBLICATION.png"
  ),
  p_scores_final,
  width = 8.5,
  height = 4.6,
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# Read fixed-style network panels
# ------------------------------------------------------------

dil_png <- readPNG(
  file.path(
    OUTDIR,
    "Figure4B_DIL_associated_network_PUBLICATION.png"
  )
)

oa_png <- readPNG(
  file.path(
    OUTDIR,
    "Figure4C_OA_associated_network_PUBLICATION.png"
  )
)

p_dil_network <- wrap_elements(
  full = rasterGrob(
    dil_png,
    interpolate = TRUE
  )
)

p_oa_network <- wrap_elements(
  full = rasterGrob(
    oa_png,
    interpolate = TRUE
  )
)

# ------------------------------------------------------------
# Combine:
# DIL boxplot is directly above DIL network.
# OA boxplot is directly above OA network.
#
# Because Panel A spans both columns and its x-axis ordering is
# DIL then OA, the vertical correspondence is unambiguous.
# ------------------------------------------------------------

Figure4_FINAL_PUBLICATION <-
  p_scores_final /
  (
    p_dil_network |
      p_oa_network
  ) +
  plot_layout(
    heights = c(
      0.82,
      1.18
    )
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      size = 20,
      face = "bold"
    )
  )

print(
  Figure4_FINAL_PUBLICATION
)

# ------------------------------------------------------------
# Save final publication PNG
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    OUTDIR,
    "Figure_4_FINAL_PUBLICATION_ALIGNED.png"
  ),
  plot = Figure4_FINAL_PUBLICATION,
  width = 12,
  height = 11,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# Save final publication PDF
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    OUTDIR,
    "Figure_4_FINAL_PUBLICATION_ALIGNED.pdf"
  ),
  plot = Figure4_FINAL_PUBLICATION,
  width = 12,
  height = 11,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

cat(
  "\n============================================================\n"
)

cat(
  "FINAL ALIGNED PUBLICATION FIGURE 4 SAVED\n"
)

cat(
  "============================================================\n"
)

cat(
  "Panel A: DIL is LEFT; OA is RIGHT\n"
)

cat(
  "Panel B: DIL-associated network is LEFT\n"
)

cat(
  "Panel C: OA-associated network is RIGHT\n"
)

cat(
  "Network label size is fixed identically in B and C.\n"
)

cat(
  "Network footer uses 'loading-derived associations', not 'interactions'.\n\n"
)

cat(
  file.path(
    OUTDIR,
    "Figure_4_FINAL_PUBLICATION_ALIGNED.png"
  ),
  "\n"
)

cat(
  file.path(
    OUTDIR,
    "Figure_4_FINAL_PUBLICATION_ALIGNED.pdf"
  ),
  "\n"
)

## ============================================================
## 18. SAVE ORIENTATION INFORMATION / SESSION INFO
## ============================================================

orientation_info <- data.frame(
  Component = 1,
  Orientation_flipped = ORIENTATION_FLIPPED,
  Positive_direction_label = "OA-associated",
  Negative_direction_label = "DIL-associated",
  Interpretation = paste(
    "Labels indicate association with the post-hoc oriented latent component;",
    "they do not imply differential enrichment or causal treatment effects."
  )
)

write.csv(
  orientation_info,
  file.path(
    OUTDIR,
    "Component1_orientation_information.csv"
  ),
  row.names = FALSE
)

sink(
  file.path(
    OUTDIR,
    "sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

cat(
  "\n============================================================\n"
)

cat(
  "MC903 OmicNetR FINAL ANALYSIS COMPLETE\n"
)

cat(
  "============================================================\n"
)

cat(
  "Output directory: ",
  normalizePath(OUTDIR),
  "\n\n",
  sep = ""
)

cat(
  "Primary publication Figure 4:\n"
)

cat(
  "  Figure_4_FINAL_Treatment_Associated_Directionality.png\n"
)

cat(
  "  Figure_4_FINAL_Treatment_Associated_Directionality.pdf\n\n"
)

cat(
  "IMPORTANT INTERPRETATION:\n"
)

cat(
  "  OA/DIL labels are based on post-hoc orientation of Component 1 using\n"
)

cat(
  "  actual treatment-group canonical sample scores. They indicate direction\n"
)

cat(
  "  along the integrated latent component and do NOT establish direct\n"
)

cat(
  "  differential enrichment, molecular interaction, or causal treatment effects.\n"
)

