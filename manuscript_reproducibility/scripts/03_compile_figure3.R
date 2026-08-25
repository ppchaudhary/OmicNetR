
## ============================================================
## FINAL PUBLICATION FIGURE 3
## Integrated MC903 multi-omics feature landscape
##
## Panel A = Circular feature-importance plot
## Panel B = Integrated gene-metabolite correlation heatmap
## Panel C = Cross-omics association network
##
## This script assumes you have already run:
## Real_data_OmicNetR_script_FINAL_ALIGNED_PUBLICATION.R
##
## and that the folder below contains:
##   Figure3A_Integrated_VIP_circle.png
##   Figure3B_Heatmap_integrated.png
##   Integrated_top20_network_edges.csv
## ============================================================

## ----------------------------
## 0. INSTALL / LOAD PACKAGES
## ----------------------------
pkgs <- c(
  "ggplot2",
  "patchwork",
  "png",
  "grid",
  "igraph"
)

missing_pkgs <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  install.packages(
    missing_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(grid)
  library(igraph)
})

## ----------------------------
## 1. SETTINGS
## ----------------------------
OUTDIR <- "OmicNetR_MC903_FINAL"

VIP_FILE <- file.path(
  OUTDIR,
  "Figure3A_Integrated_VIP_circle.png"
)

HEATMAP_FILE <- file.path(
  OUTDIR,
  "Figure3B_Heatmap_integrated.png"
)

EDGE_FILE <- file.path(
  OUTDIR,
  "Integrated_top20_network_edges.csv"
)

FINAL_PNG <- file.path(
  OUTDIR,
  "Figure_3_FINAL_PUBLICATION.png"
)

FINAL_PDF <- file.path(
  OUTDIR,
  "Figure_3_FINAL_PUBLICATION.pdf"
)

## ----------------------------
## 2. CHECK INPUT FILES
## ----------------------------
required_files <- c(
  VIP_FILE,
  HEATMAP_FILE,
  EDGE_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "These required files are missing:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

## ============================================================
## 3. REBUILD PANEL C WITH PUBLICATION-QUALITY FIXED LABELS
## ============================================================

net_top <- read.csv(
  EDGE_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Determine which weight column is available
if ("Weight_Product_raw" %in% names(net_top)) {
  raw_weight <- net_top$Weight_Product_raw
} else if ("Weight_Product" %in% names(net_top)) {
  raw_weight <- net_top$Weight_Product
} else {
  stop(
    "No Weight_Product or Weight_Product_raw column found in edge file."
  )
}

net_top$PlotWeight <- raw_weight

# ------------------------------------------------------------
# Build node table
# ------------------------------------------------------------
genes <- unique(
  as.character(
    net_top$Gene
  )
)

mets <- unique(
  as.character(
    net_top$Metabolite
  )
)

vertices <- data.frame(
  name = c(
    genes,
    mets
  ),
  type = c(
    rep(
      "Gene",
      length(genes)
    ),
    rep(
      "Metabolite",
      length(mets)
    )
  ),
  stringsAsFactors = FALSE
)

edges <- data.frame(
  from = as.character(
    net_top$Gene
  ),
  to = as.character(
    net_top$Metabolite
  ),
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(
  d = edges,
  directed = FALSE,
  vertices = vertices
)

## ------------------------------------------------------------
## Fixed node styling
## ------------------------------------------------------------
V(g)$color <- ifelse(
  V(g)$type == "Gene",
  "#1F77B4",
  "#FF7F0E"
)

V(g)$size <- ifelse(
  V(g)$type == "Gene",
  24,
  18
)

V(g)$frame.color <- "black"
V(g)$frame.width <- 0.7

# IDENTICAL label size for every node
V(g)$label.color <- "black"
V(g)$label.cex <- 0.90
V(g)$label.family <- "sans"

## ------------------------------------------------------------
## Edge styling
##
## Width = absolute loading-product magnitude
## Color = sign of loading-product score
##
## Positive = green
## Negative = red
## ------------------------------------------------------------
w <- net_top$PlotWeight

finite_w <- abs(
  w[
    is.finite(w)
  ]
)

if (length(finite_w) == 0) {
  edge_widths <- rep(
    2,
    length(w)
  )
} else {

  abs_w <- abs(w)

  abs_w[
    !is.finite(abs_w)
  ] <- min(
    finite_w
  )

  if (length(unique(abs_w)) <= 1) {

    edge_widths <- rep(
      3,
      length(abs_w)
    )

  } else {

    edge_widths <- 1.0 +
      (
        abs_w - min(abs_w)
      ) /
      (
        max(abs_w) - min(abs_w)
      ) *
      5.0
  }
}

E(g)$width <- edge_widths

E(g)$color <- ifelse(
  w >= 0,
  "#2A9D2F",
  "#D62828"
)

## ------------------------------------------------------------
## Deterministic layout
## ------------------------------------------------------------
set.seed(20260825)

lay <- layout_with_fr(
  g,
  niter = 2500
)

NETWORK_FILE <- file.path(
  OUTDIR,
  "Figure3C_CrossOmics_Association_Network_PUBLICATION.png"
)

png(
  NETWORK_FILE,
  width = 2200,
  height = 2200,
  res = 300,
  bg = "white"
)

par(
  mar = c(
    6.0,
    3.0,
    5.0,
    3.0
  ),
  xpd = NA
)

plot(
  g,
  layout = lay,
  main = "Cross-Omics Association Network",
  vertex.size = V(g)$size,
  vertex.color = V(g)$color,
  vertex.frame.color = V(g)$frame.color,
  vertex.label.cex = 0.90,
  vertex.label.color = "black",
  vertex.label.family = "sans",
  edge.width = E(g)$width,
  edge.color = E(g)$color,
  edge.curved = 0,
  asp = 1,
  rescale = TRUE
)

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

legend(
  "bottomright",
  legend = c(
    "Positive association",
    "Negative association"
  ),
  lwd = 4,
  col = c(
    "#2A9D2F",
    "#D62828"
  ),
  cex = 0.90,
  bty = "n",
  title = "Loading-product sign"
)

mtext(
  "Top 20 loading-derived associations",
  side = 1,
  line = 3.2,
  cex = 1.0
)

dev.off()

## ============================================================
## 4. LOAD PANEL A / PANEL B / PANEL C AS RASTER PANELS
## ============================================================

vip_png <- readPNG(
  VIP_FILE
)

heatmap_png <- readPNG(
  HEATMAP_FILE
)

network_png <- readPNG(
  NETWORK_FILE
)

p_vip <- wrap_elements(
  full = rasterGrob(
    vip_png,
    interpolate = TRUE
  )
)

p_heatmap <- wrap_elements(
  full = rasterGrob(
    heatmap_png,
    interpolate = TRUE
  )
)

p_network <- wrap_elements(
  full = rasterGrob(
    network_png,
    interpolate = TRUE
  )
)

## ============================================================
## 5. CREATE CONSISTENT PANEL TITLES
## ============================================================

title_theme <- theme(
  plot.title = element_text(
    size = 15,
    face = "bold",
    hjust = 0.5,
    margin = margin(
      b = 6
    )
  )
)

p_vip_final <-
  p_vip +
  plot_annotation(
    title = "Canonical Loading-Based Feature Importance"
  ) &
  title_theme

p_heatmap_final <-
  p_heatmap +
  plot_annotation(
    title = "Cross-Omics Feature Correlation Landscape"
  ) &
  title_theme

p_network_final <-
  p_network +
  plot_annotation(
    title = "Loading-Derived Cross-Omics Association Network"
  ) &
  title_theme

## ============================================================
## 6. FINAL PUBLICATION LAYOUT
##
##         A                     B
##   Feature importance       Heatmap
##
##               C
##        Association network
##
## This layout gives the network enough width and avoids
## shrinking labels.
## ============================================================

top_row <-
  p_vip_final |
  p_heatmap_final

bottom_row <-
  p_network_final

Figure3_FINAL <-
  top_row /
  bottom_row +
  plot_layout(
    heights = c(
      1,
      1.05
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

## ----------------------------
## 7. DISPLAY
## ----------------------------
print(
  Figure3_FINAL
)

## ============================================================
## 8. SAVE PUBLICATION FILES
## ============================================================

ggsave(
  filename = FINAL_PNG,
  plot = Figure3_FINAL,
  width = 12,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = FINAL_PDF,
  plot = Figure3_FINAL,
  width = 12,
  height = 12,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

## ============================================================
## 9. SAVE FIGURE 3 EDGE-TABLE INTERPRETATION
## ============================================================

edge_interpretation <- data.frame(
  Visual_property = c(
    "Blue node",
    "Orange node",
    "Edge width",
    "Green edge",
    "Red edge"
  ),
  Interpretation = c(
    "Gene",
    "Metabolite",
    "Magnitude of loading-derived association score",
    "Positive loading-product association",
    "Negative loading-product association"
  )
)

write.csv(
  edge_interpretation,
  file.path(
    OUTDIR,
    "Figure3_network_visual_encoding.csv"
  ),
  row.names = FALSE
)

cat(
  "\n============================================================\n"
)

cat(
  "FINAL PUBLICATION FIGURE 3 CREATED\n"
)

cat(
  "============================================================\n"
)

cat(
  "Panel A = canonical loading-based feature importance\n"
)

cat(
  "Panel B = cross-omics feature correlation landscape\n"
)

cat(
  "Panel C = top 20 loading-derived cross-omics associations\n\n"
)

cat(
  "Saved as:\n",
  FINAL_PNG,
  "\n",
  FINAL_PDF,
  "\n"
)
