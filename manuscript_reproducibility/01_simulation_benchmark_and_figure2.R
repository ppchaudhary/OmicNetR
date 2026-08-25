# ============================================================
# OmicNetR vs mixOmics: simulation benchmark
# Purpose: manuscript-quality comparison on identical simulated data
# Metrics: feature recovery, held-out canonical correlation, runtime,
#          recovery consistency, optional bootstrap selection stability
# ============================================================

# ----------------------------
# 0. USER SETTINGS
# ----------------------------
N_REPS <- 100          # Start with 10 for a quick test, then use 100 for paper
SAMPLE_SIZES <- c(30, 60, 120)
SIGNAL_STRENGTHS <- c(0.5, 1.0, 1.5, 2.0)
N_GENES <- 800
N_METS <- 150
N_LINKED <- 20
N_COMPONENTS <- 2
TRAIN_FRAC <- 0.70
OMIC_PENALTY_X <- 0.90
OMIC_PENALTY_Y <- 0.90
SELECTION_K_GENES <- N_LINKED
SELECTION_K_METS <- N_LINKED
OUTPUT_DIR <- "OmicNetR_benchmark_results"
MASTER_SEED <- 20260824
RUN_BOOTSTRAP_STABILITY <- TRUE
N_BOOT <- 100           # Final manuscript analysis
BOOTSTRAP_N <- 60
BOOTSTRAP_SIGNALS <- c(0.5, 1.0, 1.5)

# ----------------------------
# 1. INSTALL / LOAD PACKAGES
# ----------------------------
cran_pkgs <- c("OmicNetR", "ggplot2", "dplyr", "tidyr", "readr", "purrr")
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("mixOmics", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install("mixOmics", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(OmicNetR)
  library(mixOmics)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
})

set.seed(MASTER_SEED)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2. SIMULATE MULTI-OMICS DATA
# ----------------------------
# Two latent biological factors are used. Half of the true genes/metabolites
# load on latent factor 1 and the other half on latent factor 2.
# All remaining variables are null/noise variables.

simulate_multiomics <- function(
    n_samples,
    n_genes = 800,
    n_metabolites = 150,
    n_linked = 20,
    signal_strength = 1.5,
    noise_sd = 1,
    n_latent = 2,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  stopifnot(n_linked %% n_latent == 0)

  Z <- matrix(rnorm(n_samples * n_latent), nrow = n_samples, ncol = n_latent)

  X <- matrix(rnorm(n_samples * n_genes, sd = noise_sd),
              nrow = n_samples, ncol = n_genes)
  Y <- matrix(rnorm(n_samples * n_metabolites, sd = noise_sd),
              nrow = n_samples, ncol = n_metabolites)

  colnames(X) <- paste0("Gene_", seq_len(n_genes))
  colnames(Y) <- paste0("Met_", seq_len(n_metabolites))
  rownames(X) <- paste0("Sample_", seq_len(n_samples))
  rownames(Y) <- paste0("Sample_", seq_len(n_samples))

  per_factor <- n_linked / n_latent
  true_genes <- paste0("Gene_", seq_len(n_linked))
  true_mets <- paste0("Met_", seq_len(n_linked))

  factor_membership_gene <- rep(seq_len(n_latent), each = per_factor)
  factor_membership_met <- rep(seq_len(n_latent), each = per_factor)

  # Add shared latent structure to true features.
  # Feature-specific coefficients vary mildly to avoid identical features.
  for (j in seq_len(n_linked)) {
    f <- factor_membership_gene[j]
    coef_x <- signal_strength * runif(1, 0.8, 1.2) * sample(c(-1, 1), 1)
    X[, j] <- coef_x * Z[, f] + rnorm(n_samples, sd = noise_sd)
  }

  for (j in seq_len(n_linked)) {
    f <- factor_membership_met[j]
    coef_y <- signal_strength * runif(1, 0.8, 1.2) * sample(c(-1, 1), 1)
    Y[, j] <- coef_y * Z[, f] + rnorm(n_samples, sd = noise_sd)
  }

  list(
    X = X,
    Y = Y,
    true_genes = true_genes,
    true_mets = true_mets,
    gene_factor = setNames(factor_membership_gene, true_genes),
    met_factor = setNames(factor_membership_met, true_mets)
  )
}

# ----------------------------
# 3. TRAIN/TEST STANDARDIZATION
# ----------------------------
scale_train_test <- function(train, test) {
  mu <- colMeans(train, na.rm = TRUE)
  sdv <- apply(train, 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  train_s <- sweep(train, 2, mu, "-")
  train_s <- sweep(train_s, 2, sdv, "/")
  test_s <- sweep(test, 2, mu, "-")
  test_s <- sweep(test_s, 2, sdv, "/")

  list(train = train_s, test = test_s, center = mu, scale = sdv)
}

# ----------------------------
# 4. HELPERS FOR LOADINGS / FEATURE SELECTION
# ----------------------------
loading_importance <- function(L) {
  L <- as.matrix(L)
  if (is.null(rownames(L))) stop("Loading matrix has no feature names in rownames().")
  apply(abs(L), 1, max, na.rm = TRUE)
}

top_k_features <- function(L, k) {
  imp <- loading_importance(L)
  k <- min(k, length(imp))
  names(sort(imp, decreasing = TRUE))[seq_len(k)]
}

nonzero_features <- function(L, tol = 1e-10) {
  imp <- loading_importance(L)
  names(imp)[imp > tol]
}

feature_metrics <- function(selected, truth) {
  selected <- unique(selected)
  truth <- unique(truth)
  TP <- length(intersect(selected, truth))
  FP <- length(setdiff(selected, truth))
  FN <- length(setdiff(truth, selected))
  precision <- if ((TP + FP) == 0) 0 else TP / (TP + FP)
  recall <- if ((TP + FN) == 0) 0 else TP / (TP + FN)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)

  tibble(
    TP = TP, FP = FP, FN = FN,
    Precision = precision,
    Recall = recall,
    F1 = f1,
    N_selected = length(selected)
  )
}

safe_cor <- function(a, b) {
  if (length(a) < 3 || length(b) < 3 || sd(a) == 0 || sd(b) == 0) return(NA_real_)
  suppressWarnings(cor(a, b, use = "pairwise.complete.obs"))
}

heldout_canonical_cor <- function(X_test, Y_test, LX, LY, ncomp = 2) {
  LX <- as.matrix(LX)
  LY <- as.matrix(LY)
  ncomp <- min(ncomp, ncol(LX), ncol(LY))

  # Match test columns to loading rows.
  X_test <- X_test[, rownames(LX), drop = FALSE]
  Y_test <- Y_test[, rownames(LY), drop = FALSE]

  cors <- vapply(seq_len(ncomp), function(k) {
    tx <- as.numeric(X_test %*% LX[, k])
    ty <- as.numeric(Y_test %*% LY[, k])
    abs(safe_cor(tx, ty))
  }, numeric(1))

  mean(cors, na.rm = TRUE)
}

# ----------------------------
# 5. FIT OmicNetR
# ----------------------------
fit_omicnetr <- function(X_train, Y_train,
                         ncomp = 2,
                         penalty_x = 0.90,
                         penalty_y = 0.90) {

  aligned <- align_omics(X_train, Y_train)

  t0 <- proc.time()[[3]]
  fit <- omic_scca(
    X = aligned$X,
    Y = aligned$Y,
    n_components = ncomp,
    penalty_X = penalty_x,
    penalty_Y = penalty_y
  )
  runtime <- proc.time()[[3]] - t0

  list(
    fit = fit,
    LX = as.matrix(fit$loadings$X),
    LY = as.matrix(fit$loadings$Y),
    runtime = runtime
  )
}

# ----------------------------
# 6. FIT mixOmics sPLS IN CANONICAL MODE
# ----------------------------
fit_mixomics <- function(X_train, Y_train,
                         ncomp = 2,
                         keepX = 20,
                         keepY = 20) {

  # Same selection budget per component for a controlled benchmark.
  keepX_vec <- rep(min(keepX, ncol(X_train)), ncomp)
  keepY_vec <- rep(min(keepY, ncol(Y_train)), ncomp)

  t0 <- proc.time()[[3]]
  fit <- mixOmics::spls(
    X = X_train,
    Y = Y_train,
    ncomp = ncomp,
    mode = "canonical",
    keepX = keepX_vec,
    keepY = keepY_vec,
    scale = FALSE
  )
  runtime <- proc.time()[[3]] - t0

  list(
    fit = fit,
    LX = as.matrix(fit$loadings$X),
    LY = as.matrix(fit$loadings$Y),
    runtime = runtime
  )
}

# ----------------------------
# 7. ONE REPLICATE
# ----------------------------
run_one_replicate <- function(n_samples, signal_strength, rep_id) {

  seed_i <- MASTER_SEED + rep_id + as.integer(n_samples * 1000 + signal_strength * 100)

  sim <- simulate_multiomics(
    n_samples = n_samples,
    n_genes = N_GENES,
    n_metabolites = N_METS,
    n_linked = N_LINKED,
    signal_strength = signal_strength,
    noise_sd = 1,
    n_latent = N_COMPONENTS,
    seed = seed_i
  )

  # Train/test split shared by both methods.
  set.seed(seed_i + 999)
  n_train <- max(10, floor(TRAIN_FRAC * n_samples))
  train_idx <- sort(sample(seq_len(n_samples), size = n_train, replace = FALSE))
  test_idx <- setdiff(seq_len(n_samples), train_idx)

  X_train_raw <- sim$X[train_idx, , drop = FALSE]
  X_test_raw  <- sim$X[test_idx, , drop = FALSE]
  Y_train_raw <- sim$Y[train_idx, , drop = FALSE]
  Y_test_raw  <- sim$Y[test_idx, , drop = FALSE]

  sx <- scale_train_test(X_train_raw, X_test_raw)
  sy <- scale_train_test(Y_train_raw, Y_test_raw)

  X_train <- sx$train; X_test <- sx$test
  Y_train <- sy$train; Y_test <- sy$test

  # ---- OmicNetR ----
  om <- tryCatch(
    fit_omicnetr(X_train, Y_train,
                 ncomp = N_COMPONENTS,
                 penalty_x = OMIC_PENALTY_X,
                 penalty_y = OMIC_PENALTY_Y),
    error = function(e) e
  )

  # ---- mixOmics ----
  mx <- tryCatch(
    fit_mixomics(X_train, Y_train,
                 ncomp = N_COMPONENTS,
                 keepX = SELECTION_K_GENES,
                 keepY = SELECTION_K_METS),
    error = function(e) e
  )

  rows <- list()
  selections <- list()

  if (!inherits(om, "error")) {
    # Primary fair comparison: top-K features by max absolute loading.
    sel_g <- top_k_features(om$LX, SELECTION_K_GENES)
    sel_m <- top_k_features(om$LY, SELECTION_K_METS)
    native_g <- nonzero_features(om$LX)
    native_m <- nonzero_features(om$LY)

    gmet <- feature_metrics(sel_g, sim$true_genes)
    mmet <- feature_metrics(sel_m, sim$true_mets)
    test_cor <- heldout_canonical_cor(X_test, Y_test, om$LX, om$LY, N_COMPONENTS)

    rows[[length(rows) + 1]] <- gmet %>% mutate(
      Method = "OmicNetR", Omic = "Genes",
      n_samples = n_samples, signal_strength = signal_strength,
      Replicate = rep_id, Runtime_sec = om$runtime,
      Test_Canonical_Correlation = test_cor,
      Native_N_Selected = length(native_g), Status = "OK"
    )

    rows[[length(rows) + 1]] <- mmet %>% mutate(
      Method = "OmicNetR", Omic = "Metabolites",
      n_samples = n_samples, signal_strength = signal_strength,
      Replicate = rep_id, Runtime_sec = om$runtime,
      Test_Canonical_Correlation = test_cor,
      Native_N_Selected = length(native_m), Status = "OK"
    )

    selections[["OmicNetR_Genes"]] <- sel_g
    selections[["OmicNetR_Metabolites"]] <- sel_m
  } else {
    rows[[length(rows) + 1]] <- tibble(
      TP = NA, FP = NA, FN = NA, Precision = NA, Recall = NA, F1 = NA,
      N_selected = NA, Method = "OmicNetR", Omic = c("Genes", "Metabolites"),
      n_samples = n_samples, signal_strength = signal_strength,
      Replicate = rep_id, Runtime_sec = NA,
      Test_Canonical_Correlation = NA, Native_N_Selected = NA,
      Status = paste("ERROR:", conditionMessage(om))
    )
  }

  if (!inherits(mx, "error")) {
    sel_g <- top_k_features(mx$LX, SELECTION_K_GENES)
    sel_m <- top_k_features(mx$LY, SELECTION_K_METS)
    native_g <- nonzero_features(mx$LX)
    native_m <- nonzero_features(mx$LY)

    gmet <- feature_metrics(sel_g, sim$true_genes)
    mmet <- feature_metrics(sel_m, sim$true_mets)
    test_cor <- heldout_canonical_cor(X_test, Y_test, mx$LX, mx$LY, N_COMPONENTS)

    rows[[length(rows) + 1]] <- gmet %>% mutate(
      Method = "mixOmics", Omic = "Genes",
      n_samples = n_samples, signal_strength = signal_strength,
      Replicate = rep_id, Runtime_sec = mx$runtime,
      Test_Canonical_Correlation = test_cor,
      Native_N_Selected = length(native_g), Status = "OK"
    )

    rows[[length(rows) + 1]] <- mmet %>% mutate(
      Method = "mixOmics", Omic = "Metabolites",
      n_samples = n_samples, signal_strength = signal_strength,
      Replicate = rep_id, Runtime_sec = mx$runtime,
      Test_Canonical_Correlation = test_cor,
      Native_N_Selected = length(native_m), Status = "OK"
    )

    selections[["mixOmics_Genes"]] <- sel_g
    selections[["mixOmics_Metabolites"]] <- sel_m
  } else {
    rows[[length(rows) + 1]] <- tibble(
      TP = NA, FP = NA, FN = NA, Precision = NA, Recall = NA, F1 = NA,
      N_selected = NA, Method = "mixOmics", Omic = c("Genes", "Metabolites"),
      n_samples = n_samples, signal_strength = signal_strength,
      Replicate = rep_id, Runtime_sec = NA,
      Test_Canonical_Correlation = NA, Native_N_Selected = NA,
      Status = paste("ERROR:", conditionMessage(mx))
    )
  }

  list(metrics = bind_rows(rows), selections = selections)
}

# ----------------------------
# 8. RUN FULL SIMULATION GRID
# ----------------------------
scenario_grid <- expand.grid(
  n_samples = SAMPLE_SIZES,
  signal_strength = SIGNAL_STRENGTHS,
  Replicate = seq_len(N_REPS),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

cat("Running", nrow(scenario_grid), "simulated datasets...\n")
cat("Each dataset is evaluated by OmicNetR and mixOmics.\n\n")

all_metrics <- vector("list", nrow(scenario_grid))
all_selections <- vector("list", nrow(scenario_grid))

for (i in seq_len(nrow(scenario_grid))) {
  sc <- scenario_grid[i, ]

  if (i %% 10 == 0 || i == 1) {
    cat(sprintf("[%d/%d] n=%d, signal=%.1f, rep=%d\n",
                i, nrow(scenario_grid), sc$n_samples,
                sc$signal_strength, sc$Replicate))
  }

  ans <- run_one_replicate(
    n_samples = sc$n_samples,
    signal_strength = sc$signal_strength,
    rep_id = sc$Replicate
  )

  all_metrics[[i]] <- ans$metrics
  all_selections[[i]] <- ans$selections
}

raw_results <- bind_rows(all_metrics)
write_csv(raw_results, file.path(OUTPUT_DIR, "benchmark_raw_results.csv"))

# ----------------------------
# 9. SUMMARIZE RESULTS
# ----------------------------
summary_results <- raw_results %>%
  filter(Status == "OK") %>%
  group_by(Method, Omic, n_samples, signal_strength) %>%
  summarise(
    n_success = n(),
    Precision_mean = mean(Precision, na.rm = TRUE),
    Precision_sd = sd(Precision, na.rm = TRUE),
    Recall_mean = mean(Recall, na.rm = TRUE),
    Recall_sd = sd(Recall, na.rm = TRUE),
    F1_mean = mean(F1, na.rm = TRUE),
    F1_sd = sd(F1, na.rm = TRUE),
    TestCor_mean = mean(Test_Canonical_Correlation, na.rm = TRUE),
    TestCor_sd = sd(Test_Canonical_Correlation, na.rm = TRUE),
    Runtime_median_sec = median(Runtime_sec, na.rm = TRUE),
    Runtime_IQR_sec = IQR(Runtime_sec, na.rm = TRUE),
    Native_Selected_mean = mean(Native_N_Selected, na.rm = TRUE),
    .groups = "drop"
  )


write_csv(summary_results, file.path(OUTPUT_DIR, "benchmark_summary.csv"))

# ----------------------------
# 9B. PAIRED METHOD COMPARISONS
# ----------------------------
# Both methods are evaluated on the same simulated replicate.
# Therefore, differences in F1 and held-out canonical correlation are paired.

# ---- Paired F1 differences ----
f1_paired <- raw_results %>%
  filter(Status == "OK") %>%
  select(
    Method, Omic, n_samples, signal_strength,
    Replicate, F1
  ) %>%
  pivot_wider(
    names_from = Method,
    values_from = F1
  ) %>%
  filter(!is.na(OmicNetR), !is.na(mixOmics)) %>%
  mutate(
    Difference = OmicNetR - mixOmics
  )

f1_paired_summary <- f1_paired %>%
  group_by(Omic, n_samples, signal_strength) %>%
  summarise(
    N_pairs = n(),
    OmicNetR_mean = mean(OmicNetR, na.rm = TRUE),
    mixOmics_mean = mean(mixOmics, na.rm = TRUE),
    Mean_difference = mean(Difference, na.rm = TRUE),
    Median_difference = median(Difference, na.rm = TRUE),
    SD_difference = sd(Difference, na.rm = TRUE),
    SE_difference = SD_difference / sqrt(N_pairs),
    CI_low = Mean_difference -
      qt(0.975, df = pmax(N_pairs - 1, 1)) * SE_difference,
    CI_high = Mean_difference +
      qt(0.975, df = pmax(N_pairs - 1, 1)) * SE_difference,
    Paired_t_p = ifelse(
      N_pairs >= 2,
      t.test(Difference, mu = 0)$p.value,
      NA_real_
    ),
    Wilcoxon_p = ifelse(
      N_pairs >= 2 && any(abs(Difference) > 0, na.rm = TRUE),
      suppressWarnings(wilcox.test(
        Difference,
        mu = 0,
        exact = FALSE
      )$p.value),
      NA_real_
    ),
    .groups = "drop"
  )

write_csv(
  f1_paired,
  file.path(OUTPUT_DIR, "F1_paired_raw.csv")
)

write_csv(
  f1_paired_summary,
  file.path(OUTPUT_DIR, "F1_paired_comparison.csv")
)

# ---- Paired held-out canonical-correlation differences ----
# Test_Canonical_Correlation is duplicated across the two Omic rows,
# so keep one row per replicate/method before comparison.
cor_paired <- raw_results %>%
  filter(Status == "OK") %>%
  distinct(
    Method, n_samples, signal_strength,
    Replicate, Test_Canonical_Correlation
  ) %>%
  pivot_wider(
    names_from = Method,
    values_from = Test_Canonical_Correlation
  ) %>%
  filter(!is.na(OmicNetR), !is.na(mixOmics)) %>%
  mutate(
    Difference = OmicNetR - mixOmics
  )

cor_paired_summary <- cor_paired %>%
  group_by(n_samples, signal_strength) %>%
  summarise(
    N_pairs = n(),
    OmicNetR_mean = mean(OmicNetR, na.rm = TRUE),
    mixOmics_mean = mean(mixOmics, na.rm = TRUE),
    Mean_difference = mean(Difference, na.rm = TRUE),
    Median_difference = median(Difference, na.rm = TRUE),
    SD_difference = sd(Difference, na.rm = TRUE),
    SE_difference = SD_difference / sqrt(N_pairs),
    CI_low = Mean_difference -
      qt(0.975, df = pmax(N_pairs - 1, 1)) * SE_difference,
    CI_high = Mean_difference +
      qt(0.975, df = pmax(N_pairs - 1, 1)) * SE_difference,
    Paired_t_p = ifelse(
      N_pairs >= 2,
      t.test(Difference, mu = 0)$p.value,
      NA_real_
    ),
    Wilcoxon_p = ifelse(
      N_pairs >= 2 && any(abs(Difference) > 0, na.rm = TRUE),
      suppressWarnings(wilcox.test(
        Difference,
        mu = 0,
        exact = FALSE
      )$p.value),
      NA_real_
    ),
    .groups = "drop"
  )

write_csv(
  cor_paired,
  file.path(OUTPUT_DIR, "Heldout_correlation_paired_raw.csv")
)

write_csv(
  cor_paired_summary,
  file.path(
    OUTPUT_DIR,
    "Heldout_correlation_paired_comparison.csv"
  )
)

# Error log
error_log <- raw_results %>% filter(Status != "OK") %>% distinct(Method, n_samples, signal_strength, Replicate, Status)
write_csv(error_log, file.path(OUTPUT_DIR, "benchmark_errors.csv"))

# ----------------------------
# 10. RECOVERY CONSISTENCY / SELECTION FREQUENCY
# ----------------------------
# Because feature names and truth are held fixed across simulations, we can
# calculate how often each feature is recovered in each scenario.

selection_frequency_rows <- list()
idx <- 1
for (i in seq_along(all_selections)) {
  sc <- scenario_grid[i, ]
  sels <- all_selections[[i]]
  if (length(sels) == 0) next

  for (nm in names(sels)) {
    parts <- strsplit(nm, "_")[[1]]
    Method <- parts[1]
    Omic <- paste(parts[-1], collapse = "_")
    feats <- sels[[nm]]

    selection_frequency_rows[[idx]] <- tibble(
      Method = Method,
      Omic = Omic,
      n_samples = sc$n_samples,
      signal_strength = sc$signal_strength,
      Replicate = sc$Replicate,
      Feature = feats
    )
    idx <- idx + 1
  }
}

selection_long <- bind_rows(selection_frequency_rows)
if (nrow(selection_long) > 0) {
  selection_frequency <- selection_long %>%
    group_by(Method, Omic, n_samples, signal_strength, Feature) %>%
    summarise(Selection_Frequency = n_distinct(Replicate) / N_REPS, .groups = "drop") %>%
    mutate(
      Is_True = case_when(
        Omic == "Genes" ~ Feature %in% paste0("Gene_", seq_len(N_LINKED)),
        Omic == "Metabolites" ~ Feature %in% paste0("Met_", seq_len(N_LINKED)),
        TRUE ~ FALSE
      )
    )

  write_csv(selection_frequency, file.path(OUTPUT_DIR, "selection_frequency.csv"))
}

# ----------------------------
# 11. PUBLICATION-QUALITY FIGURES
# ----------------------------
plot_df <- raw_results %>% filter(Status == "OK")

# Figure A: F1 vs signal strength by sample size
p_f1 <- ggplot(plot_df,
               aes(x = factor(signal_strength), y = F1, fill = Method)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
  facet_grid(Omic ~ n_samples, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 12) +
  labs(
    title = "Feature recovery across signal strengths and sample sizes",
    x = "Signal strength",
    y = "F1 score",
    fill = "Method"
  )

ggsave(file.path(OUTPUT_DIR, "Figure_A_F1_feature_recovery.pdf"), p_f1, width = 11, height = 6.5)
ggsave(file.path(OUTPUT_DIR, "Figure_A_F1_feature_recovery.png"), p_f1, width = 11, height = 6.5, dpi = 300)

# Figure B: Precision and Recall
pr_df <- plot_df %>%
  select(Method, Omic, n_samples, signal_strength, Replicate, Precision, Recall) %>%
  pivot_longer(c(Precision, Recall), names_to = "Metric", values_to = "Score")

p_pr <- ggplot(pr_df,
               aes(x = factor(signal_strength), y = Score, fill = Method)) +
  geom_boxplot(outlier.size = 0.35, position = position_dodge(width = 0.8)) +
  facet_grid(Metric + Omic ~ n_samples, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 11) +
  labs(
    title = "Precision and recall of true signal feature recovery",
    x = "Signal strength",
    y = "Score",
    fill = "Method"
  )

ggsave(file.path(OUTPUT_DIR, "Figure_B_precision_recall.pdf"), p_pr, width = 12, height = 8)
ggsave(file.path(OUTPUT_DIR, "Figure_B_precision_recall.png"), p_pr, width = 12, height = 8, dpi = 300)

# Figure C: Held-out canonical correlation (one row per replicate/method; de-duplicate omics)
testcor_df <- plot_df %>%
  distinct(Method, n_samples, signal_strength, Replicate, Test_Canonical_Correlation)

p_cor <- ggplot(testcor_df,
                aes(x = factor(signal_strength), y = Test_Canonical_Correlation, fill = Method)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
  facet_wrap(~n_samples, labeller = label_both) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = 12) +
  labs(
    title = "Generalization of cross-omics canonical structure",
    x = "Signal strength",
    y = "Mean absolute held-out canonical correlation",
    fill = "Method"
  )

ggsave(file.path(OUTPUT_DIR, "Figure_C_heldout_canonical_correlation.pdf"), p_cor, width = 10, height = 4.5)
ggsave(file.path(OUTPUT_DIR, "Figure_C_heldout_canonical_correlation.png"), p_cor, width = 10, height = 4.5, dpi = 300)

# Figure D: Runtime
runtime_df <- plot_df %>%
  distinct(Method, n_samples, signal_strength, Replicate, Runtime_sec)

p_runtime <- ggplot(runtime_df,
                    aes(x = factor(n_samples), y = Runtime_sec, fill = Method)) +
  geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
  theme_classic(base_size = 12) +
  labs(
    title = "Runtime comparison",
    x = "Number of samples",
    y = "Runtime (seconds)",
    fill = "Method"
  )

ggsave(file.path(OUTPUT_DIR, "Figure_D_runtime.pdf"), p_runtime, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "Figure_D_runtime.png"), p_runtime, width = 7, height = 5, dpi = 300)

# Figure E: True vs null selection frequency for a representative scenario
if (exists("selection_frequency") && nrow(selection_frequency) > 0) {
  rep_n <- if (60 %in% SAMPLE_SIZES) 60 else SAMPLE_SIZES[1]
  rep_sig <- if (1.5 %in% SIGNAL_STRENGTHS) 1.5 else SIGNAL_STRENGTHS[1]

  freq_plot_df <- selection_frequency %>%
    filter(n_samples == rep_n, signal_strength == rep_sig) %>%
    mutate(Truth = ifelse(Is_True, "True signal", "Null feature"))

  p_freq <- ggplot(freq_plot_df,
                   aes(x = Truth, y = Selection_Frequency, fill = Method)) +
    geom_boxplot(outlier.size = 0.4, position = position_dodge(width = 0.8)) +
    facet_wrap(~Omic) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Selection consistency (n=", rep_n, ", signal=", rep_sig, ")"),
      x = NULL,
      y = "Selection frequency across simulations",
      fill = "Method"
    )

  ggsave(file.path(OUTPUT_DIR, "Figure_E_selection_frequency.pdf"), p_freq, width = 8, height = 4.5)
  ggsave(file.path(OUTPUT_DIR, "Figure_E_selection_frequency.png"), p_freq, width = 8, height = 4.5, dpi = 300)
}


# Figure F1: Paired F1 differences (OmicNetR - mixOmics)
p_f1_diff <- ggplot(
  f1_paired,
  aes(
    x = factor(signal_strength),
    y = Difference
  )
) +
  geom_boxplot(outlier.size = 0.35) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  facet_grid(
    Omic ~ n_samples,
    labeller = label_both
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Paired difference in feature-recovery F1",
    subtitle = "Positive values favor OmicNetR; negative values favor mixOmics",
    x = "Signal strength",
    y = "F1 difference (OmicNetR - mixOmics)"
  )

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_F1_paired_F1_difference.pdf"
  ),
  p_f1_diff,
  width = 11,
  height = 6.5
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_F1_paired_F1_difference.png"
  ),
  p_f1_diff,
  width = 11,
  height = 6.5,
  dpi = 300
)

# Figure F2: Paired held-out canonical-correlation differences
p_cor_diff <- ggplot(
  cor_paired,
  aes(
    x = factor(signal_strength),
    y = Difference
  )
) +
  geom_boxplot(outlier.size = 0.35) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  facet_wrap(
    ~n_samples,
    labeller = label_both
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Paired difference in held-out canonical correlation",
    subtitle = "Positive values favor OmicNetR; negative values favor mixOmics",
    x = "Signal strength",
    y = "Correlation difference (OmicNetR - mixOmics)"
  )

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_F2_paired_heldout_correlation_difference.pdf"
  ),
  p_cor_diff,
  width = 10,
  height = 4.5
)

ggsave(
  file.path(
    OUTPUT_DIR,
    "Figure_F2_paired_heldout_correlation_difference.png"
  ),
  p_cor_diff,
  width = 10,
  height = 4.5,
  dpi = 300
)

# ----------------------------
# 12. BOOTSTRAP FEATURE-SELECTION STABILITY
# ----------------------------
# Stability is evaluated at n = BOOTSTRAP_N across weak, moderate,
# and strong signal settings defined in BOOTSTRAP_SIGNALS.
#
# IMPORTANT:
# The bootstrap analysis retains the full pairwise Jaccard distribution
# rather than reducing all bootstrap comparisons to a single mean.

jaccard <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}

pairwise_jaccard <- function(selection_list) {

  if (length(selection_list) < 2) {
    return(tibble(
      Bootstrap_1 = integer(),
      Bootstrap_2 = integer(),
      Jaccard = numeric()
    ))
  }

  cmb <- combn(seq_along(selection_list), 2)

  tibble(
    Bootstrap_1 = cmb[1, ],
    Bootstrap_2 = cmb[2, ],
    Jaccard = apply(
      cmb,
      2,
      function(z) {
        jaccard(
          selection_list[[z[1]]],
          selection_list[[z[2]]]
        )
      }
    )
  )
}

bootstrap_stability_one_scenario <- function(
    n_samples = 60,
    signal_strength = 1.0,
    n_boot = 100,
    seed = 777) {

  sim <- simulate_multiomics(
    n_samples = n_samples,
    n_genes = N_GENES,
    n_metabolites = N_METS,
    n_linked = N_LINKED,
    signal_strength = signal_strength,
    noise_sd = 1,
    n_latent = N_COMPONENTS,
    seed = seed
  )

  # Standardize the full base dataset once.
  # Bootstrap resamples are then drawn from the same standardized dataset.
  Xs <- scale(sim$X)
  Ys <- scale(sim$Y)

  om_g <- list()
  om_m <- list()
  mx_g <- list()
  mx_m <- list()

  set.seed(seed + 1)

  for (b in seq_len(n_boot)) {

    idx <- sample(
      seq_len(n_samples),
      size = n_samples,
      replace = TRUE
    )

    Xb <- Xs[idx, , drop = FALSE]
    Yb <- Ys[idx, , drop = FALSE]

    # Give each bootstrap row unique matched identifiers.
    rownames(Xb) <- paste0("Boot", b, "_", seq_len(n_samples))
    rownames(Yb) <- rownames(Xb)

    om <- tryCatch(
      fit_omicnetr(
        Xb, Yb,
        N_COMPONENTS,
        OMIC_PENALTY_X,
        OMIC_PENALTY_Y
      ),
      error = function(e) NULL
    )

    mx <- tryCatch(
      fit_mixomics(
        Xb, Yb,
        N_COMPONENTS,
        SELECTION_K_GENES,
        SELECTION_K_METS
      ),
      error = function(e) NULL
    )

    if (!is.null(om)) {
      om_g[[length(om_g) + 1]] <-
        top_k_features(om$LX, SELECTION_K_GENES)

      om_m[[length(om_m) + 1]] <-
        top_k_features(om$LY, SELECTION_K_METS)
    }

    if (!is.null(mx)) {
      mx_g[[length(mx_g) + 1]] <-
        top_k_features(mx$LX, SELECTION_K_GENES)

      mx_m[[length(mx_m) + 1]] <-
        top_k_features(mx$LY, SELECTION_K_METS)
    }
  }

  bind_rows(

    pairwise_jaccard(om_g) %>%
      mutate(
        Method = "OmicNetR",
        Omic = "Genes"
      ),

    pairwise_jaccard(om_m) %>%
      mutate(
        Method = "OmicNetR",
        Omic = "Metabolites"
      ),

    pairwise_jaccard(mx_g) %>%
      mutate(
        Method = "mixOmics",
        Omic = "Genes"
      ),

    pairwise_jaccard(mx_m) %>%
      mutate(
        Method = "mixOmics",
        Omic = "Metabolites"
      )
  ) %>%
    mutate(
      n_samples = n_samples,
      signal_strength = signal_strength,
      N_boot_requested = n_boot
    )
}

if (RUN_BOOTSTRAP_STABILITY) {

  cat("\nRunning bootstrap feature-selection stability analysis...\n")

  bootstrap_all <- vector(
    "list",
    length(BOOTSTRAP_SIGNALS)
  )

  for (s in seq_along(BOOTSTRAP_SIGNALS)) {

    sig <- BOOTSTRAP_SIGNALS[s]

    cat(
      sprintf(
        "  Bootstrap scenario %d/%d: n=%d, signal=%.1f, B=%d\n",
        s,
        length(BOOTSTRAP_SIGNALS),
        BOOTSTRAP_N,
        sig,
        N_BOOT
      )
    )

    bootstrap_all[[s]] <-
      bootstrap_stability_one_scenario(
        n_samples = BOOTSTRAP_N,
        signal_strength = sig,
        n_boot = N_BOOT,
        seed = MASTER_SEED + 55555 + s * 1000
      )
  }

  stability_results <- bind_rows(bootstrap_all)

  write_csv(
    stability_results,
    file.path(
      OUTPUT_DIR,
      "bootstrap_stability_pairwise.csv"
    )
  )

  stability_summary <- stability_results %>%
    group_by(
      Method,
      Omic,
      n_samples,
      signal_strength
    ) %>%
    summarise(
      N_pairwise_comparisons =
        sum(is.finite(Jaccard)),
      Mean_Jaccard =
        mean(Jaccard, na.rm = TRUE),
      SD_Jaccard =
        sd(Jaccard, na.rm = TRUE),
      Median_Jaccard =
        median(Jaccard, na.rm = TRUE),
      Q1 =
        quantile(Jaccard, 0.25, na.rm = TRUE),
      Q3 =
        quantile(Jaccard, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(
    stability_summary,
    file.path(
      OUTPUT_DIR,
      "bootstrap_stability_summary.csv"
    )
  )

  print(stability_summary)

  # Main stability figure
  p_stability <- ggplot(
    stability_results,
    aes(
      x = Method,
      y = Jaccard,
      fill = Method
    )
  ) +
    geom_boxplot(
      width = 0.60,
      outlier.size = 0.35
    ) +
    facet_grid(
      Omic ~ signal_strength,
      labeller = label_both
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic(base_size = 12) +
    guides(fill = "none") +
    labs(
      title = paste0(
        "Bootstrap feature-selection stability (n = ",
        BOOTSTRAP_N,
        ")"
      ),
      x = NULL,
      y = "Pairwise Jaccard similarity"
    )

  ggsave(
    file.path(
      OUTPUT_DIR,
      "Figure_F_bootstrap_stability.pdf"
    ),
    p_stability,
    width = 10,
    height = 6
  )

  ggsave(
    file.path(
      OUTPUT_DIR,
      "Figure_F_bootstrap_stability.png"
    ),
    p_stability,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# ----------------------------
# 13. SAVE SESSION INFORMATION
# ----------------------------
sink(file.path(OUTPUT_DIR, "sessionInfo.txt"))
print(sessionInfo())
sink()

# ----------------------------
# 14. PRINT MAIN SUMMARY
# ----------------------------
cat("\n==============================\n")
cat("BENCHMARK COMPLETE\n")
cat("==============================\n")
cat("Results saved in:", normalizePath(OUTPUT_DIR), "\n\n")

print(summary_results %>%
        arrange(Omic, n_samples, signal_strength, Method) %>%
        select(Method, Omic, n_samples, signal_strength,
               F1_mean, F1_sd, TestCor_mean, Runtime_median_sec,
               Native_Selected_mean))

cat("\nIMPORTANT INTERPRETATION NOTES:\n")
cat("1. Primary feature-recovery comparison uses the same top-K selection budget for both methods.\n")
cat("2. Native_N_Selected is reported separately because OmicNetR penalty-based sparsity and mixOmics keepX/keepY are parameterized differently.\n")
cat("3. Held-out canonical correlation measures generalization of the learned cross-omics components.\n")
cat("4. Do not claim algorithmic superiority from a single simulation setting; interpret patterns across sample size and signal strength.\n")
cat("5. Paired method differences and 95% confidence intervals are saved for F1 and held-out canonical correlation.\n")
cat("6. Bootstrap stability is evaluated at n=", BOOTSTRAP_N,
    " for signal strengths ", paste(BOOTSTRAP_SIGNALS, collapse=", "),
    " using ", N_BOOT, " bootstrap resamples per scenario.\n", sep="")
cat("7. The loading-product OmicNetR network is not benchmarked as a direct biological interaction network here; edge recovery would be largely determined by selected feature recovery in this simulation design.\n")



# ============================================================
# 15. CREATE FINAL COMBINED MANUSCRIPT FIGURE 2
# A = Feature recovery
# B = Held-out canonical correlation
# C = Bootstrap feature-selection stability
# ============================================================

# Install/load patchwork if needed
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork", repos = "https://cloud.r-project.org")
}

library(patchwork)

# ------------------------------------------------------------
# A. Feature recovery
# ------------------------------------------------------------
p_f1_final <- p_f1 +
  labs(
    title = "Feature recovery",
    x = "Signal strength",
    y = "F1 score",
    fill = "Method"
  ) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.position = "top"
  )

# ------------------------------------------------------------
# B. Held-out canonical correlation
# ------------------------------------------------------------
p_cor_final <- p_cor +
  labs(
    title = "Held-out cross-omics correlation",
    x = "Signal strength",
    y = "Mean absolute canonical correlation",
    fill = "Method"
  ) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.position = "top"
  )

# ------------------------------------------------------------
# C. Bootstrap feature-selection stability
# ------------------------------------------------------------
p_stability_final <- p_stability +
  labs(
    title = "Bootstrap feature-selection stability",
    x = NULL,
    y = "Pairwise Jaccard similarity"
  ) +
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.position = "none"
  )

# ------------------------------------------------------------
# Combine panels vertically
# ------------------------------------------------------------
Figure2_FINAL <-
  p_f1_final /
  p_cor_final /
  p_stability_final +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      size = 18,
      face = "bold"
    )
  )

# Display figure in R/RStudio
print(Figure2_FINAL)

# ------------------------------------------------------------
# Save high-resolution PNG
# ------------------------------------------------------------
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_2_FINAL_OmicNetR_vs_mixOmics.png"
  ),
  plot = Figure2_FINAL,
  width = 12,
  height = 17,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# Save vector PDF
# ------------------------------------------------------------
ggsave(
  filename = file.path(
    OUTPUT_DIR,
    "Figure_2_FINAL_OmicNetR_vs_mixOmics.pdf"
  ),
  plot = Figure2_FINAL,
  width = 12,
  height = 17,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

cat(
  "\nFinal combined Figure 2 saved as:\n",
  file.path(
    OUTPUT_DIR,
    "Figure_2_FINAL_OmicNetR_vs_mixOmics.png"
  ),
  "\n",
  file.path(
    OUTPUT_DIR,
    "Figure_2_FINAL_OmicNetR_vs_mixOmics.pdf"
  ),
  "\n"
)







