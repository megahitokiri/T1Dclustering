# step3a_clustering_rstudio.R
# Interactive Step 3 clustering (clean restart, no GMM).
# - Algorithms: KMeans + PAM
# - Stability: KMeans bootstrap (fast, block-wise)
# - Features per dataset (A, B, C):
#       1) genetic           -> PC1..PC15
#       2) genetic_minclin   -> PCs + age + sex
#       3) genetic_fullclin  -> PCs + ALL numeric clinical variables
# - Carries 'is_outlier' flags from Step 2a to outputs
# - Robust clinical normalization + ID harmonization (uppercase & strip non-alnum)
# - Images: PC1–PC2 cluster scatter, cluster-size bars, silhouette vs K
# - Outputs under: out_three_runs/cluster_runs_min/step3_<algo>_<feature>_<dataset>/
#
# Nothing runs on source(); call run_step3(cfg) from the R console.
#
# Usage example in R:
#   source("step3a_clustering_rstudio.R")
#   cfg <- default_step3_config()
#   cfg$DATASETS <- c("A","B","C")
#   cfg$N_PCS    <- 15
#   cfg$BOOT_N   <- 200       # increase to 1000 later
#   cfg$K_RANGE  <- 2:6
#   run_step3(cfg)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(data.table); library(ggplot2); library(cluster)
})

# =================== CONFIG ================================================
default_step3_config <- function() {
  list(
    BASE_DIR    = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper",
    OUT_DIR     = NULL,
    RUN_DIR     = NULL,
    DAT_DIR     = NULL,
    METHOD      = "irlba",                 # must match Step 2a
    DATASETS    = c("A","B","C"),
    FEATURE_SETS = c("genetic","genetic_minclin","genetic_fullclin"),
    K_RANGE     = 2:6,
    N_PCS       = 15,
    SEED        = 123,
    BOOT_N      = 200,
    OUTLIER_STABILITY_CUTOFF = 0.60
  )
}

resolve_paths <- function(cfg) {
  cfg$OUT_DIR <- file.path(cfg$BASE_DIR, "out_three_runs")
  cfg$RUN_DIR <- file.path(cfg$OUT_DIR, "cluster_runs_min")
  cfg$DAT_DIR <- file.path(cfg$OUT_DIR, "dat")
  dir.create(cfg$RUN_DIR, recursive = TRUE, showWarnings = FALSE)
  cfg
}

# =================== I/O HELPERS ===========================================
scores_path   <- function(cfg, ds) file.path(cfg$RUN_DIR, paste0("step2a_pca_", tolower(cfg$METHOD), "_", ds), "pca_scores.tsv")
outliers_path <- function(cfg, ds) file.path(cfg$RUN_DIR, paste0("step2a_pca_", tolower(cfg$METHOD), "_", ds), "pca_outliers.tsv")
outdir_algo   <- function(cfg, algo, feat, ds) file.path(cfg$RUN_DIR, paste0("step3_", algo, "_", feat, "_", ds))

load_scores <- function(cfg, ds) {
  f <- scores_path(cfg, ds)
  if (!file.exists(f)) stop("Missing pca_scores.tsv for ", ds, " at ", f)
  df <- suppressMessages(readr::read_tsv(f, show_col_types = FALSE))
  if (!"IID" %in% names(df)) stop("pca_scores.tsv for ", ds, " lacks IID column")
  df
}

load_outliers <- function(cfg, ds) {
  f <- outliers_path(cfg, ds)
  if (!file.exists(f)) return(data.frame(IID=character(), is_outlier=logical()))
  dt <- suppressMessages(readr::read_tsv(f, show_col_types = FALSE))
  if (!all(c("IID","is_outlier") %in% names(dt))) return(data.frame(IID=character(), is_outlier=logical()))
  dt[, c("IID","is_outlier")]
}

ensure_sex_numeric <- function(df) {
  if (!"Cluster_sex_num" %in% names(df) || all(is.na(df$Cluster_sex_num))) {
    if ("Cluster_sex" %in% names(df)) {
      sx <- tolower(trimws(as.character(df$Cluster_sex)))
      df$Cluster_sex_num <- ifelse(grepl("^m", sx) | sx == "1", 1,
                                   ifelse(grepl("^f", sx) | sx == "0", 0, NA_real_))
    } else {
      df$Cluster_sex_num <- NA_real_
    }
  }
  df$Cluster_sex_num[!(df$Cluster_sex_num %in% c(0,1))] <- NA_real_
  df
}
# =================== CLINICAL NORMALIZATION ================================
normalize_clinical <- function(clin_df) {
  if (is.null(clin_df) || !nrow(clin_df)) return(NULL)
  df <- clin_df
  ln <- tolower(names(df)); names(df) <- ln
  
  # Subject ID candidates
  id_candidates <- c("cluster_subjectid","subjectid","subject_id","iid","id","sampleid","sample_id","participant","participant_id")
  if (!"cluster_subjectid" %in% names(df)) {
    hit <- intersect(id_candidates, names(df)); if (length(hit)) df$cluster_subjectid <- df[[hit[1]]]
  }
  # Age candidates
  age_candidates <- c("cluster_age","age","age_years","years","yrs")
  if (!"cluster_age" %in% names(df)) {
    hit <- intersect(age_candidates, names(df)); if (length(hit)) df$cluster_age <- suppressWarnings(as.numeric(df[[hit[1]]]))
  }
  # Sex candidates (raw)
  sex_candidates <- c("cluster_sex","sex","gender","biological_sex")
  if (!"cluster_sex" %in% names(df)) {
    hit <- intersect(sex_candidates, names(df)); if (length(hit)) df$cluster_sex <- df[[hit[1]]]
  }
  # Numeric sex (M=1, F=0; others -> NA)
  if (!"cluster_sex_num" %in% names(df)) {
    sx <- df$cluster_sex
    if (!is.null(sx)) {
      if (is.numeric(sx)) {
        df$cluster_sex_num <- as.numeric(sx)
      } else {
        sxs <- tolower(as.character(sx))
        df$cluster_sex_num <- ifelse(grepl("^m", sxs), 1,
                                     ifelse(grepl("^f", sxs), 0, NA_real_))
      }
    }
  }
  df$cluster_sex_num[!(df$cluster_sex_num %in% c(0,1))] <- NA_real_
  
  # Back to expected caps
  names(df) <- sub("^cluster_subjectid$", "Cluster_subjectID", names(df))
  names(df) <- sub("^cluster_age$",       "Cluster_age",       names(df))
  names(df) <- sub("^cluster_sex$",       "Cluster_sex",       names(df))
  names(df) <- sub("^cluster_sex_num$",   "Cluster_sex_num",   names(df))
  df
}

load_clinical <- function(cfg) {
  clin_path <- file.path(cfg$DAT_DIR, "cluster_clinical.csv")
  if (!file.exists(clin_path)) {
    message("No clinical file found: ", clin_path)
    return(NULL)
  }
  z <- suppressMessages(readr::read_csv(clin_path, show_col_types = FALSE))
  
  # Normalize column names or ID handling if you already have a helper for that
  if (exists("normalize_clinical")) {
    z <- normalize_clinical(z)
  }
  # Guarantee presence of Cluster_sex_num
  z <- ensure_sex_numeric(z)
  # Report quick coverage so you can see it worked
  covg <- mean(!is.na(z$Cluster_sex_num))
  message(sprintf("Clinical loaded: N=%d | Cluster_sex_num coverage=%.1f%%",
                  nrow(z), 100 * covg))
  z
}
# =================== FEATURE BUILDING ======================================
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

select_pcs <- function(scores_df, n_pcs) {
  pc_cols <- grep("^PC\\d+$", names(scores_df), value = TRUE)
  if (!length(pc_cols)) stop("No PC columns found.")
  pc_cols <- pc_cols[seq_len(min(n_pcs, length(pc_cols)))]
  list(cols = pc_cols, X = scale(as.matrix(scores_df[, pc_cols, drop = FALSE])))
}

build_features <- function(scores_df, clin_df, n_pcs, feature_set = "genetic") {
  pcs <- select_pcs(scores_df, n_pcs); X <- pcs$X; pc_cols <- pcs$cols
  # add canonical ID columns for robust join
  scores_df$IID_canon <- canon_id(scores_df$IID)
  
  # Use full PCA rows when clinical is unavailable or for genetic-only
  if (feature_set == "genetic" || is.null(clin_df) || !"Cluster_subjectID" %in% names(clin_df)) {
    idx <- seq_len(nrow(scores_df))
    return(list(X = X, used_cols = pc_cols, row_index = idx))
  }
  
  clin_df$Cluster_subjectID_canon <- canon_id(clin_df$Cluster_subjectID)
  
  merged <- dplyr::inner_join(
    scores_df %>% dplyr::select(IID, IID_canon, dplyr::all_of(pc_cols)),
    clin_df,
    by = c("IID_canon" = "Cluster_subjectID_canon")
  )
  
  # which rows (in original scores_df) are used?
  idx <- match(merged$IID, scores_df$IID)
  
  usable <- function(M) {
    is.matrix(M) && nrow(M) >= 5 && ncol(M) >= 2 && is.finite(sum(M)) && (sd(as.numeric(M)) > 0)
  }
  
  if (feature_set == "genetic_minclin") {
    want <- c("Cluster_age","Cluster_sex_num")
    have <- intersect(want, names(merged))
    if (!length(have)) {
      warning("[genetic_minclin] Age/Sex not available — using PCs only.")
      return(list(X = X, used_cols = pc_cols, row_index = seq_len(nrow(scores_df))))
    }
    mm <- merged %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(have), as.numeric)) %>%
      dplyr::select(dplyr::all_of(c(pc_cols, have))) %>%
      as.data.frame()
    keep_mask <- stats::complete.cases(mm)
    mm <- mm[keep_mask, , drop = FALSE]
    idx2 <- idx[keep_mask]
    X2 <- scale(as.matrix(mm))
    if (!usable(X2)) {
      warning(sprintf("[genetic_minclin] Too few usable rows after join (n=%d). Using PCs only.", nrow(X2)))
      return(list(X = X, used_cols = pc_cols, row_index = seq_len(nrow(scores_df))))
    }
    return(list(X = X2, used_cols = colnames(mm), row_index = idx2))
  }
  
  if (feature_set == "genetic_fullclin") {
    num_cols <- names(merged)[sapply(merged, is.numeric)]
    use_cols <- unique(c(pc_cols, num_cols))
    mm <- merged %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(setdiff(use_cols, pc_cols)), as.numeric)) %>%
      dplyr::select(dplyr::all_of(use_cols)) %>%
      as.data.frame()
    keep_mask <- stats::complete.cases(mm)
    mm <- mm[keep_mask, , drop = FALSE]
    idx2 <- idx[keep_mask]
    X2 <- scale(as.matrix(mm))

    if (!usable(X2)) {
      warning(sprintf("[genetic_fullclin] Too few usable rows after join (n=%d). Using PCs only.", nrow(X2)))
      return(list(X = X, used_cols = pc_cols, row_index = seq_len(nrow(scores_df))))
    }
    return(list(X = X2, used_cols = use_cols, row_index = idx2))
  }
  
  stop("Unknown feature_set: ", feature_set)
}

# =================== PLOTTING ==============================================
plot_pc12_clusters <- function(scores_df, clusters, out_png, title = "Clusters") {
  pcs <- grep("^PC\\d+$", names(scores_df), value = TRUE)
  if (!all(c("PC1","PC2") %in% pcs)) return(invisible(NULL))
  df <- scores_df %>% dplyr::select(IID, PC1, PC2) %>% dplyr::mutate(cluster = as.factor(clusters))
  p <- ggplot(df, aes(PC1, PC2, color = cluster)) +
    geom_point(alpha = 0.9, size = 2) +
    theme_minimal(base_size = 12) +
    labs(color = "Cluster", title = title, x = "PC1", y = "PC2")
  ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 300)
}

plot_cluster_sizes <- function(clusters, out_png, title = "Cluster sizes") {
  df <- data.frame(cluster = as.factor(clusters)) %>% count(cluster)
  p <- ggplot(df, aes(cluster, n, fill = cluster)) +
    geom_col() + theme_minimal(base_size = 12) +
    labs(title = title, x = "Cluster", y = "Count") + guides(fill = "none")
  ggsave(out_png, p, width = 6, height = 4.5, dpi = 300)
}

plot_sil_curve <- function(sil_df, out_png, title = "Silhouette vs K") {
  if (is.null(sil_df) || !nrow(sil_df)) return(invisible(NULL))
  p <- ggplot(sil_df, aes(k, avg_sil)) +
    geom_line() + geom_point() + theme_minimal(base_size = 12) +
    labs(title = title, x = "K", y = "Average silhouette")
  ggsave(out_png, p, width = 6.5, height = 4.5, dpi = 300)
}

# =================== ALGORITHMS ============================================
choose_k_silhouette <- function(X, k_range) {
  if (nrow(X) < 5) return(list(best_k = 2, sil_df = data.frame(k = integer(), avg_sil = numeric())))
  D <- try(stats::dist(X), silent = TRUE)
  if (inherits(D, "try-error")) {
    warning("dist() failed; defaulting K=2.")
    return(list(best_k = 2, sil_df = data.frame(k = integer(), avg_sil = numeric())))
  }
  best <- list(k = NA_integer_, sil = -Inf)
  sil_df <- data.frame(k = integer(), avg_sil = numeric())
  for (k in k_range) {
    km <- try(suppressWarnings(stats::kmeans(X, centers = k, nstart = 25)), silent = TRUE)
    if (inherits(km, "try-error")) next
    sil <- try(cluster::silhouette(km$cluster, D), silent = TRUE)
    if (inherits(sil, "try-error")) next
    asw <- mean(sil[, "sil_width"])
    sil_df <- rbind(sil_df, data.frame(k = k, avg_sil = asw))
    if (is.finite(asw) && asw > best$sil) best <- list(k = k, sil = asw)
  }
  if (is.na(best$k)) {
    warning("Silhouette failed to choose K; defaulting K=2.")
    return(list(best_k = 2, sil_df = sil_df))
  }
  list(best_k = best$k, sil_df = sil_df)
}

bootstrap_stability_fast <- function(X, B, k_fixed, seed) {
  set.seed(seed)
  n <- nrow(X)
  co <- matrix(0, n, n)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    Xb  <- X[idx, , drop = FALSE]
    mb  <- stats::kmeans(Xb, centers = k_fixed, nstart = 10)$cluster
    for (c in unique(mb)) {
      idxc <- idx[mb == c]
      co[idxc, idxc] <- co[idxc, idxc] + 1
    }
  }
  co / B
}

run_kmeans <- function(X, scores_df, out_dir, ds, feat_name, k_range, seed, boot_B, outlier_flag, cutoff, row_index = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(row_index)) row_index <- seq_len(nrow(scores_df))
  ssub <- scores_df[row_index, , drop = FALSE]
  oflag <- outlier_flag[row_index]
  
  ck <- choose_k_silhouette(X, k_range)
  best_k <- ck$best_k; sil_df <- ck$sil_df
  if (is.na(best_k)) stop("KMeans failed to determine K for ", ds, " [", feat_name, "]")
  set.seed(seed)
  km <- stats::kmeans(X, centers = best_k, nstart = 50)
  clusters <- km$cluster
  
  plot_pc12_clusters(ssub, clusters,
                     file.path(out_dir, "clusters_pc12.png"),
                     sprintf("KMeans (K=%d) - %s [%s]", best_k, ds, feat_name))
  plot_cluster_sizes(clusters, file.path(out_dir, "cluster_sizes.png"),
                     sprintf("KMeans cluster sizes - %s [%s]", ds, feat_name))
  plot_sil_curve(sil_df, file.path(out_dir, "kmeans_silhouette_vs_k.png"),
                 sprintf("KMeans silhouette vs K - %s [%s]", ds, feat_name))
  
  data.table::fwrite(as.data.table(sil_df), file.path(out_dir, "kmeans_silhouette_summary.tsv"), sep = "\t")
  writeLines(paste0("chosen_K\t", best_k), con = file.path(out_dir, "kmeans_chosen_K.txt"))
  
  message("    Bootstrapping KMeans (B=", boot_B, ", K=", best_k, ") …")
  co <- bootstrap_stability_fast(X, B = boot_B, k_fixed = best_k, seed = seed + 1L)
  stab <- vapply(seq_len(nrow(X)), function(i) {
    mates <- which(clusters == clusters[i]); mean(co[i, mates])
  }, numeric(1))
  
  propose <- oflag & (stab < cutoff)
  cl_tab <- data.frame(IID = ssub$IID,
                       cluster = clusters,
                       is_outlier = oflag,
                       stability = round(stab, 4),
                       propose_subgroup = propose)
  data.table::fwrite(as.data.table(cl_tab), file.path(out_dir, "clusters.tsv"), sep = "\t")
  
  #Graph Variant Stability
  # ---- Stability visuals & summaries -----------------------------------------
  # 1) Histogram of stability (overall)
  if ("stability" %in% names(cl_tab)) {
    p_hist <- ggplot(cl_tab, aes(stability)) +
      geom_histogram(bins = 30) +
      theme_minimal(base_size = 12) +
      labs(title = sprintf("KMeans bootstrap stability - %s [%s]", ds, feat_name),
           x = "Stability (0–1)", y = "Count")
    ggsave(file.path(out_dir, "kmeans_stability_hist.png"), p_hist, width = 7, height = 4.5, dpi = 300)
  }
  
  # 2) Violin/box by cluster
  if (all(c("stability","cluster") %in% names(cl_tab))) {
    p_violin <- ggplot(cl_tab, aes(x = factor(cluster), y = stability, fill = factor(cluster))) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.size = 0.6) +
      theme_minimal(base_size = 12) +
      labs(title = sprintf("KMeans stability by cluster - %s [%s]", ds, feat_name),
           x = "Cluster", y = "Bootstrap stability") +
      guides(fill = "none")
    ggsave(file.path(out_dir, "kmeans_stability_by_cluster.png"), p_violin, width = 7.5, height = 5, dpi = 300)
  }
  
  # 3) PC1–PC2 scatter colored by stability (for a quick map of weak points)
  if (all(c("PC1","PC2") %in% names(ssub)) && "stability" %in% names(cl_tab)) {
    pcstab <- cbind(ssub[, c("IID","PC1","PC2")], cl_tab[, c("stability","cluster","propose_subgroup")])
    p_grad <- ggplot(pcstab, aes(PC1, PC2, color = stability)) +
      geom_point(size = 2, alpha = 0.9) +
      theme_minimal(base_size = 12) +
      labs(title = sprintf("KMeans PC1–PC2 colored by stability - %s [%s]", ds, feat_name),
           color = "Stability")
    ggsave(file.path(out_dir, "kmeans_pc12_colored_by_stability.png"), p_grad, width = 7.5, height = 5.5, dpi = 300)
  }
  
  # 4) Summary table per cluster + overall
  if ("stability" %in% names(cl_tab)) {
    sum_all <- data.frame(
      level = "overall",
      n = nrow(cl_tab),
      mean = mean(cl_tab$stability, na.rm = TRUE),
      median = median(cl_tab$stability, na.rm = TRUE),
      sd = sd(cl_tab$stability, na.rm = TRUE),
      pct_lt_0_6 = mean(cl_tab$stability < 0.6, na.rm = TRUE)
    )
    sum_by <- cl_tab |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean = mean(stability, na.rm = TRUE),
        median = median(stability, na.rm = TRUE),
        sd = sd(stability, na.rm = TRUE),
        pct_lt_0_6 = mean(stability < 0.6, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(level = paste0("cluster_", cluster)) |>
      dplyr::select(level, dplyr::everything(), -cluster)
    
    stab_summary <- rbind(sum_all, as.data.frame(sum_by))
    data.table::fwrite(stab_summary, file.path(out_dir, "kmeans_stability_summary.tsv"), sep = "\t")
  }
  # ---------------------------------------------------------------------------
  
  if (any(propose)) {
    data.table::fwrite(as.data.table(cl_tab[cl_tab$propose_subgroup, c("IID")]),
                       file.path(dirname(out_dir), paste0("step3_outliers_proposed_subgroup_", feat_name, "_", ds, ".tsv")),
                       sep = "\t")
  }
  
  list(k = best_k, clusters = clusters, sil = sil_df, table = cl_tab)
}

run_pam <- function(X, scores_df, out_dir, ds, feat_name, k_range, outlier_flag, row_index = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(row_index)) row_index <- seq_len(nrow(scores_df))
  ssub <- scores_df[row_index, , drop = FALSE]
  oflag <- outlier_flag[row_index]
  
  D <- stats::dist(X)
  sil_summary <- list(); pam_models <- list()
  for (k in k_range) {
    pamfit <- cluster::pam(D, k = k)
    pam_models[[as.character(k)]] <- pamfit
    sil <- cluster::silhouette(pamfit$clustering, D)
    sil_summary[[as.character(k)]] <- mean(sil[, "sil_width"])
  }
  sil_df <- data.frame(k = k_range, avg_sil = unlist(sil_summary))
  best_k <- sil_df$k[which.max(sil_df$avg_sil)]
  pam_best <- pam_models[[as.character(best_k)]]
  clusters <- pam_best$clustering
  
  plot_pc12_clusters(ssub, clusters,
                     file.path(out_dir, "clusters_pc12.png"),
                     sprintf("PAM (K=%d) - %s [%s]", best_k, ds, feat_name))
  plot_cluster_sizes(clusters, file.path(out_dir, "cluster_sizes.png"),
                     sprintf("PAM cluster sizes - %s [%s]", ds, feat_name))
  plot_sil_curve(sil_df, file.path(out_dir, "pam_silhouette_vs_k.png"),
                 sprintf("PAM silhouette vs K - %s [%s]", ds, feat_name))
  
  data.table::fwrite(as.data.table(sil_df), file.path(out_dir, "pam_silhouette_summary.tsv"), sep = "\t")
  writeLines(paste0("chosen_K\t", best_k), con = file.path(out_dir, "pam_chosen_K.txt"))
  data.table::fwrite(as.data.table(data.frame(IID = ssub$IID, cluster = clusters, is_outlier = oflag)),
                     file.path(out_dir, "clusters.tsv"), sep = "\t")
  
  list(k = best_k, clusters = clusters, sil = sil_df)
}

# =================== MAIN DRIVER ===========================================
run_step3 <- function(cfg = default_step3_config()) {
  cfg <- resolve_paths(cfg)
  set.seed(cfg$SEED)
  clin <- load_clinical(cfg)
  message("Starting Step 3 (KMeans+PAM, no GMM) for: ", paste(cfg$DATASETS, collapse = ", "))
  for (ds in cfg$DATASETS) {
    message("\n=== Dataset ", ds, " ===")
    scores_df <- load_scores(cfg, ds)
    outs_df   <- load_outliers(cfg, ds)
    scores_df <- scores_df %>% dplyr::left_join(outs_df, by = "IID")
    scores_df$is_outlier[is.na(scores_df$is_outlier)] <- FALSE
    
    message("  IDs in PCA: ", length(unique(scores_df$IID)),
            " | IDs in clinical: ", if (!is.null(clin)) length(unique(clin$Cluster_subjectID)) else 0)
    
    for (feat in cfg$FEATURE_SETS) {
      message("  -- Feature set: ", feat)
      feats <- build_features(scores_df, clin, cfg$N_PCS, feat)
      X <- feats$X; idx <- feats$row_index; is_out <- scores_df$is_outlier
      
      # KMeans
      out_dir_km <- outdir_algo(cfg, "kmeans", feat, ds)
      km <- run_kmeans(X, scores_df, out_dir_km, ds, feat,
                       cfg$K_RANGE, cfg$SEED, cfg$BOOT_N, is_out, cfg$OUTLIER_STABILITY_CUTOFF,
                       row_index = idx)
      
      # PAM
      out_dir_pam <- outdir_algo(cfg, "pam", feat, ds)
      pm <- run_pam(X, scores_df, out_dir_pam, ds, feat, cfg$K_RANGE, is_out,
                    row_index = idx)
    }
  }
  invisible(cfg)
}

# End of file
