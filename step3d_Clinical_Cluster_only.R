## ==================== CLINICAL-ONLY KMEANS GRID ====================
## - Features: numeric clinical columns only (no PCs)
## - Algorithms: KMeans (+ bootstrap stability)
## - Outputs per bucket: grid_Bxxxxx/step3_kmeans_clinical_only_<DS>/
## ================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(glue)
})

## ---------- Config ----------
BASE_DIR   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper"
RUN_ROOT   <- file.path(BASE_DIR, "out_three_runs", "cluster_runs_min")
DAT_DIR    <- file.path(BASE_DIR, "out_three_runs", "dat")
METHOD     <- "irlba"  # step2a method used to get PCA scores (for IDs and optional plotting only)

BOOT_GRID  <- c(1, 10, 100, 250, 500, 1000, 1500, 2000, 2500, 3000)
DATASETS   <- c("A","B","C")
K_RANGE    <- 2:6
SEED       <- 123
OUTLIER_STABILITY_CUTOFF <- 0.60   # only used for flagging/plots

## ---------- Paths ----------
scores_path   <- function(ds) file.path(RUN_ROOT, sprintf("step2a_pca_%s_%s", tolower(METHOD), ds), "pca_scores.tsv")

## ---------- Helpers ----------
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

read_try <- function(fp, ...) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(fread(fp, ...), error = function(e) NULL)
}

# Impute numeric NA with column median (in place)
impute_median_inplace <- function(DT) {
  for (nm in names(DT)) {
    if (is.numeric(DT[[nm]]) && anyNA(DT[[nm]])) {
      med <- suppressWarnings(median(DT[[nm]], na.rm = TRUE))
      if (!is.finite(med)) med <- 0
      set(DT, which(is.na(DT[[nm]])), nm, med)
    }
  }
}

## ---------- Clinical loader & normalizer ----------
load_clinical <- function() {
  fp <- file.path(DAT_DIR, "cluster_clinical.csv")
  stopifnot(file.exists(fp))
  df <- suppressMessages(readr::read_csv(fp, show_col_types = FALSE))
  setDT(df)
  
  # Canonicalize likely ID/sex/age fields we saw in your project
  # (Keep originals; we’ll use the unified names below)
  ln <- tolower(names(df))
  setnames(df, names(df), ln)
  
  # ID -> Cluster_subjectID
  if (!"cluster_subjectid" %in% names(df)) {
    id_cands <- c("subjectid","subject_id","iid","id","sample","participant","cluster_subjectid")
    hit <- id_cands[id_cands %in% names(df)][1]
    if (!is.na(hit)) setnames(df, hit, "cluster_subjectid")
  }
  
  # age -> Cluster_age (numeric)
  if (!"cluster_age" %in% names(df)) {
    age_cands <- c("age","age_years","years","yrs","cluster_age")
    hit <- age_cands[age_cands %in% names(df)][1]
    if (!is.na(hit)) setnames(df, hit, "cluster_age")
  }
  if ("cluster_age" %in% names(df)) {
    df[, cluster_age := suppressWarnings(as.numeric(cluster_age))]
  }
  
  # sex -> Cluster_sex_num (0/1), keep original Cluster_sex if present
  if (!"cluster_sex" %in% names(df)) {
    sex_cands <- c("sex","gender","biological_sex","cluster_sex")
    hit <- sex_cands[sex_cands %in% names(df)][1]
    if (!is.na(hit)) setnames(df, hit, "cluster_sex")
  }
  if (!"cluster_sex_num" %in% names(df)) {
    if ("cluster_sex" %in% names(df)) {
      sx <- tolower(trimws(as.character(df$cluster_sex)))
      df[, cluster_sex_num := fifelse(grepl("^m", sx) | sx == "1", 1,
                                      fifelse(grepl("^f", sx) | sx == "0", 0, as.numeric(NA)))]
    } else {
      df[, cluster_sex_num := as.numeric(NA)]
    }
  }
  df[!(cluster_sex_num %in% c(0,1)), cluster_sex_num := NA_real_]
  
  # Restore expected caps and add canon ID
  setnames(df, old = "cluster_subjectid", new = "Cluster_subjectID", skip_absent = TRUE)
  setnames(df, old = "cluster_age",       new = "Cluster_age",       skip_absent = TRUE)
  setnames(df, old = "cluster_sex",       new = "Cluster_sex",       skip_absent = TRUE)
  setnames(df, old = "cluster_sex_num",   new = "Cluster_sex_num",   skip_absent = TRUE)
  
  df[, Cluster_subjectID_canon := canon_id(Cluster_subjectID)]
  df
}

## ---------- Build clinical-only matrix (no PCs) ----------
build_clinical_matrix <- function(clin_df, keep_ids = NULL, impute = TRUE) {
  if (is.null(clin_df) || !nrow(clin_df)) return(NULL)
  DT <- as.data.table(clin_df)
  
  # Optional: restrict to IDs that appear in the PCA scores (aligns sample universe)
  if (!is.null(keep_ids)) {
    DT <- DT[Cluster_subjectID_canon %in% keep_ids]
  }
  
  if (!nrow(DT)) return(NULL)
  
  # Force a few likely numeric fields to numeric if present
  force_num <- intersect(
    c("Cluster_age","Cluster_sex_num","Cluster_height","Cluster_weight","Cluster_insulin","Cluster_A1c"),
    names(DT)
  )
  for (v in force_num) set(DT, j = v, value = suppressWarnings(as.numeric(DT[[v]])))
  
  # Select numeric clinical columns only (exclude IDs)
  num_cols <- names(DT)[vapply(DT, is.numeric, logical(1))]
  use_cols <- setdiff(num_cols, c("Cluster_subjectID","Cluster_subjectID_canon"))
  
  if (!length(use_cols)) return(NULL)
  
  # Subset with .. (NO bare symbols)
  M <- DT[, ..use_cols]
  
  # Drop all-NA and zero-variance columns
  keep_col <- vapply(M, function(v) any(is.finite(v)) && (sd(v, na.rm = TRUE) > 0), logical(1))
  if (!any(keep_col)) return(NULL)
  M <- M[, which(keep_col), with = FALSE]
  use_cols <- names(M)
  
  # Impute medians (or use complete.cases)
  if (isTRUE(impute)) {
    impute_median_inplace(M)
  } else {
    cc <- stats::complete.cases(M)
    M  <- M[cc]
    DT <- DT[cc]
  }
  
  if (nrow(M) < 5 || ncol(M) < 2) return(NULL)
  
  X <- tryCatch(scale(as.matrix(M)), error = function(e) NULL)
  if (is.null(X) || !is.finite(sum(X))) return(NULL)
  
  list(
    X = X,
    row_ids = DT$Cluster_subjectID_canon,
    used_cols = names(M)
  )
}

## ---------- KMeans + stability ----------
choose_k_silhouette <- function(X, k_range) {
  if (nrow(X) < 5) return(list(best_k = 2, sil_df = data.frame(k = integer(), avg_sil = numeric())))
  D <- try(stats::dist(X), silent = TRUE)
  if (inherits(D, "try-error")) return(list(best_k = 2, sil_df = data.frame(k = integer(), avg_sil = numeric())))
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
  if (is.na(best$k)) best$k <- 2
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

## ---------- Plotters ----------
plot_pc12_clusters <- function(scores_df, row_ids, clusters, out_png, title) {
  if (!all(c("PC1","PC2") %in% names(scores_df))) return(invisible(NULL))
  ssub <- scores_df %>% mutate(IID_canon = canon_id(IID)) %>%
    filter(IID_canon %in% row_ids) %>%
    select(IID, PC1, PC2)
  if (!nrow(ssub)) return(invisible(NULL))
  df <- ssub %>% mutate(cluster = as.factor(clusters))
  p <- ggplot(df, aes(PC1, PC2, color = cluster)) +
    geom_point(alpha = 0.9, size = 2) +
    theme_minimal(base_size = 12) +
    labs(color = "Cluster", title = title, x = "PC1", y = "PC2")
  ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 300)
}

plot_cluster_sizes <- function(clusters, out_png, title) {
  df <- data.frame(cluster = as.factor(clusters)) %>% count(cluster)
  p <- ggplot(df, aes(cluster, n, fill = cluster)) +
    geom_col() + theme_minimal(base_size = 12) +
    labs(title = title, x = "Cluster", y = "Count") + guides(fill = "none")
  ggsave(out_png, p, width = 6, height = 4.5, dpi = 300)
}

plot_sil_curve <- function(sil_df, out_png, title) {
  if (is.null(sil_df) || !nrow(sil_df)) return(invisible(NULL))
  p <- ggplot(sil_df, aes(k, avg_sil)) +
    geom_line() + geom_point() + theme_minimal(base_size = 12) +
    labs(title = title, x = "K", y = "Average silhouette")
  ggsave(out_png, p, width = 6.5, height = 4.5, dpi = 300)
}

## ---------- One dataset runner ----------
run_one_dataset <- function(ds, B) {
  out_dir <- file.path(RUN_ROOT, sprintf("grid_B%05d", B), sprintf("step3_kmeans_clinical_only_%s", ds))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # load clinical
  clin <- load_clinical()
  message(sprintf("[%s] Clinical loaded: N=%d | Cluster_sex_num coverage=%.1f%%",
                  ds, nrow(clin), 100*mean(!is.na(clin$Cluster_sex_num))))
  
  # load PCA scores ONLY to align IDs + for plotting PC scatter (not as features)
  sc <- read_try(scores_path(ds))
  if (is.null(sc) || !"IID" %in% names(sc)) stop("Missing or invalid pca_scores.tsv for ", ds)
  sc[, IID_canon := canon_id(IID)]
  
  # build clinical matrix aligned to available IDs in this dataset
  cm <- build_clinical_matrix(clin, keep_ids = unique(sc$IID_canon), impute = TRUE)
  if (is.null(cm)) stop("[", ds, "] clinical-only features unusable (too few rows/cols).")
  
  X        <- cm$X
  row_ids  <- cm$row_ids
  used_cols <- cm$used_cols
  message(sprintf("[%s] clinical-only usable rows: %d; cols: %d → %s%s",
                  ds, nrow(X), ncol(X), paste(head(used_cols, 8), collapse = ", "),
                  if (length(used_cols) > 8) ", (… more)" else ""))
  
  # choose K by silhouette, fit KMeans
  set.seed(SEED)
  ck <- choose_k_silhouette(X, K_RANGE)
  best_k <- ck$best_k
  km <- stats::kmeans(X, centers = best_k, nstart = 50)
  clusters <- km$cluster
  
  # bootstrap stability
  co <- bootstrap_stability_fast(X, B = B, k_fixed = best_k, seed = SEED + 1L)
  stab <- vapply(seq_len(nrow(X)), function(i) {
    mates <- which(clusters == clusters[i]); mean(co[i, mates])
  }, numeric(1))
  
  # outputs
  fwrite(as.data.table(ck$sil_df), file.path(out_dir, "kmeans_silhouette_summary.tsv"), sep = "\t")
  writeLines(paste0("chosen_K\t", best_k), con = file.path(out_dir, "kmeans_chosen_K.txt"))
  
  cl_tab <- data.table(
    IID = row_ids,                  # canon ID here; fine for clinical-only
    cluster = clusters,
    stability = round(stab, 4)
  )
  fwrite(cl_tab, file.path(out_dir, "clusters.tsv"), sep = "\t")
  
  # visuals (reuse your style)
  plot_pc12_clusters(sc, row_ids, clusters,
                     file.path(out_dir, "clusters_pc12.png"),
                     sprintf("KMeans (K=%d) - %s [clinical_only]", best_k, ds))
  plot_cluster_sizes(clusters, file.path(out_dir, "cluster_sizes.png"),
                     sprintf("KMeans cluster sizes - %s [clinical_only]", ds))
  plot_sil_curve(ck$sil_df, file.path(out_dir, "kmeans_silhouette_vs_k.png"),
                 sprintf("KMeans silhouette vs K - %s [clinical_only]", ds))
  
  # stability hist
  p_hist <- ggplot(as.data.frame(cl_tab), aes(stability)) +
    geom_histogram(bins = 30) + theme_minimal(base_size = 12) +
    labs(title = sprintf("KMeans bootstrap stability - %s [clinical_only]", ds),
         x = "Stability (0–1)", y = "Count")
  ggsave(file.path(out_dir, "kmeans_stability_hist.png"), p_hist, width = 7, height = 4.5, dpi = 300)
  
  # violin by cluster
  p_vio <- ggplot(as.data.frame(cl_tab), aes(x = factor(cluster), y = stability, fill = factor(cluster))) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.size = 0.6) +
    theme_minimal(base_size = 12) +
    labs(title = sprintf("KMeans stability by cluster - %s [clinical_only]", ds),
         x = "Cluster", y = "Bootstrap stability") +
    guides(fill = "none")
  ggsave(file.path(out_dir, "kmeans_stability_by_cluster.png"), p_vio, width = 7.5, height = 5, dpi = 300)
  
  # summary table
  sum_all <- data.frame(
    level = "overall",
    n = nrow(cl_tab),
    mean = mean(cl_tab$stability, na.rm = TRUE),
    median = median(cl_tab$stability, na.rm = TRUE),
    sd = sd(cl_tab$stability, na.rm = TRUE),
    pct_lt_0_6 = mean(cl_tab$stability < OUTLIER_STABILITY_CUTOFF, na.rm = TRUE)
  )
  sum_by <- cl_tab |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(stability, na.rm = TRUE),
      median = median(stability, na.rm = TRUE),
      sd = sd(stability, na.rm = TRUE),
      pct_lt_0_6 = mean(stability < OUTLIER_STABILITY_CUTOFF, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(level = paste0("cluster_", cluster)) |>
    dplyr::select(level, dplyr::everything(), -cluster)
  stab_summary <- rbind(sum_all, as.data.frame(sum_by))
  fwrite(stab_summary, file.path(out_dir, "kmeans_stability_summary.tsv"), sep = "\t")
}

## ---------- Grid runner ----------
run_grid_clinical_only <- function() {
  for (B in BOOT_GRID) {
    bucket <- file.path(RUN_ROOT, sprintf("grid_B%05d", B))
    dir.create(bucket, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("\n[Clinical-only] BOOT_N=%d → %s\n", B, bucket))
    
    # (Optionally) cleanup previous clinical-only dirs in this bucket
    old <- list.dirs(bucket, full.names = TRUE, recursive = FALSE)
    old <- old[grepl("^step3_kmeans_clinical_only_", basename(old))]
    if (length(old)) unlink(old, recursive = TRUE, force = TRUE)
    cat(sprintf("  - cleaned old clinical-only dirs: %d\n", length(old)))
    
    for (ds in DATASETS) {
      cat(sprintf(">>> Dataset %s\n", ds))
      tryCatch(
        run_one_dataset(ds, B),
        error = function(e) cat(sprintf("  ! ERROR at BOOT_N=%d / DS=%s: %s\n", B, ds, conditionMessage(e)))
      )
    }
    
    # sanity: count clusters.tsv written by this script
    n_clusters <- length(list.files(bucket,
                                    pattern = "^clusters\\.tsv$",
                                    recursive = TRUE, full.names = TRUE))
    cat(sprintf("  - clusters.tsv in bucket: %d\n", n_clusters))
  }
  cat("\nDone.\n")
}

## ---------- Go! ----------
set.seed(SEED)
run_grid_clinical_only()