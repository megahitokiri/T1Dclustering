#!/usr/bin/env R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

# =================== USER SETTINGS ===================
# Root where step2a + step3 outputs live
BASE_ROOT <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"

# Pick the bootstrap bucket you want to analyze (e.g. "grid_B02500")
GRID_TAG  <- "grid_B02500"   # <-- change to grid_B05000, etc.

# Clinical file (already cleaned with Cluster_sex_num = 0/1)
CLIN_PATH <- file.path(BASE_ROOT, "..", "dat", "cluster_clinical.csv")

# Datasets & feature families to include
DATASETS  <- c("A","B","C")
FEATURES  <- c("genetic","genetic_minclin","genetic_fullclin")

# How many PCs to look for in the PCA score files (the script adapts)
N_PCS     <- 15

# Output folder for all figures/tables
OUT_DIR   <- file.path(BASE_ROOT, GRID_TAG, "viz_summary")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# =================== HELPERS =========================
p_safe <- function(...) {
  p <- try(ggsave(...), silent = TRUE)
  if (inherits(p, "try-error")) message("ggsave failed for: ", list(...)[[1]])
}

# Find a path if it exists
maybe <- function(...) {
  p <- file.path(...)
  if (file.exists(p)) p else NA_character_
}

# Load clusters (KMeans or PAM). Returns: IID, cluster, is_outlier, stability?, propose_subgroup?
load_clusters <- function(run_dir) {
  fp <- maybe(run_dir, "clusters.tsv")
  if (is.na(fp)) return(NULL)
  fread(fp, showProgress = FALSE) %>%
    mutate(IID = canon_id(IID))
}

# Load KMeans summaries if present
load_kmeans_summaries <- function(run_dir) {
  sil <- maybe(run_dir, "kmeans_silhouette_summary.tsv")
  stab <- maybe(run_dir, "kmeans_stability_summary.tsv")
  list(
    sil = if (!is.na(sil)) fread(sil) else NULL,
    stab = if (!is.na(stab)) fread(stab) else NULL
  )
}

# Load PCA scores from step2a
load_pca_scores <- function(ds) {
  # Prefer irlba; fallback to prcomp
  p1 <- maybe(BASE_ROOT, sprintf("step2a_pca_irlba_%s", ds), "pca_scores.tsv")
  p2 <- maybe(BASE_ROOT, sprintf("step2a_pca_prcomp_%s", ds), "pca_scores.tsv")
  fp <- if (!is.na(p1)) p1 else if (!is.na(p2)) p2 else NA_character_
  if (is.na(fp)) return(NULL)
  sc <- fread(fp)
  # Ensure ID present
  idcol <- intersect(c("IID","id","sample","SubjectID"), names(sc))
  stopifnot(length(idcol) >= 1)
  sc <- sc %>%
    rename(IID = !!idcol[1]) %>%
    mutate(IID = canon_id(IID))
  sc
}

# Load PCA top loadings (PC1 & PC2) if available
load_top_loadings <- function(ds) {
  p1 <- maybe(BASE_ROOT, sprintf("step2a_pca_irlba_%s", ds), "pca_top_loadings_PC1_PC2.tsv")
  p2 <- maybe(BASE_ROOT, sprintf("step2a_pca_prcomp_%s", ds), "pca_top_loadings_PC1_PC2.tsv")
  fp <- if (!is.na(p1)) p1 else if (!is.na(p2)) p2 else NA_character_
  if (is.na(fp)) return(NULL)
  fread(fp)
}

# Robust clinical loader
load_clinical <- function() {
  if (!file.exists(CLIN_PATH)) return(NULL)
  clin <- suppressMessages(readr::read_csv(CLIN_PATH, show_col_types = FALSE))
  # Standardize ID column
  idc <- intersect(c("SubjectID","Cluster_subjectID","IID","id"), names(clin))
  if (!length(idc)) return(NULL)
  clin <- clin %>% rename(SubjectID = !!idc[1])
  clin$SubjectID <- canon_id(clin$SubjectID)
  # Ensure numeric sex
  if (!("Cluster_sex_num" %in% names(clin))) {
    sx <- tolower(as.character(clin$Cluster_sex))
    clin$Cluster_sex_num <- ifelse(grepl("^m|^1$", sx), 1,
                                   ifelse(grepl("^f|^0$", sx), 0, NA_real_))
  }
  clin
}

# Build full path for step3 results
run_path <- function(alg, feat, ds) {
  file.path(BASE_ROOT, GRID_TAG, sprintf("step3_%s_%s_%s", alg, feat, ds))
}

# =================== LOAD GLOBAL DATA =================
clin <- load_clinical()

# Collect everything into one long table for convenience
rows <- list()
kmeans_sil_all <- list()
kmeans_stab_all <- list()

for (ds in DATASETS) {
  pca <- load_pca_scores(ds)
  if (is.null(pca)) { message("No PCA for ", ds); next }
  pca_cols <- grep("^PC\\d+$", names(pca), value = TRUE)
  keep_pcs <- pca_cols[seq_len(min(N_PCS, length(pca_cols)))]
  pca_small <- pca %>% select(IID, all_of(keep_pcs))
  
  for (feat in FEATURES) {
    # KMeans
    km_dir <- run_path("kmeans", feat, ds)
    km_cl  <- load_clusters(km_dir)
    # PAM
    pm_dir <- run_path("pam", feat, ds)
    pm_cl  <- load_clusters(pm_dir)
    
    if (!is.null(km_cl)) {
      ksum <- load_kmeans_summaries(km_dir)
      if (!is.null(ksum$sil)) {
        ksum$sil$dataset <- ds; ksum$sil$feature <- feat
        kmeans_sil_all[[length(kmeans_sil_all)+1]] <- ksum$sil
      }
      if (!is.null(ksum$stab)) {
        ksum$stab$dataset <- ds; ksum$stab$feature <- feat
        kmeans_stab_all[[length(kmeans_stab_all)+1]] <- ksum$stab
      }
    }
    
    # Merge PCA + clusters + clinical
    merged <- pca_small
    if (!is.null(km_cl)) merged <- merged %>% left_join(km_cl %>% select(IID, cluster_km = cluster, is_outlier, stability), by = "IID")
    if (!is.null(pm_cl)) merged <- merged %>% left_join(pm_cl %>% select(IID, cluster_pam = cluster), by = "IID")
    if (!is.null(clin))  merged <- merged %>% left_join(clin, by = c("IID" = "SubjectID"))
    
    merged$dataset <- ds
    merged$feature <- feat
    rows[[length(rows)+1]] <- merged
  }
}
big <- bind_rows(rows)

# Save a merged TSV for ad-hoc inspection
readr::write_tsv(big, file.path(OUT_DIR, "merged_pca_clusters_clinical.tsv"))

# =================== PLOTS ===================

# 1) PCA overlay: KMeans color, PAM shape; outliers circled if present
plot_pca_overlay <- function(df, ds, feat) {
  pcs <- intersect(c("PC1","PC2"), names(df))
  if (length(pcs) < 2) return(NULL)
  dd <- df %>% filter(dataset == ds, feature == feat)
  if (!nrow(dd)) return(NULL)
  dd$is_outlier <- isTRUE(dd$is_outlier) | (!is.na(dd$is_outlier) & dd$is_outlier)
  p <- ggplot(dd, aes(.data[[pcs[1]]], .data[[pcs[2]]])) +
    geom_point(aes(color = factor(cluster_km), shape = factor(cluster_pam),
                   alpha = ifelse(is_outlier, 0.9, 0.6)), size = 2) +
    scale_alpha_identity() +
    labs(title = sprintf("PCA PC1 vs PC2 — %s [%s]", ds, feat),
         x = "PC1", y = "PC2", color = "KMeans", shape = "PAM") +
    theme_minimal(base_size = 13)
  p_safe(filename = file.path(OUT_DIR, sprintf("pca_overlay_%s_%s.png", ds, feat)),
         plot = p, width = 6, height = 5, dpi = 300)
  p
}

for (ds in DATASETS) for (feat in FEATURES) plot_pca_overlay(big, ds, feat)

# 2) KMeans stability (if present): histogram + by-cluster box
if (length(kmeans_stab_all)) {
  stab_all <- bind_rows(kmeans_stab_all)
  # overall rows + cluster-specific rows are both in file; keep per-sample if present
  # But typical file has per-sample 'stability' exported elsewhere; here: plot aggregate summaries.
  # If your per-sample stability is in clusters.tsv, use that (done below).
  
  # Per-sample stability from clusters.tsv (KMeans only)
  stab_pts <- big %>% filter(!is.na(stability), !is.na(cluster_km))
  if (nrow(stab_pts)) {
    # Hist
    p <- ggplot(stab_pts, aes(stability)) +
      geom_histogram(bins = 40) +
      facet_grid(dataset ~ feature) +
      labs(title = "KMeans per-sample stability", x = "Stability", y = "Count") +
      theme_minimal(base_size = 12)
    p_safe(file.path(OUT_DIR, "kmeans_stability_hist.png"), p, 10, 8, 300)
    
    # By cluster
    p2 <- ggplot(stab_pts, aes(x = factor(cluster_km), y = stability, fill = factor(cluster_km))) +
      geom_violin(trim = FALSE, alpha = 0.6) + geom_boxplot(width = 0.2, outlier.size = 0.6) +
      facet_grid(dataset ~ feature) +
      labs(title = "KMeans stability by cluster", x = "Cluster", y = "Stability") +
      theme_minimal(base_size = 12) + theme(legend.position = "none")
    p_safe(file.path(OUT_DIR, "kmeans_stability_by_cluster.png"), p2, 10, 8, 300)
    
    # PCA colored by stability
    if (all(c("PC1","PC2") %in% names(stab_pts))) {
      p3 <- ggplot(stab_pts, aes(PC1, PC2)) +
        geom_point(aes(color = stability), size = 1.8) +
        facet_grid(dataset ~ feature) +
        labs(title = "PCA colored by KMeans stability", color = "Stability") +
        theme_minimal(base_size = 12)
      p_safe(file.path(OUT_DIR, "kmeans_pc12_colored_by_stability.png"), p3, 10, 8, 300)
    }
  }
}

# 3) Silhouette vs K (if present)
if (length(kmeans_sil_all)) {
  sil_all <- bind_rows(kmeans_sil_all)
  if (all(c("K","silhouette") %in% names(sil_all))) {
    p <- ggplot(sil_all, aes(x = K, y = silhouette, group = feature, color = feature)) +
      geom_line(size = 1.1) + geom_point(size = 2) +
      facet_grid(dataset ~ ., scales = "free_y") +
      labs(title = "KMeans average silhouette across K", x = "K", y = "Avg silhouette", color = "Feature") +
      theme_minimal(base_size = 12)
    p_safe(file.path(OUT_DIR, "kmeans_silhouette_vs_K.png"), p, 7, 9, 300)
  }
}

# 4) Cluster sizes (KMeans & PAM)
plot_cluster_sizes <- function(label_col, file_tag) {
  dd <- big %>% filter(!is.na(.data[[label_col]]))
  if (!nrow(dd)) return()
  p <- dd %>%
    count(dataset, feature, cluster = .data[[label_col]]) %>%
    ggplot(aes(x = factor(cluster), y = n, fill = factor(cluster))) +
    geom_col() +
    facet_grid(dataset ~ feature, scales = "free_y") +
    labs(title = paste("Cluster sizes —", toupper(file_tag)),
         x = "Cluster", y = "Count") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  p_safe(file.path(OUT_DIR, paste0("cluster_sizes_", file_tag, ".png")), p, 10, 8, 300)
}
plot_cluster_sizes("cluster_km", "kmeans")
plot_cluster_sizes("cluster_pam", "pam")

# 5) Clinical distributions by cluster (KMeans)
clin_vars <- c("Cluster_age","Cluster_A1c","Cluster_weight","Cluster_height","Cluster_insulin","Cluster_sex_num")
have_clin <- intersect(clin_vars, names(big))
if (length(have_clin) && any(!is.na(big$cluster_km))) {
  long <- big %>%
    select(dataset, feature, cluster_km, all_of(have_clin)) %>%
    pivot_longer(cols = all_of(have_clin), names_to = "variable", values_to = "value")
  
  # Numeric vars as boxplots
  num_vars <- setdiff(have_clin, c("Cluster_sex_num"))
  if (length(num_vars)) {
    p <- long %>%
      filter(variable %in% num_vars) %>%
      ggplot(aes(x = factor(cluster_km), y = value, fill = factor(cluster_km))) +
      geom_boxplot(outlier.size = 0.6) +
      facet_grid(variable ~ dataset + feature, scales = "free_y") +
      labs(title = "Clinical distributions by KMeans cluster", x = "Cluster", y = "Value") +
      theme_minimal(base_size = 11) + theme(legend.position = "none")
    p_safe(file.path(OUT_DIR, "clinical_boxplots_by_cluster_kmeans.png"), p, 12, 9, 300)
  }
  
  # Sex proportion
  if ("Cluster_sex_num" %in% have_clin) {
    p2 <- big %>%
      filter(!is.na(cluster_km), !is.na(Cluster_sex_num)) %>%
      mutate(sex = ifelse(Cluster_sex_num == 1, "Male", "Female")) %>%
      count(dataset, feature, cluster_km, sex) %>%
      group_by(dataset, feature, cluster_km) %>%
      mutate(prop = n / sum(n)) %>% ungroup() %>%
      ggplot(aes(x = factor(cluster_km), y = prop, fill = sex)) +
      geom_col(position = "fill") +
      facet_grid(dataset ~ feature) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(title = "Sex composition by KMeans cluster", x = "Cluster", y = "Proportion", fill = "Sex") +
      theme_minimal(base_size = 12)
    p_safe(file.path(OUT_DIR, "sex_ratio_by_cluster_kmeans.png"), p2, 10, 8, 300)
  }
}

# 6) Cluster centroids in PC space
if (all(c("PC1","PC2") %in% names(big))) {
  cent <- big %>%
    filter(!is.na(cluster_km)) %>%
    group_by(dataset, feature, cluster_km) %>%
    summarise(PC1c = mean(PC1, na.rm=TRUE), PC2c = mean(PC2, na.rm=TRUE), .groups="drop")
  p <- ggplot(big %>% filter(!is.na(cluster_km)), aes(PC1, PC2, color = factor(cluster_km))) +
    geom_point(alpha = 0.35, size = 1.5) +
    geom_point(data = cent, aes(x = PC1c, y = PC2c), color = "black", shape = 8, size = 3, inherit.aes = FALSE) +
    facet_grid(dataset ~ feature) +
    labs(title = "Cluster centroids (KMeans) in PCA space", color = "Cluster") +
    theme_minimal(base_size = 12)
  p_safe(file.path(OUT_DIR, "pca_centroids_kmeans.png"), p, 10, 8, 300)
}

# =================== "MOST INFLUENTIAL SNPs" ===================
# Idea: quantify how much PC1/PC2 separate clusters; weight SNP loadings by that separation.
# We only have top loadings for PC1/PC2 from step2a; use those if present.
snps_reports <- list()

for (ds in DATASETS) {
  # 1) Cluster separation along PCs
  dd <- big %>% filter(dataset == ds, feature == "genetic", !is.na(cluster_km))
  if (!nrow(dd) || !all(c("PC1","PC2") %in% names(dd))) next
  sep_PC1 <- abs(diff(tapply(dd$PC1, dd$cluster_km, mean, na.rm = TRUE)))
  sep_PC2 <- abs(diff(tapply(dd$PC2, dd$cluster_km, mean, na.rm = TRUE)))
  # scale to [0,1] weights (avoid both zero)
  w1 <- if ((sep_PC1 + sep_PC2) > 0) sep_PC1/(sep_PC1 + sep_PC2) else 0.5
  w2 <- if ((sep_PC1 + sep_PC2) > 0) sep_PC2/(sep_PC1 + sep_PC2) else 0.5
  
  # 2) Load top loadings PC1+PC2
  loads <- load_top_loadings(ds)
  if (is.null(loads)) next
  # Expect columns: SNP, loading, abs_loading, PC ("PC1"/"PC2")
  top1 <- loads %>% filter(PC == "PC1") %>% select(SNP, loading_PC1 = loading, abs_PC1 = abs_loading)
  top2 <- loads %>% filter(PC == "PC2") %>% select(SNP, loading_PC2 = loading, abs_PC2 = abs_loading)
  top  <- full_join(top1, top2, by = "SNP") %>% mutate(across(everything(), ~replace_na(., 0)))
  
  # 3) Score SNPs: weighted sum of absolute loadings by PC separation weights
  top <- top %>%
    mutate(score = w1 * abs(loading_PC1) + w2 * abs(loading_PC2)) %>%
    arrange(desc(score))
  
  # Save per-dataset table
  out_tsv <- file.path(OUT_DIR, sprintf("top_snps_driving_clusters_%s.tsv", ds))
  data.table::fwrite(top, out_tsv, sep = "\t")
  
  # Quick bar plot of top 20
  plot_df <- head(top, 20) %>% mutate(SNP = factor(SNP, levels = rev(SNP)))
  p <- ggplot(plot_df, aes(SNP, score)) +
    geom_col() + coord_flip() +
    labs(title = sprintf("Top SNPs driving KMeans split — %s (PC1/PC2 loadings weighted)", ds),
         x = "SNP", y = "Score") +
    theme_minimal(base_size = 12)
  p_safe(file.path(OUT_DIR, sprintf("top_snps_bar_%s.png", ds)), p, 7, 6, 300)
  
  # Collect a summary row
  snps_reports[[length(snps_reports)+1]] <- data.frame(
    dataset = ds, sep_PC1 = sep_PC1, sep_PC2 = sep_PC2,
    w_PC1 = w1, w_PC2 = w2,
    top_SNP = plot_df$SNP[1],
    top_score = plot_df$score[1]
  )
}

if (length(snps_reports)) {
  snp_sum <- bind_rows(snps_reports)
  readr::write_tsv(snp_sum, file.path(OUT_DIR, "top_snps_summary_across_datasets.tsv"))
}

# =================== DONE ===================
message("All visuals and tables written to: ", OUT_DIR)