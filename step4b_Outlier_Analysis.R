## =========================
##  Population-level viz of clusters + outliers
## =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(glue)
  library(readr)
})

## ---- Paths you already use ----
BASE <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
CLIN_PATH <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/dat/cluster_clinical.csv"

## Datasets & features to pull
DATASETS <- c("A","B","C")
FEATURES <- c("genetic","genetic_minclin","genetic_fullclin")

## (1) Choose which bucket:
##    - NULL = auto-pick the newest grid_* that actually has clusters
##    - or set e.g. "grid_B05000"
GRID_TAG <- "grid_B02500"

## ---------------- helpers ----------------
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

find_latest_grid <- function(base_dir) {
  grids <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  grids <- grids[grepl("/grid_B\\d{5}$", grids)]
  if (!length(grids)) return(NA_character_)
  mt <- file.info(grids)$mtime
  grids[order(mt, decreasing = TRUE)][1]
}

bucket_has_clusters <- function(grid_dir) {
  length(list.files(grid_dir, pattern = "^clusters\\.tsv$", recursive = TRUE)) > 0
}

load_pca <- function(base_dir, ds) {
  # prefer irlba if you have both
  f <- file.path(base_dir, sprintf("step2a_pca_irlba_%s", ds), "pca_scores.tsv")
  if (!file.exists(f)) {
    f <- file.path(base_dir, sprintf("step2a_pca_prcomp_%s", ds), "pca_scores.tsv")
  }
  stopifnot(file.exists(f))
  dt <- suppressMessages(readr::read_tsv(f, show_col_types = FALSE))
  idcol <- c("IID","id","sample","SubjectID")
  idcol <- idcol[idcol %in% names(dt)][1]
  if (is.na(idcol)) stop("No ID column found in PCA scores for dataset ", ds)
  dt <- as.data.frame(dt)
  dt$IID <- canon_id(dt[[idcol]])
  dt[, c("IID","PC1","PC2")]
}

load_clusters_one <- function(grid_dir, ds, feat) {
  # /grid_Bxxxxx/step3_kmeans_<DS>/<feature>/clusters.tsv
  d <- file.path(grid_dir, sprintf("step3_kmeans_genetic_%s", ds))
  # if your folder names are exactly "step3_kmeans_<feature>_<DS>" use this:
  d <- file.path(grid_dir, sprintf("step3_kmeans_%s_%s", feat, ds))
  f <- file.path(d, "clusters.tsv")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  # standardize expected cols
  if (!"IID" %in% names(dt)) {
    idcol <- c("IID","id","sample","SubjectID")
    idcol <- idcol[idcol %in% names(dt)][1]
    if (is.na(idcol)) return(NULL)
    dt$IID <- dt[[idcol]]
  }
  dt$IID <- canon_id(dt$IID)
  if (!"cluster" %in% names(dt)) stop("clusters.tsv missing 'cluster'")
  if (!"is_outlier" %in% names(dt)) dt$is_outlier <- FALSE
  dt$dataset <- ds
  dt$feature <- feat
  dt[, c("IID","cluster","is_outlier","dataset","feature")]
}

load_all_clusters <- function(grid_dir, datasets, features) {
  out <- list()
  for (ds in datasets) {
    for (ft in features) {
      z <- load_clusters_one(grid_dir, ds, ft)
      if (!is.null(z)) out[[paste(ds,ft,sep=":")]] <- z
    }
  }
  if (!length(out)) stop("No clusters.tsv found under: ", grid_dir)
  bind_rows(out)
}

## ---------------- pick bucket ----------------
if (is.null(GRID_TAG)) {
  cand <- find_latest_grid(BASE)
  if (is.na(cand)) stop("No grid_Bxxxxx folder found in: ", BASE)
  if (!bucket_has_clusters(cand)) {
    stop("Latest bucket found but has no clusters.tsv: ", cand)
  }
  GRID_DIR <- cand
} else {
  GRID_DIR <- file.path(BASE, GRID_TAG)
  if (!dir.exists(GRID_DIR)) stop("Bucket not found: ", GRID_DIR)
  if (!bucket_has_clusters(GRID_DIR)) stop("Bucket has no clusters.tsv: ", GRID_DIR)
}

message("Using bucket: ", GRID_DIR)

## ---------------- load data ----------------
# clinical
clin <- suppressMessages(readr::read_csv(CLIN_PATH, show_col_types = FALSE))
id_col <- c("Cluster_subjectID","IID","subject_id","id")
id_col <- id_col[id_col %in% names(clin)][1]
if (is.na(id_col)) stop("No ID column in clinical CSV")
clin <- clin %>% rename(SubjectID = !!id_col)
clin$SubjectID <- canon_id(clin$SubjectID)

# keep a few clinical variables for summaries (add/remove as you like)
CLIN_VARS <- c("Cluster_age","Cluster_A1c","Cluster_weight","Cluster_height","Cluster_sex_num")
have_vars <- intersect(CLIN_VARS, names(clin))
clin_small <- clin %>% select(SubjectID, all_of(have_vars))

# clusters
cl_all <- load_all_clusters(GRID_DIR, DATASETS, FEATURES)

# pca for all datasets
pca_all <- bind_rows(lapply(DATASETS, function(ds) {
  df <- load_pca(BASE, ds)
  df$dataset <- ds
  df
}))

# merge PCA + clusters + clinical
merged <- cl_all %>%
  left_join(pca_all, by = c("IID","dataset")) %>%
  left_join(clin_small, by = c("IID" = "SubjectID"))

# coerce types
merged$cluster <- as.factor(merged$cluster)
merged$is_outlier <- as.logical(merged$is_outlier)
merged$feature <- factor(merged$feature, levels = FEATURES)

## ---------------- output dir ----------------
OUT_DIR <- file.path(GRID_DIR, "viz_summary")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ================= PLOT 1: PCA overlay, outliers highlighted =================
p1 <- ggplot(merged, aes(PC1, PC2, color = cluster, shape = is_outlier)) +
  geom_point(alpha = 0.8, size = 1.8, stroke = 0.6) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c("Inlier","Outlier")) +
  facet_grid(dataset ~ feature) +
  labs(title = "PCA (PC1–PC2) by cluster; Outliers shown as ×",
       color = "Cluster", shape = "") +
  theme_minimal(base_size = 13)
ggsave(file.path(OUT_DIR, "pca_overlay_outliers.png"), p1, width = 12, height = 7, dpi = 200)

## ================= PLOT 2: Outlier rate per cluster =================
outlier_rate <- merged %>%
  group_by(dataset, feature, cluster) %>%
  summarise(outlier_rate = mean(is_outlier, na.rm = TRUE),
            n = dplyr::n(), .groups = "drop")

p2 <- ggplot(outlier_rate, aes(x = cluster, y = outlier_rate, fill = cluster)) +
  geom_col() +
  geom_text(aes(label = scales::percent(outlier_rate, accuracy = 0.1)), vjust = -0.3, size = 3.3) +
  facet_grid(dataset ~ feature) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Outlier fraction by cluster",
       x = "Cluster", y = "Outliers (%)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "outlier_fraction_by_cluster.png"), p2, width = 12, height = 7, dpi = 200)

## ================= PLOT 3: Median clinical values per cluster =================
if (length(have_vars)) {
  clin_long <- merged %>%
    select(dataset, feature, cluster, all_of(have_vars)) %>%
    pivot_longer(cols = all_of(have_vars), names_to = "variable", values_to = "value")
  
  clin_sum <- clin_long %>%
    group_by(dataset, feature, cluster, variable) %>%
    summarise(median_value = median(value, na.rm = TRUE),
              n_non_na = sum(!is.na(value)),
              .groups = "drop")
  
  p3 <- ggplot(clin_sum, aes(x = cluster, y = median_value, fill = cluster)) +
    geom_col(position = "dodge") +
    facet_grid(variable ~ dataset + feature, scales = "free_y") +
    labs(title = "Median clinical values per cluster",
         x = "Cluster", y = "Median") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(OUT_DIR, "clinical_medians_by_cluster.png"), p3, width = 13, height = 9, dpi = 200)
} else {
  message("No clinical variables found among: ", paste(CLIN_VARS, collapse=", "))
}

## ================= PLOT 4 (optional): Cluster sizes =================
sizes <- merged %>%
  group_by(dataset, feature, cluster) %>%
  summarise(n = dplyr::n(), .groups = "drop")

p4 <- ggplot(sizes, aes(x = cluster, y = n, fill = cluster)) +
  geom_col() +
  facet_grid(dataset ~ feature) +
  labs(title = "Cluster sizes", x = "Cluster", y = "Count") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "cluster_sizes.png"), p4, width = 12, height = 7, dpi = 200)

## ---------------- done ----------------
cat(glue("\nSaved figures in: {OUT_DIR}\n",
         " - pca_overlay_outliers.png\n",
         " - outlier_fraction_by_cluster.png\n",
         " - clinical_medians_by_cluster.png\n",
         " - cluster_sizes.png\n"))