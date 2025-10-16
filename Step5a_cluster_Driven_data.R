## ==============================================================
## B02500 summary: clinical-only cluster sizes, stability,
## and alignment with genetic / minclin / fullclin clusters
## Outputs -> grid_B02500/viz_summary/
## ==============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glue)
})

## ---- Paths ----
ROOT <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
BUCKET <- file.path(ROOT, "grid_B02500")
VIZDIR <- file.path(BUCKET, "viz_summary")
dir.create(VIZDIR, recursive = TRUE, showWarnings = FALSE)

## ---- Helpers ----
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# Add a robust JoinID (= canonized IID/SubjectID/Cluster_subjectID) and standard columns
standardize_clusters <- function(p) {
  if (!file.exists(p)) return(NULL)
  dt <- tryCatch(fread(p), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  
  # find ID column
  idc <- c("IID","SubjectID","Cluster_subjectID","id","sample")
  idc <- idc[idc %in% names(dt)][1]
  if (is.na(idc)) return(NULL)
  
  setnames(dt, idc, "ID")
  dt[, JoinID := canon_id(ID)]
  
  # feature/dataset from path (expects .../step3_kmeans_<feature>_<dataset>/clusters.tsv)
  parent <- basename(dirname(p))
  m <- regexec("^step3_kmeans_(.+)_([ABC])$", parent)
  mm <- regmatches(parent, m)[[1]]
  if (length(mm) != 3) return(NULL)
  dt[, feature := mm[2] ]
  dt[, dataset := mm[3] ]
  
  # ensure cluster & stability columns exist (stability may be absent)
  if (!"cluster" %in% names(dt)) {
    cand <- c("kmeans_cluster","pam_cluster","gmm_cluster","Cluster","clusters","label")
    hit  <- cand[cand %in% names(dt)][1]
    if (!is.na(hit)) setnames(dt, hit, "cluster")
  }
  if (!"cluster" %in% names(dt)) return(NULL)
  
  if (!"stability" %in% names(dt)) dt[, stability := NA_real_]
  
  # keep only needed
  keep <- intersect(c("JoinID","cluster","stability","dataset","feature"), names(dt))
  out <- unique(dt[, ..keep])
  out[, cluster := as.factor(cluster)]
  out
}

# Adjusted Rand Index (no external deps)
adj_rand_index <- function(x, y) {
  x <- as.integer(as.factor(x)); y <- as.integer(as.factor(y))
  n <- length(x)
  if (n == 0) return(NA_real_)
  
  tab <- table(x, y)
  a <- rowSums(tab); b <- colSums(tab)
  c2 <- function(v) sum(v * (v - 1)) / 2
  sum_ij <- c2(tab)
  sum_i  <- c2(a)
  sum_j  <- c2(b)
  tot    <- c2(n)
  
  expected <- (sum_i * sum_j) / tot
  max_idx <- 0.5 * (sum_i + sum_j)
  if (tot == expected && max_idx == expected) return(1)  # identical
  (sum_ij - expected) / (max_idx - expected)
}

## ---- Load all clusters.tsv from BUCKET and keep only desired features ----
all_files <- list.files(BUCKET, pattern = "^clusters\\.tsv$", recursive = TRUE, full.names = TRUE)
all_dt <- rbindlist(lapply(all_files, standardize_clusters), fill = TRUE)
all_dt <- all_dt[feature %in% c("genetic","genetic_minclin","genetic_fullclin","clinical_only")]

# Ensure we have something
stopifnot(nrow(all_dt) > 0)

## ---- 1) Clinical-only sizes & stability summaries ----
clin <- all_dt[feature == "clinical_only"]

# size per dataset/cluster
clin_sizes <- clin[, .(N = .N), by = .(dataset, cluster)][order(dataset, as.integer(as.character(cluster)))]
fwrite(clin_sizes, file.path(VIZDIR, "B02500_clinical_only_sizes.tsv"), sep = "\t")

# stability by cluster + overall
clin_stab_overall <- clin[, .(
  n = .N,
  mean = mean(stability, na.rm = TRUE),
  median = median(stability, na.rm = TRUE),
  sd = sd(stability, na.rm = TRUE),
  pct_lt_0_6 = mean(stability < 0.6, na.rm = TRUE)
), by = .(dataset)]

clin_stab_byclu <- clin[, .(
  n = .N,
  mean = mean(stability, na.rm = TRUE),
  median = median(stability, na.rm = TRUE),
  sd = sd(stability, na.rm = TRUE),
  pct_lt_0_6 = mean(stability < 0.6, na.rm = TRUE)
), by = .(dataset, cluster)]

fwrite(clin_stab_overall, file.path(VIZDIR, "B02500_clinical_only_stability_overall.tsv"), sep = "\t")
fwrite(clin_stab_byclu,  file.path(VIZDIR, "B02500_clinical_only_stability_bycluster.tsv"), sep = "\t")

## ---- 2) Alignment of clinical-only vs genetic/minclin/fullclin (per dataset) ----
# small availability diagnostic
avail <- all_dt[, .N, by = .(dataset, feature)][order(dataset, feature)]
fwrite(avail, file.path(VIZDIR, "B02500_available_dataset_feature_counts.tsv"), sep = "\t")

compare_feat <- function(ds, feat_ref = "clinical_only", feat_cmp = "genetic") {
  ref <- all_dt[dataset == ds & feature == feat_ref, .(JoinID, clu_ref = cluster)]
  cmp <- all_dt[dataset == ds & feature == feat_cmp, .(JoinID, clu_cmp = cluster)]
  mm  <- merge(ref, cmp, by = "JoinID", all = FALSE)
  
  # Always emit a row so downstream code has ARI column to work with
  if (!nrow(mm)) {
    return(data.table(
      dataset   = ds, ref = feat_ref, cmp = feat_cmp,
      N_overlap = 0L, ARI = NA_real_
    ))
  }
  
  # contingency (for inspection)
  tab_wide <- as.data.table(as.data.frame.matrix(table(mm$clu_ref, mm$clu_cmp)))
  setnames(tab_wide, names(tab_wide), paste0("cmp_", names(tab_wide)))
  tab_wide[, ref_cluster := rownames(as.data.frame(tab_wide))]
  setcolorder(tab_wide, c("ref_cluster", setdiff(names(tab_wide), "ref_cluster")))
  fwrite(tab_wide, file.path(VIZDIR, glue("B02500_contingency_{ds}_{feat_ref}_vs_{feat_cmp}.tsv")), sep = "\t")
  
  # ARI
  ari <- adj_rand_index(mm$clu_ref, mm$clu_cmp)
  
  # Heatmap (row-normalized) with counts printed
  ct <- as.data.frame(table(mm$clu_ref, mm$clu_cmp))
  names(ct) <- c("ref","cmp","n")
  ct <- data.table(ct)
  ct[, prop := n / sum(n), by = .(ref)]
  p <- ggplot(ct, aes(x = cmp, y = ref, fill = prop)) +
    geom_tile() +
    geom_text(aes(label = n), color = "white", size = 3) +
    scale_fill_gradient(low = "gray90", high = "steelblue") +
    labs(title = glue("{ds}: {feat_ref} vs {feat_cmp} (ARI={sprintf('%.3f', ari)})"),
         x = feat_cmp, y = feat_ref, fill = "row prop") +
    theme_minimal(base_size = 12)
  ggsave(file.path(VIZDIR, glue("B02500_heatmap_{ds}_{feat_ref}_vs_{feat_cmp}.png")),
         p, width = 6.5, height = 4.8, dpi = 160)
  
  data.table(dataset = ds, ref = feat_ref, cmp = feat_cmp, N_overlap = nrow(mm), ARI = ari)
}

datasets <- sort(unique(all_dt$dataset))
cmp_list <- lapply(datasets, function(ds) rbindlist(list(
  compare_feat(ds, "clinical_only", "genetic"),
  compare_feat(ds, "clinical_only", "genetic_minclin"),
  compare_feat(ds, "clinical_only", "genetic_fullclin")
), fill = TRUE))
pairs <- rbindlist(cmp_list, fill = TRUE)

# Write raw pairs (even if ARI is NA everywhere)
fwrite(pairs, file.path(VIZDIR, "B02500_alignment_clinical_only_vs_others.tsv"), sep = "\t")

## ---- 3) Alignment call table (skip NA ARIs, but don’t crash if all NA) ----
pairs_valid <- pairs[!is.na(ARI)]
if (nrow(pairs_valid)) {
  pairs_valid[, call := fifelse(ARI < 0.2, "clinical-distinct",
                                fifelse(ARI < 0.5, "mixed", "aligned-with-genetics"))]
  setorder(pairs_valid, dataset, cmp)
  fwrite(pairs_valid, file.path(VIZDIR, "B02500_alignment_call.tsv"), sep = "\t")
} else {
  # Emit an empty skeleton with header so downstream usage doesn't fail
  fwrite(data.table(dataset=character(), ref=character(), cmp=character(),
                    N_overlap=integer(), ARI=double(), call=character()),
         file.path(VIZDIR, "B02500_alignment_call.tsv"), sep = "\t")
}
## ---- 4) Console hints ----
cat(glue("
✔ Wrote (B02500):
  - {file.path(VIZDIR, 'B02500_clinical_only_sizes.tsv')}
  - {file.path(VIZDIR, 'B02500_clinical_only_stability_overall.tsv')}
  - {file.path(VIZDIR, 'B02500_clinical_only_stability_bycluster.tsv')}
  - {file.path(VIZDIR, 'B02500_alignment_clinical_only_vs_others.tsv')}
  - {file.path(VIZDIR, 'B02500_alignment_call.tsv')}
  - Heatmaps + contingency tables per dataset under: {VIZDIR}
"))


## ============================
## Visuals for clinical-only run
## ============================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glue)
  library(mclust)     # for adjustedRandIndex
  library(scales)
})

## ---- Paths (change only BASE/BUCKET if needed) ----
BASE   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
BUCKET <- file.path(BASE, "grid_B02500")
VIZDIR <- file.path(BUCKET, "viz_summary")
dir.create(VIZDIR, recursive = TRUE, showWarnings = FALSE)

DATASETS <- c("A","B","C")
FEATS    <- c("genetic","genetic_minclin","genetic_fullclin")
CLIN     <- "clinical_only"

canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

## Read a clusters.tsv as a tidy data.table with standardized ID/cols
read_clusters <- function(bucket, algo = "kmeans", feat, ds) {
  d <- file.path(bucket, sprintf("step3_%s_%s_%s", algo, feat, ds))
  f <- file.path(d, "clusters.tsv")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  # ID column harmonization
  idc <- c("IID","SubjectID","id","sample","Cluster_subjectID")
  idc <- idc[idc %in% names(dt)][1]
  if (is.na(idc)) return(NULL)
  setnames(dt, idc, "IID")
  dt[, IID := canon_id(IID)]
  # cluster column harmonization
  if (!"cluster" %in% names(dt)) {
    cand <- c("kmeans_cluster","pam_cluster","gmm_cluster","Cluster","clusters","label")
    hit  <- cand[cand %in% names(dt)][1]
    if (!is.na(hit)) setnames(dt, hit, "cluster")
  }
  if (!"cluster" %in% names(dt)) return(NULL)
  # stability (optional)
  if (!"stability" %in% names(dt)) dt[, stability := NA_real_]
  dt[, `:=`(dataset = ds, feature = feat, algo = "kmeans")]
  dt[]
}

## ---------- 1) Clinical-only cluster sizes (per dataset) ----------
sizes_list <- lapply(DATASETS, function(ds) read_clusters(BUCKET, "kmeans", CLIN, ds))
sizes_list <- Filter(Negate(is.null), sizes_list)

if (length(sizes_list)) {
  clin_all <- rbindlist(sizes_list, use.names = TRUE, fill = TRUE)
  p_sizes <- clin_all[, .N, by = .(dataset, cluster)][
    , cluster := factor(cluster)
  ][
    , ggplot(.SD, aes(cluster, N, fill = cluster)) +
      geom_col(width = 0.85) +
      geom_text(aes(label = N), vjust = -0.25, size = 3) +
      facet_wrap(~ dataset, scales = "free_y") +
      labs(title = "Clinical-only cluster sizes", x = "Cluster", y = "Count") +
      guides(fill = "none") +
      theme_minimal(base_size = 12)
  ]
  ggsave(file.path(VIZDIR, "clinical_only_cluster_sizes.png"), p_sizes, width = 10, height = 4.5, dpi = 160)
} else {
  message("No clinical-only clusters found under ", BUCKET)
}

## ---------- 2) Clinical-only stability by cluster ----------
if (length(sizes_list)) {
  stab_plot <- clin_all[!is.na(stability)]
  if (nrow(stab_plot)) {
    p_stab <- ggplot(stab_plot, aes(x = factor(cluster), y = stability, fill = factor(cluster))) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.size = 0.6) +
      facet_wrap(~ dataset) +
      labs(title = "Clinical-only bootstrap stability by cluster", x = "Cluster", y = "Stability") +
      guides(fill = "none") +
      theme_minimal(base_size = 12)
    ggsave(file.path(VIZDIR, "clinical_only_stability_by_cluster.png"), p_stab, width = 10, height = 4.5, dpi = 160)
  }
}

## ---------- 3) ARI: clinical_only vs each genetic feature (per dataset) ----------
compare_one <- function(ds, feat_other) {
  c1 <- read_clusters(BUCKET, "kmeans", CLIN, ds)
  c2 <- read_clusters(BUCKET, "kmeans", feat_other, ds)
  if (is.null(c1) || is.null(c2)) return(NULL)
  m <- merge(c1[, .(IID, cl_c = cluster)], c2[, .(IID, cl_o = cluster)], by = "IID")
  if (!nrow(m)) return(NULL)
  ari <- mclust::adjustedRandIndex(as.integer(as.factor(m$cl_c)),
                                   as.integer(as.factor(m$cl_o)))
  data.table(dataset = ds, compare = feat_other, ARI = ari, N_overlap = nrow(m))
}

ari_list <- lapply(DATASETS, function(ds) lapply(FEATS, function(ft) compare_one(ds, ft)))
ari_list <- Filter(Negate(is.null), unlist(ari_list, recursive = FALSE))
if (length(ari_list)) {
  ari_tab <- rbindlist(ari_list)
  fwrite(ari_tab, file.path(VIZDIR, "clinical_only_vs_genetic_ARI.tsv"), sep = "\t")
  # Heatmap-style plot
  ari_tab[, compare := factor(compare, levels = c("genetic","genetic_minclin","genetic_fullclin"))]
  p_ari <- ggplot(ari_tab, aes(compare, dataset, fill = ARI)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", ARI))) +
    scale_fill_gradient(low = "#fbe7e7", high = "#5c8ef0", limits = c(0, 1), oob = squish) +
    labs(title = "Alignment (ARI) of clinical-only clusters vs genetic feature sets",
         x = "Feature set", y = "Dataset", fill = "ARI") +
    theme_minimal(base_size = 12)
  ggsave(file.path(VIZDIR, "clinical_only_vs_genetic_ARI_heatmap.png"), p_ari, width = 8, height = 3.8, dpi = 160)
}

## ---------- 4) Confusion heatmaps: clinical_only vs genetic_fullclin ----------
conf_plot_one <- function(ds, feat_other = "genetic_fullclin") {
  c1 <- read_clusters(BUCKET, "kmeans", CLIN, ds)
  c2 <- read_clusters(BUCKET, "kmeans", feat_other, ds)
  if (is.null(c1) || is.null(c2)) return(NULL)
  m <- merge(c1[, .(IID, cl_c = factor(cluster))], c2[, .(IID, cl_o = factor(cluster))], by = "IID")
  if (!nrow(m)) return(NULL)
  tab <- as.data.table(table(m$cl_c, m$cl_o))
  setnames(tab, c("V1","V2","N"), c("ClinicalOnly","Other","N"))
  # normalize by clinical-only row
  tab[, RowTotal := sum(N), by = ClinicalOnly]
  tab[, Prop := ifelse(RowTotal > 0, N / RowTotal, 0)]
  p <- ggplot(tab, aes(Other, ClinicalOnly, fill = Prop)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(N, " (", percent(Prop, accuracy = 0.1), ")")), size = 3) +
    scale_fill_gradient(low = "#fbe7e7", high = "#5c8ef0", limits = c(0, 1), oob = squish) +
    labs(title = glue("Confusion: clinical-only vs {feat_other} — Dataset {ds}"),
         x = feat_other, y = "clinical_only", fill = "Row %") +
    theme_minimal(base_size = 12)
  out <- file.path(VIZDIR, glue("confusion_clinical_only_vs_{feat_other}_{ds}.png"))
  ggsave(out, p, width = 6.8, height = 4.8, dpi = 160)
  data.table(dataset = ds, file = out)
}

conf_done <- rbindlist(lapply(DATASETS, conf_plot_one), fill = TRUE, use.names = TRUE)

cat(glue("
✔ Saved to: {VIZDIR}
  - clinical_only_cluster_sizes.png
  - clinical_only_stability_by_cluster.png (if stability present)
  - clinical_only_vs_genetic_ARI.tsv
  - clinical_only_vs_genetic_ARI_heatmap.png
  - confusion_clinical_only_vs_genetic_fullclin_[A|B|C].png (where available)
"))