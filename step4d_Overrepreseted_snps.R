# snps_enriched_in_small_cluster.R
# Goal: identify SNPs over-represented (by direction of PCA loadings) in the
#       smaller, more stable KMeans cluster for each dataset/feature.
# Outputs:
#   - enriched_snps_<dataset>_<feature>.tsv
#   - enriched_snps_top20_<dataset>_<feature>.png

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
})

# ===================== USER CONFIG =====================
BASE_DIR   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs"
RUN_ROOT   <- file.path(BASE_DIR, "cluster_runs_min")

# Choose which bucket/run to read (point directly at a *result* folder set,
# NOT the grid; i.e. where you already have step3_kmeans_*_* folders)
RESULT_ROOT <- RUN_ROOT   # if you want a grid bucket: file.path(RUN_ROOT, "grid_B02500")

DATASETS <- c("A","B","C")
FEATURES <- c("genetic","genetic_minclin","genetic_fullclin")

# How we find the step2a PCA outputs:
PCA_METHOD <- "irlba" # change only if you ran prcomp
# Keep top N SNPs in the TSV / plot
TOP_N <- 50
# Optional stability filter: keep only members of the smaller cluster with stability >= this quantile
STABILITY_QUANTILE <- 0.50   # 0 disables filtering; 0.50 keeps >= median among small cluster
# =======================================================

canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# helper: robust read TSV/CSV
read_tab <- function(path) {
  if (!file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv","txt")) {
    suppressWarnings(tryCatch(fread(path), error = function(e) NULL))
  } else {
    suppressWarnings(tryCatch(fread(path), error = function(e) NULL))
  }
}

# Load clusters for one dataset/feature
load_clusters <- function(ds, feat) {
  # e.g.: step3_kmeans_genetic_A / clusters.tsv
  dir1 <- file.path(RESULT_ROOT, sprintf("step3_kmeans_%s_%s", feat, ds))
  f <- file.path(dir1, "clusters.tsv")
  tab <- read_tab(f)
  if (is.null(tab)) return(NULL)
  # expected: IID, cluster, is_outlier, stability (stability may be missing in PAM runs)
  if (!"IID" %in% names(tab)) {
    # sometimes saved with id or sample
    idcand <- intersect(c("id","sample","SubjectID"), names(tab))
    if (length(idcand)) tab <- tab %>% rename(IID = !!idcand[1])
  }
  tab$IID <- canon_id(tab$IID)
  tab$cluster <- as.integer(tab$cluster)
  tab$feature <- feat
  tab$dataset <- ds
  tab$dir <- dir1
  tab
}

# Load PCA scores and loadings for a dataset
load_pca_bits <- function(ds) {
  d_pca <- file.path(RUN_ROOT, sprintf("step2a_pca_%s_%s", PCA_METHOD, ds))
  scores_path <- file.path(d_pca, "pca_scores.tsv")
  # loadings combined (PC1 & PC2) if available; else read separate and bind
  loads_both <- file.path(d_pca, "pca_top_loadings_PC1_PC2.tsv")
  loads_p1   <- file.path(d_pca, "pca_top_loadings_PC1.tsv")
  loads_p2   <- file.path(d_pca, "pca_top_loadings_PC2.tsv")

  sc <- read_tab(scores_path)
  if (is.null(sc)) return(NULL)
  # normalize
  idcol <- c("IID","id","sample","SubjectID")
  idcol <- idcol[idcol %in% names(sc)][1]
  if (is.na(idcol)) return(NULL)
  sc <- sc %>% rename(IID = !!idcol)
  sc$IID <- canon_id(sc$IID)

  # Need PC1/PC2 present
  if (!all(c("PC1","PC2") %in% names(sc))) {
    return(NULL)
  }

  # Load loadings
  L <- NULL
  if (file.exists(loads_both)) {
    L <- read_tab(loads_both)
    # columns expected: SNP, loading, abs_loading, PC
    if (!all(c("SNP","loading","abs_loading","PC") %in% names(L))) L <- NULL
  }
  if (is.null(L)) {
    L1 <- if (file.exists(loads_p1)) read_tab(loads_p1) else NULL
    L2 <- if (file.exists(loads_p2)) read_tab(loads_p2) else NULL
    if (!is.null(L1)) L1$PC <- "PC1"
    if (!is.null(L2)) L2$PC <- "PC2"
    L <- dplyr::bind_rows(L1, L2)
  }
  if (is.null(L)) return(NULL)

  # Pivot loadings wide: one row per SNP with PC1/PC2 sign
  lw <- L %>%
    select(SNP, PC, loading) %>%
    tidyr::pivot_wider(names_from = PC, values_from = loading)
  list(scores = sc, loads = lw, pca_dir = d_pca)
}

# Compute which PC dominates the separation between clusters
dominant_pc <- function(scores_df, clusters_df) {
  df <- scores_df %>%
    inner_join(clusters_df %>% select(IID, cluster), by = "IID")
  if (nrow(df) < 2 || length(unique(df$cluster)) < 2) {
    return(list(pc = "PC1", diff1 = NA_real_, diff2 = NA_real_))
  }
  mu <- df %>% group_by(cluster) %>%
    summarise(mPC1 = mean(PC1), mPC2 = mean(PC2), .groups = "drop")
  # Ensure two clusters (1 and 2); if more, take the two largest
  big <- mu %>%
    inner_join(df %>% count(cluster, name="n"), by="cluster") %>%
    arrange(desc(n)) %>% head(2)

  if (nrow(big) < 2) return(list(pc = "PC1", diff1 = NA, diff2 = NA))
  # diff = cluster with larger mean - the other
  # We only need absolute magnitude to choose dominance
  d1 <- abs(diff(big$mPC1))
  d2 <- abs(diff(big$mPC2))
  pc <- if (is.na(d1) || is.na(d2)) "PC1" else if (d1 >= d2) "PC1" else "PC2"
  list(pc = pc, diff1 = d1, diff2 = d2,
       means = mu)
}

# Main worker for one dataset/feature
do_one <- function(ds, feat, top_n = TOP_N, q_stab = STABILITY_QUANTILE) {
  clus <- load_clusters(ds, feat)
  if (is.null(clus)) {
    message("No clusters for ", ds, " [", feat, "]")
    return(NULL)
  }
  pca <- load_pca_bits(ds)
  if (is.null(pca)) {
    message("No PCA bits for ", ds, " (scores/loadings missing)")
    return(NULL)
  }
  scores <- pca$scores
  loads  <- pca$loads

  # Identify smaller cluster (by size); if tie, pick the one with lower mean size order
  sizes <- clus %>% count(cluster, name = "n") %>% arrange(n)
  if (nrow(sizes) < 1) return(NULL)
  small_id <- sizes$cluster[1]
  small_n  <- sizes$n[1]

  # stability filter within smaller cluster (if stability column exists and q_stab > 0)
  keep_IIDs <- clus %>% filter(cluster == small_id) %>% pull(IID)
  if ("stability" %in% names(clus) && q_stab > 0 && small_n > 5) {
    st_sub <- clus %>% filter(cluster == small_id & !is.na(stability))
    if (nrow(st_sub) > 3) {
      thr <- quantile(st_sub$stability, probs = q_stab, na.rm = TRUE)
      keep_IIDs <- st_sub %>% filter(stability >= thr) %>% pull(IID)
    }
  }
  clus_small <- clus %>% filter(IID %in% keep_IIDs)

  # Which PC dominates separation?
  dom <- dominant_pc(scores, clus)
  dom_pc <- dom$pc

  # Which cluster has higher mean on dominant PC?
  mu <- scores %>%
    inner_join(clus %>% select(IID, cluster), by = "IID") %>%
    group_by(cluster) %>%
    summarise(mPC1 = mean(PC1), mPC2 = mean(PC2), .groups = "drop")
  # higher mean along dom_pc
  if (dom_pc == "PC1") {
    higher <- mu$cluster[which.max(mu$mPC1)]
    sep_mag <- abs(diff(sort(mu$mPC1)))
  } else {
    higher <- mu$cluster[which.max(mu$mPC2)]
    sep_mag <- abs(diff(sort(mu$mPC2)))
  }

  # Enrichment rule:
  # If small cluster == higher-on-PC, then SNPs with positive loading on that PC
  # are enriched in the small cluster; negative loading enriched in the other cluster.
  # If small cluster != higher, invert sign interpretation.
  loads$loading_dom <- if (dom_pc == "PC1") loads$PC1 else loads$PC2
  loads <- loads %>% filter(!is.na(loading_dom))

  if (length(small_id) == 0 || is.na(higher)) {
    message("Could not determine dominant direction for ", ds, " [", feat, "]")
    return(NULL)
  }

  # direction flag
  same_dir <- (small_id == higher)
  loads <- loads %>%
    mutate(direction = ifelse(loading_dom >= 0,
                              ifelse(same_dir, "enriched_in_small", "enriched_in_large"),
                              ifelse(same_dir, "enriched_in_large", "enriched_in_small")),
           abs_loading = abs(loading_dom)) %>%
    arrange(desc(abs_loading))

  out <- loads %>%
    transmute(
      dataset = ds,
      feature = feat,
      SNP,
      dominant_PC = dom_pc,
      loading = loading_dom,
      abs_loading,
      small_cluster_id = small_id,
      small_cluster_size = small_n,
      small_cluster_stab_q = q_stab,
      separation_magnitude = sep_mag,
      enriched_in_small = direction == "enriched_in_small"
    )

  # Write TSV and plot top N
  out_dir <- file.path(clus$dir[1], "snps_enriched_small")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, sprintf("enriched_snps_%s_%s.tsv", ds, feat))
  fwrite(out, out_file, sep = "\t")

  top_plot <- out %>% slice_max(order_by = abs_loading, n = top_n)
  if (nrow(top_plot) > 0) {
    p <- ggplot(top_plot,
                aes(x = reorder(SNP, abs_loading),
                    y = abs_loading,
                    fill = enriched_in_small)) +
      geom_col() +
      coord_flip() +
      labs(title = sprintf("Top %d SNPs driving %s separation (%s, %s)",
                           top_n, dom_pc, ds, feat),
           x = "SNP", y = sprintf("|loading on %s|", dom_pc),
           fill = "Enriched in small?") +
      theme_minimal(base_size = 12)
    ggsave(filename = file.path(out_dir, sprintf("enriched_snps_top%d_%s_%s.png", top_n, ds, feat)),
           plot = p, width = 9, height = 7, dpi = 300)
  }

  message(sprintf("[OK] %s/%s → small cluster %d (n=%d); dominant=%s; wrote: %s",
                  ds, feat, small_id, small_n, dom_pc, out_file))
  out
}

# ===================== RUN =========================
all_results <- list()
for (ds in DATASETS) {
  for (feat in FEATURES) {
    res <- try(do_one(ds, feat), silent = TRUE)
    if (!inherits(res, "try-error") && !is.null(res)) {
      all_results[[paste(ds, feat, sep = "_")]] <- res
    } else {
      message("Skipped: ", ds, " [", feat, "]")
    }
  }
}

# Optional: combined summary across datasets/features (top 200 by |loading|)
if (length(all_results)) {
  combo <- bind_rows(all_results)
  sum_file <- file.path(RESULT_ROOT, "enriched_snps_summary_across_sets.tsv")
  fwrite(combo %>% arrange(desc(abs_loading)) %>% head(200), sum_file, sep = "\t")
  message("Wrote combined summary: ", sum_file)
} else {
  message("No results produced.")
}



# plot_enriched_snps.R
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(patchwork)
})

# ========= USER PATHS =========
BASE_DIR <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
TSV      <- file.path(BASE_DIR, "enriched_snps_summary_across_sets.tsv")
OUT_DIR  <- file.path(BASE_DIR, "plots_enriched_snps")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ========= READ =========
stopifnot(file.exists(TSV))
dt <- fread(TSV)

# Expect columns: dataset, feature, SNP, dominant_PC, loading, abs_loading,
#                 small_cluster_id, small_cluster_size, small_cluster_stab_q,
#                 separation_magnitude, enriched_in_small (TRUE/FALSE)
req <- c("dataset","feature","SNP","dominant_PC","loading","abs_loading","enriched_in_small")
if (!all(req %in% names(dt))) stop("TSV missing required columns: ", paste(setdiff(req, names(dt)), collapse=", "))

# Clean labels for facets/legends
dt <- dt %>%
  mutate(
    set_label = paste0(dataset, " · ", feature),
    enriched_in_small = as.logical(enriched_in_small)
  )

# ========= 1) Top-N barplots per dataset/feature =========
TOP_N <- 20
plots <- list()
for (lab in unique(dt$set_label)) {
  dsub <- dt %>% filter(set_label == lab) %>% arrange(desc(abs_loading)) %>% slice_head(n = TOP_N)
  if (!nrow(dsub)) next
  p <- ggplot(dsub, aes(x = fct_reorder(SNP, abs_loading), y = abs_loading, fill = enriched_in_small)) +
    geom_col() +
    coord_flip() +
    labs(title = paste0("Top ", TOP_N, " SNPs: ", lab, " (|loading on ", dsub$dominant_PC[1], "|)"),
         x = "SNP", y = paste0("|loading on ", dsub$dominant_PC[1], "|"),
         fill = "Enriched in\nsmall cluster?") +
    theme_minimal(base_size = 12)
  plots[[lab]] <- p
  ggsave(file.path(OUT_DIR, paste0("top", TOP_N, "_", gsub("[^A-Za-z0-9]+","_", lab), ".png")),
         p, width = 9, height = 7, dpi = 300)
}
if (length(plots)) message("✓ Saved Top-", TOP_N, " barplots to: ", OUT_DIR)

# ========= 2) Heatmap across datasets/features =========
# Keep top M SNPs overall (to keep the heatmap readable)
TOP_M_OVERALL <- 50
keep_snps <- dt %>% arrange(desc(abs_loading)) %>% pull(SNP) %>% unique() %>% head(TOP_M_OVERALL)
ht <- dt %>%
  filter(SNP %in% keep_snps) %>%
  # make one column key for set
  mutate(set = set_label) %>%
  select(SNP, set, abs_loading, enriched_in_small)

# For visual cue, map enriched_in_small to shape overlay
# We'll draw a tile for |loading| and a point on top if enriched_in_small=TRUE
# Order axes by overall importance
snp_order <- ht %>% group_by(SNP) %>% summarise(max_abs = max(abs_loading)) %>% arrange(desc(max_abs)) %>% pull(SNP)
set_order <- ht %>% group_by(set) %>% summarise(max_abs = max(abs_loading)) %>% arrange(desc(max_abs)) %>% pull(set)

p_heat <- ggplot(ht, aes(x = factor(set, levels = set_order), y = factor(SNP, levels = snp_order), fill = abs_loading)) +
  geom_tile() +
  geom_point(data = subset(ht, enriched_in_small),
             aes(x = factor(set, levels = set_order), y = factor(SNP, levels = snp_order)),
             shape = 21, stroke = 0.6, size = 1.8, fill = "white", inherit.aes = FALSE) +
  scale_fill_gradient(name = "|loading|", low = "white", high = "steelblue") +
  labs(title = paste0("SNP influence heatmap (top ", TOP_M_OVERALL, " SNPs)"),
       x = "Dataset · Feature", y = "SNP") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, paste0("heatmap_top", TOP_M_OVERALL, "_snps.png")),
       p_heat, width = 12, height = 10, dpi = 300)
message("✓ Saved heatmap to: ", OUT_DIR)

# ========= 3) Lollipop plot (single set) =========
# Pick the set with the strongest separation (proxy = max |loading| median)
set_stats <- dt %>%
  group_by(set_label, dominant_PC) %>%
  summarise(med_abs = median(abs_loading, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med_abs))

if (nrow(set_stats)) {
  pick <- set_stats$set_label[1]
  dsub <- dt %>% filter(set_label == pick) %>% arrange(desc(abs_loading)) %>% slice_head(n = 30)
  p_lolli <- ggplot(dsub, aes(x = fct_reorder(SNP, abs_loading), y = abs_loading, color = enriched_in_small)) +
    geom_segment(aes(xend = SNP, y = 0, yend = abs_loading), linewidth = 0.7) +
    geom_point(size = 2) +
    coord_flip() +
    labs(title = paste0("Lollipop: ", pick, " (top 30 by |loading|)"),
         x = "SNP", y = paste0("|loading on ", dsub$dominant_PC[1], "|"),
         color = "Enriched in\nsmall cluster?") +
    theme_minimal(base_size = 12)
  ggsave(file.path(OUT_DIR, paste0("lollipop_top30_", gsub("[^A-Za-z0-9]+","_", pick), ".png")),
         p_lolli, width = 9, height = 7, dpi = 300)
  message("✓ Saved lollipop plot for: ", pick)
}

# ========= 4) Quick summary table (per set) =========
sum_tbl <- dt %>%
  group_by(set_label) %>%
  summarise(
    n_snps = n(),
    top5 = paste(head(SNP[order(-abs_loading)] ,5), collapse = ", "),
    frac_enriched_small = mean(enriched_in_small, na.rm = TRUE),
    domPC = names(which.max(table(dominant_PC))),
    .groups = "drop"
  )
fwrite(sum_tbl, file.path(OUT_DIR, "enriched_snps_quick_summary.tsv"), sep = "\t")
message("✓ Wrote quick summary TSV to: ", OUT_DIR)

# 1. Install BiocManager if you don’t have it yet
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 2. Install the SNP database package (dbSNP build 155 for GRCh38)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")

# 3. Load the package
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)

###TOP SNP
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(GenomicRanges)

# Define your query position
gr <- GRanges(seqnames = "6", ranges = IRanges(start = 32616083, end = 32616083))

# find overlapping SNPs
hits <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)
hits
