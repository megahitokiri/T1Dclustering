## =========================
## ARI: compare feature sets
## Bucket: B02500
## =========================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)      # for adjustedRandIndex
})

# ---- Paths ----
ROOT <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
BUCKET <- file.path(ROOT, "grid_B02500")   # <- change if you compare a different grid
OUTDIR <- file.path(BUCKET, "viz_summary")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Feature sets to compare (must match your folder names under BUCKET)
FEATURES <- c("genetic", "genetic_minclin", "genetic_fullclin", "clinical_only")
DATASETS <- c("A","B","C")

# ---- Helper: read clusters.tsv robustly ----
read_clusters <- function(bucket, feature, ds) {
  # expected folder: step3_kmeans_<feature>_<dataset>/clusters.tsv
  pat <- sprintf("^step3_kmeans_%s_%s$", feature, ds)
  dirs <- list.dirs(bucket, recursive = FALSE, full.names = TRUE)
  hit  <- dirs[grepl(pat, basename(dirs))]
  if (!length(hit)) return(NULL)
  f <- file.path(hit[1], "clusters.tsv")
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, na.strings = c("", "NA"))
  # standardize ID + cluster col names
  idc <- c("IID","SubjectID","id","sample")
  idc <- idc[idc %in% names(dt)][1]
  if (is.na(idc)) return(NULL)
  setnames(dt, idc, "IID")
  if (!"cluster" %in% names(dt)) {
    cand <- c("kmeans_cluster","Cluster","clusters","label")
    hitc <- cand[cand %in% names(dt)][1]
    if (!is.na(hitc)) setnames(dt, hitc, "cluster")
  }
  if (!"cluster" %in% names(dt)) return(NULL)
  dt[, .(IID = as.character(IID), cluster = as.integer(cluster))]
}

# ---- Compute ARI per dataset with PAIR-WISE overlap ----
all_rows <- list()

ari_heatmap_one_dataset <- function(ds, feats, bucket, outdir) {
  # read clusters for every available feature
  clu <- lapply(feats, function(ft) read_clusters(bucket, ft, ds))
  names(clu) <- feats
  keep <- !vapply(clu, is.null, logical(1))
  if (!any(keep)) return(invisible(NULL))
  clu <- clu[keep]
  feats_avail <- names(clu)
  
  # build ARI and N matrices using pair-wise intersections
  F <- length(feats_avail)
  M_ari <- matrix(NA_real_, F, F, dimnames = list(feats_avail, feats_avail))
  M_n   <- matrix(0L,        F, F, dimnames = list(feats_avail, feats_avail))
  
  for (i in seq_len(F)) {
    for (j in seq_len(F)) {
      ids <- intersect(clu[[i]]$IID, clu[[j]]$IID)
      if (length(ids) >= 5) {
        xi <- clu[[i]][IID %in% ids][order(IID)]$cluster
        xj <- clu[[j]][IID %in% ids][order(IID)]$cluster
        M_ari[i, j] <- mclust::adjustedRandIndex(xi, xj)
        M_n[i, j]   <- length(ids)
        if (i < j) {
          all_rows[[length(all_rows) + 1L]] <<- data.table(
            dataset = ds,
            feature1 = feats_avail[i],
            feature2 = feats_avail[j],
            n_common = length(ids),
            ARI = M_ari[i, j]
          )
        }
      } else {
        M_ari[i, j] <- NA_real_
        M_n[i, j]   <- length(ids)
      }
    }
  }
  
  # tidy for ggplot with ARI label + n
  dt <- as.data.table(as.table(M_ari))
  setnames(dt, c("feature1","feature2","ARI"))
  dt[, n_common := as.vector(M_n)]
  dt[, lab := ifelse(is.na(ARI), "–",
                     sprintf("%.2f\nn=%d", ARI, n_common))]
  
  p <- ggplot(dt, aes(feature1, feature2, fill = ARI)) +
    geom_tile(color = "white") +
    geom_text(aes(label = lab), size = 3.5, lineheight = 0.9) +
    scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#d7301f",
                         midpoint = 0, limits = c(-0.1, 1), na.value = "grey90") +
    coord_equal() +
    labs(title = sprintf("Adjusted Rand Index — dataset %s", ds),
         x = NULL, y = NULL, fill = "ARI") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave(file.path(outdir, sprintf("B02500_ARI_heatmap_%s.png", ds)),
         p, width = 7, height = 5.5, dpi = 300)
}

# run for A/B/C
for (ds in DATASETS) {
  ari_heatmap_one_dataset(ds, FEATURES, BUCKET, OUTDIR)
}

# pairwise table
pairwise_ari <- if (length(all_rows)) rbindlist(all_rows, use.names = TRUE) else data.table()
fwrite(pairwise_ari, file.path(OUTDIR, "B02500_ARI_pairwise.tsv"), sep = "\t")

cat(sprintf(
  "\n✔ ARI done. Outputs in %s\n  - B02500_ARI_heatmap_A.png (and B/C)\n  - B02500_ARI_pairwise.tsv\n",
  OUTDIR
))
