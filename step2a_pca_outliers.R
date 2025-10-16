# step2a_pca_rstudio_all_in_one.R
# Run directly in RStudio (no CLI). Outputs under:
# out_three_runs/cluster_runs_min/step2a_pca_<method>_<A|B|C>/

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(data.table); library(ggplot2); library(stringr)
})

# =================== USER CONFIG ============================================
BASE_DIR   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper"
METHOD     <- "irlba"   # "irlba" (recommended) or "prcomp"
VAR_TARGET <- 0.80
MAX_PCS    <- 15
TOP_N_LOAD <- 10
DATASETS   <- c("A","B","C")  # e.g., c("A","C") to test fewer

# =================== Paths / setup ==========================================
OUT_DIR  <- file.path(BASE_DIR, "out_three_runs")
RUN_DIR  <- file.path(OUT_DIR, "cluster_runs_min")
DAT_DIR  <- file.path(OUT_DIR, "dat")
dir.create(RUN_DIR, showWarnings = FALSE, recursive = TRUE)

RDS_FILES <- list(
  A = file.path(OUT_DIR, "A_homog_numeric.rds"),
  B = file.path(OUT_DIR, "B_homog_numeric.rds"),
  C = file.path(OUT_DIR, "C_homog_numeric.rds")
)

SEED <- 12345
set.seed(SEED)

cat("Working directory:", getwd(), "\n")
cat("BASE_DIR exists:", dir.exists(BASE_DIR), "\n")
cat("OUT_DIR exists:", dir.exists(OUT_DIR), "\n")
cat("RUN_DIR:", RUN_DIR, "\n")

# CEPH samples to drop from A/B/C
SAMPLES_DROP <- c("NA12146","NA10830","NA12236")

# =================== Helpers (robust) =======================================
extract_IIDs <- function(obj){
  if (!is.null(obj$G) && !is.null(rownames(obj$G))) return(as.character(rownames(obj$G)))
  if (!is.null(obj$fam)) {
    if ("iid" %in% names(obj$fam)) return(as.character(obj$fam$iid))
    if (ncol(obj$fam) >= 2)       return(as.character(obj$fam[[2]]))
  }
  if (!is.null(obj$G)) return(as.character(seq_len(nrow(obj$G))))
  character()
}

prep_geno <- function(G){
  G <- as.matrix(G); storage.mode(G) <- "double"
  cm <- suppressWarnings(colMeans(G, na.rm = TRUE)); cm[!is.finite(cm)] <- 0
  for (j in seq_len(ncol(G))) if (anyNA(G[, j])) G[is.na(G[, j]), j] <- cm[j]
  col_mu <- colMeans(G); G <- sweep(G, 2, col_mu, "-")
  sds <- sqrt(colSums(G^2) / pmax(1, nrow(G) - 1)); sds[!is.finite(sds) | sds == 0] <- 1
  G <- sweep(G, 2, sds, "/"); G[!is.finite(G)] <- 0
  list(G = G, miss_snp = rep(0, ncol(G)), miss_samp = rep(0, nrow(G)))
}

safe_prcomp <- function(G){
  x <- tryCatch(prcomp(G, center = FALSE, scale. = FALSE, retx = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)
  x <- tryCatch(prcomp(scale(G, center = TRUE, scale = FALSE), center = FALSE, scale. = FALSE, retx = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)
  x <- tryCatch(prcomp(scale(G), center = FALSE, scale. = FALSE, retx = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)
  X <- scale(G, center = TRUE, scale = FALSE); X[!is.finite(X)] <- 0
  sv <- tryCatch(svd(X), error = function(e) NULL); if (is.null(sv)) return(NULL)
  out <- list()
  out$x        <- sv$u %*% diag(sv$d, nrow = length(sv$d))
  colnames(out$x) <- paste0("PC", seq_len(ncol(out$x)))
  out$sdev     <- sv$d / sqrt(max(1, nrow(G) - 1))
  out$rotation <- sv$v; rownames(out$rotation) <- colnames(G); colnames(out$rotation) <- paste0("PC", seq_len(ncol(out$rotation)))
  class(out) <- "prcomp"; out
}

choose_k_from_ve <- function(eigvals, target=VAR_TARGET, max_pcs=MAX_PCS){
  total <- sum(eigvals)
  if (!is.finite(total) || total == 0) { ve <- rep(0, length(eigvals)); cumve <- cumsum(ve); return(list(n=2, ve=ve, cumve=cumve)) }
  ve <- eigvals / total; cumve <- cumsum(ve)
  n <- which(cumve >= target)[1]; if (is.na(n)) n <- length(eigvals)
  n <- min(max(n, 2), max_pcs); list(n=n, ve=ve, cumve=cumve)
}

top_loadings_dt <- function(load_mat, pc_index, row_names, top_n = TOP_N_LOAD){
  v <- load_mat[, pc_index]; ord <- order(abs(v), decreasing = TRUE); idx <- head(ord, min(top_n, length(v)))
  data.table(SNP = row_names[idx], loading = v[idx], abs_loading = abs(v[idx]), PC = paste0("PC", pc_index))
}

# ---- Outlier helpers ----
mahal_outliers <- function(scores_df, k = 5, p_cut = 0.001) {
  pcs <- grep("^PC\\d+$", names(scores_df), value = TRUE); if (!length(pcs)) return(data.frame())
  keep <- pcs[seq_len(min(k, length(pcs)))]; X <- as.matrix(scores_df[, keep, drop = FALSE]); rownames(X) <- scores_df$IID
  mu <- colMeans(X); S <- stats::cov(X); if (!all(is.finite(S))) return(data.frame())
  md <- stats::mahalanobis(X, mu, S); cutoff <- stats::qchisq(1 - p_cut, df = ncol(X))
  data.frame(IID = rownames(X), md = md, md_flag = md > cutoff, row.names = NULL)
}
zscore_outliers <- function(scores_df, pc = "PC1", z = 4) {
  if (!pc %in% names(scores_df)) return(data.frame()); x <- as.numeric(scores_df[[pc]])
  zval <- (x - mean(x)) / sd(x); data.frame(IID = scores_df$IID, z = zval, z_flag = abs(zval) > z)
}
add_labels_if <- function(plot_obj, df, aes_id = aes(label = IID)) {
  if (requireNamespace("ggrepel", quietly = TRUE)) plot_obj + ggrepel::geom_text_repel(data = df, aes_id, size = 3, max.overlaps = 20) else plot_obj
}

# Ensure irlba if requested
if (tolower(METHOD) == "irlba" && !requireNamespace("irlba", quietly = TRUE)) {
  stop("METHOD='irlba' but package 'irlba' is not installed. Run install.packages('irlba') or set METHOD <- 'prcomp'.")
} else if (tolower(METHOD) == "irlba") {
  suppressPackageStartupMessages(library(irlba))
}

# =================== Runner (per dataset) ===================================
run_one_dataset <- function(ds, method = tolower(METHOD)){
  message("\n=== PCA-only (", method, ") :: Dataset ", ds, " ===")
  
  rds_path <- RDS_FILES[[ds]]
  message("RDS path: ", rds_path, "  (exists=", file.exists(rds_path), ")")
  if (!file.exists(rds_path)) { warning("Missing RDS: ", rds_path); return(NULL) }
  
  obj <- readRDS(rds_path)
  if (!is.list(obj) || is.null(obj$G)) stop("Unexpected RDS structure: ", rds_path)
  
  IIDs <- extract_IIDs(obj); Graw <- obj$G; rownames(Graw) <- IIDs
  
  # ----- Remove specific CEPH samples (by IID) before any processing -----
  if (length(SAMPLES_DROP)) {
    keep_idx <- !(rownames(Graw) %in% SAMPLES_DROP)
    n_removed <- sum(!keep_idx)
    if (n_removed > 0) {
      message("Dataset ", ds, ": removing ", n_removed,
              " specific CEPH sample(s): ",
              paste(rownames(Graw)[!keep_idx], collapse = ", "))
      # subset genotype matrix
      Graw <- Graw[keep_idx, , drop = FALSE]
      # keep fam in sync if present
      if (!is.null(obj$fam)) {
        fam_df <- as.data.frame(obj$fam)
        # try to find IID column
        if ("iid" %in% names(fam_df)) {
          obj$fam <- fam_df[fam_df$iid %in% rownames(Graw), , drop = FALSE]
        } else if (ncol(fam_df) >= 2) {
          obj$fam <- fam_df[fam_df[[2]] %in% rownames(Graw), , drop = FALSE]
        }
      }
    } else {
      message("Dataset ", ds, ": none of the specified CEPH IIDs present.")
    }
  }
  
  ds_dir <- file.path(RUN_DIR, paste0("step2a_pca_", method, "_", ds)); dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)
  message("Output dir: ", ds_dir)
  
  prep <- prep_geno(Graw); G <- prep$G
  fwrite(data.table(SNP = colnames(G), missing_rate = prep$miss_snp), file = file.path(ds_dir, "missingness_snp.tsv"), sep = "\t")
  fwrite(data.table(IID = rownames(G),  missing_rate = prep$miss_samp), file = file.path(ds_dir, "missingness_sample.tsv"), sep = "\t")
  
  if (method == "irlba") {
    k <- min(MAX_PCS, nrow(G) - 1L, ncol(G))
    message("Running irlba with k=", k, " (n=", nrow(G), ", p=", ncol(G), ")")
    ir <- irlba::irlba(G, nv = k, nu = k, center = FALSE, scale = FALSE)
    scores <- ir$u %*% diag(ir$d, nrow = length(ir$d), ncol = length(ir$d))
    rownames(scores) <- rownames(G); colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
    eigvals <- (ir$d ^ 2); load_mat <- ir$v; row_names <- colnames(G)
  } else {
    message("Running prcomp() …")
    pca <- safe_prcomp(G); if (is.null(pca)) stop("prcomp failed even after fallbacks. Set METHOD <- 'irlba'.")
    scores <- pca$x; eigvals <- (pca$sdev ^ 2); load_mat <- pca$rotation; row_names <- colnames(G)
  }
  
  pcs_meta <- choose_k_from_ve(eigvals, target = VAR_TARGET, max_pcs = MAX_PCS); ve <- pcs_meta$ve; cumve <- pcs_meta$cumve
  
#  scores_df <- as.data.frame(scores); scores_df$IID <- rownames(scores)
#  fwrite(as.data.table(scores_df), file = file.path(ds_dir, "pca_scores.tsv"), sep = "\t")

  # ---- Add sample IDs (IIDs) to PCA scores ----
  scores_df <- as.data.frame(scores)
  scores_df$IID <- rownames(scores)                     # add IID as column
  scores_df <- dplyr::relocate(scores_df, IID)          # move IID to first column
  
  # sanity check
  stopifnot(!anyDuplicated(scores_df$IID))
  stopifnot(any(grepl("^PC\\d+$", names(scores_df))))
  
  # write PCA scores with IDs
  data.table::fwrite(
    data.table::as.data.table(scores_df),
    file = file.path(ds_dir, "pca_scores.tsv"),
    sep = "\t"
  )    
  
  fwrite(data.table(PC = seq_along(ve), variance = ve, cum_variance = cumve), file = file.path(ds_dir, "pca_variance.tsv"), sep = "\t")
  
  if (ncol(scores) >= 2) {
    p <- ggplot(scores_df, aes(PC1, PC2)) +
      geom_point(alpha = 0.8, size = 2) +
      labs(title = paste("PCA -", toupper(method), "- dataset", ds),
           x = sprintf("PC1 (%.1f%%)", ve[1]*100), y = sprintf("PC2 (%.1f%%)", ve[2]*100)) +
      theme_minimal()
    ggsave(file.path(ds_dir, "pca_plot.png"), p, width = 6, height = 5, dpi = 300)
  }
  
  scree_dt <- data.table(PC = seq_along(ve), Variance = ve, Cumulative = cumve)
  p_scree <- ggplot(scree_dt, aes(PC, Variance)) + geom_col() +
    geom_line(aes(y = Cumulative), group = 1) + geom_point(aes(y = Cumulative)) +
    theme_minimal() + labs(title = paste("Scree plot -", toupper(method), "-", ds),
                           y = "Variance (bar) / Cumulative (line)")
  ggsave(file.path(ds_dir, "pca_scree.png"), p_scree, width = 7, height = 5, dpi = 300)
  
  if (!is.null(load_mat) && ncol(load_mat) >= 1) {
    if (is.null(rownames(load_mat))) rownames(load_mat) <- row_names
    top1 <- top_loadings_dt(load_mat, 1, row_names, TOP_N_LOAD); fwrite(top1, file = file.path(ds_dir, "pca_top_loadings_PC1.tsv"), sep = "\t")
    if (ncol(load_mat) >= 2) {
      top2 <- top_loadings_dt(load_mat, 2, row_names, TOP_N_LOAD)
      fwrite(top2, file = file.path(ds_dir, "pca_top_loadings_PC2.tsv"), sep = "\t")
      fwrite(rbind(top1, top2), file = file.path(ds_dir, "pca_top_loadings_PC1_PC2.tsv"), sep = "\t")
    }
  }
  
  # Outliers (Mahalanobis + PC1 z)
  mout <- mahal_outliers(scores_df, k = 5, p_cut = 0.001)
  zout <- zscore_outliers(scores_df, pc = "PC1", z = 4)
  out_df <- merge(mout, zout, by = "IID", all = TRUE)
  out_df$md_flag[is.na(out_df$md_flag)] <- FALSE
  out_df$z_flag[is.na(out_df$z_flag)]   <- FALSE
  out_df$is_outlier <- out_df$md_flag | out_df$z_flag
  keep_cols <- intersect(c("IID","PC1","PC2"), names(scores_df))
  out_export <- merge(out_df, scores_df[, keep_cols, drop = FALSE], by = "IID", all.x = TRUE)
  fwrite(out_export, file.path(ds_dir, "pca_outliers.tsv"), sep = "\t")
  
  if (ncol(scores) >= 2) {
    plot_df <- merge(scores_df, out_df[, c("IID","is_outlier"), drop = FALSE], by = "IID", all.x = TRUE)
    plot_df$is_outlier[is.na(plot_df$is_outlier)] <- FALSE
    gp <- ggplot(plot_df, aes(PC1, PC2)) +
      geom_point(aes(color = is_outlier), alpha = 0.9, size = 2) +
      scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
      theme_minimal() +
      labs(title = paste("PCA -", toupper(method), "-", ds, "(outliers in red)"),
           x = sprintf("PC1 (%.1f%%)", ve[1]*100), y = sprintf("PC2 (%.1f%%)", ve[2]*100), color = "Outlier")
    lab_df <- subset(plot_df, is_outlier); if (nrow(lab_df)) gp <- add_labels_if(gp, lab_df, aes(label = IID))
    ggsave(file.path(ds_dir, "pca_plot_outliers.png"), gp, width = 7, height = 5, dpi = 300)
  }
  
  # OPTIONAL: PC ↔ clinical quick check (Age/Sex)
  clin_path <- file.path(DAT_DIR, "cluster_clinical.csv")
  message("Clinical CSV path:", clin_path, "  (exists=", file.exists(clin_path), ")")
  if (file.exists(clin_path)) {
    clin <- suppressMessages(readr::read_csv(clin_path, show_col_types = FALSE))
    if ("Cluster_subjectID" %in% names(clin)) {
      mini <- merge(scores_df[, c("IID","PC1","PC2")], clin[, c("Cluster_subjectID","Cluster_age","Cluster_sex")],
                    by.x = "IID", by.y = "Cluster_subjectID", all.x = FALSE)
      if (nrow(mini)) {
        if (!is.numeric(mini$Cluster_sex)) {
          sx <- tolower(as.character(mini$Cluster_sex))
          mini$Cluster_sex <- ifelse(grepl("^m", sx), 1, ifelse(grepl("^f", sx), 0, NA_real_))
        }
        cor_tbl <- data.frame(
          var = c("Age","Sex"),
          cor_PC1 = c(suppressWarnings(cor(mini$PC1, mini$Cluster_age, use="complete.obs")),
                      suppressWarnings(cor(mini$PC1, mini$Cluster_sex, use="complete.obs"))),
          cor_PC2 = c(suppressWarnings(cor(mini$PC2, mini$Cluster_age, use="complete.obs")),
                      suppressWarnings(cor(mini$PC2, mini$Cluster_sex, use="complete.obs")))
        )
        fwrite(cor_tbl, file.path(ds_dir, "pca_clinical_cor_PC1_PC2.tsv"), sep = "\t")
      }
    }
  }
  
  # OPTIONAL: batch-colored PCA if fam has 'batch'
  if (!is.null(obj$fam) && "batch" %in% names(obj$fam) && ncol(scores) >= 2) {
    fam_df <- as.data.frame(obj$fam)
    fam_df$IID <- if ("iid" %in% names(fam_df)) as.character(fam_df$iid) else as.character(fam_df[[2]])
    plot_b <- merge(scores_df[, c("IID","PC1","PC2")],
                    fam_df[, c("IID","batch")],
                    by = "IID", all.x = TRUE)
    pb <- ggplot(plot_b, aes(PC1, PC2, color = batch)) +
      geom_point(alpha = 0.9, size = 2) +
      theme_minimal() +
      labs(title = paste("PCA -", toupper(method), "-", ds, "(colored by batch)"),
           x = sprintf("PC1 (%.1f%%)", ve[1]*100),
           y = sprintf("PC2 (%.1f%%)", ve[2]*100))
    ggsave(file.path(ds_dir, "pca_plot_batch.png"), pb, width = 7, height = 5, dpi = 300)
  }
  
  # return summary row
  data.frame(
    dataset    = ds,
    n_samples  = nrow(G),
    n_snps     = ncol(G),
    var_pc1    = ve[1],
    var_pc2    = ve[2],
    cumvar_pc2 = cumve[2],
    method     = method,
    stringsAsFactors = FALSE
  )
}

# =================== EXECUTE RUN ============================================
cat("\nStarting PCA for datasets:", paste(DATASETS, collapse = ", "), "\n")
results <- lapply(DATASETS, run_one_dataset, method = tolower(METHOD))
results <- results[!vapply(results, is.null, TRUE)]

if (length(results)) {
  sumdt   <- dplyr::bind_rows(results)
  outfile <- file.path(RUN_DIR, paste0("step2a_pca_", tolower(METHOD), "_summary.tsv"))
  data.table::fwrite(sumdt, outfile, sep = "\t")
  message("\nSummary written to: ", outfile)
} else {
  warning("No datasets processed.")
}
message("\n✅ Done. Check outputs under: ", RUN_DIR)

