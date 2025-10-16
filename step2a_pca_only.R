# step2a_pca_rstudio_selfcontained.R
# Run directly in RStudio (no CLI, no source() gymnastics).
# Produces per-dataset outputs under out_three_runs/cluster_runs_min/.

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(data.table); library(ggplot2); library(stringr)
})

# =================== USER CONFIG (edit these) ===============================
BASE_DIR   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper"
METHOD     <- "irlba"   # "irlba" (recommended) or "prcomp"
VAR_TARGET <- 0.80
MAX_PCS    <- 10
TOP_N_LOAD <- 10
DATASETS   <- c("A","B","C")  # choose a subset like c("A","C") to test

# =================== Paths / setup ==========================================
OUT_DIR  <- file.path(BASE_DIR, "out_three_runs")
RUN_DIR  <- file.path(OUT_DIR, "cluster_runs_min")
dir.create(RUN_DIR, showWarnings = FALSE, recursive = TRUE)

RDS_FILES <- list(
  A = file.path(OUT_DIR, "A_homog_numeric.rds"),
  B = file.path(OUT_DIR, "B_homog_numeric.rds"),
  C = file.path(OUT_DIR, "C_homog_numeric.rds")
)

SEED <- 12345
set.seed(SEED)

# =================== Helpers (robust) =======================================
extract_IIDs <- function(obj){
  if (!is.null(obj$G) && !is.null(rownames(obj$G))) return(as.character(rownames(obj$G)))
  if (!is.null(obj$fam)) {
    if ("iid" %in% names(obj$fam)) return(as.character(obj$fam$iid))
    if (ncol(obj$fam) >= 2) return(as.character(obj$fam[[2]]))
  }
  if (!is.null(obj$G)) return(as.character(seq_len(nrow(obj$G))))
  character()
}

prep_geno <- function(G){
  G <- as.matrix(G); storage.mode(G) <- "double"
  
  # Impute by col-mean (0 if column is all NA)
  cm <- suppressWarnings(colMeans(G, na.rm = TRUE))
  cm[!is.finite(cm)] <- 0
  for (j in seq_len(ncol(G))) {
    if (anyNA(G[, j])) {
      mu <- if (is.finite(cm[j])) cm[j] else 0
      G[is.na(G[, j]), j] <- mu
    }
  }
  
  # Center/scale with guards
  col_mu <- colMeans(G)
  G <- sweep(G, 2, col_mu, "-")
  sds <- sqrt(colSums(G^2) / pmax(1, nrow(G) - 1))
  sds[!is.finite(sds) | sds == 0] <- 1
  G <- sweep(G, 2, sds, "/")
  G[!is.finite(G)] <- 0
  
  list(
    G = G,
    miss_snp  = rep(0, ncol(G)),  # imputed -> 0
    miss_samp = rep(0, nrow(G))
  )
}

safe_prcomp <- function(G){
  # Try prcomp on already standardized matrix
  x <- tryCatch(prcomp(G, center = FALSE, scale. = FALSE, retx = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)
  # Try centering-only and full scaling variants
  x <- tryCatch(prcomp(scale(G, center = TRUE, scale = FALSE), center = FALSE, scale. = FALSE, retx = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)
  x <- tryCatch(prcomp(scale(G), center = FALSE, scale. = FALSE, retx = TRUE), error = function(e) NULL)
  if (!is.null(x)) return(x)
  # Final fallback: SVD -> prcomp-like object
  X <- scale(G, center = TRUE, scale = FALSE); X[!is.finite(X)] <- 0
  sv <- tryCatch(svd(X), error = function(e) NULL)
  if (is.null(sv)) return(NULL)
  out <- list()
  out$x        <- sv$u %*% diag(sv$d, nrow = length(sv$d))
  colnames(out$x) <- paste0("PC", seq_len(ncol(out$x)))
  out$sdev     <- sv$d / sqrt(max(1, nrow(G) - 1))
  out$rotation <- sv$v
  rownames(out$rotation) <- colnames(G)
  colnames(out$rotation) <- paste0("PC", seq_len(ncol(out$rotation)))
  class(out) <- "prcomp"
  out
}

choose_k_from_ve <- function(eigvals, target=VAR_TARGET, max_pcs=MAX_PCS){
  total <- sum(eigvals)
  if (!is.finite(total) || total == 0) {
    ve <- rep(0, length(eigvals)); cumve <- cumsum(ve)
    return(list(n=2, ve=ve, cumve=cumve))
  }
  ve <- eigvals / total
  cumve <- cumsum(ve)
  n <- which(cumve >= target)[1]; if (is.na(n)) n <- length(eigvals)
  n <- min(max(n, 2), max_pcs)
  list(n=n, ve=ve, cumve=cumve)
}

top_loadings_dt <- function(load_mat, pc_index, row_names, top_n = TOP_N_LOAD){
  v <- load_mat[, pc_index]
  ord <- order(abs(v), decreasing = TRUE)
  idx <- head(ord, min(top_n, length(v)))
  data.table(SNP = row_names[idx], loading = v[idx], abs_loading = abs(v[idx]), PC = paste0("PC", pc_index))
}

# If METHOD == "irlba", ensure package is present
if (tolower(METHOD) == "irlba" && !requireNamespace("irlba", quietly = TRUE)) {
  stop("METHOD='irlba' but package 'irlba' is not installed. Run install.packages('irlba') or set METHOD <- 'prcomp'.")
}

# =================== Runner (per dataset) ===================================
run_one_dataset <- function(ds, method = tolower(METHOD)){
  message("\n=== PCA-only (", method, ") :: Dataset ", ds, " ===")
  
  rds_path <- RDS_FILES[[ds]]
  if (!file.exists(rds_path)) { warning("Missing RDS: ", rds_path); return(NULL) }
  obj <- readRDS(rds_path)
  if (!is.list(obj) || is.null(obj$G)) stop("Unexpected RDS structure: ", rds_path)
  
  IIDs <- extract_IIDs(obj)
  Graw <- obj$G; rownames(Graw) <- IIDs
  
  # prepare output dir
  ds_dir <- file.path(RUN_DIR, paste0("step2a_pca_", method, "_", ds))
  dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Prep matrix
  prep <- prep_geno(Graw); G <- prep$G
  
  # (missingness now zeros due to impute; still write files for completeness)
  fwrite(data.table(SNP = colnames(G), missing_rate = prep$miss_snp),
         file = file.path(ds_dir, "missingness_snp.tsv"), sep = "\t")
  fwrite(data.table(IID = rownames(G), missing_rate = prep$miss_samp),
         file = file.path(ds_dir, "missingness_sample.tsv"), sep = "\t")
  
  # PCA
  if (method == "irlba") {
    k <- min(MAX_PCS, nrow(G) - 1L, ncol(G))
    ir <- irlba::irlba(G, nv = k, nu = k, center = FALSE, scale = FALSE)
    scores <- ir$u %*% diag(ir$d, nrow = length(ir$d), ncol = length(ir$d))
    rownames(scores) <- rownames(G)
    colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
    eigvals <- (ir$d ^ 2)
    load_mat <- ir$v; row_names <- colnames(G)
  } else {
    pca <- safe_prcomp(G)
    if (is.null(pca)) stop("prcomp failed even after fallbacks. Set METHOD <- 'irlba'.")
    scores <- pca$x
    eigvals <- (pca$sdev ^ 2)
    load_mat <- pca$rotation; row_names <- colnames(G)
  }
  
  pcs_meta <- choose_k_from_ve(eigvals, target = VAR_TARGET, max_pcs = MAX_PCS)
  ve <- pcs_meta$ve; cumve <- pcs_meta$cumve
  
  # Save outputs
  sc_dt <- as.data.table(scores); sc_dt[, IID := rownames(scores)]
  fwrite(sc_dt, file = file.path(ds_dir, "pca_scores.tsv"), sep = "\t")
  fwrite(data.table(PC = seq_along(ve), variance = ve, cum_variance = cumve),
         file = file.path(ds_dir, "pca_variance.tsv"), sep = "\t")
  
  if (ncol(scores) >= 2) {
    p <- ggplot(as.data.frame(scores), aes(x = PC1, y = PC2)) +
      geom_point(alpha = 0.8, size = 2) +
      labs(title = paste("PCA -", toupper(method), "- dataset", ds),
           x = sprintf("PC1 (%.1f%%)", ve[1]*100),
           y = sprintf("PC2 (%.1f%%)", ve[2]*100)) +
      theme_minimal()
    ggsave(file.path(ds_dir, "pca_plot.png"), p, width = 6, height = 5, dpi = 300)
  }
  
  scree_dt <- data.table(PC = seq_along(ve), Variance = ve, Cumulative = cumve)
  p_scree <- ggplot(scree_dt, aes(x = PC, y = Variance)) +
    geom_col() +
    geom_line(aes(y = Cumulative), group = 1) +
    geom_point(aes(y = Cumulative)) +
    theme_minimal() +
    labs(title = paste("Scree plot -", toupper(method), "-", ds), y = "Variance (bar) / Cumulative (line)")
  ggsave(file.path(ds_dir, "pca_scree.png"), p_scree, width = 7, height = 5, dpi = 300)
  
  if (!is.null(load_mat) && ncol(load_mat) >= 1) {
    if (is.null(rownames(load_mat))) rownames(load_mat) <- row_names
    top1 <- top_loadings_dt(load_mat, 1, row_names, TOP_N_LOAD)
    fwrite(top1, file = file.path(ds_dir, "pca_top_loadings_PC1.tsv"), sep = "\t")
    if (ncol(load_mat) >= 2) {
      top2 <- top_loadings_dt(load_mat, 2, row_names, TOP_N_LOAD)
      fwrite(top2, file = file.path(ds_dir, "pca_top_loadings_PC2.tsv"), sep = "\t")
      fwrite(rbind(top1, top2), file = file.path(ds_dir, "pca_top_loadings_PC1_PC2.tsv"), sep = "\t")
    }
  }
  
  data.frame(
    dataset = ds,
    n_samples = nrow(G), n_snps = ncol(G),
    var_pc1 = ve[1], var_pc2 = ve[2],
    cumvar_pc2 = cumve[2],
    method = method,
    stringsAsFactors = FALSE
  )
}

# =================== Run all selected datasets ==============================
summary_rows <- lapply(DATASETS, run_one_dataset, method = tolower(METHOD))
summary_rows <- summary_rows[!vapply(summary_rows, is.null, TRUE)]
if (length(summary_rows)) {
  sumdt <- bind_rows(summary_rows)
  outfile <- file.path(RUN_DIR, paste0("step2a_pca_", tolower(METHOD), "_summary.tsv"))
  fwrite(sumdt, outfile, sep = "\t")
  message("\nSummary written to: ", outfile)
} else {
  warning("No datasets processed.")
}

message("\n✅ Done. Check outputs under: ", RUN_DIR)

