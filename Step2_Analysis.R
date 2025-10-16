#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(stringr)
})

## ======================= CONFIG (fast + minimal) ============================
BASE_DIR <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper"
OUT_DIR  <- file.path(BASE_DIR, "out_three_runs")
DAT_DIR  <- file.path(OUT_DIR, "dat")
RUN_DIR  <- file.path(OUT_DIR, "cluster_runs_min")
dir.create(RUN_DIR, showWarnings = FALSE, recursive = TRUE)

RDS_FILES <- list(
  A = file.path(OUT_DIR, "A_homog_numeric.rds"),
  B = file.path(OUT_DIR, "B_homog_numeric.rds"),
  C = file.path(OUT_DIR, "C_homog_numeric.rds")
)

EXPECTED_SNPS <- 41
SEED   <- 12345
K_RANGE <- 2:3       # tight, fast
MAX_PCS <- 10        # cap PCs to keep PCA fast
VAR_TARGET <- 0.80   # target variance for PCA

## ======================= FAST HELPERS =======================================
read_any_table <- function(path){
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) readr::read_csv(path, show_col_types = FALSE)
  else if (ext %in% c("tsv","txt")) readr::read_tsv(path, show_col_types = FALSE)
  else NULL  # skip Excel for speed/stability
}
std_id <- function(df){
  if (is.null(df) || !is.data.frame(df)) return(df)
  keys <- c("IID","iid","patient_id","subject_id","SubjectID","SampleID","sample","ID","Id","id")
  hit  <- keys[keys %in% names(df)]
  if (!length(hit)) return(df)
  df$IID <- as.character(df[[hit[1]]]); df
}
dedup_IID <- function(df){
  if (!"IID" %in% names(df)) return(tibble(IID=character()))
  df %>% filter(!is.na(IID) & IID != "") %>% distinct(IID, .keep_all = TRUE)
}
prep_geno <- function(G){
  G <- as.matrix(G); storage.mode(G) <- "double"
  m <- colMeans(G, na.rm = TRUE)
  for (j in seq_len(ncol(G))) if (anyNA(G[,j])) G[is.na(G[,j]), j] <- m[j]
  s <- apply(G, 2, sd); s[s == 0] <- 1
  G <- sweep(G, 2, colMeans(G), "-"); G <- sweep(G, 2, s, "/")
  G
}
choose_n_pcs <- function(pca, target=VAR_TARGET, max_pcs=MAX_PCS){
  ve <- (pca$sdev^2) / sum(pca$sdev^2)
  n  <- which(cumsum(ve) >= target)[1]; if (is.na(n)) n <- length(ve)
  min(max(n, 2), max_pcs)
}
best_k_elbow <- function(X, ks, seed=SEED){
  set.seed(seed)
  fits <- lapply(ks, function(k) kmeans(X, centers=k, nstart=25, iter.max=200))
  wss  <- sapply(fits, `[[`, "tot.withinss")
  drops <- -diff(wss) / wss[-length(wss)]
  idx <- which.max(c(NA, drops))
  list(k=ks[idx], km=fits[[idx]])
}

## ======================= CLINICAL (CSV/TSV only; cached) ====================
clin_cache <- file.path(RUN_DIR, "clinical_merged_dedup.csv")
if (file.exists(clin_cache)) {
  clin_df <- suppressWarnings(readr::read_csv(clin_cache, show_col_types = FALSE))
} else {
  files <- list.files(DAT_DIR, full.names = TRUE, recursive = FALSE)
  lst <- files %>% map(read_any_table) %>% map(std_id) %>% purrr::compact()
  clin_df <- if (length(lst)) Reduce(function(a,b){
    if (!"IID" %in% names(a)) a$IID <- NA_character_
    if (!"IID" %in% names(b)) b$IID <- NA_character_
    suppressMessages(full_join(a,b, by="IID"))
  }, lst) else tibble(IID=character())
  clin_df <- dedup_IID(clin_df)
  readr::write_csv(clin_df, clin_cache)
}
