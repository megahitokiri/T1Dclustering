#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(data.table)
})

## ======================= CONFIG (read-only step, v2) ========================
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

CLINICAL_PATH <- file.path(DAT_DIR, "cluster_clinical.csv")  # only this file

## ======================= helpers ===========================================
extract_IIDs <- function(obj){
  # obj is list(G, fam, map) as saved by your builder script
  # 1) prefer rownames(G) if present
  if (!is.null(obj$G) && !is.null(rownames(obj$G))) {
    return(as.character(rownames(obj$G)))
  }
  # 2) else use fam$iid (or column 2 if needed)
  if (!is.null(obj$fam)) {
    if ("iid" %in% names(obj$fam)) return(as.character(obj$fam$iid))
    if (ncol(obj$fam) >= 2) return(as.character(obj$fam[[2]]))
  }
  # 3) fallback: seq_len(nrow(G))
  if (!is.null(obj$G)) return(as.character(seq_len(nrow(obj$G))))
  character()
}

matrix_dims <- function(obj){
  if (is.null(obj$G)) return(c(NA_integer_, NA_integer_))
  c(nrow(obj$G), ncol(obj$G))
}

## ======================= 1) Read clinical (minimal) ========================
if (!file.exists(CLINICAL_PATH)) {
  stop(sprintf("Clinical file not found at: %s", CLINICAL_PATH))
}

message("[1] Loading clinical data: ", CLINICAL_PATH)
clinical <- suppressMessages(readr::read_csv(CLINICAL_PATH, show_col_types = FALSE))

message(sprintf("Clinical rows: %d; columns: %d", nrow(clinical), ncol(clinical)))
message("Clinical column names:")
message(paste(names(clinical), collapse = ", "))

preview_cols <- intersect(c("Cluster_subjectID","Cluster_age","Cluster_sex",
                            "Cluster_weight","Cluster_height","Cluster_A1c"),
                          names(clinical))
if (length(preview_cols) == 0) preview_cols <- names(clinical)[seq_len(min(10L, ncol(clinical)))]
clinical_preview <- head(clinical[, preview_cols, drop = FALSE], 10)
fwrite(as.data.table(clinical_preview), file.path(RUN_DIR, "step1v2_clinical_preview.tsv"), sep = "\t")

## ======================= 2) Read genotype RDS (A/B/C) ======================
summary_list <- list()

for (ds in names(RDS_FILES)) {
  rds_path <- RDS_FILES[[ds]]
  message("\n=== Dataset ", ds, " ===")
  if (!file.exists(rds_path)) {
    warning(sprintf("RDS for dataset %s not found: %s. Skipping.", ds, rds_path))
    next
  }
  obj <- readRDS(rds_path)
  if (!is.list(obj) || is.null(obj$G)) {
    warning(sprintf("RDS %s is not the expected list(G, fam, map).", rds_path))
  }
  dims <- matrix_dims(obj)
  IIDs <- extract_IIDs(obj)
  n_samples <- dims[1]; n_snps <- dims[2]
  
  message(sprintf("Loaded %s: samples=%s; SNPs=%s", ds, n_samples, n_snps))
  message(sprintf("IIDs source: %s",
                  if (!is.null(rownames(obj$G))) "rownames(G)" else if (!is.null(obj$fam)) "fam$iid" else "generated sequence"))
  if (length(IIDs)) {
    message(sprintf("Example IIDs (%s): %s", ds, paste(head(IIDs, 5), collapse = ", ")))
  }
  
  # Clinical alignment check (no merging, just counts)
  if (!"Cluster_subjectID" %in% names(clinical)) {
    warning("Clinical file lacks 'Cluster_subjectID'; cannot align. Skipping alignment counts.")
    n_overlap <- NA_integer_
    n_missing <- NA_integer_
  } else {
    common_ids <- intersect(IIDs, as.character(clinical$Cluster_subjectID))
    n_overlap <- length(common_ids)
    n_missing <- length(setdiff(IIDs, as.character(clinical$Cluster_subjectID)))
    # Save only the matched minimal clinical subset for this dataset preview
    if (n_overlap > 0) {
      clin_min_cols <- intersect(c("Cluster_subjectID","Cluster_age","Cluster_sex",
                                   "Cluster_weight","Cluster_height","Cluster_A1c"),
                                 names(clinical))
      clin_match <- clinical[clinical$Cluster_subjectID %in% common_ids, clin_min_cols, drop = FALSE]
      fwrite(as.data.table(head(clin_match, 20)), file.path(RUN_DIR, paste0("step1v2_", ds, "_clin_match_preview.tsv")), sep = "\t")
    }
  }
  
  summary_list[[ds]] <- data.table(
    dataset = ds,
    n_samples = n_samples,
    n_snps = n_snps,
    iids_available = length(IIDs) > 0,
    overlap_with_clinical = n_overlap,
    missing_clinical = n_missing
  )
}

## ======================= 3) Write a tiny summary ===========================
if (length(summary_list)) {
  sum_dt <- rbindlist(summary_list, fill = TRUE)
  out_file <- file.path(RUN_DIR, "step1v2_read_summary.tsv")
  fwrite(sum_dt, out_file, sep = "\t")
  message("\nSummary written to: ", out_file)
} else {
  message("\nNo genotype datasets were read.")
}

message("\nStep 1 v2 complete: data read only, no merges, no QC.")
