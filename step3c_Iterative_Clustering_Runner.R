## ==================== STEP 3 GRID RUNNER (compact & robust) ====================

## Project root + load step3 code
setwd("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/")
source("step3a_clustering_rstudio.R")

## Where step3 writes and where grid buckets live
ROOT_WRITE <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
BASE_GRID  <- ROOT_WRITE

## Grid + knobs
BOOT_GRID <- c(1, 10, 100, 250, 500, 1000, 1500, 2000,2500,3000)   # Max Bootstrap 3000
DATASETS  <- c("A","B","C")
N_PCS     <- 15
K_RANGE   <- 2:6

## ---- Helpers (minimal messages) ----
list_top_step3_dirs <- function() {
  kids <- list.dirs(ROOT_WRITE, recursive = FALSE, full.names = TRUE)
  kids[grepl("^step3_(kmeans|pam|gmm)_", basename(kids))]
}

clean_top_step3_dirs <- function() {
  dd <- list_top_step3_dirs()
  if (length(dd)) unlink(dd, recursive = TRUE, force = TRUE)
  invisible(length(dd))
}

copy_dir_contents <- function(src_dir, dst_dir) {
  if (!dir.exists(src_dir)) return(FALSE)
  if (dir.exists(dst_dir)) unlink(dst_dir, recursive = TRUE, force = TRUE)
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  items <- list.files(src_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(items)) return(FALSE)
  ok <- file.copy(items, dst_dir, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
  all(ok)
}

verify_bucket <- function(bucket) {
  hits <- list.files(bucket, pattern = "^clusters\\.tsv$", recursive = TRUE, full.names = TRUE)
  length(hits)
}

## ---- One run for a given BOOT_N ----
run_one_boot <- function(B) {
  bucket <- file.path(BASE_GRID, sprintf("grid_B%05d", B))
  dir.create(bucket, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n[Step3] BOOT_N=%d  → bucket: %s\n", B, bucket))
  
  # 1) Clean previous top-level step3_* in ROOT_WRITE
  n_del <- clean_top_step3_dirs()
  cat(sprintf("  - cleaned old result dirs: %d\n", n_del))
  
  # 2) Configure + run step3 (KMeans only)
  cfg <- default_step3_config()
  cfg$DATASETS   <- DATASETS
  cfg$N_PCS      <- N_PCS
  cfg$K_RANGE    <- K_RANGE
  cfg$BOOT_N     <- B
  cfg$RUN_KMEANS <- TRUE
  cfg$RUN_PAM    <- FALSE
  cfg$RUN_GMM    <- FALSE
  
  t0 <- Sys.time()
  run_step3(cfg)
  cat(sprintf("  - run_step3 done in %.1fs\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
  
  # 3) Collect newly written top-level step3_* dirs and copy to bucket
  produced <- list_top_step3_dirs()
  cat(sprintf("  - result dirs found: %d\n", length(produced)))
  
  moved_ok <- 0L
  for (sd in produced) {
    dst <- file.path(bucket, basename(sd))
    moved_ok <- moved_ok + as.integer(copy_dir_contents(sd, dst))
  }
  cat(sprintf("  - result dirs copied: %d\n", moved_ok))
  
  # 4) Sanity check: bucket has clusters.tsv files?
  n_clusters <- verify_bucket(bucket)
  cat(sprintf("  - clusters.tsv in bucket: %d\n", n_clusters))
  
  invisible(list(bucket = bucket, produced = produced, copied = moved_ok, n_clusters = n_clusters))
}

## ---- Run the grid (stops on error is OK; comment tryCatch if you prefer strict) ----
for (B in BOOT_GRID) {
  tryCatch(run_one_boot(B), error = function(e) {
    cat(sprintf("  ! ERROR at BOOT_N=%d: %s\n", B, conditionMessage(e)))
  })
}

## Quick peek: list created buckets
print(list.dirs(BASE_GRID, full.names = TRUE, recursive = FALSE))