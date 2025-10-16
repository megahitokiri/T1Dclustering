# =========================
# Minimal, robust merge fix
# (IDs only; safe for odd RDS)
# =========================
suppressPackageStartupMessages({
  library(data.table)
})

# ---- Paths ----
ROOT   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs"
GRID   <- "grid_B02500"
BUCKET <- file.path(ROOT, "cluster_runs_min", GRID)
VIZ    <- file.path(BUCKET, "viz_summary")
dir.create(VIZ, recursive = TRUE, showWarnings = FALSE)

# ---- Canonical join key (same as before) ----
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
drop0    <- function(x) sub("^0+", "", x)
mk_join  <- function(x) drop0(canon_id(x))

# ---- Extract IDs from many possible RDS shapes or fallback to .fam ----
load_geno_ids_ds <- function(ds) {
  rds_fp <- file.path(ROOT, paste0(ds, "_homog_numeric.rds"))
  fam_fp <- file.path(ROOT, paste0(ds, "_homog.fam"))
  
  ids <- NULL
  if (file.exists(rds_fp)) {
    X <- readRDS(rds_fp)
    
    # 1) Row names on any rectangular object
    rn <- try(rownames(X), silent = TRUE)
    if (!inherits(rn, "try-error") && !is.null(rn)) {
      ids <- rn
    }
    
    # 2) If it's a list, try common fields
    if (is.null(ids) && is.list(X)) {
      cand <- c("ids","IID","samples","rownames","row_ids")
      hit  <- cand[cand %in% names(X)][1]
      if (!is.na(hit)) {
        ids <- as.character(X[[hit]])
      }
      # some lists store a matrix/data.frame under $X / $geno with rownames
      if (is.null(ids)) {
        mat_cand <- c("X","geno","G","mat","data")
        mat_hit  <- mat_cand[mat_cand %in% names(X)][1]
        if (!is.na(mat_hit)) {
          rn2 <- try(rownames(X[[mat_hit]]), silent = TRUE)
          if (!inherits(rn2, "try-error") && !is.null(rn2)) ids <- rn2
        }
      }
    }
    
    # 3) If it’s a Matrix::dgCMatrix or similar with Dimnames
    if (is.null(ids)) {
      dn <- try(attr(X, "Dimnames"), silent = TRUE)
      if (!inherits(dn, "try-error") && !is.null(dn) && length(dn) >= 1 && !is.null(dn[[1]])) {
        ids <- dn[[1]]
      }
    }
  }
  
  # 4) Fallback to PLINK .fam (V2 = IID)
  if (is.null(ids)) {
    if (!file.exists(fam_fp)) {
      stop("Could not recover sample IDs for ", ds,
           " (no rownames/ids in RDS and no .fam at ", fam_fp, ").")
    }
    fam <- fread(fam_fp, header = FALSE)
    if (ncol(fam) < 2) stop(".fam malformed for ", ds, ": ", fam_fp)
    ids <- fam[[2]]
  }
  
  # Return a compact DT with IID_raw + JoinID
  dt <- data.table(IID_raw = as.character(ids))
  dt[, JoinID := mk_join(IID_raw)]
  unique(dt[])
}

# ---- Load cluster table and select a sensible cluster column ----
load_cluster_ds <- function(ds,
                            prefer_order = c("cluster_genetic_fullclin",
                                             "cluster_genetic",
                                             "cluster_genetic_minclin",
                                             "cluster_clinical_only")) {
  cand <- c(
    file.path(VIZ, paste0("merged_assignments_", ds, ".tsv")),
    file.path(VIZ, paste0(GRID, "_merged_assignments_", ds, ".tsv")),
    file.path(VIZ, paste0(ds, "_merged_assignments.tsv"))
  )
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("No merged_assignments file found for ", ds, " under ", VIZ)
  
  clu <- fread(hit)
  
  # If JoinID exists, normalize it; otherwise build from any IID-like column
  if ("JoinID" %in% names(clu)) {
    clu[, JoinID := mk_join(JoinID)]
  } else {
    id_cand <- c("IID","IID_raw","IID_raw.x","IID_raw.y","SubjectID","Cluster_subjectID")
    id_cand <- id_cand[id_cand %in% names(clu)]
    if (!length(id_cand)) stop("No ID-like column found in: ", hit)
    src <- id_cand[1]
    setnames(clu, src, "IID_tmp", skip_absent = TRUE)
    clu[, JoinID := mk_join(IID_tmp)]
  }
  
  # Pick a cluster_* column by preference
  cl_cols <- intersect(prefer_order, names(clu))
  if (!length(cl_cols)) stop("No cluster_* columns present in: ", hit)
  chosen <- cl_cols[1L]
  
  out <- clu[, .(JoinID, cluster = get(chosen))]
  out <- out[!is.na(cluster)]
  attr(out, "cluster_col") <- chosen
  attr(out, "file") <- hit
  out[]
}

# ---- One dataset merge (reports + writes small TSV) ----
merge_one <- function(ds) {
  message("\n>>> Dataset ", ds)
  geno_ids <- load_geno_ids_ds(ds)   # IDs only, robust to RDS shape
  message("  - IDs recovered: ", nrow(geno_ids))
  
  clu  <- load_cluster_ds(ds)
  chosen_col <- attr(clu, "cluster_col")
  message("  - Cluster file: ", attr(clu, "file"))
  message("  - Using cluster column: ", chosen_col)
  
  # Merge by JoinID
  m <- merge(geno_ids[, .(JoinID, IID_raw)], clu, by = "JoinID", all.x = FALSE, all.y = TRUE)
  message("  - Overlap n = ", nrow(m), " (cluster rows = ", nrow(clu), ")")
  
  # Per-cluster counts (in overlap)
  cnt <- m[, .N, by = cluster][order(cluster)]
  fwrite(cnt, file.path(VIZ, paste0("merge_counts_", ds, ".tsv")), sep = "\t")
  fwrite(m,   file.path(VIZ, paste0("merge_ids_", ds, ".tsv")),    sep = "\t")
  cnt[]
}

# ---- Run for A/B/C ----
invisible(lapply(c("A","B","C"), merge_one))
cat("\n✔ Merge checks done. See TSVs in:\n  - ", VIZ, "\n",
    "  * merge_counts_*.tsv (counts by cluster)\n",
    "  * merge_ids_*.tsv (overlap IDs)\n", sep = "")