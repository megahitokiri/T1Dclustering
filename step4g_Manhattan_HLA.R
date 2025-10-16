## ==========================================================
## Ancestry × Cluster: robust join, diagnostics, and plots
## Grid: B02500 (change BUCKET if you switch grids)
## ==========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glue)
})

## ---------- Paths (edit BASE/BUCKET only if needed) ----------
BASE      <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
BUCKET    <- file.path(BASE, "grid_B02500")   # change e.g. to grid_B01000 if desired
VIZDIR    <- file.path(BUCKET, "viz_summary")
ANC_RAW   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/data/annotated_ancestry_ADDAM2.csv"
ANC_CLEAN <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/data/annotated_ancestry_ADDAM2.cleaned.csv"
MERGED_IN <- file.path(VIZDIR, "merged_pca_clusters_clinical.tsv")  # used if present

dir.create(VIZDIR, showWarnings = FALSE, recursive = TRUE)

## ---------- Helpers ----------
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

read_try <- function(fp, ...) {
  if (!file.exists(fp)) return(NULL)
  tryCatch(fread(fp, ...), error = function(e) NULL)
}

normalize_cluster_col <- function(dt) {
  if (is.null(dt)) return(NULL)
  if (!("cluster" %in% names(dt))) {
    cand <- c("kmeans_cluster","pam_cluster","gmm_cluster","Cluster","clusters","label")
    hit  <- cand[cand %in% names(dt)][1]
    if (!is.na(hit)) setnames(dt, hit, "cluster")
  }
  dt
}

# ---- K-MEANS ONLY: rebuild merged table by scanning bucket ----
rebuild_merged_from_bucket <- function(bucket) {
  # 1) find only step3_kmeans_*_<A|B|C> directories anywhere under the bucket
  all_dirs <- list.dirs(bucket, recursive = TRUE, full.names = TRUE)
  keep_dirs <- all_dirs[grepl("^step3_kmeans_.+_[ABC]$", basename(all_dirs))]
  if (!length(keep_dirs)) return(NULL)
  
  # 2) from those dirs, take ONLY clusters.tsv
  cl_files <- unlist(lapply(keep_dirs, function(d)
    list.files(d, pattern = "^clusters\\.tsv$", full.names = TRUE, recursive = FALSE)
  ))
  cl_files <- cl_files[file.exists(cl_files)]
  if (!length(cl_files)) return(NULL)
  
  extract_df <- function(f) {
    # expected parent dir name: step3_kmeans_<feature>_<dataset>
    parent <- basename(dirname(f))
    m <- regexec("^step3_kmeans_(.+)_([ABC])$", parent)
    mm <- regmatches(parent, m)[[1]]
    if (length(mm) != 3) return(NULL)
    feature <- mm[2]; dataset <- mm[3]
    
    dt <- tryCatch(data.table::fread(f), error = function(e) NULL)
    if (is.null(dt)) return(NULL)
    
    # ID column → SubjectID
    idc <- c("SubjectID","IID","Cluster_subjectID","id","sample")
    idc <- idc[idc %in% names(dt)][1]
    if (is.na(idc)) return(NULL)
    data.table::setnames(dt, idc, "SubjectID")
    dt[, SubjectID := toupper(gsub("[^A-Za-z0-9]", "", as.character(SubjectID)))]
    
    # cluster column normalization
    if (!("cluster" %in% names(dt))) {
      cand <- c("kmeans_cluster","pam_cluster","gmm_cluster","Cluster","clusters","label")
      hit  <- cand[cand %in% names(dt)][1]
      if (!is.na(hit)) data.table::setnames(dt, hit, "cluster")
    }
    if (!("cluster" %in% names(dt))) return(NULL)
    
    keep <- intersect(c("SubjectID","cluster","stability","is_outlier"), names(dt))
    out  <- dt[, ..keep]
    out[, `:=`(dataset = dataset, feature = feature, algo = "kmeans")]
    out
  }
  
  pieces <- lapply(cl_files, extract_df)
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(NULL)
  
  data.table::rbindlist(pieces, fill = TRUE, use.names = TRUE)
}

## ---------- 1) Load ancestry (prefer cleaned if available) ----------
anc <- read_try(ANC_CLEAN)
if (is.null(anc)) {
  anc <- read_try(ANC_RAW)
  if (is.null(anc)) stop("Cannot read ancestry file at: ", ANC_RAW)
  
  # Standardize headers from your file: IID, Predicted_Ancestry, Pop, Color, PC1, PC2, PC3, PC4
  nm <- names(anc)
  if ("IID" %in% nm && !"SubjectID" %in% nm) setnames(anc, "IID", "SubjectID")
  if ("Predicted_Ancestry" %in% nm && !"Ancestry" %in% nm) setnames(anc, "Predicted_Ancestry", "Ancestry")
  
  need <- c("SubjectID","Ancestry")
  miss <- setdiff(need, names(anc))
  if (length(miss)) stop("Ancestry file missing columns: ", paste(miss, collapse=", "))
  
  anc[, SubjectID := canon_id(SubjectID)]
  anc[, Ancestry  := as.character(Ancestry)]
  
  fwrite(anc, ANC_CLEAN)
}

# Make a join key

drop0    <- function(x) sub("^0+", "", x)
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

## 1) Ancestry side: JoinID = k_drop0(SubjectID)
stopifnot("SubjectID" %in% names(anc))
anc[, JoinID := drop0(canon_id(SubjectID))]
setkey(anc, JoinID)

## ---------- 2) Load merged clusters (or rebuild) ----------
## ---------- 2) Load merged clusters (k-means ONLY) ----------
# Always prefer a fresh rebuild so we never accidentally include PAM
merged <- rebuild_merged_from_bucket(BUCKET)
if (is.null(merged)) stop("No k-means clusters.tsv found under: ", BUCKET)

# (Paranoid) keep only kmeans if an 'algo' column exists
if ("algo" %in% names(merged)) merged <- merged[algo == "kmeans"]

# De-duplicate: one row per (dataset, feature, SubjectID)
if ("stability" %in% names(merged)) {
  data.table::setorder(merged, dataset, feature, SubjectID, -stability)
} else {
  data.table::setorder(merged, dataset, feature, SubjectID)
}
merged <- merged[!duplicated(merged, by = c("dataset","feature","SubjectID"))]

# Standardize ID; create JoinID that matches ancestry’s drop0(canon_id(SubjectID))
drop0    <- function(x) sub("^0+", "", x)
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
merged[, SubjectID := canon_id(SubjectID)]
merged[, JoinID    := drop0(SubjectID)]
merged <- normalize_cluster_col(merged)
if (!("algo" %in% names(merged))) merged[, algo := "kmeans"]

# Basic required columns check
reqm <- c("dataset","feature","cluster","SubjectID","JoinID","algo")
miss2 <- setdiff(reqm, names(merged))
if (length(miss2)) stop("Merged table missing columns: ", paste(miss2, collapse=", "))

# Optionally cache the clean k-means-only merge for later runs
data.table::fwrite(merged, MERGED_IN, sep = "\t")
## ---------- 3) Diagnostics: header & overlap snapshot ----------
diag_tbl <- data.table(
  object = c("anc cols", "merged cols"),
  value  = c(paste(names(anc), collapse = ", "),
             paste(names(merged), collapse = ", "))
)
fwrite(diag_tbl, file.path(VIZDIR, "ancestry_merge_headers.tsv"), sep = "\t")

# Overlap counts
in_both   <- length(intersect(anc$JoinID, merged$JoinID))
only_anc  <- length(setdiff(anc$JoinID, merged$JoinID))
only_clus <- length(setdiff(merged$JoinID, anc$JoinID))
ovl_line  <- data.table(
  metric = c("N_anc","N_merged","N_overlap","N_only_anc","N_only_clusters"),
  value  = c(nrow(anc), nrow(merged), in_both, only_anc, only_clus)
)
fwrite(ovl_line, file.path(VIZDIR, "ancestry_merge_overlap_counts.tsv"), sep = "\t")

# Small samples for manual inspection
samp <- data.table(
  anc_example    = head(anc$JoinID, 20),
  merged_example = head(merged$JoinID, 20)
)
fwrite(samp, file.path(VIZDIR, "ancestry_merge_id_examples.tsv"), sep = "\t")

## ---------- 4) Smart left-join (clusters ⟵ ancestry) ----------
setDT(anc); setDT(merged)
# Keep only the columns we need from anc (avoid accidental name clashes)
anck <- anc[, .(JoinID, Ancestry, PC1, PC2)]
setkeyv(anck, "JoinID")
setkeyv(merged, "JoinID")

mf <- merge(merged, anck, by = "JoinID", all.x = TRUE)
# Keep SubjectID for readability
if (!("SubjectID" %in% names(mf))) mf[, SubjectID := JoinID]
# Make cluster a factor for nicer faceting
mf[, cluster := as.factor(cluster)]

# Save the merged table used for plotting
fwrite(mf, file.path(VIZDIR, "merged_with_ancestry.tsv"), sep = "\t")

## ---------- 5) Composition tables & plots ----------
comp <- mf[, .N, by = .(dataset, feature, algo, cluster, Ancestry)]
# enforce facet column order
comp[, feature := factor(feature, levels = c("genetic", "genetic_minclin", "genetic_fullclin"))]

# Protect against all-NA ancestry: replace NA with "Unknown" for plotting
comp[is.na(Ancestry), Ancestry := "Unknown"]
comp[, prop := N / sum(N), by = .(dataset, feature, algo, cluster)]
fwrite(comp, file.path(VIZDIR, "cluster_ancestry_counts.tsv"), sep = "\t")

# Stacked proportions
p_stack <- ggplot(comp, aes(x = cluster, y = prop, fill = Ancestry)) +
  geom_col(width = 0.85, color = "gray30") +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3, color = "black") +
  facet_grid(dataset + algo ~ feature) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Ancestry composition per cluster",
       x = "Cluster", y = "Proportion") +
  theme_minimal(base_size = 12)
ggsave(file.path(VIZDIR, "cluster_ancestry_stackedbars.png"), p_stack, width = 13, height = 8.5, dpi = 160)

# Stacked counts
p_counts <- ggplot(comp, aes(x = cluster, y = N, fill = Ancestry)) +
  geom_col(width = 0.85, color = "gray30") +
  geom_text(aes(label = N),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  facet_grid(dataset + algo ~ feature) +
  labs(title = "Ancestry counts per cluster",
       x = "Cluster", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(file.path(VIZDIR, "cluster_ancestry_counts.png"), p_counts, width = 13, height = 8.5, dpi = 160)

## ---------- 6) Association tests (Ancestry vs Cluster) ----------
run_assoc_safe <- function(df) {
  # Drop NA ancestry
  df <- df[!is.na(Ancestry)]
  out <- data.table(
    dataset = if (nrow(df)) unique(df$dataset)[1] else NA_character_,
    feature = if (nrow(df)) unique(df$feature)[1] else NA_character_,
    algo    = if (nrow(df)) unique(df$algo)[1]    else NA_character_,
    method  = NA_character_,
    p_value = as.numeric(NA),
    note    = NA_character_
  )
  if (nrow(df) == 0) { out$note <- "no rows after NA filtering"; return(out) }
  
  tab <- table(df$Ancestry, df$cluster)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    out$note <- sprintf("degenerate table: %dx%d", nrow(tab), ncol(tab))
    return(out)
  }
  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    out$note <- "zero row/col in contingency table"
    return(out)
  }
  res <- tryCatch(fisher.test(tab), error = function(e) tryCatch(chisq.test(tab), error = function(e2) NULL))
  if (inherits(res, "htest")) {
    out$method  <- res$method
    out$p_value <- res$p.value
  } else {
    out$note <- "test failed"
  }
  out
}

keys <- unique(mf[, .(dataset, feature, algo)])
assoc_list <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  sl <- mf[dataset == keys$dataset[i] & feature == keys$feature[i] & algo == keys$algo[i]]
  assoc_list[[i]] <- run_assoc_safe(sl)
}
assoc <- rbindlist(assoc_list, use.names = TRUE, fill = TRUE)
setorder(assoc, p_value)
fwrite(assoc, file.path(VIZDIR, "ancestry_cluster_association.tsv"), sep = "\t")

## ---------- 7) Optional PCA overlay if PC1/PC2 exist ----------
if (all(c("PC1","PC2") %in% names(mf))) {
  p_pca <- ggplot(mf, aes(PC1, PC2, color = Ancestry, shape = cluster)) +
    geom_point(alpha = 0.9, size = 1.8) +
    facet_grid(dataset + algo ~ feature) +
    labs(title = "PC1 vs PC2 colored by ancestry (shape = cluster)",
         shape = "Cluster") +
    theme_minimal(base_size = 12)
  ggsave(file.path(VIZDIR, "pca_ancestry_overlay.png"), p_pca, width = 13, height = 8.5, dpi = 160)
}

cat(glue("
✔ Wrote to: {VIZDIR}
  - merged_with_ancestry.tsv
  - ancestry_merge_headers.tsv
  - ancestry_merge_overlap_counts.tsv
  - ancestry_merge_id_examples.tsv
  - cluster_ancestry_counts.tsv
  - cluster_ancestry_stackedbars.png
  - cluster_ancestry_counts.png
  - ancestry_cluster_association.tsv
  {if (all(c('PC1','PC2') %in% names(mf))) ' - pca_ancestry_overlay.png' else ''}
"))

## ---------- 8) Summary table: ancestry composition per cluster ----------
summary_table <- comp[, .(
  Count = sum(N),
  Proportion = sprintf("%.1f%%", 100 * mean(prop, na.rm = TRUE))
), by = .(dataset, feature, cluster, Ancestry)]

# Ensure cluster is numeric for sorting
summary_table[, cluster := as.numeric(as.character(cluster))]

# Order the columns logically
setcolorder(summary_table, c("dataset", "feature", "cluster", "Ancestry", "Count", "Proportion"))

# Sort nicely for readability
setorder(summary_table, dataset, feature, cluster, Ancestry)

# Write it to disk
out_summary <- file.path(VIZDIR, "cluster_ancestry_summary.tsv")
fwrite(summary_table, out_summary, sep = "\t")

cat(glue("
✔ Summary table written to:
  - {out_summary}
"))