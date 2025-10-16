# ---- CONFIG ----
ROOT   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs"
RAWF   <- file.path(ROOT, "rs9271351_geno.raw")   # your PLINK --recode A .raw
SNP_POS <- "06:032616083"                         # GRCh38 position
RSID    <- "rs9271351"                            # optional, used as fallback
DS_LIST <- c("A","B","C")                         # choose one if you like, e.g. "C"

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# Normalize "06:032616083" -> "6:32616083"
pos_std <- function(x) {
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  parts <- unlist(strsplit(x, "[:_]", perl = TRUE))
  chr <- gsub("^0+", "", parts[1])
  bp  <- gsub("^0+", "", parts[2])
  paste0(chr, ":", bp)
}

# Read A1/A2 from the BIM of a dataset (e.g., C_homog.bim)
get_a1a2 <- function(ds, pos = SNP_POS) {
  bimf <- file.path(ROOT, paste0(ds, "_homog.bim"))
  stopifnot(file.exists(bimf))
  bim <- fread(bimf, col.names = c("CHR","SNP","CM","BP","A1","A2"))
  pstd <- pos_std(pos)
  # match either zero-padded "06:032616083" or normalized "6:32616083"
  hit <- bim %>% filter(SNP %in% c(pos, pstd) | paste0(CHR, ":", BP) == pstd)
  if (nrow(hit) != 1) stop("Could not uniquely find SNP in ", bimf, " (n=", nrow(hit), ")")
  list(A1 = hit$A1[1], A2 = hit$A2[1], CHR = hit$CHR[1], BP = hit$BP[1], SNP = hit$SNP[1])
}

# Find the SNP column name in .raw (handles several header styles)
find_raw_col <- function(headers, pos = SNP_POS, rsid = RSID, a1, a2) {
  pstd <- pos_std(pos)          # "6:32616083"
  pos_pat <- gsub(":", "[:_]", pstd)  # match ":" or "_" between chr/bp in headers
  
  # Common patterns:
  #  - "6:32616083_C_G" (PLINK 1.9 default)
  #  - "6:32616083_C"    (your file; dosage of A1 only)
  #  - "rs9271351_C_G"   (rsid form)
  #  - With zero padding "06:032616083_*"
  pats <- c(
    paste0("^", pos_pat, "_", a1, "_", a2, "$"),
    paste0("^", pos_pat, "_", a2, "_", a1, "$"),
    paste0("^", pos_pat, "_", a1, "$"),
    paste0("^", pos_pat, "_", a2, "$"),
    paste0("^", rsid, "_", a1, "_", a2, "$"),
    paste0("^", rsid, "_", a2, "_", a1, "$"),
    paste0("^", rsid, "_", a1, "$"),
    paste0("^", rsid, "_", a2, "$"),
    paste0("^", gsub(":", ":", pos, fixed = TRUE), "_", a1, "$"),  # "06:032616083_C"
    paste0("^", gsub(":", ":", pos, fixed = TRUE), "_", a2, "$")
  )
  
  hits <- Reduce(union, lapply(pats, function(p) grep(p, headers, value = TRUE)))
  if (!length(hits)) {
    # last resort: contains bp or rsid and ends with _A1/_A2
    hits <- grep(paste0("(32616083|", rsid, ").*(_", a1, "$|_", a2, "$|_", a1, "_", a2, "$|_", a2, "_", a1, "$)"),
                 headers, value = TRUE)
  }
  if (length(hits) == 0) {
    stop("Cannot find SNP column in .raw. Headers tried to match with A1/A2 (",
         a1, "/", a2, ").")
  }
  if (length(hits) > 1) {
    # pick the one ending with _A1 if present, else first (and warn)
    a1_last <- grep(paste0("_", a1, "$"), hits, value = TRUE)
    if (length(a1_last) == 1) return(a1_last)
    warning("Multiple header matches; using first: ", hits[1], "\nAll: ",
            paste(hits, collapse = ", "))
    return(hits[1])
  }
  hits[1]
}

# Load genotype dosage + call genotypes
# If the column is "..._C", values are 0/1/2 copies of C (A1).
# We then map: 2 -> CC, 1 -> CG, 0 -> GG (assuming A2=G).
load_raw_geno <- function(ds, rawf = RAWF, pos = SNP_POS, rsid = RSID) {
  a1a2 <- get_a1a2(ds, pos)
  A1 <- a1a2$A1; A2 <- a1a2$A2
  
  raw <- fread(rawf)
  snp_col <- find_raw_col(names(raw), pos = pos, rsid = rsid, a1 = A1, a2 = A2)
  message("Using .raw column: ", snp_col, " (interpreted as dosage of ", A1, ")")
  
  out <- raw %>%
    transmute(IID        = IID,
              dosage_A1  = .data[[snp_col]]) %>%
    mutate(genotype = dplyr::case_when(
      dosage_A1 == 2 ~ paste0(A1, A1),
      dosage_A1 == 1 ~ paste0(A1, A2),
      dosage_A1 == 0 ~ paste0(A2, A2),
      TRUE ~ NA_character_
    ),
    IID_canon = canon_id(IID))
  
  # quick sanity
  cat("\nGenotype counts (", ds, "):\n", sep = "")
  print(table(out$genotype, useNA = "ifany"))
  invisible(out)
}

# ---- EXAMPLE USAGE ----
# Pick one dataset (the BIM is just used to read A1/A2; the .raw is common)
geno_C <- load_raw_geno("C")   # will print column chosen + genotype counts
# head(geno_C)

# If you want **proportions by cluster** for one set of cluster results:
# (adjust to your bucket and feature/dataset)
merge_with_clusters <- function(ds = "C", feature = "genetic", boot = 2500) {
  tag <- sprintf("grid_B%05d", boot)
  clustf <- file.path(ROOT, "cluster_runs_min", tag,
                      paste0("step3_kmeans_", ds), feature, "clusters.tsv")
  stopifnot(file.exists(clustf))
  clust <- fread(clustf)
  idcol <- c("IID","SubjectID","id","sample")
  idcol <- idcol[idcol %in% names(clust)][1]
  clust <- clust %>% rename(ID = !!idcol) %>% mutate(IID_canon = canon_id(ID))
  
  geno <- load_raw_geno(ds)
  m <- clust %>% left_join(geno %>% select(IID_canon, genotype), by = "IID_canon")
  
  cat("\nGenotype by cluster (", ds, " / ", feature, " / B=", boot, "):\n", sep = "")
  print(with(m, prop.table(table(Cluster, genotype), 1)))
  invisible(m)
}

# Example:
# mC <- merge_with_clusters(ds="C", feature="genetic", boot=2500)

# ---------- settings ----------
BASE     <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs"
GRID_TAG <- "grid_B02500"

FEATURE_FOR_DATASET_VENN <- "genetic_minclin"
DATASET_FOR_FEATURE_VENN <- "C"

SAVE_DIR <- file.path(BASE, "cluster_runs_min", GRID_TAG, "summary_plots")
dir.create(SAVE_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(stringr); library(tidyr); library(ggplot2)
})

root_grid <- file.path(BASE, "cluster_runs_min", GRID_TAG)

# ---------- helper: discover files ----------
find_enrichment_files <- function() {
  in_grid <- c(
    Sys.glob(file.path(root_grid, "**", "enriched_snps.tsv")),
    Sys.glob(file.path(root_grid, "**", "top_snps*.tsv"))
  )
  fallbacks <- Sys.glob(file.path(BASE, "top_snps_driving_clusters_*.tsv"))
  unique(c(in_grid, fallbacks))
}

parse_ds_feat <- function(fp) {
  parts <- strsplit(fp, "/")[[1]]
  ds_idx <- which(grepl("^step3_kmeans_", parts))
  ds <- if (length(ds_idx)) sub("^step3_kmeans_", "", parts[ds_idx[1]]) else NA_character_
  feat <- if (length(ds_idx)) parts[ds_idx[1] + 1] else NA_character_
  
  if (is.na(ds) && grepl("top_snps_driving_clusters_([ABC])\\.tsv$", fp)) {
    ds <- sub(".*top_snps_driving_clusters_([ABC])\\.tsv$", "\\1", fp)
    feat <- "genetic"
  }
  list(dataset = ds, feature = feat)
}

read_enrichment_any <- function(fp) {
  x <- tryCatch(fread(fp), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  meta <- parse_ds_feat(fp)
  nm <- names(x)
  snp_col <- c("SNP","snp","variant","id","marker")
  snp_col <- snp_col[snp_col %in% nm][1]
  if (is.na(snp_col)) return(NULL)
  
  flag_col <- c("enriched_small","small_enriched","is_enriched","enriched")
  flag_col <- flag_col[flag_col %in% nm][1]
  
  if (!is.na(flag_col)) {
    y <- x %>% filter(.data[[flag_col]] %in% c(TRUE, 1, "TRUE", "yes", "Yes"))
  } else {
    eff <- c("odds_ratio","OR","logFC","effect","beta","abs_loading","loading")
    eff <- eff[eff %in% nm][1]
    if (!is.na(eff)) {
      y <- x %>% filter(is.finite(.data[[eff]])) %>% arrange(desc(.data[[eff]])) %>% slice_head(n = 50)
    } else {
      y <- x %>% slice_head(n = 50)
    }
  }
  if (!nrow(y)) return(NULL)
  y %>%
    transmute(SNP = gsub("\\s+","", .data[[snp_col]]),
              dataset = meta$dataset,
              feature = meta$feature,
              src_file = fp)
}

enrich_files <- find_enrichment_files()
if (!length(enrich_files)) stop("No enrichment/top_snps files found under grid or BASE.")

enrich <- bind_rows(lapply(enrich_files, read_enrichment_any)) %>%
  filter(!is.na(SNP) & nzchar(SNP) & !is.na(dataset) & nzchar(dataset) & !is.na(feature) & nzchar(feature)) %>%
  distinct(SNP, dataset, feature, .keep_all = TRUE)

if (!nrow(enrich)) stop("Could not parse any SNP rows from the files discovered.")

# ---------- safe print of available combinations ----------
avail <- enrich %>% distinct(dataset, feature) %>% arrange(dataset, feature)
cat("\nAvailable (dataset, feature) combos:\n")
print(as.data.frame(avail), row.names = FALSE)

# ---------- auto-fix requested comparisons if missing ----------
if (!(FEATURE_FOR_DATASET_VENN %in% avail$feature)) {
  FEATURE_FOR_DATASET_VENN <- avail$feature[1]
  message("Requested feature not found. Using feature = ", FEATURE_FOR_DATASET_VENN)
}
if (!(DATASET_FOR_FEATURE_VENN %in% avail$dataset)) {
  DATASET_FOR_FEATURE_VENN <- avail$dataset[1]
  message("Requested dataset not found. Using dataset = ", DATASET_FOR_FEATURE_VENN)
}

# ---------- build set lists ----------
make_set_list_by_dataset <- function(tbl, feature) {
  sub <- tbl %>% filter(feature == feature)
  split(unique(sub$SNP), sub$dataset)
}
make_set_list_by_feature <- function(tbl, dataset) {
  sub <- tbl %>% filter(dataset == dataset)
  split(unique(sub$SNP), sub$feature)
}

sets_ds   <- make_set_list_by_dataset(enrich, FEATURE_FOR_DATASET_VENN)
sets_feat <- make_set_list_by_feature(enrich, DATASET_FOR_FEATURE_VENN)

# remove empty sets
sets_ds   <- sets_ds[lengths(sets_ds) > 0]
sets_feat <- sets_feat[lengths(sets_feat) > 0]

if (!length(sets_ds)) stop("After filtering, no sets for the chosen feature.")
if (!length(sets_feat)) stop("After filtering, no sets for the chosen dataset.")

# ---------- save membership TSVs ----------
write_sets <- function(sets, out) {
  all <- utils::stack(sets)
  names(all) <- c("SNP","set")
  fwrite(all[order(all$set, all$SNP),], out)
}
write_sets(sets_ds,   file.path(SAVE_DIR, paste0("venn_sets_", FEATURE_FOR_DATASET_VENN, "_by_dataset.tsv")))
write_sets(sets_feat, file.path(SAVE_DIR, paste0("venn_sets_", DATASET_FOR_FEATURE_VENN, "_by_feature.tsv")))

# ---------- venn/upset ----------
have_pkg <- function(p) requireNamespace(p, quietly = TRUE)
plot_venn <- function(sets, title, outfile) {
  nset <- length(sets); if (nset < 2) { message("Need >=2 sets for: ", title); return() }
  if (have_pkg("ggVennDiagram") && nset <= 5) {
    p <- ggVennDiagram::ggVennDiagram(sets, label_alpha = 0, label = "count") +
      ggplot2::labs(title = title) + theme_minimal(base_size = 13)
    ggsave(outfile, p, width = 9, height = 7, dpi = 120)
  } else if (have_pkg("UpSetR")) {
    all_items <- sort(unique(unlist(sets)))
    mat <- sapply(sets, function(v) as.integer(all_items %in% v))
    df  <- as.data.frame(mat); rownames(df) <- all_items
    png(outfile, width = 1200, height = 800, res = 120)
    UpSetR::upset(df, sets = colnames(df), nsets = ncol(df), order.by = "freq",
                  mainbar.y.label = "Intersection size", sets.x.label = "Set size")
    mtext(title, side = 3, line = 1, font = 2, cex = 1.1)
    dev.off()
  } else {
    message("Install one plotting package:\n  install.packages('ggVennDiagram')\n  # or\n  install.packages('UpSetR')")
  }
}

plot_venn(sets_ds,
          paste0("SNPs enriched — by dataset (feature = ", FEATURE_FOR_DATASET_VENN, ")"),
          file.path(SAVE_DIR, paste0("venn_by_dataset_", FEATURE_FOR_DATASET_VENN, ".png")))
plot_venn(sets_feat,
          paste0("SNPs enriched — by feature (dataset = ", DATASET_FOR_FEATURE_VENN, ")"),
          file.path(SAVE_DIR, paste0("venn_by_feature_", DATASET_FOR_FEATURE_VENN, ".png")))

message("\nDone. Plots + TSVs saved in: ", SAVE_DIR)