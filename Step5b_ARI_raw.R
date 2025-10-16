# ===============================================================
# Merge explicit clusters.tsv by JoinID (k_drop0) + ARI w/ stability
# ===============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)  # adjustedRandIndex
})

# ---------- OUTPUT DIR ----------
OUTDIR <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------- EXPLICIT FULL PATHS ----------
PATHS <- list(
  A = list(
    genetic           = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_A/clusters.tsv",
    genetic_minclin   = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_minclin_A/clusters.tsv",
    genetic_fullclin  = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_fullclin_A/clusters.tsv",
    clinical_only     = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_clinical_only_A/clusters.tsv"
  ),
  B = list(
    genetic           = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_B/clusters.tsv",
    genetic_minclin   = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_minclin_B/clusters.tsv",
    genetic_fullclin  = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_fullclin_B/clusters.tsv",
    clinical_only     = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_clinical_only_B/clusters.tsv"
  ),
  C = list(
    genetic           = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_C/clusters.tsv",
    genetic_minclin   = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_minclin_C/clusters.tsv",
    genetic_fullclin  = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_genetic_fullclin_C/clusters.tsv",
    clinical_only     = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/step3_kmeans_clinical_only_C/clusters.tsv"
  )
)

# ---------- ID helpers ----------
canon_id <- function(x) gsub("[^A-Za-z0-9]", "", toupper(as.character(x)))
drop0    <- function(x) sub("^0+", "", x)  # k_drop0
make_join <- function(iid) drop0(canon_id(iid))

# ---------- readers ----------
read_clusters <- function(fp, feat) {
  if (is.null(fp) || !nzchar(fp) || !file.exists(fp)) {
    return(list(dt=NULL, info=data.table(feature=feat, path=fp, exists=FALSE, rows=0L,
                                         has_IID=FALSE, has_cluster=FALSE, has_stability=FALSE)))
  }
  dt <- tryCatch(fread(fp), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) {
    return(list(dt=NULL, info=data.table(feature=feat, path=fp, exists=TRUE, rows=0L,
                                         has_IID=FALSE, has_cluster=FALSE, has_stability=FALSE)))
  }
  nm <- names(dt)
  
  # ID column -> IID_raw; make JoinID = drop0(canon_id(IID_raw))
  idc <- c("IID","SubjectID","Cluster_subjectID","id","sample")
  idhit <- idc[idc %in% nm][1]
  has_IID <- !is.na(idhit)
  if (has_IID) setnames(dt, idhit, "IID_raw") else dt[, IID_raw := NA_character_]
  dt[, IID_raw := as.character(IID_raw)]
  dt[, JoinID  := make_join(IID_raw)]
  
  # cluster column -> cluster_<feat>
  clhit <- if ("cluster" %in% nm) "cluster" else {
    cand <- c("kmeans_cluster","pam_cluster","gmm_cluster","Cluster","clusters","label")
    hit <- cand[cand %in% nm][1]; ifelse(is.na(hit), NA_character_, hit)
  }
  has_cluster <- !is.na(clhit)
  clname <- paste0("cluster_", feat)
  if (has_cluster) setnames(dt, clhit, clname) else dt[, (clname) := NA_integer_]
  
  # stability -> stability_<feat>
  has_stability <- "stability" %in% nm
  stab_col <- paste0("stability_", feat)
  if (has_stability) setnames(dt, "stability", stab_col) else dt[, (stab_col) := NA_real_]
  
  keep <- c("JoinID","IID_raw", clname, stab_col, intersect("is_outlier", names(dt)))
  list(dt=unique(dt[, ..keep]),  # unique in case of duplicates
       info=data.table(feature=feat, path=fp, exists=TRUE, rows=nrow(dt),
                       has_IID=has_IID, has_cluster=has_cluster, has_stability=has_stability))
}

merge_feature_tables <- function(lst) {
  dts   <- lapply(lst, `[[`, "dt")
  infos <- rbindlist(lapply(lst, `[[`, "info"))
  dts <- Filter(Negate(is.null), dts)
  if (!length(dts)) return(list(merged=NULL, info=infos))
  merged <- Reduce(function(x,y) merge(x, y, by="JoinID", all=TRUE), dts)
  list(merged=merged, info=infos)
}

# ---------- ARI + plotting ----------
pairwise_ari <- function(v1, v2) {
  ok <- !is.na(v1) & !is.na(v2)
  if (!any(ok)) return(list(ARI=NA_real_, n=0L))
  list(ARI = adjustedRandIndex(as.integer(as.factor(v1[ok])), as.integer(as.factor(v2[ok]))),
       n = sum(ok))
}

plot_ari <- function(M, Ns, stabs_named, ds, out_png) {
  labs <- colnames(M)
  diag_labs <- labs
  for (i in seq_along(labs)) {
    mu <- stabs_named[labs[i]]
    if (!is.na(mu)) diag_labs[i] <- sprintf("%s\n(μ=%.02f)", labs[i], mu)
  }
  df <- as.data.table(as.table(M)); setnames(df, c("Var1","Var2","Freq"), c("x","y","ARI"))
  dn <- as.data.table(as.table(Ns)); setnames(dn, c("Var1","Var2","Freq"), c("x","y","n"))
  df <- merge(df, dn, by=c("x","y"), all=TRUE)
  
  df[, x := factor(x, levels=labs, labels=diag_labs)]
  df[, y := factor(y, levels=labs, labels=diag_labs)]
  
  p <- ggplot(df, aes(x, y, fill = ARI)) +
    geom_tile(color="white") +
    geom_text(aes(label = ifelse(is.na(ARI), "–", sprintf("%.02f\nn=%d", ARI, n))), size=4) +
    scale_fill_gradient(limits=c(0,1), low="white", high="firebrick", na.value="grey95") +
    coord_fixed() +
    theme_minimal(base_size=14) +
    theme(axis.text.x = element_text(angle=35, hjust=1)) +
    labs(title = sprintf("Adjusted Rand Index — dataset %s", ds),
         subtitle = "Diagonal shows mean bootstrap stability per feature (μ)",
         x=NULL, y=NULL, fill="ARI")
  ggsave(out_png, p, width=9.5, height=8, dpi=300)
}

# ---------- Runner ----------
run_ds <- function(DS, paths) {
  message(">>> Dataset ", DS)
  reads <- list(
    genetic          = read_clusters(paths$genetic,          "genetic"),
    genetic_minclin  = read_clusters(paths$genetic_minclin,  "genetic_minclin"),
    genetic_fullclin = read_clusters(paths$genetic_fullclin, "genetic_fullclin"),
    clinical_only    = read_clusters(paths$clinical_only,    "clinical_only")
  )
  info <- rbindlist(lapply(reads, `[[`, "info"))
  fwrite(info, file.path(OUTDIR, sprintf("presence_%s.tsv", DS)), sep="\t")
  
  merged_obj <- merge_feature_tables(reads)
  merged <- merged_obj$merged
  fwrite(merged, file.path(OUTDIR, sprintf("merged_assignments_%s.tsv", DS)), sep="\t", na = "")
  
  if (is.null(merged) || !nrow(merged)) return(invisible(NULL))
  
  # quick sanity — how many clinical rows do we have now?
  if ("cluster_clinical_only" %in% names(merged)) {
    cat(sprintf("[DS=%s] clinical_only non-NA rows: %d\n",
                DS, sum(!is.na(merged$cluster_clinical_only))))
  }
  
  # sizes per feature set
  feat_cols <- grep("^cluster_", names(merged), value=TRUE)
  if (length(feat_cols)) {
    sizes <- rbindlist(lapply(feat_cols, function(cn) {
      merged[!is.na(get(cn)), .N, by = get(cn)][order(get(cn))][
        , `:=`(feature = sub("^cluster_","",cn), cluster = get)]
    }), fill=TRUE)
    setnames(sizes, c("get","N","feature","cluster"), c("cluster","N","feature","cluster"))
    fwrite(sizes, file.path(OUTDIR, sprintf("sizes_%s.tsv", DS)), sep="\t")
  }
  
  # pairwise ARI (use JoinID overlap)
  feats <- sub("^cluster_", "", feat_cols)
  if (length(feat_cols) >= 2) {
    M  <- matrix(NA_real_, length(feat_cols), length(feat_cols), dimnames=list(feats,feats))
    Ns <- matrix(0L,       length(feat_cols), length(feat_cols), dimnames=list(feats,feats))
    for (i in seq_along(feat_cols)) for (j in seq_along(feat_cols)) {
      r <- pairwise_ari(merged[[feat_cols[i]]], merged[[feat_cols[j]]])
      M[i,j]  <- r$ARI; Ns[i,j] <- r$n
    }
    # mean stability per feature (stability_<feat>)
    stab_named <- setNames(rep(NA_real_, length(feats)), feats)
    for (f in feats) {
      sc <- paste0("stability_", f)
      if (sc %in% names(merged)) {
        v <- merged[[sc]]
        if (is.numeric(v)) stab_named[f] <- mean(v, na.rm = TRUE)
      }
    }
    fwrite(as.data.table(as.table(M)), file.path(OUTDIR, sprintf("ARI_pairwise_%s.tsv", DS)), sep="\t")
    plot_ari(M, Ns, stab_named, DS, file.path(OUTDIR, sprintf("ARI_heatmap_%s.png", DS)))
  }
  
  # cross-tabs clinical vs genetic
  cx <- function(lhs, rhs, tag) {
    if (!all(c(lhs, rhs) %in% names(merged))) return(invisible(NULL))
    x <- merged[!is.na(get(lhs)) & !is.na(get(rhs)), .N, by=.(get(lhs), get(rhs))]
    if (!nrow(x)) return(invisible(NULL))
    setnames(x, c("get","get.1","N"), c(lhs,rhs,"N"))
    x[, pct_of_lhs := N / sum(N), by = lhs]
    fwrite(x, file.path(OUTDIR, sprintf("crosstab_%s_%s.tsv", DS, tag)), sep="\t")
  }
  cx("cluster_clinical_only","cluster_genetic",          "clinical_vs_genetic")
  cx("cluster_clinical_only","cluster_genetic_minclin",  "clinical_vs_gen_minclin")
  cx("cluster_clinical_only","cluster_genetic_fullclin", "clinical_vs_gen_fullclin")
  
  invisible(TRUE)
}

for (ds in names(PATHS)) {
  tryCatch(run_ds(ds, PATHS[[ds]]),
           error = function(e) message("  ! ERROR DS=", ds, ": ", conditionMessage(e)))
}

cat(sprintf("\n✔ Done. Check %s for:\n- presence_*.tsv\n- merged_assignments_*.tsv (now includes JoinID)\n- sizes_*.tsv\n- ARI_pairwise_*.tsv + ARI_heatmap_*.png (diagonal shows mean stability μ)\n- crosstab_*.tsv\n", OUTDIR))

####ARI PLOT
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mclust)
  library(glue)
})

# ========== Absolute paths ==========
files <- list(
  A = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/merged_assignments_A.tsv",
  B = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/merged_assignments_B.tsv",
  C = "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/merged_assignments_C.tsv"
)

outdir <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

feature_order <- c("clinical_only", "genetic", "genetic_minclin", "genetic_fullclin")

pairwise_ari <- function(x, y) {
  ok <- complete.cases(x, y)
  if (!any(ok)) return(list(ari = NA_real_, n = 0))
  list(ari = adjustedRandIndex(x[ok], y[ok]), n = sum(ok))
}

plot_ari_heatmap <- function(tab, title, out_png) {
  p <- ggplot(tab, aes(f2, f1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 4) +
    scale_fill_gradient(low = "white", high = "#c23b2e", limits = c(0,1), na.value = "grey95", name = "ARI") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid = element_blank())
  ggsave(out_png, p, width = 7.5, height = 6, dpi = 220)
}

compute_ari <- function(dt, ds) {
  # Expected column names
  feats <- feature_order
  cols <- paste0("cluster_", feats)
  cols <- cols[cols %in% names(dt)]
  
  if (length(cols) < 2) {
    message("⚠️ Dataset ", ds, ": insufficient cluster columns found.")
    return(NULL)
  }
  
  res <- data.table()
  for (i in seq_along(cols)) {
    f1 <- feats[i]
    for (j in seq_along(cols)) {
      f2 <- feats[j]
      if (i < j) {
        a <- pairwise_ari(dt[[cols[i]]], dt[[cols[j]]])
        lab <- sprintf("%.2f\nn=%d", ifelse(is.na(a$ari), 0, a$ari), a$n)
        res <- rbind(res, data.table(f1 = f1, f2 = f2, value = a$ari, label = lab))
      } else if (i == j) {
        sc <- paste0("stability_", feats[i])
        lab <- if (sc %in% names(dt)) sprintf("stab\n%.2f", mean(dt[[sc]], na.rm = TRUE)) else "stab\nn/a"
        res <- rbind(res, data.table(f1 = f1, f2 = f2, value = NA_real_, label = lab))
      } else {
        res <- rbind(res, data.table(f1 = f1, f2 = f2, value = NA_real_, label = "–"))
      }
    }
  }
  
  res[, f1 := factor(f1, levels = feature_order)]
  res[, f2 := factor(f2, levels = feature_order)]
  
  out_tsv <- file.path(outdir, sprintf("B02500_ARI_from_merged_%s.tsv", ds))
  out_png <- file.path(outdir, sprintf("B02500_ARI_from_merged_%s.png", ds))
  fwrite(res, out_tsv, sep = "\t")
  plot_ari_heatmap(res, glue("Adjusted Rand Index — dataset {ds}"), out_png)
  message("✔ ", ds, ": wrote ", basename(out_tsv), " + ", basename(out_png))
  res[, dataset := ds]
  res
}

# ========== Run all ==========
all_tabs <- rbindlist(lapply(names(files), function(ds) {
  dt <- fread(files[[ds]])
  compute_ari(dt, ds)
}), fill = TRUE)

fwrite(all_tabs, file.path(outdir, "B02500_ARI_from_merged_ALL.tsv"), sep = "\t")
message("\n✅ Done. All ARI heatmaps and tables written to: ", outdir)
