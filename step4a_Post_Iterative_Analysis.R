# summarize_bootstrap_results.R
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# ======================= USER CONFIG ==========================
BASE <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"

BOOT_GRID <- c(1, 10, 100, 250, 500, 1000, 1500,2000, 2500, 3000)
DATASETS  <- c("A","B","C")
FEATURES  <- c("genetic","genetic_minclin","genetic_fullclin")

OUT_DIR   <- file.path(BASE, "bootstrap_grid_summary")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("BASE:   ", BASE)
message("OUTDIR: ", OUT_DIR)

# ======================= HELPERS ===============================
path_kmeans_dir <- function(B, feat, ds) {
  tag <- sprintf("grid_B%05d", B)
  file.path(BASE, tag, sprintf("step3_kmeans_%s_%s", feat, ds))
}

load_stab <- function(B, ds, feat) {
  d <- path_kmeans_dir(B, feat, ds)
  f <- file.path(d, "kmeans_stability_summary.tsv")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  dt$BOOT_N  <- B
  dt$dataset <- ds
  dt$feature <- feat
  dt
}

load_sil <- function(B, ds, feat) {
  d <- path_kmeans_dir(B, feat, ds)
  f <- file.path(d, "kmeans_silhouette_summary.tsv")
  if (!file.exists(f)) return(NULL)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)
  dt$BOOT_N  <- B
  dt$dataset <- ds
  dt$feature <- feat
  dt
}

# ======================= COLLECT ===============================
stab_bins <- list()
for (B in BOOT_GRID) {
  for (ds in DATASETS) {
    for (feat in FEATURES) {
      dt <- load_stab(B, ds, feat)
      if (!is.null(dt)) stab_bins[[length(stab_bins) + 1L]] <- as.data.table(dt)
    }
  }
}
stab_all <- data.table::rbindlist(stab_bins, fill = TRUE)
message("stab_all rows: ", if (nrow(stab_all)) nrow(stab_all) else 0)

sil_bins <- list()
for (B in BOOT_GRID) {
  for (ds in DATASETS) {
    for (feat in FEATURES) {
      dt <- load_sil(B, ds, feat)
      if (!is.null(dt)) sil_bins[[length(sil_bins) + 1L]] <- as.data.table(dt)
    }
  }
}
sil_all <- data.table::rbindlist(sil_bins, fill = TRUE)
message("sil_all rows: ", if (nrow(sil_all)) nrow(sil_all) else 0)

# Save combined tables (if any)
if (nrow(stab_all)) fwrite(stab_all, file.path(OUT_DIR, "all_kmeans_stability.tsv"), sep = "\t")
if (nrow(sil_all))  fwrite(sil_all,  file.path(OUT_DIR, "all_kmeans_silhouette.tsv"), sep = "\t")

# ======================= PLOTS ================================
plots <- list()

if (nrow(stab_all)) {
  # Keep only the "overall" rows for the BOOT_N trend
  if (!"level" %in% names(stab_all)) {
    message("WARN: 'level' column not found in stability files. Using all rows as overall.")
    stab_overall <- stab_all
  } else {
    stab_overall <- stab_all %>% filter(level == "overall")
  }
  
  # Expected columns: mean, sd (from your runner outputs)
  y_col <- if ("mean" %in% names(stab_overall)) "mean" else names(stab_overall)[which(names(stab_overall) %chin% c("mean_stability","avg_stability"))[1]]
  sd_col <- if ("sd" %in% names(stab_overall)) "sd" else names(stab_overall)[which(names(stab_overall) %chin% c("sd_stability","se_stability"))[1]]
  
  if (!is.na(y_col)) {
    p_stab <- stab_overall %>%
      group_by(BOOT_N, dataset, feature) %>%
      summarise(mean_stab = mean(.data[[y_col]], na.rm = TRUE),
                sd_stab   = if (!is.na(sd_col)) mean(.data[[sd_col]], na.rm = TRUE) else NA_real_,
                .groups = "drop") %>%
      ggplot(aes(x = BOOT_N, y = mean_stab, color = feature)) +
      geom_line(linewidth = 1.05) +
      geom_point(size = 2) +
      {if (!is.na(sd_col)) geom_errorbar(aes(ymin = mean_stab - sd_stab, ymax = mean_stab + sd_stab), width = 50)} +
      facet_wrap(~ dataset, scales = "free_y") +
      labs(title = "K-means Bootstrap Stability vs BOOT_N",
           x = "BOOT_N (resamples)", y = "Mean Stability (± SD)", color = "Feature set") +
      theme_minimal(base_size = 13)
    print(p_stab)
    ggsave(file.path(OUT_DIR, "stability_vs_boot_by_dataset.png"), p_stab, width = 10, height = 6, dpi = 150)
    plots$stability <- p_stab
  } else {
    message("Could not find a stability mean column in stability files.")
  }
} else {
  message("No stability data to plot.")
}

if (nrow(sil_all)) {
  # Expected columns: K, mean_silhouette (from your runner)
  ysil <- if ("mean_silhouette" %in% names(sil_all)) "mean_silhouette" else if ("silhouette" %in% names(sil_all)) "silhouette" else NA_character_
  
  if (!is.na(ysil) && "K" %in% names(sil_all)) {
    p_sil <- sil_all %>%
      ggplot(aes(x = K, y = .data[[ysil]], color = feature, group = interaction(feature, BOOT_N))) +
      geom_line(linewidth = 1.05) +
      geom_point(size = 2) +
      facet_grid(BOOT_N ~ dataset, labeller = label_both) +
      labs(title = "Average Silhouette vs K across BOOT_N",
           x = "K", y = "Average Silhouette", color = "Feature set") +
      theme_minimal(base_size = 13)
    print(p_sil)
    ggsave(file.path(OUT_DIR, "silhouette_vs_k_by_boot.png"), p_sil, width = 12, height = 8, dpi = 150)
    plots$silhouette <- p_sil
  } else {
    message("Silhouette files missing 'K' or a silhouette mean column.")
  }
} else {
  message("No silhouette data to plot.")
}

message("Done. Outputs (if any) saved to: ", OUT_DIR)