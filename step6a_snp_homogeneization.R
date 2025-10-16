###############################################
## Top SNPs (A/B/C) — robust read + standardize
## Outputs: cleaned A/B/C TSVs, combined TSV,
##          overlap summary, and a presence plot
###############################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ----------- PATHS (edit BASE only if needed) -----------
BASE <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary"

PATHS <- list(
  A = file.path(BASE, "top_snps_driving_clusters_A.tsv"),
  B = file.path(BASE, "top_snps_driving_clusters_B.tsv"),
  C = file.path(BASE, "top_snps_driving_clusters_C.tsv")
)

## ----------- Helpers -----------
safe_read <- function(fp, tag) {
  if (!file.exists(fp)) stop("[", tag, "] Missing file: ", fp)
  dt <- tryCatch(fread(fp), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) stop("[", tag, "] Empty or unreadable: ", fp)
  message("[", tag, "] rows=", nrow(dt), " | cols: ", paste(names(dt), collapse = ", "))
  setDT(dt)
}

# Standardize to {SNP, loading, dataset}
standardize <- function(x, dataset_tag) {
  dt <- as.data.table(x)
  nms <- names(dt)
  
  # 1) Find/rename SNP column
  idc <- c("SNP","snp","snp_id","id","rsid","RSID","variant","Variant")
  hit <- idc[idc %in% nms]
  if (length(hit)) {
    setnames(dt, hit[1], "SNP")
  } else {
    if (ncol(dt) == 1L) {
      setnames(dt, 1L, "SNP")
    } else {
      stop("[", dataset_tag, "] Could not locate an SNP column. Columns: ",
           paste(nms, collapse=", "))
    }
  }
  
  # 2) Ensure numeric 'loading' column (if present)
  ldc <- c("loading","Loading","PC1_loading","pc1_loading","abs_loading","weight")
  hitL <- ldc[ldc %in% nms]
  if (length(hitL)) {
    setnames(dt, hitL[1], "loading")
    suppressWarnings(dt[, loading := as.numeric(loading)])
  } else {
    dt[, loading := NA_real_]
  }
  
  # 3) Clean + label
  dt <- dt[!is.na(SNP) & SNP != ""]
  dt[, dataset := dataset_tag]
  dt <- dt[, .(SNP, loading, dataset)]
  unique(dt, by = c("dataset","SNP"))
}

## ----------- Read & standardize A/B/C -----------
rawA <- safe_read(PATHS$A, "A")
rawB <- safe_read(PATHS$B, "B")
rawC <- safe_read(PATHS$C, "C")

A <- standardize(rawA, "A")
B <- standardize(rawB, "B")
C <- standardize(rawC, "C")

## If a file lacks numeric loadings entirely, you can create a proxy rank (optional):
add_proxy_loading <- function(dt) {
  if (all(is.na(dt$loading))) {
    dt[, loading := rank(-seq_len(.N))/.N]  # simple decreasing proxy
  }
  dt
}
A <- add_proxy_loading(A)
B <- add_proxy_loading(B)
C <- add_proxy_loading(C)

## ----------- Save cleaned per-dataset tables -----------
fwrite(A, file.path(BASE, "top_snps_A.cleaned.tsv"), sep = "\t")
fwrite(B, file.path(BASE, "top_snps_B.cleaned.tsv"), sep = "\t")
fwrite(C, file.path(BASE, "top_snps_C.cleaned.tsv"), sep = "\t")

## ----------- Combine & save -----------
top_all <- rbindlist(list(A,B,C), use.names = TRUE, fill = TRUE)
fwrite(top_all, file.path(BASE, "top_snps_ABC_combined.tsv"), sep = "\t")

message("Combined table: ", nrow(top_all), " rows; cols: ", paste(names(top_all), collapse = ", "))

## ----------- Overlap summary (unique/shared across sets) -----------
setA <- unique(A$SNP); setB <- unique(B$SNP); setC <- unique(C$SNP)
labs <- c("A","B","C")

overlap_dt <- data.table(
  region = c("A only","B only","C only",
             "A∩B only","A∩C only","B∩C only",
             "A∩B∩C"),
  count  = c(
    length(setdiff(setA, union(setB,setC))),
    length(setdiff(setB, union(setA,setC))),
    length(setdiff(setC, union(setA,setB))),
    length(setdiff(intersect(setA,setB), setC)),
    length(setdiff(intersect(setA,setC), setB)),
    length(setdiff(intersect(setB,setC), setA)),
    length(Reduce(intersect, list(setA,setB,setC)))
  )
)

fwrite(overlap_dt, file.path(BASE, "top_snps_overlap_summary.tsv"), sep = "\t")

## ----------- Presence/overlap bar plot (quick visual) -----------
presence <- rbindlist(list(
  data.table(SNP=setA, A=1L, B=0L, C=0L),
  data.table(SNP=setB, A=0L, B=1L, C=0L),
  data.table(SNP=setC, A=0L, B=0L, C=1L)
), fill = TRUE)
presence[is.na(A), A:=0L][is.na(B), B:=0L][is.na(C), C:=0L]

presence[, sets := A + B + C]
pres_counts <- presence[, .N, by = sets][order(sets)]
fwrite(pres_counts, file.path(BASE, "top_snps_presence_counts.tsv"), sep = "\t")

p <- ggplot(pres_counts, aes(factor(sets), N)) +
  geom_col(width = 0.7, fill = "#6A5ACD") +
  geom_text(aes(label = N), vjust = -0.3, size = 4) +
  labs(title = "How many datasets each SNP appears in",
       x = "Number of datasets containing the SNP (A/B/C)", y = "Count of SNPs") +
  theme_minimal(base_size = 13)
ggsave(file.path(BASE, "top_snps_presence_counts.png"), p, width = 7, height = 4.5, dpi = 300)

cat("\n✔ Done. Wrote:\n",
    " - ", file.path(BASE, "top_snps_A.cleaned.tsv"), "\n",
    " - ", file.path(BASE, "top_snps_B.cleaned.tsv"), "\n",
    " - ", file.path(BASE, "top_snps_C.cleaned.tsv"), "\n",
    " - ", file.path(BASE, "top_snps_ABC_combined.tsv"), "\n",
    " - ", file.path(BASE, "top_snps_overlap_summary.tsv"), "\n",
    " - ", file.path(BASE, "top_snps_presence_counts.tsv"), "\n",
    " - ", file.path(BASE, "top_snps_presence_counts.png"), "\n\n")