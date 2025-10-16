## ==========================================================
## Top-SNPs → (optional) nearest-gene + KEGG + overlap visuals
## Works even if Bioconductor packages aren't installed:
## - Always: HLA plot, overlap counts, presence plot
## - Optional: KEGG (needs clusterProfiler + org.Hs.eg.db)
## - Optional: nearest gene (needs biomaRt)
## ==========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glue)
})

## ---------- EDIT ONLY THIS ----------
GRID_ID <- "grid_B02500"
VIZDIR  <- file.path("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min",
                     GRID_ID, "viz_summary")
## ------------------------------------

dir.create(VIZDIR, showWarnings = FALSE, recursive = TRUE)

## ---------- Soft deps (don’t fail if missing) ----------
have_pkg <- function(x) { requireNamespace(x, quietly = TRUE) }

have_biomaRt         <- have_pkg("biomaRt")
have_org_db          <- have_pkg("org.Hs.eg.db")
have_clusterProfiler <- have_pkg("clusterProfiler")
have_UpSetR          <- have_pkg("UpSetR")
have_ggraph          <- have_pkg("ggraph") && have_pkg("igraph")

if (!have_biomaRt)         message("Note: biomaRt not available → skipping nearest-gene lookup.")
if (!have_org_db)          message("Note: org.Hs.eg.db not available → KEGG will be skipped.")
if (!have_clusterProfiler) message("Note: clusterProfiler not available → KEGG will be skipped.")
if (!have_UpSetR)          message("Note: UpSetR not available → using overlap counts plot instead of UpSet.")
if (!have_ggraph)          message("Note: ggraph/igraph not available → skipping gene–pathway network.")

## ---------- Input files (prefer cleaned; else raw) ----------
p_clean <- list(
  A = file.path(VIZDIR, "top_snps_A.cleaned.tsv"),
  B = file.path(VIZDIR, "top_snps_B.cleaned.tsv"),
  C = file.path(VIZDIR, "top_snps_C.cleaned.tsv")
)
p_raw <- list(
  A = file.path(VIZDIR, "top_snps_driving_clusters_A.tsv"),
  B = file.path(VIZDIR, "top_snps_driving_clusters_B.tsv"),
  C = file.path(VIZDIR, "top_snps_driving_clusters_C.tsv")
)

rd <- function(fp, tag) {
  if (!file.exists(fp)) stop("[", tag, "] missing: ", fp)
  dt <- fread(fp)
  if (!nrow(dt)) stop("[", tag, "] empty: ", fp)
  setDT(dt)
}

standardize <- function(dt, dataset_tag) {
  nms <- names(dt)
  idc <- c("SNP","snp","snp_id","id","rsid","RSID","variant","Variant")
  hit <- idc[idc %in% nms]
  if (length(hit)) setnames(dt, hit[1], "SNP") else {
    # fallback: assume first col is SNP
    setnames(dt, 1L, "SNP")
  }
  ldc <- c("loading","Loading","PC1_loading","pc1_loading","abs_loading","weight")
  hitL <- ldc[ldc %in% nms]
  if (length(hitL)) setnames(dt, hitL[1], "loading") else dt[, loading := NA_real_]
  suppressWarnings(dt[, loading := as.numeric(loading)])
  dt <- dt[!is.na(SNP) & SNP != ""]
  dt[, dataset := dataset_tag]
  unique(dt[, .(SNP, loading, dataset)])
}

load_one <- function(tag) {
  fp <- if (file.exists(p_clean[[tag]])) p_clean[[tag]] else p_raw[[tag]]
  message("[", tag, "] using: ", fp)
  standardize(rd(fp, tag), tag)
}

A <- load_one("A")
B <- load_one("B")
C <- load_one("C")

## Add proxy loading if a table has no numeric loadings
add_proxy_loading <- function(dt) {
  if (all(is.na(dt$loading))) dt[, loading := rank(-seq_len(.N))/.N]
  dt
}
A <- add_proxy_loading(A); B <- add_proxy_loading(B); C <- add_proxy_loading(C)

## ---------- Combine & presence flags ----------
top_all <- rbindlist(list(A,B,C), use.names = TRUE, fill = TRUE)
fwrite(top_all, file.path(VIZDIR, "top_snps_ABC_combined.tsv"), sep = "\t")

wide_presence <- dcast(unique(top_all[, .(SNP, dataset)]),
                       SNP ~ dataset, value.var = "dataset", fun.aggregate = length)
for (d in c("A","B","C")) if (!d %in% names(wide_presence)) wide_presence[, (d) := 0L]
wide_presence[, shared_n := A + B + C]

top_all <- merge(top_all, wide_presence[, .(SNP, shared_n)], by = "SNP", all.x = TRUE)

## ---------- Parse CHR:POS (for HLA figure) ----------
parse_chr_pos <- function(x) {
  y <- gsub("^chr","", x, ignore.case = TRUE)
  ok <- grepl("^\\d+[:_]\\d+$", y)
  chr <- pos <- rep(NA_character_, length(x))
  chr[ok] <- sub("^([0-9]+)[:_](\\d+)$", "\\1", y[ok])
  pos[ok] <- sub("^([0-9]+)[:_](\\d+)$", "\\2", y[ok])
  data.table(chr = chr, pos = suppressWarnings(as.integer(pos)), is_chrpos = ok)
}
top_all <- cbind(top_all, parse_chr_pos(top_all$SNP))

## ---------- Optional: nearest gene via biomaRt ----------
if (have_biomaRt) {
  message("Annotating CHR:POS SNPs to nearest gene (±50kb)…")
  W <- 5e4
  mart <- tryCatch({
    biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  }, error = function(e) NULL)
  annotate_chrpos <- function(chr, pos) {
    if (is.null(mart) || is.na(chr) || is.na(pos)) {
      return(data.table(nearest_gene = NA_character_, dist_bp = NA_integer_))
    }
    start <- max(1L, pos - W); end <- pos + W
    g <- tryCatch(biomaRt::getBM(
      attributes = c("hgnc_symbol","start_position","end_position"),
      filters = c("chromosome_name","start","end"),
      values = list(chr, start, end), mart = mart
    ), error = function(e) NULL)
    if (is.null(g) || !nrow(g)) return(data.table(nearest_gene = NA_character_, dist_bp = NA_integer_))
    g <- as.data.table(g)
    g[, mid := (start_position + end_position)/2]
    g[, dist := abs(mid - pos)]
    g <- g[order(dist)][1]
    data.table(nearest_gene = g$hgnc_symbol, dist_bp = as.integer(g$dist))
  }
  ann <- rbindlist(lapply(seq_len(nrow(top_all)), function(i) {
    if (isTRUE(top_all$is_chrpos[i])) annotate_chrpos(top_all$chr[i], top_all$pos[i])
    else data.table(nearest_gene = NA_character_, dist_bp = NA_integer_)
  }))
  top_all <- cbind(top_all, ann)
} else {
  top_all[, `:=`(nearest_gene = NA_character_, dist_bp = NA_integer_)]
}

## ---------- Save annotated table ----------
fwrite(top_all[order(-abs(loading), SNP)], file.path(VIZDIR, "top_snps_annotated.tsv"), sep = "\t")

## ---------- Figure 1: HLA “Manhattan” (chr6 only) ----------
hla <- top_all[(chr == "6" | chr == 6) & !is.na(pos)]
if (nrow(hla)) {
  hla[, Mb := pos/1e6]
  hla[, shared := factor(shared_n >= 2, levels = c(FALSE, TRUE), labels = c("FALSE","TRUE"))]
  p_hla <- ggplot(hla, aes(Mb, abs(loading), color = shared)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(size = 2.5, alpha = 0.95) +
    ggrepel::geom_text_repel(
      aes(label = ifelse(shared == "TRUE", SNP, "")),
      size = 3, min.segment.length = 0.05, box.padding = 0.2, max.overlaps = 50
    ) +
    scale_color_manual(values = c("FALSE" = "#9e9e9e", "TRUE" = "#d7301f")) +
    labs(title = "Shared SNPs among top loadings (Chr6 / HLA region)",
         subtitle = "Red points appear in ≥2 datasets (A, B, C)",
         x = "Chr6 position (Mb)", y = "|loading|", color = "Shared A/B/C") +
    theme_minimal(base_size = 12)
  ggsave(file.path(VIZDIR, "shared_HLA_SNPs_manhattan.pdf"), p_hla, width = 9.5, height = 5.8, dpi = 300)
  ggsave(file.path(VIZDIR, "shared_HLA_SNPs_manhattan.png"), p_hla, width = 9.5, height = 5.8, dpi = 300)
}

## ---------- Figure 2: Overlap (counts by “#datasets per SNP”) ----------
presence <- rbindlist(list(
  data.table(SNP = unique(A$SNP), A = 1L, B = 0L, C = 0L),
  data.table(SNP = unique(B$SNP), A = 0L, B = 1L, C = 0L),
  data.table(SNP = unique(C$SNP), A = 0L, B = 0L, C = 1L)
), fill = TRUE)
presence[is.na(A), A := 0L][is.na(B), B := 0L][is.na(C), C := 0L]
presence[, sets := A + B + C]
pres_counts <- presence[, .N, by = sets][order(sets)]
fwrite(pres_counts, file.path(VIZDIR, "topSNP_presence_counts.tsv"), sep = "\t")

p_pres <- ggplot(pres_counts, aes(factor(sets), N)) +
  geom_col(width = 0.7, fill = "#6A5ACD") +
  geom_text(aes(label = N), vjust = -0.25, size = 4) +
  labs(title = "How many datasets each SNP appears in (A/B/C)",
       x = "# of datasets containing the SNP", y = "SNP count") +
  theme_minimal(base_size = 12)
ggsave(file.path(VIZDIR, "topSNP_presence_counts.pdf"), p_pres, width = 7.2, height = 4.6, dpi = 300)
ggsave(file.path(VIZDIR, "topSNP_presence_counts.png"), p_pres, width = 7.2, height = 4.6, dpi = 300)

## ---------- Figure 3 (robust): Intersections across A/B/C without UpSetR ----------
{
  # Build membership table (unique SNP universe)
  mem <- Reduce(function(x, y) merge(x, y, by = "SNP", all = TRUE),
                list(A[, .(SNP)], B[, .(SNP)], C[, .(SNP)]))
  # Boolean presence flags
  mem[, A := SNP %in% A$SNP]
  mem[, B := SNP %in% B$SNP]
  mem[, C := SNP %in% C$SNP]
  
  # Label each SNP by the exact combo of datasets it belongs to
  labs <- c("A","B","C")
  mem[, combo := vapply(seq_len(.N), function(i) {
    hits <- labs[c(A[i], B[i], C[i])]
    if (length(hits) == 0) "(none)" else paste(hits, collapse = "∩")
  }, character(1))]
  
  # Count per combination and order nicely
  levs <- c("A","B","C","A∩B","A∩C","B∩C","A∩B∩C")
  tab  <- mem[, .N, by = combo]
  # keep only combos that actually appear
  levs_present <- levs[levs %in% tab$combo]
  tab[, combo := factor(combo, levels = c(levs_present, setdiff(levels(factor(combo)), levs_present)))]
  setorder(tab, combo)
  
  # Save counts table
  fwrite(tab, file.path(VIZDIR, "topSNP_intersection_counts.tsv"), sep = "\t")
  
  # Plot
  p_int <- ggplot(tab[combo != "(none)"], aes(x = combo, y = N)) +
    geom_col(width = 0.7, fill = "#6A5ACD") +
    geom_text(aes(label = N), vjust = -0.25, size = 4) +
    labs(title = "SNP intersections across datasets",
         subtitle = "Counts of SNPs appearing in A, B, C and their overlaps",
         x = "Dataset combination", y = "SNP count") +
    theme_minimal(base_size = 12)
  
  ggsave(file.path(VIZDIR, "topSNP_intersections.png"), p_int, width = 8.5, height = 5.0, dpi = 300)
  ggsave(file.path(VIZDIR, "topSNP_intersections.pdf"), p_int, width = 8.5, height = 5.0, dpi = 300)
}

## ---------- Optional KEGG enrichment ----------
if (have_org_db && have_clusterProfiler) {
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(clusterProfiler)
  })
  # If nearest_gene is available (biomaRt), use it; else skip KEGG neatly
  genes <- unique(na.omit(top_all$nearest_gene))
  if (length(genes) >= 5) {
    sym2ent <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys = genes, keytype = "SYMBOL",
                                     columns = "ENTREZID")
    entrez <- unique(na.omit(sym2ent$ENTREZID))
    if (length(entrez) >= 5) {
      ek <- tryCatch(enrichKEGG(gene = entrez, organism = "hsa", pvalueCutoff = 0.1),
                     error = function(e) NULL)
      if (!is.null(ek)) {
        ekdf <- as.data.frame(ek)
        fwrite(as.data.table(ekdf), file.path(VIZDIR, "kegg_enrichment.tsv"), sep = "\t")
        if (nrow(ekdf)) {
          ekdf$Description <- factor(ekdf$Description, levels = rev(ekdf$Description))
          p_kegg <- ggplot(head(ekdf, 15), aes(Description, -log10(p.adjust))) +
            geom_col(fill = "#2b8cbe") +
            coord_flip() +
            labs(title = "KEGG enrichment (nearest genes to top SNPs)", x = "", y = "-log10(adj. p)") +
            theme_minimal(base_size = 12)
          ggsave(file.path(VIZDIR, "kegg_enrichment_top15.pdf"), p_kegg, width = 8, height = 6, dpi = 300)
          ggsave(file.path(VIZDIR, "kegg_enrichment_top15.png"), p_kegg, width = 8, height = 6, dpi = 300)
        }
      }
    } else {
      message("KEGG: <5 ENTREZ IDs → skip.")
    }
  } else {
    message("KEGG: no nearest_gene values → skip (needs biomaRt or pre-annotated genes).")
  }
}

cat(glue("
✔ Outputs → {VIZDIR}
  - top_snps_ABC_combined.tsv
  - top_snps_annotated.tsv
  - shared_HLA_SNPs_manhattan.(png|pdf)  (if CHR:POS exist on chr6)
  - topSNP_presence_counts.(png|pdf)
  - topSNP_upset.(png|pdf)                (if UpSetR installed)
  - kegg_enrichment.(tsv), kegg_enrichment_top15.(png|pdf) (if gene mapping + KEGG available)
"))