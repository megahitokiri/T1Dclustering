# ======================================================
# SNP×SNP interaction screens by dataset (A/B/C)
# - Uses *_homog_numeric.rds (rows = samples)
# - Uses viz_summary/merged_assignments_*.tsv for clusters
# - Outcome: in most-prevalent cluster (1) vs rest (0)
# - Outputs: PNG/PDF (300 dpi) + TSVs in viz_summary/
# ======================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(igraph)
  library(glue)
})

# --------- Paths ----------
ROOT   <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs"
GRID   <- "grid_B02500"
BUCKET <- file.path(ROOT, "cluster_runs_min", GRID)
VIZ    <- file.path(BUCKET, "viz_summary")
dir.create(VIZ, recursive = TRUE, showWarnings = FALSE)

# --------- ID canon (same as your pipeline) ----------
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
drop0    <- function(x) sub("^0+", "", x)
mk_join  <- function(x) drop0(canon_id(x))

# --------- SNP list (edit if you want) ----------
rsid2coord <- c(
  rs6679677="01:113761186", rs3024505="01:206766559",
  rs2111485="02:162254026", rs3087243="02:203874196",
  rs17388568="04:122408207",
  rs9500974="06:029760476", rs1233320="06:029840255", rs72848653="06:029868825", rs9259013="06:029874690",
  rs12189871="06:031284147", rs9266268="06:031358273", rs16899379="06:031375490", rs149663102="06:031376405",
  rs9268500="06:032408740", rs17840116="06:032415221", rs75658393="06:032427740",
  rs9271347="06:032615766", rs1281935="06:032616043", rs9271351="06:032616083", rs1281943="06:032630313",
  rs9405117="06:032634974", rs9469200="06:032635435", rs9273032="06:032644083", rs9273076="06:032644524",
  rs2395228="06:032655446", rs117806464="06:032658670", rs9273369="06:032658707", rs1049124="06:032659936",
  rs62406889="06:032704437", rs9275490="06:032705608", rs10947332="06:032709663",
  rs6934289="06:033077179", rs17214657="06:033079396", rs2567287="06:033081408", rs9378176="06:033081532", rs3129197="06:033103130",
  rs72928038="06:090267049", rs9388489="06:126377573", rs1738074="06:159044945",
  rs4948088="07:050959497",
  rs3842753="11:002159830"
)

SNPS <- c(
  "06:033081408",  # rs2567287
  "06:032659936",  # rs1049124
  "07:050959497",  # rs4948088
  "02:162254026",  # rs2111485
  "06:031358273",  # rs9266268
  "06:159044945",  # rs1738074
  "06:033077179",  # rs6934289
  "11:002159830",  # rs3842753
  "06:032658670",  # rs117806464
  "06:032658707",  # rs9273369
  "06:032655446",  # rs2395228
  "06:032644083",  # rs9273032
  "06:032644524",  # rs9273076
  "06:032616083",  # rs9271351
  "06:032615766",  # rs9271347
  "06:032616043",  # rs1281935
  "06:032630313",  # rs1281943
  "06:032634974",  # rs9405117
  "06:032635435",  # rs9469200
  "06:032704437"   # rs62406889
)


# --------- Robust genotype loader (IDs + selected SNPs) ----------
load_geno_matrix_ds <- function(ds, want_snps = SNPS) {
  rds_fp <- file.path(ROOT, paste0(ds, "_homog_numeric.rds"))
  if (!file.exists(rds_fp)) stop("No RDS for ", ds, ": ", rds_fp)
  X <- readRDS(rds_fp)
  
  # Try to get a numeric matrix with rownames + colnames
  get_mat <- function(obj) {
    # plain matrix / data.frame
    if (is.matrix(obj) || is.data.frame(obj)) return(as.matrix(obj))
    # sparse matrix with dimnames
    dn <- try(attr(obj, "Dimnames"), silent = TRUE)
    if (!inherits(dn, "try-error") && !is.null(dn)) {
      m <- try(as.matrix(obj), silent = TRUE)
      if (!inherits(m, "try-error")) return(m)
    }
    # list formats: common slots
    if (is.list(obj)) {
      # (1) matrix-like under common keys
      for (k in c("X","geno","G","mat","data")) {
        if (!is.null(obj[[k]])) {
          mat <- try(as.matrix(obj[[k]]), silent = TRUE)
          if (!inherits(mat, "try-error")) return(mat)
        }
      }
    }
    stop("Unsupported RDS structure for ", ds, " — could not extract a matrix.")
  }
  
  M <- get_mat(X)
  rn <- rownames(M); cn <- colnames(M)
  if (is.null(rn) || is.null(cn)) stop("Genotype matrix lacks row/col names for ", ds)
  
  # Intersect SNPs
  keep <- intersect(want_snps, cn)
  if (!length(keep)) stop("None of the requested SNPs were found in ", ds, ".")
  M2 <- M[, keep, drop = FALSE]
  
  dt <- as.data.table(M2)
  dt[, IID_raw := rn]
  dt[, JoinID := mk_join(IID_raw)]
  setcolorder(dt, c("JoinID","IID_raw", keep))
  dt[]
}

# --------- Load cluster labels (prefer fullclin) ----------
load_cluster_ds <- function(ds,
                            prefer_order = c("cluster_genetic_fullclin","cluster_genetic",
                                             "cluster_genetic_minclin","cluster_clinical_only")) {
  cand <- c(
    file.path(VIZ, paste0("merged_assignments_", ds, ".tsv")),
    file.path(VIZ, paste0(GRID, "_merged_assignments_", ds, ".tsv")),
    file.path(VIZ, paste0(ds, "_merged_assignments.tsv"))
  )
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("No merged_assignments found for ", ds)
  
  clu <- fread(hit)
  if ("JoinID" %in% names(clu)) {
    clu[, JoinID := mk_join(JoinID)]
  } else {
    id_cand <- c("IID","IID_raw","IID_raw.x","IID_raw.y","SubjectID","Cluster_subjectID")
    id_cand <- id_cand[id_cand %in% names(clu)]
    if (!length(id_cand)) stop("No ID-like column in ", hit)
    setnames(clu, id_cand[1], "IID_tmp", skip_absent = TRUE)
    clu[, JoinID := mk_join(IID_tmp)]
  }
  
  cl_cols <- intersect(prefer_order, names(clu))
  if (!length(cl_cols)) stop("No cluster_* column in ", hit)
  chosen <- cl_cols[1]
  out <- clu[, .(JoinID, cluster = get(chosen))]
  out <- out[!is.na(cluster)]
  attr(out, "cluster_col") <- chosen
  attr(out, "file") <- hit
  out[]
}

# --------- Interaction runner (logistic: in target cluster vs rest) ----------
run_interactions <- function(geno_dt, clusters_dt, alpha = 0.05) {
  # align
  m <- merge(geno_dt, clusters_dt, by = "JoinID")
  # choose target cluster = most frequent
  tgt <- m[, .N, by = cluster][order(-N)][1, cluster]
  Y   <- as.integer(m$cluster == tgt)
  
  snps <- setdiff(names(geno_dt), c("JoinID","IID_raw"))
  G    <- as.matrix(m[, ..snps])
  # per-pair p-values for interaction term in logistic regression
  pr <- t(combn(snps, 2))
  p  <- rep(NA_real_, nrow(pr))
  
  # small helper to fit safely
  get_p <- function(a, b) {
    d <- data.frame(Y = Y, A = G[, a], B = G[, b])
    d <- d[complete.cases(d), ]
    if (nrow(d) < 50) return(NA_real_)
    # guard against constant cols
    if (sd(d$A) == 0 || sd(d$B) == 0) return(NA_real_)
    fit <- try(glm(Y ~ A + B + A:B, data = d, family = binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    sm <- summary(fit)$coefficients
    ip <- try(sm["A:B","Pr(>|z|)"], silent = TRUE)
    if (inherits(ip, "try-error")) return(NA_real_) else return(ip)
  }
  
  for (k in seq_len(nrow(pr))) p[k] <- get_p(pr[k,1], pr[k,2])
  padj <- p.adjust(p, method = "BH")
  edges <- data.table(V1 = pr[,1], V2 = pr[,2], p = p, padj = padj)
  list(edges = edges, target = tgt, n = nrow(m))
}

# --------- Plots ----------
plot_tile <- function(edges, alpha = 0.05, title = "") {
  snps <- sort(unique(c(edges$V1, edges$V2)))
  M <- matrix(NA_real_, length(snps), length(snps), dimnames = list(snps, snps))
  for (i in seq_len(nrow(edges))) {
    a <- edges$V1[i]; b <- edges$V2[i]; M[a,b] <- M[b,a] <- edges$p[i]
  }
  DT <- as.data.table(as.table(M))
  setnames(DT, c("Var1","Var2","p"))
  DT <- DT[Var1 != Var2]
  DT <- DT[as.integer(factor(Var1, levels = snps)) <
             as.integer(factor(Var2, levels = snps))]
  
  DT[, sig := fifelse(is.na(p), "NA",
                      fifelse(p < alpha, "<alpha", ">=alpha"))]
  
  ggplot(DT, aes(Var2, Var1, fill = sig)) +
    geom_tile(color = "grey90", linewidth = 0.3) +
    scale_fill_manual(
      values = c("<alpha" = "black", ">=alpha" = "white", "NA" = "grey85"),
      labels = c("<alpha" = paste0("<", alpha),
                 ">=alpha" = paste0("≥", alpha),
                 "NA" = "NA")
    ) +
    coord_fixed() +
    labs(
      x = NULL, y = NULL, fill = NULL, title = title,
      subtitle = paste0("Tile = p-value of interaction; α = ", alpha)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

plot_network <- function(edges, alpha = 0.05, use_fdr = FALSE, title = "") {
  thr <- if (use_fdr) "padj" else "p"
  E <- edges[!is.na(get(thr)) & get(thr) < alpha]
  if (!nrow(E)) return(NULL)
  g <- graph_from_data_frame(E[, .(from = V1, to = V2, p, padj)], directed = FALSE)
  lay <- layout_with_fr(g)
  # edge width by -log10(p)
  ew <- scales::rescale(-log10(E(g)$p), to = c(0.5, 4))
  plot(g, layout = lay,
       vertex.size = 18, vertex.label.cex = 0.8,
       edge.width = ew, main = title)
  invisible(g)
}

# --------- Driver for A/B/C ----------
analyze_ds <- function(ds, alpha = 0.05) {
  message("\n>>> Dataset ", ds)
  geno  <- load_geno_matrix_ds(ds, SNPS)
  clus  <- load_cluster_ds(ds)
  res   <- run_interactions(geno, clus, alpha = alpha)
  
  message(glue("  Target cluster = {res$target} (n={res$n})"))
  n_sig <- res$edges[!is.na(p) & p < alpha, .N]
  message(glue("  Significant pairs at p<{alpha}: {n_sig}"))
  
  # Save TSV
  TSV <- file.path(VIZ, glue("{ds}_snp_interactions.tsv"))
  fwrite(res$edges, TSV, sep = "\t")
  
  # Tile heatmap
  P1 <- plot_tile(res$edges, alpha = alpha,
                  title = glue("{ds}: SNP×SNP interactions (logistic)"))
  ggsave(file.path(VIZ, glue("{ds}_snp_interactions_tile.png")),
         P1, width = 8.5, height = 7, dpi = 300)
  ggsave(file.path(VIZ, glue("{ds}_snp_interactions_tile.pdf")),
         P1, width = 8.5, height = 7, dpi = 300)
  
  # Network (raw p)
  png(file.path(VIZ, glue("{ds}_snp_interactions_network_p.png")),
      width = 2200, height = 1700, res = 300)
  plot_network(res$edges, alpha = alpha, use_fdr = FALSE,
               title = glue("{ds}: interaction network (p<{alpha})"))
  dev.off()
  
  # Network (FDR)
  png(file.path(VIZ, glue("{ds}_snp_interactions_network_fdr.png")),
      width = 2200, height = 1700, res = 300)
  plot_network(res$edges, alpha = alpha, use_fdr = TRUE,
               title = glue("{ds}: interaction network (FDR<{alpha})"))
  dev.off()
  
  message("  ✔ Wrote: ", basename(TSV), " + PNG/PDF plots in viz_summary/")
}

# --------- Run all ---------
invisible(lapply(c("A","B","C"), analyze_ds))
cat("\n✔ All done. Outputs in:\n  ", VIZ, "\n", sep = "")


library(ggplot2)
library(patchwork)
library(igraph)
library(scales)
library(data.table)

# --- Load interaction results for C ---
edges <- fread("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/C_snp_interactions.tsv")
edges <- edges[!is.na(p)]

# --- Panel A: Heatmap ---
edges_grid <- edges[, .(Var1 = V1, Var2 = V2, p)]
p1 <- ggplot(edges_grid, aes(Var2, Var1, fill = -log10(p))) +
  geom_tile(color = "grey90") +
  scale_fill_gradient(low = "white", high = "black", na.value = "grey90") +
  labs(title = "A: SNP×SNP −log10(p)", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Panel B: Network ---
sig <- edges[p < 0.05]
g <- graph_from_data_frame(sig[, .(from = V1, to = V2, weight = -log10(p))])
p2 <- ggraph::ggraph(g, layout = "fr") +
  ggraph::geom_edge_link(aes(width = weight), color = "grey50") +
  ggraph::geom_node_point(size = 6, color = "tomato") +
  ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  labs(title = "B: Significant SNP–SNP network (p<0.05)")

# --- Panel C: Chromosome schematic ---
chr6 <- data.table(pos = c(32635435, 32655446))
p3 <- ggplot(chr6, aes(x = pos, y = 1)) +
  geom_segment(x = 32630000, xend = 32670000, y = 1, yend = 1, color = "grey80", size = 2) +
  geom_point(aes(color = as.factor(pos)), size = 5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8")) +
  labs(title = "C: Interacting SNPs (Chr6, 32.63–32.66 Mb)",
       x = "Position (Mb)", y = NULL, color = "SNP") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_blank())

# Combine
final_plot <- (p1 | p2) / p3
ggsave("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/HLA_interaction_composite.png",
       final_plot, width = 9, height = 8, dpi = 300)

# --- Map SNPs to genes manually based on Ensembl / HLA position ---
mapped_genes <- c("HLA-DQB1", "HLA-DQA2")


suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(igraph)
  library(ggraph)
  library(data.table)
  library(scales)
})

# --- SNP coordinates and mapping ---
snp_map <- data.table(
  chrpos = c("06:032635435","06:032655446"),
  rsid = c("rs9469200","rs2395228"),
  pos = c(32635435,32655446)
)
snp_map[, label := paste0(rsid, " (", chrpos, ")")]

# --- Load interaction results for dataset C ---
edges <- fread("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/C_snp_interactions.tsv")
edges <- merge(edges, snp_map, by.x = "V1", by.y = "chrpos", all.x = TRUE)
edges <- merge(edges, snp_map, by.x = "V2", by.y = "chrpos", all.x = TRUE,
               suffixes = c("_1","_2"))

# --- Panel A: Heatmap of –log10(p) ---
p1 <- ggplot(edges, aes(V2, V1, fill = -log10(p))) +
  geom_tile(color = "grey90") +
  scale_fill_gradient(low = "white", high = "black", na.value = "grey85") +
  labs(title = "A: SNP×SNP −log10(p)", x = NULL, y = NULL, fill = "−log10(p)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# --- Panel B: Network with rsIDs and positions ---
sig <- edges[p < 0.05]
if (nrow(sig) > 0) {
  g <- graph_from_data_frame(sig[, .(from = paste0(rsid_1, " (", V1, ")"),
                                     to = paste0(rsid_2, " (", V2, ")"),
                                     weight = -log10(p))])
  p2 <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = weight), color = "grey50") +
    geom_node_point(size = 8, color = "tomato") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5) +
    theme_void() +
    labs(title = "B: Significant SNP–SNP network (p<0.05)",
         subtitle = "Showing rsIDs and CHR:POS for significant pair(s)")
} else {
  p2 <- ggplot() + theme_void() + 
    labs(title = "B: No significant interactions (p<0.05)")
}

# --- Panel C: Genomic schematic with rsID + POS labels (fixed margins) ---
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# --- Paths & SNPs ---
geno_fp <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/C_homog_numeric.rds"
snps_chrpos <- c("06:032635435","06:032655446")   # rs9469200, rs2395228
rsid_map <- c("06:032635435"="rs9469200", "06:032655446"="rs2395228")

# --- Load genotype matrix (handles matrix or list-of-matrix) ---
X <- readRDS(geno_fp)
if (is.list(X)) X <- X[[1]]
stopifnot(is.matrix(X) || is.data.frame(X))

# safety: only keep the two columns that exist
have <- intersect(snps_chrpos, colnames(X))
if (length(have) != 2) stop("Couldn't find both SNPs in C matrix. Found: ", paste(have, collapse=", "))

# --- Build a small df for plotting + counts ---
counts_dt <- rbindlist(lapply(have, function(s) {
  v <- X[, s]
  # Coerce to integer 0/1/2 if needed
  if (!is.integer(v)) v <- as.integer(round(v))
  tab <- as.integer(table(factor(v, levels = c(0,1,2))))
  data.table(chrpos = s,
             rsid   = rsid_map[[s]],
             pos    = as.integer(sub("^..:(\\d+)$","\\1", s)),
             n0 = tab[1], n1 = tab[2], n2 = tab[3])
}))
counts_dt[, label := sprintf("0/1/2 = %d/%d/%d", n0, n1, n2)]

# --- Panel C (fixed) ---
# pad the x-range to the right so labels never clip
xr <- range(counts_dt$pos)
pad <- round(0.06 * diff(xr))  # ~6% padding
xmin <- xr[1] - pad
xmax <- xr[2] + pad

p3 <- ggplot(counts_dt, aes(x = pos, y = 1, color = rsid)) +
  # genome “track”
  geom_segment(x = xmin, xend = xmax, y = 1, yend = 1, color = "grey80", linewidth = 2) +
  # points
  geom_point(size = 6) +
  # main labels: rsID + CHR:POS (stacked)
  geom_text_repel(aes(label = paste0(rsid, "\n", chrpos)),
                  nudge_y = 0.18, min.segment.length = 0,
                  box.padding = 0.3, point.padding = 0.3,
                  size = 4, fontface = "bold", lineheight = 0.95, seed = 42) +
  # genotype count labels (0/1/2) printed below points
  geom_text(aes(label = label), vjust = 1.8, size = 3.6, fontface = "plain") +
  scale_color_manual(values = c("rs9469200"="#377EB8","rs2395228"="#E41A1C")) +
  labs(title = "C: Interacting SNPs (Chr6, 32.63–32.66 Mb)",
       x = "Position (Mb)", y = NULL, color = NULL) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = scales::pretty_breaks(6),
                     labels = function(z) sprintf("%.3f", z/1e6)) +
  coord_cartesian(clip = "off") +                      # ← prevents clipping
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(t = 20, r = 80, b = 30, l = 80)  # generous margins
  ) +
  expand_limits(y = c(0.7, 1.35))  # vertical breathing room

p3 <- p3 +
  labs(
    subtitle = "Genotype counts (beneath each rsID): 0=homozygous ref, 1=heterozygous, 2=homozygous alt"
  ) +
  theme(
    plot.subtitle = element_text(size = 10, margin = margin(b = 6))
  ) +
  # optional sticky legend text inside the panel (bottom-left)
  annotate(
    "label", x = min(snp_map$pos)/1e6 + 0.003, y = 0.92,
    label = "Legend: 0 = homozygous reference\n          1 = heterozygous\n          2 = homozygous alternative",
    size = 3.2, hjust = 0, vjust = 1,
    label.size = 0.2, color = "grey20", fill = "white"
  )
# Now rebuild your composite:
# final_plot <- (p1 | p2) / p3
## =======================
## Panel D — KEGG-like pathway schematic
## =======================
library(ggplot2)
library(patchwork)
library(grid)

# --- helper to draw boxes with labels (now defined early) ---
node <- function(xmin, xmax, ymin, ymax, label,
                 fill = "#FEE0D2", col = "#CB181D") {
  list(
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = fill, color = col, linewidth = 0.5),
    annotate("text",
             x = (xmin + xmax)/2, y = (ymin + ymax)/2,
             label = label, size = 3.4, fontface = "bold")
  )
}

# --- simple schematic showing where HLA-DQA2/DQB1 sit in antigen presentation ---
p4 <- ggplot() +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 7,
           fill = "#FAFAFA", color = "#DDDDDD") +
  
  node(0.5, 3.2, 5.6, 6.6, "APC / Endosome\n(MHC-II loading)",
       fill = "#EFF3FF", col = "#6BAED6") +
  node(3.5, 4.8, 5.6, 6.6, "HLA-DM/DO",
       fill = "#EFF3FF", col = "#6BAED6") +
  node(5.1, 8.8, 5.6, 6.6, "HLA-DQA2 / HLA-DQB1\n(MHC-II α/β)",
       fill = "#FEE0D2", col = "#CB181D") +
  node(5.1, 8.8, 3.9, 4.9, "Peptide–MHC-II complex",
       fill = "#FFF7BC", col = "#D95F0E") +
  node(5.1, 6.8, 2.2, 3.2, "TCR",
       fill = "#ECECEC", col = "#9E9E9E") +
  node(7.1, 8.8, 2.2, 3.2, "CD4",
       fill = "#ECECEC", col = "#9E9E9E") +
  
  geom_curve(aes(x = 3.2, y = 6.1, xend = 3.5, yend = 6.1),
             curvature = 0.05, arrow = arrow(length = unit(3, "pt")),
             color = "#6BAED6") +
  geom_curve(aes(x = 4.8, y = 6.1, xend = 5.1, yend = 6.1),
             curvature = 0.05, arrow = arrow(length = unit(3, "pt")),
             color = "#6BAED6") +
  geom_curve(aes(x = 7.0, y = 5.6, xend = 7.0, yend = 4.9),
             curvature = 0, arrow = arrow(length = unit(3, "pt")),
             color = "#D95F0E") +
  geom_curve(aes(x = 6.0, y = 3.9, xend = 6.0, yend = 3.2),
             curvature = 0, arrow = arrow(length = unit(3, "pt")),
             color = "#9E9E9E") +
  geom_curve(aes(x = 7.9, y = 3.9, xend = 7.9, yend = 3.2),
             curvature = 0, arrow = arrow(length = unit(3, "pt")),
             color = "#9E9E9E") +
  
  # highlight our two interacting SNPs
  annotate("label",
           x = 6.95, y = 6.85,
           label = "rs9469200 (06:032635435)\nrs2395228 (06:032655446)",
           size = 3, label.size = 0.2, fill = "white") +
  
  labs(
    title = "D: Antigen processing & presentation (KEGG hsa04612)",
    subtitle = "Highlighted: MHC-II α/β (HLA-DQA2 / HLA-DQB1) containing the interacting SNPs",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold")) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 8), expand = FALSE)

## =======================
## Combine panels and save
## =======================
final_plot <- (p1 | p2) / (p3 | p4)

ggsave(
  "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/HLA_interaction_composite_ABCD.png",
  final_plot, width = 14, height = 10, dpi = 300
)
ggsave(
  "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/grid_B02500/viz_summary/HLA_interaction_composite_ABCD.pdf",
  final_plot, width = 14, height = 10
)
cat("\n✅ Saved composite figure with rsIDs + positions:\n  - HLA_interaction_composite.png\n  - HLA_interaction_composite.pdf\n")