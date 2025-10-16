# --- CONFIG ---
BASE <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
GRID <- "grid_B02500"               # bootstrap folder
TOPN <- 20                          # top SNPs per dataset
OUTDIR <- file.path(BASE, GRID, "qc_topSNP_overlap")

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(glue)
  library(ggVennDiagram)
})

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Helper: path to top loadings (PC1)
topfile <- function(ds) file.path(BASE, sprintf("step2a_pca_irlba_%s", ds), "pca_top_loadings_PC1.tsv")

# Read top-N SNPs by |loading|
read_topN <- function(ds) {
  f <- topfile(ds)
  if (!file.exists(f)) stop(glue("Missing file for {ds}: {f}"))
  dt <- fread(f)
  # find columns
  snp_col <- intersect(c("SNP","snp","variant"), names(dt))[1]
  load_col <- intersect(c("abs_loading","abs_PC1_loading","loading","PC1_loading"), names(dt))[1]
  if (is.na(snp_col)) stop(glue("No SNP column found in {f}"))
  if (is.na(load_col)) stop(glue("No loading column found in {f}"))
  dt %>%
    mutate(abs_loading = abs(.data[[load_col]])) %>%
    arrange(desc(abs_loading)) %>%
    slice_head(n = TOPN) %>%
    transmute(dataset = ds, SNP = .data[[snp_col]], abs_loading)
}

topA <- read_topN("A")
topB <- read_topN("B")
topC <- read_topN("C")

# --- Merge union + overlaps ---
union_tbl <- bind_rows(topA, topB, topC) %>%
  mutate(in_A = dataset=="A",
         in_B = dataset=="B",
         in_C = dataset=="C") %>%
  group_by(SNP) %>%
  summarise(in_A = any(in_A), in_B = any(in_B), in_C = any(in_C),
            max_abs_loading = max(abs_loading, na.rm=TRUE),
            .groups="drop") %>%
  arrange(desc(max_abs_loading))

fwrite(union_tbl, file.path(OUTDIR, "top20_union_membership_genetic.tsv"))

# Overlap counts
overlap_counts <- tibble(
  region = c("A_only","B_only","C_only","A∩B_only","A∩C_only","B∩C_only","A∩B∩C"),
  count  = c(
    sum(union_tbl$in_A & !union_tbl$in_B & !union_tbl$in_C),
    sum(union_tbl$in_B & !union_tbl$in_A & !union_tbl$in_C),
    sum(union_tbl$in_C & !union_tbl$in_A & !union_tbl$in_B),
    sum(union_tbl$in_A & union_tbl$in_B & !union_tbl$in_C),
    sum(union_tbl$in_A & union_tbl$in_C & !union_tbl$in_B),
    sum(union_tbl$in_B & union_tbl$in_C & !union_tbl$in_A),
    sum(union_tbl$in_A & union_tbl$in_B & union_tbl$in_C)
  )
)
fwrite(overlap_counts, file.path(OUTDIR, "top20_overlap_counts_genetic.tsv"))

# Venn diagram (A/B/C)
sets <- list(
  A = topA$SNP,
  B = topB$SNP,
  C = topC$SNP
)

p <- ggVennDiagram(sets, label_alpha = 0, label_size = 5, edge_size = 1.1) +
  scale_fill_gradient(low = "#cfe2f3", high = "#08519c") +
  ggtitle(glue("Top-{TOPN} SNPs by |PC1 loading| — feature = genetic")) +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(OUTDIR, "top20_overlap_venn_genetic.jpeg"), p,
       width = 8, height = 7, dpi = 300)

# Also save raw source table
source_tbl <- bind_rows(
  topA %>% transmute(SNP, source_dataset="A", abs_loading),
  topB %>% transmute(SNP, source_dataset="B", abs_loading),
  topC %>% transmute(SNP, source_dataset="C", abs_loading)
) %>%
  arrange(desc(abs_loading))
fwrite(source_tbl, file.path(OUTDIR, "topSNP_loading_sources_genetic.tsv"))

message("\n✅ Done!")
message("Outputs written to: ", OUTDIR)

# ===========================
# Top-20 PC1 SNP overlap + Fisher tests
# ===========================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2)
  library(ggVennDiagram)
})

# --- Paths ---
BASE <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min"
OUTDIR <- file.path(BASE, "grid_B02500", "viz_summary")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# --- Helper: load top-20 by |loading| for one dataset (feature=genetic) ---
load_top20 <- function(ds) {
  # Use the PCA loadings produced in step2a (irlba) for this dataset
  f <- file.path(BASE, sprintf("step2a_pca_irlba_%s", ds), "pca_top_loadings_PC1.tsv")
  stopifnot(file.exists(f))
  dt <- fread(f)
  # Try to be robust to column naming:
  # common formats: "SNP", "variant", "marker"; loading column: "loading", "PC1_loading", "abs_loading"
  snp_col <- intersect(c("SNP","variant","marker","id"), names(dt))[1]
  load_col <- intersect(c("loading","PC1_loading","abs_loading","value"), names(dt))[1]
  if (is.na(snp_col) || is.na(load_col)) {
    stop(sprintf("Cannot find SNP/loading columns in %s; have: %s", f, paste(names(dt), collapse=", ")))
  }
  dt <- dt %>%
    mutate(abs_loading = abs(.data[[load_col]])) %>%
    arrange(desc(abs_loading)) %>%
    transmute(SNP = .data[[snp_col]], abs_loading) %>%
    distinct(SNP, .keep_all = TRUE) %>%
    slice_head(n = 20)
  dt
}

topA <- load_top20("A")
topB <- load_top20("B")
topC <- load_top20("C")

sets <- list(
  A = topA$SNP,
  B = topB$SNP,
  C = topC$SNP
)

# --- Build region membership (no poking at ggVenn internals) ---
A <- sets$A; B <- sets$B; C <- sets$C
AB  <- intersect(A,B); AC <- intersect(A,C); BC <- intersect(B,C)
ABC <- intersect(AB, C)

regions <- list(
  `A_only`  = setdiff(A, union(B, C)),
  `B_only`  = setdiff(B, union(A, C)),
  `C_only`  = setdiff(C, union(A, B)),
  `A&B`     = setdiff(AB, C),
  `A&C`     = setdiff(AC, B),
  `B&C`     = setdiff(BC, A),
  `A&B&C`   = ABC
)

# Save region table with SNP names
region_tbl <- tibble(
  region = names(regions),
  count  = sapply(regions, length),
  SNPs   = sapply(regions, function(x) paste(x, collapse = ", "))
)
fwrite(region_tbl, file.path(OUTDIR, "top20_overlap_snps_by_region.tsv"))

# --- Fisher's exact tests (pairwise) ---
# Universe = union of all top-20 lists
U <- union(A, union(B, C))

fisher_pair <- function(X, Y, nameX, nameY) {
  a <- length(intersect(X, Y))               # in both
  b <- length(setdiff(X, Y))                 # in X only
  c <- length(setdiff(Y, X))                 # in Y only
  d <- length(setdiff(U, union(X, Y)))       # in neither (within U)
  mat <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE,
                dimnames = list(c("inY","not_inY"), c("inX","not_inX")))
  ft <- fisher.test(mat, alternative = "greater")
  tibble(pair = paste0(nameX, " vs ", nameY),
         a_both = a, b_Xonly = b, c_Yonly = c, d_neither = d,
         odds_ratio = unname(ft$estimate),
         p_value = ft$p.value)
}

fisher_results <- bind_rows(
  fisher_pair(A, B, "A", "B"),
  fisher_pair(A, C, "A", "C"),
  fisher_pair(B, C, "B", "C")
)
fwrite(fisher_results, file.path(OUTDIR, "top20_pairwise_fisher.tsv"))

# --- Venn plot with counts *and* SNP labels in a caption file ---
# (Directly writing all SNPs into the shapes gets cluttered; we keep the figure clean
#  and store full labels in the TSV above.)
venn_plot <- ggVennDiagram(sets, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "#d6eaf8", high = "#1b4f72") +
  ggtitle("Top-20 SNPs by |PC1 loading| — feature = genetic") +
  theme(plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(OUTDIR, "top20_overlap_venn_genetic.jpeg"),
       venn_plot, width = 9, height = 8, dpi = 300)

# --- Also print fisher results to console for quick view ---
print(fisher_results)
message("✅ Wrote:",
        "\n  • ", file.path(OUTDIR, "top20_overlap_venn_genetic.jpeg"),
        "\n  • ", file.path(OUTDIR, "top20_overlap_snps_by_region.tsv"),
        "\n  • ", file.path(OUTDIR, "top20_pairwise_fisher.tsv"))