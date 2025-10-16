setwd("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/")

# 1) Load (nothing runs yet)
source("step3a_clustering_rstudio.R")

##Sanity Checks
# Make a quick overlap report for each dataset
cfg <- default_step3_config(); cfg <- resolve_paths(cfg)
clin <- load_clinical(cfg)

canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))
clin$Cluster_subjectID_canon <- canon_id(clin$Cluster_subjectID)

for (ds in c("A","B","C")) {
  sc  <- load_scores(cfg, ds)
  sc$IID_canon <- canon_id(sc$IID)
  
  overlap <- intersect(sc$IID_canon, clin$Cluster_subjectID_canon)
  only_pca  <- setdiff(sc$IID_canon, clin$Cluster_subjectID_canon)
  only_clin <- setdiff(clin$Cluster_subjectID_canon, sc$IID_canon)
  
  cat(sprintf("\n[%s] PCA IDs=%d | Clinical IDs=%d | Overlap=%d | PCA-only=%d | Clin-only=%d\n",
              ds, length(unique(sc$IID_canon)), length(unique(clin$Cluster_subjectID_canon)),
              length(overlap), length(only_pca), length(only_clin)))
  
  # write a small report for manual inspection
  out <- data.frame(
    PCA_IID_raw   = head(sc$IID, 30),
    PCA_IID_canon = head(sc$IID_canon, 30)
  )
  out2 <- data.frame(
    Clin_ID_raw   = head(clin$Cluster_subjectID, 30),
    Clin_ID_canon = head(clin$Cluster_subjectID_canon, 30)
  )
  data.table::fwrite(out,  file.path(cfg$RUN_DIR, sprintf("step3_id_preview_pca_%s.tsv", ds)), sep="\t")
  data.table::fwrite(out2, file.path(cfg$RUN_DIR, sprintf("step3_id_preview_clin_%s.tsv", ds)), sep="\t")
  
  data.table::fwrite(data.frame(only_pca=head(only_pca,100)),
                     file.path(cfg$RUN_DIR, sprintf("step3_id_only_pca_%s.tsv", ds)), sep="\t")
  data.table::fwrite(data.frame(only_clin=head(only_clin,100)),
                     file.path(cfg$RUN_DIR, sprintf("step3_id_only_clin_%s.tsv", ds)), sep="\t")
}

cfg <- default_step3_config(); cfg <- resolve_paths(cfg)
clin <- load_clinical(cfg)

table(is.na(clin$Cluster_age))
table(is.na(clin$Cluster_sex_num))

mean(complete.cases(clin[, c("Cluster_age", "Cluster_sex_num")]))
# 2) Configure
cfg <- default_step3_config()
cfg$DATASETS <- c("A","B","C")  # pick any
cfg$N_PCS    <- 15
cfg$BOOT_N   <- 1500             # raise to 1000 later
cfg$K_RANGE  <- 2:6
# 3) Go
run_step3(cfg)