setwd("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/")
# ------------------ USER PATHS ------------------
CLIN_PATH <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/dat/cluster_clinical.csv"
PCA_SCORES_GLOB <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/cluster_runs_min/step2a_pca_*_*/pca_scores.tsv"
# ------------------------------------------------

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(glue)
})

canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

# ---------- 1) Read clinical + basic cleaning ----------
stopifnot(file.exists(CLIN_PATH))
clin <- suppressMessages(readr::read_csv(CLIN_PATH, show_col_types = FALSE))

# standardize ID column name
id_col <- c("Cluster_subjectID","IID","subject_id","id")
id_col <- id_col[id_col %in% names(clin)][1]
if (is.na(id_col)) stop("Could not find an ID column in clinical file.")
clin <- clin %>% rename(SubjectID = !!id_col)

# trim/normalize IDs
clin$SubjectID <- canon_id(clin$SubjectID)

# try to create Cluster_sex_num if not present
if (!"Cluster_sex_num" %in% names(clin)) {
  if ("Cluster_sex" %in% names(clin)) {
    sx <- tolower(as.character(clin$Cluster_sex))
    clin$Cluster_sex_num <- ifelse(grepl("^m", sx), 1,
                                   ifelse(grepl("^f", sx), 0, NA_real_))
  } else {
    clin$Cluster_sex_num <- NA_real_
  }
}

# ---------- 2) Column completeness table ----------
# % non-missing for every column
comp <- clin %>%
  summarize(across(everything(), ~ mean(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "coverage") %>%
  arrange(desc(coverage))

# show a quick summary
cat("\n=== Clinical completeness (top 25 by coverage) ===\n")
print(head(comp, 25), n = 25)
cat("\nLowest-coverage columns:\n")
print(tail(comp, 10), n = 10)

# highlight key fields you care about
key_fields <- c("Cluster_age","Cluster_sex","Cluster_sex_num")
key_fields <- key_fields[key_fields %in% names(clin)]
if (length(key_fields)) {
  cat("\nKey fields coverage:\n")
  print(dplyr::filter(comp, column %in% key_fields))
}

# duplicates / empty IDs
dup_ids <- clin %>%
  count(SubjectID, name = "n") %>% filter(n > 1)
n_empty <- sum(!nzchar(clin$SubjectID))
cat(glue("\nDuplicates in SubjectID: {nrow(dup_ids)} rows\n"))
if (nrow(dup_ids)) print(dup_ids, n = nrow(dup_ids))
cat(glue("Empty/blank SubjectID: {n_empty}\n"))

# ---------- 3) Row completeness for a chosen set ----------
# (a) minimal set = age + sex_num
need_min <- intersect(c("Cluster_age","Cluster_sex_num"), names(clin))
n_min_complete <- if (length(need_min)) sum(complete.cases(clin[, need_min, drop=FALSE])) else NA

# (b) “full clinical” = everything except ID-like columns
exclude <- c("SubjectID","Cluster_subjectID","IID")
have_full <- setdiff(names(clin), exclude)
n_full_complete <- sum(complete.cases(clin[, have_full, drop=FALSE]))

cat(glue("\nRow completeness:\n  N total clinical rows: {nrow(clin)}\n",
         if (!is.na(n_min_complete)) glue("  Complete for [age+sex_num]: {n_min_complete}\n") else "",
         "  Complete for [ALL selected clinical columns]: {n_full_complete}\n"))

# Which columns are driving the missingness most?
miss_impact <- clin %>%
  mutate(.row_ok = complete.cases(across(all_of(have_full)))) %>%
  filter(!.row_ok) %>%
  summarize(across(all_of(have_full), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "n_missing") %>%
  arrange(desc(n_missing))
cat("\nColumns most responsible for rows being dropped in FULL set:\n")
print(head(miss_impact, 20), n = 20)

# ---------- 4) Compare with PCA IDs (A/B/C) ----------
pca_files <- Sys.glob(PCA_SCORES_GLOB)
if (length(pca_files)) {
  cat("\n=== Overlap with PCA scores ===\n")
  for (ff in pca_files) {
    ds <- sub(".*step2a_pca_([a-z]+)_([ABC]).*", "\\2", ff, ignore.case = TRUE)
    pcs <- suppressMessages(readr::read_tsv(ff, show_col_types = FALSE))
    idcol <- c("IID","id","sample","SubjectID")
    idcol <- idcol[idcol %in% names(pcs)][1]
    if (is.na(idcol)) next
    pcs_ids <- canon_id(pcs[[idcol]])
    ovl <- length(intersect(pcs_ids, clin$SubjectID))
    cat(glue("  [{ds}] PCA N={length(unique(pcs_ids))} | Clinical N={nrow(clin)} | Overlap={ovl}\n"))
  }
} else {
  cat("\n(No PCA score files found with the current glob.)\n")
}

# ---------- 5) Optional: write a TSV report ----------
out_dir <- dirname(CLIN_PATH)
readr::write_tsv(comp, file.path(out_dir, "clinical_column_coverage.tsv"))
readr::write_tsv(miss_impact, file.path(out_dir, "clinical_missingness_drivers.tsv"))

####Cleaner

CLIN_PATH <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/dat/cluster_clinical.csv"

library(readr); library(dplyr); library(stringr)

clin <- suppressMessages(readr::read_csv(CLIN_PATH, show_col_types = FALSE))

# Robust derivation from Cluster_sex → Cluster_sex_num (M/male/1 → 1, F/female/0 → 0, else NA)
to_sex_num <- function(x) {
  sx <- tolower(trimws(as.character(x)))
  out <- ifelse(grepl("^m", sx) | sx == "1", 1,
                ifelse(grepl("^f", sx) | sx == "0", 0, NA_real_))
  as.numeric(out)
}

# If Cluster_sex_num is missing OR all NA, (re)create it from Cluster_sex
if (!"Cluster_sex_num" %in% names(clin) || all(is.na(clin$Cluster_sex_num))) {
  message("Creating/overwriting Cluster_sex_num from Cluster_sex …")
  clin$Cluster_sex_num <- to_sex_num(clin$Cluster_sex)
} else {
  # Fill only the NA entries using Cluster_sex
  idx <- is.na(clin$Cluster_sex_num)
  if (any(idx)) {
    message("Filling NA values in Cluster_sex_num from Cluster_sex …")
    clin$Cluster_sex_num[idx] <- to_sex_num(clin$Cluster_sex[idx])
  }
}

# Sanity: force invalid values to NA (only 0/1 allowed)
clin$Cluster_sex_num[!(clin$Cluster_sex_num %in% c(0,1))] <- NA_real_

# Quick check
cat("Non-missing counts:\n")
print(colSums(!is.na(clin[, c("Cluster_sex","Cluster_sex_num","Cluster_age")])), row.names = FALSE)

# Save back (overwrite the clinical CSV)
readr::write_csv(clin, CLIN_PATH)
message("Wrote cleaned clinical file: ", CLIN_PATH)




#!/usr/bin/env Rscript

# ================== CONFIG ==================
setwd("/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/")
RAW_PATH  <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/dat/cluster_clinical.csv"
CLEAN_PATH <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/dat/cluster_clinical_clean.csv"
# ============================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(glue)
})

# ------------------ HELPERS ------------------
canon_id <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

to_sex_num <- function(x) {
  sx <- tolower(trimws(as.character(x)))
  out <- ifelse(grepl("^m", sx) | sx == "1", 1,
                ifelse(grepl("^f", sx) | sx == "0", 0, NA_real_))
  as.numeric(out)
}

# ------------------ LOAD ------------------
stopifnot(file.exists(RAW_PATH))
clin <- suppressMessages(read_csv(RAW_PATH, show_col_types = FALSE))

message(glue("Loaded raw clinical file: {RAW_PATH}"))
message(glue("Rows: {nrow(clin)} | Columns: {ncol(clin)}"))

# Standardize ID column name
id_col <- c("Cluster_subjectID","IID","subject_id","id")
id_col <- id_col[id_col %in% names(clin)][1]
if (is.na(id_col)) stop("Could not find an ID column in clinical file.")
clin <- clin %>% rename(SubjectID = !!id_col)
clin$SubjectID <- canon_id(clin$SubjectID)

# ------------------ FIX SEX_NUM ------------------
if (!"Cluster_sex_num" %in% names(clin) || all(is.na(clin$Cluster_sex_num))) {
  message("Creating new 'Cluster_sex_num' from 'Cluster_sex'...")
  clin$Cluster_sex_num <- to_sex_num(clin$Cluster_sex)
} else {
  idx <- is.na(clin$Cluster_sex_num)
  if (any(idx)) {
    message(glue("Filling {sum(idx)} NA values in 'Cluster_sex_num' from 'Cluster_sex'..."))
    clin$Cluster_sex_num[idx] <- to_sex_num(clin$Cluster_sex[idx])
  }
}
clin$Cluster_sex_num[!(clin$Cluster_sex_num %in% c(0,1))] <- NA_real_

# ------------------ QUICK IMPUTATION ------------------
# optional — ensures no row loss in clustering
if ("Cluster_age" %in% names(clin)) {
  med_age <- median(clin$Cluster_age, na.rm = TRUE)
  clin$Cluster_age[is.na(clin$Cluster_age)] <- med_age
}

if ("Cluster_sex_num" %in% names(clin)) {
  mode_sex <- as.numeric(names(sort(table(clin$Cluster_sex_num), decreasing = TRUE)[1]))
  clin$Cluster_sex_num[is.na(clin$Cluster_sex_num)] <- mode_sex
}

# ------------------ CHECKS ------------------
coverage <- clin %>%
  summarize(across(everything(), ~ mean(!is.na(.)))) %>%
  pivot_longer(everything(), names_to = "column", values_to = "coverage") %>%
  arrange(desc(coverage))

key_fields <- c("Cluster_age","Cluster_sex","Cluster_sex_num")
key_fields <- key_fields[key_fields %in% names(clin)]

cat("\n=== Key Fields Coverage ===\n")
print(filter(coverage, column %in% key_fields))

# ------------------ SAVE CLEAN FILE ------------------
write_csv(clin, CLEAN_PATH)
message(glue("\n✅ Wrote cleaned clinical file: {CLEAN_PATH}"))
message(glue("Rows: {nrow(clin)} | Columns: {ncol(clin)}"))

# ------------------ VALIDATION SUMMARY ------------------
cat("\n=== Clean File Summary ===\n")
cat(glue("Cluster_sex_num unique values: {toString(sort(unique(clin$Cluster_sex_num)))}\n"))
cat(glue("Cluster_age range: {range(clin$Cluster_age, na.rm = TRUE)}\n"))

cov_sex <- mean(!is.na(clin$Cluster_sex_num))
cov_age <- mean(!is.na(clin$Cluster_age))
cat(glue("\nCoverage →  sex_num: {round(100*cov_sex,1)}%   |   age: {round(100*cov_age,1)}%\n"))
if (cov_sex < 0.95) warning("⚠️ Coverage for Cluster_sex_num is still low. Check source file.")

cat("\nDone.\n")


library(readr)
library(dplyr)

# path to your file
clin_path <- "/Users/jose/Documents/PHD experimental Medicine/ISPAD_Paper/out_three_runs/dat/cluster_clinical.csv"

clin <- read_csv(clin_path, show_col_types = FALSE)

# Fix Cluster_sex_num:  -1 = NA → then impute by mode
clin <- clin %>%
  mutate(
    Cluster_sex_num = case_when(
      Cluster_sex == 1 ~ 1,     # male
      Cluster_sex == 2 ~ 0,     # female
      Cluster_sex == -1 ~ NA_real_,  # unknown
      TRUE ~ Cluster_sex_num
    )
  )

# Replace NA by most common value (mode)
mode_sex <- as.numeric(names(sort(table(clin$Cluster_sex_num), decreasing = TRUE)[1]))
clin$Cluster_sex_num[is.na(clin$Cluster_sex_num)] <- mode_sex

# Optional: fill missing ages by median
if ("Cluster_age" %in% names(clin)) {
  clin$Cluster_age[is.na(clin$Cluster_age)] <- median(clin$Cluster_age, na.rm = TRUE)
}

# Save back same name so step3 picks it up automatically
write_csv(clin, clin_path)

cat("\n✅ Clinical file fixed and saved in place.\n")
table(clin$Cluster_sex_num, useNA="ifany")
