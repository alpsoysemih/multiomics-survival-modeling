############################################################
# Meta-analysis of Cox models: METABRIC + TCGA-BRCA
# - Sheet-level models preserved
# - E / C / E_C / E_C_int tables
# - Filtered sheets by pooled CI_lo > 1
# - Remove rows where rounded pooled CI_lo == 1
# - Sort filtered sheets:
#     gene (ascending) + pooled HR (descending)
# - 3 significant digits for all numerics
# - Column sets obey E / C / int rules
############################################################

library(tidyverse)
library(readxl)
library(openxlsx)
library(stringr)

############################################################
# Paths
############################################################

tcga_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Results/TCGA-BRCA"
metabric_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/METABRIC"
results_dir  <- "/Volumes/Expansion/Prognostic_Analysis/Results/Meta-Analysis"

if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

tcga_file     <- file.path(tcga_dir,     "TCGA_BRCA_OS_all_Cox_all_models_by_gene_with_LRT.xlsx")
metabric_file <- file.path(metabric_dir, "METABRIC_OS_all_Cox_all_models_by_gene_with_LRT.xlsx")
out_file      <- file.path(results_dir,  "METABRIC_TCGA-BRCA_OS_Gene_Meta-Analysis.xlsx")

############################################################
# Helper functions
############################################################

# Read all sheets from a cohort file, excluding unwanted sheet names.
read_all_sheets <- function(path, cohort_name) {
  exclude_pat <- c("_filter", "_filter2", "LR_expr", "LR_cna", "LR_baseline")
  sheets <- excel_sheets(path)
  sheets <- sheets[!str_detect(sheets, paste(exclude_pat, collapse = "|"))]
  
  out <- lapply(sheets, function(sh) {
    df <- read_excel(path, sheet = sh) %>% as_tibble()
    # Make sure gene and model exist
    if (!"gene" %in% names(df)) stop(paste("Sheet has no 'gene' column:", sh))
    if (!"model" %in% names(df)) df$model <- sh
    df$sheet  <- sh
    df$cohort <- cohort_name
    df
  })
  names(out) <- sheets
  out
}

# Coerce numeric-like columns (including possible comma decimals) to numeric.
coerce_numeric_cols <- function(df) {
  num_candidates <- c(
    "n", "events", "c_index",
    "HR_expr", "CI_lo_expr", "CI_hi_expr",
    "HR_cna",  "CI_lo_cna",  "CI_hi_cna",
    "HR_expr_cna", "CI_lo_expr_cna", "CI_hi_expr_cna"
  )
  for (nm in intersect(num_candidates, names(df))) {
    df[[nm]] <- df[[nm]] %>%
      as.character() %>%
      gsub(",", ".", ., fixed = FALSE) %>%
      as.numeric()
  }
  df
}

# Ensure all needed columns exist (missing ones as NA_real_)
ensure_needed_cols <- function(df) {
  needed <- c(
    "n", "events", "c_index",
    "HR_expr", "CI_lo_expr", "CI_hi_expr",
    "HR_cna",  "CI_lo_cna",  "CI_hi_cna",
    "HR_expr_cna", "CI_lo_expr_cna", "CI_hi_expr_cna"
  )
  missing <- setdiff(needed, names(df))
  for (m in missing) df[[m]] <- NA_real_
  df
}

# Standard error from confidence interval (log-scale).
se_from_ci <- function(ci_lo, ci_hi) {
  (log(ci_hi) - log(ci_lo)) / (2 * 1.96)
}

# Format numeric columns to 3 significant digits.
fmt_all_numeric <- function(df) {
  df %>% mutate(across(where(is.numeric), ~ signif(., digits = 3)))
}

############################################################
# Load METABRIC and TCGA-BRCA, sheet-wise
############################################################

met_ls  <- read_all_sheets(metabric_file, "METABRIC")
tcga_ls <- read_all_sheets(tcga_file,     "TCGA-BRCA")

met_ls  <- lapply(met_ls,  coerce_numeric_cols)
tcga_ls <- lapply(tcga_ls, coerce_numeric_cols)

met_ls  <- lapply(met_ls,  ensure_needed_cols)
tcga_ls <- lapply(tcga_ls, ensure_needed_cols)

############################################################
# Join METABRIC and TCGA per sheet/model
############################################################

common_models <- intersect(names(met_ls), names(tcga_ls))

joined <- bind_rows(lapply(common_models, function(m) {
  df_m <- met_ls[[m]]
  df_t <- tcga_ls[[m]]
  
  # Ensure both have 'model' column and same model name
  if (!"model" %in% names(df_m)) df_m$model <- m
  if (!"model" %in% names(df_t)) df_t$model <- m
  
  inner_join(
    df_m,
    df_t,
    by = c("gene", "model"),
    suffix = c("_METABRIC", "_TCGA")
  )
}))

############################################################
# Classify model type based on available HR columns
############################################################

joined <- joined %>%
  mutate(
    has_expr = !is.na(HR_expr_METABRIC) & !is.na(HR_expr_TCGA),
    has_cna  = !is.na(HR_cna_METABRIC)  & !is.na(HR_cna_TCGA),
    has_int  = !is.na(HR_expr_cna_METABRIC) & !is.na(HR_expr_cna_TCGA),
    model_type = case_when(
      has_expr & !has_cna & !has_int ~ "E_only",
      has_cna  & !has_expr & !has_int ~ "C_only",
      has_expr & has_cna  & !has_int ~ "E_C",
      has_expr & has_cna  & has_int  ~ "E_C_int",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(model_type))

############################################################
# Split into model-type-specific tables
############################################################

E_only_raw  <- joined %>% filter(model_type == "E_only")
C_only_raw  <- joined %>% filter(model_type == "C_only")
E_C_raw     <- joined %>% filter(model_type == "E_C")
E_C_int_raw <- joined %>% filter(model_type == "E_C_int")

############################################################
# Meta-analysis functions
############################################################

# Expression-only meta-analysis (HR_expr)
meta_expr <- function(df) {
  if (nrow(df) == 0) return(df)
  df %>%
    mutate(
      logE1 = log(HR_expr_METABRIC),
      logE2 = log(HR_expr_TCGA),
      seE1  = se_from_ci(CI_lo_expr_METABRIC, CI_hi_expr_METABRIC),
      seE2  = se_from_ci(CI_lo_expr_TCGA,     CI_hi_expr_TCGA),
      wE1   = 1 / seE1^2,
      wE2   = 1 / seE2^2,
      logHR_E_pool = (wE1 * logE1 + wE2 * logE2) / (wE1 + wE2),
      seEp  = sqrt(1 / (wE1 + wE2)),
      HR_expr_pooled      = exp(logHR_E_pool),
      CI_lo_expr_pooled   = exp(logHR_E_pool - 1.96 * seEp),
      CI_hi_expr_pooled   = exp(logHR_E_pool + 1.96 * seEp)
    )
}

# CNA-only meta-analysis (HR_cna)
meta_cna <- function(df) {
  if (nrow(df) == 0) return(df)
  df %>%
    mutate(
      logC1 = log(HR_cna_METABRIC),
      logC2 = log(HR_cna_TCGA),
      seC1  = se_from_ci(CI_lo_cna_METABRIC, CI_hi_cna_METABRIC),
      seC2  = se_from_ci(CI_lo_cna_TCGA,     CI_hi_cna_TCGA),
      wC1   = 1 / seC1^2,
      wC2   = 1 / seC2^2,
      logHR_C_pool = (wC1 * logC1 + wC2 * logC2) / (wC1 + wC2),
      seCp  = sqrt(1 / (wC1 + wC2)),
      HR_cna_pooled      = exp(logHR_C_pool),
      CI_lo_cna_pooled   = exp(logHR_C_pool - 1.96 * seCp),
      CI_hi_cna_pooled   = exp(logHR_C_pool + 1.96 * seCp)
    )
}

# Interaction meta-analysis (HR_expr_cna)
meta_int <- function(df) {
  if (nrow(df) == 0) return(df)
  df %>%
    mutate(
      logI1 = log(HR_expr_cna_METABRIC),
      logI2 = log(HR_expr_cna_TCGA),
      seI1  = se_from_ci(CI_lo_expr_cna_METABRIC, CI_hi_expr_cna_METABRIC),
      seI2  = se_from_ci(CI_lo_expr_cna_TCGA,     CI_hi_expr_cna_TCGA),
      wI1   = 1 / seI1^2,
      wI2   = 1 / seI2^2,
      logHR_I_pool = (wI1 * logI1 + wI2 * logI2) / (wI1 + wI2),
      seIp  = sqrt(1 / (wI1 + wI2)),
      HR_expr_cna_pooled      = exp(logHR_I_pool),
      CI_lo_expr_cna_pooled   = exp(logHR_I_pool - 1.96 * seIp),
      CI_hi_expr_cna_pooled   = exp(logHR_I_pool + 1.96 * seIp)
    )
}

############################################################
# Run meta-analysis per model type
############################################################

E_only   <- E_only_raw   %>% meta_expr()
C_only   <- C_only_raw   %>% meta_cna()
E_C      <- E_C_raw      %>% meta_expr() %>% meta_cna()
E_C_int  <- E_C_int_raw  %>% meta_expr() %>% meta_cna() %>% meta_int()

############################################################
# Filters based on pooled CI_lo > 1 (numeric, pre-rounding)
############################################################

E_filter             <- E_only   %>% filter(CI_lo_expr_pooled     > 1)
C_filter             <- C_only   %>% filter(CI_lo_cna_pooled      > 1)
E_C_filter_by_E      <- E_C      %>% filter(CI_lo_expr_pooled     > 1)
E_C_filter_by_C      <- E_C      %>% filter(CI_lo_cna_pooled      > 1)
E_C_int_filter_by_E   <- E_C_int %>% filter(CI_lo_expr_pooled     > 1)
E_C_int_filter_by_C   <- E_C_int %>% filter(CI_lo_cna_pooled      > 1)
E_C_int_filter_by_int <- E_C_int %>% filter(CI_lo_expr_cna_pooled > 1)

############################################################
# Column sets (E / C / E_C / E_C_int) — obey rules
############################################################

# Always keep these core columns
safe_cols <- c(
  "gene", "model",
  "n_METABRIC", "n_TCGA",
  "events_METABRIC", "events_TCGA",
  "c_index_METABRIC", "c_index_TCGA"
)

# Expression-only sheet columns (no 'cna' or 'expr_cna' columns)
keep_E <- c(
  safe_cols,
  "HR_expr_METABRIC",  "HR_expr_TCGA",
  "CI_lo_expr_METABRIC", "CI_lo_expr_TCGA",
  "CI_hi_expr_METABRIC", "CI_hi_expr_TCGA",
  "HR_expr_pooled",
  "CI_lo_expr_pooled", "CI_hi_expr_pooled"
)

# CNA-only sheet columns (no 'expr' or 'expr_cna' columns)
keep_C <- c(
  safe_cols,
  "HR_cna_METABRIC",  "HR_cna_TCGA",
  "CI_lo_cna_METABRIC", "CI_lo_cna_TCGA",
  "CI_hi_cna_METABRIC", "CI_hi_cna_TCGA",
  "HR_cna_pooled",
  "CI_lo_cna_pooled", "CI_hi_cna_pooled"
)

# Expression + CNA (no interaction columns here)
keep_EC <- c(
  safe_cols,
  # Expression part
  "HR_expr_METABRIC",  "HR_expr_TCGA",
  "CI_lo_expr_METABRIC", "CI_lo_expr_TCGA",
  "CI_hi_expr_METABRIC", "CI_hi_expr_TCGA",
  "HR_expr_pooled",
  "CI_lo_expr_pooled", "CI_hi_expr_pooled",
  # CNA part
  "HR_cna_METABRIC",  "HR_cna_TCGA",
  "CI_lo_cna_METABRIC", "CI_lo_cna_TCGA",
  "CI_hi_cna_METABRIC", "CI_hi_cna_TCGA",
  "HR_cna_pooled",
  "CI_lo_cna_pooled", "CI_hi_cna_pooled"
)

# Expression + CNA + Interaction
keep_EC_int <- c(
  keep_EC,
  "HR_expr_cna_METABRIC",  "HR_expr_cna_TCGA",
  "CI_lo_expr_cna_METABRIC", "CI_lo_expr_cna_TCGA",
  "CI_hi_expr_cna_METABRIC", "CI_hi_expr_cna_TCGA",
  "HR_expr_cna_pooled",
  "CI_lo_expr_cna_pooled", "CI_hi_expr_cna_pooled"
)

# Apply column pruning
E_only  <- E_only  %>% dplyr::select(any_of(keep_E))
E_filter<- E_filter%>% dplyr::select(any_of(keep_E))

C_only  <- C_only  %>% dplyr::select(any_of(keep_C))
C_filter<- C_filter%>% dplyr::select(any_of(keep_C))

E_C     <- E_C     %>% dplyr::select(any_of(keep_EC))
E_C_filter_by_E <- E_C_filter_by_E %>% dplyr::select(any_of(keep_EC))
E_C_filter_by_C <- E_C_filter_by_C %>% dplyr::select(any_of(keep_EC))

E_C_int <- E_C_int %>% dplyr::select(any_of(keep_EC_int))
E_C_int_filter_by_E   <- E_C_int_filter_by_E   %>% dplyr::select(any_of(keep_EC_int))
E_C_int_filter_by_C   <- E_C_int_filter_by_C   %>% dplyr::select(any_of(keep_EC_int))
E_C_int_filter_by_int <- E_C_int_filter_by_int %>% dplyr::select(any_of(keep_EC_int))

############################################################
# Apply 3 significant digits for numeric columns
############################################################

E_only  <- fmt_all_numeric(E_only)
E_filter<- fmt_all_numeric(E_filter)

C_only  <- fmt_all_numeric(C_only)
C_filter<- fmt_all_numeric(C_filter)

E_C     <- fmt_all_numeric(E_C)
E_C_filter_by_E <- fmt_all_numeric(E_C_filter_by_E)
E_C_filter_by_C <- fmt_all_numeric(E_C_filter_by_C)

E_C_int <- fmt_all_numeric(E_C_int)
E_C_int_filter_by_E   <- fmt_all_numeric(E_C_int_filter_by_E)
E_C_int_filter_by_C   <- fmt_all_numeric(E_C_int_filter_by_C)
E_C_int_filter_by_int <- fmt_all_numeric(E_C_int_filter_by_int)

############################################################
# For filtered sheets:
# - drop rows where rounded pooled CI_lo == 1
# - sort by gene (asc), then relevant HR_pooled (desc)
############################################################

clean_filtered <- function(df, ci_col, hr_col) {
  if (nrow(df) == 0) return(df)
  df %>%
    # Remove rows whose (rounded) pooled CI_lo is exactly 1
    filter(.data[[ci_col]] != 1) %>%
    arrange(gene, desc(.data[[hr_col]]))
}

# E_filter: expression pooled CI and HR
E_filter             <- clean_filtered(E_filter,
                                       ci_col = "CI_lo_expr_pooled",
                                       hr_col = "HR_expr_pooled")

# C_filter: CNA pooled CI and HR
C_filter             <- clean_filtered(C_filter,
                                       ci_col = "CI_lo_cna_pooled",
                                       hr_col = "HR_cna_pooled")

# E_C filters
E_C_filter_by_E      <- clean_filtered(E_C_filter_by_E,
                                       ci_col = "CI_lo_expr_pooled",
                                       hr_col = "HR_expr_pooled")

E_C_filter_by_C      <- clean_filtered(E_C_filter_by_C,
                                       ci_col = "CI_lo_cna_pooled",
                                       hr_col = "HR_cna_pooled")

# E_C_int filters
E_C_int_filter_by_E   <- clean_filtered(E_C_int_filter_by_E,
                                        ci_col = "CI_lo_expr_pooled",
                                        hr_col = "HR_expr_pooled")

E_C_int_filter_by_C   <- clean_filtered(E_C_int_filter_by_C,
                                        ci_col = "CI_lo_cna_pooled",
                                        hr_col = "HR_cna_pooled")

E_C_int_filter_by_int <- clean_filtered(E_C_int_filter_by_int,
                                        ci_col = "CI_lo_expr_cna_pooled",
                                        hr_col = "HR_expr_cna_pooled")

############################################################
# Write Excel workbook
############################################################

wb <- createWorkbook()

addWorksheet(wb, "E_only")
writeData(wb, "E_only", E_only)

addWorksheet(wb, "E_filter")
writeData(wb, "E_filter", E_filter)

addWorksheet(wb, "C_only")
writeData(wb, "C_only", C_only)

addWorksheet(wb, "C_filter")
writeData(wb, "C_filter", C_filter)

addWorksheet(wb, "E_C")
writeData(wb, "E_C", E_C)

addWorksheet(wb, "E_C_filter_by_E")
writeData(wb, "E_C_filter_by_E", E_C_filter_by_E)

addWorksheet(wb, "E_C_filter_by_C")
writeData(wb, "E_C_filter_by_C", E_C_filter_by_C)

addWorksheet(wb, "E_C_int")
writeData(wb, "E_C_int", E_C_int)

addWorksheet(wb, "E_C_int_filter_by_E")
writeData(wb, "E_C_int_filter_by_E", E_C_int_filter_by_E)

addWorksheet(wb, "E_C_int_filter_by_C")
writeData(wb, "E_C_int_filter_by_C", E_C_int_filter_by_C)

addWorksheet(wb, "E_C_int_filter_by_int")
writeData(wb, "E_C_int_filter_by_int", E_C_int_filter_by_int)

saveWorkbook(wb, out_file, overwrite = TRUE)

cat("Meta-analysis finished. File written to:\n", out_file, "\n")