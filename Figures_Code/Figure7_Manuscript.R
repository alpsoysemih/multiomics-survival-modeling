# Load required libraries quietly and suppress startup messages
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(survival); library(broom)
  library(tibble); library(ggplot2)
  library(org.Hs.eg.db)
  library(forcats)
  library(cowplot)
  library(scales)
  library(survRM2)
  library(timeROC)
  library(purrr)
})

# Define a null-coalescing operator for safe default values
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ===================================================================
## 0) PATHS & INPUT FILES
## ===================================================================

# Define dataset and output directories
metabric_dir <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
gdsc_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
out_dir      <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"

# Create output directory if it does not exist and report the path
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(">>> Output directory:", out_dir, "\n")

# Define input file paths for GDSC CNA and drug response
gdsc_cna_file  <- file.path(gdsc_dir, "GDSC_CNA.Paclitaxel.tsv")
gdsc_resp_file <- file.path(gdsc_dir, "GDSC_response.Paclitaxel.tsv")

# Define input file paths for METABRIC expression, CNA, and clinical tables
metabric_expr_file <- file.path(
  metabric_dir,
  "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt"
)
metabric_cna_file  <- file.path(
  metabric_dir,
  "data_cna.txt"
)
metabric_survival_file <- file.path(
  metabric_dir,
  "brca_metabric_clinical_data.tsv"
)

# Define input file paths for TCGA-BRCA expression, CNA, and survival/clinical covariates
tcga_expr_file <- file.path(
  tcga_dir,
  "TCGA-BRCA_exprs.z.tsv"
)
tcga_cna_file  <- file.path(
  tcga_dir,
  "TCGA-BRCA_CNA.tsv"
)
tcga_surv_file <- file.path(
  tcga_dir,
  "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv"
)

# Define signature genes (one expression gene + three CNA genes)
sig_expr <- "SERPINE1"
sig_cna  <- c("CLDN1","LAMB2","LYPD6B")
sig_all  <- c(sig_expr, sig_cna)

## ===================================================================
## 1) GDSC CNA direction via Fisher
## ===================================================================

# Load GDSC CNA and response labels for CNA-direction inference
cat(">>> Loading GDSC ...\n")

gdsc_cna  <- read_tsv(gdsc_cna_file,  show_col_types = FALSE)
gdsc_resp <- read_tsv(gdsc_resp_file, show_col_types = FALSE)

# Map Entrez gene IDs to HGNC gene symbols for interpretability
gdsc_cna$Gene_Symbol <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(gdsc_cna$gene_id),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Reorder columns to keep Gene_Symbol first and drop the original Entrez column
gdsc_cna <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, everything(), -gene_id)

# Binarize CNA calls (0 = no event, 1 = any non-zero CNA event)
gdsc_cna[, -1] <- ifelse(gdsc_cna[, -1] != 0, 1, 0)

# Build a sample-level response table for downstream joining
resp <- gdsc_resp %>%
  transmute(
    sample = as.character(sample_name),
    status = response
  )

# Convert CNA matrix to long format and merge with drug response status
cna_long <- gdsc_cna %>%
  pivot_longer(-Gene_Symbol, names_to = "sample", values_to = "CNA") %>%
  inner_join(resp, by = "sample")

# Compute per-gene Fisher exact test (CNA vs response) to obtain OR and p-value
fisher_or <- cna_long %>%
  group_by(Gene_Symbol) %>%
  summarise({
    tab <- table(CNA, status)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      tibble(OR = NA_real_, p = NA_real_)
    } else {
      ft <- fisher.test(tab)
      tibble(OR = unname(ft$estimate), p = ft$p.value)
    }
  }) %>%
  ungroup()

## ===================================================================
## 2) METABRIC
## ===================================================================

# Load METABRIC expression, CNA, and clinical covariates
cat(">>> Loading METABRIC ...\n")

met_expr <- read_tsv(metabric_expr_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol) %>%
  dplyr::select(-Entrez_Gene_Id)

met_cna <- read_tsv(metabric_cna_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol) %>%
  dplyr::select(-Entrez_Gene_Id)

# Binarize METABRIC CNA calls (0 = no event, 1 = any non-zero CNA event)
met_cna_bin <- met_cna
met_cna_bin[, -1] <- ifelse(met_cna_bin[, -1] != 0, 1, 0)

# Parse clinical variables and harmonize OS, age, stage, and PAM50 labels
metabric_survival <- read_tsv(metabric_survival_file, show_col_types = FALSE) 
  
metabric_surv_small <- metabric_survival %>%
  dplyr::rename(
    sample      = `Sample ID`,
    Age         = `Age at Diagnosis`,
    Tumor_Stage = `Tumor Stage`,
    OS_TIME_m   = `Overall Survival (Months)`,
    OS_Status   = `Overall Survival Status`,
    PAM50_raw   = `Pam50 + Claudin-low subtype`
  ) %>%
  dplyr::mutate(
    OS_TIME  = as.numeric(OS_TIME_m),
    OS_EVENT = ifelse(grepl("^1", OS_Status), 1, 0),
    Age      = as.numeric(Age),
    Stage_simple = dplyr::case_when(
      Tumor_Stage %in% c("0", "Stage 0", "Stage0")        ~ NA_character_,
      Tumor_Stage %in% c("1","I","Stage I")               ~ "StageI",
      Tumor_Stage %in% c("2","II","Stage II")             ~ "StageII",
      Tumor_Stage %in% c("3","III","Stage III",
                         "4","IV","Stage IV")             ~ "StageIII_IV",
      TRUE                                                ~ "StageX"
    ),
    PAM50 = {
      f <- factor(PAM50_raw)
      # drop Claudin-low (and any spelling variants if present)
      f <- droplevels(f[!(f %in% c("Claudin-low", "Claudin low", "Claudin_low"))])
      if ("LumA" %in% levels(f)) stats::relevel(f, ref = "LumA") else f
    },
    Stage_simple = {
      f <- factor(
        Stage_simple,
        levels = c("StageI", "StageII", "StageIII_IV", "StageX")
      )
      if ("StageI" %in% levels(f)) stats::relevel(f, ref = "StageI") else f
    }
  ) %>%
  dplyr::filter(
    !is.na(OS_TIME),
    !is.na(OS_EVENT)
  ) %>%
  dplyr::select(sample, OS_TIME, OS_EVENT, Age, Stage_simple, PAM50)

# Convert METABRIC expression matrix to long format for the signature genes
expr_long <- met_expr %>%
  filter(Gene_Symbol %in% sig_all) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

# Convert METABRIC CNA to long format and align CNA direction using GDSC-derived OR
cna_long_met <- met_cna_bin %>%
  filter(Gene_Symbol %in% sig_cna) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  ) %>%
  left_join(fisher_or, by = "Gene_Symbol") %>%
  mutate(
    dir         = case_when(
      OR > 1  ~  1,
      OR < 1  ~ -1,
      TRUE    ~ NA_real_
    ),
    cna_aligned = ifelse(dir == 1, cna, 1 - cna)
  )

# Extract per-sample signature variables (expression + aligned CNA)
serp  <- expr_long    %>% filter(Gene_Symbol == sig_expr) %>%
  transmute(sample, SERPINE1_expr = expr)
cl1   <- cna_long_met %>% filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = cna_aligned)
lamb2 <- cna_long_met %>% filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = cna_aligned)
ly6   <- cna_long_met %>% filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = cna_aligned)

# Merge clinical covariates with signature variables and keep complete cases only
met_sig <- metabric_surv_small %>%
  left_join(serp,  by = "sample") %>%
  left_join(cl1,   by = "sample") %>%
  left_join(lamb2, by = "sample") %>%
  left_join(ly6,   by = "sample") %>%
  drop_na()

# Report METABRIC sample size after merging and filtering
cat(">>> METABRIC signature N:", nrow(met_sig), "\n")

## ===================================================================
## 3) TCGA-BRCA
## ===================================================================

# Load TCGA-BRCA expression, CNA, and survival/clinical covariates
cat(">>> Loading TCGA-BRCA ...\n")

tcga_expr     <- read_tsv(tcga_expr_file, show_col_types = FALSE)
tcga_cna      <- read_tsv(tcga_cna_file,  show_col_types = FALSE)
tcga_survival <- read_tsv(tcga_surv_file, show_col_types = FALSE)

# Map Entrez gene IDs to HGNC gene symbols for both expression and CNA tables
tcga_expr$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_expr$ENTREZID),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

tcga_cna$gene <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_cna$gene),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Standardize gene identifier column names across TCGA tables
colnames(tcga_expr)[colnames(tcga_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(tcga_cna)[colnames(tcga_cna)   == "gene"]     <- "Gene_Symbol"

# Binarize TCGA CNA calls (0 = no event, 1 = any non-zero CNA event)
tcga_cna_bin <- tcga_cna
tcga_cna_bin[, -1] <- ifelse(tcga_cna_bin[, -1] != 0, 1, 0)

# Build a harmonized TCGA clinical table (OS, age, stage, PAM50)
tcga_clin_small <- tcga_survival %>%
  dplyr::transmute(
    sample,
    OS_TIME_d = as.numeric(OS.time),  # days
    OS_EVENT  = as.numeric(OS),       # 0 / 1
    Age       = as.numeric(Age),
    Stage     = Stage,
    PAM50_raw = BRCA_Subtype_PAM50
  ) %>%
  mutate(
    ## OS time in months
    OS_TIME  = OS_TIME_d / 30.44,
    OS_EVENT = OS_EVENT,
    Stage_simple = case_when(
      Stage %in% c("Stage I", "Stage IA", "Stage IB")                                 ~ "StageI",
      Stage %in% c("Stage II", "Stage IIA", "Stage IIB")                              ~ "StageII",
      Stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV") ~ "StageIII_IV",
      TRUE                                                                            ~ "StageX"
    ),
    PAM50 = {
      f <- factor(PAM50_raw)
      if ("LumA" %in% levels(f)) stats::relevel(f, ref = "LumA") else f
    },
    Stage_simple = {
      f <- factor(Stage_simple, levels = c("StageI", "StageII", "StageIII_IV", "StageX"))
      if ("StageI" %in% levels(f)) stats::relevel(f, ref = "StageI") else f
    }
  ) %>%
  filter(
    !is.na(OS_TIME),
    !is.na(OS_EVENT)
  ) %>%
  dplyr::select(sample, OS_TIME, OS_EVENT, Age, Stage_simple, PAM50)

# Drop extremely rare stage levels to improve model stability
stage_counts_tcga <- table(tcga_clin_small$Stage_simple)
rare_stages_tcga  <- names(stage_counts_tcga[stage_counts_tcga < 5])
if (length(rare_stages_tcga) > 0) {
  tcga_clin_small <- tcga_clin_small %>%
    filter(!Stage_simple %in% rare_stages_tcga) %>%
    mutate(Stage_simple = droplevels(Stage_simple))
}

# Report TCGA sample size after clinical filtering
cat(">> TCGA-BRCA OS sample size:", nrow(tcga_clin_small), "\n")

# Convert TCGA expression matrix to long format and harmonize barcodes to sample level
expr_tcga <- tcga_expr %>%
  filter(Gene_Symbol %in% sig_all) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "bc",
    values_to = "expr"
  ) %>%
  mutate(sample = substr(bc, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, expr)

# Convert TCGA CNA matrix to long format and align CNA direction using GDSC-derived OR
cna_tcga <- tcga_cna_bin %>%
  filter(Gene_Symbol %in% sig_cna) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "bc",
    values_to = "cna"
  ) %>%
  mutate(sample = substr(bc, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, cna) %>%
  left_join(fisher_or, by = "Gene_Symbol") %>%
  mutate(
    dir         = case_when(
      OR > 1  ~  1,
      OR < 1  ~ -1,
      TRUE    ~ NA_real_
    ),
    cna_aligned = ifelse(dir == 1, cna, 1 - cna)
  )

# Extract per-sample signature variables for TCGA (expression + aligned CNA)
serp_tcga <- expr_tcga %>%
  filter(Gene_Symbol == sig_expr) %>%
  transmute(sample, SERPINE1_expr = expr)

cl1_tcga  <- cna_tcga %>%
  filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = cna_aligned)

lamb2_tcga <- cna_tcga %>%
  filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = cna_aligned)

ly6_tcga  <- cna_tcga %>%
  filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = cna_aligned)

# Merge clinical covariates with signature variables and keep complete cases only
tcga_sig <- tcga_clin_small %>%
  left_join(serp_tcga,  by = "sample") %>%
  left_join(cl1_tcga,   by = "sample") %>%
  left_join(lamb2_tcga, by = "sample") %>%
  left_join(ly6_tcga,   by = "sample") %>%
  drop_na()

# Report TCGA sample size after merging and filtering
cat(">>> TCGA-BRCA signature N:", nrow(tcga_sig), "\n")

## ===================================================================
## 4) Train signature-only model on METABRIC and derive risk score
## ===================================================================

# Fit a Cox model on METABRIC using signature variables and compute risk scores
cat(">>> Fitting signature-only model on METABRIC ...\n")

cox_sig_train <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LAMB2_CNA + LYPD6B_CNA,
  data = met_sig
)

# Predict linear predictors (risk scores) in training and external cohort
met_sig$score  <- predict(cox_sig_train, newdata = met_sig,  type = "lp")
tcga_sig$score <- predict(cox_sig_train, newdata = tcga_sig, type = "lp")

# Define median cutoff using training cohort scores and assign risk groups
cut_median <- median(met_sig$score, na.rm = TRUE)
cat(">>> METABRIC signature score median:", cut_median, "\n")

met_sig <- met_sig %>%
  mutate(
    risk_group = ifelse(score >= cut_median, "High", "Low"),
    risk_group = factor(risk_group, levels = c("Low","High"))
  )

tcga_sig <- tcga_sig %>%
  mutate(
    risk_group = ifelse(score >= cut_median, "High", "Low"),
    risk_group = factor(risk_group, levels = c("Low","High"))
  )

## ===================================================================
## 5) Forest plots (Figure 6a–b)
## ===================================================================

# Fit Cox models for forest plots (score-only and score+clinical) and format HR tables
cat(">>> Building forest plots ...\n")

cox_met_score_only  <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ score,
  data = met_sig
)

cox_tcga_score_only <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ score,
  data = tcga_sig
)

# Helper: extract score HR and confidence interval from a Cox model
tidy_score_only <- function(cox_fit, cohort_label) {
  broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(term == "score") %>%
    dplyr::mutate(
      HR         = estimate,
      HR_low     = conf.low,
      HR_high    = conf.high,
      term_clean = "Risk score",
      Cohort     = cohort_label
    )
}

# Build forest input table for score-only models across cohorts (Figure 6a)
df_forest_a <- dplyr::bind_rows(
  tidy_score_only(cox_met_score_only,  "METABRIC"),
  tidy_score_only(cox_tcga_score_only, "TCGA-BRCA")
) %>%
  dplyr::mutate(
    term_factor = factor(term_clean, levels = "Risk score"),
    Cohort      = factor(Cohort, levels = c("METABRIC","TCGA-BRCA"))
  )

# Fit a multivariable Cox model in TCGA including risk score + clinical covariates
cox_tcga_full <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ score + Age + Stage_simple + PAM50,
  data = tcga_sig
)

# Helper: format HR table for a multivariable Cox model with readable labels
tidy_forest_full <- function(cox_fit, cohort_label) {
  broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(!is.na(estimate)) %>%
    dplyr::mutate(
      HR      = estimate,
      HR_low  = conf.low,
      HR_high = conf.high,
      term_clean = dplyr::case_when(
        term == "score" ~ "Risk score",
        term == "Age"   ~ "Age (years)",
        stringr::str_detect(term, "Stage_simple") ~
          paste0("Stage: ", stringr::str_replace(term, "Stage_simple", "")),
        stringr::str_detect(term, "PAM50") ~
          stringr::str_replace(term, "PAM50", "PAM50: "),
        TRUE ~ term
      ),
      Cohort = cohort_label
    )
}

# Build forest input table for the TCGA multivariable model (Figure 6b)
df_forest_b <- tidy_forest_full(cox_tcga_full, "TCGA-BRCA") %>%
  dplyr::mutate(
    term_factor = forcats::fct_rev(
      factor(term_clean, levels = unique(term_clean))
    )
  )

# Format x-axis tick labels with fixed decimal places
lab_hr_3dec <- function(x) formatC(x, digits = 2, format = "f")

## Figure 6a
# Create forest plot for risk score only (METABRIC and TCGA)
p_forest_a <- ggplot(df_forest_a,
                     aes(x = HR, y = term_factor)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = HR_low, xmax = HR_high),
                  size = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", HR)),
    vjust = -1.0, hjust = 0.5,
    fontface = "bold",
    size = 3.5
  ) +
  scale_x_log10(
    breaks = c(1, 3, 10, 30),
    labels = lab_hr_3dec
  ) +
  facet_wrap(~Cohort, nrow = 1) +
  theme_bw() +
  labs(
    x = "HR",
    y = "Risk score",
    title = "Risk scores without clinical covariates"
  ) +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin  = margin(t = 10, r = 5, b = 5, l = 12)
  )

## Figure 6b
# Create forest plot for TCGA multivariable model (risk score + clinical covariates)
p_forest_b <- ggplot(df_forest_b,
                     aes(x = HR, y = term_factor)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = HR_low, xmax = HR_high),
                  size = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", HR)),
    vjust = -1.0, hjust = 0.5,
    fontface = "bold",
    size = 3.5
  ) +
  scale_x_log10(
    breaks = c(1, 3, 10, 30),
    labels = lab_hr_3dec
  ) +
  theme_bw() +
  labs(
    x = "HR",
    y = "",
    title = "Risk score and clinical covariates (TCGA-BRCA)"
  ) +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(face = "bold", size = 11),
    legend.position = "none",
    plot.margin  = margin(t = 10, r = 5, b = 5, l = 12)
  )

## (Optional) Save Figure 6 here if needed
## fig6 <- plot_grid(p_forest_a, p_forest_b, ncol = 2,
##                   labels = c("a.", "b."), label_size = 16,
##                   label_fontface = "bold")
## ggsave(file.path(out_dir, "Figure_6.jpeg"), fig6, width = 12, height = 5, dpi = 300)

## ===================================================================
## Figure 7 – RMST difference & time-dependent AUC (up to 200 months)
## ===================================================================

# Build Figure 7 components: RMST difference curve and time-dependent AUC curve
cat(">>> Building Figure 7 (RMST & time-dependent AUC up to 200 months) ...\n")

## -------------------------------------------------------------------
## 1) Helper: RMST curve (High – Low) for a given cohort
## -------------------------------------------------------------------

# Compute RMST(High) - RMST(Low) across a grid of truncation times (tau)
compute_rmst_curve <- function(dat, cohort_label, t_grid, max_tau_cap = 200) {
  ## dat: data.frame with OS_TIME, OS_EVENT, risk_group (Low/High)
  ## t_grid: candidate tau values (months)
  
  t_low_max  <- max(dat$OS_TIME[dat$risk_group == "Low"],  na.rm = TRUE)
  t_high_max <- max(dat$OS_TIME[dat$risk_group == "High"], na.rm = TRUE)
  tau_max_by_data <- min(t_low_max, t_high_max)
  
  tau_max <- min(max_tau_cap, tau_max_by_data)
  t_grid_use <- t_grid[t_grid <= tau_max]
  
  cat(sprintf("   %s: using tau up to %.1f months\n", cohort_label, tau_max))
  
  # Robust extractor for RMST estimates from rmst2 outputs across object types
  get_rmst_val <- function(x) {
    if (is.list(x) && !is.data.frame(x) && !is.matrix(x)) {
      nm <- names(x)
      pref <- intersect(c("RMST", "rmst", "est", "estimate"), nm)
      if (length(pref) > 0) {
        return(get_rmst_val(x[[pref[1]]]))
      } else {
        return(get_rmst_val(x[[1]]))
      }
    }
    if (is.matrix(x)) {
      if ("RMST" %in% rownames(x) && "Est." %in% colnames(x)) {
        return(as.numeric(x["RMST", "Est."]))
      } else {
        return(as.numeric(x[1]))
      }
    }
    if (is.data.frame(x)) {
      rn <- rownames(x)
      cn <- colnames(x)
      if ("RMST" %in% rn && "Est." %in% cn) {
        return(as.numeric(x["RMST", "Est."]))
      }
      if ("term" %in% cn && "RMST" %in% x$term) {
        est_col <- intersect(c("Est.","est","estimate"), cn)
        if (length(est_col) > 0L) {
          return(as.numeric(x[x$term == "RMST", est_col[1]][1]))
        }
      }
      return(as.numeric(x[[1]][1]))
    }
    if (is.numeric(x)) {
      if (!is.null(names(x)) && "RMST" %in% names(x)) {
        return(as.numeric(x["RMST"]))
      } else {
        return(as.numeric(x[1]))
      }
    }
    as.numeric(x)[1]
  }
  
  # Loop over tau values and compute RMST difference per cohort
  purrr::map_dfr(t_grid_use, function(tt) {
    fit <- survRM2::rmst2(
      time   = dat$OS_TIME,
      status = dat$OS_EVENT,
      arm    = ifelse(dat$risk_group == "High", 1, 0),
      tau    = tt
    )
    
    rmst_high <- get_rmst_val(fit$RMST.arm1)
    rmst_low  <- get_rmst_val(fit$RMST.arm0)
    
    tibble::tibble(
      Cohort    = cohort_label,
      tau       = tt,
      rmst_diff = as.numeric(rmst_high - rmst_low)
    )
  })
}

## -------------------------------------------------------------------
## 2) RMST curves (up to 200 months)
## -------------------------------------------------------------------

# Define truncation-time grid (tau) and compute RMST difference curves per cohort
t_grid_rmst <- seq(12, 200, by = 12)

cat("   Computing RMST curves ...\n")
rmst_met  <- compute_rmst_curve(met_sig,  "METABRIC",  t_grid_rmst)
rmst_tcga <- compute_rmst_curve(tcga_sig, "TCGA-BRCA", t_grid_rmst)

# Combine RMST results across cohorts for plotting
rmst_df <- bind_rows(rmst_met, rmst_tcga)

# Plot RMST difference curves (High - Low) up to 200 months
p_rmst <- ggplot(rmst_df,
                 aes(x = tau, y = rmst_diff,
                     colour = Cohort, group = Cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(
    limits = c(0, 200),
    breaks = seq(0, 200, by = 50)
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(min(rmst_df$rmst_diff) / 3) * 3,
      ceiling(max(rmst_df$rmst_diff) / 3) * 3,
      by = 3
    )
  ) +
  scale_colour_manual(values = c("METABRIC" = "#E41A1C",
                                 "TCGA-BRCA" = "#377EB8")) +
  theme_bw() +
  labs(
    x = "Time (months)",
    y = "RMST difference",
    title = "RMST difference"
  ) +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(face = "bold", size = 11)
  )

## -------------------------------------------------------------------
## 3) Time-dependent AUC of risk score (up to 200 months)
## -------------------------------------------------------------------

# Compute time-dependent AUC of the continuous risk score up to a common follow-up horizon
cat("   Computing time-dependent AUC curves ...\n")

t_max_met      <- min(200, max(met_sig$OS_TIME,  na.rm = TRUE))
t_max_tcga     <- min(200, max(tcga_sig$OS_TIME, na.rm = TRUE))
t_max_common   <- floor(min(t_max_met, t_max_tcga))
t_grid_auc     <- seq(12, t_max_common, by = 12)

cat(sprintf("   Time-dependent AUC evaluated up to %d months\n", t_max_common))

# Estimate time-dependent ROC/AUC for METABRIC
time_roc_met <- timeROC::timeROC(
  T      = met_sig$OS_TIME,
  delta  = met_sig$OS_EVENT,
  marker = met_sig$score,
  cause  = 1,
  times  = t_grid_auc,
  iid    = FALSE
)

# Estimate time-dependent ROC/AUC for TCGA-BRCA
time_roc_tcga <- timeROC::timeROC(
  T      = tcga_sig$OS_TIME,
  delta  = tcga_sig$OS_EVENT,
  marker = tcga_sig$score,
  cause  = 1,
  times  = t_grid_auc,
  iid    = FALSE
)

# Build a tidy table of AUC values for plotting
auc_df <- tibble(
  time   = rep(t_grid_auc, 2),
  AUC    = c(time_roc_met$AUC, time_roc_tcga$AUC),
  Cohort = rep(c("METABRIC", "TCGA-BRCA"),
               each = length(t_grid_auc))
)

# Plot time-dependent AUC curves up to 200 months
p_auc <- ggplot(auc_df,
                aes(x = time, y = AUC,
                    colour = Cohort, group = Cohort)) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(
    limits = c(0, 200),
    breaks = seq(0, 200, by = 50)
  ) +
  scale_y_continuous(
    limits = c(0.40, 1.00),
    labels = number_format(accuracy = 0.01),
    breaks = seq(0.40, 1.00, by = 0.10),
    expand = c(0, 0)
  ) +
  scale_colour_manual(values = c("METABRIC" = "#E41A1C",
                                 "TCGA-BRCA" = "#377EB8")) +
  theme_bw() +
  labs(
    x = "Time (months)",
    y = "Time-dependent AUC",
    title = "Time-dependent AUC of risk score"
  ) +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold", size = 11),
    legend.text  = element_text(face = "bold", size = 11)
  )

## -------------------------------------------------------------------
## 4) Final Figure 7 (a: RMST, b: AUC)
## -------------------------------------------------------------------

# Assemble Figure 7 panels (a: RMST difference, b: time-dependent AUC) and save to disk
fig7 <- cowplot::plot_grid(
  p_rmst, p_auc,
  ncol           = 2,
  rel_widths     = c(1, 1),
  labels         = c("A.", "B."),
  label_size     = 16,
  label_fontface = "bold",
  hjust          = -0.1
)

ggsave(
  file.path(out_dir, "Figure 7.jpeg"),
  fig7,
  width  = 14,
  height = 6,
  dpi    = 300
)

# Print completion messages
cat(">>> Figure 7 saved (up to 200 months).\n")
cat(">>> Script finished successfully.\n")