# Load required libraries quietly (I/O, wrangling, survival, plotting, annotation)
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(survival); library(broom)
  library(tibble); library(ggplot2)
  library(org.Hs.eg.db)
  library(forcats)
  library(cowplot)
})

# Null-coalescing helper: return a if not NULL, otherwise b
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ===================================================================
## 0) PATHS
## ===================================================================

# Define dataset and output directories
metabric_dir <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
gdsc_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
out_dir      <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"

# Create output directory if it does not exist
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Define signature genes (one expression gene + three CNA genes)
sig_expr <- "SERPINE1"
sig_cna  <- c("CLDN1","LYPD6B","LAMB2")
sig_all  <- c(sig_expr, sig_cna)

# Print output directory for logging
cat(">>> Output directory:", out_dir, "\n")

## ===================================================================
## 1) GDSC CNA direction via Fisher
## ===================================================================

# Announce GDSC processing
cat(">>> Loading GDSC ...\n")

# Define GDSC CNA and response file paths
gdsc_cna_file  <- file.path(gdsc_dir,"GDSC_CNA.Paclitaxel.tsv")
gdsc_resp_file <- file.path(gdsc_dir,"GDSC_response.Paclitaxel.tsv")

# Load GDSC CNA matrix and response labels
gdsc_cna  <- read_tsv(gdsc_cna_file, show_col_types = FALSE)
gdsc_resp <- read_tsv(gdsc_resp_file, show_col_types = FALSE)

## Map Entrez → Symbol
# Convert Entrez gene IDs to gene symbols for interpretability
gdsc_cna$Gene_Symbol <- mapIds(
  org.Hs.eg.db,
  as.character(gdsc_cna$gene_id),
  "SYMBOL","ENTREZID"
)
# Move Gene_Symbol to the first column and drop the original gene_id column
gdsc_cna <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, everything(), -gene_id)

## Binary CNA
# Binarize CNA calls (0 = no event, 1 = any non-zero CNA event)
gdsc_cna[,-1] <- ifelse(gdsc_cna[,-1] != 0, 1, 0)

# Build sample-level response table for joining (sample name + response status)
resp <- gdsc_resp %>%
  transmute(
    sample = as.character(sample_name),
    status = response
  )

## Long CNA for Fisher
# Convert CNA matrix to long format and merge with response status
cna_long <- gdsc_cna %>%
  pivot_longer(-Gene_Symbol, names_to = "sample", values_to = "CNA") %>%
  inner_join(resp, by = "sample")

# Compute per-gene Fisher exact test (CNA vs response) to get OR and p-value
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

# Announce METABRIC processing
cat(">>> Loading METABRIC ...\n")

# Define METABRIC expression, CNA, and clinical file paths
expr_file <- file.path(metabric_dir,"data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
cna_file  <- file.path(metabric_dir,"data_cna.txt")
clin_file <- file.path(metabric_dir,"brca_metabric_clinical_data.tsv")

# Load METABRIC expression and standardize gene identifier column name
met_expr <- read_tsv(expr_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = `Hugo_Symbol`) %>%
  dplyr::select(-Entrez_Gene_Id)

# Load METABRIC CNA and standardize gene identifier column name
met_cna <- read_tsv(cna_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = `Hugo_Symbol`) %>%
  dplyr::select(-Entrez_Gene_Id)

## Binary CNA
# Binarize METABRIC CNA calls (0 = no event, 1 = any non-zero CNA event)
met_cna_bin <- met_cna
met_cna_bin[,-1] <- ifelse(met_cna_bin[,-1] != 0, 1, 0)

## Clinical cleanup (harmonized with TCGA-style Stage_simple + PAM50)
# Load and harmonize clinical variables (OS time/event, age, stage, PAM50)
clin <- read_tsv(clin_file, show_col_types = FALSE) %>%
  dplyr::rename(
    sample    = `Sample ID`,
    Age       = `Age at Diagnosis`,
    Stage_raw = `Tumor Stage`,
    OS_TM     = `Overall Survival (Months)`,
    OS_ST     = `Overall Survival Status`,
    PAM50_raw = `Pam50 + Claudin-low subtype`
  ) %>%
  mutate(
    OS_TIME  = as.numeric(OS_TM),
    OS_EVENT = ifelse(str_detect(OS_ST, "^1"), 1, 0),
    Age      = as.numeric(Age),
    Stage_simple = dplyr::case_when(
      Stage_raw %in% c("0", "Stage 0", "Stage0")        ~ NA_character_,
      Stage_raw %in% c("1","I","Stage I")               ~ "StageI",
      Stage_raw %in% c("2","II","Stage II")             ~ "StageII",
      Stage_raw %in% c("3","III","Stage III",
                       "4","IV","Stage IV")             ~ "StageIII_IV",
      TRUE                                                ~ "StageX"
    ),
    PAM50 = {
      f <- factor(PAM50_raw)
      # drop Claudin-low (and any spelling variants if present)
      f <- droplevels(f[!(f %in% c("Claudin-low", "Claudin low", "Claudin_low"))])
      if ("LumA" %in% levels(f)) stats::relevel(f, ref = "LumA") else f
    }
  ) %>%
  filter(!is.na(OS_TIME)) %>%
  mutate(
    Stage_simple = {
      f <- factor(
        Stage_simple,
        levels = c("StageI", "StageII", "StageIII_IV", "StageX")
      )
      if ("StageI" %in% levels(f)) stats::relevel(f, ref = "StageI") else f
    }
  ) %>%
  dplyr::select(sample, OS_TIME, OS_EVENT, Age, Stage_simple, PAM50) %>%
  droplevels()

## Expression long
# Convert expression matrix to long format for signature genes only
expr_long <- met_expr %>%
  filter(Gene_Symbol %in% sig_all) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

## CNA aligned using GDSC Fisher OR
# Convert CNA matrix to long format for CNA signature genes and align direction using GDSC OR
cna_long_met <- met_cna_bin %>%
  filter(Gene_Symbol %in% sig_cna) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  ) %>%
  left_join(fisher_or, by = "Gene_Symbol") %>%
  mutate(
    dir        = case_when(OR > 1 ~  1,
                           OR < 1 ~ -1,
                           TRUE   ~ NA_real_),
    cna_aligned = ifelse(dir == 1, cna, 1 - cna)
  )

## Join signature variables
# Extract SERPINE1 expression per sample
serp  <- expr_long    %>% filter(Gene_Symbol == sig_expr) %>%
  transmute(sample, SERPINE1_expr = expr)
# Extract aligned CNA indicators per sample (one column per gene)
cl1   <- cna_long_met %>% filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = cna_aligned)
lamb2 <- cna_long_met %>% filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = cna_aligned)
ly6   <- cna_long_met %>% filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = cna_aligned)

# Merge clinical covariates with signature variables and drop incomplete cases
met_sig <- clin %>%
  left_join(serp,  by = "sample") %>%
  left_join(cl1,   by = "sample") %>%
  left_join(ly6,   by = "sample") %>%
  left_join(lamb2, by = "sample") %>%
  drop_na()

# Report METABRIC usable sample size after merging
cat(">>> METABRIC signature N:", nrow(met_sig), "\n")

## ===================================================================
## 3) TCGA-BRCA
## ===================================================================

# Announce TCGA processing
cat(">>> Loading TCGA-BRCA ...\n")

# Load TCGA expression, CNA, and survival/clinical data
tcga_expr     <- read_tsv(file.path(tcga_dir,"TCGA-BRCA_exprs.z.tsv"), show_col_types = FALSE)
tcga_cna      <- read_tsv(file.path(tcga_dir,"TCGA-BRCA_CNA.tsv"),      show_col_types = FALSE)
tcga_survival <- read_tsv(file.path(tcga_dir,"TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv"),
                          show_col_types = FALSE)

## Map genes to SYMBOL
# Map expression Entrez IDs to gene symbols (in-place replacement)
tcga_expr$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_expr$ENTREZID),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Map CNA Entrez IDs to gene symbols (in-place replacement)
tcga_cna$gene <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_cna$gene),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Standardize gene column names across TCGA tables
colnames(tcga_expr)[colnames(tcga_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(tcga_cna)[colnames(tcga_cna)   == "gene"]     <- "Gene_Symbol"

# Binarize TCGA CNA calls (0 = no event, 1 = any non-zero CNA event)
tcga_cna_bin <- tcga_cna
tcga_cna_bin[, -1] <- ifelse(tcga_cna_bin[, -1] != 0, 1, 0)

## Clinical cleanup
# Select and harmonize TCGA clinical variables (OS, age, stage, PAM50)
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

## Drop ultra-rare stages (<5 samples)
# Identify stage levels with very small sample size and remove them to stabilize modeling
stage_counts_tcga <- table(tcga_clin_small$Stage_simple)
rare_stages_tcga  <- names(stage_counts_tcga[stage_counts_tcga < 5])
if (length(rare_stages_tcga) > 0) {
  tcga_clin_small <- tcga_clin_small %>%
    filter(!Stage_simple %in% rare_stages_tcga) %>%
    mutate(Stage_simple = droplevels(Stage_simple))
}

# Report TCGA clinical sample size after filtering
cat(">> TCGA-BRCA OS sample size:", nrow(tcga_clin_small), "\n")

## Expression long (TCGA)
# Convert TCGA expression matrix to long format for signature genes and standardize barcodes
expr_tcga <- tcga_expr %>%
  filter(Gene_Symbol %in% sig_all) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "bc",
    values_to = "expr"
  ) %>%
  mutate(sample = substr(bc, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, expr)

## CNA aligned using GDSC Fisher OR (TCGA)
# Convert TCGA CNA matrix to long format and align CNA direction using GDSC OR
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
    dir        = case_when(OR > 1 ~  1,
                           OR < 1 ~ -1,
                           TRUE   ~ NA_real_),
    cna_aligned = ifelse(dir == 1, cna, 1 - cna)
  )

# Extract SERPINE1 expression per TCGA sample
serp_tcga <- expr_tcga %>%
  filter(Gene_Symbol == sig_expr) %>%
  transmute(sample, SERPINE1_expr = expr)

# Extract aligned CNA indicators per TCGA sample (one column per gene)
cl1_tcga  <- cna_tcga %>%
  filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = cna_aligned)

lamb2_tcga <- cna_tcga %>%
  filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = cna_aligned)

ly6_tcga  <- cna_tcga %>%
  filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = cna_aligned)

# Merge TCGA clinical covariates with signature variables and drop incomplete cases
tcga_sig <- tcga_clin_small %>%
  left_join(serp_tcga,  by = "sample") %>%
  left_join(cl1_tcga,   by = "sample") %>%
  left_join(ly6_tcga,   by = "sample") %>%
  left_join(lamb2_tcga, by = "sample") %>%
  drop_na()

# Report TCGA usable sample size after merging
cat(">>> TCGA-BRCA signature N:", nrow(tcga_sig), "\n")

## ===================================================================
## 4) Train signature-only model on METABRIC and derive risk score
## ===================================================================

# Announce training step for the signature-only Cox model
cat(">>> Fitting signature-only model on METABRIC ...\n")

# Fit signature-only Cox model on METABRIC (training cohort)
cox_sig_train <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LYPD6B_CNA + LAMB2_CNA,
  data = met_sig
)

met_sig$score  <- predict(cox_sig_train, newdata = met_sig,  type = "lp")
tcga_sig$score <- predict(cox_sig_train, newdata = tcga_sig, type = "lp")

# Compute linear predictor (risk score) for METABRIC and project to TCGA
met_sig$score  <- predict(cox_sig_train, newdata = met_sig,  type = "lp")
tcga_sig$score <- predict(cox_sig_train, newdata = tcga_sig, type = "lp")

# Define median cutoff using training cohort scores
cut_median <- median(met_sig$score, na.rm = TRUE)
cat(">>> METABRIC signature score median:", cut_median, "\n")

# Assign Low/High risk groups in METABRIC based on training median cutoff
met_sig <- met_sig %>%
  mutate(
    risk_group = ifelse(score >= cut_median, "High", "Low"),
    risk_group = factor(risk_group, levels = c("Low","High"))
  )

# Assign Low/High risk groups in TCGA using the same training median cutoff
tcga_sig <- tcga_sig %>%
  mutate(
    risk_group = ifelse(score >= cut_median, "High", "Low"),
    risk_group = factor(risk_group, levels = c("Low","High"))
  )

## ===================================================================
## 5) Forest plots
## ===================================================================

# Announce forest plot construction
cat(">>> Building forest plots ...\n")

# Fit score-only Cox model on METABRIC (HR of risk score)
cox_met_score_only  <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ score,
  data = met_sig
)

# Fit score-only Cox model on TCGA (HR of risk score)
cox_tcga_score_only <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ score,
  data = tcga_sig
)

# Helper: extract HR/CI for the risk score term and add cohort label
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

# Build forest table for score-only models (METABRIC + TCGA)
df_forest_a <- dplyr::bind_rows(
  tidy_score_only(cox_met_score_only,  "METABRIC"),
  tidy_score_only(cox_tcga_score_only, "TCGA-BRCA")
) %>%
  dplyr::mutate(
    term_factor = factor(term_clean, levels = "Risk score"),
    Cohort      = factor(Cohort, levels = c("METABRIC","TCGA-BRCA"))
  )

# Fit full Cox model on TCGA including score + clinical covariates
cox_tcga_full <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ score + Age + Stage_simple + PAM50,
  data = tcga_sig
)

# Helper: tidy full Cox model and generate readable term labels
tidy_forest_full <- function(cox_fit, cohort_label) {
  broom::tidy(cox_fit, exponentiate = TRUE, conf.int = TRUE) %>%
    dplyr::filter(!is.na(estimate)) %>%
    dplyr::mutate(
      HR      = estimate,
      HR_low  = conf.low,
      HR_high = conf.high,
      
      ## ---- LABEL CLEANUP ----
      term_clean = dplyr::case_when(
        term == "score" ~ "Risk score",
        term == "Age"   ~ "Age (years)",
        
        ## ---- STAGE LABELS ----
        stringr::str_detect(term, "Stage_simple") ~
          paste0("Stage: ", stringr::str_replace(term, "Stage_simple", "")),
        
        ## ---- PAM50 LABELS ----
        stringr::str_detect(term, "PAM50") ~
          stringr::str_replace(term, "PAM50", "PAM50: "),
        
        TRUE ~ term
      ),
      Cohort = cohort_label
    )
}

## ---- Desired labels for plot terms (used only for plotting) ----
term_labels_fig6 <- c(
  age                     = "Age",
  stage_simpleStageII     = "Stage II",
  stage_simpleStageIII_IV = "Stage III/IV",
  
  pam50LumB               = "PAM50 (LumB)",
  pam50Basal              = "PAM50 (Basal)",
  pam50Her2               = "PAM50 (HER2)",
  pam50Normal             = "PAM50 (Normal)",
  
  risk_score              = "Gene signature score",
  ras_cna_score           = "Ras signaling pathway\n(CNA score)"
)

## ---- Plot-only relabeling for Tag B (do NOT change model fitting) ----
df_forest_b <- tidy_forest_full(cox_tcga_full, "TCGA-BRCA") %>%
  # Remove StageX ONLY from the plot
  dplyr::filter(term != "Stage_simpleStageX") %>%
  # Convert to plotting keys that match term_labels_fig6
  dplyr::mutate(
    term_key = dplyr::case_when(
      term == "Age"                  ~ "age",
      term == "score"                ~ "risk_score",
      term == "ras_cna_score"        ~ "ras_cna_score",
      term == "Stage_simpleStageII"  ~ "stage_simpleStageII",
      term == "Stage_simpleStageIII_IV" ~ "stage_simpleStageIII_IV",
      term == "PAM50LumB"            ~ "pam50LumB",
      term == "PAM50Basal"           ~ "pam50Basal",
      term == "PAM50Her2"            ~ "pam50Her2",
      term == "PAM50Normal"          ~ "pam50Normal",
      TRUE                           ~ NA_character_
    ),
    term_clean = dplyr::if_else(
      is.na(term_key),
      term_clean,                            # fallback (won't affect shown terms if all mapped)
      unname(term_labels_fig6[term_key])     # apply your exact labels
    )
  ) %>%
  # Keep only terms you defined in term_labels_fig6 (plot-only)
  dplyr::filter(term_key %in% names(term_labels_fig6)) %>%
  # Order y-axis exactly as in term_labels_fig6 (and reverse for forest style)
  dplyr::mutate(
    term_factor = forcats::fct_rev(
      factor(term_key, levels = names(term_labels_fig6), labels = unname(term_labels_fig6))
    )
  )

## --- Shared x-axis label formatter (fixed decimals) ---
lab_hr_3dec <- function(x) formatC(x, digits = 2, format = "f")

## a plot
# Forest plot (a): risk score HR in METABRIC and TCGA (score-only models)
p_forest_a <- ggplot(df_forest_a,
                     aes(x = HR, y = term_factor)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = HR_low, xmax = HR_high),
                  size = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", HR)),  # numeric HR label above point
    vjust = -1.0, hjust = 0.5,
    fontface = "bold",
    size = 3.5
  ) +
  scale_x_log10(
    breaks = c(1, 3, 10, 30),
    labels = lab_hr_3dec          # fixed decimals on axis ticks
  ) +
  facet_wrap(~Cohort, nrow = 1) +
  theme_bw() +
  labs(
    x = "HR",
    y = "",
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

## b plot
# Forest plot (b): TCGA full model HRs (risk score + clinical covariates)
p_forest_b <- ggplot(
  data = df_forest_b %>%
    dplyr::filter(!stringr::str_detect(as.character(term_factor), "Stage: StageX")),
  aes(x = HR, y = term_factor)
) +
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
    title = "Risk score with clinical covariates (TCGA-BRCA)"
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

## ===================================================================
## 6) KM curves and Number at risk tables (survminer + patchwork)
## ===================================================================

# Announce Kaplan–Meier and risk table generation
cat(">>> Building KM curves and risk tables (survminer-style) ...\n")

## ---------- METABRIC ----------

# Fit KM curves stratified by risk group (METABRIC)
fit_km_met6 <- survfit(Surv(OS_TIME, OS_EVENT) ~ risk_group, data = met_sig)

# Set METABRIC x-axis maximum time (months)
tmax_met <- 400

# Compute log-rank test p-value (METABRIC)
lr_met   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ risk_group, data = met_sig)
p_met_lr <- 1 - pchisq(lr_met$chisq, length(lr_met$n) - 1)

# Build KM plot + number-at-risk table (survminer output object)
g_met <- survminer::ggsurvplot(
  fit_km_met6,
  data               = met_sig,
  risk.table         = TRUE,
  pval               = FALSE,   # add p-value manually
  legend.labs        = c("Low", "High"),
  legend.title       = "Risk group",
  xlab               = "Time (months)",
  ylab               = "OS probability",
  xlim               = c(0, tmax_met),
  ggtheme            = theme_bw(),
  palette            = c("#E41A1C", "#377EB8"),
  risk.table.height  = 0.28
)

# Create annotation text for KM plot (median cutoff + log-rank p)
ann_met_label <- sprintf("median score = %.3g\nlog-rank p = %.3g",
                         cut_median, p_met_lr)

# Style KM plot and add annotation (METABRIC)
p_km_met_main <- g_met$plot +
  annotate("text",
           x = tmax_met * 0.55,
           y = 0.95,
           label = ann_met_label,
           hjust = 0, vjust = 1,
           size  = 4.5,
           fontface = "bold") +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    legend.text  = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold", size = 11)
  ) +
  labs(title = "METABRIC")

# Style number-at-risk table (METABRIC)
p_risk_met <- g_met$table +
  theme_bw() +
  labs(x = NULL) +          # remove x-axis label from the risk table
  theme(
    axis.title.y = element_blank(),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    panel.grid   = element_blank(),
    plot.margin  = margin(t = 0, r = 5, b = 5, l = 5)
  )

# Combine KM plot and risk table vertically (METABRIC)
p_km_met_panel <- p_km_met_main / p_risk_met +
  patchwork::plot_layout(heights = c(3, 1)) &
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))


## ---------- TCGA-BRCA ----------

# Fit KM curves stratified by risk group (TCGA-BRCA)
fit_km_tcga6 <- survfit(Surv(OS_TIME, OS_EVENT) ~ risk_group, data = tcga_sig)

# Set TCGA x-axis maximum time (months)
tmax_tcga <- 300

# Compute log-rank test p-value (TCGA-BRCA)
lr_tcga   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ risk_group, data = tcga_sig)
p_tcga_lr <- 1 - pchisq(lr_tcga$chisq, length(lr_tcga$n) - 1)

# Build KM plot + number-at-risk table (TCGA-BRCA)
g_tcga <- survminer::ggsurvplot(
  fit_km_tcga6,
  data               = tcga_sig,
  risk.table         = TRUE,
  pval               = FALSE,
  legend.labs        = c("Low", "High"),
  legend.title       = "Risk group",
  xlab               = "Time (months)",
  ylab               = "OS probability",
  xlim               = c(0, tmax_tcga),
  ggtheme            = theme_bw(),
  palette            = c("#E41A1C", "#377EB8"),
  risk.table.height  = 0.28
)

# Create annotation text for KM plot (median cutoff + log-rank p)
ann_tcga_label <- sprintf("median score = %.3g\nlog-rank p = %.3g",
                          cut_median, p_tcga_lr)

# Style KM plot and add annotation (TCGA-BRCA)
p_km_tcga_main <- g_tcga$plot +
  annotate("text",
           x = tmax_tcga * 0.55,
           y = 0.95,
           label = ann_tcga_label,
           hjust = 0, vjust = 1,
           size  = 4.5,
           fontface = "bold") +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12,
                                margin = margin(r = 8)),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    legend.text  = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold", size = 11)
  ) +
  labs(title = "TCGA-BRCA")

# Style number-at-risk table (TCGA-BRCA)
p_risk_tcga <- g_tcga$table +
  theme_bw() +
  labs(x = NULL) +          # remove x-axis label from the risk table
  theme(
    axis.title.y = element_blank(),
    axis.text.x  = element_text(face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    panel.grid   = element_blank(),
    plot.margin  = margin(t = 0, r = 5, b = 5, l = 5)
  )

# Combine KM plot and risk table vertically (TCGA-BRCA)
p_km_tcga_panel <- p_km_tcga_main / p_risk_tcga +
  patchwork::plot_layout(heights = c(3, 1)) &
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))


## ===================================================================
## FINAL FIGURE 6  (forest top: a., b.; bottom: c. as KM panels)
## ===================================================================

## ÜST SATIR: a ve b (small spacer; tags only on this row)
# Build top row layout with spacer in the middle and panel labels a/b
top_row <- cowplot::plot_grid(
  p_forest_a,
  NULL,
  p_forest_b,
  ncol           = 3,
  rel_widths     = c(1, 0.08, 2),
  labels         = c("A.", "", "B."),
  label_size     = 16,
  label_fontface = "bold",
  hjust          = -0.1
)

# Build bottom row layout for KM+risk tables with a single panel label c
bottom_row <- cowplot::plot_grid(
  p_km_met_panel,
  NULL,
  p_km_tcga_panel,
  ncol           = 3,
  rel_widths     = c(1, 0.03, 1),
  labels         = c("C.", "", ""),
  label_size     = 16,
  label_fontface = "bold",
  hjust          = -0.1
) +
  theme(plot.margin = margin(t = 20, r = 0, b = 0, l = 0))   # add a small top margin

# Stack top and bottom rows into the final Figure 6 layout
fig6 <- cowplot::plot_grid(
  top_row,
  bottom_row,
  ncol        = 1,
  rel_heights = c(1, 1.3),
  align       = "v"
)

# Save Figure 6 as a high-resolution JPEG
ggsave(
  file.path(out_dir, "Figure 6.jpeg"),
  fig6 + theme(plot.margin = margin(t = 10, r = 40, b = 10, l = 40)),
  width  = 14,
  height = 11,
  dpi    = 300
)

# Print completion message
cat(">>> Figure 6 saved.\n")