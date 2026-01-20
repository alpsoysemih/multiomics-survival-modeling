# Load required libraries quietly
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(survival); library(broom)
  library(tibble); library(ggplot2); library(patchwork)
  library(org.Hs.eg.db); library(cowplot)
})

# Null-coalescing helper: return a unless NULL, else b
`%||%` <- function(a, b) if (!is.null(a)) a else b

## ===================================================================
## 0) PATHS
## ===================================================================

# Define input/output directories
metabric_dir <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
gdsc_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
out_dir      <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"

# Create output directory if needed
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Define signature genes (expression + CNA) and combine
sig_expr <- "SERPINE1"
sig_cna  <- c("CLDN1","CXCR4","LYPD6B","LAMB2")   
sig_all  <- c(sig_expr, sig_cna)

# Print output directory
cat(">>> Output directory:", out_dir, "\n")

# Define a bold theme helper for consistent styling
theme_bold_all <- theme(
  axis.title.x  = element_text(face = "bold"),
  axis.title.y  = element_text(face = "bold"),
  axis.text.x   = element_text(face = "bold"),
  axis.text.y   = element_text(face = "bold"),
  legend.title  = element_text(face = "bold"),
  legend.text   = element_text(face = "bold"),
  strip.text    = element_text(face = "bold")
)

## ===================================================================
## 1) GDSC CNA direction via Fisher
## ===================================================================

# Announce GDSC loading
cat(">>> Loading GDSC ...\n")

# Define GDSC input files
gdsc_cna_file  <- file.path(gdsc_dir,"GDSC_CNA.Paclitaxel.tsv")
gdsc_resp_file <- file.path(gdsc_dir,"GDSC_response.Paclitaxel.tsv")

# Read GDSC CNA and response tables
gdsc_cna  <- read_tsv(gdsc_cna_file, show_col_types = FALSE)
gdsc_resp <- read_tsv(gdsc_resp_file, show_col_types = FALSE)

# Map Entrez IDs to gene symbols for CNA table
gdsc_cna$Gene_Symbol <- mapIds(
  org.Hs.eg.db,
  as.character(gdsc_cna$gene_id),
  "SYMBOL","ENTREZID"
)

# Reorder columns and drop Entrez gene_id
gdsc_cna <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, everything(), -gene_id)

# Binarize CNA values (0 vs non-zero)
gdsc_cna[,-1] <- ifelse(gdsc_cna[,-1] != 0, 1, 0)

# Build a compact response table (sample, status)
resp <- gdsc_resp %>%
  transmute(
    sample = as.character(sample_name),
    status = response
  )

# Convert CNA matrix to long format and join with response
cna_long <- gdsc_cna %>%
  pivot_longer(-Gene_Symbol, names_to = "sample", values_to = "CNA") %>%
  inner_join(resp, by = "sample")

# Compute Fisher exact test OR and p-value per gene (CNA vs response)
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

# Announce METABRIC loading
cat(">>> Loading METABRIC ...\n")

# Define METABRIC input files
expr_file <- file.path(metabric_dir,"data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
cna_file  <- file.path(metabric_dir,"data_cna.txt")
clin_file <- file.path(metabric_dir,"brca_metabric_clinical_data.tsv")

# Read METABRIC expression and standardize gene column name
met_expr <- read_tsv(expr_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = `Hugo_Symbol`) %>%
  dplyr::select(-Entrez_Gene_Id)

# Read METABRIC CNA and standardize gene column name
met_cna <- read_tsv(cna_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = `Hugo_Symbol`) %>%
  dplyr::select(-Entrez_Gene_Id)

# Binarize METABRIC CNA values (0 vs non-zero)
met_cna_bin <- met_cna
met_cna_bin[,-1] <- ifelse(met_cna_bin[,-1] != 0, 1, 0)

# Read and harmonize clinical data (OS, Age, Stage_simple, PAM50)
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
    },
    PAM50 = factor(PAM50)
  ) %>%
  dplyr::select(sample, OS_TIME, OS_EVENT, Age, Stage_simple, PAM50) %>%
  droplevels()

# Convert METABRIC expression to long format for signature genes
expr_long <- met_expr %>%
  filter(Gene_Symbol %in% sig_all) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

# Convert METABRIC CNA to long format and align direction using GDSC Fisher OR
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
    cna_aligned = case_when(
      dir ==  1 ~ cna,
      dir == -1 ~ 1 - cna,
      TRUE      ~ cna   # OR/dir NA ise dokunma
    )
  )

# Extract signature columns per sample
serp  <- expr_long    %>% filter(Gene_Symbol == sig_expr) %>%
  transmute(sample, SERPINE1_expr = expr)
cl1   <- cna_long_met %>% filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = cna_aligned)
lamb2 <- cna_long_met %>% filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = cna_aligned)
cx4   <- cna_long_met %>% filter(Gene_Symbol == "CXCR4") %>%
  transmute(sample, CXCR4_CNA = cna_aligned)
ly6   <- cna_long_met %>% filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = cna_aligned)

# Merge clinical + signature variables and drop missing values
met_sig <- clin %>%
  left_join(serp,  by = "sample") %>%
  left_join(cl1,   by = "sample") %>%
  left_join(cx4,   by = "sample") %>%
  left_join(ly6,   by = "sample") %>%
  left_join(lamb2, by = "sample") %>%   # <-- eklendi
  drop_na()

# Report METABRIC signature sample size
cat(">>> METABRIC signature N:", nrow(met_sig), "\n")

## ===================================================================
## 3) TCGA
## ===================================================================

# Announce TCGA loading
cat(">>> Loading TCGA-BRCA ...\n")

# Read TCGA expression, CNA, and clinical tables
tcga_expr <- read_tsv(file.path(tcga_dir,"TCGA-BRCA_exprs.z.tsv"), show_col_types = FALSE)
tcga_cna  <- read_tsv(file.path(tcga_dir,"TCGA-BRCA_CNA.tsv"),      show_col_types = FALSE)
tcga_clin <- read_tsv(file.path(tcga_dir,"TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv"),
                      show_col_types = FALSE)

# Map Entrez IDs to gene symbols for TCGA expression
tcga_expr$Gene_Symbol <- mapIds(
  org.Hs.eg.db,
  as.character(tcga_expr$ENTREZID),
  "SYMBOL","ENTREZID"
)

# Reorder columns and drop ENTREZID
tcga_expr <- tcga_expr %>%
  dplyr::select(Gene_Symbol, everything(), -ENTREZID)

# Map Entrez IDs to gene symbols for TCGA CNA
tcga_cna$Gene_Symbol <- mapIds(
  org.Hs.eg.db,
  as.character(tcga_cna$gene),
  "SYMBOL","ENTREZID"
)

# Reorder columns and drop Entrez gene column
tcga_cna <- tcga_cna %>%
  dplyr::select(Gene_Symbol, everything(), -gene)

# Binarize TCGA CNA values (0 vs non-zero)
tcga_cna_bin <- tcga_cna
tcga_cna_bin[,-1] <- ifelse(tcga_cna_bin[,-1] != 0, 1, 0)

# Clean TCGA clinical data (drop Normal-like, harmonize OS/Stage/PAM50)
tcga_clin_small <- tcga_clin %>%
  filter(!BRCA_Subtype_PAM50 %in% c("Normal","Normal-like")) %>%
  mutate(
    OS_TIME  = as.numeric(OS.time) / 30.44,
    OS_EVENT = as.numeric(OS),
    Age      = as.numeric(Age),
    Stage_simple = case_when(
      Stage %in% c("Stage I", "Stage IA", "Stage IB")                                 ~ "StageI",
      Stage %in% c("Stage II", "Stage IIA", "Stage IIB")                              ~ "StageII",
      Stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV") ~ "StageIII_IV",
      TRUE                                                                            ~ "StageX"
    ),
    PAM50 = case_when(
      str_detect(BRCA_Subtype_PAM50,"LumA")  ~ "LumA",
      str_detect(BRCA_Subtype_PAM50,"LumB")  ~ "LumB",
      str_detect(BRCA_Subtype_PAM50,"Basal") ~ "Basal",
      str_detect(BRCA_Subtype_PAM50,"Her2")  ~ "Her2",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(OS_TIME)) %>%
  mutate(
    Stage_simple = {
      f <- factor(Stage_simple, levels = c("StageI", "StageII", "StageIII_IV", "StageX"))
      if ("StageI" %in% levels(f)) stats::relevel(f, ref = "StageI") else f
    },
    PAM50 = factor(PAM50)
  ) %>%
  dplyr::select(sample, OS_TIME, OS_EVENT, Age, Stage_simple, PAM50) %>%
  droplevels()

# Convert TCGA expression to long format and standardize sample barcode
expr_tcga <- tcga_expr %>%
  filter(Gene_Symbol %in% sig_all) %>%
  pivot_longer(
    -Gene_Symbol,
    names_to  = "bc",
    values_to = "expr"
  ) %>%
  mutate(sample = substr(bc, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, expr)

# Convert TCGA CNA to long format, standardize sample barcode, and align direction using GDSC Fisher OR
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

# Extract TCGA signature columns per sample
serp_tcga <- expr_tcga %>%
  filter(Gene_Symbol == sig_expr) %>%
  transmute(sample, SERPINE1_expr = expr)

cl1_tcga  <- cna_tcga %>%
  filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = cna_aligned)

lamb2_tcga <- cna_tcga %>%
  filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = cna_aligned)

cx4_tcga  <- cna_tcga %>%
  filter(Gene_Symbol == "CXCR4") %>%
  transmute(sample, CXCR4_CNA = cna_aligned)

ly6_tcga  <- cna_tcga %>%
  filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = cna_aligned)

# Merge TCGA clinical + signature variables and drop missing values
tcga_sig <- tcga_clin_small %>%
  left_join(serp_tcga,  by = "sample") %>%
  left_join(cl1_tcga,   by = "sample") %>%
  left_join(cx4_tcga,   by = "sample") %>%
  left_join(ly6_tcga,   by = "sample") %>%
  left_join(lamb2_tcga, by = "sample") %>%  # <-- eklendi
  drop_na()

# Report TCGA signature sample size
cat(">>> TCGA-BRCA signature N:", nrow(tcga_sig), "\n")

## ===================================================================
## 4) PANEL A — Correlation heatmap (Expr–CNA: point-biserial + Wilcoxon,
##                                  CNA–CNA: phi + Fisher)
## ===================================================================

corr_vars <- c("SERPINE1_expr","CLDN1_CNA","CXCR4_CNA","LYPD6B_CNA","LAMB2_CNA")

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# phi coefficient (binary-binary) from a 2x2 table
phi_from_2x2 <- function(tab) {
  if (!all(dim(tab) == c(2, 2))) return(NA_real_)
  a <- as.numeric(tab[1,1]); b <- as.numeric(tab[1,2])
  c <- as.numeric(tab[2,1]); d <- as.numeric(tab[2,2])
  denom <- (a + b) * (c + d) * (a + c) * (b + d)
  if (is.na(denom) || denom == 0) return(NA_real_)
  (a*d - b*c) / sqrt(denom)
}

build_corr_df <- function(df, cohort_label) {
  
  vars <- corr_vars
  out  <- vector("list", length(vars) * length(vars))
  k <- 1
  
  for (v1 in vars) {
    for (v2 in vars) {
      
      r <- NA_real_
      p <- NA_real_
      
      if (v1 == v2) {
        
        r <- 1
        p <- NA_real_
        
      } else {
        
        is_expr1 <- grepl("_expr$", v1)
        is_expr2 <- grepl("_expr$", v2)
        is_cna1  <- grepl("_CNA$",  v1)
        is_cna2  <- grepl("_CNA$",  v2)
        
        # Expr–CNA (either direction): point-biserial correlation + Wilcoxon p
        if ((is_expr1 && is_cna2) || (is_cna1 && is_expr2)) {
          
          # ensure x is expression, y is CNA
          if (is_expr1) { x <- df[[v1]]; y <- df[[v2]] } else { x <- df[[v2]]; y <- df[[v1]] }
          
          # force CNA to clean 0/1
          y01 <- dplyr::if_else(is.na(y), NA_real_, as.numeric(y != 0))
          
          # need exactly 2 levels to test
          y_levels <- sort(unique(y01[!is.na(y01)]))
          if (length(y_levels) < 2) {
            r <- NA_real_
            p <- NA_real_
          } else {
            r <- suppressWarnings(cor(x, y01, use = "pairwise.complete.obs"))
            p <- tryCatch(
              wilcox.test(x ~ factor(y01))$p.value,
              error = function(e) NA_real_
            )
          }
          
          # CNA–CNA: phi + Fisher
        } else if (is_cna1 && is_cna2) {
          
          y1 <- dplyr::if_else(is.na(df[[v1]]), NA_real_, as.numeric(df[[v1]] != 0))
          y2 <- dplyr::if_else(is.na(df[[v2]]), NA_real_, as.numeric(df[[v2]] != 0))
          
          tab <- table(y1, y2, useNA = "no")
          
          if (all(dim(tab) == c(2, 2))) {
            r <- phi_from_2x2(tab)
            p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
          } else {
            r <- NA_real_
            p <- NA_real_
          }
        }
      }
      
      out[[k]] <- tibble::tibble(
        var1   = v1,
        var2   = v2,
        r      = r,
        p      = p,
        star   = p_to_stars(p),
        Cohort = cohort_label
      )
      k <- k + 1
    }
  }
  
  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      var1 = factor(var1, levels = vars),
      var2 = factor(var2, levels = rev(vars))
    )
}

corr_df <- dplyr::bind_rows(
  build_corr_df(met_sig,  "METABRIC"),
  build_corr_df(tcga_sig, "TCGA-BRCA")
)

corr_df <- corr_df %>%
  dplyr::mutate(
    star = dplyr::case_when(
      is.na(p)     ~ "",
      p < 0.001    ~ "***",
      p < 0.01     ~ "**",
      p < 0.05     ~ "*",
      TRUE         ~ ""
    ),
    lab = dplyr::if_else(
      is.na(r),
      "",
      paste0(sprintf("%.2f", r), star)
    )
  )

p_cor <- ggplot(corr_df, aes(x = var1, y = var2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = lab), size = 4, fontface = "bold") +
  facet_wrap(~Cohort, nrow = 1) +
  scale_fill_gradient2(
    low  = "#4575B4",
    mid  = "white",
    high = "#D73027",
    limits = c(-1, 1),
    name   = "Correlation",
    labels = function(x) sprintf("%.2f", x)
  ) +
  labs(x = NULL, y = NULL, title = "Correlation") +
  theme_bw() +
  theme(
    plot.title   = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.text.x  = element_text(size = 12, angle = 40, hjust = 1, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 12, face = "bold"),
    strip.text   = element_text(size = 12, face = "bold")
  )

## ===================================================================
## 5) Cox models (Signature-only vs +clinical)
## ===================================================================

# Define Cox model formulas (use model strings as labels)
model_spec <- list(
  "SERPINE1+CLDN1+LAMB2+CXCR4"  =
    as.formula("Surv(OS_TIME,OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LAMB2_CNA + CXCR4_CNA"),
  "SERPINE1+CLDN1+LAMB2+LYPD6B" =
    as.formula("Surv(OS_TIME,OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LAMB2_CNA + LYPD6B_CNA")
)

# Fit signature-only Cox models on METABRIC
cox_sig_met  <- lapply(model_spec, function(f) coxph(f, data = met_sig))

# Fit signature+clinical Cox models on METABRIC
cox_full_met <- lapply(model_spec, function(f) {
  coxph(update(f, . ~ . + Age + Stage_simple + PAM50), data = met_sig)
})

# Fit signature-only Cox models on TCGA
cox_sig_tcga  <- lapply(model_spec, function(f) coxph(f, data = tcga_sig))

# Fit signature+clinical Cox models on TCGA
cox_full_tcga <- lapply(model_spec, function(f) {
  coxph(update(f, . ~ . + Age + Stage_simple + PAM50), data = tcga_sig)
})

# Helper: compute C-index and its SE for a fitted Cox model
cindex_fun <- function(fit, data) {
  lp <- predict(fit, newdata = data)
  s  <- summary(coxph(Surv(OS_TIME, OS_EVENT) ~ lp, data = data))
  c(s$concordance[1], s$concordance[2])
}

## ===================================================================
## PANEL B — C-index (Signature-only vs Signature+clinical, METABRIC & TCGA)
## ===================================================================

# Collect C-index results across models/cohorts/groups
ci_rows <- list()

# Loop over models and compute C-index for both cohorts and both model groups
for (m in names(model_spec)) {
  cm <- cindex_fun(cox_sig_met[[m]], met_sig)
  ci_rows[[length(ci_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "METABRIC",
    Group  = "Signature-only",
    C      = cm[1],
    SE     = cm[2]
  )
  ct <- cindex_fun(cox_sig_tcga[[m]], tcga_sig)
  ci_rows[[length(ci_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "TCGA-BRCA",
    Group  = "Signature-only",
    C      = ct[1],
    SE     = ct[2]
  )
  cmf <- cindex_fun(cox_full_met[[m]], met_sig)
  ci_rows[[length(ci_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "METABRIC",
    Group  = "Signature + Clinical",
    C      = cmf[1],
    SE     = cmf[2]
  )
  ctf <- cindex_fun(cox_full_tcga[[m]], tcga_sig)
  ci_rows[[length(ci_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "TCGA-BRCA",
    Group  = "Signature + Clinical",
    C      = ctf[1],
    SE     = ctf[2]
  )
}

# Combine C-index rows into a single table
df_ci <- bind_rows(ci_rows)

# Order facets left-to-right
df_ci$Group <- factor(df_ci$Group, levels = c("Signature-only","Signature + Clinical"))

# Plot C-index bars with error bars and labels
p_ci <- ggplot(df_ci, aes(x = Model, y = C, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  geom_errorbar(
    aes(ymin = C - SE, ymax = C + SE),
    position = position_dodge(width = 0.6),
    width = 0.2
  ) +
  geom_text(
    aes(label = sprintf("%.3f", C), y = C + SE + 0.025),
    position = position_dodge(width = 0.6),
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(~Group, nrow = 1) +
  labs(
    x     = "",
    y     = "C-index",
    title = "Model accuracy (C-index)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", hjust = .5, size = 15),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(face = "bold", size = 12, angle = 30, hjust = 1),
    axis.text.y  = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 12),
    legend.text  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12)
  )

## ===================================================================
## PANEL C — AIC (METABRIC & TCGA, Signature-only vs Signature+clinical)
## ===================================================================

# (Optional) placeholder list if you want loop-built AIC rows
aic_rows <- list()

# Compute AIC values per model/cohort/group (loop-based)
for (m in names(model_spec)) {
  aic_rows[[length(aic_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "METABRIC",
    Group  = "Signature-only",
    AIC    = AIC(cox_sig_met[[m]])
  )
  aic_rows[[length(aic_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "METABRIC",
    Group  = "Signature + Clinical",
    AIC    = AIC(cox_full_met[[m]])
  )
  aic_rows[[length(aic_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "TCGA-BRCA",
    Group  = "Signature-only",
    AIC    = AIC(cox_sig_tcga[[m]])
  )
  aic_rows[[length(aic_rows) + 1]] <- tibble(
    Model  = m,
    Cohort = "TCGA-BRCA",
    Group  = "Signature + Clinical",
    AIC    = AIC(cox_full_tcga[[m]])
  )
}

# Build AIC table (explicitly enumerated) and compute per-facet label offsets
df_aic <- bind_rows(aic_rows) %>%
  group_by(Group) %>%
  mutate(
    aic_offset = 0.02 * (max(AIC) - min(AIC)),
    label_y    = AIC + aic_offset
  ) %>%
  ungroup()

df_aic$Group <- factor(df_aic$Group,
                       levels = c("Signature-only","Signature + Clinical"))

# Order facets left-to-right
df_aic$Group <- factor(df_aic$Group, levels = c("Signature-only","Signature + Clinical"))

# Plot AIC bars with labels above bars
p_aic <- ggplot(df_aic, aes(x = Model, y = AIC, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  geom_text(
    aes(y = label_y, label = round(AIC, 0)),
    position = position_dodge(width = 0.6),
    vjust = 0,
    fontface = "bold",
    size = 4
  ) +
  facet_wrap(~Group, nrow = 1) +
  theme_bw() +
  labs(
    x     = "",
    y     = "AIC",
    title = "Model complexity (AIC)"
  ) +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 15),
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 5)),
    strip.text   = element_text(face = "bold", size = 12),
    legend.text  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12)
  )

## ===================================================================
## PANEL D — LRT vs SERPINE1+CLDN1 (Signature-only & Signature+Clinical)
## ===================================================================

# Helper: LRT p-value for nested Cox models
lrt_compare <- function(form_small, form_big, data) {
  broom::tidy(
    anova(
      coxph(form_small, data = data),
      coxph(form_big,  data = data),
      test = "LRT"
    )
  )$p.value[2]
}

# Define baseline signature-only model
form_sig0       <- Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LAMB2_CNA

# Define signature-only extensions (add one CNA feature)
form_sig_cxcr4  <- update(form_sig0, . ~ . + CXCR4_CNA)
form_sig_lypd6b <- update(form_sig0, . ~ . + LYPD6B_CNA)

# Define clinical-only model (MAKE SURE THIS EXISTS)
form_clin <- Surv(OS_TIME, OS_EVENT) ~ Age + Stage_simple + PAM50
form_clin_sig0      <- update(form_clin, . ~ . + SERPINE1_expr + CLDN1_CNA + LAMB2_CNA)

# Define clinical + extended signature models
form_clin_sig_cxcr4 <- update(form_clin_sig0, . ~ . + CXCR4_CNA)
form_clin_sig_lypd6b<- update(form_clin_sig0, . ~ . + LYPD6B_CNA)

# Convert p-values to significance stars
p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Collect LRT results across cohorts and model groups
rows <- list()

# Loop over cohorts to compare nested models with LRT
for (cohort_name in c("METABRIC","TCGA-BRCA")) {
  
  # Select cohort-specific signature table
  dat <- if (cohort_name == "METABRIC") met_sig else tcga_sig
  
  # LRT: add CXCR4 vs baseline (signature-only)
  p_cx  <- lrt_compare(form_sig0,      form_sig_cxcr4,  dat)
  
  # LRT: add LYPD6B vs baseline (signature-only)
  p_ly  <- lrt_compare(form_sig0,      form_sig_lypd6b, dat)
  
  # Store signature-only LRT p-values
  rows[[length(rows) + 1]] <- tibble(
    Model  = c("SERPINE1+CLDN1+LAMB2+CXCR4",
               "SERPINE1+CLDN1+LAMB2+LYPD6B"),
    Cohort = cohort_name,
    Group  = "Signature-only",
    p      = c(p_cx, p_ly)
  )
  
  # LRT: add CXCR4 vs baseline (signature+clinical)
  p_cx_cl <- lrt_compare(form_clin_sig0, form_clin_sig_cxcr4,  dat)
  
  # LRT: add LYPD6B vs baseline (signature+clinical)
  p_ly_cl <- lrt_compare(form_clin_sig0, form_clin_sig_lypd6b, dat)
  
  # Store signature+clinical LRT p-values
  rows[[length(rows) + 1]] <- tibble(
    Model  = c("SERPINE1+CLDN1+LAMB2+CXCR4",
               "SERPINE1+CLDN1+LAMB2+LYPD6B"),
    Cohort = cohort_name,
    Group  = "Signature + Clinical",
    p      = c(p_cx_cl, p_ly_cl)
  )
}

# Format LRT results for plotting (-log10(p), stars, label positions)
df_lrt <- bind_rows(rows) %>%
  mutate(
    Group = factor(Group,
                   levels = c("Signature-only","Signature + Clinical")),
    Model = factor(Model,
                   levels = c("SERPINE1+CLDN1+LAMB2+CXCR4",
                              "SERPINE1+CLDN1+LAMB2+LYPD6B")),
    neg   = -log10(p),
    star  = vapply(p, p_to_stars, character(1)),
    label_y  = neg + 0.03,
    label_txt = ifelse(
      star == "",
      sprintf("%.2f", neg),
      paste0(sprintf("%.2f", neg), star)
    )
  )

# Plot LRT bars with significance threshold line and labels
p_lrt <- ggplot(df_lrt, aes(x = Model, y = neg, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.6) +
  geom_text(
    aes(y = label_y, label = label_txt),
    position = position_dodge(width = 0.6),
    vjust = 0,
    fontface = "bold",
    size = 4
  ) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "red") +
  facet_wrap(~Group, nrow = 1) +
  theme_bw() +
  labs(
    x     = "",
    y     = expression(bold(-log10(p))),
    title = "Feature contribution (LRT)"
  ) +
  theme(
    plot.title   = element_text(face = "bold", hjust = .5, size = 15),
    axis.text.x  = element_text(angle = 30, hjust = 1, face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    strip.text   = element_text(face = "bold", size = 12),
    legend.text  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12)
  )

## ===================================================================
## FINAL FIGURE — assemble panels and save
## ===================================================================

# Add left margin to all panels (prevents left-side clipping)
left_margin_theme <- theme(
  plot.margin = margin(t = 5, r = 5, b = 5, l = 65)  
)

# Apply left margin to all panels
p_aic <- p_aic + left_margin_theme

# Build the top row: Panel (a) correlation heatmap + Panel (b) C-index barplot
row_top <- cowplot::plot_grid(
  p_cor, p_ci,
  labels         = c("A.", "B."),     # panel tags
  label_size     = 18,               # tag font size
  label_fontface = "bold",           # tag font style
  ncol           = 2,                # two plots side-by-side
  align          = "h",              # align horizontally within the row
  axis           = "tb"              # align top/bottom axes
)

# Build the bottom row: Panel (c) AIC barplot + Panel (d) LRT barplot
row_bottom <- cowplot::plot_grid(
  p_aic, p_lrt,
  labels         = c("C.", "D."),     # panel tags
  label_size     = 18,               # tag font size
  label_fontface = "bold",           # tag font style
  ncol           = 2,                # two plots side-by-side
  align          = "h",              # align horizontally within the row
  axis           = "tb"              # align top/bottom axes
)

# Stack the two rows to form the final 2x2 multi-panel figure
fig <- cowplot::plot_grid(
  row_top,
  row_bottom,
  ncol        = 1,                   # rows stacked vertically
  rel_heights = c(1, 1)              # equal row heights
)

# Save the final figure as a high-resolution JPEG
ggsave(
  file.path(out_dir, "Figure 5.jpeg"),
  fig,
  width  = 20,                       # output width (inches)
  height = 11,                       # output height (inches)
  dpi    = 300                       # resolution
)

# Print completion message
cat(">>> DONE.\n")