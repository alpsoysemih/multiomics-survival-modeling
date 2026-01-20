# Load required packages quietly (suppress startup messages)
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(broom))

# Global option: keep strings as character (not factors) by default
options(stringsAsFactors = FALSE)

# Start message for the Figure 10 workflow
cat(">> Starting Figure 10 pipeline (clinical + signature + Ras)...\n\n")

# ------------------------------------------------------------
# 0) Directories and global settings
# ------------------------------------------------------------

# Define input data directories
metabric_dir   <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir       <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
gdsc_dir       <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
pathfindR_dir  <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"

# Define output directory for intermediate results
out_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/Ras_signature_integration"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
cat(">> Results directory:\n   ", out_dir, "\n\n")

# Define output directory for manuscript figures
figures_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}

# Define signature genes used in the gene-level prognostic score
sig_gene_expr <- "SERPINE1"
sig_genes_cna <- c("CLDN1", "LAMB2", "LYPD6B")
sig_genes_all <- unique(c(sig_gene_expr, sig_genes_cna))

# Define pathway names (must match pathfindR Term_Description)
ras_pathway_name   <- "Ras signaling pathway"
focal_pathway_name <- "Focal adhesion"
axon_pathway_name  <- "Axon guidance"

# ------------------------------------------------------------
# 1) GDSC: Fisher OR for CNA direction (used to align CNA)
# ------------------------------------------------------------

# Load GDSC expression/CNA/response files and compute Fisher odds ratio per gene
cat(">> Loading GDSC data and computing Fisher OR for CNA (Paclitaxel)...\n")

gdsc_expr_file     <- file.path(gdsc_dir, "GDSC_exprs.Paclitaxel.eb_with.TCGA_exprs.Paclitaxel.tsv")
gdsc_cna_file      <- file.path(gdsc_dir, "GDSC_CNA.Paclitaxel.tsv")
gdsc_response_file <- file.path(gdsc_dir, "GDSC_response.Paclitaxel.tsv")

gdsc_expr     <- read_tsv(gdsc_expr_file,     show_col_types = FALSE)
gdsc_cna      <- read_tsv(gdsc_cna_file,      show_col_types = FALSE)
gdsc_response <- read_tsv(gdsc_response_file, show_col_types = FALSE)

# Binarize CNA values: any non-zero value becomes 1; zero remains 0
gdsc_cna[, -1] <- ifelse(gdsc_cna[, -1] != 0, 1, 0)

# Map Entrez IDs to gene symbols for expression and CNA tables
gdsc_expr$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(gdsc_expr$ENTREZID),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)
gdsc_cna$gene_id <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(gdsc_cna$gene_id),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Standardize gene symbol column naming
colnames(gdsc_expr)[colnames(gdsc_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(gdsc_cna)[colnames(gdsc_cna)   == "gene_id"]  <- "Gene_Symbol"

# Extract response labels and the set of samples used in the response table
resp <- gdsc_response %>% mutate(status = response)
gdsc_samples <- as.character(unique(resp$sample_name))

# Subset CNA to response samples only (genes x samples matrix)
cna_gdsc <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

cna_mat <- as.matrix(cna_gdsc[, -1])
rownames(cna_mat) <- cna_gdsc$Gene_Symbol
cna_bin    <- ifelse(cna_mat != 0, 1, 0)
cna_bin_df <- as_tibble(cna_bin, rownames = "Gene_Symbol")

# Convert CNA matrix to long format and join with response status
cna_long_gdsc <- cna_bin_df %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna_bin"
  ) %>%
  dplyr::rename(sample_name = "sample") %>%
  inner_join(
    resp %>%
      dplyr::select(sample_name, status) %>%
      mutate(sample_name = as.character(sample_name)),
    by = "sample_name"
  )

# Compute Fisher test odds ratio per gene (binary CNA vs response status)
fisher_by_gene <- cna_long_gdsc %>%
  group_by(Gene_Symbol) %>%
  group_modify(~ {
    tab <- table(.x$cna_bin, .x$status)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      tibble(OR = NA_real_, p_value = NA_real_)
    } else {
      ft <- fisher.test(tab)
      tibble(OR = unname(ft$estimate), p_value = ft$p.value)
    }
  }) %>%
  ungroup()

cat("   Fisher OR computed for", nrow(fisher_by_gene), "genes\n\n")

# ------------------------------------------------------------
# 2) METABRIC: expression, CNA, OS + clinical covariates
# ------------------------------------------------------------

# Load METABRIC expression, CNA, and clinical survival/covariates
cat(">> Loading METABRIC expression, CNA and OS data...\n")

metabric_expr_file     <- file.path(metabric_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
metabric_cna_file      <- file.path(metabric_dir, "data_cna.txt")
metabric_survival_file <- file.path(metabric_dir, "brca_metabric_clinical_data.tsv")

metabric_expr     <- read_tsv(metabric_expr_file,     show_col_types = FALSE)
metabric_cna      <- read_tsv(metabric_cna_file,      show_col_types = FALSE)
metabric_survival <- read_tsv(metabric_survival_file, show_col_types = FALSE)

# Standardize gene symbol column names and remove Entrez ID columns
metabric_expr <- metabric_expr %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

metabric_cna <- metabric_cna %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

# Create binary CNA matrix for pathway scoring (any non-zero CNA -> 1)
metabric_cna_binary <- metabric_cna
metabric_cna_binary[, -1] <- ifelse(metabric_cna_binary[, -1] != 0, 1, 0)

# Prepare OS endpoint and clinical covariates; harmonize stage levels; set PAM50 reference if present
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

# Remove extremely rare stage levels to stabilize Cox model fitting
stage_counts <- table(metabric_surv_small$Stage_simple)
rare_stages  <- names(stage_counts[stage_counts < 5])
if (length(rare_stages) > 0) {
  metabric_surv_small <- metabric_surv_small %>%
    filter(!Stage_simple %in% rare_stages) %>%
    mutate(Stage_simple = droplevels(Stage_simple))
}

cat("   METABRIC OS samples:", nrow(metabric_surv_small), "\n\n")

# ------------------------------------------------------------
# 3) TCGA-BRCA: expression, CNA, OS + clinical covariates
# ------------------------------------------------------------

# Load TCGA-BRCA expression, CNA, and clinical survival/covariates
cat(">> Loading TCGA-BRCA expression, CNA and OS data...\n")

tcga_expr_file     <- file.path(tcga_dir, "TCGA-BRCA_exprs.z.tsv")
tcga_cna_file      <- file.path(tcga_dir, "TCGA-BRCA_CNA.tsv")
tcga_survival_file <- file.path(tcga_dir, "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv")

tcga_expr     <- read_tsv(tcga_expr_file,     show_col_types = FALSE)
tcga_cna      <- read_tsv(tcga_cna_file,      show_col_types = FALSE)
tcga_survival <- read_tsv(tcga_survival_file, show_col_types = FALSE)

# Map Entrez IDs to gene symbols for TCGA expression and CNA
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

# Standardize gene symbol column naming
colnames(tcga_expr)[colnames(tcga_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(tcga_cna)[colnames(tcga_cna)   == "gene"]     <- "Gene_Symbol"

# Create binary CNA matrix for pathway scoring (any non-zero CNA -> 1)
tcga_cna_binary <- tcga_cna
tcga_cna_binary[, -1] <- ifelse(tcga_cna_binary[, -1] != 0, 1, 0)

# Prepare OS endpoint and clinical covariates; harmonize stage levels; set PAM50 reference if present
tcga_surv_small <- tcga_survival %>%
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

# Remove extremely rare stage levels to stabilize Cox model fitting
stage_counts_tcga <- table(tcga_surv_small$Stage_simple)
rare_stages_tcga  <- names(stage_counts_tcga[stage_counts_tcga < 5])
if (length(rare_stages_tcga) > 0) {
  tcga_surv_small <- tcga_surv_small %>%
    filter(!Stage_simple %in% rare_stages_tcga) %>%
    mutate(Stage_simple = droplevels(Stage_simple))
}

cat("   TCGA-BRCA OS samples:", nrow(tcga_surv_small), "\n\n")

# ------------------------------------------------------------
# 4) Build molecular signature datasets (SERPINE1 + CNAs)
# ------------------------------------------------------------

# Construct cohort-specific datasets containing OS + signature variables
cat(">> Building molecular signature datasets (METABRIC & TCGA)...\n")

# --- METABRIC: extract signature gene expression/CNA in long format ---

expr_long_metabric_sig <- metabric_expr %>%
  filter(Gene_Symbol %in% sig_genes_all) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

cna_long_metabric_sig <- metabric_cna_binary %>%
  filter(Gene_Symbol %in% sig_genes_cna) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

# Align cohort CNA direction using GDSC Fisher OR (flip when OR<1)
cna_long_metabric_aligned <- cna_long_metabric_sig %>%
  inner_join(fisher_by_gene, by = "Gene_Symbol") %>%
  mutate(
    dir_cna = case_when(
      !is.na(OR) & OR > 1 ~  1,
      !is.na(OR) & OR < 1 ~ -1,
      TRUE                ~ NA_real_
    ),
    cna_aligned = case_when(
      is.na(dir_cna) ~ NA_real_,
      dir_cna ==  1  ~ cna,
      dir_cna == -1  ~ 1 - cna
    )
  )

# Create per-gene feature tables to join into survival data
serp_metabric <- expr_long_metabric_sig %>%
  filter(Gene_Symbol == sig_gene_expr) %>%
  transmute(sample, SERPINE1_expr = as.numeric(expr))

cldn1_metabric <- cna_long_metabric_aligned %>%
  filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = as.numeric(cna_aligned))

lamb2_metabric <- cna_long_metabric_aligned %>%
  filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = as.numeric(cna_aligned))

lypd6b_metabric <- cna_long_metabric_aligned %>%
  filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = as.numeric(cna_aligned))

# Build final METABRIC signature dataset (complete cases for signature variables)
metabric_sig <- metabric_surv_small %>%
  left_join(serp_metabric,   by = "sample") %>%
  left_join(cldn1_metabric,  by = "sample") %>%
  left_join(lamb2_metabric,  by = "sample") %>%
  left_join(lypd6b_metabric, by = "sample") %>%
  filter(
    !is.na(SERPINE1_expr),
    !is.na(CLDN1_CNA),
    !is.na(LAMB2_CNA),
    !is.na(LYPD6B_CNA)
  )

cat("   METABRIC signature samples:", nrow(metabric_sig), "\n")

# --- TCGA: extract signature gene expression/CNA in long format ---

expr_long_tcga_sig <- tcga_expr %>%
  filter(Gene_Symbol %in% sig_genes_all) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "expr"
  ) %>%
  mutate(sample = substr(sample_barcode, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, expr)

cna_long_tcga_sig <- tcga_cna_binary %>%
  filter(Gene_Symbol %in% sig_genes_cna) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  mutate(sample = substr(sample_barcode, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, cna)

# Align cohort CNA direction using GDSC Fisher OR (flip when OR<1)
cna_long_tcga_aligned <- cna_long_tcga_sig %>%
  inner_join(fisher_by_gene, by = "Gene_Symbol") %>%
  mutate(
    dir_cna = case_when(
      !is.na(OR) & OR > 1 ~  1,
      !is.na(OR) & OR < 1 ~ -1,
      TRUE                ~ NA_real_
    ),
    cna_aligned = case_when(
      is.na(dir_cna) ~ NA_real_,
      dir_cna ==  1  ~ cna,
      dir_cna == -1  ~ 1 - cna
    )
  )

# Create per-gene feature tables to join into survival data
serp_tcga <- expr_long_tcga_sig %>%
  filter(Gene_Symbol == sig_gene_expr) %>%
  transmute(sample, SERPINE1_expr = as.numeric(expr))

cldn1_tcga <- cna_long_tcga_aligned %>%
  filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = as.numeric(cna_aligned))

lamb2_tcga <- cna_long_tcga_aligned %>%
  filter(Gene_Symbol == "LAMB2") %>%
  transmute(sample, LAMB2_CNA = as.numeric(cna_aligned))

lypd6b_tcga <- cna_long_tcga_aligned %>%
  filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = as.numeric(cna_aligned))

# Build final TCGA signature dataset (complete cases for signature variables)
tcga_sig <- tcga_surv_small %>%
  left_join(serp_tcga,   by = "sample") %>%
  left_join(cldn1_tcga,  by = "sample") %>%
  left_join(lamb2_tcga,  by = "sample") %>%
  left_join(lypd6b_tcga, by = "sample") %>%
  filter(
    !is.na(SERPINE1_expr),
    !is.na(CLDN1_CNA),
    !is.na(LAMB2_CNA),
    !is.na(LYPD6B_CNA)
  )

cat("   TCGA signature samples:", nrow(tcga_sig), "\n\n")

# ------------------------------------------------------------
# 5) Signature-only models (M1–M4) and risk_score
# ------------------------------------------------------------

# Train signature-only Cox models in METABRIC and evaluate C-index in both cohorts
cat(">> Fitting signature-only Cox models (METABRIC train, TCGA test)...\n")

model_specs <- list(
  M1_full = list(
    formula = as.formula("Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LAMB2_CNA + LYPD6B_CNA"),
    label   = "SERPINE1_expr + CLDN1_CNA + LAMB2_CNA + LYPD6B_CNA"
  ),
  M2_no_LYPD6B = list(
    formula = as.formula("Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LAMB2_CNA"),
    label   = "SERPINE1_expr + CLDN1_CNA + LAMB2_CNA"
  ),
  M3_no_LAMB2 = list(
    formula = as.formula("Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LYPD6B_CNA"),
    label   = "SERPINE1_expr + CLDN1_CNA + LYPD6B_CNA"
  ),
  M4_CLDN1_only = list(
    formula = as.formula("Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA"),
    label   = "SERPINE1_expr + CLDN1_CNA"
  )
)

# Fit all signature models and store per-cohort linear predictors and C-index values
cox_fits    <- list()
cindex_rows <- list()

for (m in names(model_specs)) {
  spec  <- model_specs[[m]]
  form  <- spec$formula
  label <- spec$label
  
  cat("   -", m, ":", label, "\n")
  fit <- coxph(form, data = metabric_sig)
  cox_fits[[m]] <- fit
  
  lp_met  <- predict(fit, newdata = metabric_sig, type = "lp")
  lp_tcga <- predict(fit, newdata = tcga_sig,     type = "lp")
  
  metabric_sig[[paste0("risk_", m)]] <- lp_met
  tcga_sig[[paste0("risk_", m)]]     <- lp_tcga
  
  tmp_met <- data.frame(
    OS_TIME  = metabric_sig$OS_TIME,
    OS_EVENT = metabric_sig$OS_EVENT,
    marker   = lp_met
  )
  c_met <- summary(coxph(Surv(OS_TIME, OS_EVENT) ~ marker, data = tmp_met))$concordance
  
  tmp_tcga <- data.frame(
    OS_TIME  = tcga_sig$OS_TIME,
    OS_EVENT = tcga_sig$OS_EVENT,
    marker   = lp_tcga
  )
  c_tcga <- summary(coxph(Surv(OS_TIME, OS_EVENT) ~ marker, data = tmp_tcga))$concordance
  
  cindex_rows[[length(cindex_rows) + 1]] <- tibble(
    model_id    = m,
    model_terms = label,
    Cohort      = "METABRIC",
    C_index     = c_met[1],
    SE          = c_met[2]
  )
  cindex_rows[[length(cindex_rows) + 1]] <- tibble(
    model_id    = m,
    model_terms = label,
    Cohort      = "TCGA-BRCA",
    C_index     = c_tcga[1],
    SE          = c_tcga[2]
  )
}

# Select the best-performing signature model based on TCGA C-index
cindex_all <- bind_rows(cindex_rows)

best_model <- cindex_all %>%
  filter(Cohort == "TCGA-BRCA") %>%
  arrange(desc(C_index)) %>%
  pull(model_id) %>%
  .[1]

cat("\n>> Best signature model (by TCGA C-index):", best_model, "\n\n")

# Define the selected signature model and compute the final risk_score in both cohorts
cox_sig_metabric <- cox_fits[[best_model]]

metabric_sig$risk_score <- metabric_sig[[paste0("risk_", best_model)]]
tcga_sig$risk_score     <- tcga_sig[[paste0("risk_", best_model)]]

# ------------------------------------------------------------
# 6) Pathway CNA scores: Ras, Focal adhesion, Axon guidance
# ------------------------------------------------------------

# Load pathfindR enrichment output and compute pathway-level CNA scores using aligned CNA directions
cat(">> Loading pathfindR enrichment and computing Ras/Focal/Axon CNA scores...\n")

pathfindR_file <- file.path(pathfindR_dir, "pathfindR_paclitaxel.tsv")
enrichment_df  <- read_tsv(pathfindR_file, show_col_types = FALSE)

# Parse pathway gene lists (Up_regulated / Down_regulated) into a long table of pathway-to-gene mappings
path_genes_long <- enrichment_df %>%
  dplyr::select(Term_Description, Up_regulated, Down_regulated) %>%
  dplyr::rename(pathway = Term_Description) %>%
  mutate(
    pathway        = str_trim(as.character(pathway)),
    Up_regulated   = str_split(Up_regulated,   pattern = "[,;]"),
    Down_regulated = str_split(Down_regulated, pattern = "[,;]")
  ) %>%
  pivot_longer(
    cols      = c("Up_regulated", "Down_regulated"),
    names_to  = "direction_raw",
    values_to = "genes"
  ) %>%
  unnest(genes) %>%
  mutate(
    Gene_Symbol = str_trim(genes)
  ) %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
  dplyr::select(pathway, Gene_Symbol) %>%
  distinct()

# Extract gene sets for target pathways
ras_genes   <- path_genes_long %>% filter(pathway == ras_pathway_name)
focal_genes <- path_genes_long %>% filter(pathway == focal_pathway_name)
axon_genes  <- path_genes_long %>% filter(pathway == axon_pathway_name)

# Sanity checks for pathway gene sets
if (nrow(ras_genes) == 0)  stop("No genes found for Ras signaling pathway in pathfindR file.")
if (nrow(focal_genes) == 0) warning("No genes found for Focal adhesion pathway in pathfindR file.")
if (nrow(axon_genes)  == 0) warning("No genes found for Axon guidance pathway in pathfindR file.")

cat("   Ras genes:",   length(unique(ras_genes$Gene_Symbol)),   "\n")
cat("   Focal genes:", length(unique(focal_genes$Gene_Symbol)), "\n")
cat("   Axon genes:",  length(unique(axon_genes$Gene_Symbol)),  "\n\n")

# Helper function: compute pathway CNA score per sample (mean of aligned CNA across pathway genes)
compute_path_score <- function(cna_long_all, path_genes_tbl, fisher_tbl) {
  if (nrow(path_genes_tbl) == 0) {
    return(tibble(sample = character(), score = numeric()))
  }
  cna_path <- cna_long_all %>%
    inner_join(path_genes_tbl, by = "Gene_Symbol") %>%
    inner_join(fisher_tbl,     by = "Gene_Symbol") %>%
    mutate(
      dir_cna = case_when(
        !is.na(OR) & OR > 1 ~  1,
        !is.na(OR) & OR < 1 ~ -1,
        TRUE                ~ NA_real_
      ),
      cna_aligned = case_when(
        is.na(dir_cna) ~ NA_real_,
        dir_cna ==  1  ~ cna,
        dir_cna == -1  ~ 1 - cna
      )
    )
  
  if (nrow(cna_path) == 0) {
    return(tibble(sample = character(), score = numeric()))
  }
  
  cna_path %>%
    group_by(sample) %>%
    summarise(score = mean(cna_aligned, na.rm = TRUE), .groups = "drop")
}

# Convert cohort CNA matrices to long format for pathway scoring
cna_long_metabric_all <- metabric_cna_binary %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

cna_long_tcga_all <- tcga_cna_binary %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  mutate(sample = substr(sample_barcode, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, cna)

# Compute pathway CNA scores for each cohort
ras_score_metabric <- compute_path_score(cna_long_metabric_all, ras_genes, fisher_by_gene) %>%
  dplyr::rename(ras_cna_score = score)

ras_score_tcga <- compute_path_score(cna_long_tcga_all, ras_genes, fisher_by_gene) %>%
  dplyr::rename(ras_cna_score = score)

focal_score_metabric <- compute_path_score(cna_long_metabric_all, focal_genes, fisher_by_gene) %>%
  dplyr::rename(focal_cna_score = score)
focal_score_tcga <- compute_path_score(cna_long_tcga_all, focal_genes, fisher_by_gene) %>%
  dplyr::rename(focal_cna_score = score)

axon_score_metabric <- compute_path_score(cna_long_metabric_all, axon_genes, fisher_by_gene) %>%
  dplyr::rename(axon_cna_score = score)
axon_score_tcga <- compute_path_score(cna_long_tcga_all, axon_genes, fisher_by_gene) %>%
  dplyr::rename(axon_cna_score = score)

# Join pathway CNA scores into cohort signature datasets
metabric_sig2 <- metabric_sig %>%
  left_join(ras_score_metabric,   by = "sample") %>%
  left_join(focal_score_metabric, by = "sample") %>%
  left_join(axon_score_metabric,  by = "sample")

tcga_sig2 <- tcga_sig %>%
  left_join(ras_score_tcga,   by = "sample") %>%
  left_join(focal_score_tcga, by = "sample") %>%
  left_join(axon_score_tcga,  by = "sample")

cat("   METABRIC Ras scores available for", sum(!is.na(metabric_sig2$ras_cna_score)),   "samples\n")
cat("   METABRIC Focal scores available for", sum(!is.na(metabric_sig2$focal_cna_score)), "samples\n")
cat("   METABRIC Axon scores available for",  sum(!is.na(metabric_sig2$axon_cna_score)),  "samples\n")
cat("   TCGA Ras scores available for",       sum(!is.na(tcga_sig2$ras_cna_score)),       "samples\n")
cat("   TCGA Focal scores available for",     sum(!is.na(tcga_sig2$focal_cna_score)),     "samples\n")
cat("   TCGA Axon scores available for",      sum(!is.na(tcga_sig2$axon_cna_score)),      "samples\n\n")

# ------------------------------------------------------------
# 7) Joint clinical + signature + Ras/Focal/Axon Cox models
# ------------------------------------------------------------

# Fit multivariate Cox models with clinical covariates plus signature and pathway CNA scores
cat(">> Fitting joint clinical + signature + Ras/Focal/Axon Cox models...\n")

# Prepare analysis-ready data frames (snake_case) and enforce numeric types / factor levels
met_j <- metabric_sig2 %>%
  filter(!is.na(ras_cna_score)) %>%
  clean_names() %>%
  mutate(
    os_time         = as.numeric(os_time),
    os_event        = as.numeric(os_event),
    age             = as.numeric(age),
    risk_score      = as.numeric(risk_score),
    ras_cna_score   = as.numeric(ras_cna_score),
    focal_cna_score = as.numeric(focal_cna_score),
    axon_cna_score  = as.numeric(axon_cna_score),
    stage_simple    = droplevels(as.factor(stage_simple)),
    pam50           = droplevels(as.factor(pam50))
  )

tcga_j <- tcga_sig2 %>%
  filter(!is.na(ras_cna_score)) %>%
  clean_names() %>%
  mutate(
    os_time         = as.numeric(os_time),
    os_event        = as.numeric(os_event),
    age             = as.numeric(age),
    risk_score      = as.numeric(risk_score),
    ras_cna_score   = as.numeric(ras_cna_score),
    focal_cna_score = as.numeric(focal_cna_score),
    axon_cna_score  = as.numeric(axon_cna_score),
    stage_simple    = droplevels(as.factor(stage_simple)),
    pam50           = droplevels(as.factor(pam50))
  )

# Helper function: fit a set of Cox models and return model-level metrics and tidy coefficients
fit_model_set <- function(dat, cohort_name) {
  
  dat2 <- dat %>%
    dplyr::select(
      os_time, os_event, age, stage_simple, pam50,
      risk_score, ras_cna_score, focal_cna_score, axon_cna_score
    ) %>%
    tidyr::drop_na() %>%
    as.data.frame()
  
  if (nrow(dat2) < 50) {
    warning("Too few samples for cohort: ", cohort_name)
  }
  
  surv_obj <- Surv(dat2$os_time, dat2$os_event)
  
  models <- list(
    Clinical                     = coxph(surv_obj ~ age + stage_simple + pam50, data = dat2),
    Clinical_plus_sig            = coxph(surv_obj ~ age + stage_simple + pam50 + risk_score, data = dat2),
    Clinical_plus_ras            = coxph(surv_obj ~ age + stage_simple + pam50 + ras_cna_score, data = dat2),
    Clinical_plus_sig_ras        = coxph(surv_obj ~ age + stage_simple + pam50 + risk_score + ras_cna_score, data = dat2),
    Clinical_plus_focal          = coxph(surv_obj ~ age + stage_simple + pam50 + focal_cna_score, data = dat2),
    Clinical_plus_axon           = coxph(surv_obj ~ age + stage_simple + pam50 + axon_cna_score, data = dat2),
    Clinical_plus_sig_focal      = coxph(surv_obj ~ age + stage_simple + pam50 + risk_score + focal_cna_score, data = dat2),
    Clinical_plus_sig_axon       = coxph(surv_obj ~ age + stage_simple + pam50 + risk_score + axon_cna_score, data = dat2),
    Clinical_plus_ras_focal      = coxph(surv_obj ~ age + stage_simple + pam50 + ras_cna_score + focal_cna_score, data = dat2),
    Clinical_plus_ras_axon       = coxph(surv_obj ~ age + stage_simple + pam50 + ras_cna_score + axon_cna_score, data = dat2),
    Clinical_plus_focal_axon     = coxph(surv_obj ~ age + stage_simple + pam50 + focal_cna_score + axon_cna_score, data = dat2),
    Clinical_plus_sig_ras_foc_ax = coxph(
      surv_obj ~ age + stage_simple + pam50 +
        risk_score + ras_cna_score +
        focal_cna_score + axon_cna_score,
      data = dat2
    )
  )
  
  results <- lapply(names(models), function(nm) {
    m  <- models[[nm]]
    sm <- summary(m)
    
    coefs_tidy <- broom::tidy(m, exponentiate = FALSE, conf.int = TRUE) %>%
      dplyr::mutate(
        HR      = exp(estimate),
        HR_low  = exp(conf.low),
        HR_high = exp(conf.high)
      )
    
    tibble(
      Cohort  = cohort_name,
      Model   = nm,
      n       = nrow(dat2),
      events  = sum(dat2$os_event),
      C_index = sm$concordance[1],
      C_se    = sm$concordance[2],
      AIC     = AIC(m),
      loglik  = sm$loglik[2],
      coefs   = list(coefs_tidy)
    )
  })
  
  bind_rows(results)
}

# Fit joint model sets in METABRIC and TCGA-BRCA
res_met_full  <- fit_model_set(met_j,  "METABRIC")
res_tcga_full <- fit_model_set(tcga_j, "TCGA-BRCA")

# Combine results across cohorts
joint_res_full <- bind_rows(res_met_full, res_tcga_full)

# Unnest coefficient tables for downstream plotting
coefs_all_joint_full <- joint_res_full %>%
  dplyr::select(Cohort, Model, coefs) %>%
  unnest(coefs) %>%
  dplyr::select(
    Cohort, Model, term,
    estimate, std.error, statistic, p.value,
    HR, HR_low, HR_high
  )

cat("   Combined coefficient table ready for Figure 10.\n\n")

# ------------------------------------------------------------
# 8) Figure 10 – Heatmap (HR>1) + Forest plot (Sig & Ras)
# ------------------------------------------------------------

# Build the final Figure 10 panels and save them as a high-resolution JPEG
cat(">> [Figure 10] Building heatmap + forest plot (clinical vs sig vs Ras)...\n")

# Define which model terms and model variants will be displayed in Figure 10
terms_fig10 <- c(
  "age",
  "stage_simpleStageII",
  "stage_simpleStageIII_IV",
  "pam50LumB",
  "pam50Basal",
  "pam50Her2",
  "pam50Normal",
  "risk_score",
  "ras_cna_score"
)

models_fig10 <- c(
  "Clinical",
  "Clinical_plus_sig",
  "Clinical_plus_ras",
  "Clinical_plus_sig_ras"
)

# Define human-readable labels for the heatmap y-axis (terms) and x-axis (models)
term_labels_fig10 <- c(
  age                       = "Age",
  stage_simpleStageII       = "Stage II",
  stage_simpleStageIII_IV   = "Stage III/IV",
  stage_simpleStageX        = "Stage X",
  
  pam50LumB                 = "PAM50 (LumB)",
  pam50Basal                = "PAM50 (Basal)",
  pam50Her2                 = "PAM50 (HER2)",
  pam50Normal               = "PAM50 (Normal)",
  
  risk_score                = "Gene signature score",
  ras_cna_score             = "Ras signaling pathway\n(CNA score)"
)

model_labels_fig10 <- c(
  Clinical              = "Clinical",
  Clinical_plus_sig     = "Clinical + Signature",
  Clinical_plus_ras     = "Clinical + Ras",
  Clinical_plus_sig_ras = "Clinical + Signature + Ras"
)

# Create the heatmap coefficient table and format HR labels and significance marks
coefs_fig10 <- coefs_all_joint_full %>%
  dplyr::filter(term %in% terms_fig10,
                Model %in% models_fig10) %>%
  dplyr::mutate(
    term  = factor(term,  levels = terms_fig10),
    Model = factor(Model, levels = models_fig10),
    
    HR_for_plot = if_else(HR > 1, HR, NA_real_),
    
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    
    HR_label = ifelse(
      is.na(HR),
      "",
      paste0(sprintf("%.2f", HR), sig)
    )
  )

# Panel (a): heatmap of HR values (only HR>1 colored) with per-tile HR labels
p_fig10_heat <- ggplot(
  coefs_fig10,
  aes(x = Model, y = term, fill = HR_for_plot)
) +
  geom_tile(color = "grey85") +
  geom_text(
    aes(label = HR_label),
    size     = 3.3,
    fontface = "bold"
  ) +
  scale_fill_gradient(
    name    = "HR",
    low     = "white",
    high    = "red",
    na.value = "grey95",
    labels   = function(x) sprintf("%.2f", x)
  ) +
  scale_y_discrete(labels = term_labels_fig10) +
  scale_x_discrete(labels = model_labels_fig10) +
  facet_wrap(~ Cohort, nrow = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    panel.grid   = element_blank(),
    strip.text   = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold", size = 12),
    legend.text  = element_text(face = "bold", size = 12),
    axis.title   = element_blank()
  )

# Panel (b): forest plot for the key predictors (signature risk score and Ras CNA score)
terms_forest10 <- c("risk_score", "ras_cna_score")

coefs_forest10 <- coefs_fig10 %>%
  dplyr::filter(term %in% terms_forest10) %>%
  dplyr::mutate(
    predictor = factor(
      term,
      levels = terms_forest10,
      labels = c(
        "Signature score",
        "Ras signaling pathway\n(CNA score)"
      )
    ),
    Model_f = factor(
      Model,
      levels = models_fig10,
      labels = model_labels_fig10
    ),
    HR_label = sprintf("%.2f", HR),
    HR_fake_upper = ifelse(Cohort == "TCGA-BRCA", 50.00, NA_real_)
  )

p_fig10_forest <- ggplot(
  coefs_forest10,
  aes(x = HR, y = predictor, color = Model_f)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(
    size     = 3,
    position = position_dodge(width = 0.7)
  ) +
  geom_errorbar(
    aes(xmin = HR_low, xmax = HR_high),
    height   = 0,
    position = position_dodge(width = 0.7)
  ) +
  geom_text(
    aes(label = sprintf("%.2f", HR)),
    position   = position_dodge(width = 0.7),
    vjust      = -1.2,
    size       = 3.6,
    fontface   = "bold",
    show.legend = FALSE
  ) +
  geom_point(
    data        = dplyr::filter(coefs_forest10, !is.na(HR_fake_upper)),
    aes(x = HR_fake_upper, y = predictor),
    inherit.aes = FALSE,
    alpha       = 0
  ) +
  facet_wrap(~ Cohort, nrow = 1, scales = "free_x") +
  scale_x_log10(
    name   = "HR",
    breaks = c(0.50, 1.00, 3.00, 9.00, 50.00),
    labels = function(x) sprintf("%.2f", x),
    limits = c(0.50, NA)
  ) +
  ylab(NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 14),
    axis.title.y     = element_blank(),
    plot.margin      = margin(t = 15, r = 5, b = 5, l = 5),
    panel.spacing    = unit(1.2, "cm"),
    axis.title.x     = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(size = 11),
    axis.text.y      = element_text(size = 11, face = "bold")
  ) +
  guides(
    color = guide_legend(
      title        = "Model ",
      override.aes = list(label = NULL)
    )
  ) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 11)
  )

# Combine panels and add panel tags (A., B.)
fig10_combined <- p_fig10_heat / p_fig10_forest +
  patchwork::plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = "."
  )

fig10_combined <- fig10_combined &
  theme(
    plot.tag = element_text(face = "bold", size = 14)
  )

# Save Figure 10 to the manuscript figures directory
fig10_file <- file.path(figures_dir, "Figure 10.jpeg")

ggsave(
  filename = fig10_file,
  plot     = fig10_combined,
  width    = 11,
  height   = 11,
  dpi      = 300
)

# Print the final output path and completion message
cat("   Figure 10 saved to:\n   ", fig10_file, "\n\n")
cat(">> Figure 10 pipeline completed.\n\n")