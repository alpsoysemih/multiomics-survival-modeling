# Load core tidyverse-style packages for I/O, wrangling, and reshaping
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

# Load annotation and survival modeling packages
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(survival))

# Load plotting and layout utilities
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

# Load model tidying and tibble helpers
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(tibble))

# Keep character columns as characters by default
options(stringsAsFactors = FALSE)

# Define input directories used throughout the script
metabric_dir   <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir       <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
gdsc_dir       <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
pathfindR_dir  <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"

# Create output directory for manuscript figures if it does not exist
figures_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}

# Print the main output location
cat(">> All outputs will be written under:\n   ", figures_dir, "\n\n")

# Define signature genes used for expression and CNA-based features
sig_gene_expr <- "SERPINE1"
sig_genes_cna <- c("CLDN1", "CXCR4", "LYPD6B")
sig_genes_all <- unique(c(sig_gene_expr, sig_genes_cna))

# Define target pathway names (must match pathfindR Term_Description)
ras_pathway_name    <- "Ras signaling pathway"
focal_pathway_name  <- "Focal adhesion"
axon_pathway_name   <- "Axon guidance"

# Load GDSC data and compute per-gene Fisher ORs to determine CNA directionality
cat(">> [1] Loading GDSC data and computing Fisher OR for CNA...\n")

# Define GDSC input files (expression, CNA, and response)
gdsc_expr_file     <- file.path(gdsc_dir, "GDSC_exprs.Paclitaxel.eb_with.TCGA_exprs.Paclitaxel.tsv")
gdsc_cna_file      <- file.path(gdsc_dir, "GDSC_CNA.Paclitaxel.tsv")
gdsc_response_file <- file.path(gdsc_dir, "GDSC_response.Paclitaxel.tsv")

# Read GDSC input tables
gdsc_expr     <- read_tsv(gdsc_expr_file,     show_col_types = FALSE)
gdsc_cna      <- read_tsv(gdsc_cna_file,      show_col_types = FALSE)
gdsc_response <- read_tsv(gdsc_response_file, show_col_types = FALSE)

# Binarize CNA calls (any non-zero CNA becomes 1)
gdsc_cna[, -1] <- ifelse(gdsc_cna[, -1] != 0, 1, 0)

# Map Entrez IDs to gene symbols for expression and CNA matrices
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

# Harmonize column names to a shared Gene_Symbol field
colnames(gdsc_expr)[colnames(gdsc_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(gdsc_cna)[colnames(gdsc_cna)   == "gene_id"]  <- "Gene_Symbol"

# Restrict CNA matrix to samples that have response labels
resp <- gdsc_response %>% mutate(status = response)
gdsc_samples <- as.character(unique(resp$sample_name))

# Subset CNA matrix to response-annotated samples only
cna_gdsc <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

# Convert CNA table into a binary matrix with gene symbols as rownames
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

# Compute Fisher odds ratio and p-value per gene (CNA vs response)
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

# Report the number of genes with computed Fisher ORs
cat("   Fisher OR computed for", nrow(fisher_by_gene), "genes\n\n")

# Load METABRIC expression, CNA, and overall survival tables
cat(">> [2] Loading METABRIC expression, CNA and OS data...\n")

# Define METABRIC input files
metabric_expr_file     <- file.path(metabric_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
metabric_cna_file      <- file.path(metabric_dir, "data_cna.txt")
metabric_survival_file <- file.path(metabric_dir, "brca_metabric_clinical_data.tsv")

# Read METABRIC expression, CNA, and clinical survival data
metabric_expr     <- read_tsv(metabric_expr_file,     show_col_types = FALSE)
metabric_cna      <- read_tsv(metabric_cna_file,      show_col_types = FALSE)
metabric_survival <- read_tsv(metabric_survival_file, show_col_types = FALSE)

# Standardize gene identifier column names for expression and CNA tables
metabric_expr <- metabric_expr %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

metabric_cna <- metabric_cna %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

# Create a binarized CNA matrix for downstream pathway scoring
metabric_cna_binary <- metabric_cna
metabric_cna_binary[, -1] <- ifelse(metabric_cna_binary[, -1] != 0, 1, 0)

# Extract METABRIC overall survival time (months) and event indicator
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

# Print METABRIC OS sample size
cat("   METABRIC OS samples:", nrow(metabric_surv_small), "\n\n")

# Load TCGA-BRCA expression, CNA, and overall survival tables
cat(">> [3] Loading TCGA-BRCA expression, CNA and OS data...\n")

# Define TCGA-BRCA input files
tcga_expr_file     <- file.path(tcga_dir, "TCGA-BRCA_exprs.z.tsv")
tcga_cna_file      <- file.path(tcga_dir, "TCGA-BRCA_CNA.tsv")
tcga_survival_file <- file.path(tcga_dir, "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv")

# Read TCGA-BRCA expression, CNA, and clinical survival data
tcga_expr     <- read_tsv(tcga_expr_file,     show_col_types = FALSE)
tcga_cna      <- read_tsv(tcga_cna_file,      show_col_types = FALSE)
tcga_survival <- read_tsv(tcga_survival_file, show_col_types = FALSE)

# Map Entrez identifiers to gene symbols for TCGA expression and CNA matrices
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

# Harmonize column names to a shared Gene_Symbol field
colnames(tcga_expr)[colnames(tcga_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(tcga_cna)[colnames(tcga_cna)   == "gene"]     <- "Gene_Symbol"

# Create a binarized CNA matrix for downstream pathway scoring
tcga_cna_binary <- tcga_cna
tcga_cna_binary[, -1] <- ifelse(tcga_cna_binary[, -1] != 0, 1, 0)

# Extract TCGA overall survival time (months) and event indicator
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

# Print TCGA OS sample size
cat("   TCGA-BRCA OS samples:", nrow(tcga_surv_small), "\n\n")

# Build cohort-specific signature datasets (expression + aligned CNA features)
cat(">> [4] Building molecular signature datasets (METABRIC & TCGA)...\n")

# Convert METABRIC signature gene expression matrix to long format
expr_long_metabric_sig <- metabric_expr %>%
  filter(Gene_Symbol %in% sig_genes_all) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

# Convert METABRIC signature gene CNA matrix to long format
cna_long_metabric_sig <- metabric_cna_binary %>%
  filter(Gene_Symbol %in% sig_genes_cna) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

# Align METABRIC CNA direction per gene using GDSC-derived Fisher ORs
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

# Extract SERPINE1 expression and aligned CNA features into wide format (METABRIC)
serp_metabric <- expr_long_metabric_sig %>%
  filter(Gene_Symbol == sig_gene_expr) %>%
  transmute(sample, SERPINE1_expr = as.numeric(expr))

cldn1_metabric <- cna_long_metabric_aligned %>%
  filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = as.numeric(cna_aligned))

cxcr4_metabric <- cna_long_metabric_aligned %>%
  filter(Gene_Symbol == "CXCR4") %>%
  transmute(sample, CXCR4_CNA = as.numeric(cna_aligned))

lypd6b_metabric <- cna_long_metabric_aligned %>%
  filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = as.numeric(cna_aligned))

# Join METABRIC survival with signature expression/CNA features and keep complete cases
metabric_sig <- metabric_surv_small %>%
  left_join(serp_metabric,  by = "sample") %>%
  left_join(cldn1_metabric, by = "sample") %>%
  left_join(cxcr4_metabric, by = "sample") %>%
  left_join(lypd6b_metabric, by = "sample") %>%
  filter(
    !is.na(SERPINE1_expr),
    !is.na(CLDN1_CNA),
    !is.na(CXCR4_CNA),
    !is.na(LYPD6B_CNA)
  )

# Report METABRIC signature dataset size
cat("   METABRIC signature samples:", nrow(metabric_sig), "\n")

# Convert TCGA signature gene expression matrix to long format (barcode -> first 15 chars)
expr_long_tcga_sig <- tcga_expr %>%
  filter(Gene_Symbol %in% sig_genes_all) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "expr"
  ) %>%
  mutate(sample = substr(sample_barcode, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, expr)

# Convert TCGA signature gene CNA matrix to long format (barcode -> first 15 chars)
cna_long_tcga_sig <- tcga_cna_binary %>%
  filter(Gene_Symbol %in% sig_genes_cna) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  mutate(sample = substr(sample_barcode, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, cna)

# Align TCGA CNA direction per gene using GDSC-derived Fisher ORs
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

# Extract SERPINE1 expression and aligned CNA features into wide format (TCGA)
serp_tcga <- expr_long_tcga_sig %>%
  filter(Gene_Symbol == sig_gene_expr) %>%
  transmute(sample, SERPINE1_expr = as.numeric(expr))

cldn1_tcga <- cna_long_tcga_aligned %>%
  filter(Gene_Symbol == "CLDN1") %>%
  transmute(sample, CLDN1_CNA = as.numeric(cna_aligned))

cxcr4_tcga <- cna_long_tcga_aligned %>%
  filter(Gene_Symbol == "CXCR4") %>%
  transmute(sample, CXCR4_CNA = as.numeric(cna_aligned))

lypd6b_tcga <- cna_long_tcga_aligned %>%
  filter(Gene_Symbol == "LYPD6B") %>%
  transmute(sample, LYPD6B_CNA = as.numeric(cna_aligned))

# Join TCGA survival with signature expression/CNA features and keep complete cases
tcga_sig <- tcga_surv_small %>%
  left_join(serp_tcga,  by = "sample") %>%
  left_join(cldn1_tcga, by = "sample") %>%
  left_join(cxcr4_tcga, by = "sample") %>%
  left_join(lypd6b_tcga, by = "sample") %>%
  filter(
    !is.na(SERPINE1_expr),
    !is.na(CLDN1_CNA),
    !is.na(CXCR4_CNA),
    !is.na(LYPD6B_CNA)
  )

# Report TCGA signature dataset size
cat("   TCGA signature samples:", nrow(tcga_sig), "\n\n")

# Load pathfindR results and compute aligned CNA pathway scores for Ras/Focal/Axon
cat(">> [5] Loading pathfindR enrichment and computing pathway CNA scores...\n")

# Read the pathfindR enrichment file used for pathway gene definitions
pathfindR_file <- file.path(pathfindR_dir, "pathfindR_paclitaxel.tsv")
enrichment_df  <- read_tsv(pathfindR_file, show_col_types = FALSE)

# Parse pathway gene members from up- and down-regulated gene lists
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

# Subset the pathway-to-gene map to the three target pathways
ras_genes <- path_genes_long %>% filter(pathway == ras_pathway_name)
focal_genes <- path_genes_long %>% filter(pathway == focal_pathway_name)
axon_genes  <- path_genes_long %>% filter(pathway == axon_pathway_name)

# Validate that the requested pathways exist in the pathfindR file
if (nrow(ras_genes) == 0)  stop("No genes found for Ras signaling pathway in pathfindR file.")
if (nrow(focal_genes) == 0) warning("No genes found for Focal adhesion pathway in pathfindR file.")
if (nrow(axon_genes)  == 0) warning("No genes found for Axon guidance pathway in pathfindR file.")

# Print the number of unique genes used per pathway
cat("   Ras genes:",   length(unique(ras_genes$Gene_Symbol)),   "\n")
cat("   Focal genes:", length(unique(focal_genes$Gene_Symbol)), "\n")
cat("   Axon genes:",  length(unique(axon_genes$Gene_Symbol)),  "\n\n")

# Compute a sample-level pathway CNA score by averaging aligned CNA across pathway genes
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

# Convert METABRIC CNA matrix to long format for pathway scoring
cna_long_metabric_all <- metabric_cna_binary %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

# Convert TCGA CNA matrix to long format for pathway scoring (barcode -> first 15 chars)
cna_long_tcga_all <- tcga_cna_binary %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  mutate(sample = substr(sample_barcode, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, cna)

# Compute Ras pathway CNA scores for METABRIC and TCGA
ras_score_metabric <- compute_path_score(cna_long_metabric_all, ras_genes, fisher_by_gene) %>%
  dplyr::rename(ras_cna_score = score)
ras_score_tcga <- compute_path_score(cna_long_tcga_all, ras_genes, fisher_by_gene) %>%
  dplyr::rename(ras_cna_score = score)

# Compute Focal adhesion pathway CNA scores for METABRIC and TCGA
focal_score_metabric <- compute_path_score(cna_long_metabric_all, focal_genes, fisher_by_gene) %>%
  dplyr::rename(focal_cna_score = score)
focal_score_tcga <- compute_path_score(cna_long_tcga_all, focal_genes, fisher_by_gene) %>%
  dplyr::rename(focal_cna_score = score)

# Compute Axon guidance pathway CNA scores for METABRIC and TCGA
axon_score_metabric <- compute_path_score(cna_long_metabric_all, axon_genes, fisher_by_gene) %>%
  dplyr::rename(axon_cna_score = score)
axon_score_tcga <- compute_path_score(cna_long_tcga_all, axon_genes, fisher_by_gene) %>%
  dplyr::rename(axon_cna_score = score)

# Attach pathway CNA scores to the signature datasets
metabric_sig2 <- metabric_sig %>%
  left_join(ras_score_metabric,   by = "sample") %>%
  left_join(focal_score_metabric, by = "sample") %>%
  left_join(axon_score_metabric,  by = "sample")

tcga_sig2 <- tcga_sig %>%
  left_join(ras_score_tcga,   by = "sample") %>%
  left_join(focal_score_tcga, by = "sample") %>%
  left_join(axon_score_tcga,  by = "sample")

# Report availability of pathway CNA scores in each cohort
cat("   METABRIC Ras scores available for", sum(!is.na(metabric_sig2$ras_cna_score)),   "samples\n")
cat("   METABRIC Focal scores available for", sum(!is.na(metabric_sig2$focal_cna_score)), "samples\n")
cat("   METABRIC Axon scores available for",  sum(!is.na(metabric_sig2$axon_cna_score)),  "samples\n")
cat("   TCGA Ras scores available for",       sum(!is.na(tcga_sig2$ras_cna_score)),       "samples\n")
cat("   TCGA Focal scores available for",     sum(!is.na(tcga_sig2$focal_cna_score)),     "samples\n")
cat("   TCGA Axon scores available for",      sum(!is.na(tcga_sig2$axon_cna_score)),      "samples\n\n")

# Build Figure 9: HR heatmap and forest plot for signature/pathway CNA models (no clinical covariates)
cat(">> [6] Building Figure 9 heatmap + forest plot (signature & Ras/Focal/Axon CNA, no clinical covariates)...\n")

# Fit the 3-feature signature Cox model in METABRIC to obtain coefficients for risk scoring
sig_fit_3genes <- coxph(
  Surv(OS_TIME, OS_EVENT) ~ SERPINE1_expr + CLDN1_CNA + LYPD6B_CNA,
  data = metabric_sig2
)

# Print fitted signature coefficients for reference
cat("   Signature (3-gene) coefficients (METABRIC):\n")
print(summary(sig_fit_3genes)$coefficients[, c("coef", "exp(coef)")])

# Compute signature linear predictor scores in METABRIC and project them to TCGA
metabric_sig2$sig_score_3g <- predict(
  sig_fit_3genes,
  newdata = metabric_sig2,
  type    = "lp"
)
tcga_sig2$sig_score_3g <- predict(
  sig_fit_3genes,
  newdata = tcga_sig2,
  type    = "lp"
)

# Prepare complete-case datasets for Cox models without clinical covariates
prep_minimal_dat <- function(dat, cohort_name) {
  dat2 <- dat %>%
    dplyr::select(
      OS_TIME, OS_EVENT,
      sig_score_3g,
      ras_cna_score,
      focal_cna_score,
      axon_cna_score
    ) %>%
    tidyr::drop_na() %>%
    as.data.frame()
  
  if (nrow(dat2) < 30) {
    warning("Too few samples for minimal models in cohort: ", cohort_name)
  }
  
  dat2
}

# Build minimal modeling datasets for METABRIC and TCGA-BRCA
met_min  <- prep_minimal_dat(metabric_sig2, "METABRIC")
tcga_min <- prep_minimal_dat(tcga_sig2,    "TCGA-BRCA")

# Fit multiple Cox models (signature and pathway CNA scores) without clinical covariates and extract HRs
fit_noclin_sig_models <- function(dat2, cohort_name) {
  
  surv_obj <- Surv(dat2$OS_TIME, dat2$OS_EVENT)
  
  models <- list(
    Sig_only                = coxph(surv_obj ~ sig_score_3g, data = dat2),
    Ras_only                = coxph(surv_obj ~ ras_cna_score, data = dat2),
    Axon_only               = coxph(surv_obj ~ axon_cna_score, data = dat2),
    Focal_only              = coxph(surv_obj ~ focal_cna_score, data = dat2),
    Sig_plus_Ras            = coxph(surv_obj ~ sig_score_3g + ras_cna_score, data = dat2),
    Sig_plus_Axon           = coxph(surv_obj ~ sig_score_3g + axon_cna_score, data = dat2),
    Sig_plus_Focal          = coxph(surv_obj ~ sig_score_3g + focal_cna_score, data = dat2),
    Sig_plus_Ras_Axon       = coxph(surv_obj ~ sig_score_3g + ras_cna_score + axon_cna_score,  data = dat2),
    Sig_plus_Ras_Focal      = coxph(surv_obj ~ sig_score_3g + ras_cna_score + focal_cna_score, data = dat2),
    Sig_plus_Axon_Focal     = coxph(surv_obj ~ sig_score_3g + axon_cna_score + focal_cna_score, data = dat2),
    Sig_plus_Ras_Focal_Axon = coxph(
      surv_obj ~ sig_score_3g + ras_cna_score +
        focal_cna_score + axon_cna_score,
      data = dat2
    )
  )
  
  par_terms <- c("sig_score_3g", "ras_cna_score", "axon_cna_score", "focal_cna_score")
  
  res_list <- lapply(names(models), function(nm) {
    fit <- models[[nm]]
    tt  <- broom::tidy(fit, exponentiate = FALSE, conf.int = TRUE)
    
    tt2 <- tt %>%
      dplyr::filter(term %in% par_terms) %>%
      dplyr::mutate(
        HR      = exp(estimate),
        HR_low  = exp(conf.low),
        HR_high = exp(conf.high)
      ) %>%
      dplyr::select(
        term, estimate, std.error, statistic, p.value,
        HR, HR_low, HR_high
      )
    
    if (nrow(tt2) == 0) {
      return(NULL)
    }
    
    tt2 %>%
      dplyr::mutate(
        Cohort = cohort_name,
        Model  = nm
      )
  })
  
  dplyr::bind_rows(res_list)
}

# Fit and combine model coefficients across cohorts for Figure 9
coefs_met_fig9  <- fit_noclin_sig_models(met_min,  "METABRIC")
coefs_tcga_fig9 <- fit_noclin_sig_models(tcga_min, "TCGA-BRCA")

coefs_fig9_all <- dplyr::bind_rows(coefs_met_fig9, coefs_tcga_fig9)

# Define terms, model order, and display labels used in Figure 9
terms_fig9 <- c(
  "sig_score_3g",
  "ras_cna_score",
  "axon_cna_score",
  "focal_cna_score"
)

models_fig9 <- c(
  "Sig_only",
  "Ras_only",
  "Axon_only",
  "Focal_only",
  "Sig_plus_Ras",
  "Sig_plus_Axon",
  "Sig_plus_Focal",
  "Sig_plus_Ras_Axon",
  "Sig_plus_Ras_Focal",
  "Sig_plus_Axon_Focal",
  "Sig_plus_Ras_Focal_Axon"
)

term_labels_fig9 <- c(
  sig_score_3g   = "Gene signature score",
  ras_cna_score  = "Ras signaling pathway\n(CNA score)",
  axon_cna_score = "Axon guidance (CNA score)",
  focal_cna_score= "Focal adhesion (CNA score)"
)

model_labels_fig9 <- c(
  Sig_only                = "Signature only",
  Ras_only                = "Ras only",
  Axon_only               = "Axon only",
  Focal_only              = "Focal only",
  Sig_plus_Ras            = "Signature + Ras",
  Sig_plus_Axon           = "Signature + Axon",
  Sig_plus_Focal          = "Signature + Focal",
  Sig_plus_Ras_Axon       = "Signature + Ras + Axon",
  Sig_plus_Ras_Focal      = "Signature + Ras + Focal",
  Sig_plus_Axon_Focal     = "Signature + Axon + Focal",
  Sig_plus_Ras_Focal_Axon = "Signature + Ras + Focal + Axon"
)

# Format coefficients for the HR heatmap panel (color only for HR > 1 and add significance stars)
coefs_fig9_heat <- coefs_fig9_all %>%
  dplyr::filter(term %in% terms_fig9,
                Model %in% models_fig9) %>%
  dplyr::mutate(
    Cohort = factor(Cohort, levels = c("METABRIC", "TCGA-BRCA")),
    term   = factor(term,  levels = terms_fig9),
    Model  = factor(Model, levels = models_fig9),
    HR_for_plot = dplyr::if_else(HR > 1, HR, NA_real_),
    sig = dplyr::case_when(
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

# Build panel (a): heatmap of HR values with in-tile labels
p_fig9_heat <- ggplot(
  coefs_fig9_heat,
  aes(x = Model, y = term, fill = HR_for_plot)
) +
  geom_tile(color = "grey85") +
  geom_text(aes(label = HR_label),
            size = 3.5,
            fontface = "bold") +
  scale_fill_gradient(
    name    = "HR",
    low     = "white",
    high    = "red",
    na.value = "grey95",
    labels = function(x) sprintf("%.2f", x)
  ) +
  scale_y_discrete(labels = term_labels_fig9) +
  scale_x_discrete(labels = model_labels_fig9) +
  facet_wrap(~ Cohort, nrow = 1) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, face = "bold", size = 11),
    axis.text.y  = element_text(face = "bold", size = 11),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    strip.text   = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text  = element_text(face = "bold", size = 11)
  )

# Prepare coefficients for panel (b): forest plot of HR and 95% CI
coefs_fig9_forest <- coefs_fig9_all %>%
  dplyr::filter(term %in% terms_fig9,
                Model %in% models_fig9) %>%
  dplyr::mutate(
    Cohort = factor(Cohort, levels = c("METABRIC", "TCGA-BRCA")),
    predictor = factor(
      term,
      levels = terms_fig9,
      labels = term_labels_fig9
    ),
    Model = factor(
      Model,
      levels = models_fig9,
      labels = model_labels_fig9
    )
  )

# Build panel (b): forest plot with model-specific colors and log-scale x-axis
p_fig9_forest <- ggplot(
  coefs_fig9_forest,
  aes(x = HR, y = predictor, color = Model)
) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_point(
    size     = 2.4,
    position = position_dodge(width = 0.7)
  ) +
  geom_errorbar(
    aes(xmin = HR_low, xmax = HR_high),
    height   = 0,
    position = position_dodge(width = 0.7)
  ) +
  scale_x_log10(
    name   = "HR",
    labels = function(x) sprintf("%.2f", x)
  ) +
  facet_wrap(~ Cohort, nrow = 1) +
  ylab(NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold", size = 12),
    legend.text      = element_text(size = 11),
    strip.text       = element_text(face = "bold", size = 14),
    axis.text.x      = element_text(size = 11),
    axis.text.y      = element_text(face = "bold", size = 11),
    axis.title.x     = element_text(face = "bold")
  )

# Combine heatmap and forest panels into a single Figure 9 with bold panel tags
combined_fig9 <- (p_fig9_heat / p_fig9_forest) +
  patchwork::plot_layout(heights = c(1, 1.1)) +
  patchwork::plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  ) &
  labs(fill = "Models") &
  theme(
    plot.tag = element_text(face = "bold", size = 16)
  )

# Save the combined Figure 9 as a high-resolution JPEG
fig9_file <- file.path(figures_dir, "Figure 9.jpeg")

ggsave(
  filename = fig9_file,
  plot     = combined_fig9,
  width    = 16,
  height   = 10,
  dpi      = 300
)

# Print the output file path for Figure 9
cat("   Figure 9 (panels a: heatmap, b: forest) saved to:\n   ", fig9_file, "\n\n")
cat(">> [6] Done.\n\n")