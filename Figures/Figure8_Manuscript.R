# Load required packages quietly (one by one) to keep logs clean
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))

# Use character strings by default instead of factors
options(stringsAsFactors = FALSE)

## ------------------------------------------------------------
## 0) Directories
## ------------------------------------------------------------

# Define dataset, intermediate results, and figure output directories
metabric_data_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_data_dir         <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
gdsc_dir              <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
pathfindR_dir         <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"

# Define cohort-specific results directories for pathway-level Cox outputs
metabric_results_dir  <- "/Volumes/Expansion/Prognostic_Analysis/Results/METABRIC"
tcga_results_dir      <- "/Volumes/Expansion/Prognostic_Analysis/Results/TCGA-BRCA"

# Create the manuscript figures directory if missing and print the save location
figures_dir           <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}
cat(">> Figures will be saved to:\n   ", figures_dir, "\n\n")

## ------------------------------------------------------------
## 1) Clinical OS data (METABRIC & TCGA-BRCA)
## ------------------------------------------------------------

# Load and preprocess METABRIC overall survival and clinical covariates
cat(">> Loading METABRIC clinical OS data...\n")

metabric_survival_file <- file.path(metabric_data_dir, "brca_metabric_clinical_data.tsv")
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

# Drop very rare stage categories to avoid unstable modeling strata
stage_counts <- table(metabric_surv_small$Stage_simple)
rare_stages  <- names(stage_counts[stage_counts < 5])
if (length(rare_stages) > 0) {
  metabric_surv_small <- metabric_surv_small %>%
    dplyr::filter(!Stage_simple %in% rare_stages) %>%
    dplyr::mutate(Stage_simple = droplevels(Stage_simple))
}

cat("   METABRIC OS samples:", nrow(metabric_surv_small), "\n\n")

# Load and preprocess TCGA-BRCA overall survival and clinical covariates
cat(">> Loading TCGA-BRCA clinical OS data...\n")

tcga_survival_file <- file.path(tcga_data_dir, "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv")
tcga_survival <- read_tsv(tcga_survival_file, show_col_types = FALSE)

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


# Drop very rare stage categories to avoid unstable modeling strata
stage_counts_tcga <- table(tcga_surv_small$Stage_simple)
rare_stages_tcga  <- names(stage_counts_tcga[stage_counts_tcga < 5])
if (length(rare_stages_tcga) > 0) {
  tcga_surv_small <- tcga_surv_small %>%
    dplyr::filter(!Stage_simple %in% rare_stages_tcga) %>%
    dplyr::mutate(Stage_simple = droplevels(Stage_simple))
}

cat("   TCGA-BRCA OS samples:", nrow(tcga_surv_small), "\n\n")

## ------------------------------------------------------------
## 2) GDSC: Fisher OR for CNA direction (fisher_by_gene)
## ------------------------------------------------------------

# Load GDSC data and compute gene-wise CNA direction using Fisher's exact test
cat(">> Loading GDSC data and computing Fisher OR for CNA...\n")

gdsc_expr_file     <- file.path(gdsc_dir, "GDSC_exprs.Paclitaxel.eb_with.TCGA_exprs.Paclitaxel.tsv")
gdsc_cna_file      <- file.path(gdsc_dir, "GDSC_CNA.Paclitaxel.tsv")
gdsc_response_file <- file.path(gdsc_dir, "GDSC_response.Paclitaxel.tsv")

gdsc_expr     <- read_tsv(gdsc_expr_file,     show_col_types = FALSE)
gdsc_cna      <- read_tsv(gdsc_cna_file,      show_col_types = FALSE)
gdsc_response <- read_tsv(gdsc_response_file, show_col_types = FALSE)

# Binarize CNA calls (0 = no event, 1 = any non-zero CNA event)
gdsc_cna[, -1] <- ifelse(gdsc_cna[, -1] != 0, 1, 0)

# Map Entrez gene IDs to HGNC symbols for expression and CNA matrices
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

# Standardize gene identifier column names across GDSC tables
colnames(gdsc_expr)[colnames(gdsc_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(gdsc_cna)[colnames(gdsc_cna)   == "gene_id"]  <- "Gene_Symbol"

# Prepare sample-level response status for CNA-response association tests
resp <- gdsc_response %>%
  mutate(status = response)

gdsc_samples <- as.character(unique(resp$sample_name))

# Subset CNA matrix to samples present in the response table
cna_gdsc <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

# Build a binary CNA matrix with gene symbols as rownames
cna_mat <- as.matrix(cna_gdsc[, -1])
rownames(cna_mat) <- cna_gdsc$Gene_Symbol
cna_bin    <- ifelse(cna_mat != 0, 1, 0)
cna_bin_df <- as_tibble(cna_bin, rownames = "Gene_Symbol")

# Convert binary CNA matrix to long format and join response labels
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

# Compute Fisher OR and p-value per gene (CNA vs response) to define CNA direction
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

cat("   GDSC Fisher OR computed for", nrow(fisher_by_gene), "genes\n\n")

## ------------------------------------------------------------
## 3) METABRIC & TCGA CNA matrices (binary, long format)
## ------------------------------------------------------------

# Load and binarize METABRIC CNA matrix, then convert to long format
cat(">> Loading METABRIC CNA data...\n")

metabric_cna_file <- file.path(metabric_data_dir, "data_cna.txt")
metabric_cna <- read_tsv(metabric_cna_file, show_col_types = FALSE)

metabric_cna <- metabric_cna %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

metabric_cna_binary <- metabric_cna
metabric_cna_binary[, -1] <- ifelse(metabric_cna_binary[, -1] != 0, 1, 0)

cna_long_metabric_all <- metabric_cna_binary %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

cat("   METABRIC CNA long rows:", nrow(cna_long_metabric_all), "\n\n")

# Load and binarize TCGA-BRCA CNA matrix, map Entrez IDs, and convert to long format
cat(">> Loading TCGA-BRCA CNA data...\n")

tcga_cna_file <- file.path(tcga_data_dir, "TCGA-BRCA_CNA.tsv")
tcga_cna <- read_tsv(tcga_cna_file, show_col_types = FALSE)

tcga_cna$gene <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_cna$gene),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)
colnames(tcga_cna)[colnames(tcga_cna) == "gene"] <- "Gene_Symbol"

tcga_cna_binary <- tcga_cna
tcga_cna_binary[, -1] <- ifelse(tcga_cna_binary[, -1] != 0, 1, 0)

cna_long_tcga_all <- tcga_cna_binary %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  mutate(sample = substr(sample_barcode, 1, 15)) %>%
  dplyr::select(Gene_Symbol, sample, cna)

cat("   TCGA-BRCA CNA long rows:", nrow(cna_long_tcga_all), "\n\n")

## ------------------------------------------------------------
## 4) Pathway ↔ gene mapping from pathfindR
## ------------------------------------------------------------

# Load pathfindR enrichment output and build a pathway-to-gene membership table
cat(">> Loading pathfindR enrichment and building pathway–gene map...\n")

pathfindR_file <- file.path(pathfindR_dir, "pathfindR_paclitaxel.tsv")  ## adjust if needed
enrichment_df <- read_tsv(pathfindR_file, show_col_types = FALSE)

path_genes_long <- enrichment_df %>%
  dplyr::select(Term_Description, Up_regulated, Down_regulated) %>%
  dplyr::rename(pathway = Term_Description) %>%
  mutate(
    pathway        = stringr::str_trim(as.character(pathway)),
    Up_regulated   = str_split(Up_regulated,   pattern = "[,;]"),
    Down_regulated = str_split(Down_regulated, pattern = "[,;]")
  ) %>%
  tidyr::pivot_longer(
    cols      = c("Up_regulated", "Down_regulated"),
    names_to  = "direction_raw",
    values_to = "genes"
  ) %>%
  tidyr::unnest(genes) %>%
  mutate(
    Gene_Symbol = str_trim(genes)
  ) %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
  dplyr::select(pathway, Gene_Symbol) %>%
  distinct()

cat("   Number of unique pathways in pathfindR:", length(unique(path_genes_long$pathway)), "\n\n")

## ------------------------------------------------------------
## 5) Pathway-level Cox results (E_C_filter) and common pathways
## ------------------------------------------------------------

# Read pathway-level Cox results and identify pathways shared across cohorts
cat(">> Reading pathway-level Cox results (E_C_filter)...\n")

metabric_cox_file <- file.path(
  metabric_results_dir,
  "METABRIC_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx"
)
tcga_cox_file <- file.path(
  tcga_results_dir,
  "TCGA_BRCA_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx"
)

metabric_cox <- read_excel(
  metabric_cox_file,
  sheet = "E_C_filter"
) %>%
  clean_names() %>%
  mutate(pathway = stringr::str_trim(as.character(pathway)))

tcga_cox <- read_excel(
  tcga_cox_file,
  sheet = "E_C_filter"
) %>%
  clean_names() %>%
  mutate(pathway = stringr::str_trim(as.character(pathway)))

# Validate required columns for downstream HR/CI annotation in KM panels
required_cols <- c("pathway", "hr_cna", "ci_lo_cna", "ci_hi_cna")

if (!all(required_cols %in% colnames(metabric_cox))) {
  stop(
    "METABRIC E_C_filter sheet is missing required columns: ",
    paste(setdiff(required_cols, colnames(metabric_cox)), collapse = ", ")
  )
}
if (!all(required_cols %in% colnames(tcga_cox))) {
  stop(
    "TCGA E_C_filter sheet is missing required columns: ",
    paste(setdiff(required_cols, colnames(tcga_cox)), collapse = ", ")
  )
}

# Keep only the HR/CI columns needed for plotting and remove duplicates
metabric_cox_small <- metabric_cox %>%
  dplyr::select(all_of(required_cols)) %>%
  distinct()

tcga_cox_small <- tcga_cox %>%
  dplyr::select(all_of(required_cols)) %>%
  distinct()

# Identify pathways present in both METABRIC and TCGA results
common_pathways <- intersect(
  metabric_cox_small$pathway,
  tcga_cox_small$pathway
)

cat("   Number of common pathways:", length(common_pathways), "\n\n")
if (length(common_pathways) == 0) {
  stop("No common pathways found between METABRIC and TCGA-BRCA E_C_filter sheets.")
}

# Restrict pathway–gene membership table to pathways shared across cohorts
path_genes_common <- path_genes_long %>%
  filter(pathway %in% common_pathways)

## ------------------------------------------------------------
## 6) Compute sample-level CNA pathway scores (METABRIC & TCGA)
## ------------------------------------------------------------

# Compute per-sample CNA pathway scores in METABRIC using aligned CNA directions
cat(">> Computing CNA-based pathway scores for METABRIC...\n")

cna_path_metabric <- cna_long_metabric_all %>%
  inner_join(path_genes_common, by = "Gene_Symbol") %>%
  inner_join(fisher_by_gene,    by = "Gene_Symbol") %>%
  mutate(
    dir_cna = dplyr::case_when(
      !is.na(OR) & OR > 1 ~  1,
      !is.na(OR) & OR < 1 ~ -1,
      TRUE                ~ NA_real_
    ),
    cna_aligned = dplyr::case_when(
      is.na(dir_cna) ~ NA_real_,
      dir_cna ==  1  ~ cna,
      dir_cna == -1  ~ 1 - cna
    )
  )

path_cna_score_metabric <- cna_path_metabric %>%
  group_by(sample, pathway) %>%
  summarise(
    cna_score = mean(cna_aligned, na.rm = TRUE),
    .groups   = "drop"
  )

cat("   METABRIC sample–pathway CNA scores:", nrow(path_cna_score_metabric), "rows\n\n")

# Compute per-sample CNA pathway scores in TCGA-BRCA using aligned CNA directions
cat(">> Computing CNA-based pathway scores for TCGA-BRCA...\n")

cna_path_tcga <- cna_long_tcga_all %>%
  inner_join(path_genes_common, by = "Gene_Symbol") %>%
  inner_join(fisher_by_gene,    by = "Gene_Symbol") %>%
  mutate(
    dir_cna = dplyr::case_when(
      !is.na(OR) & OR > 1 ~  1,
      !is.na(OR) & OR < 1 ~ -1,
      TRUE                ~ NA_real_
    ),
    cna_aligned = dplyr::case_when(
      is.na(dir_cna) ~ NA_real_,
      dir_cna ==  1  ~ cna,
      dir_cna == -1  ~ 1 - cna
    )
  )

path_cna_score_tcga <- cna_path_tcga %>%
  group_by(sample, pathway) %>%
  summarise(
    cna_score = mean(cna_aligned, na.rm = TRUE),
    .groups   = "drop"
  )

cat("   TCGA-BRCA sample–pathway CNA scores:", nrow(path_cna_score_tcga), "rows\n\n")

# Join pathway scores with survival tables and keep valid OS records only
metabric_scores_surv <- path_cna_score_metabric %>%
  left_join(metabric_surv_small, by = "sample") %>%
  filter(!is.na(OS_TIME), !is.na(OS_EVENT))

tcga_scores_surv <- path_cna_score_tcga %>%
  left_join(tcga_surv_small, by = "sample") %>%
  filter(!is.na(OS_TIME), !is.na(OS_EVENT))

cat("   METABRIC CNA scores with OS:", nrow(metabric_scores_surv), "rows\n")
cat("   TCGA    CNA scores with OS:", nrow(tcga_scores_surv), "rows\n\n")

## ------------------------------------------------------------
## 7) Helper: safe filenames and KM panel builder
## ------------------------------------------------------------

# Convert pathway names into filesystem-safe strings for output naming
make_safe_filename <- function(pathway) {
  out <- gsub("[^A-Za-z0-9]+", "_", pathway)
  out <- gsub("_+", "_", out)
  out <- gsub("^_|_$", "", out)
  if (nchar(out) == 0) out <- "pathway"
  out
}

# Build a Kaplan–Meier panel with a risk table and HR/CI subtitle for Figure 8
make_km_panel_fig8 <- function(dat,
                               cohort_name,
                               cut_val,
                               tmax,
                               pathway_name,
                               hr_cna,
                               ci_lo_cna,
                               ci_hi_cna) {
  
  dat <- dat %>%
    mutate(
      group = ifelse(cna_score > cut_val, "High", "Low"),
      group = factor(group, levels = c("Low", "High"))
    )
  
  ## Log-rank p-value
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ group, data = dat)
  p_lr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  
  fit <- survfit(Surv(OS_TIME, OS_EVENT) ~ group, data = dat)
  
  g <- survminer::ggsurvplot(
    fit,
    data               = dat,
    risk.table         = TRUE,
    risk.table.fontsize = 5,
    risk.table.height  = 0.28,
    pval               = FALSE,
    conf.int           = FALSE,
    legend.title       = "CNA burden ",
    legend.labs        = c("Low CNA burden", "High CNA burden"),
    xlab               = "Time (months)",
    ylab               = "OS probability",
    xlim               = c(0, tmax),
    ggtheme            = theme_bw(),
    palette            = c("#E41A1C", "#377EB8")
  )
  
  ## Risk table y-axis: Low / High
  strata_levels <- levels(g$table$data$strata)
  if (length(strata_levels) == 2) {
    g$table <- g$table +
      scale_y_discrete(
        breaks = strata_levels,
        labels = c("Low", "High")
      )
  }
  
  
  ann_label <- sprintf("median CNA score = %.3f\nlog-rank p = %.3g",
                       cut_val, p_lr)
  
  ## KM plot
  p_km <- g$plot +
    annotate(
      "text",
      x        = tmax * 0.55,
      y        = 0.95,
      label    = ann_label,
      hjust    = 0,
      vjust    = 1,
      size     = 4.5,
      fontface = "bold"
    )  +
    scale_x_continuous(
      limits = c(0, tmax),
      breaks = seq(0, tmax, by = 100),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title    = paste0(pathway_name, " (", cohort_name, ", OS)"),
      subtitle = sprintf("HR = %.2f, 95%% CI: [%.2f - %.2f]",
                         hr_cna, ci_lo_cna, ci_hi_cna)
    ) +
    theme(
      plot.title    = element_text(face = "bold", hjust = .5, size = 20),
      plot.subtitle = element_text(hjust = .5, size = 16),
      axis.title.x  = element_text(face = "bold", size = 16),
      axis.title.y  = element_text(face = "bold", size = 16),
      axis.text.x   = element_text(face = "bold", size = 16),
      axis.text.y   = element_text(face = "bold", size = 16),
      legend.title  = element_text(face = "bold", size = 16),
      legend.text   = element_text(size = 16)
    )
  
  ## Risk table
  p_risk <- g$table +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x  = element_text(face = "bold", size = 16),
      axis.text.y  = element_text(face = "bold", size = 16),
      axis.title.x = element_text(size = 16),
      panel.grid   = element_blank()
    )
  
  ## Combine: KM on top, risk table at the bottom
  p_panel <- p_km / p_risk +
    patchwork::plot_layout(heights = c(3, 1))
  
  p_panel
}

## ------------------------------------------------------------
## 8) KM plots for 3 selected pathways → Figure 8.jpeg
##     a) Ras signaling pathway
##     b) Axon guidance
##     c) Focal adhesion
## ------------------------------------------------------------

# Generate Figure 8 by plotting KM curves for three predefined pathways in both cohorts
cat(">> Building Figure 8 (KM curves for Ras / Axon guidance / Focal adhesion)...\n")

## Axis limits for METABRIC and TCGA
tmax_met  <- 400
tmax_tcga <- 300

## Selected pathways
path_vec <- c("Ras signaling pathway",
              "Axon guidance",
              "Focal adhesion")

row_panels <- list()

## Spacer for horizontal gap between METABRIC and TCGA panels
spacer_h <- ggplot2::ggplot() + ggplot2::theme_void()

for (pw in path_vec) {
  cat("   - Pathway:", pw, "\n")
  
  ## Extract HR and CI for METABRIC and TCGA from Cox tables
  met_row <- metabric_cox_small %>%
    filter(pathway == pw)
  tcga_row <- tcga_cox_small %>%
    filter(pathway == pw)
  
  if (nrow(met_row) == 0 || nrow(tcga_row) == 0) {
    next
  }
  
  ## Take the first entry if duplicated
  met_row  <- met_row[1, , drop = FALSE]
  tcga_row <- tcga_row[1, , drop = FALSE]
  
  hr_met  <- as.numeric(met_row$hr_cna)
  lo_met  <- as.numeric(met_row$ci_lo_cna)
  hi_met  <- as.numeric(met_row$ci_hi_cna)
  
  hr_tcga <- as.numeric(tcga_row$hr_cna)
  lo_tcga <- as.numeric(tcga_row$ci_lo_cna)
  hi_tcga <- as.numeric(tcga_row$ci_hi_cna)
  
  ## CNA pathway score + OS for the selected pathway
  met_dat <- metabric_scores_surv %>%
    filter(pathway == pw) %>%
    dplyr::select(sample, pathway, cna_score, OS_TIME, OS_EVENT)
  
  tcga_dat <- tcga_scores_surv %>%
    filter(pathway == pw) %>%
    dplyr::select(sample, pathway, cna_score, OS_TIME, OS_EVENT)
  
  if (nrow(met_dat) < 10 || nrow(tcga_dat) < 10) {
    next
  }
  
  ## METABRIC median – the same cut value is reused for TCGA
  cut_val <- median(met_dat$cna_score, na.rm = TRUE)
  if (is.na(cut_val)) {
    next
  }
  
  ## KM + risk table panels (with HR & CI subtitle)
  p_met_panel  <- make_km_panel_fig8(
    met_dat,  "METABRIC",  cut_val, tmax_met,
    pw, hr_met, lo_met, hi_met
  )
  p_tcga_panel <- make_km_panel_fig8(
    tcga_dat, "TCGA-BRCA", cut_val, tmax_tcga,
    pw, hr_tcga, lo_tcga, hi_tcga
  )
  
  ## Horizontal row: METABRIC | spacer | TCGA
  row_panels[[pw]] <- patchwork::wrap_plots(
    p_met_panel,
    spacer_h,
    p_tcga_panel,
    ncol   = 3,
    widths = c(1, 0.18, 1),
    guides = "collect"
  ) &
    theme(
      legend.position = "bottom",
      legend.title    = element_text(face = "bold", size = 18),
      legend.text     = element_text(size = 18)
    )
}

## Order rows: a) Ras, b) Axon, c) Focal
row_ras   <- row_panels[["Ras signaling pathway"]]
row_axon  <- row_panels[["Axon guidance"]]
row_focal <- row_panels[["Focal adhesion"]]

## Hide legends in the first two rows so that only the last row keeps a legend
if (!is.null(row_ras)) {
  row_ras  <- row_ras  & theme(legend.position = "none")
}
if (!is.null(row_axon)) {
  row_axon <- row_axon & theme(legend.position = "none")
}

## Vertical spacer between rows
spacer <- ggplot2::ggplot() + ggplot2::theme_void()

## Stack three rows using cowplot, with labels a., b., c.
fig8 <- cowplot::plot_grid(
  row_ras,
  spacer,                     # gap between row 1 and row 2
  row_axon,
  spacer,                     # gap between row 2 and row 3
  row_focal,
  ncol           = 1,
  rel_heights    = c(1, 0.08, 1, 0.08, 1),
  labels         = c("A.", "", "B.", "", "C."),
  label_size     = 24,
  label_fontface = "bold",
  hjust          = -0.1
)

## Save as Figure 8
ggsave(
  filename = file.path(figures_dir, "Figure 8.jpeg"),
  plot     = fig8 + theme(plot.margin = margin(t = 10, r = 20, b = 10, l = 20)),
  width    = 14,
  height   = 20,
  dpi      = 300
)

cat(">> Figure 8 saved to:\n   ", file.path(figures_dir, "Figure 8.jpeg"), "\n")