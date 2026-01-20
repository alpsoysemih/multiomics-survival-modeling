# Load required packages quietly
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(org.Hs.eg.db))

# Disable automatic string-to-factor conversion
options(stringsAsFactors = FALSE)

## ------------------------------------------------------------
## Common KM theme (slightly larger axis labels/ticks)
## ------------------------------------------------------------
# Define a consistent theme for Kaplan–Meier plots
km_theme <- theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title.x  = element_text(size = 12),
    axis.title.y  = element_text(size = 12),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 12),
    legend.title  = element_text(face = "bold", size = 13),  
    legend.text   = element_text(size = 13)
  )

## ------------------------------------------------------------
## 0) Directories
## ------------------------------------------------------------

# Define input and output directories for datasets and results
metabric_dir         <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir             <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
results_root         <- "/Volumes/Expansion/Prognostic_Analysis/Results"

# Define cohort-specific results directories
metabric_results_dir <- file.path(results_root, "METABRIC")
tcga_results_dir     <- file.path(results_root, "TCGA-BRCA")

# Create the manuscript figures directory if it does not exist
figures_dir <- file.path(results_root, "Figures_Manuscript")
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}

# Print the figure output path
cat(">> Figures will be saved in:\n   ", figures_dir, "\n\n")

## ------------------------------------------------------------
## 1) METABRIC: expression, CNA, OS survival
## ------------------------------------------------------------

# Define METABRIC file paths for expression, CNA, and clinical survival data
metabric_expr_file     <- file.path(metabric_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
metabric_cna_file      <- file.path(metabric_dir, "data_cna.txt")
metabric_survival_file <- file.path(metabric_dir, "brca_metabric_clinical_data.tsv")

# Read METABRIC expression, CNA, and clinical survival tables
cat(">> Reading METABRIC data...\n")

metabric_expr     <- read_tsv(metabric_expr_file,     show_col_types = FALSE)
metabric_cna      <- read_tsv(metabric_cna_file,      show_col_types = FALSE)
metabric_survival <- read_tsv(metabric_survival_file, show_col_types = FALSE)

## Drop Entrez ID, keep gene symbol
# Keep gene symbols and remove Entrez identifiers from expression matrix
metabric_expr <- metabric_expr %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

# Keep gene symbols and remove Entrez identifiers from CNA matrix
metabric_cna <- metabric_cna %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

## Binary CNA (0 vs != 0)
# Convert CNA values into a binary indicator (0 vs non-zero)
metabric_cna_binary <- metabric_cna
metabric_cna_binary[, -1] <- ifelse(metabric_cna_binary[, -1] != 0, 1, 0)

## METABRIC OS (months) + covariates
# Build a cleaned METABRIC survival table with OS and clinical covariates
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

# Report the number of METABRIC samples with OS data
cat("   Number of METABRIC samples with OS data: ",
    nrow(metabric_surv_small), "\n\n")

## Long format expression / CNA
# Convert METABRIC expression matrix to long format (one row per gene-sample)
expr_long_metabric <- metabric_expr %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

# Convert METABRIC binary CNA matrix to long format (one row per gene-sample)
cna_long_metabric <- metabric_cna_binary %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

## Helper: per-gene METABRIC data
# Merge expression, CNA, and survival covariates for a single gene (METABRIC)
get_gene_surv_data_met <- function(gene_name) {
  expr_long_metabric %>%
    dplyr::filter(Gene_Symbol == gene_name) %>%
    dplyr::select(sample, expr) %>%
    dplyr::inner_join(
      cna_long_metabric %>%
        dplyr::filter(Gene_Symbol == gene_name) %>%
        dplyr::select(sample, cna),
      by = "sample"
    ) %>%
    dplyr::inner_join(metabric_surv_small, by = "sample") %>%
    dplyr::mutate(
      expr = as.numeric(expr),
      cna  = as.numeric(cna)
    ) %>%
    dplyr::filter(
      !is.na(OS_TIME),
      !is.na(OS_EVENT)
    ) %>%
    dplyr::distinct()
}

## ------------------------------------------------------------
## 2) TCGA-BRCA: expression, CNA, OS survival
## ------------------------------------------------------------

# Define TCGA-BRCA file paths for expression, CNA, and clinical survival data
tcga_expr_file     <- file.path(tcga_dir, "TCGA-BRCA_exprs.z.tsv")
tcga_cna_file      <- file.path(tcga_dir, "TCGA-BRCA_CNA.tsv")
tcga_survival_file <- file.path(tcga_dir, "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv")

# Read TCGA-BRCA expression, CNA, and clinical survival tables
cat(">> Reading TCGA-BRCA data...\n")

tcga_expr     <- read_tsv(tcga_expr_file,     show_col_types = FALSE)
tcga_cna      <- read_tsv(tcga_cna_file,      show_col_types = FALSE)
tcga_survival <- read_tsv(tcga_survival_file, show_col_types = FALSE)

## Binary CNA (0 vs != 0)
# Convert CNA values into a binary indicator (0 vs non-zero)
tcga_cna_binary <- tcga_cna
tcga_cna_binary[, -1] <- ifelse(tcga_cna_binary[, -1] != 0, 1, 0)

## Map Entrez IDs → gene symbols
# Map TCGA expression Entrez IDs to gene symbols
tcga_expr$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_expr$ENTREZID),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Map TCGA CNA Entrez IDs to gene symbols
tcga_cna_binary$gene <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_cna_binary$gene),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Rename mapped gene symbol columns to a shared name
colnames(tcga_expr)[colnames(tcga_expr) == "ENTREZID"]         <- "Gene_Symbol"
colnames(tcga_cna_binary)[colnames(tcga_cna_binary) == "gene"] <- "Gene_Symbol"

## TCGA OS (days) + covariates → months
# Build a cleaned TCGA survival table and convert OS time from days to months
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

# Report the number of TCGA-BRCA samples with OS data
cat("   Number of TCGA-BRCA samples with OS data: ",
    nrow(tcga_surv_small), "\n\n")

## Long format expression / CNA
# Convert TCGA expression matrix to long format and harmonize sample barcodes
expr_long_tcga <- tcga_expr %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "expr"
  ) %>%
  dplyr::mutate(
    sample = substr(sample_barcode, 1, 15)
  )

# Convert TCGA binary CNA matrix to long format and harmonize sample barcodes
cna_long_tcga <- tcga_cna_binary %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  dplyr::mutate(
    sample = substr(sample_barcode, 1, 15)
  )

## Helper: per-gene TCGA data
# Merge expression, CNA, and survival covariates for a single gene (TCGA-BRCA)
get_gene_surv_data_tcga <- function(gene_name) {
  expr_long_tcga %>%
    dplyr::filter(Gene_Symbol == gene_name) %>%
    dplyr::select(sample, expr) %>%
    dplyr::inner_join(
      cna_long_tcga %>%
        dplyr::filter(Gene_Symbol == gene_name) %>%
        dplyr::select(sample, cna),
      by = "sample"
    ) %>%
    dplyr::inner_join(tcga_surv_small, by = "sample") %>%
    dplyr::mutate(
      expr = as.numeric(expr),
      cna  = as.numeric(cna)
    ) %>%
    dplyr::filter(
      !is.na(OS_TIME),
      !is.na(OS_EVENT)
    ) %>%
    dplyr::distinct()
}

## ------------------------------------------------------------
## 3) Read Cox Excel (E_C_Age_Stg_PAM50) and extract HRs
## ------------------------------------------------------------

# Define file paths to Cox regression result workbooks
cox_met_file  <- file.path(
  metabric_results_dir,
  "METABRIC_OS_all_Cox_all_models_by_gene_with_LRT.xlsx"
)
cox_tcga_file <- file.path(
  tcga_results_dir,
  "TCGA_BRCA_OS_all_Cox_all_models_by_gene_with_LRT.xlsx"
)

# Specify the sheet that contains the multivariable Cox model results
sheet_name <- "E_C_Age_Stg_PAM50"

# Read Cox regression results for METABRIC and TCGA-BRCA
cat(">> Reading Cox results (sheet: ", sheet_name, ")...\n", sep = "")

cox_met  <- readxl::read_excel(cox_met_file,  sheet = sheet_name)
cox_tcga <- readxl::read_excel(cox_tcga_file, sheet = sheet_name)

## Helper to extract HR / CI / p for a given gene and effect type
# Extract HR, 95% CI, and p-value for expression or CNA for a single gene
get_cox_stats <- function(df, gene_name, effect = c("expr", "cna")) {
  effect <- match.arg(effect)
  
  if (!("gene" %in% names(df))) {
    stop("Cox table does not have a 'gene' column.")
  }
  
  row <- df %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::slice(1)
  
  if (nrow(row) == 0L) {
    stop("Gene '", gene_name, "' not found in Cox table.")
  }
  
  if (effect == "expr") {
    if (!all(c("HR_expr", "CI_lo_expr", "CI_hi_expr", "p_expr") %in% names(row))) {
      stop("Expression columns missing in Cox table for gene ", gene_name)
    }
    return(list(
      HR   = as.numeric(row$HR_expr),
      LCL  = as.numeric(row$CI_lo_expr),
      UCL  = as.numeric(row$CI_hi_expr),
      PVAL = as.numeric(row$p_expr)
    ))
  } else {
    if (!all(c("HR_cna", "CI_lo_cna", "CI_hi_cna", "p_cna") %in% names(row))) {
      stop("CNA columns missing in Cox table for gene ", gene_name)
    }
    return(list(
      HR   = as.numeric(row$HR_cna),
      LCL  = as.numeric(row$CI_lo_cna),
      UCL  = as.numeric(row$CI_hi_cna),
      PVAL = as.numeric(row$p_cna)
    ))
  }
}

## ------------------------------------------------------------
## 4) KM helper functions
## ------------------------------------------------------------

## SERPINE1 expression – METABRIC (optimal cutoff)
# Create a METABRIC KM plot for expression using an optimal surv_cutpoint threshold
make_km_expr_opt_met <- function(gene_name, HR, LCL, UCL, PVAL, tag = NULL) {  
  d <- get_gene_surv_data_met(gene_name) %>%
    dplyr::filter(!is.na(expr))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    stop("Not enough data for METABRIC expression KM: ", gene_name)
  }
  
  cut_obj <- survminer::surv_cutpoint(
    d %>% dplyr::select(OS_TIME, OS_EVENT, expr),
    time      = "OS_TIME",
    event     = "OS_EVENT",
    variables = "expr",
    minprop   = 0.1
  )
  cut_val <- cut_obj$cutpoint$cutpoint[1]
  
  d <- d %>%
    dplyr::mutate(
      expr_group = factor(
        ifelse(expr > cut_val, "High", "Low"),
        levels = c("Low", "High")
      )
    )
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  subtitle_txt <- sprintf("HR = %.2f, 95%% CI: [%.2f – %.2f]", HR, LCL, UCL)
  pval_txt     <- paste0("log-rank p = ", signif(lr_p, 3))
  
  g <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    pval                 = pval_txt,
    legend.labs          = levels(d$expr_group),
    title                = paste0(gene_name, " Expression (METABRIC, OS)"),
    legend.title         = "Expression ",
    xlab                 = "Time (months)",
    ylab                 = "OS probability",
    ggtheme              = theme_bw(),
    risk.table.fontsize  = 3.7
  )
  
  plot_part  <- g$plot  +
    labs(subtitle = subtitle_txt, tag = tag) +
    km_theme +
    theme(
      plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
      axis.title.x = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold", size = 13)
    )
  
  table_part <- g$table + theme_bw() +
    theme(
      axis.title.y = element_blank(),   # remove "Strata"
      axis.title.x = element_blank(),   # remove x axis label
      axis.text.x  = element_blank(),   # remove x tick labels
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12, face = "bold")
    )
  
  plot_part / table_part + patchwork::plot_layout(heights = c(3, 1))
}

## SERPINE1 expression – TCGA (optimal cutoff, OS in months)
# Create a TCGA-BRCA KM plot for expression using an optimal surv_cutpoint threshold
make_km_expr_opt_tcga <- function(gene_name, HR, LCL, UCL, PVAL, tag = NULL) {
  d <- get_gene_surv_data_tcga(gene_name) %>%
    dplyr::filter(!is.na(expr))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    stop("Not enough data for TCGA-BRCA expression KM: ", gene_name)
  }
  
  cut_obj <- survminer::surv_cutpoint(
    d %>% dplyr::select(OS_TIME, OS_EVENT, expr),
    time      = "OS_TIME",
    event     = "OS_EVENT",
    variables = "expr",
    minprop   = 0.1
  )
  cut_val <- cut_obj$cutpoint$cutpoint[1]
  
  d <- d %>%
    dplyr::mutate(
      expr_group = factor(
        ifelse(expr > cut_val, "High", "Low"),
        levels = c("Low", "High")
      )
    )
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  subtitle_txt <- sprintf("HR = %.2f, 95%% CI: [%.2f – %.2f]", HR, LCL, UCL)
  pval_txt     <- paste0("log-rank p = ", signif(lr_p, 3))
  
  g <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    pval                 = pval_txt,
    legend.labs          = levels(d$expr_group),
    title                = paste0(gene_name, " Expression (TCGA-BRCA, OS)"),
    xlab                 = "Time (months)",
    ylab                 = "OS probability",
    ggtheme              = theme_bw(),
    legend.title         = "Expression ",
    risk.table.fontsize  = 4
  )
  
  plot_part  <- g$plot  +
    labs(subtitle = subtitle_txt, tag = tag) +
    km_theme +
    theme(
      plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
      axis.title.x = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold", size = 13)
    )
  
  table_part <- g$table + theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12, face = "bold")
    )
  
  plot_part / table_part + patchwork::plot_layout(heights = c(3, 1))
}

## CNA – METABRIC (0 vs 1), with optional legend
# Create a METABRIC KM plot for binary CNA status (0 vs 1)
make_km_cna_met <- function(gene_name, HR, LCL, UCL, PVAL,
                            show_legend = FALSE, tag = NULL) {
  d <- get_gene_surv_data_met(gene_name) %>%
    dplyr::filter(!is.na(cna))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    stop("Not enough data for METABRIC CNA KM: ", gene_name)
  }
  
  d <- d %>%
    dplyr::mutate(
      cna_group = factor(
        ifelse(cna == 0, "CNA = 0", "CNA = 1"),
        levels = c("CNA = 0", "CNA = 1")
      )
    )
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  subtitle_txt <- sprintf("HR = %.2f, 95%% CI: [%.2f – %.2f]", HR, LCL, UCL)
  pval_txt     <- paste0("log-rank p = ", signif(lr_p, 3))
  
  g <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    pval                 = pval_txt,
    legend.labs          = levels(d$cna_group),
    title                = paste0(gene_name, " CNA (METABRIC, OS)"),
    xlab                 = "Time (months)",
    ylab                 = "",          # no y-label for CNA panels
    ggtheme              = theme_bw(),
    legend.title         = "CNA ",
    risk.table.fontsize    = 4
  )
  
  plot_part <- g$plot +
    labs(subtitle = subtitle_txt, tag = tag) +
    km_theme +
    theme(
      plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
      axis.title.x = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold", size = 13)
    )
  if (!show_legend) {
    plot_part <- plot_part + theme(legend.position = "none")
  }
  
  table_part <- g$table + theme_bw() +
    theme(
      axis.title.y = element_blank(),   # remove "Strata"
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12, face = "bold")
    )
  
  plot_part / table_part + patchwork::plot_layout(heights = c(3, 1))
}

# Create a TCGA-BRCA KM plot for binary CNA status (0 vs 1)
make_km_cna_tcga <- function(gene_name, HR, LCL, UCL, PVAL,
                             show_legend = FALSE, tag = NULL) {
  d <- get_gene_surv_data_tcga(gene_name) %>%
    dplyr::filter(!is.na(cna))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    stop("Not enough data for TCGA-BRCA CNA KM: ", gene_name)
  }
  
  d <- d %>%
    dplyr::mutate(
      cna_group = factor(
        ifelse(cna == 0, "CNA = 0", "CNA = 1"),
        levels = c("CNA = 0", "CNA = 1")
      )
    )
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  subtitle_txt <- sprintf("HR = %.2f, 95%% CI: [%.2f – %.2f]", HR, LCL, UCL)
  pval_txt     <- paste0("log-rank p = ", signif(lr_p, 3))
  
  g <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    pval                 = pval_txt,
    legend.labs          = levels(d$cna_group),
    title                = paste0(gene_name, " CNA (TCGA-BRCA, OS)"),
    xlab                 = "Time (months)",
    ylab                 = "",          # no y-label for CNA panels
    ggtheme              = theme_bw(),
    legend.title         = "CNA ",
    risk.table.fontsize    = 4
  )
  
  plot_part <- g$plot +
    labs(subtitle = subtitle_txt, tag = tag) +
    km_theme +
    theme(
      plot.tag = element_text(face = "bold", size = 16, hjust = 0, vjust = 1),
      axis.title.x = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold", size = 13)
    )
  if (!show_legend) {
    plot_part <- plot_part + theme(legend.position = "none")
  }
  
  table_part <- g$table + theme_bw() +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 12, face = "bold")
    )
  
  plot_part / table_part + patchwork::plot_layout(heights = c(3, 1))
}
 

## ------------------------------------------------------------
## 5) Build KM panels for all genes  (UPDATED: add LAMB2; retag CNA panels)
## ------------------------------------------------------------

## SERPINE1 expression
cat(">> Building KM panels for SERPINE1 expression...\n")

serp_met_stats  <- get_cox_stats(cox_met,  "SERPINE1", effect = "expr")
serp_tcga_stats <- get_cox_stats(cox_tcga, "SERPINE1", effect = "expr")

p_serp_met  <- make_km_expr_opt_met(
  gene_name = "SERPINE1",
  HR        = serp_met_stats$HR,
  LCL       = serp_met_stats$LCL,
  UCL       = serp_met_stats$UCL,
  PVAL      = serp_met_stats$PVAL,
  tag       = "A."
)

p_serp_tcga <- make_km_expr_opt_tcga(
  gene_name = "SERPINE1",
  HR        = serp_tcga_stats$HR,
  LCL       = serp_tcga_stats$LCL,
  UCL       = serp_tcga_stats$UCL,
  PVAL      = serp_tcga_stats$PVAL,
  tag       = NULL
)

## CXCR4 CNA – no legend (tag B.)
cat(">> Building KM panels for CXCR4 CNA...\n")

cxcr4_met_stats  <- get_cox_stats(cox_met,  "CXCR4", effect = "cna")
cxcr4_tcga_stats <- get_cox_stats(cox_tcga, "CXCR4", effect = "cna")

p_cxcr4_met <- make_km_cna_met(
  gene_name   = "CXCR4",
  HR          = cxcr4_met_stats$HR,
  LCL         = cxcr4_met_stats$LCL,
  UCL         = cxcr4_met_stats$UCL,
  PVAL        = cxcr4_met_stats$PVAL,
  show_legend = FALSE,
  tag         = "B."
)

p_cxcr4_tcga <- make_km_cna_tcga(
  gene_name   = "CXCR4",
  HR          = cxcr4_tcga_stats$HR,
  LCL         = cxcr4_tcga_stats$LCL,
  UCL         = cxcr4_tcga_stats$UCL,
  PVAL        = cxcr4_tcga_stats$PVAL,
  show_legend = FALSE,
  tag         = NULL
)

## LAMB2 CNA – no legend (NEW, tag C.)
cat(">> Building KM panels for LAMB2 CNA...\n")

lamb2_met_stats  <- get_cox_stats(cox_met,  "LAMB2", effect = "cna")
lamb2_tcga_stats <- get_cox_stats(cox_tcga, "LAMB2", effect = "cna")

p_lamb2_met <- make_km_cna_met(
  gene_name   = "LAMB2",
  HR          = lamb2_met_stats$HR,
  LCL         = lamb2_met_stats$LCL,
  UCL         = lamb2_met_stats$UCL,
  PVAL        = lamb2_met_stats$PVAL,
  show_legend = FALSE,
  tag         = "C."
)

p_lamb2_tcga <- make_km_cna_tcga(
  gene_name   = "LAMB2",
  HR          = lamb2_tcga_stats$HR,
  LCL         = lamb2_tcga_stats$LCL,
  UCL         = lamb2_tcga_stats$UCL,
  PVAL        = lamb2_tcga_stats$PVAL,
  show_legend = FALSE,
  tag         = NULL
)

## CLDN1 CNA – no legend (MOVE: tag D.)
cat(">> Building KM panels for CLDN1 CNA...\n")

cldn1_met_stats  <- get_cox_stats(cox_met,  "CLDN1", effect = "cna")
cldn1_tcga_stats <- get_cox_stats(cox_tcga, "CLDN1", effect = "cna")

p_cldn1_met <- make_km_cna_met(
  gene_name   = "CLDN1",
  HR          = cldn1_met_stats$HR,
  LCL         = cldn1_met_stats$LCL,
  UCL         = cldn1_met_stats$UCL,
  PVAL        = cldn1_met_stats$PVAL,
  show_legend = FALSE,
  tag         = "D."
)

p_cldn1_tcga <- make_km_cna_tcga(
  gene_name   = "CLDN1",
  HR          = cldn1_tcga_stats$HR,
  LCL         = cldn1_tcga_stats$LCL,
  UCL         = cldn1_tcga_stats$UCL,
  PVAL        = cldn1_tcga_stats$PVAL,
  show_legend = FALSE,
  tag         = NULL
)

## LYPD6B CNA – legend ONLY here (MOVE: tag E.)
cat(">> Building KM panels for LYPD6B CNA...\n")

lypd6b_met_stats  <- get_cox_stats(cox_met,  "LYPD6B", effect = "cna")
lypd6b_tcga_stats <- get_cox_stats(cox_tcga, "LYPD6B", effect = "cna")

p_lypd6b_met <- make_km_cna_met(
  gene_name   = "LYPD6B",
  HR          = lypd6b_met_stats$HR,
  LCL         = lypd6b_met_stats$LCL,
  UCL         = lypd6b_met_stats$UCL,
  PVAL        = lypd6b_met_stats$PVAL,
  show_legend = TRUE,
  tag         = "E."
)

p_lypd6b_tcga <- make_km_cna_tcga(
  gene_name   = "LYPD6B",
  HR          = lypd6b_tcga_stats$HR,
  LCL         = lypd6b_tcga_stats$LCL,
  UCL         = lypd6b_tcga_stats$UCL,
  PVAL        = lypd6b_tcga_stats$PVAL,
  show_legend = FALSE,
  tag         = NULL
)

## ------------------------------------------------------------
## 6) 2 x 5 grid: top row METABRIC, bottom row TCGA-BRCA
##    A column (SERPINE1) has its own legend at the bottom
##    B–E (CNA columns) share a common legend at the bottom
## ------------------------------------------------------------

row_spacer <- patchwork::plot_spacer()

col_a <- patchwork::wrap_plots(
  p_serp_met,
  row_spacer,
  p_serp_tcga,
  nrow    = 3,
  ncol    = 1,
  heights = c(1, 0.05, 1),
  guides  = "collect"
) & theme(
  legend.position = "bottom",
  legend.text     = element_text( size = 13),
  legend.key.size = unit(0.8, "cm"),
  legend.spacing  = unit(0.4, "cm")
)

## B–E CNA grid: CXCR4 (B), LAMB2 (C), CLDN1 (D), LYPD6B (E)
col_bcde <- patchwork::wrap_plots(
  # top row: METABRIC
  p_cxcr4_met, p_lamb2_met, p_cldn1_met, p_lypd6b_met,
  # middle row: spacer
  row_spacer,  row_spacer,  row_spacer,  row_spacer,
  # bottom row: TCGA-BRCA
  p_cxcr4_tcga, p_lamb2_tcga, p_cldn1_tcga, p_lypd6b_tcga,
  nrow    = 3,
  ncol    = 4,
  heights = c(1, 0.05, 1),
  guides  = "collect"
) & theme(
  legend.position = "bottom",
  legend.text     = element_text(size = 13),
  legend.key.size = unit(0.8, "cm"),
  legend.spacing  = unit(0.4, "cm")
)

combined_km <- (col_a | col_bcde) +
  patchwork::plot_layout(widths = c(1, 4))

fig_file_jpeg <- file.path(figures_dir, "Figure 3.jpeg")

ggplot2::ggsave(
  filename = fig_file_jpeg,
  plot     = combined_km,
  width    = 26,   # widened for 5 columns
  height   = 14,
  dpi      = 300
)

cat(">> 2 x 5 KM figure with risk tables saved to:\n")
cat("   ", fig_file_jpeg, "\n\n")