# Load required libraries quietly to keep the console output clean
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ggplotify))


## --------------------------------------------------
## 0) Helper functions
## --------------------------------------------------

# Sanitize strings for safe filenames by replacing non-alphanumerics with underscores
clean_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^[:alnum:]]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# Convert TCGA survival time from days to months
tcga_to_months <- function(x) x / 30.4375

# Convert a model sheet name into a readable formula-like label
pretty_model_label <- function(sheet_name) {
  base <- gsub("_filter2$", "", sheet_name)
  base <- gsub("_filter$",  "", base)
  tokens <- unlist(strsplit(base, "_"))
  map <- c(
    "E"     = "Expr",
    "C"     = "CNA",
    "EC"    = "Expr×CNA",
    "Age"   = "Age",
    "Stg"   = "Stage",
    "PAM50" = "PAM50"
  )
  pretty <- ifelse(tokens %in% names(map), map[tokens], tokens)
  paste(pretty, collapse = " + ")
}

## --------------------------------------------------
## 1) Directories and output
## --------------------------------------------------

# Define input dataset directories and output result locations
metabric_dir         <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir             <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
pathfindR_dir        <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"
gdsc_dir             <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
metabric_results_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/METABRIC"
tcga_results_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Results/TCGA-BRCA"

# Create the manuscript figures directory if it does not exist
figures_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

## --------------------------------------------------
## 2) METABRIC expression, CNA, survival
## --------------------------------------------------

# Set file paths for METABRIC expression, CNA, and clinical survival data
metabric_expr_file     <- file.path(metabric_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
metabric_cna_file      <- file.path(metabric_dir, "data_cna.txt")
metabric_survival_file <- file.path(metabric_dir, "brca_metabric_clinical_data.tsv")

# Read METABRIC expression, CNA, and survival tables
metabric_expr     <- read_tsv(metabric_expr_file,     show_col_types = FALSE)
metabric_cna      <- read_tsv(metabric_cna_file,      show_col_types = FALSE)
metabric_survival <- read_tsv(metabric_survival_file, show_col_types = FALSE)

## Drop Entrez ID, keep gene symbol
metabric_expr <- metabric_expr %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

metabric_cna <- metabric_cna %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol)

## Make CNA binary (0 vs !=0)
metabric_cna_binary <- metabric_cna
metabric_cna_binary[, -1] <- ifelse(metabric_cna_binary[, -1] != 0, 1, 0)

## METABRIC OS + covariates (months)
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

# Report how many METABRIC samples have valid OS information after filtering
message("Number of METABRIC samples with OS data: ", nrow(metabric_surv_small))

## --------------------------------------------------
## 3) METABRIC expression & CNA long format
## --------------------------------------------------

# Convert METABRIC expression matrix into long format (sample-level rows per gene)
expr_long_metabric <- metabric_expr %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

# Convert METABRIC binary CNA matrix into long format (sample-level rows per gene)
cna_long_metabric <- metabric_cna_binary %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

## --------------------------------------------------
## 4) GDSC CNA vs response: Fisher OR per gene
## --------------------------------------------------

# Set file paths for GDSC expression, CNA, and Paclitaxel response data
gdsc_expr_file     <- file.path(gdsc_dir, "GDSC_exprs.Paclitaxel.eb_with.TCGA_exprs.Paclitaxel.tsv")
gdsc_cna_file      <- file.path(gdsc_dir, "GDSC_CNA.Paclitaxel.tsv")
gdsc_response_file <- file.path(gdsc_dir, "GDSC_response.Paclitaxel.tsv")

# Read GDSC expression, CNA, and response tables
gdsc_expr     <- read_tsv(gdsc_expr_file,     show_col_types = FALSE)
gdsc_cna      <- read_tsv(gdsc_cna_file,      show_col_types = FALSE)
gdsc_response <- read_tsv(gdsc_response_file, show_col_types = FALSE)

# Binarize GDSC CNA values into 0 vs 1
gdsc_cna[, -1] <- ifelse(gdsc_cna[, -1] != 0, 1, 0)

# Map Entrez IDs to gene symbols for the GDSC expression table
gdsc_expr$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(gdsc_expr$ENTREZID),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Map Entrez IDs to gene symbols for the GDSC CNA table
gdsc_cna$gene_id <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(gdsc_cna$gene_id),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# Standardize gene symbol column names across expression and CNA tables
colnames(gdsc_expr)[colnames(gdsc_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(gdsc_cna)[colnames(gdsc_cna)   == "gene_id"]  <- "Gene_Symbol"

# Create a response table with a consistent status column
resp <- gdsc_response %>%
  dplyr::mutate(status = response)

# Identify samples present in the response table for subsetting matrices
gdsc_samples <- as.character(unique(resp$sample_name))

# Subset GDSC expression matrix to the samples of interest
expr_gdsc <- gdsc_expr %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

# Z-score normalize expression per gene (row-wise) and replace NAs with zero
expr_mat <- as.matrix(expr_gdsc[, -1])
rownames(expr_mat) <- expr_gdsc$Gene_Symbol
expr_z <- t(scale(t(expr_mat)))
expr_z[is.na(expr_z)] <- 0
expr_z_df <- as_tibble(expr_z, rownames = "Gene_Symbol")

# Subset GDSC CNA matrix to the samples of interest
cna_gdsc <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

# Build a binary CNA matrix and convert to tibble with gene symbols
cna_mat <- as.matrix(cna_gdsc[, -1])
rownames(cna_mat) <- cna_gdsc$Gene_Symbol
cna_bin <- ifelse(cna_mat != 0, 1, 0)
cna_bin_df <- as_tibble(cna_bin, rownames = "Gene_Symbol")

# Convert GDSC expression Z-scores to long format and join with response status
expr_long_gdsc <- expr_z_df %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr_z"
  ) %>%
  dplyr::rename(sample_name = "sample") %>%
  dplyr::inner_join(
    resp %>%
      dplyr::select(sample_name, status) %>%
      dplyr::mutate(sample_name = as.character(sample_name)),
    by = "sample_name"
  )

# Convert GDSC CNA binary calls to long format and join with response status
cna_long_gdsc <- cna_bin_df %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna_bin"
  ) %>%
  dplyr::rename(sample_name = "sample") %>%
  dplyr::inner_join(
    resp %>%
      dplyr::select(sample_name, status) %>%
      dplyr::mutate(sample_name = as.character(sample_name)),
    by = "sample_name"
  )

# Compute Fisher's exact test odds ratio and p-value per gene (CNA vs response)
fisher_by_gene <- cna_long_gdsc %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::group_modify(~ {
    tab <- table(.x$cna_bin, .x$status)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      tibble(OR = NA_real_, p_value = NA_real_)
    } else {
      ft <- fisher.test(tab)
      tibble(OR = unname(ft$estimate), p_value = ft$p.value)
    }
  }) %>%
  dplyr::ungroup()

## --------------------------------------------------
## 5) PathfindR enrichment and pathway–gene mapping
## --------------------------------------------------

# Read pathfindR enrichment results for Paclitaxel and expand pathway-to-gene mappings
enrichment_file <- file.path(pathfindR_dir, "pathfindR_paclitaxel.tsv")
enrichment_df   <- read_tsv(enrichment_file, show_col_types = FALSE)

# Build a long table of pathway, gene, and direction (Up/Down) from enrichment output
path_genes_long <- enrichment_df %>%
  dplyr::mutate(
    up_genes   = Up_regulated,
    down_genes = Down_regulated
  ) %>%
  tidyr::unite("genes_str", Up_regulated, Down_regulated, sep = ",") %>%
  dplyr::rename(pathway = "Term_Description") %>%
  tidyr::separate_rows(up_genes,   sep = "[,;]") %>%
  tidyr::separate_rows(down_genes, sep = "[,;]") %>%
  tidyr::separate_rows(genes_str,  sep = "[,;]") %>%
  dplyr::mutate(
    Gene_Symbol = stringr::str_trim(genes_str),
    direction   = dplyr::case_when(
      genes_str %in% up_genes   ~ "Up",
      genes_str %in% down_genes ~ "Down",
      TRUE                      ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
  dplyr::select(pathway, Gene_Symbol, direction) %>%
  dplyr::distinct() %>%
  tidyr::drop_na()


## --------------------------------------------------
## 6) TCGA-BRCA expression, CNA, survival
## --------------------------------------------------

# Define TCGA input files for expression, CNA, and survival
tcga_expr_file     <- file.path(tcga_dir, "TCGA-BRCA_exprs.z.tsv")
tcga_cna_file      <- file.path(tcga_dir, "TCGA-BRCA_CNA.tsv")
tcga_survival_file <- file.path(tcga_dir, "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv")

# Read TCGA-BRCA expression, CNA, and survival tables
tcga_expr     <- read_tsv(tcga_expr_file,     show_col_types = FALSE)
tcga_cna      <- read_tsv(tcga_cna_file,      show_col_types = FALSE)
tcga_survival <- read_tsv(tcga_survival_file, show_col_types = FALSE)

# Binarize TCGA CNA into 0 vs 1
tcga_cna[, -1] <- ifelse(tcga_cna[, -1] != 0, 1, 0)

# Map Entrez IDs to gene symbols for TCGA expression and CNA tables
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

# Standardize gene symbol column names across TCGA expression and CNA tables
colnames(tcga_expr)[colnames(tcga_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(tcga_cna)[colnames(tcga_cna)   == "gene"]     <- "Gene_Symbol"

# Build a compact TCGA survival table (OS in days; converted later when plotting)
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

## --------------------------------------------------
## 7) TCGA-BRCA expression & CNA long format
## --------------------------------------------------

# Convert TCGA CNA matrix to long format (sample-level rows per gene)
cna_long_tcga <- tcga_cna %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  dplyr::mutate(
    sample = substr(sample_barcode, 1, 15)
  ) %>%
  dplyr::select(Gene_Symbol, sample, cna)

# Compute TCGA pathway CNA score per sample using the same OR-driven direction rule
cna_path_long_tcga <- cna_long_tcga %>%
  dplyr::inner_join(path_genes_long, by = "Gene_Symbol", relationship = "many-to-many") %>%
  dplyr::inner_join(fisher_by_gene,  by = "Gene_Symbol", relationship = "many-to-many") %>%
  dplyr::mutate(
    dir_cna = dplyr::case_when(
      !is.na(OR) & OR > 1 ~  1,
      !is.na(OR) & OR < 1 ~ -1,
      TRUE               ~ NA_real_
    ),
    cna_contrib = dplyr::case_when(
      is.na(dir_cna) ~ NA_real_,
      dir_cna ==  1  ~ cna,
      dir_cna == -1  ~ 1 - cna
    )
  )

path_cna_score_tcga_sample <- cna_path_long_tcga %>%
  dplyr::group_by(sample, pathway) %>%
  dplyr::summarise(
    cna     = mean(cna_contrib, na.rm = TRUE),
    n_genes = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_genes > 1)

# Return METABRIC survival + pathway CNA score for a single pathway
get_pathway_surv_data <- function(pathway_name) {
  path_cna_score_metabric_sample %>%
    dplyr::filter(.data$pathway == !!pathway_name) %>%
    dplyr::select(sample, cna) %>%
    dplyr::inner_join(metabric_surv_small, by = "sample") %>%
    dplyr::mutate(
      cna     = as.numeric(cna),
      OS_TIME = as.numeric(OS_TIME),
      OS_EVENT = as.numeric(OS_EVENT)
    ) %>%
    dplyr::filter(!is.na(OS_TIME), !is.na(OS_EVENT)) %>%
    dplyr::distinct()
}

# Return TCGA survival + pathway CNA score for a single pathway (OS in days here)
get_pathway_surv_data_tcga_OS <- function(pathway_name) {
  path_cna_score_tcga_sample %>%
    dplyr::filter(.data$pathway == !!pathway_name) %>%
    dplyr::select(sample, cna) %>%
    dplyr::inner_join(tcga_surv_small, by = "sample") %>%
    dplyr::mutate(
      cna      = as.numeric(cna),
      OS_TIME  = as.numeric(OS_TIME),
      OS_EVENT = as.numeric(OS_EVENT)
    ) %>%
    dplyr::filter(!is.na(OS_TIME), !is.na(OS_EVENT)) %>%
    dplyr::distinct()
}

## --------------------------------------------------
## 8) METABRIC pathway scores (expression & CNA)
## --------------------------------------------------

# Compute signed expression contributions per pathway using gene direction (Up/Down)
expr_path_metabric <- expr_long_metabric %>%
  dplyr::inner_join(path_genes_long, by = "Gene_Symbol", relationship = "many-to-many") %>%
  dplyr::mutate(
    dir_factor = dplyr::case_when(
      direction == "Up"   ~  1,
      direction == "Down" ~ -1,
      TRUE                ~ NA_real_
    ),
    expr_contrib = abs(expr) * dir_factor
  )

# Aggregate expression pathway scores per sample as the mean contribution across pathway genes
path_expr_score_metabric_sample <- expr_path_metabric %>%
  dplyr::group_by(sample, pathway) %>%
  dplyr::summarise(
    expr    = mean(expr_contrib, na.rm = TRUE),
    n_genes = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_genes > 1)

# Compute CNA pathway contributions using Fisher OR directionality as the sign definition
cna_path_long_metabric <- cna_long_metabric %>%
  dplyr::inner_join(path_genes_long, by = "Gene_Symbol", relationship = "many-to-many") %>%
  dplyr::inner_join(fisher_by_gene,  by = "Gene_Symbol", relationship = "many-to-many") %>%
  dplyr::mutate(
    dir_cna = dplyr::case_when(
      !is.na(OR) & OR > 1 ~  1,
      !is.na(OR) & OR < 1 ~ -1,
      TRUE               ~ NA_real_
    ),
    cna_contrib = dplyr::case_when(
      is.na(dir_cna) ~ NA_real_,
      dir_cna ==  1  ~ cna,
      dir_cna == -1  ~ 1 - cna
    )
  )

# Aggregate CNA pathway scores per sample as the mean CNA contribution across pathway genes
path_cna_score_metabric_sample <- cna_path_long_metabric %>%
  dplyr::group_by(sample, pathway) %>%
  dplyr::summarise(
    cna = mean(cna_contrib, na.rm = TRUE),
    .groups = "drop"
  )

############################################################
## 9) Forest helpers: per-dataset + meta HRs for a pathway
############################################################

# Load patchwork for combining ggplots into multi-panel layouts
suppressPackageStartupMessages(library(patchwork))

# Extract per-dataset HRs and compute fixed-effect meta-analysis HRs for a given pathway
collect_hr_and_meta_for_pathway <- function(pathway) {
  
  ## Per-sheet rows for METABRIC
  rows_met <- purrr::map_dfr(met_filter_sheets, function(sh) {
    df <- readxl::read_xlsx(cox_file_met, sheet = sh)
    if (!"pathway" %in% names(df) || !"HR_cna" %in% names(df)) return(NULL)
    row <- df %>% dplyr::filter(.data$pathway == !!pathway)
    if (nrow(row) == 0 || is.na(row$HR_cna[1])) return(NULL)
    tibble::tibble(
      dataset = "METABRIC",
      sheet   = sh,
      model   = pretty_model_label(sh),
      HR      = as.numeric(row$HR_cna[1]),
      LCL     = as.numeric(row$CI_lo_cna[1]),
      UCL     = as.numeric(row$CI_hi_cna[1])
    )
  })
  
  ## Per-sheet rows for TCGA-BRCA
  rows_tcga <- purrr::map_dfr(tcga_filter_sheets, function(sh) {
    df <- readxl::read_xlsx(cox_file_tcga, sheet = sh)
    if (!"pathway" %in% names(df) || !"HR_cna" %in% names(df)) return(NULL)
    row <- df %>% dplyr::filter(.data$pathway == !!pathway)
    if (nrow(row) == 0 || is.na(row$HR_cna[1])) return(NULL)
    tibble::tibble(
      dataset = "TCGA-BRCA",
      sheet   = sh,
      model   = pretty_model_label(sh),
      HR      = as.numeric(row$HR_cna[1]),
      LCL     = as.numeric(row$CI_lo_cna[1]),
      UCL     = as.numeric(row$CI_hi_cna[1])
    )
  })
  
  # Merge per-dataset rows and stop early if nothing is available
  all_rows <- dplyr::bind_rows(rows_met, rows_tcga)
  if (nrow(all_rows) == 0) return(NULL)
  
  ## log(HR) and SE for meta-analysis
  all_rows <- all_rows %>%
    dplyr::mutate(
      yi  = log(HR),
      sei = (log(UCL) - log(LCL)) / (2 * 1.96)
    )
  
  ## Fixed-effect meta-analysis per model
  meta_df <- all_rows %>%
    dplyr::group_by(model) %>%
    dplyr::group_modify(function(d, key) {
      if (nrow(d) == 1 || any(is.na(d$sei))) {
        tibble::tibble(
          HR_meta  = d$HR[1],
          LCL_meta = d$LCL[1],
          UCL_meta = d$UCL[1]
        )
      } else {
        fit <- tryCatch(
          metafor::rma.uni(yi = d$yi, sei = d$sei, method = "FE"),
          error = function(e) NULL
        )
        if (is.null(fit)) {
          tibble::tibble(
            HR_meta  = d$HR[1],
            LCL_meta = d$LCL[1],
            UCL_meta = d$UCL[1]
          )
        } else {
          tibble::tibble(
            HR_meta  = exp(as.numeric(fit$b)),
            LCL_meta = exp(fit$ci.lb),
            UCL_meta = exp(fit$ci.ub)
          )
        }
      }
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      dataset = "Meta Analysis"
    ) %>%
    dplyr::rename(
      HR  = HR_meta,
      LCL = LCL_meta,
      UCL = UCL_meta
    )
  
  ## Combine per-dataset + meta rows
  out_df <- all_rows %>%
    dplyr::select(dataset, model, HR, LCL, UCL) %>%
    dplyr::bind_rows(meta_df %>% dplyr::select(dataset, model, HR, LCL, UCL))
  
  out_df
}

############################################################
## 10) Forest plot panel (b): facet by dataset
############################################################

# Build a faceted forest plot showing HRs per dataset (including a meta-analysis facet)
forest_facet_for_pathway <- function(pathway, hr_df) {
  if (is.null(hr_df) || nrow(hr_df) == 0) return(NULL)
  
  df <- hr_df %>%
    dplyr::mutate(
      model   = factor(model, levels = rev(unique(model))),
      dataset = factor(
        dataset,
        levels = c("METABRIC", "TCGA-BRCA", "Meta Analysis")
      ),
      ## Text label: always 3 decimal places
      HR_label = formatC(HR, format = "f", digits = 2)
    )
  
  ## For x-axis limits and breaks
  x_min <- min(df$LCL, na.rm = TRUE)
  x_max <- max(df$UCL, na.rm = TRUE)
  
  ggplot(df, aes(x = HR, y = model)) +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    
    ## HR labels: right above each point
    geom_text(
      aes(x = HR, y = model, label = HR_label),
      vjust    = -0.8,
      hjust    = 0.5,
      fontface = "bold",
      size     = 4
    ) +
    
    ## x-axis: fixed scale, not log
    scale_x_continuous(
      limits = c(x_min, x_max * 1.25),              # extra room on the right
      breaks = scales::pretty_breaks(n = 5),
      labels = function(x) formatC(x, format = "f", digits = 2)
    ) +
    
    xlab("HR") +
    ylab("") +
    facet_wrap(~dataset, nrow = 1) +
    ggtitle("Ras signaling pathway - CNA (CoxPH regression models)") +
    theme_bw(base_size = 10) +
    theme(
      strip.text   = element_text(face = "bold", size = 12),
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y  = element_text(size = 12, face = "bold"),
      axis.text.x  = element_text(size = 12),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.ticks.x = element_line(linewidth = 0.8),
      plot.margin  = margin(5, 10, 5, 5),
      panel.spacing.x = grid::unit(1.2, "cm"),  
      panel.spacing.y = grid::unit(0.4, "cm")    
    )
}

############################################################
## 11) KM helpers: METABRIC & TCGA-BRCA (CNA, median split)
##    - No hard p-value filter
##    - Only minimal sample / event checks
##    - Legend extracted later, not inside ggsurvplot
############################################################

# Generate a METABRIC KM plot for pathway CNA score using a median split
km_metabric_cna_median <- function(pathway, HR, LCL, UCL) {
  d <- get_pathway_surv_data(pathway) %>%
    dplyr::filter(!is.na(cna))
  
  if (nrow(d) < 10 || sum(d$OS_EVENT, na.rm = TRUE) < 3) return(NULL)
  
  med <- stats::median(d$cna, na.rm = TRUE)
  
  d <- d %>%
    dplyr::mutate(
      cna_group = factor(
        ifelse(cna > med, "High CNA Burden", "Low CNA Burden"),
        levels = c("Low CNA Burden", "High CNA Burden")
      )
    )
  
  if (length(unique(d$cna_group)) < 2) return(NULL)
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  subtitle_txt <- sprintf("HR = %.2f, 95%% CI: [%.2f – %.2f]", HR, LCL, UCL)
  
  g <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    risk.table.y.text    = TRUE,
    risk.table.col       = NULL,
    risk.table.title     = "Number at risk",
    risk.table.height    = 0.45,
    risk.table.fontsize  = 4.5,
    pval                 = paste0("log-rank p = ", signif(lr_p, 3)),
    legend.title         = "",
    legend.labs          = levels(d$cna_group),
    title                = "Ras signaling pathway - CNA (METABRIC, OS)",
    xlab                 = "Time (months)",
    ylab                 = "OS probability",
    ## >>> All axis settings go here <<<
    ggtheme = theme_bw() +
      theme(
        axis.title.x      = element_text(size = 13, face = "bold"),
        axis.title.y      = element_text(size = 13, face = "bold"),
        axis.text.x       = element_text(size = 12),
        axis.text.y       = element_text(size = 12),
        axis.ticks.length = grid::unit(0.3, "lines"),
        plot.margin       = margin(-2, 5, -2, 5),
        panel.spacing     = grid::unit(0, "pt")
      ),
    tables.theme = theme_bw(base_size = 12)
  )
  
  ## Subsequent theme only for title/subtitle; no axis.text.* here
  g$plot <- g$plot +
    labs(subtitle = subtitle_txt) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )
  
  if (!is.null(g$table)) {
    
    ## First set axis fonts
    g$table <- g$table +
      theme(
        axis.text.x  = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
    
    ## Remove colour mappings from all layers and force black
    for (i in seq_along(g$table$layers)) {
      # drop colour mappings coming from aesthetics
      g$table$layers[[i]]$mapping$colour <- NULL
      g$table$layers[[i]]$mapping$fill   <- NULL
      
      # force black for lines / text
      g$table$layers[[i]]$aes_params$colour <- "black"
      g$table$layers[[i]]$aes_params$fill   <- "black"
    }
  }
  
  list(obj = g, p_lr = lr_p)
}

# Generate a TCGA-BRCA KM plot for pathway CNA score using a median split (OS converted to months)
km_tcga_cna_median <- function(pathway, HR, LCL, UCL) {
  d <- get_pathway_surv_data_tcga_OS(pathway) %>%
    dplyr::filter(!is.na(cna)) %>%
    dplyr::mutate(OS_TIME_m = OS_TIME) 
  
  if (nrow(d) < 10 || sum(d$OS_EVENT, na.rm = TRUE) < 3) return(NULL)
  
  med <- stats::median(d$cna, na.rm = TRUE)
  
  d <- d %>%
    dplyr::mutate(
      cna_group = factor(
        ifelse(cna > med, "High CNA Burden", "Low CNA Burden"),
        levels = c("Low CNA Burden", "High CNA Burden")
      )
    )
  
  if (length(unique(d$cna_group)) < 2) return(NULL)
  
  fit  <- survfit(Surv(OS_TIME_m, OS_EVENT) ~ cna_group, data = d)
  lr   <- survdiff(Surv(OS_TIME_m, OS_EVENT) ~ cna_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  subtitle_txt <- sprintf("HR = %.2f, 95%% CI: [%.2f – %.2f]", HR, LCL, UCL)
  
  g <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    risk.table.y.text    = TRUE,
    risk.table.col       = NULL,
    risk.table.title     = "Number at risk",
    risk.table.height    = 0.45,
    risk.table.fontsize  = 4.5,
    pval                 = paste0("log-rank p = ", signif(lr_p, 3)),
    legend.title         = "",
    legend.labs          = levels(d$cna_group),
    title                = "Ras signaling pathway - CNA (TCGA-BRCA, OS)",
    xlab                 = "Time (months)",
    ylab                 = "OS probability",
    ggtheme = theme_bw() +
      theme(
        axis.title.x      = element_text(size = 13, face = "bold"),
        axis.title.y      = element_text(size = 13, face = "bold"),
        axis.text.x       = element_text(size = 12),
        axis.text.y       = element_text(size = 12),
        axis.ticks.length = grid::unit(0.3, "lines"),
        plot.margin       = margin(-2, 5, -2, 5),
        panel.spacing     = grid::unit(0, "pt")
      ),
    tables.theme = theme_bw(base_size = 12)
  )
  
  g$plot <- g$plot +
    labs(subtitle = subtitle_txt) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "bottom"
    )
  
  if (!is.null(g$table)) {
    
    ## First set axis fonts
    g$table <- g$table +
      theme(
        axis.text.x  = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
    
    ## Remove colour mappings from all layers and force black
    for (i in seq_along(g$table$layers)) {
      # drop colour mappings coming from aesthetics
      g$table$layers[[i]]$mapping$colour <- NULL
      g$table$layers[[i]]$mapping$fill   <- NULL
      
      # force black for lines / text
      g$table$layers[[i]]$aes_params$colour <- "black"
      g$table$layers[[i]]$aes_params$fill   <- "black"
    }
  }
  
  list(obj = g, p_lr = lr_p)
}

############################################################
## 12) Build Figure 4 (KM left + forest right) and save
############################################################

# Define Cox result file paths used for both joint HR extraction and forest/meta helpers
cox_file_met  <- file.path(metabric_results_dir, "METABRIC_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx")
cox_file_tcga <- file.path(tcga_results_dir,     "TCGA_BRCA_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx")

# Define which sheets to scan for forest/meta panels
met_filter_sheets  <- excel_sheets(cox_file_met)
tcga_filter_sheets <- excel_sheets(cox_file_tcga)

# Pick the full-model sheet used to populate KM subtitles (single model per cohort)
full_model_sheet <- "E_C_Age_Stg_PAM50_filter"

# Build a joint pathway table (METABRIC + TCGA) with the exact column names used downstream
cna_stats_joint <- {
  df_met  <- readxl::read_xlsx(cox_file_met,  sheet = full_model_sheet)
  df_tcga <- readxl::read_xlsx(cox_file_tcga, sheet = full_model_sheet)
  
  if (!"pathway" %in% names(df_met))  stop("METABRIC Cox sheet lacks 'pathway' column: ", full_model_sheet)
  if (!"pathway" %in% names(df_tcga)) stop("TCGA Cox sheet lacks 'pathway' column: ", full_model_sheet)
  
  req_cols <- c("HR_cna", "CI_lo_cna", "CI_hi_cna")
  if (!all(req_cols %in% names(df_met)))  stop("METABRIC Cox sheet lacks CNA columns in: ", full_model_sheet)
  if (!all(req_cols %in% names(df_tcga))) stop("TCGA Cox sheet lacks CNA columns in: ", full_model_sheet)
  
  dplyr::full_join(
    df_met  %>%
      dplyr::select(pathway, HR_cna, CI_lo_cna, CI_hi_cna) %>%
      dplyr::rename(
        HR_cna_met = HR_cna,
        LCL_met    = CI_lo_cna,
        UCL_met    = CI_hi_cna
      ),
    df_tcga %>%
      dplyr::select(pathway, HR_cna, CI_lo_cna, CI_hi_cna) %>%
      dplyr::rename(
        HR_cna_tcga = HR_cna,
        LCL_tcga    = CI_lo_cna,
        UCL_tcga    = CI_hi_cna
      ),
    by = "pathway"
  )
}

# Define the target pathway to visualize across cohorts and models
target_pathway <- "Ras signaling pathway"

## Joint HR row from full model (E_C_Age_Stg_PAM50_filter)
stats_row <- cna_stats_joint %>%
  dplyr::filter(.data$pathway == target_pathway) %>%
  dplyr::slice(1)

## KM panels
km_met  <- km_metabric_cna_median(
  pathway = target_pathway,
  HR      = stats_row$HR_cna_met,
  LCL     = stats_row$LCL_met,
  UCL     = stats_row$UCL_met
)
km_tcga <- km_tcga_cna_median(
  pathway = target_pathway,
  HR      = stats_row$HR_cna_tcga,
  LCL     = stats_row$LCL_tcga,
  UCL     = stats_row$UCL_tcga
)

## If one cohort still fails, just warn and skip
if (is.null(km_met) || is.null(km_tcga)) {
  warning("KM plot could not be generated for one of the cohorts (METABRIC or TCGA-BRCA).")
} else {
  ## 1) Remove legends from internal KM plots
  km_met$obj$plot  <- km_met$obj$plot  + theme(legend.position = "none")
  km_tcga$obj$plot <- km_tcga$obj$plot + theme(legend.position = "none")
  
  km_met$obj$plot <- km_met$obj$plot +
    theme(
      plot.margin = margin(t = 1, r = 0, b = 0, l = 0)
    )
  
  km_tcga$obj$plot <- km_tcga$obj$plot +
    theme(
      plot.margin = margin(t = 0, r = 0, b = 1, l = 0)
    )
  
  ## 2) Extract combined legend from TCGA plot
  tcga_legend <- cowplot::get_legend(
    km_tcga$obj$plot +
      theme(
        legend.position = "bottom",
        legend.title    = element_text(face = "bold", size = 12),
        legend.text     = element_text(face = "bold", size = 12),
        legend.key.size = unit(0.9, "cm"),          # increase symbol size
        legend.spacing  = unit(0.4, "cm")           # increase spacing between items
      )
  )
  
  ## 3) Construct left panel in vertical order:
  ##    METABRIC (plot + table) -> spacer -> TCGA (plot + table) -> legend at bottom
  ## 1) Turn KM objects into a single grob each
  km_to_grob <- function(km_obj){
    gridExtra::arrangeGrob(
      km_obj$plot,
      km_obj$table,
      ncol = 1,
      heights = c(0.75, 0.25)
    )
  }
  
  ## 2) Spacer grob
  blank_mid <- grid::grobTree(
    grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
  )
  
  ## 3) Left panel (a)
  left_panel <- gridExtra::arrangeGrob(
    km_to_grob(km_met$obj),
    blank_mid,
    km_to_grob(km_tcga$obj),
    tcga_legend,
    ncol = 1,
    heights = c(0.78, 0.04, 0.78, 0.10)
  )
  
  ## 4) Full-bleed style for panel a
  left_gg <- ggplot() +
    annotation_custom(left_panel) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)  # remove all outer margins
    )
  
  ## 5) Forest data and panel b (side-by-side facets)
  hr_df    <- collect_hr_and_meta_for_pathway(target_pathway)
  right_gg <- forest_facet_for_pathway(target_pathway, hr_df) +
    scale_x_continuous(
      labels = function(x) formatC(x, format = "f", digits = 2)
    ) +
    ggtitle("Ras signaling pathway - CNA (CoxPH regression models)") +
    theme(
      axis.text.y   = element_text(size = 12, face = "bold"),
      strip.text    = element_text(size = 12, face = "bold"),
      axis.title.x  = element_text(size = 12, face = "bold"),
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x   = element_text(size = 12),
      
      ## Forest panel margin — mainly top–bottom
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ## 6) Convert left grob to ggplot (for patchwork)
  left_gg <- ggplot() +
    annotation_custom(left_panel) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0)
    )
  
  ## 7) Combine panels with patchwork: a (left, 25%) | b (right, 75%)
  combined <- (left_gg | right_gg) +
    patchwork::plot_layout(widths = c(1, 1.8)) +
    patchwork::plot_annotation(
      tag_levels = "A",
      tag_suffix = "."
    )
  
  ## Apply tag style at the very end (highest priority)
  combined <- combined & theme(
    plot.tag = element_text(face = "bold", size = 16)
  )
  
  ## 8) Save as JPEG
  out_file <- file.path(figures_dir, "Figure 4.jpeg")
  ggsave(out_file, combined, width = 20, height = 11, dpi = 300)
  message("Combined KM + forest figure saved to: ", out_file)
}