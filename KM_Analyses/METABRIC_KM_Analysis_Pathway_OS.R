## --------------------------------------------------
## 0) Packages
## --------------------------------------------------

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(writexl))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(AnnotationDbi))

## --------------------------------------------------
## 1) Directories
## --------------------------------------------------

metabric_dir         <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
pathfindR_dir        <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"
gdsc_dir             <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
metabric_results_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/METABRIC"

if (!dir.exists(metabric_results_dir)) {
  dir.create(metabric_results_dir, recursive = TRUE)
}

## Output folder for KM plots (OS, all METABRIC, PATHWAY-level)
km_dir <- file.path(metabric_results_dir, "METABRIC_KM_Pathway_OS")
if (!dir.exists(km_dir)) {
  dir.create(km_dir, recursive = TRUE)
}

## Helper to clean pathway names for filenames
clean_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^[:alnum:]]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

## --------------------------------------------------
## 2) Read METABRIC expression, CNA, survival
## --------------------------------------------------

metabric_expr_file     <- file.path(metabric_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
metabric_cna_file      <- file.path(metabric_dir, "data_cna.txt")
metabric_survival_file <- file.path(metabric_dir, "brca_metabric_clinical_data.tsv")

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

## Make CNA binary (0 vs !=0) → 0 / 1
metabric_cna_binary <- metabric_cna
metabric_cna_binary[, -1] <- ifelse(metabric_cna_binary[, -1] != 0, 1, 0)

## --------------------------------------------------
## 3) METABRIC OS (ALL patients) + covariates (Age, Stage_simple, PAM50)
## --------------------------------------------------

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

message("Number of METABRIC samples with OS data: ", nrow(metabric_surv_small))

## --------------------------------------------------
## 4) Expression & CNA in long format (per gene; used for pathway scoring)
## --------------------------------------------------

expr_long_metabric <- metabric_expr %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

cna_long_metabric <- metabric_cna_binary %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

## --------------------------------------------------
## 5) GDSC data for Fisher OR (gene-level CNA vs response)
## --------------------------------------------------

gdsc_expr_file     <- file.path(gdsc_dir, "GDSC_exprs.Paclitaxel.eb_with.TCGA_exprs.Paclitaxel.tsv")
gdsc_cna_file      <- file.path(gdsc_dir, "GDSC_CNA.Paclitaxel.tsv")
gdsc_response_file <- file.path(gdsc_dir, "GDSC_response.Paclitaxel.tsv")

gdsc_expr     <- read_tsv(gdsc_expr_file,     show_col_types = FALSE)
gdsc_cna      <- read_tsv(gdsc_cna_file,      show_col_types = FALSE)
gdsc_response <- read_tsv(gdsc_response_file, show_col_types = FALSE)

## Make GDSC CNA binary (0 vs !=0)
gdsc_cna[, -1] <- ifelse(gdsc_cna[, -1] != 0, 1, 0)

## Map Entrez IDs to gene symbols in GDSC
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

colnames(gdsc_expr)[colnames(gdsc_expr) == "ENTREZID"] <- "Gene_Symbol"
colnames(gdsc_cna)[colnames(gdsc_cna) == "gene_id"]   <- "Gene_Symbol"

## GDSC samples and response
resp <- gdsc_response %>%
  dplyr::mutate(status = response)

gdsc_samples <- as.character(unique(resp$sample_name))

## GDSC expression: subset to GDSC samples and z-score by gene
expr_gdsc <- gdsc_expr %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

expr_mat <- as.matrix(expr_gdsc[, -1])
rownames(expr_mat) <- expr_gdsc$Gene_Symbol
expr_z <- t(scale(t(expr_mat)))
expr_z[is.na(expr_z)] <- 0
expr_z_df <- as_tibble(expr_z, rownames = "Gene_Symbol")

## GDSC CNA: subset to GDSC samples, already binary
cna_gdsc <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

cna_mat <- as.matrix(cna_gdsc[, -1])
rownames(cna_mat) <- cna_gdsc$Gene_Symbol
cna_bin <- ifelse(cna_mat != 0, 1, 0)
cna_bin_df <- as_tibble(cna_bin, rownames = "Gene_Symbol")

## GDSC expression in long format with response
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

## GDSC CNA in long format with response
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

## Fisher OR per gene (GDSC CNA vs response)
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
## 6) PathfindR enrichment and pathway–gene mapping
## --------------------------------------------------

enrichment_file <- file.path(pathfindR_dir, "pathfindR_paclitaxel.tsv")
enrichment_df   <- read_tsv(enrichment_file, show_col_types = FALSE)

selected_terms <- c(
  "hsa04520", "hsa05146", "hsa04371", "hsa04360", "hsa04020",
  "hsa04218", "hsa05230", "hsa05207", "hsa05231", "hsa04820",
  "hsa04512", "hsa01521", "hsa04915", "hsa04510", "hsa04068",
  "hsa04540", "hsa04066", "hsa04390", "hsa04010", "hsa04916",
  "hsa05206", "hsa04814", "hsa04650", "hsa04080", "hsa04921",
  "hsa05235", "hsa04151", "hsa05200", "hsa05205", "hsa00620",
  "hsa04015", "hsa04014", "hsa04810", "hsa04926", "hsa04924",
  "hsa04668", "hsa04530", "hsa05202", "hsa04310", "hsa04024",
  "hsa04022", "hsa04115"
)

enrichment_df <- enrichment_df %>%
  dplyr::filter(ID %in% selected_terms)

## Build pathway–gene-direction table
path_genes_long <- enrichment_df %>%
  dplyr::mutate(
    up_genes   = Up_regulated,
    down_genes = Down_regulated
  ) %>%
  tidyr::unite("genes_str", Up_regulated, Down_regulated, sep = ",") %>%
  dplyr::rename(
    pathway = "Term_Description"
  ) %>%
  tidyr::separate_rows(up_genes,   sep = "[,;]") %>%
  tidyr::separate_rows(down_genes, sep = "[,;]") %>%
  tidyr::separate_rows(genes_str,  sep = "[,;]") %>%
  dplyr::mutate(
    Gene_Symbol = stringr::str_trim(genes_str),
    direction   = dplyr::case_when(
      genes_str %in% up_genes   ~ "Up",
      genes_str %in% down_genes ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
  dplyr::select(pathway, Gene_Symbol, direction) %>%
  dplyr::distinct() %>%
  tidyr::drop_na()

## --------------------------------------------------
## 7) Pathway scores in METABRIC (expression & CNA)
## --------------------------------------------------

## 7a) Pathway expression scores (sample × pathway)
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

path_expr_score_metabric <- expr_path_metabric %>%
  dplyr::group_by(sample, pathway) %>%
  dplyr::summarise(
    expr    = mean(expr_contrib, na.rm = TRUE),
    n_genes = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_genes > 1)

## 7b) Pathway CNA scores (sample × pathway), using GDSC Fisher OR for direction
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
      dir_cna ==  1  ~ cna,     ## CNA presence associated with worse outcome
      dir_cna == -1  ~ 1 - cna  ## CNA absence associated with worse outcome
    )
  )

path_cna_score_metabric <- cna_path_long_metabric %>%
  dplyr::group_by(sample, pathway) %>%
  dplyr::summarise(
    cna = mean(cna_contrib, na.rm = TRUE),
    .groups = "drop"
  )

## --------------------------------------------------
## 8) Helper: get per-pathway survival + pathway scores (OS, ALL patients)
## --------------------------------------------------
## NOTE: Unlike gene-level get_gene_surv_data, we do NOT require both expr and cna
## in joins. We left_join expression and CNA pathway scores to survival and then
## each KM function filters on the columns it needs.

get_pathway_surv_data <- function(pathway_name) {
  metabric_surv_small %>%
    dplyr::left_join(
      path_expr_score_metabric %>%
        dplyr::filter(pathway == pathway_name) %>%
        dplyr::select(sample, expr_score = expr),
      by = "sample"
    ) %>%
    dplyr::left_join(
      path_cna_score_metabric %>%
        dplyr::filter(pathway == pathway_name) %>%
        dplyr::select(sample, cna_score = cna),
      by = "sample"
    ) %>%
    dplyr::mutate(
      expr = as.numeric(expr_score),
      cna  = as.numeric(cna_score)
    ) %>%
    dplyr::filter(
      !is.na(OS_TIME),
      !is.na(OS_EVENT)
    ) %>%
    dplyr::distinct()
}

## --------------------------------------------------
## 9) KM plotting helpers (PATHWAY-LEVEL, OS ALL PATIENTS)
## --------------------------------------------------

## 9a) Pathway CNA-based KM (2 groups, median cut on CNA score)
plot_km_pathway_cna <- function(pathway_name,
                                HR, LCL, UCL, PVAL,
                                out_dir,
                                model_suffix,
                                require_sig_p = TRUE) {
  d <- get_pathway_surv_data(pathway_name) %>%
    dplyr::filter(!is.na(cna))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("Pathway CNA KM (OS all): insufficient data (n<30 or <5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  med <- stats::median(d$cna, na.rm = TRUE)
  
  d <- d %>%
    dplyr::mutate(
      cna_group = factor(
        ifelse(cna > med, "High CNA burden", "Low CNA burden"),
        levels = c("Low CNA burden", "High CNA burden")
      )
    )
  
  if (length(unique(d$cna_group)) < 2) {
    message("Pathway CNA KM (OS all): only one CNA group present, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  if (require_sig_p && (is.na(lr_p) || lr_p >= 0.05)) {
    message("Pathway CNA KM (OS all): log-rank p >= 0.05, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  subtitle_txt <- sprintf(
    "HR = %.2f, 95%% CI: [%.2f – %.2f]",
    HR, LCL, UCL
  )
  
  p <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    risk.table.y.text    = TRUE,
    risk.table.col       = NULL,
    risk.table.title     = "Number at risk",
    risk.table.height    = 0.25,
    risk.table.fontsize  = 3,
    pval                 = paste0("log-rank p = ", signif(lr_p, 3)),
    legend.title         = "",
    legend.labs          = levels(d$cna_group),
    title                = paste0(pathway_name, " - CNA (OS, METABRIC pathways)"),
    xlab                 = "Time (months)",
    ylab                 = "Probability",
    ggtheme              = theme_bw(),
    tables.theme         = theme_bw()
  )
  
  p$plot <- p$plot +
    labs(subtitle = subtitle_txt) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  if (!is.null(p$table)) {
    if ("strata" %in% names(p$table$data)) {
      old_levels <- levels(p$table$data$strata)
      p$table$data$strata <- factor(
        p$table$data$strata,
        levels = old_levels,
        labels = levels(d$cna_group)[seq_along(old_levels)]
      )
    }
    
    p$table <- p$table +
      theme(
        strip.background = element_blank(),
        strip.text       = element_blank(),
        axis.title.y     = element_blank(),
        axis.title.x     = element_blank()
      )
  }
  
  file_tag  <- paste0("KM_OS_pathway_CNA_", clean_filename(pathway_name), model_suffix, ".pdf")
  out_file  <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 9, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved METABRIC OS PATHWAY CNA KM plot: ", out_file)
  invisible(out_file)
}

## 9b) Pathway expression – median cut (2 groups on expr score)
plot_km_pathway_expr_median <- function(pathway_name,
                                        HR, LCL, UCL, PVAL,
                                        out_dir,
                                        model_suffix,
                                        require_sig_p = TRUE) {
  d <- get_pathway_surv_data(pathway_name) %>%
    dplyr::filter(!is.na(expr))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("Pathway expr-median KM (OS all): insufficient data (n<30 or <5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  med <- stats::median(d$expr, na.rm = TRUE)
  
  d <- d %>%
    dplyr::mutate(
      expr_group = factor(
        ifelse(expr > med, "Active", "Suppressed"),
        levels = c("Suppressed", "Active")
      )
    )
  
  if (length(unique(d$expr_group)) < 2) {
    message("Pathway expr-median KM (OS all): only one expression group present, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  if (require_sig_p && (is.na(lr_p) || lr_p >= 0.05)) {
    message("Pathway expr-median KM (OS all): log-rank p >= 0.05, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  subtitle_txt <- sprintf(
    "HR = %.2f, 95%% CI: [%.2f – %.2f]",
    HR, LCL, UCL
  )
  
  p <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    risk.table.y.text    = TRUE,
    risk.table.col       = NULL,
    risk.table.title     = "Number at risk",
    risk.table.height    = 0.25,
    risk.table.fontsize  = 3,
    pval                 = paste0("log-rank p = ", signif(lr_p, 3)),
    legend.title         = "",
    legend.labs          = levels(d$expr_group),
    title                = paste0(pathway_name, " - Expression median (OS, METABRIC pathways)"),
    xlab                 = "Time (months)",
    ylab                 = "Probability",
    ggtheme              = theme_bw(),
    tables.theme         = theme_bw()
  )
  
  p$plot <- p$plot +
    labs(subtitle = subtitle_txt) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  if (!is.null(p$table)) {
    if ("strata" %in% names(p$table$data)) {
      old_levels <- levels(p$table$data$strata)
      p$table$data$strata <- factor(
        p$table$data$strata,
        levels = old_levels,
        labels = levels(d$expr_group)[seq_along(old_levels)]
      )
    }
    
    p$table <- p$table +
      theme(
        strip.background = element_blank(),
        strip.text       = element_blank(),
        axis.title.y     = element_blank(),
        axis.title.x     = element_blank()
      )
  }
  
  file_tag <- paste0("KM_OS_pathway_expr_median_", clean_filename(pathway_name), model_suffix, ".pdf")
  out_file <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 9, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved METABRIC OS PATHWAY expression-median KM plot: ", out_file)
  invisible(out_file)
}

## 9c) Pathway expression – optimal cut (maxstat, 2 groups on expr score)
plot_km_pathway_expr_optimal <- function(pathway_name,
                                         HR, LCL, UCL, PVAL,
                                         out_dir,
                                         model_suffix,
                                         require_sig_p = TRUE) {
  d <- get_pathway_surv_data(pathway_name) %>%
    dplyr::filter(!is.na(expr))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("Pathway expr-optimal KM (OS all): insufficient data (n<30 or <5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  tmp <- d %>%
    dplyr::select(OS_TIME, OS_EVENT, expr)
  
  cut_res <- tryCatch(
    {
      surv_cutpoint(
        tmp,
        time      = "OS_TIME",
        event     = "OS_EVENT",
        variables = "expr",
        minprop   = 0.2
      )
    },
    error = function(e) {
      message("Pathway expr-optimal KM (OS all): surv_cutpoint failed for pathway ", pathway_name, " – skipping. Error: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(cut_res)) return(NULL)
  
  cat_res <- surv_categorize(cut_res)
  d$expr_group_opt_raw <- cat_res$expr
  
  d <- d %>%
    dplyr::mutate(
      expr_group_opt = dplyr::recode(
        expr_group_opt_raw,
        "low"  = "Suppressed",
        "high" = "Active",
        .default = NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(expr_group_opt)) %>%
    dplyr::mutate(
      expr_group_opt = factor(expr_group_opt, levels = c("Suppressed", "Active"))
    )
  
  if (length(unique(d$expr_group_opt)) < 2) {
    message("Pathway expr-optimal KM (OS all): only one optimal expression group present, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_group_opt, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_group_opt, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  if (require_sig_p && (is.na(lr_p) || lr_p >= 0.05)) {
    message("Pathway expr-optimal KM (OS all): log-rank p >= 0.05, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  subtitle_txt <- sprintf(
    "HR = %.2f, 95%% CI: [%.2f – %.2f]",
    HR, LCL, UCL
  )
  
  pval_txt <- paste0("log-rank p = ", signif(lr_p, 3))
  
  p <- suppressWarnings(
    ggsurvplot(
      fit,
      data                 = d,
      risk.table           = TRUE,
      risk.table.y.text    = TRUE,
      risk.table.col       = NULL,
      risk.table.title     = "Number at risk",
      risk.table.height    = 0.25,
      risk.table.fontsize  = 3,
      pval                 = pval_txt,
      legend.title         = "",
      legend.labs          = levels(d$expr_group_opt),
      title                = paste0(pathway_name, " - Expression optimal (OS, METABRIC pathways)"),
      xlab                 = "Time (months)",
      ylab                 = "Probability",
      ggtheme              = theme_bw(),
      tables.theme         = theme_bw()
    )
  )
  
  p$plot <- p$plot +
    labs(subtitle = subtitle_txt) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  if (!is.null(p$table)) {
    if ("strata" %in% names(p$table$data)) {
      old_levels <- levels(p$table$data$strata)
      p$table$data$strata <- factor(
        p$table$data$strata,
        levels = old_levels,
        labels = levels(d$expr_group_opt)[seq_along(old_levels)]
      )
    }
    
    p$table <- p$table +
      theme(
        strip.background = element_blank(),
        strip.text       = element_blank(),
        axis.title.y     = element_blank(),
        axis.title.x     = element_blank()
      )
  }
  
  file_tag <- paste0("KM_OS_pathway_expr_optimal_", clean_filename(pathway_name), model_suffix, ".pdf")
  out_file <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 9, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved METABRIC OS PATHWAY expression-optimal KM plot: ", out_file)
  invisible(out_file)
}

## 9d) Pathway expr:cna interaction KM on pathway scores
plot_km_pathway_expr_cna_interaction <- function(pathway_name,
                                                 HR, LCL, UCL, PVAL,
                                                 out_dir,
                                                 model_suffix,
                                                 four_groups = TRUE) {
  d <- get_pathway_surv_data(pathway_name) %>%
    dplyr::filter(!is.na(expr), !is.na(cna))
  
  ## Minimum sample / event check
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("Pathway Expr:CNA KM (OS all): insufficient data (n<30 or <5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  ## ==========================================
  ## 4-GROUP VERSION (expr High/Low x CNA High/Low)
  ## ==========================================
  if (four_groups) {
    
    ## FIXED: median tanımları eksikti
    med_expr <- stats::median(d$expr, na.rm = TRUE)
    
    ## CNA binary olduğu için median gerekmez, doğrudan 0/1:
    d <- d %>%
      dplyr::mutate(
        expr_group = factor(
          ifelse(expr > med_expr, "Active", "Suppressed"),
          levels = c("Suppressed", "Active")
        ),
        cna_group = factor(
          ifelse(cna == 1, "High CNA Burden", "Low CNA Burden"),
          levels = c("Low CNA Burden", "High CNA Burden")
        ),
        int_group = factor(
          paste(expr_group, cna_group, sep = " & "),
          levels = c(
            "Suppressed & Low CNA Burden",
            "Suppressed & High CNA Burden",
            "Active & Low CNA Burden",
            "Active & High CNA Burden"
          )
        )
      )
    
    ## ==========================================
    ## 2-GROUP VERSION (expr*cna interaction score median)
    ## ==========================================
  } else {
    
    d <- d %>%
      dplyr::mutate(int_score = expr * cna)
    
    if (length(unique(na.omit(d$int_score))) < 2) {
      message("Pathway Expr:CNA KM (OS all): interaction score <2 values, skipping pathway: ", pathway_name)
      return(NULL)
    }
    
    med_int <- stats::median(d$int_score, na.rm = TRUE)
    
    d <- d %>%
      dplyr::mutate(
        int_group = factor(
          ifelse(int_score > med_int, "High interaction", "Low interaction"),
          levels = c("Low interaction", "High interaction")
        )
      )
  }
  
  ## Not enough groups
  if (length(na.omit(unique(d$int_group))) < 2) {
    message("Pathway Expr:CNA KM (OS all): fewer than 2 interaction groups present, skipping: ", pathway_name)
    return(NULL)
  }
  
  ## ==========================================
  ## Survival fit
  ## ==========================================
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ int_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ int_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  subtitle_txt <- sprintf(
    "HR = %.2f, 95%% CI: [%.2f – %.2f]",
    HR, LCL, UCL
  )
  
  ## ==========================================
  ## KM Plot
  ## ==========================================
  p <- ggsurvplot(
    fit,
    data                 = d,
    risk.table           = TRUE,
    risk.table.y.text    = TRUE,
    risk.table.col       = NULL,
    risk.table.title     = "Number at risk",
    risk.table.height    = if (four_groups) 0.30 else 0.25,
    risk.table.fontsize  = 3,
    pval                 = paste0("log-rank p = ", signif(lr_p, 3)),
    legend.title         = "",
    legend.labs          = levels(d$int_group),
    title                = paste0(pathway_name, " - Expression × CNA (OS, METABRIC pathways)"),
    xlab                 = "Time (months)",
    ylab                 = "OS probability",
    ggtheme              = theme_bw(),
    tables.theme         = theme_bw()
  )
  
  p$plot <- p$plot +
    labs(subtitle = subtitle_txt) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ## Risk table düzeltme
  if (!is.null(p$table)) {
    if ("strata" %in% names(p$table$data)) {
      p$table$data$strata <- factor(
        p$table$data$strata,
        labels = levels(d$int_group)
      )
    }
    p$table <- p$table +
      theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
      )
  }
  
  ## ==========================================
  ## SAVE FILE
  ## ==========================================
  file_tag <- paste0("KM_OS_pathway_exprCNA_interaction_",
                     clean_filename(pathway_name),
                     model_suffix,
                     ".pdf")
  out_file <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 9, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved METABRIC OS PATHWAY expr:CNA interaction KM plot: ", out_file)
  invisible(out_file)
}

## 9e) Extra diagnostics for pathway expr:CNA interaction
extra_pathway_expr_cna_diagnostics <- function(pathway_name,
                                               out_dir,
                                               model_suffix,
                                               four_groups = TRUE) {
  d <- get_pathway_surv_data(pathway_name) %>%
    dplyr::filter(!is.na(expr), !is.na(cna))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("Extra PATHWAY expr:CNA diagnostics (OS all): insufficient data (n<30 or <5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  summary_lines <- c()
  add_line <- function(...) {
    summary_lines <<- c(summary_lines, paste0(...))
  }
  
  add_line("Pathway: ", pathway_name)
  add_line("Model suffix: ", model_suffix)
  add_line("N total: ", nrow(d), "; Events: ", sum(d$OS_EVENT, na.rm = TRUE))
  add_line("----")
  
  ## 1) Four-group KM summary (Low/High expr x Low/High CNA)
  if (four_groups) {
    med_expr <- stats::median(d$expr, na.rm = TRUE)
    med_cna  <- stats::median(d$cna,  na.rm = TRUE)
    
    d4 <- d %>%
      dplyr::mutate(
        expr_group = factor(
          ifelse(expr > med_expr, "High expr", "Low expr"),
          levels = c("Low expr", "High expr")
        ),
        cna_group = factor(
          ifelse(cna > med_cna, "High CNA", "Low CNA"),
          levels = c("Low CNA", "High CNA")
        ),
        int_group = factor(
          paste(expr_group, cna_group, sep = " & "),
          levels = c(
            "Low expr & Low CNA",
            "Low expr & High CNA",
            "High expr & Low CNA",
            "High expr & High CNA"
          )
        )
      ) %>%
      dplyr::filter(!is.na(int_group))
    
    if (length(unique(na.omit(d4$int_group))) >= 2) {
      tab4 <- d4 %>%
        dplyr::group_by(int_group) %>%
        dplyr::summarise(
          n      = dplyr::n(),
          events = sum(OS_EVENT, na.rm = TRUE),
          .groups = "drop"
        )
      
      fit4  <- survfit(Surv(OS_TIME, OS_EVENT) ~ int_group, data = d4)
      lr4   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ int_group, data = d4)
      lr4_p <- 1 - pchisq(lr4$chisq, df = length(lr4$n) - 1)
      
      add_line("Four-group KM (Low/High expr x Low/High CNA):")
      add_line(capture.output(print(tab4)))
      add_line("Four-group log-rank p: ", signif(lr4_p, 4))
      add_line("----")
    } else {
      add_line("Four-group KM: fewer than 2 groups with data, skipping.")
      add_line("----")
    }
  } else {
    add_line("Four-group KM: not requested (four_groups = FALSE).")
    add_line("----")
  }
  
  ## 2) Two-group KM using interaction score (expr*cna), median split
  d_int <- d %>%
    dplyr::mutate(
      int_score = expr * cna
    )
  
  if (length(unique(na.omit(d_int$int_score))) >= 2) {
    med_int <- stats::median(d_int$int_score, na.rm = TRUE)
    
    d_int <- d_int %>%
      dplyr::mutate(
        int_group_median = factor(
          ifelse(int_score > med_int, "High interaction", "Low interaction"),
          levels = c("Low interaction", "High interaction")
        )
      )
    
    if (length(unique(na.omit(d_int$int_group_median))) >= 2) {
      fit2  <- survfit(Surv(OS_TIME, OS_EVENT) ~ int_group_median, data = d_int)
      lr2   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ int_group_median, data = d_int)
      lr2_p <- 1 - pchisq(lr2$chisq, df = length(lr2$n) - 1)
      
      tab2 <- d_int %>%
        dplyr::group_by(int_group_median) %>%
        dplyr::summarise(
          n      = dplyr::n(),
          events = sum(OS_EVENT, na.rm = TRUE),
          .groups = "drop"
        )
      
      add_line("Two-group KM (interaction score = expr*cna, median cut):")
      add_line("Median interaction score: ", signif(med_int, 4))
      add_line(capture.output(print(tab2)))
      add_line("Two-group log-rank p: ", signif(lr2_p, 4))
      add_line("----")
    } else {
      add_line("Two-group KM (interaction score): only one group present, skipping.")
      add_line("----")
    }
  } else {
    add_line("Interaction score has <2 unique values, skipping 2-group KM.")
    add_line("----")
  }
  
  ## 3) Optimal cutpoint for interaction score (expr*cna)
  if ("int_score" %in% names(d_int) && length(unique(na.omit(d_int$int_score))) >= 3) {
    tmp_int <- d_int %>%
      dplyr::select(OS_TIME, OS_EVENT, int_score)
    
    cut_res_int <- tryCatch(
      {
        surv_cutpoint(
          tmp_int,
          time      = "OS_TIME",
          event     = "OS_EVENT",
          variables = "int_score",
          minprop   = 0.2
        )
      },
      error = function(e) {
        NULL
      }
    )
    
    if (!is.null(cut_res_int)) {
      cat_int <- surv_categorize(cut_res_int)
      d_int_opt <- d_int %>%
        dplyr::mutate(
          int_group_opt_raw = cat_int$int_score,
          int_group_opt = dplyr::recode(
            int_group_opt_raw,
            "low"  = "Low interaction",
            "high" = "High interaction",
            .default = NA_character_
          )
        ) %>%
        dplyr::filter(!is.na(int_group_opt)) %>%
        dplyr::mutate(
          int_group_opt = factor(
            int_group_opt,
            levels = c("Low interaction", "High interaction")
          )
        )
      
      if (length(unique(na.omit(d_int_opt$int_group_opt))) >= 2) {
        fit_opt  <- survfit(Surv(OS_TIME, OS_EVENT) ~ int_group_opt, data = d_int_opt)
        lr_opt   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ int_group_opt, data = d_int_opt)
        lr_opt_p <- 1 - pchisq(lr_opt$chisq, df = length(lr_opt$n) - 1)
        
        tab_opt <- d_int_opt %>%
          dplyr::group_by(int_group_opt) %>%
          dplyr::summarise(
            n      = dplyr::n(),
            events = sum(OS_EVENT, na.rm = TRUE),
            .groups = "drop"
          )
        
        cp_df <- cut_res_int$cutpoint
        if ("cutpoint" %in% colnames(cp_df)) {
          cut_value <- as.numeric(cp_df$cutpoint[1])
        } else if ("int_score" %in% colnames(cp_df)) {
          cut_value <- as.numeric(cp_df$int_score[1])
        } else {
          cut_value <- NA_real_
        }
        
        add_line("Optimal cutpoint KM (interaction score):")
        if (!is.na(cut_value)) {
          add_line("Cutpoint on interaction score: ", signif(cut_value, 4))
        } else {
          add_line("Cutpoint on interaction score: (could not extract numeric value)")
        }
        add_line(capture.output(print(tab_opt)))
        add_line("Optimal-cut log-rank p: ", signif(lr_opt_p, 4))
        add_line("----")
      } else {
        add_line("Optimal cutpoint: only one group after categorization, skipping.")
        add_line("----")
      }
    } else {
      add_line("Optimal cutpoint (interaction score): surv_cutpoint failed, skipping.")
      add_line("----")
    }
  } else {
    add_line("Optimal cutpoint (interaction score): not enough unique values, skipping.")
    add_line("----")
  }
  
  ## 4) Quartile-based KMs (on expr pathway score, with CNA score)
  expr_q <- tryCatch(
    {
      stats::quantile(d$expr, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
    },
    error = function(e) {
      NA
    }
  )
  
  if (any(is.na(expr_q))) {
    add_line("Quartiles: could not compute expression quartiles, skipping quartile-based KMs.")
    add_line("----")
  } else {
    d_q <- d %>%
      dplyr::mutate(
        expr_quartile = cut(
          expr,
          breaks = c(-Inf, expr_q[1], expr_q[2], expr_q[3], Inf),
          labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"),
          include.lowest = TRUE,
          right = TRUE
        ),
        cna_group_bin = factor(
          ifelse(cna > stats::median(cna, na.rm = TRUE), "High CNA", "Low CNA"),
          levels = c("Low CNA", "High CNA")
        )
      ) %>%
      dplyr::filter(!is.na(expr_quartile), !is.na(cna_group_bin))
    
    ## 4a) 8-group KM: expr_quartile x CNA group
    d8 <- d_q %>%
      dplyr::mutate(
        q_cna_group = interaction(expr_quartile, cna_group_bin, sep = " & ", drop = TRUE)
      )
    
    if (length(unique(na.omit(d8$q_cna_group))) >= 2) {
      fit8  <- survfit(Surv(OS_TIME, OS_EVENT) ~ q_cna_group, data = d8)
      lr8   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ q_cna_group, data = d8)
      lr8_p <- 1 - pchisq(lr8$chisq, df = length(lr8$n) - 1)
      
      tab8 <- d8 %>%
        dplyr::group_by(q_cna_group) %>%
        dplyr::summarise(
          n      = dplyr::n(),
          events = sum(OS_EVENT, na.rm = TRUE),
          .groups = "drop"
        )
      
      add_line("Quartile-based KM (8 groups: expr quartile x CNA group):")
      add_line(capture.output(print(tab8)))
      add_line("8-group log-rank p: ", signif(lr8_p, 4))
      add_line("----")
      
      p8 <- ggsurvplot(
        fit8,
        data                 = d8,
        risk.table           = TRUE,
        risk.table.y.text    = TRUE,
        risk.table.col       = NULL,
        risk.table.title     = "Number at risk",
        risk.table.height    = 0.30,
        risk.table.fontsize  = 3,
        pval                 = paste0("log-rank p = ", signif(lr8_p, 3)),
        legend.title         = "",
        legend.labs          = levels(d8$q_cna_group),
        title                = paste0(pathway_name, " - Expr quartiles x CNA (8 groups, OS, METABRIC pathways)"),
        xlab                 = "Time (months)",
        ylab                 = "Probability",
        ggtheme              = theme_bw(),
        tables.theme         = theme_bw()
      )
      
      p8$plot <- p8$plot +
        labs(subtitle = "Quartiles of pathway expression (Q1–Q4) crossed with CNA burden (Low/High)") +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
      if (!is.null(p8$table)) {
        if ("strata" %in% names(p8$table$data)) {
          old_levels <- levels(p8$table$data$strata)
          p8$table$data$strata <- factor(
            p8$table$data$strata,
            levels = old_levels,
            labels = levels(d8$q_cna_group)[seq_along(old_levels)]
          )
        }
        
        p8$table <- p8$table +
          theme(
            strip.background = element_blank(),
            strip.text       = element_blank(),
            axis.title.y     = element_blank(),
            axis.title.x     = element_blank()
          )
      }
      
      quart8_file <- file.path(
        out_dir,
        paste0("KM_OS_pathway_exprCNA_quartiles8_", clean_filename(pathway_name), model_suffix, ".pdf")
      )
      
      pdf(quart8_file, width = 9, height = 7)
      gridExtra::grid.arrange(p8$plot, p8$table, heights = c(0.7, 0.3))
      dev.off()
      
      message("Saved METABRIC OS PATHWAY quartile 8-group KM plot: ", quart8_file)
    } else {
      add_line("Quartile 8-group KM: fewer than 2 groups present, skipping.")
      add_line("----")
    }
    
    ## 4b) CNA=Low subset
    d0 <- d_q %>%
      dplyr::filter(cna_group_bin == "Low CNA")
    
    if (nrow(d0) >= 30 &&
        sum(d0$OS_EVENT, na.rm = TRUE) >= 5 &&
        length(unique(na.omit(d0$expr_quartile))) >= 2) {
      
      fit0  <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_quartile, data = d0)
      lr0   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_quartile, data = d0)
      lr0_p <- 1 - pchisq(lr0$chisq, df = length(lr0$n) - 1)
      
      tab0 <- d0 %>%
        dplyr::group_by(expr_quartile) %>%
        dplyr::summarise(
          n      = dplyr::n(),
          events = sum(OS_EVENT, na.rm = TRUE),
          .groups = "drop"
        )
      
      add_line("Quartile KM in Low CNA subset (4 groups: expr quartiles):")
      add_line(capture.output(print(tab0)))
      add_line("Low CNA quartile log-rank p: ", signif(lr0_p, 4))
      add_line("----")
      
      p0 <- ggsurvplot(
        fit0,
        data                 = d0,
        risk.table           = TRUE,
        risk.table.y.text    = TRUE,
        risk.table.col       = NULL,
        risk.table.title     = "Number at risk",
        risk.table.height    = 0.30,
        risk.table.fontsize  = 3,
        pval                 = paste0("log-rank p = ", signif(lr0_p, 3)),
        legend.title         = "",
        legend.labs          = levels(d0$expr_quartile),
        title                = paste0(pathway_name, " - Expr quartiles (Low CNA, OS, METABRIC pathways)"),
        xlab                 = "Time (months)",
        ylab                 = "Probability",
        ggtheme              = theme_bw(),
        tables.theme         = theme_bw()
      )
      
      p0$plot <- p0$plot +
        labs(subtitle = "Quartiles of pathway expression (Q1–Q4) within Low CNA burden") +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
      if (!is.null(p0$table)) {
        if ("strata" %in% names(p0$table$data)) {
          old_levels <- levels(p0$table$data$strata)
          p0$table$data$strata <- factor(
            p0$table$data$strata,
            levels = old_levels,
            labels = levels(d0$expr_quartile)[seq_along(old_levels)]
          )
        }
        
        p0$table <- p0$table +
          theme(
            strip.background = element_blank(),
            strip.text       = element_blank(),
            axis.title.y     = element_blank(),
            axis.title.x     = element_blank()
          )
      }
      
      quart0_file <- file.path(
        out_dir,
        paste0("KM_OS_pathway_exprCNA_quartiles_LowCNA_", clean_filename(pathway_name), model_suffix, ".pdf")
      )
      
      pdf(quart0_file, width = 9, height = 6)
      gridExtra::grid.arrange(p0$plot, p0$table, heights = c(0.7, 0.3))
      dev.off()
      
      message("Saved METABRIC OS PATHWAY quartile KM plot (Low CNA): ", quart0_file)
    } else {
      add_line("Quartile KM in Low CNA: insufficient data or <2 quartiles present, skipping.")
      add_line("----")
    }
    
    ## 4c) CNA=High subset
    d1 <- d_q %>%
      dplyr::filter(cna_group_bin == "High CNA")
    
    if (nrow(d1) >= 30 &&
        sum(d1$OS_EVENT, na.rm = TRUE) >= 5 &&
        length(unique(na.omit(d1$expr_quartile))) >= 2) {
      
      fit1  <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_quartile, data = d1)
      lr1   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_quartile, data = d1)
      lr1_p <- 1 - pchisq(lr1$chisq, df = length(lr1$n) - 1)
      
      tab1 <- d1 %>%
        dplyr::group_by(expr_quartile) %>%
        dplyr::summarise(
          n      = dplyr::n(),
          events = sum(OS_EVENT, na.rm = TRUE),
          .groups = "drop"
        )
      
      add_line("Quartile KM in High CNA subset (4 groups: expr quartiles):")
      add_line(capture.output(print(tab1)))
      add_line("High CNA quartile log-rank p: ", signif(lr1_p, 4))
      add_line("----")
      
      p1 <- ggsurvplot(
        fit1,
        data                 = d1,
        risk.table           = TRUE,
        risk.table.y.text    = TRUE,
        risk.table.col       = NULL,
        risk.table.title     = "Number at risk",
        risk.table.height    = 0.30,
        risk.table.fontsize  = 3,
        pval                 = paste0("log-rank p = ", signif(lr1_p, 3)),
        legend.title         = "",
        legend.labs          = levels(d1$expr_quartile),
        title                = paste0(pathway_name, " - Expr quartiles (High CNA, OS, METABRIC pathways)"),
        xlab                 = "Time (months)",
        ylab                 = "Probability",
        ggtheme              = theme_bw(),
        tables.theme         = theme_bw()
      )
      
      p1$plot <- p1$plot +
        labs(subtitle = "Quartiles of pathway expression (Q1–Q4) within High CNA burden") +
        theme(
          plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)
        )
      
      if (!is.null(p1$table)) {
        if ("strata" %in% names(p1$table$data)) {
          old_levels <- levels(p1$table$data$strata)
          p1$table$data$strata <- factor(
            p1$table$data$strata,
            levels = old_levels,
            labels = levels(d1$expr_quartile)[seq_along(old_levels)]
          )
        }
        
        p1$table <- p1$table +
          theme(
            strip.background = element_blank(),
            strip.text       = element_blank(),
            axis.title.y     = element_blank(),
            axis.title.x     = element_blank()
          )
      }
      
      quart1_file <- file.path(
        out_dir,
        paste0("KM_OS_pathway_exprCNA_quartiles_HighCNA_", clean_filename(pathway_name), model_suffix, ".pdf")
      )
      
      pdf(quart1_file, width = 9, height = 6)
      gridExtra::grid.arrange(p1$plot, p1$table, heights = c(0.7, 0.3))
      dev.off()
      
      message("Saved METABRIC OS PATHWAY quartile KM plot (High CNA): ", quart1_file)
    } else {
      add_line("Quartile KM in High CNA: insufficient data or <2 quartiles present, skipping.")
      add_line("----")
    }
  }
  
  ## 5) Refit Cox model on pathway scores: Surv ~ expr + cna + expr:cna
  cox_int <- tryCatch(
    {
      coxph(Surv(OS_TIME, OS_EVENT) ~ expr + cna + expr:cna, data = d)
    },
    error = function(e) {
      NULL
    }
  )
  
  if (!is.null(cox_int)) {
    cox_sum <- summary(cox_int)
    
    if ("expr:cna" %in% rownames(cox_sum$coefficients)) {
      beta_int <- cox_sum$coefficients["expr:cna", "coef"]
      HR_int   <- cox_sum$coefficients["expr:cna", "exp(coef)"]
      p_int    <- cox_sum$coefficients["expr:cna", "Pr(>|z|)"]
      se_int   <- cox_sum$coefficients["expr:cna", "se(coef)"]
      LCL_int  <- exp(beta_int - 1.96 * se_int)
      UCL_int  <- exp(beta_int + 1.96 * se_int)
      
      add_line("Refitted Cox (Surv ~ expr + cna + expr:cna) on PATHWAY scores:")
      add_line(sprintf("  expr:cna: HR=%.3f, 95%% CI [%.3f, %.3f], p=%.4g",
                       HR_int, LCL_int, UCL_int, p_int))
    } else {
      add_line("Refitted Cox: expr:cna coefficient not found in summary.")
    }
    
    zph_res <- tryCatch(
      {
        cox.zph(cox_int)
      },
      error = function(e) {
        NULL
      }
    )
    
    if (!is.null(zph_res)) {
      zph_tab <- zph_res$table
      
      if ("expr:cna" %in% rownames(zph_tab)) {
        p_ph_int <- zph_tab["expr:cna", "p"]
        add_line(sprintf("PH test (cox.zph) for expr:cna: p=%.4g", p_ph_int))
      } else {
        add_line("PH test: expr:cna row not found in zph table.")
      }
      
      if ("GLOBAL" %in% rownames(zph_tab)) {
        p_ph_global <- zph_tab["GLOBAL", "p"]
        add_line(sprintf("PH test (cox.zph) GLOBAL p=%.4g", p_ph_global))
      }
    } else {
      add_line("PH test (cox.zph): failed or not computed.")
    }
    
    add_line("----")
  } else {
    add_line("Refitted Cox (expr + cna + expr:cna): model failed, skipping.")
    add_line("----")
  }
  
  ## 6) Spline Cox model on expr score with interaction by cna
  spline_pdf_file <- file.path(
    out_dir,
    paste0("Cox_spline_OS_pathway_exprCNA_", clean_filename(pathway_name), model_suffix, ".pdf")
  )
  
  spline_fit <- tryCatch(
    {
      coxph(Surv(OS_TIME, OS_EVENT) ~ splines::ns(expr, df = 3) * cna, data = d)
    },
    error = function(e) {
      NULL
    }
  )
  
  if (!is.null(spline_fit)) {
    expr_q   <- stats::quantile(d$expr, probs = c(0.05, 0.95), na.rm = TRUE)
    expr_seq <- seq(expr_q[1], expr_q[2], length.out = 50)
    
    ## Use median expr in high CNA as reference; if no high CNA, use all
    med_expr_cna_high <- stats::median(d$expr[d$cna > stats::median(d$cna, na.rm = TRUE)], na.rm = TRUE)
    if (is.na(med_expr_cna_high)) {
      med_expr_cna_high <- stats::median(d$expr, na.rm = TRUE)
    }
    
    newdata_seq <- data.frame(expr = expr_seq, cna = 1)
    newdata_ref <- data.frame(expr = med_expr_cna_high, cna = 1)
    
    lp_seq <- as.numeric(predict(spline_fit, newdata = newdata_seq, type = "lp"))
    lp_ref <- as.numeric(predict(spline_fit, newdata = newdata_ref, type = "lp"))[1]
    
    HR_seq <- exp(lp_seq - lp_ref)
    
    df_plot <- data.frame(
      expr = expr_seq,
      HR   = HR_seq
    )
    
    g <- ggplot(df_plot, aes(x = expr, y = HR)) +
      geom_line() +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(
        title = paste0(pathway_name, " - Spline Cox for pathway expr (cna=1, OS)"),
        subtitle = "Model: Surv ~ ns(expr, 3) * cna (relative hazard vs median expr at high CNA)",
        x = "Pathway expression score",
        y = "Relative hazard (cna = 1)"
      ) +
      theme_bw() +
      theme(
        plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    pdf(spline_pdf_file, width = 7, height = 5)
    print(g)
    dev.off()
    
    add_line("Spline Cox plot saved to: ", spline_pdf_file)
  } else {
    add_line("Spline Cox model (ns(expr,3)*cna): fitting failed, no spline PDF produced.")
  }
  
  add_line("----")
  
  txt_file <- file.path(
    out_dir,
    paste0("KM_OS_pathway_exprCNA_interaction_details_", clean_filename(pathway_name), model_suffix, ".txt")
  )
  
  writeLines(summary_lines, con = txt_file)
  message("Extra PATHWAY expr:CNA diagnostics saved for ", pathway_name, " to: ", txt_file)
  
  invisible(txt_file)
}

## --------------------------------------------------
## 10) Read OS Cox Excel (PATHWAY-level) and loop over sheets/pathways
## --------------------------------------------------

cox_file <- file.path(
  metabric_results_dir,
  "METABRIC_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx"
)

all_sheets <- readxl::excel_sheets(cox_file)

## Exclude sheets whose names contain:
## "_filter", "_filter1", "_filter2", "LR_"
valid_sheets <- all_sheets[
  !grepl("_filter",  all_sheets) &
    !grepl("_filter1", all_sheets) &
    !grepl("_filter2", all_sheets) &
    !grepl("LR_",      all_sheets) &
    all_sheets != "Pathway_scores"
]

message("PATHWAY sheets to be processed (OS all): ", paste(valid_sheets, collapse = ", "))

for (sh in valid_sheets) {
  message("\n==============================")
  message("Processing PATHWAY sheet (OS all): ", sh)
  
  df <- readxl::read_xlsx(cox_file, sheet = sh)
  
  ## Feature column: can be "pathway" or (fallback) "gene"
  feat_col <- NULL
  if ("pathway" %in% names(df)) {
    feat_col <- "pathway"
  } else if ("gene" %in% names(df)) {
    feat_col <- "gene"
  } else {
    message("Sheet ", sh, " does not contain 'pathway' or 'gene' column. Skipping.")
    next
  }
  
  ## Create a folder for this sheet
  sheet_dir <- file.path(km_dir, sh)
  if (!dir.exists(sheet_dir)) {
    dir.create(sheet_dir, recursive = TRUE)
  }
  
  ## Does the sheet name contain "E", "C", "EC"?
  tokens      <- unlist(strsplit(sh, "_"))
  cov_tokens  <- tokens[tokens %in% c("E", "C", "EC")]
  has_EC_name <- "EC" %in% cov_tokens
  
  model_suffix <- if (length(cov_tokens) > 0) {
    paste0("_", paste(cov_tokens, collapse = "_"))
  } else {
    ""
  }
  
  ## CI columns
  has_ci_expr    <- "CI_lo_expr"     %in% names(df)
  has_ci_cna     <- "CI_lo_cna"      %in% names(df)
  has_ci_exprcna <- "CI_lo_expr_cna" %in% names(df)
  
  if (!has_ci_expr && !has_ci_cna && !has_ci_exprcna) {
    message("Sheet ", sh, " has no CI_lo columns. Skipping.")
    next
  }
  
  for (i in seq_len(nrow(df))) {
    pathway_name <- df[[feat_col]][i]
    if (is.na(pathway_name) || pathway_name == "") next
    
    HR_expr <- if ("HR_expr" %in% names(df)) as.numeric(df$HR_expr[i]) else NA_real_
    HR_cna  <- if ("HR_cna"  %in% names(df)) as.numeric(df$HR_cna[i])  else NA_real_
    HR_int  <- if ("HR_expr_cna" %in% names(df)) as.numeric(df$HR_expr_cna[i]) else NA_real_
    
    cond_expr_raw    <- has_ci_expr    && !is.na(df$CI_lo_expr[i])     && df$CI_lo_expr[i]     > 1.00
    cond_cna_raw     <- has_ci_cna     && !is.na(df$CI_lo_cna[i])      && df$CI_lo_cna[i]      > 1.00
    cond_exprcna_raw <- has_ci_exprcna && !is.na(df$CI_lo_expr_cna[i]) && df$CI_lo_expr_cna[i] > 1.00
    
    cond_expr    <- cond_expr_raw    && !is.na(HR_expr) && HR_expr > 1.00
    cond_cna     <- cond_cna_raw     && !is.na(HR_cna)  && HR_cna  > 1.00
    cond_exprcna <- cond_exprcna_raw && !is.na(HR_int)  && HR_int  > 1.00
    
    if (!cond_expr && !cond_cna && !cond_exprcna) next
    
    message("  Pathway: ", pathway_name,
            " | expr? ", cond_expr,
            " | cna? ", cond_cna,
            " | exprCNA? ", cond_exprcna,
            " | EC_in_name? ", has_EC_name)
    
    ## Expression KM (median + optimal)
    if (cond_expr && all(c("HR_expr", "CI_lo_expr", "CI_hi_expr", "p_expr") %in% names(df))) {
      LCL_expr  <- as.numeric(df$CI_lo_expr[i])
      UCL_expr  <- as.numeric(df$CI_hi_expr[i])
      PVAL_expr <- as.numeric(df$p_expr[i])
      
      plot_km_pathway_expr_median(
        pathway_name  = pathway_name,
        HR            = HR_expr,
        LCL           = LCL_expr,
        UCL           = UCL_expr,
        PVAL          = PVAL_expr,
        out_dir       = sheet_dir,
        model_suffix  = model_suffix,
        require_sig_p = TRUE
      )
      
      plot_km_pathway_expr_optimal(
        pathway_name  = pathway_name,
        HR            = HR_expr,
        LCL           = LCL_expr,
        UCL           = UCL_expr,
        PVAL          = PVAL_expr,
        out_dir       = sheet_dir,
        model_suffix  = model_suffix,
        require_sig_p = TRUE
      )
    } else if (cond_expr) {
      message("    (expr) HR_expr/CI/p columns missing for ", pathway_name, " in sheet ", sh, " – skipping expr plots.")
    }
    
    ## CNA KM
    if (cond_cna && all(c("HR_cna", "CI_lo_cna", "CI_hi_cna", "p_cna") %in% names(df))) {
      LCL_cna  <- as.numeric(df$CI_lo_cna[i])
      UCL_cna  <- as.numeric(df$CI_hi_cna[i])
      PVAL_cna <- as.numeric(df$p_cna[i])
      
      plot_km_pathway_cna(
        pathway_name  = pathway_name,
        HR            = HR_cna,
        LCL           = LCL_cna,
        UCL           = UCL_cna,
        PVAL          = PVAL_cna,
        out_dir       = sheet_dir,
        model_suffix  = model_suffix,
        require_sig_p = TRUE
      )
    } else if (cond_cna) {
      message("    (cna) HR_cna/CI/p columns missing for ", pathway_name, " in sheet ", sh, " – skipping CNA plots.")
    }
    
    ## Expr:CNA interaction KM + diagnostics
    if (cond_exprcna && all(c("HR_expr_cna", "CI_lo_expr_cna", "CI_hi_expr_cna", "p_expr_cna") %in% names(df))) {
      LCL_int  <- as.numeric(df$CI_lo_expr_cna[i])
      UCL_int  <- as.numeric(df$CI_hi_expr_cna[i])
      PVAL_int <- as.numeric(df$p_expr_cna[i])
      
      plot_km_pathway_expr_cna_interaction(
        pathway_name = pathway_name,
        HR           = HR_int,
        LCL          = LCL_int,
        UCL          = UCL_int,
        PVAL         = PVAL_int,
        out_dir      = sheet_dir,
        model_suffix = model_suffix,
        four_groups  = has_EC_name
      )
      
      extra_pathway_expr_cna_diagnostics(
        pathway_name = pathway_name,
        out_dir      = sheet_dir,
        model_suffix = model_suffix,
        four_groups  = has_EC_name
      )
    } else if (cond_exprcna) {
      message("    (exprCNA) HR_expr_cna/CI/p columns missing for ", pathway_name, " in sheet ", sh, " – skipping interaction plots.")
    }
  }
}

message("\nAll METABRIC OS PATHWAY KM plots and interaction diagnostics have been saved under: ", km_dir)