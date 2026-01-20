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

## --------------------------------------------------
## 1) Directories
## --------------------------------------------------

tcga_dir         <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
pathfindR_dir    <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"
gdsc_dir         <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
tcga_results_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/TCGA-BRCA"

if (!dir.exists(tcga_results_dir)) {
  dir.create(tcga_results_dir, recursive = TRUE)
}

## Output folder for KM plots (OS, all TCGA-BRCA, PATHWAY)
km_dir <- file.path(tcga_results_dir, "TCGA_KM_Pathway_OS")
if (!dir.exists(km_dir)) {
  dir.create(km_dir, recursive = TRUE)
}

## --------------------------------------------------
## 2) Read TCGA-BRCA expression, CNA, survival
## --------------------------------------------------

tcga_expr_file     <- file.path(tcga_dir, "TCGA-BRCA_exprs.z.tsv")
tcga_cna_file      <- file.path(tcga_dir, "TCGA-BRCA_CNA.tsv")
tcga_survival_file <- file.path(tcga_dir, "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv")

tcga_expr     <- read_tsv(tcga_expr_file,     show_col_types = FALSE)
tcga_cna      <- read_tsv(tcga_cna_file,      show_col_types = FALSE)
tcga_survival <- read_tsv(tcga_survival_file, show_col_types = FALSE)

## Make CNA binary (0 vs !=0) → 0 / 1
tcga_cna_binary <- tcga_cna
tcga_cna_binary[, -1] <- ifelse(tcga_cna_binary[, -1] != 0, 1, 0)

## Entrez ID → Gene Symbol for TCGA expr & CNA
tcga_expr$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_expr$ENTREZID),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

tcga_cna_binary$gene <- mapIds(
  org.Hs.eg.db,
  keys      = as.character(tcga_cna_binary$gene),
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

colnames(tcga_expr)[colnames(tcga_expr) == "ENTREZID"]         <- "Gene_Symbol"
colnames(tcga_cna_binary)[colnames(tcga_cna_binary) == "gene"] <- "Gene_Symbol"

## --------------------------------------------------
## 3) TCGA-BRCA OS + covariates (Age, Stage_simple, PAM50) – ALL patients
##    (same time scale as gene-level KM: days)
## --------------------------------------------------

tcga_surv_small <- tcga_survival %>%
  dplyr::transmute(
    sample,
    OS_TIME_d = as.numeric(OS.time),   ## days
    OS_EVENT  = as.numeric(OS),        ## 0 / 1
    Age        = as.numeric(Age),
    Stage      = factor(Stage),
    PAM50      = factor(BRCA_Subtype_PAM50)
  ) %>%
  dplyr::filter(
    !is.na(OS_TIME_d),
    !is.na(OS_EVENT)
  )

tcga_surv_small <- tcga_surv_small %>%
  dplyr::mutate(
    PAM50 = {
      f <- factor(PAM50)
      if ("LumA" %in% levels(f)) stats::relevel(f, ref = "LumA") else f
    },
    Stage = {
      f <- factor(Stage)
      if ("Stage I" %in% levels(f)) stats::relevel(f, ref = "Stage I") else f
    }
  )

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

message("Number of TCGA-BRCA samples with OS data (ALL patients): ", nrow(tcga_surv_small))

## --------------------------------------------------
## 4) GDSC data (for Fisher OR direction) – same as pathway Cox
## --------------------------------------------------

gdsc_expr_file     <- file.path(gdsc_dir, "GDSC_exprs.Paclitaxel.eb_with.TCGA_exprs.Paclitaxel.tsv")
gdsc_cna_file      <- file.path(gdsc_dir, "GDSC_CNA.Paclitaxel.tsv")
gdsc_response_file <- file.path(gdsc_dir, "GDSC_response.Paclitaxel.tsv")

gdsc_expr     <- read_tsv(gdsc_expr_file,     show_col_types = FALSE)
gdsc_cna      <- read_tsv(gdsc_cna_file,      show_col_types = FALSE)
gdsc_response <- read_tsv(gdsc_response_file, show_col_types = FALSE)

gdsc_cna[, -1] <- ifelse(gdsc_cna[, -1] != 0, 1, 0)

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
colnames(gdsc_cna)[colnames(gdsc_cna)   == "gene_id"]  <- "Gene_Symbol"

resp <- gdsc_response %>%
  mutate(status = response)

gdsc_samples <- as.character(unique(resp$sample_name))

expr_gdsc <- gdsc_expr %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

expr_mat <- as.matrix(expr_gdsc[, -1])
rownames(expr_mat) <- expr_gdsc$Gene_Symbol
expr_z   <- t(scale(t(expr_mat)))
expr_z[is.na(expr_z)] <- 0
expr_z_df <- as_tibble(expr_z, rownames = "Gene_Symbol")

cna_gdsc <- gdsc_cna %>%
  dplyr::select(Gene_Symbol, all_of(gdsc_samples))

cna_mat <- as.matrix(cna_gdsc[, -1])
rownames(cna_mat) <- cna_gdsc$Gene_Symbol
cna_bin    <- ifelse(cna_mat != 0, 1, 0)
cna_bin_df <- as_tibble(cna_bin, rownames = "Gene_Symbol")

expr_long_gdsc <- expr_z_df %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr_z"
  ) %>%
  dplyr::rename(sample_name = "sample") %>%
  inner_join(
    resp %>%
      dplyr::select(sample_name, status) %>%
      mutate(sample_name = as.character(sample_name)),
    by = "sample_name"
  )

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

## --------------------------------------------------
## 5) pathfindR enrichment & pathway–gene mapping – same selection as Cox
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

path_genes_long <- enrichment_df %>%
  mutate(
    up_genes   = Up_regulated,
    down_genes = Down_regulated
  ) %>%
  unite("genes_str", Up_regulated, Down_regulated, sep = ",") %>%
  dplyr::rename(pathway = "Term_Description") %>%
  separate_rows(up_genes,   sep = "[,;]") %>%
  separate_rows(down_genes, sep = "[,;]") %>%
  separate_rows(genes_str,  sep = "[,;]") %>%
  mutate(
    Gene_Symbol = str_trim(genes_str),
    direction   = case_when(
      genes_str %in% up_genes   ~ "Up",
      genes_str %in% down_genes ~ "Down",
      TRUE                      ~ NA_character_
    )
  ) %>%
  filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
  dplyr::select(pathway, Gene_Symbol, direction) %>%
  distinct() %>%
  na.omit()

## --------------------------------------------------
## 6) TCGA expression & CNA in long format (for pathway scores)
## --------------------------------------------------

expr_long_tcga <- tcga_expr %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "expr"
  ) %>%
  dplyr::mutate(
    sample = substr(sample_barcode, 1, 15)
  )

cna_long_tcga <- tcga_cna_binary %>%
  tidyr::pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample_barcode",
    values_to = "cna"
  ) %>%
  dplyr::mutate(
    sample = substr(sample_barcode, 1, 15)
  )

## --------------------------------------------------
## 7) Pathway expression scores per sample (TCGA-BRCA)
## --------------------------------------------------

expr_path_tcga <- expr_long_tcga %>%
  inner_join(path_genes_long, by = "Gene_Symbol", relationship = "many-to-many") %>%
  mutate(
    dir_factor   = case_when(
      direction == "Up"   ~  1,
      direction == "Down" ~ -1,
      TRUE                ~ NA_real_
    ),
    expr_contrib = abs(expr) * dir_factor
  )

path_expr_score_tcga_sample <- expr_path_tcga %>%
  group_by(sample, pathway) %>%
  summarise(
    expr    = mean(expr_contrib, na.rm = TRUE),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  filter(n_genes > 1)

## --------------------------------------------------
## 8) Pathway CNA scores per sample (aligned with GDSC Fisher OR)
## --------------------------------------------------

cna_path_long_tcga <- cna_long_tcga %>%
  inner_join(path_genes_long, by = "Gene_Symbol", relationship = "many-to-many") %>%
  inner_join(fisher_by_gene,  by = "Gene_Symbol", relationship = "many-to-many") %>%
  mutate(
    dir_cna = case_when(
      !is.na(OR) & OR > 1 ~  1,
      !is.na(OR) & OR < 1 ~ -1,
      TRUE                ~ NA_real_
    ),
    cna_contrib = case_when(
      is.na(dir_cna) ~ NA_real_,
      dir_cna ==  1  ~ cna,
      dir_cna == -1  ~ 1 - cna
    )
  )

path_cna_score_tcga_sample <- cna_path_long_tcga %>%
  group_by(sample, pathway) %>%
  summarise(
    cna = mean(cna_contrib, na.rm = TRUE),
    .groups   = "drop"
  )

## --------------------------------------------------
## 9) Helper: get per-pathway survival + omics data (OS, ALL TCGA-BRCA)
## --------------------------------------------------

get_pathway_surv_data_tcga_OS <- function(pathway_name) {
  path_expr_score_tcga_sample %>%
    dplyr::filter(pathway == pathway_name) %>%
    dplyr::select(sample, pathway, expr) %>%
    dplyr::inner_join(
      path_cna_score_tcga_sample %>%
        dplyr::filter(pathway == pathway_name) %>%
        dplyr::select(sample, pathway, cna),
      by = c("sample", "pathway")
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

## --------------------------------------------------
## 10) KM plotting helpers (PATHWAY, TCGA-BRCA, OS)
## --------------------------------------------------

## 10a) CNA-based KM (2 groups)
plot_km_cna_tcga_OS_pathway <- function(pathway_name,
                                         HR, LCL, UCL, PVAL,
                                         out_dir,
                                         model_suffix,
                                         require_sig_p = TRUE) {
  d <- get_pathway_surv_data_tcga_OS(pathway_name) %>%
    dplyr::filter(!is.na(cna))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("PATHWAY CNA KM (OS TCGA-BRCA): insufficient data (n < 30 or < 5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  d <- d %>%
    dplyr::mutate(
      cna_group = factor(
        ifelse(cna == 0, "CNA = 0", "CNA = 1"),
        levels = c("CNA = 0", "CNA = 1")
      )
    )
  
  if (length(unique(d$cna_group)) < 2) {
    message("PATHWAY CNA KM (OS TCGA-BRCA): only one CNA group present, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ cna_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  if (require_sig_p && (is.na(lr_p) || lr_p >= 0.05)) {
    message("PATHWAY CNA KM (OS TCGA-BRCA): log-rank p >= 0.05, skipping pathway: ", pathway_name)
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
    title                = paste0(pathway_name, " - CNA (OS, TCGA-BRCA)"),
    xlab                 = "Time (days)",
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
  
  file_tag <- paste0("KM_OS_CNA_", make.names(pathway_name), model_suffix, ".pdf")
  out_file <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 8, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved TCGA-BRCA OS PATHWAY CNA KM plot: ", out_file)
  invisible(out_file)
}

## 10b) Expression – median cut (2 groups)
plot_km_expr_median_tcga_OS_pathway <- function(pathway_name,
                                                 HR, LCL, UCL, PVAL,
                                                 out_dir,
                                                 model_suffix,
                                                 require_sig_p = TRUE) {
  d <- get_pathway_surv_data_tcga_OS(pathway_name) %>%
    dplyr::filter(!is.na(expr))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("PATHWAY Expr-median KM (OS TCGA-BRCA): insufficient data (n < 30 or < 5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  med <- stats::median(d$expr, na.rm = TRUE)
  
  d <- d %>%
    dplyr::mutate(
      expr_group = factor(
        ifelse(expr > med, "High", "Low"),
        levels = c("Low", "High")
      )
    )
  
  if (length(unique(d$expr_group)) < 2) {
    message("PATHWAY Expr-median KM (OS TCGA-BRCA): only one expression group present, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  if (require_sig_p && (is.na(lr_p) || lr_p >= 0.05)) {
    message("PATHWAY Expr-median KM (OS TCGA-BRCA): log-rank p >= 0.05, skipping pathway: ", pathway_name)
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
    title                = paste0(pathway_name, " - Expression median (OS, TCGA-BRCA)"),
    xlab                 = "Time (days)",
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
  
  file_tag <- paste0("KM_OS_expr_median_", make.names(pathway_name), model_suffix, ".pdf")
  out_file <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 8, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved TCGA-BRCA OS PATHWAY expression-median KM plot: ", out_file)
  invisible(out_file)
}

## 10c) Expression – optimal cut (maxstat, 2 groups)
plot_km_expr_optimal_tcga_OS_pathway <- function(pathway_name,
                                                  HR, LCL, UCL, PVAL,
                                                  out_dir,
                                                  model_suffix,
                                                  require_sig_p = TRUE) {
  d <- get_pathway_surv_data_tcga_OS(pathway_name) %>%
    dplyr::filter(!is.na(expr))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("PATHWAY Expr-optimal KM (OS TCGA-BRCA): insufficient data (n < 30 or < 5 events), skipping pathway: ", pathway_name)
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
      message("PATHWAY Expr-optimal KM (OS TCGA-BRCA): surv_cutpoint failed for ", pathway_name, " – skipping. Error: ", e$message)
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
        "low"  = "Low",
        "high" = "High",
        .default = NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(expr_group_opt)) %>%
    dplyr::mutate(
      expr_group_opt = factor(expr_group_opt, levels = c("Low", "High"))
    )
  
  if (length(unique(d$expr_group_opt)) < 2) {
    message("PATHWAY Expr-optimal KM (OS TCGA-BRCA): only one optimal expression group present, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ expr_group_opt, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  if (require_sig_p && (is.na(lr_p) || lr_p >= 0.05)) {
    message("PATHWAY Expr-optimal KM (OS TCGA-BRCA): log-rank p >= 0.05, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  fit <- survfit(Surv(OS_TIME, OS_EVENT) ~ expr_group_opt, data = d)
  
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
      title                = paste0(pathway_name, " - Expression optimal (OS, TCGA-BRCA)"),
      xlab                 = "Time (days)",
      ylab                 = "OS probability",
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
  
  file_tag <- paste0("KM_OS_expr_optimal_", make.names(pathway_name), model_suffix, ".pdf")
  out_file <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 8, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved TCGA-BRCA OS PATHWAY expression-optimal KM plot: ", out_file)
  invisible(out_file)
}

## 10d) expr:cna interaction KM (PATHWAY)
plot_km_expr_cna_interaction_tcga_OS_pathway <- function(pathway_name,
                                                          HR, LCL, UCL, PVAL,
                                                          out_dir,
                                                          model_suffix,
                                                          four_groups = TRUE) {
  d <- get_pathway_surv_data_tcga_OS(pathway_name) %>%
    dplyr::filter(!is.na(expr), !is.na(cna))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("PATHWAY Expr:CNA KM (OS TCGA-BRCA): insufficient data (n < 30 or < 5 events), skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  if (four_groups) {
    med <- stats::median(d$expr, na.rm = TRUE)
    
    d <- d %>%
      dplyr::mutate(
        expr_group = factor(
          ifelse(expr > med, "High exp", "Low exp"),
          levels = c("Low exp", "High exp")
        ),
        cna_group = factor(
          ifelse(cna == 0, "CNA = 0", "CNA = 1"),
          levels = c("CNA = 0", "CNA = 1")
        ),
        int_group = factor(
          paste(expr_group, cna_group, sep = " & "),
          levels = c("Low exp & CNA = 0",
                     "Low exp & CNA = 1",
                     "High exp & CNA = 0",
                     "High exp & CNA = 1")
        )
      )
  } else {
    d <- d %>%
      dplyr::mutate(
        int_score = expr * cna
      )
    
    if (length(unique(na.omit(d$int_score))) < 2) {
      message("PATHWAY Expr:CNA KM (OS TCGA-BRCA): interaction score has < 2 unique values, skipping pathway: ", pathway_name)
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
  
  if (length(na.omit(unique(d$int_group))) < 2) {
    message("PATHWAY Expr:CNA KM (OS TCGA-BRCA): fewer than 2 interaction groups present, skipping pathway: ", pathway_name)
    return(NULL)
  }
  
  fit  <- survfit(Surv(OS_TIME, OS_EVENT) ~ int_group, data = d)
  lr   <- survdiff(Surv(OS_TIME, OS_EVENT) ~ int_group, data = d)
  lr_p <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)
  
  ## For interaction: always plot, regardless of p
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
    risk.table.height    = if (four_groups) 0.30 else 0.25,
    risk.table.fontsize  = 3,
    pval                 = paste0("log-rank p = ", signif(lr_p, 3)),
    legend.title         = "",
    legend.labs          = levels(d$int_group),
    title                = paste0(pathway_name, " - Expression x CNA (OS, TCGA-BRCA)"),
    xlab                 = "Time (days)",
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
  
  if (!is.null(p$table)) {
    if ("strata" %in% names(p$table$data)) {
      old_levels <- levels(p$table$data$strata)
      p$table$data$strata <- factor(
        p$table$data$strata,
        levels = old_levels,
        labels = levels(d$int_group)[seq_along(old_levels)]
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
  
  file_tag <- paste0("KM_OS_exprCNA_interaction_", make.names(pathway_name), model_suffix, ".pdf")
  out_file <- file.path(out_dir, file_tag)
  
  pdf(out_file, width = 8, height = 6)
  gridExtra::grid.arrange(p$plot, p$table, heights = c(0.7, 0.3))
  dev.off()
  
  message("Saved TCGA-BRCA OS PATHWAY expr:CNA interaction KM plot: ", out_file)
  invisible(out_file)
}

## 10e) Extra diagnostics for expr:CNA interaction pathways (same as gene version, but pathway_name)
extra_expr_cna_diagnostics_tcga_OS_pathway <- function(pathway_name,
                                                        out_dir,
                                                        model_suffix,
                                                        four_groups = TRUE) {
  d <- get_pathway_surv_data_tcga_OS(pathway_name) %>%
    dplyr::filter(!is.na(expr), !is.na(cna))
  
  if (nrow(d) < 30 || sum(d$OS_EVENT, na.rm = TRUE) < 5) {
    message("Extra expr:CNA diagnostics (TCGA-BRCA OS PATHWAY): insufficient data (n < 30 or < 5 events), skipping pathway: ", pathway_name)
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
  
  ## 1) Four-group summary (Low/High x CNA=0/1) if requested
  if (four_groups) {
    med_expr <- stats::median(d$expr, na.rm = TRUE)
    
    d4 <- d %>%
      dplyr::mutate(
        expr_group = factor(
          ifelse(expr > med_expr, "High exp", "Low exp"),
          levels = c("Low exp", "High exp")
        ),
        cna_group = factor(
          ifelse(cna == 0, "CNA = 0", "CNA = 1"),
          levels = c("CNA = 0", "CNA = 1")
        ),
        int_group = factor(
          paste(expr_group, cna_group, sep = " & "),
          levels = c(
            "Low exp & CNA = 0",
            "Low exp & CNA = 1",
            "High exp & CNA = 0",
            "High exp & CNA = 1"
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
      
      add_line("Four-group KM (Low/High x CNA=0/1):")
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
    add_line("Interaction score has < 2 unique values, skipping 2-group KM.")
    add_line("----")
  }
  
  ## 3) Optimal cutpoint for interaction score (if feasible)
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
  
  ## 4) Quartile-based KMs (expression quartiles, with CNA)
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
        cna_group = factor(
          ifelse(cna == 0, "CNA = 0", "CNA = 1"),
          levels = c("CNA = 0", "CNA = 1")
        )
      ) %>%
      dplyr::filter(!is.na(expr_quartile), !is.na(cna_group))
    
    ## 4a) 8-group KM: expr_quartile x CNA
    d8 <- d_q %>%
      dplyr::mutate(
        q_cna_group = interaction(expr_quartile, cna_group, sep = " & ", drop = TRUE)
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
      
      add_line("Quartile-based KM (8 groups: expr quartile x CNA):")
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
        title                = paste0(pathway_name, " - Expr quartiles x CNA (8 groups, OS, TCGA-BRCA)"),
        xlab                 = "Time (days)",
        ylab                 = "OS probability",
        ggtheme              = theme_bw(),
        tables.theme         = theme_bw()
      )
      
      p8$plot <- p8$plot +
        labs(subtitle = "Quartiles of expression (Q1–Q4) crossed with CNA (0/1)") +
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
        paste0("KM_OS_exprCNA_quartiles8_", make.names(pathway_name), model_suffix, ".pdf")
      )
      
      pdf(quart8_file, width = 9, height = 7)
      gridExtra::grid.arrange(p8$plot, p8$table, heights = c(0.7, 0.3))
      dev.off()
      
      message("Saved TCGA-BRCA OS PATHWAY quartile 8-group KM plot: ", quart8_file)
    } else {
      add_line("Quartile 8-group KM: fewer than 2 groups present, skipping.")
      add_line("----")
    }
    
    ## 4b) 4-group KM: expr quartiles in CNA=0 subset
    d0 <- d_q %>%
      dplyr::filter(cna_group == "CNA = 0")
    
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
      
      add_line("Quartile KM in CNA=0 subset (4 groups: expr quartiles):")
      add_line(capture.output(print(tab0)))
      add_line("CNA=0 quartile log-rank p: ", signif(lr0_p, 4))
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
        title                = paste0(pathway_name, " - Expr quartiles (CNA = 0, OS, TCGA-BRCA)"),
        xlab                 = "Time (days)",
        ylab                 = "OS probability",
        ggtheme              = theme_bw(),
        tables.theme         = theme_bw()
      )
      
      p0$plot <- p0$plot +
        labs(subtitle = "Quartiles of expression (Q1–Q4) within CNA = 0") +
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
        paste0("KM_OS_exprCNA_quartiles_CNA0_", make.names(pathway_name), model_suffix, ".pdf")
      )
      
      pdf(quart0_file, width = 8, height = 6)
      gridExtra::grid.arrange(p0$plot, p0$table, heights = c(0.7, 0.3))
      dev.off()
      
      message("Saved TCGA-BRCA OS PATHWAY quartile KM plot (CNA=0): ", quart0_file)
    } else {
      add_line("Quartile KM in CNA=0: insufficient data or < 2 quartiles present, skipping.")
      add_line("----")
    }
    
    ## 4c) 4-group KM: expr quartiles in CNA=1 subset
    d1 <- d_q %>%
      dplyr::filter(cna_group == "CNA = 1")
    
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
      
      add_line("Quartile KM in CNA=1 subset (4 groups: expr quartiles):")
      add_line(capture.output(print(tab1)))
      add_line("CNA=1 quartile log-rank p: ", signif(lr1_p, 4))
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
        title                = paste0(pathway_name, " - Expr quartiles (CNA = 1, OS, TCGA-BRCA)"),
        xlab                 = "Time (days)",
        ylab                 = "OS probability",
        ggtheme              = theme_bw(),
        tables.theme         = theme_bw()
      )
      
      p1$plot <- p1$plot +
        labs(subtitle = "Quartiles of expression (Q1–Q4) within CNA = 1") +
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
        paste0("KM_OS_exprCNA_quartiles_CNA1_", make.names(pathway_name), model_suffix, ".pdf")
      )
      
      pdf(quart1_file, width = 8, height = 6)
      gridExtra::grid.arrange(p1$plot, p1$table, heights = c(0.7, 0.3))
      dev.off()
      
      message("Saved TCGA-BRCA OS PATHWAY quartile KM plot (CNA=1): ", quart1_file)
    } else {
      add_line("Quartile KM in CNA=1: insufficient data or < 2 quartiles present, skipping.")
      add_line("----")
    }
  }
  
  ## 5) Refit Cox model: Surv ~ expr + cna + expr:cna (PATHWAY-level scores)
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
      
      add_line("Refitted Cox (Surv ~ expr + cna + expr:cna):")
      add_line(sprintf("  expr:cna: HR=%.3f, 95%% CI [%.3f, %.3f], p=%.4g",
                       HR_int, LCL_int, UCL_int, p_int))
    } else {
      add_line("Refitted Cox: expr:cna coefficient not found in summary.")
    }
    
    ## PH test (cox.zph)
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
  
  ## 6) Spline Cox model: Surv ~ ns(expr, 3) * cna; plot HR vs expr for cna=1
  spline_pdf_file <- file.path(
    out_dir,
    paste0("Cox_spline_exprCNA_", make.names(pathway_name), model_suffix, ".pdf")
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
    expr_q2   <- stats::quantile(d$expr, probs = c(0.05, 0.95), na.rm = TRUE)
    expr_seq  <- seq(expr_q2[1], expr_q2[2], length.out = 50)
    
    med_expr_cna1 <- stats::median(d$expr[d$cna == 1], na.rm = TRUE)
    if (is.na(med_expr_cna1)) {
      med_expr_cna1 <- stats::median(d$expr, na.rm = TRUE)
    }
    
    newdata_seq <- data.frame(expr = expr_seq, cna = 1)
    newdata_ref <- data.frame(expr = med_expr_cna1, cna = 1)
    
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
        title    = paste0(pathway_name, " - Spline Cox for expr (cna = 1, OS, TCGA-BRCA)"),
        subtitle = "Model: Surv ~ ns(expr, 3) * cna (relative hazard vs median expr, cna = 1)",
        x        = "Pathway expression score (z-like)",
        y        = "Relative hazard (cna = 1)"
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
  
  ## Write summary to text file
  txt_file <- file.path(
    out_dir,
    paste0("KM_OS_exprCNA_interaction_details_", make.names(pathway_name), model_suffix, ".txt")
  )
  
  writeLines(summary_lines, con = txt_file)
  message("Extra expr:CNA diagnostics saved for pathway ", pathway_name, " to: ", txt_file)
  
  invisible(txt_file)
}

## --------------------------------------------------
## 11) Read TCGA-BRCA OS PATHWAY Cox Excel (NEW sheet names) and loop
## --------------------------------------------------

cox_file_tcga_OS_pathway <- file.path(
  tcga_results_dir,
  "TCGA_BRCA_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx"
)

all_sheets_tcga_OS_pathway <- readxl::excel_sheets(cox_file_tcga_OS_pathway)

## Exclude sheets whose names contain:
## "_filter", "_filter1", "_filter2", "LR_"
valid_sheets_tcga_OS_pathway <- all_sheets_tcga_OS_pathway[
  !grepl("_filter",  all_sheets_tcga_OS_pathway) &
    !grepl("_filter1", all_sheets_tcga_OS_pathway) &
    !grepl("_filter2", all_sheets_tcga_OS_pathway) &
    !grepl("LR_",      all_sheets_tcga_OS_pathway)
]

message("PATHWAY sheets to be processed (TCGA-BRCA OS): ", paste(valid_sheets_tcga_OS_pathway, collapse = ", "))

for (sh in valid_sheets_tcga_OS_pathway) {
  message("\n==============================")
  message("Processing PATHWAY sheet (TCGA-BRCA OS): ", sh)
  
  df <- readxl::read_xlsx(cox_file_tcga_OS_pathway, sheet = sh)
  
  ## Feature column: prefer "pathway", fallback to "gene" if needed
  feat_col <- NULL
  if ("pathway" %in% names(df)) {
    feat_col <- "pathway"
  } else if ("gene" %in% names(df)) {
    feat_col <- "gene"
  } else {
    message("Sheet ", sh, " does not contain 'pathway' or 'gene' column. Skipping.")
    next
  }
  
  ## Skip folder creation for Pathway_scores sheet
  if (sh == "Pathway_scores") {
    sheet_dir <- km_dir  ## use parent folder, do not create subfolder
  } else {
    sheet_dir <- file.path(km_dir, sh)
    if (!dir.exists(sheet_dir)) {
      dir.create(sheet_dir, recursive = TRUE)
    }
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
    
    ## Extract HR values (may be NA)
    HR_expr <- if ("HR_expr" %in% names(df)) as.numeric(df$HR_expr[i]) else NA_real_
    HR_cna  <- if ("HR_cna"  %in% names(df)) as.numeric(df$HR_cna[i])  else NA_real_
    HR_int  <- if ("HR_expr_cna" %in% names(df)) as.numeric(df$HR_expr_cna[i]) else NA_real_
    
    ## CI conditions
    cond_expr_raw    <- has_ci_expr    && !is.na(df$CI_lo_expr[i])     && df$CI_lo_expr[i]     > 1.00
    cond_cna_raw     <- has_ci_cna     && !is.na(df$CI_lo_cna[i])      && df$CI_lo_cna[i]      > 1.00
    cond_exprcna_raw <- has_ci_exprcna && !is.na(df$CI_lo_expr_cna[i]) && df$CI_lo_expr_cna[i] > 1.00
    
    ## Add HR > 1.00 condition
    cond_expr    <- cond_expr_raw    && !is.na(HR_expr) && HR_expr > 1.00
    cond_cna     <- cond_cna_raw     && !is.na(HR_cna)  && HR_cna  > 1.00
    cond_exprcna <- cond_exprcna_raw && !is.na(HR_int)  && HR_int  > 1.00
    
    ## If none of the conditions are satisfied, skip this pathway
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
      
      ## Median cut (log-rank p < 0.05 enforced inside)
      plot_km_expr_median_tcga_OS_pathway(
        pathway_name  = pathway_name,
        HR            = HR_expr,
        LCL           = LCL_expr,
        UCL           = UCL_expr,
        PVAL          = PVAL_expr,
        out_dir       = sheet_dir,
        model_suffix  = model_suffix,
        require_sig_p = TRUE
      )
      
      ## Optimal cut (log-rank p < 0.05 enforced inside)
      plot_km_expr_optimal_tcga_OS_pathway(
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
    
    ## CNA KM (2 groups; log-rank p < 0.05 inside)
    if (cond_cna && all(c("HR_cna", "CI_lo_cna", "CI_hi_cna", "p_cna") %in% names(df))) {
      LCL_cna  <- as.numeric(df$CI_lo_cna[i])
      UCL_cna  <- as.numeric(df$CI_hi_cna[i])
      PVAL_cna <- as.numeric(df$p_cna[i])
      
      plot_km_cna_tcga_OS_pathway(
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
    
    ## Expr:CNA interaction KM + extra diagnostics
    if (cond_exprcna && all(c("HR_expr_cna", "CI_lo_expr_cna", "CI_hi_expr_cna", "p_expr_cna") %in% names(df))) {
      LCL_int  <- as.numeric(df$CI_lo_expr_cna[i])
      UCL_int  <- as.numeric(df$CI_hi_expr_cna[i])
      PVAL_int <- as.numeric(df$p_expr_cna[i])
      
      ## Main interaction KM plot
      plot_km_expr_cna_interaction_tcga_OS_pathway(
        pathway_name = pathway_name,
        HR           = HR_int,
        LCL          = LCL_int,
        UCL          = UCL_int,
        PVAL         = PVAL_int,
        out_dir      = sheet_dir,
        model_suffix = model_suffix,
        four_groups  = has_EC_name
      )
      
      ## Extra diagnostics + quartile KMs for ALL interaction pathways
      extra_expr_cna_diagnostics_tcga_OS_pathway(
        pathway_name = pathway_name,
        out_dir      = sheet_dir,
        model_suffix = model_suffix,
        four_groups  = has_EC_name
      )
    } else if (cond_exprcna) {
      message("    (exprCNA) HR_expr_cna/CI/p columns missing for ", pathway_name, " in sheet ", sh, " – skipping interaction plots.")
    }
  } ## end for each row/pathway
}

message("\nAll TCGA-BRCA OS PATHWAY KM plots and interaction diagnostics have been saved under: ", km_dir)