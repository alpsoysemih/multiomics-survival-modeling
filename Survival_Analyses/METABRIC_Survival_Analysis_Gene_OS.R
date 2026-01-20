## ==================================================
## METABRIC OS (all patients): Gene-based Cox + LRT (with expr:CNA interaction + filtered sheets)
## ==================================================

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(writexl))

## --------------------------------------------------
## 0) Helper functions
## --------------------------------------------------

get_val <- function(tf, term, col) {
  x <- tf[[col]][tf$term == term]
  if (length(x) == 0) NA_real_ else x
}

clean_name <- function(x) {
  x %>%
    str_replace_all("[^[:alnum:]_]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

## Map full covariate tokens to short ones in model names
shorten_model_name <- function(model_name) {
  # If this is a filtered sheet name, keep the suffix as is
  is_filtered <- grepl("_filtered$", model_name)
  base_name   <- if (is_filtered) sub("_filtered$", "", model_name) else model_name
  
  parts <- strsplit(base_name, "_")[[1]]
  parts_short <- vapply(parts, function(tok) {
    if (tok == "expr") {
      "E"
    } else if (tok == "cna") {
      "C"
    } else if (tok == "exprCNA") {
      "EC"
    } else if (tok == "Stage") {
      "Stg"
    } else {
      tok
    }
  }, character(1))
  
  new_base <- paste(parts_short, collapse = "_")
  if (is_filtered) {
    paste0(new_base, "_filtered")
  } else {
    new_base
  }
}

## --------------------------------------------------
## 1) Directories
## --------------------------------------------------

metabric_dir  <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
pathfindR_dir <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"
gdsc_dir      <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/GDSC"
metabric_results_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/METABRIC"

if (!dir.exists(metabric_results_dir)) {
  dir.create(metabric_results_dir, recursive = TRUE)
}

## --------------------------------------------------
## 2) METABRIC data (expression, CNA, survival)
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

## Make CNA binary (0 vs != 0)
metabric_cna_binary <- metabric_cna
metabric_cna_binary[, -1] <- ifelse(metabric_cna_binary[, -1] != 0, 1, 0)

## --------------------------------------------------
## 3) METABRIC OS (all patients): survival + covariates (Age, Stage_simple, PAM50)
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

## --------------------------------------------------
## 3b) Drop stage categories with < 5 patients in OS subset
## --------------------------------------------------

stage_counts <- table(metabric_surv_small$Stage_simple)
rare_stages  <- names(stage_counts[stage_counts < 5])

if (length(rare_stages) > 0) {
  metabric_surv_small <- metabric_surv_small %>%
    dplyr::filter(!Stage_simple %in% rare_stages) %>%
    dplyr::mutate(Stage_simple = droplevels(Stage_simple))
}

## --------------------------------------------------
## 4) GDSC data (for Fisher OR direction at gene level)
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
## 5) pathfindR enrichment: get gene set to analyze
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

up_down_genes <- enrichment_df %>%
  dplyr::select(Up_regulated, Down_regulated) %>%
  tidyr::pivot_longer(cols = dplyr::everything(), values_to = "genes") %>%
  dplyr::filter(!is.na(genes), genes != "")

genes_union <- up_down_genes %>%
  dplyr::mutate(genes = stringr::str_split(genes, ",")) %>%
  tidyr::unnest(genes) %>%
  dplyr::mutate(genes = stringr::str_trim(genes)) %>%
  dplyr::filter(genes != "")

genes <- unique(genes_union$genes)
message("There are ", length(genes), " up- or down-regulated genes in enrichment result")

genes_metabric <- intersect(
  genes,
  intersect(
    intersect(metabric_expr$Gene_Symbol, metabric_cna_binary$Gene_Symbol),
    fisher_by_gene$Gene_Symbol
  )
)

message("Number of genes used for METABRIC OS gene-based Cox models: ",
        length(genes_metabric))

## --------------------------------------------------
## 6) METABRIC gene-level expression & CNA (long format)
## --------------------------------------------------

expr_long_metabric <- metabric_expr %>%
  filter(Gene_Symbol %in% genes_metabric) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "expr"
  )

cna_long_metabric <- metabric_cna_binary %>%
  filter(Gene_Symbol %in% genes_metabric) %>%
  pivot_longer(
    cols      = -Gene_Symbol,
    names_to  = "sample",
    values_to = "cna"
  )

cna_long_metabric_aligned <- cna_long_metabric %>%
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

## --------------------------------------------------
## 7) Model combinations (expr, cna, exprCNA, Age, Stage, PAM50)
##   - Now we create BOTH:
##     * expr + cna + covars   (no interaction)
##     * expr + cna + exprCNA + covars (with interaction)
## --------------------------------------------------

all_covars <- c("expr", "cna", "Age", "Stage", "PAM50")

model_list <- list()

for (k in 1:length(all_covars)) {
  cmb <- combn(all_covars, k, simplify = FALSE)
  for (s in cmb) {
    ## We only keep models that contain at least expr or cna
    if (!("expr" %in% s || "cna" %in% s)) next
    
    s_sorted <- intersect(all_covars, s)
    
    if ("expr" %in% s_sorted && "cna" %in% s_sorted) {
      ## --- 7a) model WITHOUT interaction: expr + cna + other covariates ---
      other_covars_no_int <- setdiff(s_sorted, c("expr", "cna"))
      model_name_no_int   <- paste(c("expr", "cna", other_covars_no_int), collapse = "_")
      covars_vec_no_int   <- c("expr", "cna", other_covars_no_int)
      model_list[[model_name_no_int]] <- covars_vec_no_int
      
      ## --- 7b) model WITH interaction: expr + cna + exprCNA + other covariates ---
      other_covars <- setdiff(s_sorted, c("expr", "cna"))
      model_name   <- paste(c("expr", "cna", "exprCNA", other_covars), collapse = "_")
      covars_vec   <- c("expr", "cna", "exprCNA", other_covars)
      model_list[[model_name]] <- covars_vec
      
    } else {
      ## Models that have only expr or only cna (no interaction)
      model_name <- paste(s_sorted, collapse = "_")
      covars_vec <- s_sorted
      model_list[[model_name]] <- covars_vec
    }
  }
}

## Baseline model without expression/CNA
model_list[["Age_Stage_PAM50"]] <- c("Age", "Stage", "PAM50")

model_specs <- tibble(
  model  = names(model_list),
  covars = unname(model_list)
)

## --------------------------------------------------
## 8) Cox function for a single gene + single model (expr, cna, expr:cna)
## --------------------------------------------------

fit_one_gene_one_model_metabric <- function(gene_name, covars, model_name) {
  d <- expr_long_metabric %>%
    filter(Gene_Symbol == gene_name) %>%
    dplyr::rename(expr = expr) %>%
    inner_join(
      cna_long_metabric_aligned %>%
        filter(Gene_Symbol == gene_name) %>%
        dplyr::select(Gene_Symbol, sample, cna_aligned),
      by = c("Gene_Symbol", "sample")
    ) %>%
    inner_join(metabric_surv_small, by = "sample") %>%
    transmute(
      gene         = Gene_Symbol,
      sample       = sample,
      expr         = as.numeric(expr),
      cna          = as.numeric(cna_aligned),
      OS_TIME      = OS_TIME,
      OS_EVENT     = OS_EVENT,
      Age          = as.numeric(Age),
      Stage_simple = Stage_simple,
      PAM50        = PAM50
    ) %>%
    distinct()
  
  n_events <- sum(d$OS_EVENT, na.rm = TRUE)
  if (nrow(d) < 30) return(NULL)
  if (n_events < 5) return(NULL)
  
  if ("expr" %in% covars && length(unique(na.omit(d$expr))) < 2) return(NULL)
  if ("cna"  %in% covars && length(unique(na.omit(d$cna)))  < 2) return(NULL)
  
  if ("expr" %in% covars) {
    d$expr <- scale(d$expr)[, 1]
  }
  
  ## Build RHS terms for the Cox model
  rhs_terms <- c()
  if ("expr"   %in% covars) rhs_terms <- c(rhs_terms, "expr")
  if ("cna"    %in% covars) rhs_terms <- c(rhs_terms, "cna")
  if ("Age"    %in% covars) rhs_terms <- c(rhs_terms, "Age")
  if ("Stage"  %in% covars) rhs_terms <- c(rhs_terms, "Stage_simple")
  if ("PAM50"  %in% covars) rhs_terms <- c(rhs_terms, "PAM50")
  
  ## If this model is flagged as having an expr–CNA interaction, add expr:cna
  if ("exprCNA" %in% covars) {
    rhs_terms <- c(rhs_terms, "expr:cna")
  }
  
  rhs <- paste(rhs_terms, collapse = " + ")
  fml <- as.formula(paste("Surv(OS_TIME, OS_EVENT) ~", rhs))
  
  fit <- try(coxph(fml, data = d, ties = "efron"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  tf      <- tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  cox_sum <- summary(fit)
  
  c_index    <- unname(cox_sum$concordance[1])
  c_index_se <- unname(cox_sum$concordance[2])
  logtest_p  <- cox_sum$logtest["pvalue"]
  wald_p     <- cox_sum$waldtest["pvalue"]
  score_p    <- cox_sum$sctest["pvalue"]
  
  row_list <- list(
    gene       = gene_name,
    model      = model_name,
    n          = nrow(d),
    events     = n_events,
    c_index    = c_index,
    c_index_se = c_index_se,
    logtest_p  = logtest_p,
    wald_p     = wald_p,
    score_p    = score_p
  )
  
  if ("expr" %in% covars) {
    row_list$HR_expr    <- get_val(tf, "expr", "estimate")
    row_list$CI_lo_expr <- get_val(tf, "expr", "conf.low")
    row_list$CI_hi_expr <- get_val(tf, "expr", "conf.high")
    row_list$p_expr     <- get_val(tf, "expr", "p.value")
  }
  
  if ("cna" %in% covars) {
    row_list$HR_cna    <- get_val(tf, "cna", "estimate")
    row_list$CI_lo_cna <- get_val(tf, "cna", "conf.low")
    row_list$CI_hi_cna <- get_val(tf, "cna", "conf.high")
    row_list$p_cna     <- get_val(tf, "cna", "p.value")
  }
  
  if ("expr" %in% covars && "cna" %in% covars && "exprCNA" %in% covars) {
    row_list$HR_expr_cna    <- get_val(tf, "expr:cna", "estimate")
    row_list$CI_lo_expr_cna <- get_val(tf, "expr:cna", "conf.low")
    row_list$CI_hi_expr_cna <- get_val(tf, "expr:cna", "conf.high")
    row_list$p_expr_cna     <- get_val(tf, "expr:cna", "p.value")
  }
  
  if ("Age" %in% covars) {
    row_list$HR_Age    <- get_val(tf, "Age", "estimate")
    row_list$CI_lo_Age <- get_val(tf, "Age", "conf.low")
    row_list$CI_hi_Age <- get_val(tf, "Age", "conf.high")
    row_list$p_Age     <- get_val(tf, "Age", "p.value")
  }
  
  if ("Stage" %in% covars) {
    stage_rows <- tf %>%
      dplyr::filter(stringr::str_starts(term, "Stage_simple"))
    
    if (nrow(stage_rows) > 0) {
      for (i in seq_len(nrow(stage_rows))) {
        term_i      <- stage_rows$term[i]
        stage_label <- sub("^Stage_simple", "Stage", term_i)
        safe_nm     <- clean_name(stage_label)
        
        hr_name <- paste0("HR_",    safe_nm)
        lo_name <- paste0("CI_lo_", safe_nm)
        hi_name <- paste0("CI_hi_", safe_nm)
        p_name  <- paste0("p_",     safe_nm)
        
        row_list[[hr_name]] <- stage_rows$estimate[i]
        row_list[[lo_name]] <- stage_rows$conf.low[i]
        row_list[[hi_name]] <- stage_rows$conf.high[i]
        row_list[[p_name]]  <- stage_rows$p.value[i]
      }
    }
  }
  
  if ("PAM50" %in% covars) {
    pam_rows <- tf %>%
      dplyr::filter(str_starts(term, "PAM50"))
    
    if (nrow(pam_rows) > 0) {
      for (i in seq_len(nrow(pam_rows))) {
        term_i  <- pam_rows$term[i]
        safe_nm <- clean_name(term_i)
        
        hr_name <- paste0("HR_",    safe_nm)
        lo_name <- paste0("CI_lo_", safe_nm)
        hi_name <- paste0("CI_hi_", safe_nm)
        p_name  <- paste0("p_",     safe_nm)
        
        row_list[[hr_name]] <- pam_rows$estimate[i]
        row_list[[lo_name]] <- pam_rows$conf.low[i]
        row_list[[hi_name]] <- pam_rows$conf.high[i]
        row_list[[p_name]]  <- pam_rows$p.value[i]
      }
    }
  }
  
  as_tibble(row_list)
}

## --------------------------------------------------
## 9) Fit all models for all genes (METABRIC OS, all)
## --------------------------------------------------

results_by_model_metabric <- list()

for (i in seq_len(nrow(model_specs))) {
  model_name <- model_specs$model[i]
  covars     <- model_specs$covars[[i]]
  
  message("=== METABRIC OS (all) | Fitting model: ", model_name,
          " | covariates: ", paste(covars, collapse = ", "), " ===")
  
  model_rows <- map(
    genes_metabric,
    ~ fit_one_gene_one_model_metabric(.x, covars, model_name)
  )
  
  model_df <- model_rows %>%
    compact() %>%
    bind_rows()
  
  if (nrow(model_df) > 0) {
    results_by_model_metabric[[model_name]] <- model_df
    message(">>> METABRIC OS (all) | Finished model: ", model_name,
            " | number of fitted genes: ",
            length(model_rows) - sum(sapply(model_rows, is.null)))
  } else {
    message("!!! METABRIC OS (all) | No suitable gene for model: ", model_name)
  }
}

## --------------------------------------------------
## 10) LRT helper functions (with expr:cna included when applicable)
## --------------------------------------------------

get_gene_data_metabric <- function(gene_name) {
  expr_long_metabric %>%
    filter(Gene_Symbol == gene_name) %>%
    dplyr::rename(expr = expr) %>%
    inner_join(
      cna_long_metabric_aligned %>%
        filter(Gene_Symbol == gene_name) %>%
        dplyr::select(Gene_Symbol, sample, cna_aligned),
      by = c("Gene_Symbol", "sample")
    ) %>%
    inner_join(metabric_surv_small, by = "sample") %>%
    transmute(
      gene         = Gene_Symbol,
      sample       = sample,
      expr         = as.numeric(expr),
      cna          = as.numeric(cna_aligned),
      OS_TIME      = OS_TIME,
      OS_EVENT     = OS_EVENT,
      Age          = as.numeric(Age),
      Stage_simple = Stage_simple,
      PAM50        = PAM50
    ) %>%
    distinct()
}

fit_pair_lrt_gene_metabric <- function(gene_name, d, vars_small, vars_big,
                                       name_small, name_big) {
  ## Map "Stage" → "Stage_simple" but keep the exprCNA flag
  vars_small_data <- ifelse(vars_small == "Stage", "Stage_simple", vars_small)
  vars_big_data   <- ifelse(vars_big   == "Stage", "Stage_simple", vars_big)
  
  ## exprCNA is NOT a column in d, it is only a flag for adding expr:cna term
  ## So we remove it from the list of columns we will subset from d
  vars_small_no_int <- setdiff(vars_small_data, "exprCNA")
  vars_big_no_int   <- setdiff(vars_big_data,   "exprCNA")
  
  ## All real variables that must be present in data
  all_vars <- unique(c(vars_small_no_int, vars_big_no_int))
  needed   <- c("OS_TIME", "OS_EVENT", all_vars)
  
  d2 <- d[, needed, drop = FALSE]
  d2 <- d2[complete.cases(d2), , drop = FALSE]
  if (nrow(d2) < 30) return(NULL)
  n_events <- sum(d2$OS_EVENT, na.rm = TRUE)
  if (n_events < 5) return(NULL)
  
  ## Scale expression if present
  if ("expr" %in% all_vars) {
    d2$expr <- scale(d2$expr)[, 1]
  }
  
  ## Build RHS for small model
  rhs_small_terms <- character(0)
  if ("expr"         %in% vars_small_no_int) rhs_small_terms <- c(rhs_small_terms, "expr")
  if ("cna"          %in% vars_small_no_int) rhs_small_terms <- c(rhs_small_terms, "cna")
  if ("Age"          %in% vars_small_no_int) rhs_small_terms <- c(rhs_small_terms, "Age")
  if ("Stage_simple" %in% vars_small_no_int) rhs_small_terms <- c(rhs_small_terms, "Stage_simple")
  if ("PAM50"        %in% vars_small_no_int) rhs_small_terms <- c(rhs_small_terms, "PAM50")
  ## If exprCNA flag is present in the small model, add expr:cna
  if ("exprCNA" %in% vars_small_data) {
    rhs_small_terms <- c(rhs_small_terms, "expr:cna")
  }
  rhs_small <- paste(rhs_small_terms, collapse = " + ")
  fml_small <- as.formula(paste("Surv(OS_TIME, OS_EVENT) ~", rhs_small))
  
  ## Build RHS for big model
  rhs_big_terms <- character(0)
  if ("expr"         %in% vars_big_no_int) rhs_big_terms <- c(rhs_big_terms, "expr")
  if ("cna"          %in% vars_big_no_int) rhs_big_terms <- c(rhs_big_terms, "cna")
  if ("Age"          %in% vars_big_no_int) rhs_big_terms <- c(rhs_big_terms, "Age")
  if ("Stage_simple" %in% vars_big_no_int) rhs_big_terms <- c(rhs_big_terms, "Stage_simple")
  if ("PAM50"        %in% vars_big_no_int) rhs_big_terms <- c(rhs_big_terms, "PAM50")
  ## If exprCNA flag is present in the big model, add expr:cna
  if ("exprCNA" %in% vars_big_data) {
    rhs_big_terms <- c(rhs_big_terms, "expr:cna")
  }
  rhs_big <- paste(rhs_big_terms, collapse = " + ")
  fml_big <- as.formula(paste("Surv(OS_TIME, OS_EVENT) ~", rhs_big))
  
  ## Fit the two nested models
  fit_small <- try(coxph(fml_small, data = d2, ties = "efron"), silent = TRUE)
  fit_big   <- try(coxph(fml_big,   data = d2, ties = "efron"), silent = TRUE)
  if (inherits(fit_small, "try-error") || inherits(fit_big, "try-error")) {
    return(NULL)
  }
  
  ## Likelihood ratio test between small and big
  a <- anova(fit_small, fit_big, test = "LRT")
  
  tibble(
    gene         = gene_name,
    small_model  = name_small,
    big_model    = name_big,
    df_small     = length(unique(vars_small)),
    df_big       = length(unique(vars_big)),
    logLik_small = as.numeric(a$loglik[1]),
    logLik_big   = as.numeric(a$loglik[2]),
    LRT_stat     = as.numeric(a$Chisq[2]),
    LRT_df       = as.numeric(a$Df[2]),
    LRT_p        = as.numeric(a$`Pr(>|Chi|)`[2])
  )
}

lr_tests_for_gene_metabric <- function(gene_name) {
  d <- get_gene_data_metabric(gene_name)
  n_events <- sum(d$OS_EVENT, na.rm = TRUE)
  if (nrow(d) < 30 || n_events < 5) return(NULL)
  
  lr_expr_list <- list()
  lr_cna_list  <- list()
  lr_base_list <- list()
  
  model_names <- names(model_list)
  
  for (i in seq_along(model_names)) {
    for (j in seq_along(model_names)) {
      if (j <= i) next
      
      name_i <- model_names[i]
      name_j <- model_names[j]
      cov_i  <- model_list[[name_i]]
      cov_j  <- model_list[[name_j]]
      
      if (all(cov_i %in% cov_j)) {
        vars_small <- cov_i
        vars_big   <- cov_j
        name_small <- name_i
        name_big   <- name_j
      } else if (all(cov_j %in% cov_i)) {
        vars_small <- cov_j
        vars_big   <- cov_i
        name_small <- name_j
        name_big   <- name_i
      } else {
        next
      }
      
      if (setequal(vars_small, vars_big)) next
      
      row <- fit_pair_lrt_gene_metabric(
        gene_name, d, vars_small, vars_big,
        name_small, name_big
      )
      if (is.null(row)) next
      
      ## classification for LRT outputs (expr, cna, baseline),
      ## now using the explicit "exprCNA" flag in the covariate sets
      
      vars_all <- union(vars_small, vars_big)
      
      has_expr    <- "expr"    %in% vars_all
      has_cna     <- "cna"     %in% vars_all
      has_exprCNA <- "exprCNA" %in% vars_all
      
      ## 1) Any comparison that involves expression OR exprCNA → LR_expr sheet
      if (has_expr || has_exprCNA) {
        lr_expr_list[[length(lr_expr_list) + 1L]] <- row
      }
      
      ## 2) Any comparison that involves CNA OR exprCNA → LR_cna sheet
      if (has_cna || has_exprCNA) {
        lr_cna_list[[length(lr_cna_list) + 1L]]  <- row
      }
      
      ## 3) Baseline sheet:
      ##    small model = Age + Stage + PAM50
      ##    big model = baseline + any of expr / cna / exprCNA
      if (setequal(vars_small, c("Age", "Stage", "PAM50"))) {
        lr_base_list[[length(lr_base_list) + 1L]] <- row
      }
    }
  }
  
  list(
    expr = lr_expr_list %>% compact() %>% bind_rows(),
    cna  = lr_cna_list  %>% compact() %>% bind_rows(),
    base = lr_base_list %>% compact() %>% bind_rows()
  )
}

## --------------------------------------------------
## 11) Collect LRT results for all genes (METABRIC OS all)
## --------------------------------------------------

lr_expr_list_all_metabric <- list()
lr_cna_list_all_metabric  <- list()
lr_base_list_all_metabric <- list()

for (g in genes_metabric) {
  lr_res <- lr_tests_for_gene_metabric(g)
  if (is.null(lr_res)) next
  
  if (!is.null(lr_res$expr) && nrow(lr_res$expr) > 0) {
    lr_expr_list_all_metabric[[g]] <- lr_res$expr
  }
  if (!is.null(lr_res$cna)  && nrow(lr_res$cna)  > 0) {
    lr_cna_list_all_metabric[[g]]  <- lr_res$cna
  }
  if (!is.null(lr_res$base) && nrow(lr_res$base) > 0) {
    lr_base_list_all_metabric[[g]] <- lr_res$base
  }
}

lr_expr_df_metabric <- lr_expr_list_all_metabric %>%
  compact() %>%
  bind_rows() %>%
  distinct(gene, small_model, big_model, .keep_all = TRUE)

lr_cna_df_metabric <- lr_cna_list_all_metabric %>%
  compact() %>%
  bind_rows() %>%
  distinct(gene, small_model, big_model, .keep_all = TRUE)

lr_base_df_metabric <- lr_base_list_all_metabric %>%
  compact() %>%
  bind_rows() %>%
  distinct(gene, small_model, big_model, .keep_all = TRUE)

## --------------------------------------------------
## 12) p-adjust (BH) for Cox models, rounding, sorting, filtered sheets
## --------------------------------------------------

total_genes <- length(genes_metabric)
message("Minimum number of hypotheses used for BH correction (total genes): ",
        total_genes)

process_model_df <- function(df) {
  if (nrow(df) == 0) return(df)
  
  ## BH-adjusted p-values for all p-columns
  p_cols_raw <- names(df)[grepl("(^p_|_p$)", names(df))]
  for (pc in p_cols_raw) {
    adj_col <- paste0(pc, "_adjusted")
    p_vec   <- suppressWarnings(as.numeric(df[[pc]]))
    lp      <- sum(!is.na(p_vec))
    if (lp == 0) {
      df[[adj_col]] <- NA_real_
    } else {
      df[[adj_col]] <- p.adjust(
        p_vec,
        method = "BH",
        n = max(total_genes, lp)
      )
    }
    name_vec <- names(df)
    name_vec_no_adj <- name_vec[name_vec != adj_col]
    idx <- match(pc, name_vec_no_adj)
    new_order <- append(name_vec_no_adj, adj_col, after = idx)
    df <- df[, new_order]
  }
  
  ## Round HR, CI and c_index/c_index_se
  hr_ci_cols <- names(df)[grepl("^(HR|CI_lo|CI_hi|c_index|c_index_se)", names(df))]
  if (length(hr_ci_cols) > 0) {
    df[hr_ci_cols] <- lapply(df[hr_ci_cols], function(x) {
      xnum <- suppressWarnings(as.numeric(x))
      ifelse(is.na(xnum), NA_real_, signif(xnum, 3))
    })
  }
  
  ## Round all p and adjusted p columns
  p_cols_all <- names(df)[grepl("(^p_|_p$|_adjusted$)", names(df))]
  p_cols_all <- unique(p_cols_all[p_cols_all %in% names(df)])
  
  if (length(p_cols_all) > 0) {
    df[p_cols_all] <- lapply(df[p_cols_all], function(x) {
      xnum <- suppressWarnings(as.numeric(x))
      ifelse(is.na(xnum), NA_real_, signif(xnum, 3))
    })
  }
  
  ## Sort rows: first by CI_lo_cna (if present), else by CI_lo_expr (if present), both descending
  if ("CI_lo_cna" %in% names(df)) {
    df <- df %>% arrange(desc(CI_lo_cna))
  } else if ("CI_lo_expr" %in% names(df)) {
    df <- df %>% arrange(desc(CI_lo_expr))
  }
  
  df
}

results_by_model_metabric_proc <- lapply(results_by_model_metabric, process_model_df)

## Filtered sheets:
##  - keep genes with CI_lo_expr > 1 OR CI_lo_cna > 1 OR (if present) CI_lo_expr_cna > 1
##  - Do NOT create a filtered sheet for the baseline model "Age_Stage_PAM50"
filtered_by_model_metabric <- list()

for (nm in names(results_by_model_metabric_proc)) {
  df <- results_by_model_metabric_proc[[nm]]
  
  if (nm == "Age_Stage_PAM50" || nrow(df) == 0) {
    next
  }
  
  mask_expr <- rep(FALSE, nrow(df))
  mask_cna  <- rep(FALSE, nrow(df))
  mask_int  <- rep(FALSE, nrow(df))
  
  if ("CI_lo_expr" %in% names(df)) {
    mask_expr <- !is.na(df$CI_lo_expr) & df$CI_lo_expr > 1.00
  }
  
  if ("CI_lo_cna" %in% names(df)) {
    mask_cna <- !is.na(df$CI_lo_cna) & df$CI_lo_cna > 1.00
  }
  
  if ("CI_lo_expr_cna" %in% names(df)) {
    mask_int <- !is.na(df$CI_lo_expr_cna) & df$CI_lo_expr_cna > 1.00
  }
  
  mask <- mask_expr | mask_cna | mask_int
  df_filt <- df[mask, , drop = FALSE]
  filtered_by_model_metabric[[paste0(nm, "_filtered")]] <- df_filt
}

## LRT: NO multiple-testing correction, only rounding
process_lrt_df <- function(df) {
  if (!("LRT_p" %in% names(df))) return(df)
  if (nrow(df) == 0) return(df)
  
  num_cols_lrt <- intersect(
    c("logLik_small", "logLik_big", "LRT_stat", "LRT_p"),
    names(df)
  )
  if (length(num_cols_lrt) > 0) {
    df[num_cols_lrt] <- lapply(df[num_cols_lrt], function(x) {
      xnum <- suppressWarnings(as.numeric(x))
      ifelse(is.na(xnum), NA_real_, signif(xnum, 3))
    })
  }
  
  df
}

lr_expr_df_metabric <- process_lrt_df(lr_expr_df_metabric)
lr_cna_df_metabric  <- process_lrt_df(lr_cna_df_metabric)
lr_base_df_metabric <- process_lrt_df(lr_base_df_metabric)

## Abbreviate model names inside LRT result data frames
if (nrow(lr_expr_df_metabric) > 0) {
  lr_expr_df_metabric <- lr_expr_df_metabric %>%
    mutate(
      small_model = vapply(small_model, shorten_model_name, character(1)),
      big_model   = vapply(big_model,   shorten_model_name, character(1))
    )
}

if (nrow(lr_cna_df_metabric) > 0) {
  lr_cna_df_metabric <- lr_cna_df_metabric %>%
    mutate(
      small_model = vapply(small_model, shorten_model_name, character(1)),
      big_model   = vapply(big_model,   shorten_model_name, character(1))
    )
}

if (nrow(lr_base_df_metabric) > 0) {
  lr_base_df_metabric <- lr_base_df_metabric %>%
    mutate(
      small_model = vapply(small_model, shorten_model_name, character(1)),
      big_model   = vapply(big_model,   shorten_model_name, character(1))
    )
}

## --------------------------------------------------
## 13) Excel export:
##   - main models
##   - _filter versions (CI_lo_expr / CI_lo_cna / CI_lo_expr_cna > 1.00)
##   - for expr:CNA models: extra ..._filter2 with only CI_lo_expr_cna > 1.00
## --------------------------------------------------

xls_models <- list()

for (nm in names(results_by_model_metabric_proc)) {
  df <- results_by_model_metabric_proc[[nm]]
  
  ## --- shorten the model name (E / C / EC / Age / Stg / PAM50) ---
  nm_short <- shorten_model_name(nm)
  
  ## --- main sheet ---
  xls_models[[nm_short]] <- df
  
  ## --- filtered sheet (union of expr / cna / expr_cna) ---
  if (nm != "Age_Stage_PAM50") {
    filt_name <- paste0(nm, "_filtered")
    
    df_filt <- filtered_by_model_metabric[[filt_name]]
    
    if (!is.null(df_filt)) {
      ## rename "..._filtered" → "..._filter"
      nm_filt_short <- paste0(nm_short, "_filter")
      xls_models[[nm_filt_short]] <- df_filt
    }
    
    ## --- filter2 sheet: only CI_lo_expr_cna > 1.00 for interaction models ---
    if ("CI_lo_expr_cna" %in% names(df)) {
      df_filt2 <- df %>%
        dplyr::filter(!is.na(CI_lo_expr_cna) & CI_lo_expr_cna > 1.00) %>%
        dplyr::arrange(desc(CI_lo_expr_cna))
      
      nm_filt2_short <- paste0(nm_short, "_filter2")
      xls_models[[nm_filt2_short]] <- df_filt2
    }
  }
}

out_file_metabric <- file.path(
  metabric_results_dir,
  "METABRIC_OS_all_Cox_all_models_by_gene_with_LRT.xlsx"
)

write_xlsx(
  c(
    xls_models,
    list(
      LR_expr     = lr_expr_df_metabric,
      LR_cna      = lr_cna_df_metabric,
      LR_baseline = lr_base_df_metabric
    )
  ),
  out_file_metabric
)

message(
  "All gene-level Cox models for METABRIC (OS, all patients) with BOTH expr+cna and expr+cna+expr:CNA models, ",
  "BH-adjusted Cox p-values, CI_lo_expr/CI_lo_cna/CI_lo_expr_cna > 1.00 filtered (_filter) sheets, ",
  "extra _filter2 sheets for CI_lo_expr_cna > 1.00, ",
  "and raw (uncorrected) LRT results (including models with expr:CNA interaction and non-interaction) ",
  "have been written to: ",
  out_file_metabric
)