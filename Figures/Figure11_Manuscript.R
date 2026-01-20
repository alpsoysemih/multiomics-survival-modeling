# Load required packages quietly
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggnewscale))

# Global option: keep strings as character (not factors) by default
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Directories
# ------------------------------------------------------------

# Define input directories
pathfindR_dir        <- "/Volumes/Expansion/Prognostic_Analysis/Pathway Enrichment Result"
metabric_results_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/METABRIC"
tcga_results_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Results/TCGA-BRCA"

# Define output directory for the final figure
figures_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}

# Define the final output file path (only this file will be written)
network_file <- file.path(figures_dir, "Figure 11.jpeg")

# Print the output file path for the network figure
cat(">> Network figure will be saved as:\n   ", network_file, "\n\n")

# ------------------------------------------------------------
# Read pathfindR enrichment results (Paclitaxel)
# ------------------------------------------------------------

# Define the pathfindR enrichment file and stop if it is missing
pathfindR_file <- file.path(pathfindR_dir, "pathfindR_paclitaxel.tsv")
if (!file.exists(pathfindR_file)) {
  stop("Cannot find pathfindR file: ", pathfindR_file)
}

# Load enrichment results
enrichment_df <- readr::read_tsv(pathfindR_file, show_col_types = FALSE)

# ------------------------------------------------------------
# Define key pathway and gene sets
# ------------------------------------------------------------

# Define signature genes (expression + CNA)
sig_genes_all <- c("SERPINE1", "CLDN1", "CXCR4", "LAMB2", "LYPD6B")

# Define pathway names of interest (must match Term_Description)
ras_pathway_name   <- "Ras signaling pathway"
focal_pathway_name <- "Focal adhesion"
axon_pathway_name  <- "Axon guidance"

# Define KEGG pathway IDs selected for the network
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

# ------------------------------------------------------------
# Build pathway–gene table for network
# ------------------------------------------------------------

# Build a long pathway-to-gene table with expression direction (up/down) from selected KEGG terms
path_genes_net <- enrichment_df %>%
  dplyr::filter(ID %in% selected_terms) %>%
  dplyr::select(ID, Term_Description, Up_regulated, Down_regulated) %>%
  dplyr::rename(pathway = Term_Description) %>%
  dplyr::mutate(
    pathway        = stringr::str_trim(as.character(pathway)),
    Up_regulated   = stringr::str_split(Up_regulated,   pattern = "[,;]"),
    Down_regulated = stringr::str_split(Down_regulated, pattern = "[,;]")
  ) %>%
  tidyr::pivot_longer(
    cols      = c("Up_regulated", "Down_regulated"),
    names_to  = "direction_raw",
    values_to = "genes"
  ) %>%
  tidyr::unnest(genes) %>%
  dplyr::mutate(
    Gene_Symbol = stringr::str_trim(genes),
    expr_dir = dplyr::case_when(
      direction_raw == "Up_regulated"   ~ "up",
      direction_raw == "Down_regulated" ~ "down",
      TRUE                              ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Gene_Symbol), Gene_Symbol != "") %>%
  dplyr::distinct(ID, pathway, Gene_Symbol, expr_dir)

# Identify pathways that contain at least one signature gene
path_with_sig <- path_genes_net %>%
  dplyr::filter(Gene_Symbol %in% sig_genes_all) %>%
  dplyr::pull(pathway) %>%
  unique()

# Define the final set of pathways to keep (target pathways + pathways containing signature genes)
pathways_keep <- union(
  c(ras_pathway_name, focal_pathway_name, axon_pathway_name),
  path_with_sig
)

# Filter network table to retained pathways
path_genes_net_filt <- path_genes_net %>%
  dplyr::filter(pathway %in% pathways_keep)

# Keep all genes in the selected pathways (no additional gene filtering)
path_genes_net_filt2 <- path_genes_net_filt

# Report the number of pathways and genes retained for the network
cat("   Pathways kept: ", length(unique(path_genes_net_filt2$pathway)), "\n")
cat("   Genes kept (all genes in selected pathways): ",
    length(unique(path_genes_net_filt2$Gene_Symbol)), "\n\n")

# ------------------------------------------------------------
# Summarize expression direction per gene (up / down / mixed)
# ------------------------------------------------------------

# Summarize gene-level direction across all retained pathways
gene_expr_summary <- path_genes_net_filt2 %>%
  dplyr::group_by(Gene_Symbol) %>%
  dplyr::summarise(
    any_up   = any(expr_dir == "up",   na.rm = TRUE),
    any_down = any(expr_dir == "down", na.rm = TRUE),
    expr_dir_gene = dplyr::case_when(
      any_up & !any_down ~ "up",
      any_down & !any_up ~ "down",
      any_up & any_down  ~ "mixed",
      TRUE               ~ "unknown"
    ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Network edges (pathway–gene)
# ------------------------------------------------------------

# Define genes to highlight as signature nodes in the network
highlight_genes <- c("SERPINE1", "CLDN1", "LYPD6B", "LAMB2", "CXCR4")

# Define genes to highlight as signature nodes in the network
highlight_genes <- c("SERPINE1", "CLDN1", "LYPD6B", "LAMB2", "CXCR4")

# ---- NEW: genes that are prognostic-associated but NOT independent (NOT in signature)
# (sig_gene_union is computed later in your code; so we will define this AFTER sig_gene_union is created)
# We'll create a placeholder here; the real vector will be created after sig_gene_union.

edges_net <- path_genes_net_filt2 %>%
  dplyr::select(from = pathway, to = Gene_Symbol) %>%
  dplyr::mutate(edge_type = "other", edge_alpha = 0.6)

# ------------------------------------------------------------
# Helper functions: extract significant HR info and pool effects across cohorts
# ------------------------------------------------------------

# Extract significant nodes from an Excel sheet using lower CI bounds and compute log(HR) and SE
extract_sig_table <- function(file, sheet, name_col) {
  empty_tbl <- tibble::tibble(
    name       = character(),
    level      = character(),
    sig_type   = character(),
    HR_use     = numeric(),
    logHR_use  = numeric(),
    se_use     = numeric()
  )
  if (!file.exists(file)) {
    warning("File not found: ", file)
    return(empty_tbl)
  }
  df <- readxl::read_xlsx(file, sheet = sheet)
  if (!(name_col %in% colnames(df))) {
    warning("Column ", name_col, " not found in ", file, " / sheet ", sheet)
    return(empty_tbl)
  }
  df[[name_col]] <- df[[name_col]] %>%
    as.character() %>%
    stringr::str_trim() %>%
    sub("^[^:]+:\\s*", "", .)
  
  has_cna     <- all(c("CI_lo_cna",  "HR_cna")  %in% names(df))
  has_expr    <- all(c("CI_lo_expr", "HR_expr") %in% names(df))
  has_hi_cna  <- "CI_hi_cna"  %in% names(df)
  has_hi_expr <- "CI_hi_expr" %in% names(df)
  
  if (!has_cna & !has_expr) {
    warning("No CI_lo_cna/HR_cna or CI_lo_expr/HR_expr in ", file, " / sheet ", sheet)
    return(empty_tbl)
  }
  
  ## At least one of CNA or expression CI_lo > 1
  cond <- rep(FALSE, nrow(df))
  if (has_cna)  cond <- cond | (df$CI_lo_cna  > 1)
  if (has_expr) cond <- cond | (df$CI_lo_expr > 1)
  
  df2 <- df[cond & !is.na(df[[name_col]]), , drop = FALSE]
  if (nrow(df2) == 0) return(empty_tbl)
  
  sig_type <- dplyr::case_when(
    has_cna  & df2$CI_lo_cna  > 1 & !(has_expr & df2$CI_lo_expr > 1) ~ "cna",
    has_expr & df2$CI_lo_expr > 1 & !(has_cna  & df2$CI_lo_cna  > 1) ~ "expr",
    has_cna  & has_expr &
      df2$CI_lo_cna > 1 & df2$CI_lo_expr > 1 ~ "expr",  # <- BURAYI DEĞİŞTİR
    TRUE ~ NA_character_
  )
  
  HR_use <- dplyr::case_when(
    sig_type == "cna"  ~ df2$HR_cna,
    sig_type == "expr" ~ df2$HR_expr,
    TRUE               ~ NA_real_
  )
  logHR_use <- log(HR_use)
  
  se_cna  <- se_expr <- rep(NA_real_, nrow(df2))
  if (has_cna && has_hi_cna) {
    se_cna <- (log(df2$CI_hi_cna) - log(df2$CI_lo_cna)) / (2 * 1.96)
  }
  if (has_expr && has_hi_expr) {
    se_expr <- (log(df2$CI_hi_expr) - log(df2$CI_lo_expr)) / (2 * 1.96)
  }
  se_use <- dplyr::case_when(
    sig_type == "cna"  ~ se_cna,
    sig_type == "expr" ~ se_expr,
    TRUE               ~ NA_real_
  )
  
  out <- tibble::tibble(
    name       = df2[[name_col]],
    level      = ifelse(name_col == "gene", "gene", "pathway"),
    sig_type   = sig_type,
    HR_use     = HR_use,
    logHR_use  = logHR_use,
    se_use     = se_use
  ) %>%
    dplyr::filter(!is.na(sig_type), !is.na(HR_use), !is.na(logHR_use))
  
  if (nrow(out) == 0) return(empty_tbl)
  out
}

# Fixed-effect pooling of log(HR) across two cohorts (fallback to mean log(HR) if SEs are invalid)
meta_pool <- function(logHR1, se1, logHR2, se2, HR1, HR2) {
  w1 <- 1 / (se1^2)
  w2 <- 1 / (se2^2)
  valid <- is.finite(w1) & is.finite(w2) & (w1 > 0) & (w2 > 0)
  pooled_logHR <- ifelse(
    valid,
    (w1 * logHR1 + w2 * logHR2) / (w1 + w2),
    (log(HR1) + log(HR2)) / 2
  )
  pooled_logHR
}

# ------------------------------------------------------------
# Read METABRIC Cox results (gene and pathway)
# ------------------------------------------------------------

# Define METABRIC Cox result files and extract significant gene/pathway signals
met_gene_file  <- file.path(metabric_results_dir, "METABRIC_OS_all_Cox_all_models_by_gene_with_LRT.xlsx")
met_path_file  <- file.path(metabric_results_dir, "METABRIC_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx")

met_gene_EC    <- extract_sig_table(met_gene_file, "E_C_filter",               "gene")
met_gene_EC_AS <- extract_sig_table(met_gene_file, "E_C_Age_Stg_PAM50_filter", "gene")
met_path_EC    <- extract_sig_table(met_path_file, "E_C_filter",               "pathway")
met_path_EC_AS <- extract_sig_table(met_path_file, "E_C_Age_Stg_PAM50_filter", "pathway")

met_EC_all    <- dplyr::bind_rows(met_gene_EC,    met_path_EC)
met_EC_AS_all <- dplyr::bind_rows(met_gene_EC_AS, met_path_EC_AS)

# ------------------------------------------------------------
# Read TCGA-BRCA Cox results (gene and pathway)
# ------------------------------------------------------------

# Define TCGA Cox result files and extract significant gene/pathway signals
tcga_gene_file  <- file.path(tcga_results_dir, "TCGA_BRCA_OS_all_Cox_all_models_by_gene_with_LRT.xlsx")
tcga_path_file  <- file.path(tcga_results_dir, "TCGA_BRCA_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx")

tcga_gene_EC    <- extract_sig_table(tcga_gene_file, "E_C_filter",               "gene")
tcga_gene_EC_AS <- extract_sig_table(tcga_gene_file, "E_C_Age_Stg_PAM50_filter", "gene")
tcga_path_EC    <- extract_sig_table(tcga_path_file, "E_C_filter",               "pathway")
tcga_path_EC_AS <- extract_sig_table(tcga_path_file, "E_C_Age_Stg_PAM50_filter", "pathway")

tcga_EC_all    <- dplyr::bind_rows(tcga_gene_EC,    tcga_path_EC)
tcga_EC_AS_all <- dplyr::bind_rows(tcga_gene_EC_AS, tcga_path_EC_AS)

# ------------------------------------------------------------
# Build meta-analytic node-level HR information (two layers)
# ------------------------------------------------------------

## E_C_filter
common_EC <- dplyr::inner_join(
  met_EC_all  %>%
    dplyr::select(name, level,
                  sig_type_met  = sig_type,
                  HR_met        = HR_use,
                  logHR_met     = logHR_use,
                  se_met        = se_use),
  tcga_EC_all %>%
    dplyr::select(name, level,
                  sig_type_tcga = sig_type,
                  HR_tcga       = HR_use,
                  logHR_tcga    = logHR_use,
                  se_tcga       = se_use),
  by = c("name", "level")
)

common_EC_cna <- common_EC %>%
  dplyr::filter(sig_type_met == "cna", sig_type_tcga == "cna") %>%
  dplyr::mutate(
    shape_EC  = "cna",
    logHR_EC  = meta_pool(logHR_met, se_met, logHR_tcga, se_tcga,
                          HR_met, HR_tcga),
    HR_EC     = exp(logHR_EC)
  ) %>%
  dplyr::select(name, level, shape_EC, logHR_EC, HR_EC)

common_EC_expr <- common_EC %>%
  dplyr::filter(sig_type_met == "expr", sig_type_tcga == "expr") %>%
  dplyr::mutate(
    shape_EC  = "expr",
    logHR_EC  = meta_pool(logHR_met, se_met, logHR_tcga, se_tcga,
                          HR_met, HR_tcga),
    HR_EC     = exp(logHR_EC)
  ) %>%
  dplyr::select(name, level, shape_EC, logHR_EC, HR_EC)

common_EC_nodes <- dplyr::bind_rows(common_EC_cna, common_EC_expr) %>%
  dplyr::group_by(name, level) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

## E_C_Age_Stg_PAM50_filter
common_EC_AS <- dplyr::inner_join(
  met_EC_AS_all  %>%
    dplyr::select(name, level,
                  sig_type_met  = sig_type,
                  HR_met        = HR_use,
                  logHR_met     = logHR_use,
                  se_met        = se_use),
  tcga_EC_AS_all %>%
    dplyr::select(name, level,
                  sig_type_tcga = sig_type,
                  HR_tcga       = HR_use,
                  logHR_tcga    = logHR_use,
                  se_tcga       = se_use),
  by = c("name", "level")
)

common_EC_AS_cna <- common_EC_AS %>%
  dplyr::filter(sig_type_met == "cna", sig_type_tcga == "cna") %>%
  dplyr::mutate(
    shape_EC_AS  = "cna",
    logHR_EC_AS  = meta_pool(logHR_met, se_met, logHR_tcga, se_tcga,
                             HR_met, HR_tcga),
    HR_EC_AS     = exp(logHR_EC_AS)
  ) %>%
  dplyr::select(name, level, shape_EC_AS, logHR_EC_AS, HR_EC_AS)

common_EC_AS_expr <- common_EC_AS %>%
  dplyr::filter(sig_type_met == "expr", sig_type_tcga == "expr") %>%
  dplyr::mutate(
    shape_EC_AS  = "expr",
    logHR_EC_AS  = meta_pool(logHR_met, se_met, logHR_tcga, se_tcga,
                             HR_met, HR_tcga),
    HR_EC_AS     = exp(logHR_EC_AS)
  ) %>%
  dplyr::select(name, level, shape_EC_AS, logHR_EC_AS, HR_EC_AS)

common_EC_AS_nodes <- dplyr::bind_rows(common_EC_AS_cna, common_EC_AS_expr) %>%
  dplyr::group_by(name, level) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

# Merge node-level information across the two layers
node_sig_info <- dplyr::full_join(
  common_EC_nodes,
  common_EC_AS_nodes,
  by = c("name", "level")
)

# ------------------------------------------------------------
# Compute the union of significant genes across layers (gene-level only)
# ------------------------------------------------------------
sig_gene_union <- union(
  common_EC_nodes    %>% dplyr::filter(level == "gene") %>% dplyr::pull(name),
  common_EC_AS_nodes %>% dplyr::filter(level == "gene") %>% dplyr::pull(name)
)

# ------------------------------------------------------------
# Edges: recolor edges for "prognostic-associated but not independent" genes
# ------------------------------------------------------------

# These are prognostic-associated genes (appear in sig_gene_union)
# but NOT among the 5 independent signature genes
assoc_not_indep_genes <- setdiff(sig_gene_union, sig_genes_all)

cat(">> Prognostic-associated but NOT independent genes found:", 
    length(assoc_not_indep_genes), "\n")
print(sort(assoc_not_indep_genes))

# Update edge_type for these 4 genes (or however many)
edges_net <- edges_net %>%
  dplyr::mutate(
    edge_type = dplyr::case_when(
      to %in% highlight_genes        ~ "signature",          # keep your original highlight
      to %in% assoc_not_indep_genes  ~ "assoc_not_indep",    # NEW red edges
      TRUE                           ~ "other"
    ),
    edge_alpha = dplyr::case_when(
      edge_type == "signature"       ~ 1,
      edge_type == "assoc_not_indep" ~ 1,
      TRUE                           ~ 0.6
    )
  )

# ============================================================
# HYPERGEOMETRIC TEST PER PATHWAY (NO BH CORRECTION)
# Background = genes that appear in the drawn network
# ============================================================

# 1) Background genes = all unique genes in the network
genes_in_network <- unique(path_genes_net_filt2$Gene_Symbol)
N_bg <- length(genes_in_network)

# 2) "Success" genes in background = prognostic genes that also exist in network
prog_genes_in_network <- intersect(sig_gene_union, genes_in_network)
M_success <- length(prog_genes_in_network)

cat("=== Hypergeometric background ===\n")
cat("Background genes (N):", N_bg, "\n")
cat("Prognostic genes in background (M):", M_success, "\n\n")

# Helper: enrichment p-value = P(X >= k)
# X ~ Hypergeometric(N_bg, M_success, n_draw)
hypergeom_p_enrich <- function(k, n_draw, M, N) {
  # If no successes possible or nothing drawn, return NA
  if (is.na(k) || is.na(n_draw) || is.na(M) || is.na(N)) return(NA_real_)
  if (N <= 0 || n_draw <= 0) return(NA_real_)
  if (M < 0 || M > N) return(NA_real_)
  if (k < 0) return(NA_real_)
  if (k > n_draw) return(0)
  1 - phyper(q = k - 1, m = M, n = N - M, k = n_draw)
}

# 3) Per-pathway counts + p-values
per_pathway_hyper <- path_genes_net_filt2 %>%
  dplyr::distinct(pathway, Gene_Symbol) %>%
  dplyr::mutate(is_prognostic = Gene_Symbol %in% prog_genes_in_network) %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(
    n_genes_pathway      = dplyr::n_distinct(Gene_Symbol),  # n_draw
    n_prog_genes_pathway = sum(is_prognostic, na.rm = TRUE),# k
    prog_genes_list      = paste(sort(unique(Gene_Symbol[is_prognostic])), collapse = ", "),
    p_hyper_enrich       = hypergeom_p_enrich(
      k      = n_prog_genes_pathway,
      n_draw = n_genes_pathway,
      M      = M_success,
      N      = N_bg
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(p_hyper_enrich, dplyr::desc(n_prog_genes_pathway), pathway)

cat("=== Per-pathway hypergeometric enrichment (NO BH) ===\n")
print(per_pathway_hyper, n = Inf)

cat("\n=== Nominal p < 0.05 (NO BH) ===\n")
print(dplyr::filter(per_pathway_hyper, p_hyper_enrich < 0.05), n = Inf)

# Stars for figure labels
p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

pathway_labels_with_stars <- per_pathway_hyper %>%
  dplyr::mutate(
    star = vapply(p_hyper_enrich, p_to_stars, character(1)),
    plot_label_star = ifelse(
      star == "",
      paste0(pathway, " (", n_prog_genes_pathway, " out of ", n_genes_pathway, ")"),
      paste0(pathway, " (", n_prog_genes_pathway, " out of ", n_genes_pathway, ") ", star)
    )
  ) %>%
  dplyr::select(pathway, plot_label_star)

# ------------------------------------------------------------
# Node table (pathways + genes) with all annotations
# ------------------------------------------------------------

# Create pathway nodes and attach star labels (from hypergeometric p-values)
nodes_pathways <- tibble::tibble(
  name = unique(path_genes_net_filt2$pathway) %>% stringr::str_trim(),
  type = "pathway"
) %>%
  dplyr::left_join(
    pathway_labels_with_stars,
    by = c("name" = "pathway")
  ) %>%
  dplyr::mutate(
    plot_label = plot_label_star
  )

# Create gene nodes and attach expression-direction summary and signature membership
nodes_genes <- tibble::tibble(
  name = unique(path_genes_net_filt2$Gene_Symbol) %>% stringr::str_trim(),
  type = "gene"
) %>%
  dplyr::left_join(gene_expr_summary, by = c("name" = "Gene_Symbol")) %>%
  dplyr::mutate(
    is_signature = name %in% sig_genes_all
  )

# Combine nodes and attach pooled HR information, shapes, and plot labels
nodes_net <- dplyr::bind_rows(nodes_pathways, nodes_genes) %>%
  dplyr::left_join(
    node_sig_info,
    by = c("name" = "name", "type" = "level")
  ) %>%
  dplyr::mutate(
    expr_group = dplyr::case_when(
      type == "pathway"         ~ "pathway",
      expr_dir_gene == "up"     ~ "gene_up",
      expr_dir_gene == "down"   ~ "gene_down",
      expr_dir_gene == "mixed"  ~ "gene_mixed",
      TRUE                      ~ "gene_unknown"
    ),
    shape_EC_factor    = factor(shape_EC,    levels = c("cna", "expr")),
    shape_EC_AS_factor = factor(shape_EC_AS, levels = c("cna", "expr")),
    plot_label = dplyr::case_when(
      type == "pathway" & !is.na(plot_label) ~ plot_label,
      TRUE                                   ~ name
    )
  )

# ------------------------------------------------------------
# Layout: pathways on a circle, genes in the center
# ------------------------------------------------------------

# Position pathway nodes on a circle
pathway_names <- nodes_net$name[nodes_net$type == "pathway"]
n_p   <- length(pathway_names)
theta <- seq(0, 2*pi, length.out = n_p + 1)[1:n_p]

df_path <- tibble::tibble(
  name = pathway_names,
  x = 1 * cos(theta),
  y = 1 * sin(theta)
)

# Position gene nodes near the origin with a fixed random seed for reproducibility
genes <- nodes_net$name[nodes_net$type == "gene"]
set.seed(123)
df_genes <- tibble::tibble(
  name = genes,
  x = rnorm(length(genes), 0, 0.25),
  y = rnorm(length(genes), 0, 0.25)
)

# Combine node coordinates into a single layout table
nodes_pos <- dplyr::bind_rows(df_path, df_genes)

# Attach coordinates to the node annotation table
nodes_net_layout <- nodes_net %>%
  dplyr::left_join(nodes_pos, by = "name")

# Build an igraph object from pathway–gene edges and node attributes
g_net <- igraph::graph_from_data_frame(
  d = edges_net,
  vertices = nodes_net_layout,
  directed = FALSE
)

nodes_df <- nodes_net_layout

# ------------------------------------------------------------
# Build the network plot (no outer box)
# ------------------------------------------------------------

p_net <- ggraph(
  g_net,
  layout = "manual",
  x = x,
  y = y
) +
  ## Edges
  geom_edge_link(
    aes(colour = edge_type, alpha = edge_alpha),
    linewidth = 0.4
  ) +
  scale_edge_alpha_identity() +
  scale_edge_colour_manual(
    values = c(
      other           = "grey80",
      signature       = "blue",
      assoc_not_indep = "red"
    ),
    guide = "none"
  ) +
  
  ## Background nodes
  geom_node_point(
    data = nodes_df,
    aes(x = x, y = y),
    shape  = 21,
    size   = 3,
    fill   = "white",
    colour = "grey60",
    stroke = 0.4
  ) +
  
  ## ===== RED LAYER: Expr + CNA (E_C_filter) =====
ggnewscale::new_scale_fill() +
  geom_node_point(
    data = nodes_df %>% dplyr::filter(!is.na(shape_EC_factor)),
    aes(
      x     = x,
      y     = y,
      fill  = HR_EC,
      shape = shape_EC_factor
    ),
    colour = "red",
    size   = 4,
    stroke = 0.7
  ) +
  scale_shape_manual(
    name   = "Effect type",
    values = c(cna = 21, expr = 22),
    labels = c(cna = "CNA", expr = "Expression"),
    limits = c("cna", "expr"),
    drop   = FALSE,
    guide  = ggplot2::guide_legend(
      override.aes = list(
        shape  = c(21, 22),
        colour = "black",
        fill   = "white",
        size   = 4,
        stroke = 0.9
      )
    )
  ) +
  scale_fill_gradient2(
    name     = "HR (Expr + CNA)",
    low      = "#FFECEC",
    mid      = "#FFB2B2",
    high     = "red",
    midpoint = 1,
    na.value = "transparent",
    labels   = scales::number_format(accuracy = 0.01)
  ) +
  
  ## ===== BLUE LAYER: Expr + CNA + Age + Stage + PAM50 =====
ggnewscale::new_scale_fill() +
  geom_node_point(
    data = nodes_df %>%
      dplyr::filter(!is.na(HR_EC_AS), !is.na(shape_EC_AS_factor)),
    aes(
      x     = x,
      y     = y,
      fill  = HR_EC_AS,
      shape = shape_EC_AS_factor
    ),
    colour      = "blue",
    size        = 3.5,
    stroke      = 0.7,
    show.legend = TRUE
  ) +
  scale_fill_gradient2(
    name     = "HR (Expr + CNA + Age + Stage + PAM50)",
    low      = "#EDF2FF",
    mid      = "#BFD6FF",
    high     = "blue",
    midpoint = 1,
    na.value = "transparent",
    labels   = scales::number_format(accuracy = 0.01)
  ) +
  
  ## Node labels (pathways now include stars)
  ggrepel::geom_text_repel(
    data = nodes_df,
    aes(x = x, y = y, label = plot_label),
    size          = 3,
    box.padding   = 0.4,
    segment.alpha = 0.4,
    max.overlaps  = Inf
  ) +
  
  theme_void() +
  ggtitle("Network plot with prognostic genes and pathways\n") +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "right",
    legend.title    = element_text(size = 9,  face = "bold"),
    legend.text     = element_text(size = 8,  face = "bold"),
    panel.border    = element_blank(),
    plot.background = element_blank(),
    panel.background= element_blank()
  )

## ------------------------------------------------------------
## Pathway p-value labels (drawn separately in RED on the figure)
## ------------------------------------------------------------

pathway_p_df <- per_pathway_hyper %>%
  dplyr::transmute(
    name    = stringr::str_trim(pathway),
    p_value = p_hyper_enrich,
    p_label = ifelse(is.na(p_value), "", paste0("p=", formatC(p_value, format = "e", digits = 2)))
  ) %>%
  dplyr::left_join(
    nodes_df %>% dplyr::select(name, type, x, y),
    by = "name"
  ) %>%
  dplyr::filter(type == "pathway", p_label != "")


# Output file path
out_file <- "/Volumes/Expansion/Prognostic_Analysis/Results/Hypergeometric test results.xlsx"

# Write to Excel
openxlsx::write.xlsx(
  per_pathway_hyper,
  file = out_file,
  sheetName = "Hypergeometric_Test",
  overwrite = TRUE
)

cat(">> Hypergeometric results written to:\n", out_file, "\n")


# ------------------------------------------------------------
# Save the figure
# ------------------------------------------------------------
ggplot2::ggsave(
  filename = network_file,
  plot     = p_net,
  width    = 11,
  height   = 6,
  dpi      = 300
)

cat(">> Network figure saved to:\n   ", network_file, "\n")