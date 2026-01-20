# Load required packages quietly
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(cowplot))

# Set global options to avoid automatic factor conversion
options(stringsAsFactors = FALSE)

# Define input/output directories and create the figure output folder if needed
metabric_results_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/METABRIC"
tcga_results_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Results/TCGA-BRCA"

figures_dir <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"
if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
}

cat(">> Figures will be saved in:\n   ", figures_dir, "\n\n")

# Convert selected columns to numeric if they exist in the input table
to_numeric_if_exists <- function(df, cols) {
  common <- intersect(cols, colnames(df))
  if (length(common) == 0L) return(df)
  df %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(common),
      ~ suppressWarnings(as.numeric(.x))
    ))
}

# Standardize the identifier column to a common name (ID) for gene or pathway tables
standardize_id_column <- function(df, prefer = c("gene", "pathway")) {
  prefer <- match.arg(prefer)
  cn <- colnames(df)
  cn_lower <- tolower(cn)
  
  if (prefer == "gene") {
    targets <- c("gene", "genes", "gene_symbol", "symbol", "hgnc_symbol")
  } else {
    targets <- c("pathway", "term", "pathway_name", "kegg_pathway")
  }
  
  idx <- which(cn_lower %in% targets)
  if (length(idx) == 0L) {
    stop(
      "Could not find ID column for prefer='", prefer,
      "'. Columns are: ", paste(cn, collapse = ", ")
    )
  }
  id_col <- cn[idx[1]]
  df %>% dplyr::rename(ID = !!id_col)
}

# Define colors for effect types (Expression vs CNA)
effect_colors <- c(
  "Expression" = "blue",
  "CNA"        = "red"
)

# Define a shared ggplot theme for all panels
theme_base <- ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid        = ggplot2::element_blank(),
    axis.line         = ggplot2::element_line(color = "black"),
    axis.ticks        = ggplot2::element_line(color = "black"),
    plot.title        = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.title.x      = ggplot2::element_text(face = "bold"),
    axis.title.y      = ggplot2::element_text(face = "bold"),
    axis.text         = ggplot2::element_text(color = "black"),
    strip.text        = ggplot2::element_text(face = "bold"),
    legend.title      = ggplot2::element_text(face = "bold")
  )

# Build a volcano plot using the significant effect (Expression or CNA) per ID
make_volcano <- function(df_filter,
                         cohort_label,
                         id_type = c("Gene","Pathway")) {
  
  id_type <- match.arg(id_type)
  
  df <- df_filter %>%
    to_numeric_if_exists(c(
      "HR_expr","CI_lo_expr","CI_hi_expr","p_expr",
      "HR_cna","CI_lo_cna","CI_hi_cna","p_cna"
    )) %>%
    dplyr::mutate(
      expr_sig = !is.na(CI_lo_expr) & CI_lo_expr > 1,
      cna_sig  = !is.na(CI_lo_cna)  & CI_lo_cna  > 1,
      effect_type = dplyr::case_when(
        expr_sig & !cna_sig ~ "Expression",
        !expr_sig & cna_sig ~ "CNA",
        expr_sig & cna_sig  ~ "Expression",  ## treat “both” as expression
        TRUE                ~ "None"
      ),
      plot_HR = dplyr::case_when(
        effect_type == "Expression" ~ HR_expr,
        effect_type == "CNA"        ~ HR_cna,
        TRUE                        ~ NA_real_
      ),
      plot_p = dplyr::case_when(
        effect_type == "Expression" ~ p_expr,
        effect_type == "CNA"        ~ p_cna,
        TRUE                        ~ NA_real_
      )
    ) %>%
    dplyr::filter(
      effect_type != "None",
      !is.na(plot_HR), plot_HR > 0,
      !is.na(plot_p),  plot_p  > 0
    ) %>%
    dplyr::mutate(
      neglog10p = -log10(plot_p),
      Effect    = factor(effect_type,
                         levels = c("Expression","CNA"))
    )
  
  if (nrow(df) == 0L) {
    warning("No significant entries for volcano in ", cohort_label, " (", id_type, ")")
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_size = 13) +
        ggplot2::ggtitle(paste(id_type, "volcano plot"))
    )
  }
  
  # Define the p-value threshold line for plotting
  thr_p <- -log10(0.05)
  
  # Define x-axis breaks and force inclusion of HR=1
  hr_range   <- range(df$plot_HR, na.rm = TRUE)
  base_break <- scales::pretty_breaks(n = 4)(hr_range)
  x_breaks   <- sort(unique(c(base_break, 1)))
  
  # Define y-axis breaks and force inclusion of the p-value threshold
  y_range  <- range(df$neglog10p, thr_p, na.rm = TRUE)
  y_breaks <- sort(unique(c(scales::pretty_breaks(n = 4)(y_range), thr_p)))
  
  ggplot2::ggplot(df, ggplot2::aes(x = plot_HR,
                                   y = neglog10p,
                                   color = Effect,
                                   label = ID)) +
    ggplot2::geom_hline(yintercept = thr_p,
                        linetype = "dashed",
                        size = 0.4) +
    ggplot2::geom_vline(xintercept = 1,
                        linetype = "dashed",
                        size = 0.4) +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggrepel::geom_text_repel(
      size         = 3.8,
      fontface     = "bold",
      max.overlaps = Inf,
      box.padding  = 0.3,
      show.legend  = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = function(v) sprintf("%.2f", v)
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = function(v) sprintf("%.2f", v)
    ) +
    ggplot2::scale_color_manual(
      name   = "Effect type ",
      values = effect_colors,
      breaks = c("Expression","CNA"),
      labels = c("Expression","CNA")
    ) +
    ggplot2::labs(
      x = "HR",
      y = "-log10(p)",
      title = paste0("Volcano plot (", cohort_label, ", ", id_type, ")")
    ) +
    theme_base +
    ggplot2::theme(
      axis.title.y = ggplot2::element_text(face = "bold", size = 12)
    )
}

# Extract IDs that are significant by expression or CNA (CI_lo > 1 criterion)
get_significant_ids <- function(df_filter) {
  df_filter %>%
    to_numeric_if_exists(c("CI_lo_expr","CI_lo_cna")) %>%
    dplyr::filter(
      (!is.na(CI_lo_expr) & CI_lo_expr > 1) |
        (!is.na(CI_lo_cna) & CI_lo_cna > 1)
    ) %>%
    dplyr::pull(ID) %>%
    unique()
}

# Convert a wide table into a long forest-plot table for Expression and CNA effects
prepare_forest_long <- function(df_filter, cohort_label) {
  df_num <- df_filter %>%
    to_numeric_if_exists(c(
      "HR_expr","CI_lo_expr","CI_hi_expr","p_expr",
      "HR_cna","CI_lo_cna","CI_hi_cna","p_cna"
    ))
  
  expr_long <- df_num %>%
    dplyr::select(ID, HR_expr, CI_lo_expr, CI_hi_expr, p_expr) %>%
    dplyr::filter(!is.na(HR_expr),
                  !is.na(CI_lo_expr),
                  !is.na(CI_hi_expr),
                  HR_expr > 0) %>%
    dplyr::mutate(
      Effect = "Expression",
      HR     = HR_expr,
      CI_lo  = CI_lo_expr,
      CI_hi  = CI_hi_expr,
      p      = p_expr
    ) %>%
    dplyr::select(ID, Effect, HR, CI_lo, CI_hi, p)
  
  cna_long <- df_num %>%
    dplyr::select(ID, HR_cna, CI_lo_cna, CI_hi_cna, p_cna) %>%
    dplyr::filter(!is.na(HR_cna),
                  !is.na(CI_lo_cna),
                  !is.na(CI_hi_cna),
                  HR_cna > 0) %>%
    dplyr::mutate(
      Effect = "CNA",
      HR     = HR_cna,
      CI_lo  = CI_lo_cna,
      CI_hi  = CI_hi_cna,
      p      = p_cna
    ) %>%
    dplyr::select(ID, Effect, HR, CI_lo, CI_hi, p)
  
  dplyr::bind_rows(expr_long, cna_long) %>%
    dplyr::mutate(
      cohort = cohort_label
    )
}

# Build a forest plot comparing two cohorts plus fixed-effect meta-analysis
make_forest_two_cohorts <- function(df_met_long,
                                    df_tcga_long,
                                    common_ids,
                                    id_type = c("Gene","Pathway")) {
  id_type <- match.arg(id_type)
  
  df_all <- dplyr::bind_rows(df_met_long, df_tcga_long) %>%
    dplyr::filter(ID %in% common_ids)
  
  if (nrow(df_all) == 0L) {
    warning("No entries for forest plot (", id_type, ") after intersection.")
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_size = 13) +
        ggplot2::ggtitle(paste(id_type, "forest (no overlap)"))
    )
  }
  
  df_all <- df_all %>%
    dplyr::filter(!is.na(CI_lo), CI_lo > 1, HR > 0)
  
  if (nrow(df_all) == 0L) {
    warning("All entries had CI_lo <= 1 for ", id_type, " forest.")
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_size = 13) +
        ggplot2::ggtitle(paste(id_type, "forest (no CI_lo > 1.0)"))
    )
  }
  
  df_all <- df_all %>%
    dplyr::group_by(ID, Effect) %>%
    dplyr::mutate(n_cohort = dplyr::n_distinct(cohort)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_cohort == 2)
  
  if (nrow(df_all) == 0L) {
    warning("No (ID,Effect) pairs present in both cohorts for ", id_type, " forest.")
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_size = 13) +
        ggplot2::ggtitle(paste0("Forest plot", " (" ,id_type, ")"))
    )
  }
  
  df_all <- df_all %>%
    dplyr::mutate(
      logHR = log(HR),
      se    = (log(CI_hi) - log(CI_lo)) / (2 * 1.96)
    )
  
  meta_df <- df_all %>%
    dplyr::group_by(ID, Effect) %>%
    dplyr::summarise(
      n_cohort = dplyr::n_distinct(cohort),
      logHR_meta = {
        w <- 1 / (se^2)
        sum(w * logHR, na.rm = TRUE) / sum(w, na.rm = TRUE)
      },
      se_meta = sqrt(1 / sum(1 / (se^2), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      HR    = exp(logHR_meta),
      CI_lo = exp(logHR_meta - 1.96 * se_meta),
      CI_hi = exp(logHR_meta + 1.96 * se_meta),
      z     = logHR_meta / se_meta,
      p     = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
      cohort = "Meta Analysis"
    ) %>%
    dplyr::select(ID, Effect, HR, CI_lo, CI_hi, p, cohort)
  
  plot_df <- df_all %>%
    dplyr::select(ID, Effect, HR, CI_lo, CI_hi, p, cohort) %>%
    dplyr::bind_rows(meta_df) %>%
    dplyr::mutate(
      cohort = factor(cohort,
                      levels = c("METABRIC", "TCGA-BRCA", "Meta Analysis")),
      Effect = factor(Effect, levels = c("Expression", "CNA")),
      ID_Effect = paste0(ID, " (", Effect, ")")
    )
  
  plot_df <- plot_df %>%
    dplyr::arrange(ID, Effect, cohort) %>%
    dplyr::mutate(
      ID_Effect = factor(ID_Effect, levels = rev(unique(ID_Effect)))
    )
  
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = HR, y = ID_Effect, xmin = CI_lo, xmax = CI_hi)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",
                        color = "grey40", size = 0.4) +
    ggplot2::geom_errorbarh(height = 0.25, color = "grey30") +
    ggplot2::geom_point(
      ggplot2::aes(color = Effect),
      size = 2.8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", HR)),
      vjust = -1.0,
      hjust = 0.5,
      size  = 3.4,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_x_log10(
      labels = function(v) sprintf("%.2f", v)
    ) +
    ggplot2::scale_color_manual(
      name   = "Effect type",
      values = effect_colors,
      breaks = c("Expression","CNA"),
      labels = c("Expression","CNA")
    ) +
    ggplot2::facet_wrap(~ cohort, ncol = 3) +
    ggplot2::labs(
      x = "HR",
      y = NULL,
      title = paste0("Forest plot", " (" ,id_type, ")")
    ) +
    theme_base +
    ggplot2::theme(
      strip.text   = ggplot2::element_text(face = "bold"),
      axis.text.y  = ggplot2::element_text(face = "bold"),
      
      # >>> ADD THESE <<<
      panel.spacing.x = grid::unit(1.2, "cm"),   # space between facet columns
      panel.spacing.y = grid::unit(0.6, "cm"),   # (optional) vertical spacing
      plot.margin     = ggplot2::margin(t = 5, r = 15, b = 5, l = 5),
    )
}

# Define input Excel files for gene-level results and load the specified worksheet
met_gene_file <- file.path(
  metabric_results_dir,
  "METABRIC_OS_all_Cox_all_models_by_gene_with_LRT.xlsx"
)
tcga_gene_file <- file.path(
  tcga_results_dir,
  "TCGA_BRCA_OS_all_Cox_all_models_by_gene_with_LRT.xlsx"
)

cat(">> [Gene-level] Reading E_C_Age_Stg_PAM50_filter sheets...\n")

met_gene_raw <- readxl::read_excel(
  met_gene_file,
  sheet = "E_C_Age_Stg_PAM50_filter"
)
tcga_gene_raw <- readxl::read_excel(
  tcga_gene_file,
  sheet = "E_C_Age_Stg_PAM50_filter"
)

# Harmonize gene identifier columns across cohorts
met_gene_filter <- met_gene_raw %>%
  standardize_id_column(prefer = "gene")
tcga_gene_filter <- tcga_gene_raw %>%
  standardize_id_column(prefer = "gene")

# Create gene-level volcano plots for METABRIC and TCGA-BRCA
p_gene_vol_met  <- make_volcano(
  met_gene_filter,
  cohort_label = "METABRIC",
  id_type      = "Gene"
)
p_gene_vol_tcga <- make_volcano(
  tcga_gene_filter,
  cohort_label = "TCGA-BRCA",
  id_type      = "Gene"
)

# Remove redundant axis labels in the top row panels
p_gene_vol_met  <- p_gene_vol_met  + ggplot2::labs(x = NULL)
p_gene_vol_tcga <- p_gene_vol_tcga + ggplot2::labs(x = NULL, y = NULL)

# Identify significant genes per cohort and compute intersection
sig_genes_met  <- get_significant_ids(met_gene_filter)
sig_genes_tcga <- get_significant_ids(tcga_gene_filter)
common_genes   <- intersect(sig_genes_met, sig_genes_tcga)

cat("   Significant genes (METABRIC): ", length(sig_genes_met),  "\n")
cat("   Significant genes (TCGA):     ", length(sig_genes_tcga), "\n")
cat("   Overlapping genes:            ", length(common_genes),   "\n\n")

# Prepare long-format tables for forest plotting and build gene forest plot
met_gene_long  <- prepare_forest_long(met_gene_filter, cohort_label = "METABRIC")
tcga_gene_long <- prepare_forest_long(tcga_gene_filter, cohort_label = "TCGA-BRCA")

p_gene_forest <- make_forest_two_cohorts(
  df_met_long  = met_gene_long,
  df_tcga_long = tcga_gene_long,
  common_ids   = common_genes,
  id_type      = "Gene"
)

# Remove redundant x-axis label in the gene forest panel
p_gene_forest <- p_gene_forest + ggplot2::labs(x = NULL)

# Define input Excel files for pathway-level results and load the specified worksheet
met_path_file <- file.path(
  metabric_results_dir,
  "METABRIC_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx"
)
tcga_path_file <- file.path(
  tcga_results_dir,
  "TCGA_BRCA_OS_all_Cox_all_models_by_pathway_with_LRT.xlsx"
)

cat(">> [Pathway-level] Reading E_C_Age_Stg_PAM50_filter sheets...\n")

met_path_raw <- readxl::read_excel(
  met_path_file,
  sheet = "E_C_Age_Stg_PAM50_filter"
)
tcga_path_raw <- readxl::read_excel(
  tcga_path_file,
  sheet = "E_C_Age_Stg_PAM50_filter"
)

# Harmonize pathway identifier columns across cohorts
met_path_filter <- met_path_raw %>%
  standardize_id_column(prefer = "pathway")
tcga_path_filter <- tcga_path_raw %>%
  standardize_id_column(prefer = "pathway")

# Create pathway-level volcano plots for METABRIC and TCGA-BRCA
p_path_vol_met  <- make_volcano(
  met_path_filter,
  cohort_label = "METABRIC",
  id_type      = "Pathway"
)
p_path_vol_tcga <- make_volcano(
  tcga_path_filter,
  cohort_label = "TCGA-BRCA",
  id_type      = "Pathway"
)

# Remove redundant y-axis label in the right pathway volcano panel
p_path_vol_tcga <- p_path_vol_tcga + ggplot2::labs(y = NULL)

# Identify significant pathways per cohort and compute intersection
sig_path_met  <- get_significant_ids(met_path_filter)
sig_path_tcga <- get_significant_ids(tcga_path_filter)
common_paths  <- intersect(sig_path_met, sig_path_tcga)

cat("   Significant pathways (METABRIC): ", length(sig_path_met),  "\n")
cat("   Significant pathways (TCGA):     ", length(sig_path_tcga), "\n")
cat("   Overlapping pathways:            ", length(common_paths),  "\n\n")

# Prepare long-format tables for forest plotting and build pathway forest plot
met_path_long  <- prepare_forest_long(met_path_filter, cohort_label = "METABRIC")
tcga_path_long <- prepare_forest_long(tcga_path_filter, cohort_label = "TCGA-BRCA")

p_path_forest <- make_forest_two_cohorts(
  df_met_long  = met_path_long,
  df_tcga_long = tcga_path_long,
  common_ids   = common_paths,
  id_type      = "Pathway"
)

# Configure a single shared legend from panel a and disable legends elsewhere
p_gene_vol_met <- p_gene_vol_met + ggplot2::theme(legend.position = "bottom")

p_gene_vol_tcga <- p_gene_vol_tcga + ggplot2::theme(legend.position = "none")
p_gene_forest   <- p_gene_forest   + ggplot2::theme(legend.position = "none")
p_path_vol_met  <- p_path_vol_met  + ggplot2::theme(legend.position = "none")
p_path_vol_tcga <- p_path_vol_tcga + ggplot2::theme(legend.position = "none")
p_path_forest   <- p_path_forest   + ggplot2::theme(legend.position = "none")

# Extract the shared legend from the METABRIC gene volcano panel
legend_shared <- cowplot::get_legend(
  p_gene_vol_met +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title    = ggplot2::element_text(face = "bold"),
      legend.text     = ggplot2::element_text(size = 14),
      plot.margin     = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
    )
)

# Build legend-free panel variants to assemble the 2x3 grid
p_gene_vol_met_noleg  <- p_gene_vol_met  + ggplot2::theme(legend.position = "none")
p_gene_vol_tcga_noleg <- p_gene_vol_tcga + ggplot2::theme(legend.position = "none")
p_gene_forest_noleg   <- p_gene_forest   + ggplot2::theme(legend.position = "none")
p_path_vol_met_noleg  <- p_path_vol_met  + ggplot2::theme(legend.position = "none")
p_path_vol_tcga_noleg <- p_path_vol_tcga + ggplot2::theme(legend.position = "none")
p_path_forest_noleg   <- p_path_forest   + ggplot2::theme(legend.position = "none")

# Assemble the 2x3 grid layout using patchwork
grid_main <-
  (p_gene_vol_met_noleg | p_gene_vol_tcga_noleg | p_gene_forest_noleg) /
  (p_path_vol_met_noleg | p_path_vol_tcga_noleg | p_path_forest_noleg)

# Add panel tags (a–f) in the grid
grid_main <-
  grid_main +
  patchwork::plot_annotation(
    tag_levels = "A",
    tag_suffix = "."
  )

# Style the panel tags and tighten plot margins
grid_main <-
  grid_main &
  ggplot2::theme(
    plot.tag    = ggplot2::element_text(face = "bold", size = 16),
    plot.margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5)
  )

# Combine the grid and the shared legend into one figure with cowplot
combined_plot <- cowplot::plot_grid(
  grid_main,
  legend_shared,
  ncol        = 1,
  rel_heights = c(1, 0.08)
)

# Define the output file path for the final figure
fig_file <- file.path(
  figures_dir,
  "Figure 2.jpeg"
)

# Save the combined figure to disk with specified dimensions and resolution
ggplot2::ggsave(
  filename = fig_file,
  plot     = combined_plot,
  width    = 20,
  height   = 11,
  dpi      = 300
)

cat(">> 2 x 3 grid figure saved to:\n   ", fig_file, "\n\n")
cat(">> Script finished.\n")