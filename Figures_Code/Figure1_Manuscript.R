# Load required packages quietly
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(tibble); library(ggplot2)
  library(cowplot)
})

# Define a null-coalescing helper for safer defaults
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Define input/output directories and create output folder if needed
metabric_dir <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"
tcga_dir     <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
out_dir      <- "/Volumes/Expansion/Prognostic_Analysis/Results/Figures_Manuscript"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Set output figure path and print it
out_file <- file.path(out_dir, "Figure 1.jpeg")
cat(">>> Output:", out_file, "\n")

# Normalize METABRIC sample IDs to a consistent uppercase trimmed format
norm_met <- function(x) toupper(str_trim(as.character(x)))

# Normalize TCGA barcodes to patient-level identifiers (first 12 characters)
norm_tcga_patient <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  substr(x, 1, 12)
}

# Extract sample IDs from wide tables where genes are rows and samples are columns
get_wide_sample_ids <- function(df, gene_col = "Gene_Symbol") {
  setdiff(colnames(df), gene_col)
}

# Define a shared ggplot theme for boxed panels
theme_fig1 <- theme_minimal(base_size = 13) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.background = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# Set manual colors for cohorts
cohort_cols <- c("METABRIC" = "red", "TCGA-BRCA" = "blue")

# Load METABRIC expression, CNA, and clinical data
cat(">>> Loading METABRIC...\n")

met_expr_file <- file.path(metabric_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
met_cna_file  <- file.path(metabric_dir, "data_cna.txt")
met_clin_file <- file.path(metabric_dir, "brca_metabric_clinical_data.tsv")

# Read METABRIC expression data and standardize gene column naming
met_expr <- read_tsv(met_expr_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol) %>%
  dplyr::select(-Entrez_Gene_Id)

# Read METABRIC CNA data and standardize gene column naming
met_cna <- read_tsv(met_cna_file, show_col_types = FALSE) %>%
  dplyr::rename(Gene_Symbol = Hugo_Symbol) %>%
  dplyr::select(-Entrez_Gene_Id)

# Read and curate METABRIC clinical data, derive OS time/event and simplify stage labels
met_clin <- read_tsv(met_clin_file, show_col_types = FALSE) %>%
  dplyr::rename(
    sample    = `Sample ID`,
    Age       = `Age at Diagnosis`,
    Stage_raw = `Tumor Stage`,
    OS_TM     = `Overall Survival (Months)`,
    OS_ST     = `Overall Survival Status`,
    PAM50_raw = `Pam50 + Claudin-low subtype`
  ) %>%
  mutate(
    sample   = norm_met(sample),
    OS_TIME  = suppressWarnings(as.numeric(OS_TM)),
    OS_EVENT = ifelse(str_detect(OS_ST %||% "", "^1"), 1, 0),
    Age      = suppressWarnings(as.numeric(Age)),
    Stage_simple = case_when(
      Stage_raw %in% c("0","Stage 0","Stage0")  ~ "Stage0",
      Stage_raw %in% c("1","I","Stage I")       ~ "StageI",
      Stage_raw %in% c("2","II","Stage II")     ~ "StageII",
      Stage_raw %in% c("3","III","Stage III")   ~ "StageIII",
      Stage_raw %in% c("4","IV","Stage IV")     ~ "StageIV",
      is.na(Stage_raw) | Stage_raw == ""        ~ NA_character_,
      TRUE                                      ~ NA_character_
    ),
    PAM50 = as.character(PAM50_raw)
  ) %>%
  filter(!is.na(OS_TIME), !is.na(OS_EVENT)) %>%
  mutate(
    Stage_simple = factor(
      Stage_simple,
      levels = c("Stage0", "StageI", "StageII", "StageIII", "StageIV")
    )
  ) %>%
  dplyr::select(sample, OS_TIME, OS_EVENT, Age, Stage_simple, PAM50)

# Load TCGA-BRCA expression, CNA, and survival/clinical covariates
cat(">>> Loading TCGA-BRCA...\n")

tcga_expr_file <- file.path(tcga_dir, "TCGA-BRCA_exprs.z.tsv")
tcga_cna_file  <- file.path(tcga_dir, "TCGA-BRCA_CNA.tsv")
tcga_surv_file <- file.path(tcga_dir, "TCGA_BRCA_survival_with_PAM50_Age_Stage.tsv")

# Read TCGA expression and CNA matrices plus survival table
tcga_expr <- read_tsv(tcga_expr_file, show_col_types = FALSE)
tcga_cna  <- read_tsv(tcga_cna_file,  show_col_types = FALSE)
tcga_surv <- read_tsv(tcga_surv_file, show_col_types = FALSE)

# Build TCGA clinical table, convert OS time to months, and harmonize stage and PAM50 labels
tcga_clin <- tcga_surv %>%
  transmute(
    sample     = as.character(sample),
    patient    = norm_tcga_patient(sample),
    OS_EVENT   = suppressWarnings(as.numeric(OS)),
    OS_TIME_d  = suppressWarnings(as.numeric(OS.time)),
    Age        = suppressWarnings(as.numeric(Age)),
    Stage_raw  = as.character(Stage),
    PAM50_raw  = as.character(BRCA_Subtype_PAM50)
  ) %>%
  mutate(
    OS_TIME = OS_TIME_d / 30.44,
    Stage_simple = case_when(
      is.na(Stage_raw)                                                      ~ NA_character_,
      Stage_raw %in% c("Stage I", "Stage IA", "Stage IB")                    ~ "StageI",
      Stage_raw %in% c("Stage II", "Stage IIA", "Stage IIB")                 ~ "StageII",
      Stage_raw %in% c("Stage III", "Stage IIIA", "Stage IIIB","Stage IIIC") ~ "StageIII",
      Stage_raw %in% c("Stage IV")                                           ~ "StageIV",
      Stage_raw %in% c("Stage X")                                            ~ "StageX",
      TRUE                                                               ~ NA_character_
    ),
    PAM50 = as.character(PAM50_raw),
    Stage_simple = factor(
      Stage_simple,
      levels = c("StageI", "StageII", "StageIII", "StageIV", "StageX")
    )
  ) %>%
  filter(!is.na(OS_TIME), !is.na(OS_EVENT)) %>%
  dplyr::select(sample, patient, OS_TIME, OS_EVENT, Age, Stage_simple, PAM50)

# Compute Expr ∩ CNA ∩ Survival intersections for both cohorts
cat(">>> Computing All3 intersections...\n")

# Get METABRIC sample IDs present in expression, CNA, and survival tables
met_ids_expr <- norm_met(get_wide_sample_ids(met_expr, "Gene_Symbol"))
met_ids_cna  <- norm_met(get_wide_sample_ids(met_cna,  "Gene_Symbol"))
met_ids_surv <- unique(met_clin$sample)

# Keep METABRIC samples available across all three data types
met_all3 <- Reduce(intersect, list(unique(met_ids_expr), unique(met_ids_cna), unique(met_ids_surv)))
met_all3_clin <- met_clin %>% filter(sample %in% met_all3)
N_met_all3 <- nrow(met_all3_clin)

# Get TCGA patient IDs present in expression, CNA, and survival tables
tcga_ids_expr <- norm_tcga_patient(get_wide_sample_ids(tcga_expr, "Gene_Symbol"))
tcga_ids_cna  <- norm_tcga_patient(get_wide_sample_ids(tcga_cna,  "Gene_Symbol"))
tcga_ids_surv <- unique(tcga_clin$patient)

# Keep TCGA patients available across all three data types
tcga_all3_patient <- Reduce(intersect, list(unique(tcga_ids_expr), unique(tcga_ids_cna), unique(tcga_ids_surv)))
tcga_all3_clin <- tcga_clin %>% filter(patient %in% tcga_all3_patient)
N_tcga_all3 <- nrow(tcga_all3_clin)

cat(">>> METABRIC All3 N:", N_met_all3, "\n")
cat(">>> TCGA-BRCA All3 N:", N_tcga_all3, "\n")

# Create plotting-ready clinical variables with explicit NA categories and cohort-specific rules
met_plot <- met_all3_clin %>%
  mutate(
    Cohort = "METABRIC",
    Age_plot = Age,
    Stage_plot = case_when(
      is.na(Stage_simple) | Stage_simple == "" ~ "NA",
      Stage_simple == "Stage0"   ~ "Stage 0",
      Stage_simple == "StageI"   ~ "Stage I",
      Stage_simple == "StageII"  ~ "Stage II",
      Stage_simple == "StageIII" ~ "Stage III",
      Stage_simple == "StageIV"  ~ "Stage IV",
      TRUE ~ "NA"
    ),
    PAM50_plot = case_when(
      is.na(PAM50) | PAM50 == "" ~ "NA",
      toupper(PAM50) == "NC" ~ "NA",          # NC -> NA (METABRIC only)
      str_to_lower(PAM50) == "claudin-low" ~ "Claudin-low",
      TRUE ~ PAM50
    )
  ) %>%
  mutate(
    Stage_plot = factor(Stage_plot, levels = c("Stage 0","Stage I","Stage II","Stage III","Stage IV","NA")),
    PAM50_plot = factor(
      PAM50_plot,
      levels = c("LumA","LumB","Basal","Her2","Claudin-low","Normal","NA")
    )
  )

# Create plotting-ready clinical variables for TCGA with Stage X handling and no Claudin-low
tcga_plot <- tcga_all3_clin %>%
  mutate(
    Cohort = "TCGA-BRCA",
    Age_plot = Age,
    Stage_plot = case_when(
      is.na(Stage_simple) | Stage_simple == "" ~ "NA",
      Stage_simple == "StageI"   ~ "Stage I",
      Stage_simple == "StageII"  ~ "Stage II",
      Stage_simple == "StageIII" ~ "Stage III",
      Stage_simple == "StageIV"  ~ "Stage IV",
      Stage_simple == "StageX"   ~ "Stage X",
      TRUE ~ "Stage X"
    ),
    PAM50_plot = case_when(
      is.na(PAM50) | PAM50 == "" ~ "NA",
      toupper(PAM50) == "NC" ~ "NA",
      str_to_lower(PAM50) == "claudin-low" ~ "NA",
      TRUE ~ PAM50
    )
  ) %>%
  mutate(
    Stage_plot = factor(Stage_plot, levels = c("Stage I","Stage II","Stage III","Stage IV","Stage X","NA")),
    PAM50_plot = factor(PAM50_plot, levels = c("LumA","LumB","Basal","Her2","Normal","NA"))
  )

# Compute cohort-specific Age summary statistics for facet labels
age_stats <- bind_rows(met_plot, tcga_plot) %>%
  mutate(Age_use = ifelse(!is.na(Age_plot) & Age_plot >= 20, Age_plot, NA_real_)) %>%
  group_by(Cohort) %>%
  summarise(
    N_all3 = ifelse(Cohort[1] == "METABRIC", N_met_all3, N_tcga_all3),
    mean_age = mean(Age_use, na.rm = TRUE),
    median_age = median(Age_use, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Cohort_label = paste0(
      Cohort, " (N = ", N_all3,
      ", mean = ", round(mean_age, 0),
      ", median = ", round(median_age, 0), ")"
    )
  )

# Map cohort names to Age facet labels containing N/mean/median
facet_labeller_age <- function(x) {
  lbl <- age_stats$Cohort_label[match(x, age_stats$Cohort)]
  ifelse(is.na(lbl), x, lbl)
}

# Map cohort names to facet labels containing only N
facet_labeller_n <- function(x) {
  ifelse(
    x == "METABRIC",
    paste0("METABRIC (N = ", N_met_all3, ")"),
    ifelse(
      x == "TCGA-BRCA",
      paste0("TCGA-BRCA (N = ", N_tcga_all3, ")"),
      x
    )
  )
}

# Combine cohorts into a single plotting table with ordered cohort factor
plot_data_all <- bind_rows(met_plot, tcga_plot) %>%
  mutate(Cohort = factor(Cohort, levels = c("METABRIC","TCGA-BRCA")))

# Subset valid ages for histogram plotting
plot_age <- plot_data_all %>%
  filter(!is.na(Age_plot), Age_plot >= 20)

# Plot Age distributions faceted by cohort
p_age <- ggplot(plot_age, aes(x = Age_plot, fill = Cohort)) +
  geom_histogram(bins = 30, color = NA) +
  facet_wrap(~Cohort, ncol = 2, labeller = as_labeller(facet_labeller_age)) +
  scale_fill_manual(values = cohort_cols) +
  scale_x_continuous(breaks = seq(20, 100, 10), limits = c(20, 100)) +
  labs(title = "Age", x = "Age", y = "Count", fill = "Cohort") +
  theme_fig1 +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", size = 16)
  )

# Tabulate stage counts per cohort for bar plotting
stage_counts_df <- plot_data_all %>%
  mutate(
    Stage_plot2 = as.character(Stage_plot),
    Stage_plot2 = ifelse(is.na(Stage_plot2) | Stage_plot2 == "", "NA", Stage_plot2),
    Stage_plot2 = factor(Stage_plot2, levels = c("Stage 0","Stage I","Stage II","Stage III","Stage IV","Stage X","NA"))
  ) %>%
  count(Cohort, Stage_plot2, name = "n")

# Plot tumor stage distributions faceted by cohort
p_stage <- ggplot(stage_counts_df, aes(x = Stage_plot2, y = n, fill = Cohort)) +
  geom_col() +
  geom_text(aes(label = n), fontface = "bold", vjust = -0.25, size = 4) +
  facet_wrap(~Cohort, ncol = 2, scales = "free_x", labeller = as_labeller(facet_labeller_n)) +
  scale_fill_manual(values = cohort_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    title = "Tumor Stage",
    x = "Stage",
    y = "Count",
    fill = "Cohort"
  ) +
  theme_fig1 +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", size = 16)
  )

# Tabulate PAM50 subtype counts per cohort for bar plotting
pam_counts_df <- plot_data_all %>%
  mutate(
    PAM50_plot2 = as.character(PAM50_plot),
    PAM50_plot2 = ifelse(is.na(PAM50_plot2) | PAM50_plot2 == "", "NA", PAM50_plot2)
  ) %>%
  count(Cohort, PAM50_plot2, name = "n") %>%
  mutate(
    PAM50_plot2 = factor(PAM50_plot2, levels = c("LumA","LumB","Basal","Her2","Claudin-low","Normal","NC","NA"))
  )

# Plot PAM50 subtype distributions faceted by cohort
p_pam50 <- ggplot(pam_counts_df, aes(x = PAM50_plot2, y = n, fill = Cohort)) +
  geom_col() +
  geom_text(aes(label = n), fontface = "bold", vjust = -0.25, size = 4) +
  facet_wrap(~Cohort, ncol = 2, scales = "free_x", labeller = as_labeller(facet_labeller_n)) +
  scale_fill_manual(values = cohort_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    title = "PAM50 Subtype",
    x = "PAM50",
    y = "Count",
    fill = "Cohort"
  ) +
  theme_fig1 +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", size = 16)
  )

# Add panel tags for multi-panel figure labeling
p_age_tag   <- ggdraw(p_age)   + draw_label("A.", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 16)
p_stage_tag <- ggdraw(p_stage) + draw_label("B.", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 16)
p_pam_tag   <- ggdraw(p_pam50) + draw_label("C.", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 16)

# Extract a shared legend to place once at the bottom of the combined figure
legend_shared <- cowplot::get_legend(
  (p_age +
     theme(
       legend.position = "bottom",
       legend.title = element_text(face = "bold"),
       legend.text  = element_text(size = 12)
     ))
)

# Assemble the final stacked figure with a shared legend
final_fig <- cowplot::plot_grid(
  cowplot::plot_grid(p_age_tag, p_stage_tag, p_pam_tag, ncol = 1, rel_heights = c(1, 1, 1)),
  legend_shared,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Save the final figure to disk
ggsave(out_file, final_fig, width = 12, height = 14, dpi = 300)
cat(">>> Saved:", out_file, "\n")