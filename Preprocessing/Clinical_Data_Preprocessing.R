# Load packages
suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))


# Define directories
tcga_dir <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/TCGA-BRCA"
metabric_dir <- "/Volumes/Expansion/Prognostic_Analysis/Datasets/METABRIC"


###################################################################################################
###################################### 1. TCGA-BRCA dataset #######################################
###################################################################################################

# Read TCGA-BRCA survival data (OS only will be retained)
tcga_surv_brca <- readr::read_tsv(
  file.path(tcga_dir, "TCGA_BRCA_survival.txt"),
  show_col_types = FALSE
)

# Convert overall survival time from days to years
tcga_surv_brca <- tcga_surv_brca %>%
  mutate(
    OS.time.years = round(as.numeric(OS.time) / 365.25)
  ) %>%
  relocate(OS.time.years, .after = OS.time)

# Retrieve PAM50 molecular subtype annotations
pam50 <- TCGAquery_subtype(tumor = "BRCA")

pam50_subtypes <- pam50 %>%
  dplyr::select(patient, BRCA_Subtype_PAM50)

# Merge OS survival data with PAM50 subtypes
tcga_surv_pam50 <- tcga_surv_brca %>%
  left_join(pam50_subtypes, by = c("_PATIENT" = "patient"))

# Download clinical metadata for TCGA-BRCA
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# Harmonize patient identifier column name
if ("_PATIENT" %in% colnames(tcga_surv_pam50)) {
  tcga_surv_pam50 <- dplyr::rename(
    tcga_surv_pam50,
    PATIENT = `_PATIENT`
  )
}

# Extract age (years) and pathological stage
clinical_age_stage <- clinical %>%
  transmute(
    PATIENT = bcr_patient_barcode,
    Age = case_when(
      "age_at_diagnosis" %in% names(clinical) ~ round(as.numeric(age_at_diagnosis) / 365.25),
      "days_to_birth" %in% names(clinical) ~ round(-as.numeric(days_to_birth) / 365.25),
      TRUE ~ NA_real_
    ),
    Stage = if ("ajcc_pathologic_stage" %in% names(clinical))
      as.character(ajcc_pathologic_stage) else NA_character_
  )

# Merge OS survival data with clinical covariates
tcga_surv_os <- tcga_surv_pam50 %>%
  left_join(clinical_age_stage, by = "PATIENT")

# Save OS-only survival dataset
out_file <- file.path(tcga_dir, "TCGA_BRCA_OS_with_PAM50_Age_Stage.tsv")
write_tsv(tcga_surv_os, out_file)
message("Saved OS-only TCGA survival data to: ", out_file)


###################################################################################################
####################################### 2. METABRIC dataset #######################################
###################################################################################################

# Read METABRIC expression, CNA, and clinical survival datasets
metabric_expr_file <- file.path(metabric_dir, "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt")
metabric_cna_file  <- file.path(metabric_dir, "data_cna.txt")
metabric_survival_file <- file.path(metabric_dir, "data_clinical_patient.txt")

metabric_expr <- read_tsv(metabric_expr_file, show_col_types = FALSE)
metabric_cna  <- read_tsv(metabric_cna_file, show_col_types = FALSE)

# Binarize CNA values
metabric_cna[, -1] <- ifelse(metabric_cna[, -1] != 0, 1, 0)

metabric_survival <- read_tsv(metabric_survival_file, show_col_types = FALSE)
metabric_survival <- metabric_survival[-c(1:4),]

# Rename OS-related columns only
metabric_survival <- metabric_survival %>% 
  dplyr::rename(
    OS_STATUS = `Overall Survival Status`,
    OS_MONTHS = `Overall Survival (Months)`
  )

# Convert OS status to binary event indicator
metabric_survival$OS_STATUS <- ifelse(
  grepl("^0:", metabric_survival$OS_STATUS), 0,
  ifelse(grepl("^1:", metabric_survival$OS_STATUS), 1, NA)
)

# Convert OS time from months to years
metabric_survival <- metabric_survival %>%
  mutate(
    OS_MONTHS = as.numeric(OS_MONTHS),
    OS_YEARS  = as.integer(round(OS_MONTHS / 12))
  )

# Print dataset dimensions
dim(metabric_expr)
print(paste("There are", ncol(metabric_expr), "patients in expression data"))
print(paste("There are", nrow(metabric_expr), "genes in expression data"))

dim(metabric_cna)
print(paste("There are", ncol(metabric_cna), "patients in CNA data"))
print(paste("There are", nrow(metabric_cna), "genes in CNA data"))

dim(metabric_survival)
print(paste("There are", nrow(metabric_survival), "patients in OS survival data"))