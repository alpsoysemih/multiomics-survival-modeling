# Preprocessing

This directory contains preprocessing steps required to run the survival modeling pipeline in **METABRIC** and **TCGA-BRCA**.

> **Important note (TCGA-BRCA molecular preprocessing):**  
> TCGA-BRCA **gene expression** and **CNA** preprocessing can be reproduced by running the scripts available in my GitHub repository:  
> **alpsoysemih/multiomics-drug-response-resistance-tl**.  
> Once TCGA-BRCA expression and CNA matrices are generated in processed form, the remaining steps required for survival modeling—such as **z-score normalization**, **CNA binarization**, and computation of **pathway CNA scores**—are implemented directly inside the downstream analysis scripts in this repository.

---

## Clinical data preprocessing

The script in this folder focuses on preparing **overall survival (OS)** and key **clinical covariates** (age, tumor stage, PAM50 subtype) for both cohorts.

### 1) TCGA-BRCA clinical preprocessing

**Inputs**
- `TCGA_BRCA_survival.txt` (survival table; OS endpoint retained)
- PAM50 subtype annotations via `TCGAquery_subtype(tumor = "BRCA")`
- Clinical metadata via `GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")`

**Steps**
- Read TCGA-BRCA survival data and retain OS variables.
- Convert OS time from **days to years** (`OS.time.years = round(OS.time / 365.25)`).
- Retrieve **PAM50 subtypes** (TCGAquery_subtype) and merge with survival data.
- Download clinical metadata and extract:
  - **Age (years)** (derived from `age_at_diagnosis` if present; otherwise from `days_to_birth`)
  - **Pathological tumor stage** (`ajcc_pathologic_stage`, if available)
- Harmonize patient identifiers (rename `_PATIENT` → `PATIENT`) and merge all covariates into a single table.
- Save the OS-only clinical dataset as:
  - `TCGA_BRCA_OS_with_PAM50_Age_Stage.tsv`

**Output**
- `TCGA_BRCA_OS_with_PAM50_Age_Stage.tsv`

---

### 2) METABRIC clinical preprocessing (plus CNA binarization)

**Inputs**
- Expression: `data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt`
- CNA: `data_cna.txt`
- Clinical survival: `data_clinical_patient.txt`

**Steps**
- Load METABRIC expression and CNA matrices.
- **Binarize CNA values**:
  - Any non-zero GISTIC value → `1` (altered)
  - Copy-neutral (`0`) → `0`
- Load clinical survival table and remove non-data header rows.
- Rename OS columns:
  - `Overall Survival Status` → `OS_STATUS`
  - `Overall Survival (Months)` → `OS_MONTHS`
- Convert OS status to a binary event indicator:
  - `"0:"` → `0` (censored)
  - `"1:"` → `1` (event)
- Convert OS time from **months to years**:
  - `OS_YEARS = round(OS_MONTHS / 12)`
- Print dataset dimensions for sanity checks (expression, CNA, and survival).

**Outputs**
- A clinical OS table with harmonized `OS_STATUS`, `OS_MONTHS`, and `OS_YEARS` fields (kept in-memory for downstream scripts).
- A binarized CNA matrix (altered vs. copy-neutral), used in downstream survival analyses.

---

## Notes

- This preprocessing stage prepares **clinical endpoints and covariates** used throughout CoxPH and KM analyses.
- Molecular preprocessing (especially for TCGA-BRCA expression/CNA matrices) can be reproduced externally as noted above, while key transformations (e.g., z-scoring and CNA binarization used for modeling) are also available within downstream scripts to ensure end-to-end reproducibility.
