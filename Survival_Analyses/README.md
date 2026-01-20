# Survival Analyses (Cox Proportional Hazards Models)

This directory contains R scripts implementing Cox proportional hazards (CoxPH) survival analyses for gene-level and pathway-level features in the METABRIC (training) and TCGA-BRCA (validation) cohorts. These analyses quantify prognostic effects using hazard ratios (HRs) and form the statistical foundation for downstream Kaplan–Meier stratification, meta-analysis, and signature construction.

All analyses use **overall survival (OS)** as the endpoint.

---

## METABRIC_Survival_Analysis_Gene_OS.R

This script performs gene-level CoxPH survival analyses in the METABRIC training cohort to identify genes independently associated with OS.

### Key steps:
- Loads preprocessed METABRIC gene expression, CNA, and clinical survival data
- Fits CoxPH regression models for each candidate gene using:
  - Expression only
  - CNA only
  - Expression + CNA
  - Expression × CNA interaction
- Evaluates models:
  - Without clinical covariates
  - With adjustment for age, tumor stage, and PAM50 subtype
- Extracts HRs, 95% confidence intervals (CIs), and Wald test p-values
- Identifies independently prognostic genes based on statistical significance and effect size
- Outputs structured result tables for downstream visualization (volcano plots, forest plots, and KM analyses)

This script defines the gene-level prognostic landscape in the training cohort.

---

## METABRIC_Survival_Analysis_Pathway_OS.R

This script conducts pathway-level CoxPH survival analyses in METABRIC using pathway activity scores derived from paclitaxel resistance–associated pathways.

### Key steps:
- Loads pathway-level activity scores derived from:
  - Gene expression
  - Copy number alterations (CNAs)
- Fits CoxPH models for each pathway using:
  - Expression-based pathway scores
  - CNA-based pathway scores
- Adjusts all models for clinical covariates:
  - Age  
  - Tumor stage  
  - PAM50 subtype  
- Estimates HRs, 95% CIs, and p-values for each pathway
- Identifies pathways that remain independently prognostic after clinical adjustment
- Produces pathway-level survival result tables for downstream volcano plots and meta-analysis

This script establishes pathway-level prognostic relevance in the training cohort.

---

## TCGA-BRCA_Survival_Analysis_Gene_OS.R

This script validates gene-level CoxPH survival associations in the independent TCGA-BRCA cohort using the same modeling framework applied in METABRIC.

### Key steps:
- Loads preprocessed TCGA-BRCA gene expression, CNA, and OS data
- Fits CoxPH models for each candidate gene using identical model specifications:
  - Expression only
  - CNA only
  - Expression + CNA
  - Expression × CNA interaction
- Adjusts models for:
  - Age  
  - Tumor stage  
  - PAM50 subtype  
- Extracts HRs, 95% CIs, and statistical significance metrics
- Enables direct cross-cohort comparison of gene-level prognostic effects
- Outputs results for meta-analysis and validation figures

This script provides independent validation of gene-level prognostic signals.

---

## TCGA-BRCA_Survival_Analysis_Pathway_OS.R

This script evaluates pathway-level CoxPH survival associations in TCGA-BRCA, serving as an independent validation of pathway-level findings from METABRIC.

### Key steps:
- Loads TCGA-BRCA pathway activity scores derived from gene expression and CNAs
- Fits CoxPH models using the same pathway definitions as in METABRIC
- Adjusts all models for age, tumor stage, and PAM50 subtype
- Estimates HRs and 95% CIs for each pathway
- Assesses cross-cohort consistency of pathway-level prognostic effects
- Generates pathway-level result tables for meta-analysis and figure generation

This script confirms which resistance-associated pathways retain prognostic relevance in an independent cohort.

---

## Notes
- All CoxPH models use **overall survival (OS)** as the outcome.
- Clinical covariates are harmonized across cohorts.
- Results from these scripts are used directly in:
  - Volcano plots
  - Forest plots
  - Kaplan–Meier stratification
  - Meta-analyses
  - Prognostic signature construction
