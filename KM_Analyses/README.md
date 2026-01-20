# Kaplan–Meier Survival Analyses

This directory contains R scripts implementing Kaplan–Meier (KM) survival analyses for gene-level and pathway-level features in the METABRIC (training) and TCGA-BRCA (validation) cohorts. All analyses focus on overall survival (OS) and are aligned with the resistance-informed prognostic framework described in the manuscript.

---

## METABRIC_KM_Analysis_Gene_OS.R

This script performs gene-level Kaplan–Meier survival analyses in the METABRIC cohort for genes identified as independently prognostic in Cox proportional hazards (CoxPH) modeling.

### Key steps:
- Loads preprocessed METABRIC gene expression, CNA, and clinical OS data
- Stratifies samples based on:
  - **Gene expression** using optimal cutoffs derived from maximally selected rank statistics (maxstat)
  - **Gene-level CNA status** using cohort-specific median CNA values
- Generates KM survival curves for each selected gene
- Computes log-rank test p-values to assess survival differences between risk groups
- Estimates hazard ratios (HRs) and 95% confidence intervals (CIs) using multivariate CoxPH models adjusted for:
  - Age  
  - Tumor stage  
  - PAM50 subtype  
- Outputs publication-ready KM plots with:
  - Survival time shown in months
  - Numbers at risk displayed below each curve

This script supports the training-cohort evaluation of gene-level prognostic markers.

---

## METABRIC_KM_Analysis_Pathway_OS.R

This script conducts pathway-level Kaplan–Meier survival analyses in the METABRIC cohort using pathway-level CNA burden scores derived from paclitaxel resistance–associated pathways.

### Key steps:
- Loads pathway-level CNA scores and METABRIC clinical survival data
- Stratifies patients into low and high CNA burden groups using cohort-specific median cutoffs
- Performs KM survival analyses for each pathway
- Evaluates statistical significance using log-rank tests
- Estimates HRs and 95% CIs using CoxPH models (without clinical covariates to isolate molecular effects)
- Generates KM plots illustrating survival stratification by pathway-level CNA burden

This script establishes pathway-level prognostic effects in the training cohort.

---

## TCGA-BRCA_KM_Analysis_Gene_OS.R

This script replicates gene-level KM survival analyses in the independent TCGA-BRCA validation cohort using cutoffs and features defined in METABRIC.

### Key steps:
- Loads preprocessed TCGA-BRCA gene expression, CNA, and OS data
- Applies:
  - Cohort-specific maxstat-derived expression cutoffs
  - Cohort-specific median CNA cutoffs for gene-level CNAs
- Generates KM survival curves for the same set of prognostic genes identified in METABRIC
- Performs log-rank tests to assess survival differences
- Estimates HRs and 95% CIs using multivariate CoxPH models adjusted for:
  - Age  
  - Tumor stage  
  - PAM50 subtype  
- Produces KM plots consistent in style and structure with METABRIC analyses

This script enables independent validation of gene-level prognostic markers.

---

## TCGA-BRCA_KM_Analysis_Pathway_OS.R

This script evaluates pathway-level CNA burden in TCGA-BRCA using Kaplan–Meier analyses with cutoffs defined in the METABRIC training cohort.

### Key steps:
- Loads TCGA-BRCA pathway-level CNA scores and clinical survival data
- Applies **median CNA burden cutoffs derived exclusively from METABRIC**, without re-optimization
- Performs KM survival analyses for each pathway
- Computes log-rank test statistics
- Estimates HRs and 95% CIs using CoxPH models without clinical covariates
- Generates KM plots directly comparable to METABRIC pathway-level analyses

This script provides strict cross-cohort validation of pathway-level CNA-based survival stratification.

---

## Notes
- All analyses focus on **overall survival (OS)**.
- Survival time is expressed in **months** in all KM plots.
- Numbers at risk are displayed below each KM curve.
