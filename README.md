# multiomics-survival-modeling
Multi-omics survival modeling in breast cancer integrating gene expression and copy number alterations to identify paclitaxel resistance–informed prognostic genes, pathways, and a validated gene signature. Includes CoxPH, KM, RMST, time-dependent AUC, pathway scoring, network analysis, and cross-cohort validation in METABRIC and TCGA-BRCA.



# Paclitaxel Resistance–Informed Multi-Omics Survival Modeling in Breast Cancer
This repository contains the full analysis pipeline accompanying the manuscript:

“Multi-omics survival modeling identifies paclitaxel resistance–informed prognostic genes and pathways and derives a validated gene signature in breast cancer.”

The code implements a resistance-informed survival modeling framework that integrates gene expression, copy number alterations (CNAs), and clinical covariates to identify robust prognostic genes, pathways, and a parsimonious gene-level prognostic signature validated across independent cohorts.

⸻

Overview

Most prognostic biomarkers in breast cancer are identified independently of treatment-related biology. In this study, we explicitly anchor survival analyses to paclitaxel resistance–associated genes and pathways, previously identified using deep learning–based drug response modeling and pathfindR enrichment analysis.

Using two independent breast cancer cohorts (METABRIC and TCGA-BRCA), we perform:

- Gene-level and pathway-level Cox proportional hazards (CoxPH) modeling  
- Integration of gene expression and copy number alteration (CNA) data  
- Clinical covariate adjustment (age, tumor stage, PAM50 subtype)  
- Cross-cohort validation with fixed, training-derived cutoffs  
- Kaplan–Meier (KM) analysis, log-rank testing, RMST, and time-dependent AUC analyses  
- Network-based contextualization and hypergeometric over-representation testing  

The repository provides fully reproducible code for all analyses and figures reported in the manuscript.

⸻

Data Sources

This project uses publicly available datasets:

- **METABRIC**
  - Gene expression data (Illumina microarrays, z-scores)
  - Copy number alteration (CNA) data (GISTIC scores)
  - Clinical survival data
  - Source: Synapse (syn1688369)

- **TCGA-BRCA**
  - RNA-seq gene expression data (RSEM → TPM → log2 transformation → z-score normalization)
  - CNA profiles (Affymetrix SNP 6.0 arrays; CBS segmentation)
  - Clinical data retrieved using the `TCGAbiolinks` R package

All analyses use overall survival (OS) as the endpoint.

⚠️ Raw data are not redistributed in this repository. Scripts assume the user has appropriate access to the original datasets.

⸻

Analysis Workflow

## Data Preprocessing
- Expression normalization and z-score standardization
- CNA binarization (altered vs. copy-neutral)
- Clinical covariate harmonization (tumor stage and PAM50 subtype)

## Gene-Level Survival Modeling
- Cox proportional hazards (CoxPH) models using:
  - Expression only
  - CNA only
  - Expression + CNA
  - Expression × CNA interaction
- Models evaluated with and without clinical covariates
- Cross-cohort meta-analysis of hazard ratios

## Pathway-Level Survival Modeling
- Pathway activity scores derived from paclitaxel resistance–associated pathways
- Separate expression-based and CNA-based pathway scores
- Identification of independently prognostic pathways

## Prognostic Signature Construction
- Selection of non-redundant prognostic genes
- Multicollinearity assessment using:
  - Point-biserial correlation
  - Phi coefficient
- Model selection based on:
  - Concordance index (C-index)
  - Akaike information criterion (AIC)
  - Likelihood ratio tests (LRT)
- Final gene signature trained in METABRIC and validated in TCGA-BRCA

## Survival Stratification and Temporal Evaluation
- Kaplan–Meier analyses with fixed cutoffs
- Restricted Mean Survival Time (RMST) analysis
- Time-dependent AUC analyses up to 200 months

## Network Analysis
- Integration of prognostic genes and pathways
- Visualization of gene–pathway relationships
- Hypergeometric testing for non-random accumulation of prognostic genes

⸻

## Key Findings Implemented in This Repository
- Identification of **SERPINE1** (expression) and **CLDN1, CXCR4, LAMB2, and LYPD6B** (CNAs) as reproducible prognostic genes across independent cohorts
- Construction of a parsimonious gene-level prognostic signature with robust cross-cohort validation
- Identification of **Ras signaling** as the only pathway retaining independent prognostic significance after clinical adjustment
- Demonstration of a non-random accumulation of prognostic genes within the **Axon guidance** pathway
- Evidence that gene-level prognostic signatures outperform pathway-level CNA scores in survival modeling

⸻

Software & Dependencies

## Software and Package Versions
All analyses were performed in **R (version ≥ 4.4.2)**.

### Core R Packages
- **survival** (v3.8.3): Cox proportional hazards modeling, C-index, AIC, and likelihood ratio tests
- **survminer** (v0.5.1): Kaplan–Meier visualization and log-rank testing
- **survRM2** (v1.0.4): Restricted mean survival time (RMST) analysis
- **timeROC** (v0.4): Time-dependent AUC analysis
- **TCGAbiolinks** (v2.34.1): Retrieval of TCGA-BRCA clinical and molecular data
- **dplyr** (v1.1.4) and **tidyr** (v1.3.1): Data manipulation and preprocessing
- **ggplot2** (v4.0.1): Figure generation and visualization
- **igraph** (v2.2.1): Network construction and hypergeometric over-representation testing
- **stats** (base R): Statistical testing and meta-analysis utilities

### Notes
- Meta-analyses were performed using inverse-variance weighting on log hazard ratios.
- All figures in the repository were generated using **ggplot2** and related visualization utilities.
- Exact package versions and methodological details are reported in the associated manuscript.

⸻

## Reproducibility
- **METABRIC** is used as the training cohort.
- **TCGA-BRCA** serves as an independent validation cohort.

- **Gene-level analyses**
  - Expression-based gene stratification uses **cohort-specific cutoffs**, determined by maximally selected rank statistics (maxstat).
  - CNA-based gene stratification uses **cohort-specific median cutoffs**.

- **Pathway-level and gene-signature analyses**
  - All cutoffs are derived **exclusively from METABRIC**.
  - The same predefined cutoffs are applied **unchanged** to TCGA-BRCA.

- No cohort-specific tuning, re-optimization, or threshold adjustment is performed for pathway-level or signature-level validation.


⸻

Contact

Corresponding author:
Osman Ugur Sezerman
📧 ugur.sezerman@acibadem.edu.tr
