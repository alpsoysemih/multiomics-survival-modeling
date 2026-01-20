# Meta-Analyses of Prognostic Associations

This directory contains R scripts implementing fixed-effect meta-analyses to evaluate the cross-cohort consistency of gene-level and pathway-level prognostic associations identified in METABRIC (training cohort) and TCGA-BRCA (independent validation cohort).

Meta-analyses combine hazard ratio (HR) estimates derived separately from each cohort to quantify shared prognostic effects while accounting for cohort-specific variance.

---

## OS_Gene_Meta-Analysis.R

This script performs gene-level fixed-effect meta-analyses of Cox proportional hazards (CoxPH) results obtained from METABRIC and TCGA-BRCA.

### Key steps:
- Loads gene-level CoxPH results generated independently in:
  - METABRIC_Survival_Analysis_Gene_OS.R
  - TCGA-BRCA_Survival_Analysis_Gene_OS.R
- Selects genes evaluated in both cohorts using identical model specifications
- Extracts:
  - Log hazard ratios (log HRs)
  - Standard errors
  - 95% confidence intervals (CIs)
- Performs fixed-effect meta-analysis using inverse-variance weighting:
  - Pooled log HR = weighted mean of cohort-specific log HRs
  - Weights defined by the inverse of the variance
- Computes pooled HRs and corresponding 95% CIs
- Assesses statistical significance of combined effects
- Outputs meta-analysis tables used for:
  - Forest plots
  - Cross-cohort reproducibility assessment
  - Selection of robust, reproducible prognostic genes

This script identifies gene-level prognostic signals that are consistent across independent datasets.

---

## OS_Pathway_Meta-Analysis.R

This script performs pathway-level fixed-effect meta-analyses of prognostic associations derived from pathway activity scores.

### Key steps:
- Loads pathway-level CoxPH results from:
  - METABRIC_Survival_Analysis_Pathway_OS.R
  - TCGA-BRCA_Survival_Analysis_Pathway_OS.R
- Separately processes:
  - Expression-based pathway scores
  - CNA-based pathway scores
- Extracts log HRs, standard errors, and CIs for each pathway
- Applies inverse-variance weighted fixed-effect meta-analysis
- Estimates pooled HRs and 95% CIs for each pathway
- Identifies pathways with consistent prognostic effects across cohorts
- Produces summary tables for:
  - Forest plots
  - Comparative pathway-level evaluation
  - Downstream integration with gene-signature analyses

This script establishes which paclitaxel resistance–associated pathways retain prognostic relevance across cohorts.

---

## Notes
- All meta-analyses use fixed-effect models to emphasize shared biological signals rather than cohort heterogeneity.
- Meta-analyses are performed on the log HR scale.
- Results are used directly in manuscript figures and tables.
- No additional model tuning or cohort-specific optimization is performed during meta-analysis.
