# PILCRS: Pesticide-Informed Liver Cancer Risk Score

Analysis code accompanying:

> **Urinary Pesticide Profiles and Liver Disease Risk in Thailand: A Machine Learning Risk Prediction Model**
> Patel DP, Loffredo CA, et al., on behalf of the TIGER-LC Consortium.
> *Cancer Research Communications* (under review).

> **Important:** This repository contains analysis code only. No participant-level data are included.
> The scripts require the TIGER-LC study datasets, which are available on request (see **Data availability**).

## Overview

This repository provides the R code used to generate all main and supplementary tables and figures in the
manuscript. The analyses evaluate creatinine-normalized urinary pesticide burden, summarized as Pesticide
Load Scores (PLS11, including glyphosate-pathway analytes; PLS8, excluding them) and urinary cypermethrin,
in relation to chronic liver disease (CLD) and hepatocellular carcinoma (HCC) in a multicenter hospital-based
case-control study in Thailand (228 CLD, 116 HCC, 249 hospital-based controls). The code covers:

- Descriptive statistics and group comparisons
- Quartile-based logistic regression, adjusted and stratified models, and sensitivity analyses
- Multiple imputation by chained equations for missing cypermethrin
- Development of the Pesticide-Informed Liver Cancer Risk Score (PILCRS) and gradient-boosted (XGBoost) classifiers
- Internal validation by 10-fold cross-validation and 1,000-iteration bootstrap optimism correction
- Model interpretation (SHAP), reclassification (NRI, IDI), E-values, and population attributable fractions

## Repository contents

| Script | Generates |
|---|---|
| `PILCRS_Tables.R` | Tables 1–3; Supplementary Tables S2–S17 |
| `PILCRS_Figures.R` | Figures 1–3; Supplementary Figures S1–S6; values reported in Supplementary Tables S18–S19 |

Supplementary Table S1 (LC-MS/MS assay parameters) is descriptive and has no associated code.

## Requirements

- R version 4.6.1 (tested) and RStudio
- R packages (installed automatically if missing):
  - Tables: `readxl`, `dplyr`, `tidyr`, `tibble`, `gridExtra`, `logistf`, `pROC`, `mice`
  - Figures: `readxl`, `ggplot2`, `ggpubr`, `dplyr`, `tibble`, `rlang`, `fmsb`, `ggsignif`, `rstatix`,
    `readr`, `xgboost`, `FSA`, `patchwork`, `pROC`, `boot`, `SHAPforxgboost`, `data.table`, `tidyverse`, `broom`

## How to run

1. Open `PILCRS_Tables.R` or `PILCRS_Figures.R` in RStudio and click **Source**.
2. For each dataset not already loaded, a file window opens; select the corresponding Excel file
   (dataset names are listed at the top of each script).
3. Results are displayed in RStudio and saved to an `output` folder:
   - Tables: printed to the console, kept in the R environment (`Table_1`, `Table_2`, ..., `Table_S17`),
     and saved as PDF and CSV files.
   - Figures: displayed in the Plots pane and saved as 600 dpi TIFF files.
4. A summary at the end of the tables script reports the status of each analysis step.

Approximate run time on a standard desktop: tables, a few minutes; figures, about 6 minutes
(dominated by the 1,000-iteration bootstrap).

## Reproducibility

- Random seeds are fixed (`set.seed(2025)`) for cross-validation, bootstrap resampling, and multiple
  imputation, so repeated runs reproduce identical results.
- Quartile cutpoints are fixed to those reported in Table 2 so that all quartile-based analyses use
  the same exposure categories.
- No results are entered manually; all values are computed from the input datasets at run time.

## Data availability

Individual-level data from the TIGER-LC study are not publicly available because of participant privacy
and TIGER-LC Consortium governance requirements. De-identified data and data dictionaries are available on
reasonable request from the corresponding authors, following approval of a methodologically sound proposal
and completion of a data-use agreement, subject to institutional and ethical restrictions.

## Citation

If you use this code, please cite the article above. Full citation details will be added on publication.

## License

Released under the MIT License. See `LICENSE`.
