# PILCRS: Pesticide-Informed Liver Cancer Risk Score

Analysis code accompanying:

> **Urinary Pesticide Profiles and Liver Disease Risk in Thailand: A Machine Learning Risk Prediction Model**
> Patel DP, Loffredo CA, et al., on behalf of the TIGER-LC Consortium. *Cancer Research Communications* (under review).

> **Important:** This repository contains analysis code only. No participant-level data are included. The scripts require the TIGER-LC study datasets, which are available on request (see **Data availability**).

## Overview

This repository provides the R code used to generate all main and supplementary tables and figures in the manuscript. The analyses evaluate creatinine-normalized urinary pesticide burden, summarized as Pesticide Load Scores (PLS11, including glyphosate-pathway analytes; PLS8, excluding them) and urinary cypermethrin, in relation to chronic liver disease (CLD) and hepatocellular carcinoma (HCC) in a multicenter hospital-based case-control study in Thailand (228 CLD, 116 HCC, 249 hospital-based controls).

Analyses include:

- Descriptive statistics and group comparisons
- Quartile-based logistic regression, adjusted and stratified models, and sensitivity analyses
- Multiple imputation by chained equations for missing cypermethrin
- Development of the Pesticide-Informed Liver Cancer Risk Score (PILCRS) and gradient-boosted (XGBoost) classifiers
- Internal validation by 10-fold cross-validation and 1,000-iteration bootstrap optimism correction
- Model interpretation (SHAP), reclassification (NRI, IDI), E-values, and population attributable fractions

## Model specification

| Component | Specification |
|---|---|
| Outcomes | CLD vs hospital-based controls; HCC vs hospital-based controls |
| Exposure metrics | PLS11, PLS8, creatinine-normalized cypermethrin |
| Base covariates | Age, sex, alcohol use, agricultural occupation |
| Extended covariates | Base covariates plus HBV and HCV serostatus |
| Regression model | Multivariable logistic regression (PILCRS) |
| Machine learning model | XGBoost (maximum depth 3, learning rate 0.1, 300 rounds) |
| Internal validation | 10-fold cross-validation; 1,000 bootstrap resamples with optimism correction |
| Performance measures | AUC, Brier score, calibration slope, calibration plots, NRI, IDI |

## Repository contents

| Script | Generates |
|---|---|
| `PILCRS_Tables.R` | Tables 1–3; Supplementary Tables S2–S17 |
| `PILCRS_Figures.R` | Figures 1–3; Supplementary Figures S1–S6; values reported in Supplementary Tables S18–S19 |

Supplementary Table S1 (LC-MS/MS assay parameters) is descriptive and has no associated code.

## Requirements

- R version 4.6.1
- R packages:
  - Tables: `readxl`, `dplyr`, `tidyr`, `tibble`, `gridExtra`, `logistf`, `pROC`, `mice`
  - Figures: `readxl`, `ggplot2`, `ggpubr`, `dplyr`, `tibble`, `rlang`, `fmsb`, `ggsignif`, `rstatix`, `readr`, `xgboost`, `FSA`, `patchwork`, `pROC`, `boot`, `SHAPforxgboost`, `data.table`, `tidyverse`, `broom`

## Usage

1. Open `PILCRS_Tables.R` or `PILCRS_Figures.R` in RStudio and run the script.
2. When prompted, select the Excel file for each required dataset. Dataset names are listed at the top of each script.
3. Tables are returned to the R environment (`Table_1` to `Table_S17`) and saved as PDF and CSV files in an `output` folder. Figures are displayed in RStudio and saved as 600 dpi TIFF files in the same folder.

## Reproducibility

- Random seeds are fixed (`set.seed(2025)`) for cross-validation, bootstrap resampling, and multiple imputation.
- Quartile cutpoints are fixed to those reported in Table 2.
- All reported values are computed from the input datasets at run time.

## Reporting standards

The study is reported in accordance with the STROBE statement for observational studies and the TRIPOD+AI statement for clinical prediction models using regression or machine learning methods.

## Intended use

The code is provided to support transparency and reproducibility of the published analyses. PILCRS has been internally validated only and is not intended for clinical decision-making. External validation in independent, population-based cohorts is required before any clinical application.

## Data availability

Individual-level data from the TIGER-LC study are not publicly available because of participant privacy and TIGER-LC Consortium governance requirements. De-identified data and data dictionaries are available on reasonable request from the corresponding authors, following approval of a methodologically sound proposal and completion of a data-use agreement, subject to institutional and ethical restrictions.

## Funding

Supported in part by the Intramural Research Program of the National Institutes of Health, National Cancer Institute, Center for Cancer Research, and by the Chulabhorn Research Institute and Thailand Science Research and Innovation.

## Contact

Correspondence regarding the study and data access should be addressed to the corresponding authors listed in the article.

## Citation

Please cite the article above. Full citation details will be added on publication.

## License

Released under the MIT License. See `LICENSE`.
