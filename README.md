# PILCRS-pesticide-liver-Thailand
Analysis code for "Urinary Pesticide Profiles and Liver Disease Risk in Thailand: A Machine Learning Risk Prediction Model" (Cancer Research Communications). Code only; data available on request.

# PILCRS-pesticide-liver-Thailand

Analysis code for:

**Urinary Pesticide Profiles and Liver Disease Risk in Thailand: A Machine Learning Risk Prediction Model**
*Cancer Research Communications* (under review).

## Scripts

| Script | Produces |
|---|---|
| `R/PILCRS_Tables_11222026_42.R` | Tables 1–3 and Supplementary Tables S2–S17 |
| `R/PILCRS_Figures_11222026_41.R` | Figures 1–3 and Supplementary Figures S1–S6 (with Supplementary Tables S18–S19 values) |

## How to run

1. Open a script in RStudio and click **Source**.
2. A file window opens for each dataset that is not already loaded; select the matching Excel file.
3. Tables and figures are shown in RStudio and saved to an `output` folder (tables as PDF and CSV; figures as 600 dpi TIFF).

Required R packages are listed at the top of each script.

## Data availability

This repository contains analysis code only. No participant data are included.

Individual-level data from the TIGER-LC study are not publicly available because of participant
privacy and TIGER-LC Consortium governance requirements. De-identified data and data dictionaries
are available on reasonable request from the corresponding authors, following approval of a
methodologically sound proposal and completion of a data-use agreement, subject to institutional
and ethical restrictions.

## License

MIT License. See `LICENSE`.
