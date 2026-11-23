###############################################################################
# PILCRS_Figures_11222026_41.R
# All main and supplementary FIGURES for:
#   "Urinary Pesticide Profiles and Liver Disease Risk in Thailand: A Machine
#    Learning Risk Prediction Model" (Cancer Research Communications)
#
# Order: Figure 1A, 1B-1D, 2A, 2B, 2C, 3 (with Supplementary Figures S3-S4),
#        then Supplementary Figures S1, S2, S5, S6.
#
# HOW TO RUN
#   When you click "Source", a file window opens for each dataset that is not
#   already loaded; pick the matching Excel file (first sheet is read).
#   Datasets (object names used in the script):
#     Pesticides_Northern_hospital_samples_Madja_LHC_1, TIGERL_2025_Madja_IJC_Dax01,
#     Copy_of_TIGERL_2025_Madja_IJC_Dax01, PILCRS_1, PILCRS_1_MASTER,
#     PILCRS_1_MASTER1, Risk_Per_Calc
#   Then click "Source" (Ctrl+Shift+S). Every plot is printed to the RStudio
#   Plots pane (use the arrows to page through them).
#   Plots appear in the RStudio Plots pane; 600 dpi TIFF files are saved to "output".
#
# Required packages: ggplot2, ggpubr, dplyr, tibble, rlang, fmsb, ggsignif, rstatix,
#   readr, xgboost, FSA, patchwork, pROC, boot, SHAPforxgboost, data.table,
#   tidyverse, broom
###############################################################################

###dir.create("output", showWarnings = FALSE)

## ---------------- Output folder ----------------
OUT_DIR <- "H:/Sep2023_Desktop Backup/TGLC_Chris_01112021/Pesticide Load Project/2025_Madja_Dax_Combine/Manuscript/Submissions/Cancer Research Communication/Major Revision_1/Data Analysis Final/Github/R/Figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Cannot reach the Figures folder: ", OUT_DIR)

# Send every file.path("output", ...) in this script to OUT_DIR
file.path <- function(...) {
  a <- list(...)
  if (length(a) && identical(a[[1]], "output")) a[[1]] <- OUT_DIR
  do.call(base::file.path, a)
}

# Saves a ggplot or patchwork figure to OUT_DIR as a 600 dpi TIFF file
save_fig <- function(name, plot, width = 7, height = 6) {
  f <- file.path(OUT_DIR, paste0(name, ".tiff"))
  ggplot2::ggsave(f, plot = plot, width = width, height = height, units = "in",
                  dpi = 600, device = "tiff", compression = "lzw")
  cat("Saved:", f, "\n")
}

## ---------------- Load the Excel datasets (pick each file manually) ----------------
# Datasets already in the RStudio environment are kept; for any that are missing,
# a window opens so you can choose the Excel file yourself.
load_excel_manually <- function(objects) {
  if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
  for (obj in objects) {
    if (exists(obj, envir = globalenv())) {
      cat("Already loaded:", obj, "\n")
      next
    }
    if (!interactive()) stop("Dataset not loaded: ", obj, ". Load it in RStudio before running.")
    cat("\nSelect the Excel file for:", obj, "\n")
    f <- file.choose()
    assign(obj, readxl::read_excel(f, sheet = 1, .name_repair = "unique_quiet"), envir = globalenv())
    cat("Loaded", obj, "(", nrow(get(obj, envir = globalenv())), "rows ) from", basename(f), "\n")
  }
}

load_excel_manually(c(
  "Pesticides_Northern_hospital_samples_Madja_LHC_1",  # Pesticides_Northern_hospital_samples_Madja-LHC_1.xlsx
  "TIGERL_2025_Madja_IJC_Dax01",                       # TIGERL_2025_Madja_IJC_Dax01.xlsx
  "Copy_of_TIGERL_2025_Madja_IJC_Dax01",               # Copy_of_TIGERL_2025_Madja_IJC_Dax01.xlsx
  "PILCRS_1",                                          # PILCRS_1.xlsx
  "PILCRS_1_MASTER",                                   # PILCRS_1_MASTER.xlsx
  "PILCRS_1_MASTER1",                                  # PILCRS_1_MASTER1.xlsx
  "Risk_Per_Calc"                                      # Risk_Per_Calc.xlsx
))


###############################################################################
###############################################################################
# FIGURE 1A: RADAR PLOT OF MEDIAN URINARY PESTICIDE CONCENTRATIONS
###############################################################################
###############################################################################

################## Radar Plot (Final, 600 DPI Export) ########################

suppressPackageStartupMessages({
  library(dplyr)
  library(fmsb)
  library(tibble)
})

#=========================#
# Load and Prepare Data
#=========================#
data_file <- Pesticides_Northern_hospital_samples_Madja_LHC_1

# Locate tissue/group column
candidate_cols <- c("tissuetype","tissue_type","TissueType","Tissue","Group","group","cohort","Cohort")
tt_col <- intersect(candidate_cols, names(data_file))
if (length(tt_col) == 0) {
  stop(paste0("Could not find a tissue/group column. Columns: ", paste(names(data_file), collapse = ", ")))
}
tt_col <- tt_col[1]

# Clean & map to HBC/CLD/HCC
tt_raw <- toupper(trimws(as.character(data_file[[tt_col]])))
tt_raw <- gsub("^PC$", "HBC", tt_raw)

map_to_group <- function(x) {
  is_hbc <- grepl("\\b(HBC|HEALTHY|CONTROL|NON[- ]?CANCER|CTRL)\\b", x)
  is_cld <- grepl("\\b(CLD|CLF|CHRONIC LIVER|CIRRHOSIS|FIBROSIS|CHRONIC\\s*LD)\\b", x)
  is_hcc <- grepl("\\b(HCC|HEPATOCELLULAR|CARCINOMA)\\b", x)
  ifelse(is_hbc, "HBC", ifelse(is_cld, "CLD", ifelse(is_hcc, "HCC", NA)))
}
data_file$group_std <- map_to_group(tt_raw)

cat("\n=== Mapping check (Original -> Mapped) ===\n")
print(table(Original = tt_raw, Mapped = data_file$group_std, useNA = "ifany"))

#=========================#
# Define pesticides + rename mapping
#=========================#
pesticide_cols <- c(
  "Pendimethalin_Concnormcreatinine",
  "metsulfuronmethyl_Conc.normcreatinine",
  "dichlorophenoxyaceticacid_conc_normcreatininecolourimetric",
  "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine",
  "bromodiolon_calcresponse",
  "butachlor_calcresponse",
  "oxadiazon_calcresponse",
  "Cypermethrin (nM)",
  "Glyphosate (nM)_M",
  "Phosphoric acid(nM)_M",
  "AMPA (nM)_M"
)

pesticide_rename <- c(
  "Pendimethalin_Concnormcreatinine" = "Pendimethalin (nM)",
  "metsulfuronmethyl_Conc.normcreatinine" = "Metsulfuron-methyl (nM)",
  "dichlorophenoxyaceticacid_conc_normcreatininecolourimetric" = "2,4-Dichlorophenoxyacetic acid (nM)",
  "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine" = "Flocoumafen (nM)",
  "bromodiolon_calcresponse" = "Bromadiolone (pM)",
  "butachlor_calcresponse" = "Butachlor (pM)",
  "oxadiazon_calcresponse" = "Oxadiazon (pM)",
  "Cypermethrin (nM)" = "Cypermethrin (nM)",
  "Glyphosate (nM)_M" = "Glyphosate (nM)",
  "Phosphoric acid(nM)_M" = "PPA (nM)",
  "AMPA (nM)_M" = "AMPA (nM)"
)

#=========================#
# Compute Median per Group
#=========================#
df_medians <- data_file %>%
  filter(group_std %in% c("HBC", "CLD", "HCC")) %>%
  select(group_std, dplyr::all_of(pesticide_cols)) %>%
  group_by(group_std) %>%
  summarise(across(dplyr::all_of(pesticide_cols), ~ median(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  column_to_rownames("group_std") %>%
  t() %>% as.data.frame()

if (nrow(df_medians) == 0) {
  stop("After mapping, no rows found for HBC/CLD/HCC. Check the 'Mapping check' table above and adjust patterns if needed.")
}

# Apply display names
rownames(df_medians) <- pesticide_rename[rownames(df_medians)]

#=========================#
# Normalize Each Pesticide (Row-wise 0–1)
#=========================#
# Save group names BEFORE normalization (they are the *columns* right now)
group_names <- colnames(df_medians)

df_scaled <- as.data.frame(t(apply(df_medians, 1, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    rep(0.5, length(x))   # flat/all-NA -> midline
  } else {
    (x - rng[1]) / (rng[2] - rng[1])
  }
})))

# Transpose back so rows = groups, and RESTORE the group names explicitly
df_scaled <- as.data.frame(t(df_scaled))
rownames(df_scaled) <- group_names

#=========================#
# Enforce Group Order + Aesthetics
#=========================#
group_order    <- c("HBC", "CLD", "HCC")
groups_present <- intersect(group_order, rownames(df_scaled))
if (length(groups_present) == 0) {
  cat("\nRow names after normalization:\n"); print(rownames(df_scaled))
  stop("No groups found among HBC/CLD/HCC after normalization. (Row names printed above.)")
}
df_scaled <- df_scaled[groups_present, , drop = FALSE]

# Stable palette
pal <- c(HBC = "#1f77b4", CLD = "#2ca02c", HCC = "#d62728")
pcols  <- unname(pal[groups_present])
pfcols <- adjustcolor(pcols, alpha.f = 0.25)
line_widths <- rep(1.8, length(groups_present))
line_types  <- rep(1,   length(groups_present))

# fmsb data
radar_data <- rbind(
  rep(1, ncol(df_scaled)),
  rep(0, ncol(df_scaled)),
  df_scaled
)
colnames(radar_data) <- colnames(df_scaled)

#=========================#
# Export (600 DPI TIFF)
#=========================#
output_path <- file.path("output", "Radar_Plot_Pesticide_Median_1.tiff")
tiff(output_path, width = 17, height = 12.5, units = "in", res = 600, compression = "lzw")

#=========================#
# Plot Radar Chart
#=========================#
par(cex.axis = 3.5, font.axis = 2)
par(cex = 1.2, text.font = 2)

radarchart(
  radar_data,
  axistype    = 1,
  pcol        = pcols,
  pfcol       = pfcols,
  plwd        = line_widths,
  plty        = line_types,
  cglcol      = "black",
  cglty       = 2,
  axislabcol  = "black",
  caxislabels = seq(0, 1, 0.25),
  cglwd       = 1.5,
  vlcex       = 1.6,
  title       = ""
)

legend("topright",
       legend   = groups_present,
       bty      = "n",
       pch      = 20,
       col      = pcols,
       text.col = "black",
       text.font= 2,
       pt.cex   = 5,
       cex      = 2.5
)

legend("bottomleft",
       title = "Exposure Level (normalized)",
       legend = c("High (1.0)", "↑", "Medium (0.5)", "↑", "Low (0.0)"),
       text.col = "black",
       text.font = 2,
       cex = 1.1,
       bty = "n"
)

dev.off()






























####################################################
########################################Final###########################################
#########################################09152025#####################################












################### Radar Plot (Final, 600 DPI Export) ########################

suppressPackageStartupMessages({
  library(dplyr)
  library(fmsb)
  library(tibble)
})

#=========================#
# Load and Prepare Data
#=========================#
data_file <- Pesticides_Northern_hospital_samples_Madja_LHC_1

# Locate tissue/group column
candidate_cols <- c("tissuetype","tissue_type","TissueType","Tissue","Group","group","cohort","Cohort")
tt_col <- intersect(candidate_cols, names(data_file))
if (length(tt_col) == 0) {
  stop(paste0("Could not find a tissue/group column. Columns: ", paste(names(data_file), collapse = ", ")))
}
tt_col <- tt_col[1]

# Clean & map to HBC/CLD/HCC
tt_raw <- toupper(trimws(as.character(data_file[[tt_col]])))
tt_raw <- gsub("^PC$", "HBC", tt_raw)

map_to_group <- function(x) {
  is_hbc <- grepl("\\b(HBC|HEALTHY|CONTROL|NON[- ]?CANCER|CTRL)\\b", x)
  is_cld <- grepl("\\b(CLD|CLF|CHRONIC LIVER|CIRRHOSIS|FIBROSIS|CHRONIC\\s*LD)\\b", x)
  is_hcc <- grepl("\\b(HCC|HEPATOCELLULAR|CARCINOMA)\\b", x)
  ifelse(is_hbc, "HBC", ifelse(is_cld, "CLD", ifelse(is_hcc, "HCC", NA)))
}
data_file$group_std <- map_to_group(tt_raw)

cat("\n=== Mapping check (Original -> Mapped) ===\n")
print(table(Original = tt_raw, Mapped = data_file$group_std, useNA = "ifany"))

#=========================#
# Define pesticides + rename mapping
#=========================#
pesticide_cols <- c(
  "Pendimethalin_Concnormcreatinine",
  "metsulfuronmethyl_Conc.normcreatinine",
  "dichlorophenoxyaceticacid_conc_normcreatininecolourimetric",
  "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine",
  "bromodiolon_calcresponse",
  "butachlor_calcresponse",
  "oxadiazon_calcresponse",
  "Cypermethrin (nM)",
  "Glyphosate (nM)_M",
  "Phosphoric acid(nM)_M",
  "AMPA (nM)_M"
)

pesticide_rename <- c(
  "Pendimethalin_Concnormcreatinine" = "Pendimethalin (nM)",
  "metsulfuronmethyl_Conc.normcreatinine" = "Metsulfuron-methyl (nM)",
  "dichlorophenoxyaceticacid_conc_normcreatininecolourimetric" = "2,4-Dichlorophenoxyacetic acid (nM)",
  "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine" = "Flocoumafen (nM)",
  "bromodiolon_calcresponse" = "Bromadiolone (pM)",
  "butachlor_calcresponse" = "Butachlor (pM)",
  "oxadiazon_calcresponse" = "Oxadiazon (pM)",
  "Cypermethrin (nM)" = "Cypermethrin (nM)",
  "Glyphosate (nM)_M" = "Glyphosate (nM)",
  "Phosphoric acid(nM)_M" = "PPA (nM)",
  "AMPA (nM)_M" = "AMPA (nM)"
)

#=========================#
# Compute Median per Group
#=========================#
df_medians <- data_file %>%
  filter(group_std %in% c("HBC", "CLD", "HCC")) %>%
  select(group_std, dplyr::all_of(pesticide_cols)) %>%
  group_by(group_std) %>%
  summarise(across(dplyr::all_of(pesticide_cols), ~ median(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  column_to_rownames("group_std") %>%
  t() %>% as.data.frame()

if (nrow(df_medians) == 0) {
  stop("After mapping, no rows found for HBC/CLD/HCC. Check the 'Mapping check' table above and adjust patterns if needed.")
}

# Apply display names
rownames(df_medians) <- pesticide_rename[rownames(df_medians)]

#=========================#
# Normalize Each Pesticide (Row-wise 0–1)
#=========================#
# Save group names BEFORE normalization (they are the *columns* right now)
group_names <- colnames(df_medians)

df_scaled <- as.data.frame(t(apply(df_medians, 1, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    rep(0.5, length(x))   # flat/all-NA -> midline
  } else {
    (x - rng[1]) / (rng[2] - rng[1])
  }
})))

# Transpose back so rows = groups, and RESTORE the group names explicitly
df_scaled <- as.data.frame(t(df_scaled))
rownames(df_scaled) <- group_names

#=========================#
# Enforce Group Order + Aesthetics
#=========================#
group_order    <- c("HBC", "CLD", "HCC")
groups_present <- intersect(group_order, rownames(df_scaled))
if (length(groups_present) == 0) {
  cat("\nRow names after normalization:\n"); print(rownames(df_scaled))
  stop("No groups found among HBC/CLD/HCC after normalization. (Row names printed above.)")
}
df_scaled <- df_scaled[groups_present, , drop = FALSE]

#=========================#
# Colors (YOUR requested palette)
#=========================#
pal <- c(HBC = "#1F5A99", CLD = "#33a02c", HCC = "#e31a1c")
pcols  <- unname(pal[groups_present])
pfcols <- adjustcolor(pcols, alpha.f = 0.25)
line_widths <- rep(1.8, length(groups_present))
line_types  <- rep(1,   length(groups_present))

#=========================#
# fmsb data
#=========================#
radar_data <- rbind(
  rep(1, ncol(df_scaled)),
  rep(0, ncol(df_scaled)),
  df_scaled
)
colnames(radar_data) <- colnames(df_scaled)

#=========================#
# Small helper to draw the chart (used for preview + export)
#=========================#
.draw_radar <- function() {
  par(cex.axis = 3.5, font.axis = 2)
  par(cex = 1.2, text.font = 2)
  
  radarchart(
    radar_data,
    axistype    = 1,
    pcol        = pcols,
    pfcol       = pfcols,
    plwd        = line_widths,
    plty        = line_types,
    cglcol      = "black",
    cglty       = 2,
    axislabcol  = "black",
    caxislabels = seq(0, 1, 0.25),
    cglwd       = 1.5,
    vlcex       = 1.6,
    title       = ""
  )
  
  legend("topright",
         legend   = groups_present,
         bty      = "n",
         pch      = 20,
         col      = pcols,
         text.col = "black",
         text.font= 2,
         pt.cex   = 5,
         cex      = 2.5
  )
  
  legend("bottomleft",
         title = "Exposure Level (normalized)",
         legend = c("High (1.0)", "↑", "Medium (0.5)", "↑", "Low (0.0)"),
         text.col = "black",
         text.font = 2,
         cex = 1.1,
         bty = "n"
  )
}

#=========================#
# Preview in RStudio (Plots pane)
#=========================#
.draw_radar()   # <-- you’ll see this immediately in RStudio

#=========================#
# Export (600 DPI TIFF)
#=========================#
output_path <- file.path("output", "Radar_Plot_Pesticide_Median_1.tiff")
tiff(output_path, width = 17, height = 12.5, units = "in", res = 600, compression = "lzw")
.draw_radar()
dev.off()


###############################################################################
###############################################################################
# FIGURE 1B-1D: PLS11 FOLD CHANGE, PLS8 FOLD CHANGE, CYPERMETHRIN
###############################################################################
###############################################################################

############################# FINAL FOR THE MANUSCRIPT ##################
################# PLS FCA with standard approach (x-axis updated with n) #################

# Load required libraries
library(ggplot2)
library(ggpubr)

# Customize variables
data_file <- "TIGERL_2025_Madja_IJC_Dax01"
group_column <- "tissuetypeORG"
value_column <- "PLSFC_GroupAM"
y_axis_label <- expression(PLS[11]~"FC")
x_axis_label <- ""

# Load the data
filtered_data <- get(data_file)
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])
filtered_data <- filtered_data[!is.na(filtered_data[[group_column]]) & !is.na(filtered_data[[value_column]]), ]
filtered_data <- filtered_data[filtered_data[[group_column]] %in% c("HBC", "CLD", "HCC"), ]
filtered_data[[group_column]] <- factor(filtered_data[[group_column]], levels = c("HBC", "CLD", "HCC"))

# Count sample sizes
n_counts <- as.data.frame(table(filtered_data[[group_column]]))
colnames(n_counts) <- c("Group", "n")
n_counts$Label <- paste0(n_counts$Group, " (n=", n_counts$n, ")")
group_levels <- c("HBC", "CLD", "HCC")
n_counts <- n_counts[match(group_levels, n_counts$Group), ]

# Add updated group labels to data
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = n_counts$Group,
                                    labels = n_counts$Label)

# Remove near-zero values and add pseudo-count
filtered_data <- filtered_data[filtered_data[[value_column]] >= 0.1, ]
filtered_data[[value_column]] <- log2(filtered_data[[value_column]] + 0.1 + 1)

# Calculate medians for labeling
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", round(medians$Median, 2))

median_positions <- aggregate(filtered_data[[value_column]],
                              by = list(Group = filtered_data$Group_Label),
                              FUN = function(x) quantile(x, 0.5, na.rm = TRUE))
colnames(median_positions) <- c("Group_Label", "Position")
medians$Position <- median_positions$Position

# Define group colors
group_colors <- c("HBC" = "#1f78b4", "CLD" = "#33a02c", "HCC" = "#e31a1c")
names(group_colors) <- n_counts$Label

# Create the plot
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label),
            size = 7, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c(n_counts$Label[1], n_counts$Label[2]),
    c(n_counts$Label[2], n_counts$Label[3]),
    c(n_counts$Label[1], n_counts$Label[3])
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = x_axis_label, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 1.5, color = "black"),
    axis.line.x.bottom = element_line(size = 1.5, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

#============================#
# Save the plot (600 DPI)
#============================#
ggsave(
  filename = file.path("output", "PLSFC_GroupAM_Boxplot_600dpi.tiff"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 600,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

# Show the plot in RStudio
print(p)








############# PLS FCB, REMOVE VALUE ABOVE 5 IN Y AXIS + ADD n IN X LABELS + SAVE 600 DPI #############

# Load required libraries
library(ggplot2)
library(ggpubr)

# Customize these variables as needed
data_file <- "Copy_of_TIGERL_2025_Madja_IJC_Dax01"
group_column <- "tissuetypeORG"
value_column <- "PLSFC_GroupBM"
y_axis_label <- expression(PLS[8]~"FC")  # Corrected subscript for PLS₈
x_axis_label <- ""

# Load the data
filtered_data <- get(data_file)

# Replace "PC" with "HBC" if needed
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

# Remove NA values
filtered_data <- filtered_data[!is.na(filtered_data[[group_column]]) & !is.na(filtered_data[[value_column]]), ]

# Filter only HBC, CLD, HCC
group_levels <- c("HBC", "CLD", "HCC")
filtered_data <- filtered_data[filtered_data[[group_column]] %in% group_levels, ]
filtered_data[[group_column]] <- factor(filtered_data[[group_column]], levels = group_levels)

# Count sample sizes
n_counts <- as.data.frame(table(filtered_data[[group_column]]))
colnames(n_counts) <- c("Group", "n")
n_counts$Label <- paste0(n_counts$Group, " (n=", n_counts$n, ")")
n_counts <- n_counts[match(group_levels, n_counts$Group), ]

# Add group labels with n
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = n_counts$Group,
                                    labels = n_counts$Label)

# Remove near-zero and extreme values
filtered_data <- filtered_data[filtered_data[[value_column]] >= 0.1, ]
filtered_data[[value_column]] <- log2(filtered_data[[value_column]] + 0.1 + 1)
filtered_data <- filtered_data[filtered_data[[value_column]] <= 5, ]

# Median labels
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", round(medians$Median, 2))

median_positions <- aggregate(filtered_data[[value_column]],
                              by = list(Group_Label = filtered_data$Group_Label),
                              FUN = function(x) quantile(x, 0.5, na.rm = TRUE))
colnames(median_positions) <- c("Group_Label", "Position")
medians$Position <- median_positions$Position

# Group colors
group_colors <- c("HBC" = "#1f78b4", "CLD" = "#33a02c", "HCC" = "#e31a1c")
names(group_colors) <- n_counts$Label

# Create plot
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label),
            size = 7, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c(n_counts$Label[1], n_counts$Label[2]),
    c(n_counts$Label[2], n_counts$Label[3]),
    c(n_counts$Label[1], n_counts$Label[3])
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = x_axis_label, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 1.5, color = "black"),
    axis.line.x.bottom = element_line(size = 1.5, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

# Save 600 DPI TIFF
ggsave(
  filename = file.path("output", "PLSFC_GroupBM_Boxplot_600dpi.tiff"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 600,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

# Optionally display the plot
print(p)




















############ FINAL REVISED: Cypermethrin_nM | Fixed n for HBC/HCC | Auto n for CLD | 600 DPI EXPORT ############

# Load required libraries
library(ggplot2)
library(ggpubr)

# Set parameters
data_file <- "TIGERL_2025_Madja_IJC_Dax01"
group_column <- "tissuetypeORG"
value_column <- "Cypermethrin_nM"
y_axis_label <- "Cypermethrin (nM)"

# Load data
filtered_data <- get(data_file)
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

# Filter and order groups
group_levels <- c("HBC", "CLD", "HCC")
filtered_data <- filtered_data[filtered_data[[group_column]] %in% group_levels, ]
filtered_data[[group_column]] <- factor(filtered_data[[group_column]], levels = group_levels)

# Set sample sizes (manual for HBC/HCC, auto for CLD)
n_counts <- data.frame(
  Group = group_levels,
  n = c(249,
        sum(filtered_data[[group_column]] == "CLD"),
        116)
)
n_counts$Label <- paste0(n_counts$Group, " (n=", n_counts$n, ")")

# Add x-axis group labels
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = n_counts$Group,
                                    labels = n_counts$Label)

# Log2-transform Cypermethrin values
filtered_data[[value_column]] <- log2(filtered_data[[value_column]] + 1)

# Compute medians for labeling
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", round(medians$Median, 2))
medians$Position <- medians$Median

# Group colors
group_colors <- c("HBC" = "#1f78b4", "CLD" = "#33a02c", "HCC" = "#e31a1c")
names(group_colors) <- n_counts$Label

# Create the plot
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label),
            size = 7, fontface = "bold", color = "black", vjust = -1.2) +
  stat_compare_means(comparisons = list(
    c(n_counts$Label[1], n_counts$Label[2]),
    c(n_counts$Label[2], n_counts$Label[3]),
    c(n_counts$Label[1], n_counts$Label[3])
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 1.5, color = "black"),
    axis.line.x.bottom = element_line(size = 1.5, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

# Save plot as 600 DPI TIFF
ggsave(
  filename = file.path("output", "Cypermethrin_Boxplot_600dpi.tiff"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 600,
  units = "in",
  device = "tiff",
  compression = "lzw"
)

# Optionally display in RStudio
print(p)


###############################################################################
###############################################################################
# FIGURE 2A: PILCRS11, PILCRS_CYP, PILCRS8 BOXPLOTS
###############################################################################
###############################################################################

################################## PILCRS₁₁ — Manuscript Boxplot ##############################################

# Load required libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rlang)

# Load your dataset (already in environment)
PILCRS_1 <- PILCRS_1  # Replace with actual file path if needed

# Define column names
group_column <- "tissuetypeORG"
value_column <- "PILCRSA"  # This is the log2-transformed input column
y_axis_label <- expression("PILCRS"[11])  # Subscripted y-axis label

# Prepare dataset
filtered_data <- PILCRS_1 %>%
  mutate(
    !!group_column := gsub("PC", "HBC", !!sym(group_column)),
    !!value_column := log2(as.numeric(!!sym(value_column)))
  ) %>%
  filter(!is.na(!!sym(group_column)) & !is.na(!!sym(value_column))) %>%
  filter(!!sym(group_column) %in% c("HBC", "CLD", "HCC")) %>%
  mutate(Group_Label = factor(
    case_when(
      !!sym(group_column) == "HBC" ~ "HBC (n=249)",
      !!sym(group_column) == "CLD" ~ "CLD (n=228)",
      !!sym(group_column) == "HCC" ~ "HCC (n=116)"
    ),
    levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")
  ))

# Calculate median values
medians <- filtered_data %>%
  group_by(Group_Label) %>%
  summarise(
    Median = median(!!sym(value_column), na.rm = TRUE),
    Position = quantile(!!sym(value_column), 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Label = paste0("Med: ", round(Median, 4)))

# Define custom color scheme
group_colors <- c(
  "HBC (n=249)" = "#1f78b4",
  "CLD (n=228)" = "#33a02c",
  "HCC (n=116)" = "#e31a1c"
)

# Create the plot
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column), fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5, inherit.aes = FALSE) +
  stat_compare_means(
    comparisons = list(
      c("HBC (n=249)", "CLD (n=228)"),
      c("CLD (n=228)", "HCC (n=116)"),
      c("HBC (n=249)", "HCC (n=116)")
    ),
    method = "wilcox.test",
    label = "p.signif",
    size = 7,
    fontface = "bold",
    color = "black",
    step.increase = 0.08,
    bracket.size = 1.2
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +
  coord_cartesian(ylim = c(
    min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
    max(filtered_data[[value_column]], na.rm = TRUE) * 1.3
  )) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 1.5, color = "black"),
    axis.line.x.bottom = element_line(size = 1.5, color = "black")
  )

# Save high-resolution TIFF (600 DPI)
ggsave(
  filename = file.path("output", "Figure2_PILCRS11.tiff"),
  plot = p,
  device = "tiff",
  dpi = 600,
  width = 7,
  height = 6,
  units = "in",
  compression = "lzw"
)

# Also print the plot to RStudio Viewer
print(p)















####################################### CYPERMETHRIN: PILCRS–CYP (log2) ##########################################

# Load required libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rlang)

# Load dataset
PILCRS_1 <- PILCRS_1  # Replace if loading externally

# Define variables
group_column <- "tissuetypeORG"
value_column <- "PILCRSCYP"

# Data preparation
filtered_data <- PILCRS_1 %>%
  mutate(
    !!group_column := gsub("PC", "HBC", !!sym(group_column)),
    !!value_column := log2(as.numeric(!!sym(value_column)))
  ) %>%
  filter(!is.na(!!sym(group_column)) & !is.na(!!sym(value_column))) %>%
  filter(!!sym(group_column) %in% c("HBC", "CLD", "HCC")) %>%
  mutate(Group_Label = factor(
    case_when(
      !!sym(group_column) == "HBC" ~ "HBC (n=249)",
      !!sym(group_column) == "CLD" ~ "CLD (n=228)",
      !!sym(group_column) == "HCC" ~ "HCC (n=116)"
    ),
    levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")
  ))

# Median annotation
medians <- filtered_data %>%
  group_by(Group_Label) %>%
  summarise(
    Median = median(!!sym(value_column), na.rm = TRUE),
    Position = quantile(!!sym(value_column), 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Label = paste0("Med: ", round(Median, 4)))

# Define color palette
group_colors <- c(
  "HBC (n=249)" = "#1f78b4",
  "CLD (n=228)" = "#33a02c",
  "HCC (n=116)" = "#e31a1c"
)

# Generate the boxplot
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5, inherit.aes = FALSE) +
  stat_compare_means(
    comparisons = list(
      c("HBC (n=249)", "CLD (n=228)"),
      c("CLD (n=228)", "HCC (n=116)"),
      c("HBC (n=249)", "HCC (n=116)")
    ),
    method = "wilcox.test",
    label = "p.signif",
    size = 7,
    fontface = "bold",
    color = "black",
    step.increase = 0.08,
    bracket.size = 1.2
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = NULL,
    y = expression("PILCRS"["CYP"]),
    title = ""
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 1.5, color = "black"),
    axis.line.x.bottom = element_line(size = 1.5, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

# Save high-resolution figure (optional)
ggsave(
  filename = file.path("output", "Figure2_PILCRS_CYP.tiff"),
  plot = p,
  device = "tiff",
  dpi = 600,
  width = 7,
  height = 6,
  units = "in",
  compression = "lzw"
)

# Display the plot
print(p)

















#################################################################################################
# Manuscript Figure: PILCRS₈ (log₂ + shift) Boxplot with Group Comparisons
#################################################################################################

# Load required libraries
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rlang)

# Define variables
group_column <- "tissuetypeORG"
value_column <- "PILCRSB_NEW"
y_axis_label <- expression(PILCRS[8])  # Subscript 8 for figure label

# Log-transform with shift to ensure positivity
shift <- abs(min(log2(as.numeric(PILCRS_1[[value_column]])), na.rm = TRUE)) + 1

# Data preparation
filtered_data <- PILCRS_1 %>%
  mutate(
    !!group_column := gsub("PC", "HBC", !!sym(group_column)),
    !!value_column := log2(as.numeric(!!sym(value_column))) + shift
  ) %>%
  filter(!is.na(!!sym(group_column)) & !is.na(!!sym(value_column))) %>%
  filter(!!sym(group_column) %in% c("HBC", "CLD", "HCC")) %>%
  mutate(Group_Label = factor(
    case_when(
      !!sym(group_column) == "HBC" ~ "HBC (n=249)",
      !!sym(group_column) == "CLD" ~ "CLD (n=228)",
      !!sym(group_column) == "HCC" ~ "HCC (n=116)"
    ),
    levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")
  ))

# Median text labels
medians <- filtered_data %>%
  group_by(Group_Label) %>%
  summarise(
    Median = median(!!sym(value_column), na.rm = TRUE),
    Position = quantile(!!sym(value_column), 0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Label = paste0("Med: ", round(Median, 4)))

# Group-specific colors
group_colors <- c(
  "HBC (n=249)" = "#1f78b4",
  "CLD (n=228)" = "#33a02c",
  "HCC (n=116)" = "#e31a1c"
)

# Generate the plot
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5, inherit.aes = FALSE) +
  stat_compare_means(
    comparisons = list(
      c("HBC (n=249)", "CLD (n=228)"),
      c("CLD (n=228)", "HCC (n=116)"),
      c("HBC (n=249)", "HCC (n=116)")
    ),
    method = "wilcox.test",
    label = "p.signif",
    size = 7,
    fontface = "bold",
    color = "black",
    step.increase = 0.08,
    bracket.size = 1.2
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = NULL,
    y = y_axis_label,
    title = NULL
  ) +
  coord_cartesian(ylim = c(3, 11.5)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 22, face = "bold", color = "black"),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 1.5, color = "black"),
    axis.line.x.bottom = element_line(size = 1.5, color = "black")
  )

# Print the plot
print(p)

# Optional: Save as high-resolution TIFF
ggsave(
  filename = file.path("output", "Figure2_PILCRS8.tiff"),
  plot = p,
  device = "tiff",
  dpi = 600,
  width = 7,
  height = 6,
  units = "in",
  compression = "lzw"
)


###############################################################################
###############################################################################
# FIGURE 2B: PREDICTED RISK SCORE (XGBOOST) BY GROUP
###############################################################################
###############################################################################

########################################################
# XGBoost Risk Prediction Using PLS₁₁ (PLS_GroupA) + Covariates
# Author: [Your Name]
# Purpose: Generate covariate-adjusted disease risk predictions 
#          using PLS₁₁ as main exposure variable.
########################################################

########################################################
# XGBoost Risk Prediction Using PLS₁₁ (PLS_GroupA) + Covariates
# Purpose: Generate covariate-adjusted disease risk predictions 
#          using PLS₁₁ as main exposure variable.
########################################################

# Load required libraries
library(readr)
library(ggplot2)
library(xgboost)
library(dplyr)
library(ggsignif)
library(FSA)
library(rstatix)

#-----------------------
# Load and prepare data
#-----------------------
df <- PILCRS_1_MASTER

# Filter for relevant disease groups
df <- df[df$tissuetypeORG %in% c("PC", "CLD", "HCC"), ]

# Recode group labels
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

#-----------------------
# Define covariates (PLS₁₁ = PLS_GroupA)
#-----------------------
covariates <- c("PLS_GroupA", "Gender_Cat", "Occupation new", "Alcoholcat", "Ageorg")
df <- df[complete.cases(df[, covariates]), ]

#-----------------------
# Create binary outcome
#-----------------------
df$RiskLabel <- ifelse(df$tissuetypeORG == "HBC (n=249)", 0, 1)

#-----------------------
# Train XGBoost model
#-----------------------
X <- as.matrix(df[, covariates])
y <- df$RiskLabel

dtrain <- xgb.DMatrix(data = X, label = y)

model <- xgb.train(
  params = list(
    objective     = "binary:logistic",
    max_depth     = 3,
    learning_rate = 0.1
  ),
  data    = dtrain,
  nrounds = 300,
  verbose = 0
)

# Predict probabilities
df$Prob_CLD_HCC <- predict(model, dtrain)

#-----------------------
# Statistical tests
#-----------------------
kruskal_result <- kruskal.test(Prob_CLD_HCC ~ tissuetypeORG, data = df)
print(kruskal_result)

dunn <- dunnTest(Prob_CLD_HCC ~ tissuetypeORG, data = df, method = "bonferroni")
dunn_res <- dunn$res
dunn_res$signif <- cut(dunn_res$P.adj,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns"))
print(dunn_res)

#-----------------------
# Generate plot
#-----------------------
p <- ggplot(df, aes(x = tissuetypeORG, y = Prob_CLD_HCC, fill = tissuetypeORG)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.75, size = 1.1) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.4, color = "black") +
  geom_signif(
    comparisons = list(
      c("HBC (n=249)", "CLD (n=228)"),
      c("CLD (n=228)", "HCC (n=116)"),
      c("HBC (n=249)", "HCC (n=116)")
    ),
    map_signif_level = function(p) {
      if (p >= 0.05) return("ns")
      else if (p < 0.001) return("***")
      else if (p < 0.01) return("**")
      else return("*")
    },
    test = "wilcox.test",
    step_increase = 0.08,
    textsize = 7
  ) +
  labs(
    x = NULL,
    y = "Predicted Risk Score",
    title = expression(bold("") * bold(PLS[11]) * bold(" + Covariates"))
  ) +
  scale_y_continuous(
    limits = c(0, 1.24),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x     = element_text(size = 20, face = "bold", color = "black"),
    axis.text.y     = element_text(size = 20, face = "bold", color = "black"),
    axis.line.x     = element_line(size = 1.5, color = "black"),
    axis.line.y     = element_line(size = 1.5, color = "black"),
    panel.grid.major.y = element_line(color = "gray80", size = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", size = 0.4),
    panel.grid.minor.x = element_blank(),
    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Print plot
print(p)
save_fig("Figure2B_XGBoost_PLS11", p, width = 7, height = 6)


#-----------------------
# End of script
#-----------------------





















########################################################
# XGBoost Risk Prediction Using PLS₈ (PLS_GroupB) + Covariates
# Purpose: Generate covariate-adjusted disease risk predictions 
#          using PLS₈ as main exposure variable.
########################################################

# Load required libraries
library(readr)
library(ggplot2)
library(xgboost)
library(dplyr)
library(ggsignif)
library(FSA)
library(rstatix)

#-----------------------
# Load and prepare data
#-----------------------
df <- PILCRS_1_MASTER

# Filter for relevant disease groups
df <- df[df$tissuetypeORG %in% c("PC", "CLD", "HCC"), ]

# Recode group labels
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

#-----------------------
# Define covariates (PLS₈ = PLS_GroupB)
#-----------------------
covariates <- c("PLS_GroupB", "Gender_Cat", "Occupation new", "Alcoholcat", "Ageorg")
df <- df[complete.cases(df[, covariates]), ]

#-----------------------
# Create binary outcome
#-----------------------
df$RiskLabel <- ifelse(df$tissuetypeORG == "HBC (n=249)", 0, 1)

#-----------------------
# Train XGBoost model
#-----------------------
X <- as.matrix(df[, covariates])
y <- df$RiskLabel

dtrain <- xgb.DMatrix(data = X, label = y)

model <- xgb.train(
  params = list(
    objective     = "binary:logistic",
    max_depth     = 3,
    learning_rate = 0.1
  ),
  data    = dtrain,
  nrounds = 300,
  verbose = 0
)

# Predict probabilities
df$Prob_CLD_HCC <- predict(model, dtrain)

#-----------------------
# Statistical tests
#-----------------------
kruskal_result <- kruskal.test(Prob_CLD_HCC ~ tissuetypeORG, data = df)
print(kruskal_result)

dunn <- dunnTest(Prob_CLD_HCC ~ tissuetypeORG, data = df, method = "bonferroni")
dunn_res <- dunn$res
dunn_res$signif <- cut(dunn_res$P.adj,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns"))
print(dunn_res)

#-----------------------
# Generate plot
#-----------------------
p <- ggplot(df, aes(x = tissuetypeORG, y = Prob_CLD_HCC, fill = tissuetypeORG)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.75, size = 1.1) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.4, color = "black") +
  geom_signif(
    comparisons = list(
      c("HBC (n=249)", "CLD (n=228)"),
      c("CLD (n=228)", "HCC (n=116)"),
      c("HBC (n=249)", "HCC (n=116)")
    ),
    map_signif_level = function(p) {
      if (p >= 0.05) return("ns")
      else if (p < 0.001) return("***")
      else if (p < 0.01) return("**")
      else return("*")
    },
    test = "wilcox.test",
    step_increase = 0.08,
    textsize = 7
  ) +
  labs(
    x = NULL,
    y = "Predicted Risk Score",
    title = expression(bold("") * bold(PLS[8]) * bold(" + Covariates"))
  ) +
  scale_y_continuous(
    limits = c(0, 1.24),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x     = element_text(size = 20, face = "bold", color = "black"),
    axis.text.y     = element_text(size = 20, face = "bold", color = "black"),
    axis.line.x     = element_line(size = 1.5, color = "black"),
    axis.line.y     = element_line(size = 1.5, color = "black"),
    panel.grid.major.y = element_line(color = "gray80", size = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", size = 0.4),
    panel.grid.minor.x = element_blank(),
    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Print plot
print(p)
save_fig("Figure2B_XGBoost_PLS8", p, width = 7, height = 6)


#-----------------------
# End of script
#-----------------------






















########################################################
# XGBoost Risk Prediction Using Cypermethrin_nM + Covariates
# Purpose: Predict CLD/HCC risk using Cypermethrin_nM exposure
#          adjusted for demographic and occupational covariates
########################################################

# Load required libraries
library(readr)
library(ggplot2)
library(xgboost)
library(dplyr)
library(ggsignif)
library(FSA)
library(rstatix)

#-----------------------
# Load and prepare data
#-----------------------
df <- PILCRS_1_MASTER

# Filter to PC, CLD, HCC only
df <- df[df$tissuetypeORG %in% c("PC", "CLD", "HCC"), ]

# Recode group labels
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

#-----------------------
# Define covariates
#-----------------------
covariates <- c("Cypermethrin_nM", "Gender_Cat", "Occupation new", "Alcoholcat", "Ageorg")
df <- df[complete.cases(df[, covariates]), ]

#-----------------------
# Create binary outcome
#-----------------------
df$RiskLabel <- ifelse(df$tissuetypeORG == "HBC (n=249)", 0, 1)

#-----------------------
# Train XGBoost model
#-----------------------
X <- as.matrix(df[, covariates])
y <- df$RiskLabel

dtrain <- xgb.DMatrix(data = X, label = y)

model <- xgb.train(
  params = list(
    objective     = "binary:logistic",
    max_depth     = 3,
    learning_rate = 0.1
  ),
  data    = dtrain,
  nrounds = 300,
  verbose = 0
)

# Predict probabilities
df$Prob_CLD_HCC <- predict(model, dtrain)

#-----------------------
# Statistical testing
#-----------------------
kruskal_result <- kruskal.test(Prob_CLD_HCC ~ tissuetypeORG, data = df)
print(kruskal_result)

dunn <- dunnTest(Prob_CLD_HCC ~ tissuetypeORG, data = df, method = "bonferroni")
dunn_res <- dunn$res
dunn_res$signif <- cut(dunn_res$P.adj,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns"))
print(dunn_res)

#-----------------------
# Generate plot
#-----------------------
p <- ggplot(df, aes(x = tissuetypeORG, y = Prob_CLD_HCC, fill = tissuetypeORG)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "black") +
  geom_signif(
    comparisons = list(
      c("HBC (n=249)", "CLD (n=228)"),
      c("CLD (n=228)", "HCC (n=116)"),
      c("HBC (n=249)", "HCC (n=116)")
    ),
    map_signif_level = function(p) {
      if (p >= 0.05) return("ns")
      else if (p < 0.001) return("***")
      else if (p < 0.01) return("**")
      else return("*")
    },
    test = "wilcox.test",
    step_increase = 0.08,
    textsize = 7
  ) +
  labs(
    x = NULL,
    y = "Predicted Risk Score",
    title = "Cypermethrin + Covariates"
  ) +
  scale_y_continuous(
    limits = c(0, 1.24),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray85", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 1.5, color = "black"),
    axis.line.y = element_line(size = 1.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 20, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 20, face = "bold", color = "black"),
    plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

print(p)
save_fig("Figure2B_XGBoost_Cypermethrin", p, width = 7, height = 6)


###############################################################################
###############################################################################
# FIGURE 2C: QUARTILE-STRATIFIED PREDICTED RISK PROBABILITY
###############################################################################
###############################################################################

#########################################################################
######################### For the Main Figure ###########################
#########################################################################

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})

# ---- Data prep ----
df <- Risk_Per_Calc %>%
  mutate(Status = factor(tissuetypeORG_DUPLICATE, levels = c("PC","CLD","HCC"))) %>%
  select(where(is.numeric), Status) %>%
  filter(complete.cases(.))

# ---- Colors (matching curves) ----
COLORS <- c("CLD" = "#33a02c",  # green
            "HCC" = "#e31a1c")  # red

# ---- Smoothing helper (monotone cubic) ----
smooth_q <- function(x, y, n = 200) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  s <- splinefun(x, y, method = "monoH.FC")
  tibble(
    x_s = seq(min(x), max(x), length.out = n),
    y_s = s(seq(min(x), max(x), length.out = n))
  )
}

# ---- Custom theme ----
theme_ref <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "gray85", size = 0.5),
      panel.grid.minor = element_blank(),
      axis.line = element_line(size = 1.2, color = "black"),
      axis.title.x = element_text(size = 22, face = "bold", color = "black"),
      axis.title.y = element_text(size = 22, face = "bold", color = "black"),
      axis.text.x  = element_text(size = 20, face = "bold", color = "black"),
      axis.text.y  = element_text(size = 20, face = "bold", color = "black"),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
      legend.position = "none"
    )
}

# ---- Plotting function with embedded legend ----
panel_pilcrs_curvy_CI <- function(data, exposure_var, exposure_label) {
  # Data prep
  dfx <- data %>%
    filter(!is.na(.data[[exposure_var]])) %>%
    mutate(SCORE = as.numeric(.data[[exposure_var]]),
           SCORE_quart = factor(ntile(SCORE, 4), labels = c("Q1","Q2","Q3","Q4")))
  
  newgrid <- tibble(SCORE_quart = factor(c("Q1","Q2","Q3","Q4"),
                                         levels = c("Q1","Q2","Q3","Q4")))
  
  # ---- CLD vs PC ----
  cld_dat <- dfx %>% filter(Status %in% c("PC","CLD")) %>%
    mutate(Case = as.integer(Status == "CLD"))
  fit_cld <- glm(Case ~ SCORE_quart, data = cld_dat, family = binomial())
  prc <- predict(fit_cld, newdata = newgrid, type = "link", se.fit = TRUE)
  cld_q <- newgrid %>%
    mutate(x = as.numeric(SCORE_quart),
           p_hat = plogis(prc$fit) * 100,
           lo = plogis(prc$fit - 1.96*prc$se.fit) * 100,
           hi = plogis(prc$fit + 1.96*prc$se.fit) * 100,
           Outcome = "CLD")
  
  # ---- HCC vs PC ----
  hcc_dat <- dfx %>% filter(Status %in% c("PC","HCC")) %>%
    mutate(Case = as.integer(Status == "HCC"))
  fit_hcc <- glm(Case ~ SCORE_quart, data = hcc_dat, family = binomial())
  prh <- predict(fit_hcc, newdata = newgrid, type = "link", se.fit = TRUE)
  hcc_q <- newgrid %>%
    mutate(x = as.numeric(SCORE_quart),
           p_hat = plogis(prh$fit) * 100,
           lo = plogis(prh$fit - 1.96*prh$se.fit) * 100,
           hi = plogis(prh$fit + 1.96*prh$se.fit) * 100,
           Outcome = "HCC")
  
  # Combine & enforce monotone
  qpts <- bind_rows(cld_q, hcc_q) %>%
    group_by(Outcome) %>%
    mutate(p_hat = cummax(p_hat),
           lo = cummax(lo),
           hi = cummax(hi)) %>%
    ungroup()
  
  # Smooth curves
  smoothed <- qpts %>%
    group_by(Outcome) %>%
    group_modify(~{
      tibble(
        x_s = smooth_q(.x$x, .x$p_hat)$x_s,
        p_s = smooth_q(.x$x, .x$p_hat)$y_s,
        lo_s = smooth_q(.x$x, .x$lo)$y_s,
        hi_s = smooth_q(.x$x, .x$hi)$y_s
      )
    }) %>% ungroup()
  
  # ---- Build plot ----
  p <- ggplot(smoothed, aes(x = x_s, y = p_s, color = Outcome, fill = Outcome)) +
    geom_ribbon(aes(ymin = lo_s, ymax = hi_s), alpha = 0.25) +
    geom_line(size = 1.8, lineend = "round") +
    scale_x_continuous(breaks = 1:4, labels = c("Q1","Q2","Q3","Q4")) +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
    scale_color_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    labs(x = "Quartile", y = "Predicted Risk Probability (%)") +
    ggtitle(exposure_label) +
    theme_ref() +
    
    # ---- Embedded legend with clear space ----
  annotate("rect", xmin = 1.05, xmax = 1.3, ymin = 95, ymax = 96.5,
           fill = alpha(COLORS["CLD"], 0.2), color = NA) +
    annotate("segment", x = 1.07, xend = 1.28, y = 95.75, yend = 95.75,
             color = COLORS["CLD"], size = 3, lineend = "round") +
    annotate("text", x = 1.45, y = 95.75, label = "CLD vs PC",
             hjust = 0, size = 6, fontface = "bold", color = "black") +
    
    annotate("rect", xmin = 1.05, xmax = 1.3, ymin = 90, ymax = 91.5,
             fill = alpha(COLORS["HCC"], 0.2), color = NA) +
    annotate("segment", x = 1.07, xend = 1.28, y = 90.75, yend = 90.75,
             color = COLORS["HCC"], size = 3, lineend = "round") +
    annotate("text", x = 1.45, y = 90.75, label = "HCC vs PC",
             hjust = 0, size = 6, fontface = "bold", color = "black")
  
  return(p)
}

# ---- Generate Panels ----
pA <- panel_pilcrs_curvy_CI(df, "PILCRSA", expression(bold(PILCRS[11])))
pB <- panel_pilcrs_curvy_CI(df, "PILCRSB_NEW", expression(bold(PILCRS[8])))
pC <- panel_pilcrs_curvy_CI(df, "PILCRSCYP", expression(bold(PILCRS[CYP])))

# Print panels
print(pA)
print(pB)
print(pC)
save_fig("Figure2C_PILCRS11", pA, width = 7, height = 6)
save_fig("Figure2C_PILCRS8", pB, width = 7, height = 6)
save_fig("Figure2C_PILCRS_CYP", pC, width = 7, height = 6)

# Combined vertical layout
final_plot <- pA / pB / pC
print(final_plot)


###############################################################################
###############################################################################
# FIGURE 3 + SUPPLEMENTARY FIGURES S3 AND S4: BOOTSTRAP INTERNAL VALIDATION AND CALIBRATION
###############################################################################
###############################################################################

while (!is.null(dev.list())) dev.off()   # close any leftover file device so plots show in RStudio

library(xgboost)
# Works with both old (1.x) and new (3.x) xgboost versions
fit_xgb <- function(X, y) {
  xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1),
            data = xgb.DMatrix(X, label = y), nrounds = 300, verbose = 0)
}


# ==========================#######################
# INTERNAL VALIDATION: PILCRS - PLS11 + Covariates (CLD and HCC separately)
# ==========================#######################

# Load libraries
library(xgboost)
library(pROC)
library(boot)
library(ggplot2)
library(dplyr)
library(patchwork)

# -------------------------
# STEP 1: Prepare Data
# -------------------------
df <- PILCRS_1_MASTER
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

# -------------------------
# STEP 2: Define Bootstrap Validation Function
# -------------------------
run_validation <- function(df, outcome_label, save_prefix) {
  df <- df[df$tissuetypeORG %in% c("HBC (n=249)", outcome_label), ]
  covariates <- c("PLS_GroupA", "Gender_Cat", "Alcoholcat", "Ageorg", "Occupation new")
  df <- df[complete.cases(df[, covariates]), ]
  df$RiskLabel <- ifelse(df$tissuetypeORG == "HBC (n=249)", 0, 1)
  X_full <- data.matrix(df[, covariates])
  y_full <- df$RiskLabel
  
  # STEP 3: Bootstrap
  bootstrap_metrics <- function(data, indices) {
    boot_data <- data[indices, ]
    X_boot <- data.matrix(boot_data[, covariates])
    y_boot <- boot_data$RiskLabel
    model_boot <- fit_xgb(X_boot, y_boot)
    pred_boot <- predict(model_boot, X_boot)
    pred_test <- predict(model_boot, X_full)
    auc_app <- auc(y_boot, pred_boot)
    auc_test <- auc(y_full, pred_test)
    logit_test <- log(pred_test / (1 - pred_test))
    calib_model <- glm(y_full ~ logit_test, family = binomial)
    calib_slope <- coef(calib_model)[2]
    return(c(auc_app, auc_test, calib_slope))
  }
  
  set.seed(2025)
  results <- boot(data = df, statistic = bootstrap_metrics, R = 1000)
  auc_app <- results$t[, 1]
  auc_test <- results$t[, 2]
  calib_slopes <- results$t[, 3]
  
  # STEP 4: Final Model
  model_orig <- fit_xgb(X_full, y_full)
  pred_orig <- predict(model_orig, X_full)
  auc_orig <- auc(y_full, pred_orig)
  auc_corrected <- auc_orig - mean(auc_app - auc_test)
  slope_corrected <- mean(calib_slopes)
  ci_slope <- quantile(calib_slopes, probs = c(0.025, 0.975))
  
  # STEP 5: Output
  cat(paste0("\n==== RESULTS FOR MANUSCRIPT: ", outcome_label, " ===="))
  cat("\nOriginal AUC: ", round(auc_orig, 3))
  cat("\nOptimism-corrected AUC: ", round(auc_corrected, 3))
  cat("\nCalibration slope: ", round(slope_corrected, 3),
      " (95% CI: ", round(ci_slope[1], 3), "–", round(ci_slope[2], 3), ")")
  cat("\nValidation: 1000x Bootstrap\n")
  
  # STEP 6A: AUC Histogram
  fig_auc <- ggplot(data.frame(AUC = auc_test), aes(x = AUC)) +
    geom_histogram(fill = "skyblue", color = "black", bins = 30) +
    geom_vline(xintercept = mean(auc_test), linetype = "dashed", color = "red", size = 1.2) +
    labs(title = paste("AUC Distribution:", outcome_label), x = "AUC", y = "Frequency") +
    theme_minimal(base_size = 14)
  
  print(fig_auc)
  save_fig(paste0("FigureS4_AUC_", save_prefix), fig_auc, width = 7, height = 5)
  
  # STEP 6B: Slope Histogram
  fig_slope <- ggplot(data.frame(Slope = calib_slopes), aes(x = Slope)) +
    geom_histogram(fill = "lightgreen", color = "black", bins = 30) +
    geom_vline(xintercept = slope_corrected, linetype = "dashed", color = "darkgreen", size = 1.2) +
    labs(title = paste("Calibration Slope:", outcome_label), x = "Calibration Slope", y = "Frequency") +
    theme_minimal(base_size = 14)
  
  print(fig_slope)
  save_fig(paste0("FigureS3_CalibrationSlope_", save_prefix), fig_slope, width = 7, height = 5)
  
  # STEP 7: Calibration Plot
  df$PredictedProb <- pred_orig
  df_cal <- df %>% mutate(bin = ntile(PredictedProb, 10)) %>%
    group_by(bin) %>% summarise(mean_pred = mean(PredictedProb), mean_obs = mean(RiskLabel), n = n())
  
  annot_text <- paste0("AUC (bootstrap): ", round(auc_corrected, 3),
                       "\nSlope: ", round(slope_corrected, 2), " (95% CI: ",
                       round(ci_slope[1], 2), "–", round(ci_slope[2], 2), ")\nValidation: 1000x bootstrap")
  
  fig_calib <- ggplot(df_cal, aes(x = mean_pred, y = mean_obs)) +
    geom_point(aes(shape = "Observed"), size = 3, color = "black") +
    geom_smooth(aes(linetype = "LOESS"), method = "loess", se = TRUE, span = 1, color = "blue", linewidth = 1) +
    geom_abline(aes(linetype = "Ideal"), slope = 1, intercept = 0, color = "gray30") +
    annotate("text", x = 0.05, y = 0.95, label = annot_text, hjust = 0, vjust = 1, size = 5) +
    scale_linetype_manual(values = c("LOESS" = "solid", "Ideal" = "dashed")) +
    scale_shape_manual(values = c("Observed" = 16)) +
    labs(title = paste("Calibration Plot: PILCRS11 + Covariates -", outcome_label),
         x = "Mean predicted probability", y = "Observed probability", linetype = NULL, shape = NULL) +
    theme_minimal(base_size = 15) +
    theme(axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "right")
  
  print(fig_calib)
  save_fig(paste0("Figure3_Calibration_", save_prefix), fig_calib, width = 8, height = 6)
  
  return(fig_calib)
}

# -------------------------
# FINAL: Run Both Models and Save Combined Panel
# -------------------------
p_cld <- run_validation(df, "CLD (n=228)", "PILCRS11_CLD")
p_hcc <- run_validation(df, "HCC (n=116)", "PILCRS11_HCC")

# Show combined panel
combined_plot <- p_cld + p_hcc
print(combined_plot)
















# ==========================#######################
# INTERNAL VALIDATION: PILCRS - PLS₈ + Covariates (CLD and HCC separately)
# ==========================#######################

# Load libraries
library(xgboost)
library(pROC)
library(boot)
library(ggplot2)
library(dplyr)
library(patchwork)

# -------------------------
# STEP 1: Prepare Data
# -------------------------
df <- PILCRS_1_MASTER
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

# -------------------------
# STEP 2: Define Bootstrap Validation Function
# -------------------------
run_validation <- function(df, outcome_label, save_prefix) {
  df <- df[df$tissuetypeORG %in% c("HBC (n=249)", outcome_label), ]
  covariates <- c("PLS_GroupB", "Gender_Cat", "Alcoholcat", "Ageorg", "Occupation new")
  df <- df[complete.cases(df[, covariates]), ]
  df$RiskLabel <- ifelse(df$tissuetypeORG == "HBC (n=249)", 0, 1)
  X_full <- data.matrix(df[, covariates])
  y_full <- df$RiskLabel
  
  # STEP 3: Bootstrap
  bootstrap_metrics <- function(data, indices) {
    boot_data <- data[indices, ]
    X_boot <- data.matrix(boot_data[, covariates])
    y_boot <- boot_data$RiskLabel
    model_boot <- fit_xgb(X_boot, y_boot)
    pred_boot <- predict(model_boot, X_boot)
    pred_test <- predict(model_boot, X_full)
    auc_app <- auc(y_boot, pred_boot)
    auc_test <- auc(y_full, pred_test)
    logit_test <- log(pred_test / (1 - pred_test))
    calib_model <- glm(y_full ~ logit_test, family = binomial)
    calib_slope <- coef(calib_model)[2]
    return(c(auc_app, auc_test, calib_slope))
  }
  
  set.seed(2025)
  results <- boot(data = df, statistic = bootstrap_metrics, R = 1000)
  auc_app <- results$t[, 1]
  auc_test <- results$t[, 2]
  calib_slopes <- results$t[, 3]
  
  # STEP 4: Final Model
  model_orig <- fit_xgb(X_full, y_full)
  pred_orig <- predict(model_orig, X_full)
  auc_orig <- auc(y_full, pred_orig)
  auc_corrected <- auc_orig - mean(auc_app - auc_test)
  slope_corrected <- mean(calib_slopes)
  ci_slope <- quantile(calib_slopes, probs = c(0.025, 0.975))
  
  # STEP 5: Output
  cat(paste0("\n==== RESULTS FOR MANUSCRIPT (PLS₈ + Covariates): ", outcome_label, " ===="))
  cat("\nOriginal AUC: ", round(auc_orig, 3))
  cat("\nOptimism-corrected AUC: ", round(auc_corrected, 3))
  cat("\nCalibration slope: ", round(slope_corrected, 3),
      " (95% CI: ", round(ci_slope[1], 3), "–", round(ci_slope[2], 3), ")")
  cat("\nValidation: 1000x Bootstrap\n")
  
  # STEP 6A: AUC Histogram
  fig_auc <- ggplot(data.frame(AUC = auc_test), aes(x = AUC)) +
    geom_histogram(fill = "skyblue", color = "black", bins = 30) +
    geom_vline(xintercept = mean(auc_test), linetype = "dashed", color = "red", size = 1.2) +
    labs(title = paste("AUC Distribution:", outcome_label), x = "AUC", y = "Frequency") +
    theme_minimal(base_size = 14)
  
  print(fig_auc)
  save_fig(paste0("FigureS4_AUC_", save_prefix), fig_auc, width = 7, height = 5)
  
  # STEP 6B: Slope Histogram
  fig_slope <- ggplot(data.frame(Slope = calib_slopes), aes(x = Slope)) +
    geom_histogram(fill = "lightgreen", color = "black", bins = 30) +
    geom_vline(xintercept = slope_corrected, linetype = "dashed", color = "darkgreen", size = 1.2) +
    labs(title = paste("Calibration Slope:", outcome_label), x = "Calibration Slope", y = "Frequency") +
    theme_minimal(base_size = 14)
  
  print(fig_slope)
  save_fig(paste0("FigureS3_CalibrationSlope_", save_prefix), fig_slope, width = 7, height = 5)
  
  # STEP 7: Calibration Plot
  df$PredictedProb <- pred_orig
  df_cal <- df %>% mutate(bin = ntile(PredictedProb, 10)) %>%
    group_by(bin) %>% summarise(mean_pred = mean(PredictedProb), mean_obs = mean(RiskLabel), n = n())
  
  annot_text <- paste0("AUC (bootstrap): ", round(auc_corrected, 3),
                       "\nSlope: ", round(slope_corrected, 2), " (95% CI: ",
                       round(ci_slope[1], 2), "–", round(ci_slope[2], 2), ")\nValidation: 1000x bootstrap")
  
  fig_calib <- ggplot(df_cal, aes(x = mean_pred, y = mean_obs)) +
    geom_point(aes(shape = "Observed"), size = 3, color = "black") +
    geom_smooth(aes(linetype = "LOESS"), method = "loess", se = TRUE, span = 1, color = "blue", linewidth = 1) +
    geom_abline(aes(linetype = "Ideal"), slope = 1, intercept = 0, color = "gray30") +
    annotate("text", x = 0.05, y = 0.95, label = annot_text, hjust = 0, vjust = 1, size = 5) +
    scale_linetype_manual(values = c("LOESS" = "solid", "Ideal" = "dashed")) +
    scale_shape_manual(values = c("Observed" = 16)) +
    labs(title = paste("Calibration Plot: PILCRS (PLS₈ + Covariates) -", outcome_label),
         x = "Mean predicted probability", y = "Observed probability", linetype = NULL, shape = NULL) +
    theme_minimal(base_size = 15) +
    theme(axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "right")
  
  print(fig_calib)
  save_fig(paste0("Figure3_Calibration_", save_prefix), fig_calib, width = 8, height = 6)
  
  return(fig_calib)
}

# -------------------------
# FINAL: Run Both Models and Save Combined Panel
# -------------------------
p_cld <- run_validation(df, "CLD (n=228)", "PILCRS8_CLD")
p_hcc <- run_validation(df, "HCC (n=116)", "PILCRS8_HCC")

# Show combined panel
print(p_cld + p_hcc)



















# ==========================#######################
# INTERNAL VALIDATION: Cypermethrin_nM + Covariates (CLD and HCC separately)
# ==========================#######################

# Load libraries
library(xgboost)
library(pROC)
library(boot)
library(ggplot2)
library(dplyr)
library(patchwork)

# -------------------------
# STEP 1: Prepare Data
# -------------------------
df <- PILCRS_1_MASTER
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

# -------------------------
# STEP 2: Define Bootstrap Validation Function
# -------------------------
run_validation <- function(df, outcome_label, save_prefix) {
  df <- df[df$tissuetypeORG %in% c("HBC (n=249)", outcome_label), ]
  covariates <- c("Cypermethrin_nM", "Gender_Cat", "Alcoholcat", "Ageorg", "Occupation new")
  df <- df[complete.cases(df[, covariates]), ]
  df$RiskLabel <- ifelse(df$tissuetypeORG == "HBC (n=249)", 0, 1)
  X_full <- data.matrix(df[, covariates])
  y_full <- df$RiskLabel
  
  # STEP 3: Bootstrap
  bootstrap_metrics <- function(data, indices) {
    boot_data <- data[indices, ]
    X_boot <- data.matrix(boot_data[, covariates])
    y_boot <- boot_data$RiskLabel
    model_boot <- fit_xgb(X_boot, y_boot)
    pred_boot <- predict(model_boot, X_boot)
    pred_test <- predict(model_boot, X_full)
    auc_app <- auc(y_boot, pred_boot)
    auc_test <- auc(y_full, pred_test)
    logit_test <- log(pred_test / (1 - pred_test))
    calib_model <- glm(y_full ~ logit_test, family = binomial)
    calib_slope <- coef(calib_model)[2]
    return(c(auc_app, auc_test, calib_slope))
  }
  
  set.seed(2025)
  results <- boot(data = df, statistic = bootstrap_metrics, R = 1000)
  auc_app <- results$t[, 1]
  auc_test <- results$t[, 2]
  calib_slopes <- results$t[, 3]
  
  # STEP 4: Final Model
  model_orig <- fit_xgb(X_full, y_full)
  pred_orig <- predict(model_orig, X_full)
  auc_orig <- auc(y_full, pred_orig)
  auc_corrected <- auc_orig - mean(auc_app - auc_test)
  slope_corrected <- mean(calib_slopes)
  ci_slope <- quantile(calib_slopes, probs = c(0.025, 0.975))
  
  # STEP 5: Output
  cat(paste0("\n==== RESULTS FOR MANUSCRIPT (Cypermethrin + Covariates): ", outcome_label, " ===="))
  cat("\nOriginal AUC: ", round(auc_orig, 3))
  cat("\nOptimism-corrected AUC: ", round(auc_corrected, 3))
  cat("\nCalibration slope: ", round(slope_corrected, 3),
      " (95% CI: ", round(ci_slope[1], 3), "–", round(ci_slope[2], 3), ")")
  cat("\nValidation: 1000x Bootstrap\n")
  
  # STEP 6A: AUC Histogram
  fig_auc <- ggplot(data.frame(AUC = auc_test), aes(x = AUC)) +
    geom_histogram(fill = "skyblue", color = "black", bins = 30) +
    geom_vline(xintercept = mean(auc_test), linetype = "dashed", color = "red", size = 1.2) +
    labs(title = paste("AUC Distribution:", outcome_label), x = "AUC", y = "Frequency") +
    theme_minimal(base_size = 14)
  
  print(fig_auc)
  save_fig(paste0("FigureS4_AUC_", save_prefix), fig_auc, width = 7, height = 5)
  
  # STEP 6B: Slope Histogram
  fig_slope <- ggplot(data.frame(Slope = calib_slopes), aes(x = Slope)) +
    geom_histogram(fill = "lightgreen", color = "black", bins = 30) +
    geom_vline(xintercept = slope_corrected, linetype = "dashed", color = "darkgreen", size = 1.2) +
    labs(title = paste("Calibration Slope:", outcome_label), x = "Calibration Slope", y = "Frequency") +
    theme_minimal(base_size = 14)
  
  print(fig_slope)
  save_fig(paste0("FigureS3_CalibrationSlope_", save_prefix), fig_slope, width = 7, height = 5)
  
  # STEP 7: Calibration Plot
  df$PredictedProb <- pred_orig
  df_cal <- df %>% mutate(bin = ntile(PredictedProb, 10)) %>%
    group_by(bin) %>% summarise(mean_pred = mean(PredictedProb), mean_obs = mean(RiskLabel), n = n())
  
  annot_text <- paste0("AUC (bootstrap): ", round(auc_corrected, 3),
                       "\nSlope: ", round(slope_corrected, 2), " (95% CI: ",
                       round(ci_slope[1], 2), "–", round(ci_slope[2], 2), ")\nValidation: 1000x bootstrap")
  
  fig_calib <- ggplot(df_cal, aes(x = mean_pred, y = mean_obs)) +
    geom_point(aes(shape = "Observed"), size = 3, color = "black") +
    geom_smooth(aes(linetype = "LOESS"), method = "loess", se = TRUE, span = 1, color = "blue", linewidth = 1) +
    geom_abline(aes(linetype = "Ideal"), slope = 1, intercept = 0, color = "gray30") +
    annotate("text", x = 0.05, y = 0.95, label = annot_text, hjust = 0, vjust = 1, size = 5) +
    scale_linetype_manual(values = c("LOESS" = "solid", "Ideal" = "dashed")) +
    scale_shape_manual(values = c("Observed" = 16)) +
    labs(title = paste("Calibration Plot: Cypermethrin + Covariates -", outcome_label),
         x = "Mean predicted probability", y = "Observed probability", linetype = NULL, shape = NULL) +
    theme_minimal(base_size = 15) +
    theme(axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          legend.position = "right")
  
  print(fig_calib)
  save_fig(paste0("Figure3_Calibration_", save_prefix), fig_calib, width = 8, height = 6)
  
  return(fig_calib)
}

# -------------------------
# FINAL: Run Both Models and Save Combined Panel
# -------------------------
p_cyp_cld <- run_validation(df, "CLD (n=228)", "Cypermethrin_CLD")
p_cyp_hcc <- run_validation(df, "HCC (n=116)", "Cypermethrin_HCC")

# Show combined panel
print(p_cyp_cld + p_cyp_hcc)


###############################################################################
###############################################################################
# SUPPLEMENTARY FIGURE S1: INDIVIDUAL URINARY PESTICIDE BOXPLOTS
###############################################################################
###############################################################################

############ FINAL REVISED: Pendimethalin (nM) with Fixed Sample Sizes ############

# Load required libraries
library(ggplot2)
library(ggpubr)

# Parameters
data_file <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column <- "tissuetype"
value_column <- "Pendimethalin_Concnormcreatinine"
y_axis_label <- "Pendimethalin (nM)"

# Load data
filtered_data <- get(data_file)

# Replace PC with HBC
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

# Clean and filter
filtered_data <- filtered_data[!is.na(filtered_data[[group_column]]) & !is.na(filtered_data[[value_column]]), ]
group_levels <- c("HBC", "CLD", "HCC")
filtered_data <- filtered_data[filtered_data[[group_column]] %in% group_levels, ]
filtered_data[[group_column]] <- factor(filtered_data[[group_column]], levels = group_levels)

# Fixed group labels with manual sample sizes
n_counts <- data.frame(
  Group = group_levels,
  n = c(249, 228, 116)
)
n_counts$Label <- paste0(n_counts$Group, " (n=", n_counts$n, ")")
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = n_counts$Group,
                                    labels = n_counts$Label)

# Calculate medians
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", round(medians$Median, 4))
medians$Position <- medians$Median

# Set group colors
group_colors <- c("HBC" = "#1f78b4", "CLD" = "#33a02c", "HCC" = "#e31a1c")
names(group_colors) <- n_counts$Label

# Create plot
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label), 
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c(n_counts$Label[1], n_counts$Label[2]),
    c(n_counts$Label[2], n_counts$Label[3]),
    c(n_counts$Label[1], n_counts$Label[3])
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

# Print plot
print(p)
save_fig("FigureS1_D_Pendimethalin", p, width = 7, height = 5)

















########################################################
##### Metsulfuron-methyl (nM) Plot — No Transformation #
################2A#######################################################
############ FINAL REVISED: Metsulfuron-methyl (nM) — No Transformation ############

# Load required libraries
library(ggplot2)
library(ggpubr)

#=========================#
# Define data and columns
#=========================#
data_file <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column  <- "tissuetype"
value_column  <- "metsulfuronmethyl_Conc.normcreatinine"
y_axis_label  <- "Metsulfuron-methyl (nM)"

#=========================#
# Load and Prepare Data
#=========================#
filtered_data <- get(data_file)
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

# Filter to valid tissue types and non-NA values
filtered_data <- filtered_data[
  !is.na(filtered_data[[group_column]]) &
    !is.na(filtered_data[[value_column]]) &
    filtered_data[[group_column]] %in% c("HBC", "CLD", "HCC"), ]

# Set factor levels
group_levels <- c("HBC", "CLD", "HCC")
filtered_data[[group_column]] <- factor(filtered_data[[group_column]], levels = group_levels)

#=========================#
# Manually Set Sample Sizes
#=========================#
n_counts <- data.frame(
  Group = group_levels,
  n = c(249, 228, 116)
)
n_counts$Label <- paste0(n_counts$Group, " (n=", n_counts$n, ")")
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = n_counts$Group,
                                    labels = n_counts$Label)

#=========================#
# Calculate Medians
#=========================#
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", format(round(medians$Median, 4), nsmall = 4))
medians$Position <- medians$Median

#=========================#
# Define group colors
#=========================#
group_colors <- c("HBC" = "#1f78b4", "CLD" = "#33a02c", "HCC" = "#e31a1c")
names(group_colors) <- n_counts$Label

#=========================#
# Create the plot
#=========================#
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column), 
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label), 
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c(n_counts$Label[1], n_counts$Label[2]),
    c(n_counts$Label[2], n_counts$Label[3]),
    c(n_counts$Label[1], n_counts$Label[3])
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

#=========================#
# Print the plot
#=========================#
print(p)







#########2B##############################################
########################################################
##### Metsulfuron-methyl (nM) — FINAL PLOT w/ N ########
########################################################

# Load required libraries
library(ggplot2)
library(ggsignif)
library(dplyr)
library(rstatix)

#=========================#
# Define Parameters
#=========================#
data_file       <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column    <- "tissuetype"
value_column    <- "metsulfuronmethyl_Conc.normcreatinine"
y_axis_label    <- "Metsulfuron-methyl (nM)"
y_axis_limits   <- c(0, 0.35)

#=========================#
# Load and Prepare Data
#=========================#
df <- get(data_file)
df[[group_column]] <- gsub("PC", "HBC", df[[group_column]])

df <- df %>%
  filter(!is.na(.data[[group_column]]),
         !is.na(.data[[value_column]]),
         .data[[value_column]] > 0,
         .data[[group_column]] %in% c("HBC", "CLD", "HCC")) %>%
  mutate(
    group = factor(.data[[group_column]],
                   levels = c("HBC", "CLD", "HCC"),
                   labels = c("HBC (n = 249)", "CLD (n = 228)", "HCC (n = 116)"))
  )

#=========================#
# Compute Medians
#=========================#
medians <- df %>%
  group_by(group) %>%
  summarise(Median = median(.data[[value_column]], na.rm = TRUE)) %>%
  mutate(Label = paste0("Med: ", formatC(Median, format = "f", digits = 4)),
         Position = Median)

#=========================#
# Wilcoxon + Significance
#=========================#
stat_test <- wilcox_test(df, as.formula(paste(value_column, "~ group"))) %>%
  add_significance("p") %>%
  arrange(group1, group2) %>%
  mutate(
    p_label = p.signif,
    y.position = c(0.29, 0.25, 0.33)
  )

colnames(stat_test)[colnames(stat_test) == "group1"] <- "group"
colnames(stat_test)[colnames(stat_test) == "group2"] <- "group_to"

#=========================#
# Colors for Plot
#=========================#
group_colors <- c("HBC (n = 249)" = "#1f78b4",
                  "CLD (n = 228)" = "#33a02c",
                  "HCC (n = 116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(df, aes(x = group, y = .data[[value_column]], fill = group, color = group)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = group, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  geom_signif(data = stat_test,
              aes(xmin = group,
                  xmax = group_to,
                  annotations = p_label,
                  y_position = y.position),
              manual = TRUE,
              textsize = 6,
              size = 1.2,
              tip_length = 0.012,
              color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = y_axis_limits) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 22, face = "bold", color = "black"),
    axis.text       = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  )

#=========================#
# Display the Plot
#=========================#
print(p)
save_fig("FigureS1_A_Metsulfuron_methyl", p, width = 7, height = 5)









########## 3 ####################################################
##############################################################
## FINAL: 2,4-Dichlorophenoxyacetic acid (nM) Boxplot Script
##############################################################

# Load required libraries
library(ggplot2)
library(ggpubr)

#=========================#
# Define Parameters
#=========================#
data_file <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column  <- "tissuetype"
value_column  <- "dichlorophenoxyaceticacid_conc_normcreatininecolourimetric"
y_axis_label  <- "2,4-Dichlorophenoxyacetic acid (nM)"

#=========================#
# Load and Prepare Data
#=========================#
filtered_data <- get(data_file)
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

# Filter relevant data
filtered_data <- filtered_data[
  !is.na(filtered_data[[group_column]]) &
    !is.na(filtered_data[[value_column]]) &
    filtered_data[[group_column]] %in% c("HBC", "CLD", "HCC"), ]

#=========================#
# Set custom group labels
#=========================#
group_levels <- c("HBC", "CLD", "HCC")
group_labels <- c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = group_levels,
                                    labels = group_labels)

#=========================#
# Compute Medians
#=========================#
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", format(round(medians$Median, 4), nsmall = 4))
medians$Position <- medians$Median

#=========================#
# Define group colors
#=========================#
group_colors <- c("HBC (n=249)" = "#1f78b4",
                  "CLD (n=228)" = "#33a02c",
                  "HCC (n=116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians, aes(x = Group_Label, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c("HBC (n=249)", "CLD (n=228)"),
    c("CLD (n=228)", "HCC (n=116)"),
    c("HBC (n=249)", "HCC (n=116)")
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +  # <---- x-axis label removed
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),  # <---- suppress x-axis label
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

#=========================#
# Print the Plot
#=========================#
print(p)
save_fig("FigureS1_E_2_4_D", p, width = 7, height = 5)










###########4A##############################################

########################################################
## FINAL: Flocoumafen (nM) Plot with Sample Size Labels
########################################################

# Load required libraries
library(ggplot2)
library(ggpubr)

#=========================#
# Define Parameters
#=========================#
filtered_data <- Pesticides_Northern_hospital_samples_Madja_LHC_1
group_column  <- "tissuetype"
value_column  <- "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine"
y_axis_label  <- "Flocoumafen (nM)"

#=========================#
# Prepare and Filter Data
#=========================#
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

filtered_data <- filtered_data[
  !is.na(filtered_data[[group_column]]) &
    !is.na(filtered_data[[value_column]]) &
    filtered_data[[group_column]] %in% c("HBC", "CLD", "HCC"), ]

#=========================#
# Set custom group labels
#=========================#
group_levels <- c("HBC", "CLD", "HCC")
group_labels <- c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = group_levels,
                                    labels = group_labels)

#=========================#
# Compute Medians
#=========================#
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", format(round(medians$Median, 4), nsmall = 4))
medians$Position <- medians$Median

#=========================#
# Define Colors
#=========================#
group_colors <- c("HBC (n=249)" = "#1f78b4",
                  "CLD (n=228)" = "#33a02c",
                  "HCC (n=116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = Group_Label, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c("HBC (n=249)", "CLD (n=228)"),
    c("CLD (n=228)", "HCC (n=116)"),
    c("HBC (n=249)", "HCC (n=116)")
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +  # Removed x-axis label
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),  # Suppress x-axis title
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  )

#=========================#
# Print the Plot
#=========================#
print(p)












##########4B########################################


########################################################
## FINAL: Flocoumafen (nM) Plot with Group Sample Labels
########################################################

# Load required libraries
library(ggplot2)
library(ggpubr)

#=========================#
# Define and prepare data
#=========================#
filtered_data <- Pesticides_Northern_hospital_samples_Madja_LHC_1
group_column  <- "tissuetype"
value_column  <- "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine"
y_axis_label  <- "Flocoumafen (nM)"

# Standardize labels
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

# Keep only valid values and tissue groups
filtered_data <- filtered_data[
  !is.na(filtered_data[[group_column]]) &
    !is.na(filtered_data[[value_column]]) &
    filtered_data[[group_column]] %in% c("HBC", "CLD", "HCC"), ]

# Define group label with sample size
group_levels <- c("HBC", "CLD", "HCC")
group_labels <- c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = group_levels,
                                    labels = group_labels)

#=========================#
# Compute Medians
#=========================#
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", format(round(medians$Median, 4), nsmall = 4))
medians$Position <- medians$Median

#=========================#
# Define Colors
#=========================#
group_colors <- c("HBC (n=249)" = "#1f78b4",
                  "CLD (n=228)" = "#33a02c",
                  "HCC (n=116)" = "#e31a1c")

#=========================#
# Create Plot
#=========================#
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = Group_Label, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c("HBC (n=249)", "CLD (n=228)"),
    c("CLD (n=228)", "HCC (n=116)"),
    c("HBC (n=249)", "HCC (n=116)")
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +  # Removed x-axis label
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),  # Explicitly remove x-axis title
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  )

#=========================#
# Print the Plot
#=========================#
print(p)
save_fig("FigureS1_F_Flocoumafen", p, width = 7, height = 5)












##############5#################################
########################################################
## FINAL: Bromadiolone (pM) Plot with Sample Size Labels
########################################################

# Load required libraries
library(ggplot2)
library(ggpubr)

#=========================#
# Define Parameters
#=========================#
data_file     <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column  <- "tissuetype"
value_column  <- "bromodiolon_calcresponse"
y_axis_label  <- "Bromadiolone (pM)"

#=========================#
# Load and Prepare Data
#=========================#
filtered_data <- get(data_file)
filtered_data[[group_column]] <- gsub("PC", "HBC", filtered_data[[group_column]])

# Filter non-zero and valid groups
filtered_data <- filtered_data[
  !is.na(filtered_data[[group_column]]) &
    !is.na(filtered_data[[value_column]]) &
    filtered_data[[value_column]] > 0 &
    filtered_data[[group_column]] %in% c("HBC", "CLD", "HCC"), ]

# Convert from nM to pM
filtered_data[[value_column]] <- filtered_data[[value_column]] * 1000

#=========================#
# Set custom group labels with sample size
#=========================#
group_levels <- c("HBC", "CLD", "HCC")
group_labels <- c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")
filtered_data$Group_Label <- factor(filtered_data[[group_column]],
                                    levels = group_levels,
                                    labels = group_labels)

#=========================#
# Calculate Medians
#=========================#
medians <- aggregate(filtered_data[[value_column]],
                     by = list(Group_Label = filtered_data$Group_Label),
                     FUN = median, na.rm = TRUE)
colnames(medians) <- c("Group_Label", "Median")
medians$Label <- paste0("Med: ", format(round(medians$Median, 2), nsmall = 2))
medians$Position <- medians$Median

#=========================#
# Define group colors
#=========================#
group_colors <- c("HBC (n=249)" = "#1f78b4",
                  "CLD (n=228)" = "#33a02c",
                  "HCC (n=116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(filtered_data, aes(x = Group_Label, y = !!sym(value_column),
                               fill = Group_Label, color = Group_Label)) +
  geom_boxplot(outlier.shape = NA, width = 0.8, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = Group_Label, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  stat_compare_means(comparisons = list(
    c("HBC (n=249)", "CLD (n=228)"),
    c("CLD (n=228)", "HCC (n=116)"),
    c("HBC (n=249)", "HCC (n=116)")
  ),
  method = "wilcox.test",
  label = "p.signif",
  size = 6,
  fontface = "bold",
  color = "black",
  step.increase = 0.08,
  bracket.size = 1.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(x = NULL, y = y_axis_label, title = "") +  # Removed x-axis label
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  ) +
  ylim(min(filtered_data[[value_column]], na.rm = TRUE) * 0.9,
       max(filtered_data[[value_column]], na.rm = TRUE) * 1.3)

#=========================#
# Print the Plot
#=========================#
print(p)
save_fig("FigureS1_G_Bromadiolone", p, width = 7, height = 5)







########################################################
##### Butachlor (nM) — FINAL Plot with N in Labels #####
########################################################

# Load required libraries
library(ggplot2)
library(ggsignif)
library(dplyr)
library(rstatix)

#=========================#
# Define Parameters
#=========================#
data_file       <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column    <- "tissuetype"
value_column    <- "butachlor_calcresponse"
y_axis_label    <- "Butachlor (nM)"
y_axis_limits   <- c(0, 0.004)

#=========================#
# Load and Prepare Data
#=========================#
df <- get(data_file)
df[[group_column]] <- gsub("PC", "HBC", df[[group_column]])

df <- df %>%
  filter(!is.na(.data[[group_column]]),
         !is.na(.data[[value_column]]),
         .data[[value_column]] > 0,
         .data[[group_column]] %in% c("HBC", "CLD", "HCC")) %>%
  mutate(
    group = factor(.data[[group_column]],
                   levels = c("HBC", "CLD", "HCC"),
                   labels = c("HBC (n = 249)", "CLD (n = 228)", "HCC (n = 116)"))
  )

#=========================#
# Compute Medians
#=========================#
medians <- df %>%
  group_by(group) %>%
  summarise(Median = median(.data[[value_column]], na.rm = TRUE)) %>%
  mutate(Label = paste0("Med: ", formatC(Median, format = "f", digits = 4)),
         Position = Median)

#=========================#
# Wilcoxon + Significance Labels
#=========================#
stat_test <- wilcox_test(df, as.formula(paste(value_column, "~ group"))) %>%
  add_significance("p") %>%
  arrange(group1, group2) %>%
  mutate(
    p_label = p.signif,
    y.position = c(0.0035, 0.0033, 0.0037)
  )

colnames(stat_test)[colnames(stat_test) == "group1"] <- "group"
colnames(stat_test)[colnames(stat_test) == "group2"] <- "group_to"

#=========================#
# Group Colors
#=========================#
group_colors <- c("HBC (n = 249)" = "#1f78b4",
                  "CLD (n = 228)" = "#33a02c",
                  "HCC (n = 116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(df, aes(x = group, y = .data[[value_column]],
                    fill = group, color = group)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = group, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  geom_signif(data = stat_test,
              aes(xmin = group,
                  xmax = group_to,
                  annotations = p_label,
                  y_position = y.position),
              manual = TRUE,
              textsize = 6,
              size = 1.2,
              tip_length = 0.0012,
              color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = y_axis_limits) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 22, face = "bold", color = "black"),
    axis.text       = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  )

#=========================#
# Print the Plot
#=========================#
print(p)






########################################################
##### Butachlor (pM) Plot — Final with N in Labels #####
########################################################

# Load required libraries
library(ggplot2)
library(ggsignif)
library(dplyr)
library(rstatix)

#=========================#
# Define Parameters
#=========================#
data_file       <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column    <- "tissuetype"
value_column    <- "butachlor_calcresponse"
y_axis_label    <- "Butachlor (pM)"
y_axis_limits   <- c(0, 4)

#=========================#
# Load and Prepare Data
#=========================#
df <- get(data_file)

df[[group_column]] <- gsub("PC", "HBC", df[[group_column]])

df <- df %>%
  filter(!is.na(.data[[group_column]]),
         !is.na(.data[[value_column]]),
         .data[[value_column]] > 0,
         .data[[group_column]] %in% c("HBC", "CLD", "HCC")) %>%
  mutate(
    group = factor(.data[[group_column]],
                   levels = c("HBC", "CLD", "HCC"),
                   labels = c("HBC (n = 249)", "CLD (n = 228)", "HCC (n = 116)")),
    value_pM = .data[[value_column]] * 1000
  )

#=========================#
# Compute Medians
#=========================#
medians <- df %>%
  group_by(group) %>%
  summarise(Median = median(value_pM, na.rm = TRUE)) %>%
  mutate(
    Label = paste0("Med: ", formatC(Median, format = "f", digits = 4)),
    Position = Median
  )

#=========================#
# Compute Wilcoxon + Significance Labels
#=========================#
stat_test <- wilcox_test(df, value_pM ~ group) %>%
  add_significance("p") %>%
  arrange(group1, group2) %>%
  mutate(
    p_label = p.signif,
    y.position = c(3.5, 3.3, 3.7)
  )

colnames(stat_test)[colnames(stat_test) == "group1"] <- "group"
colnames(stat_test)[colnames(stat_test) == "group2"] <- "group_to"

#=========================#
# Define Group Colors
#=========================#
group_colors <- c("HBC (n = 249)" = "#1f78b4",
                  "CLD (n = 228)" = "#33a02c",
                  "HCC (n = 116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(df, aes(x = group, y = value_pM, fill = group, color = group)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = group, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  geom_signif(data = stat_test,
              aes(xmin = group,
                  xmax = group_to,
                  annotations = p_label,
                  y_position = y.position),
              manual = TRUE,
              textsize = 6,
              size = 1.2,
              tip_length = 0.0012,
              color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = y_axis_limits) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 22, face = "bold", color = "black"),
    axis.text       = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  )

#=========================#
# Print the Plot
#=========================#
print(p)
save_fig("FigureS1_B_Butachlor_pM", p, width = 7, height = 5)











#####7#################################################Oringianl Standrad for all others##############################################################
##7######################################################
##### Oxadiazon (nM) Plot — Signif vs NS Labels ########
########################################################
########################################################
##### Oxadiazon (nM) Plot — With Sample Sizes #########
########################################################

# Load required libraries
library(ggplot2)
library(ggsignif)
library(dplyr)
library(rstatix)

#=========================#
# Define Parameters
#=========================#
data_file       <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column    <- "tissuetype"
value_column    <- "oxadiazon_calcresponse"
y_axis_label    <- "Oxadiazon (nM)"
y_axis_limits   <- c(0, 0.01)

#=========================#
# Load and Prepare Data
#=========================#
df <- get(data_file)
df[[group_column]] <- gsub("PC", "HBC", df[[group_column]])

df <- df %>%
  filter(!is.na(.data[[group_column]]),
         !is.na(.data[[value_column]]),
         .data[[value_column]] > 0,
         .data[[group_column]] %in% c("HBC", "CLD", "HCC")) %>%
  mutate(Group = factor(.data[[group_column]],
                        levels = c("HBC", "CLD", "HCC"),
                        labels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)")))

#=========================#
# Compute Medians
#=========================#
medians <- df %>%
  group_by(Group) %>%
  summarise(Median = median(.data[[value_column]], na.rm = TRUE)) %>%
  mutate(Label = paste0("Med: ", formatC(Median, format = "f", digits = 4)),
         Position = Median)

#=========================#
# Compute Wilcoxon + Significance Labels
#=========================#
stat_test <- wilcox_test(df, as.formula(paste(value_column, "~ Group"))) %>%
  add_significance("p") %>%
  arrange(group1, group2) %>%
  mutate(p_label = p.signif,
         y.position = c(0.0087, 0.0081, 0.0093))

colnames(stat_test)[colnames(stat_test) == "group1"] <- "Group"
colnames(stat_test)[colnames(stat_test) == "group2"] <- "Group_to"

#=========================#
# Define Group Colors
#=========================#
group_colors <- c("HBC (n=249)" = "#1f78b4",
                  "CLD (n=228)" = "#33a02c",
                  "HCC (n=116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(df, aes(x = Group, y = .data[[value_column]],
                    fill = Group, color = Group)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = Group, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  geom_signif(data = stat_test,
              aes(xmin = Group, xmax = Group_to,
                  annotations = p_label, y_position = y.position),
              manual = TRUE,
              textsize = 6, size = 1.2, tip_length = 0.012, color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = y_axis_limits) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    axis.title.y    = element_text(size = 22, face = "bold", color = "black"),
    axis.text       = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  )

#=========================#
# Print the Plot
#=========================#
print(p)



########################################################
##### Oxadiazon (pM) Plot — Final Presentation #########
########################################################

# Load required libraries
library(ggplot2)
library(ggsignif)
library(dplyr)
library(rstatix)

#=========================#
# Define Parameters
#=========================#
data_file       <- "Pesticides_Northern_hospital_samples_Madja_LHC_1"
group_column    <- "tissuetype"
value_column    <- "oxadiazon_calcresponse"
y_axis_label    <- "Oxadiazon (pM)"
y_axis_limits   <- c(0, 10)

#=========================#
# Load and Prepare Data
#=========================#
df <- get(data_file)

df[[group_column]] <- gsub("PC", "HBC", df[[group_column]])

df <- df %>%
  filter(!is.na(.data[[group_column]]),
         !is.na(.data[[value_column]]),
         .data[[value_column]] > 0,
         .data[[group_column]] %in% c("HBC", "CLD", "HCC")) %>%
  mutate(
    group = factor(.data[[group_column]],
                   levels = c("HBC", "CLD", "HCC"),
                   labels = c("HBC (n = 249)", "CLD (n = 228)", "HCC (n = 116)")),
    value_pM = .data[[value_column]] * 1000
  )

#=========================#
# Compute Medians
#=========================#
medians <- df %>%
  group_by(group) %>%
  summarise(Median = median(value_pM, na.rm = TRUE)) %>%
  mutate(Label = paste0("Med: ", formatC(Median, format = "f", digits = 4)),
         Position = Median)

#=========================#
# Compute Wilcoxon Tests
#=========================#
stat_test <- wilcox_test(df, value_pM ~ group) %>%
  add_significance("p") %>%
  arrange(group1, group2) %>%
  mutate(
    p_label = p.signif,
    y.position = c(8.7, 8.1, 9.3)
  )

colnames(stat_test)[colnames(stat_test) == "group1"] <- "group"
colnames(stat_test)[colnames(stat_test) == "group2"] <- "group_to"

#=========================#
# Define Group Colors
#=========================#
group_colors <- c("HBC (n = 249)" = "#1f78b4",
                  "CLD (n = 228)" = "#33a02c",
                  "HCC (n = 116)" = "#e31a1c")

#=========================#
# Create the Plot
#=========================#
p <- ggplot(df, aes(x = group, y = value_pM, fill = group, color = group)) +
  geom_boxplot(width = 0.8, outlier.shape = NA, alpha = 0.8, size = 1.2) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  geom_text(data = medians,
            aes(x = group, y = Position, label = Label),
            size = 6, fontface = "bold", color = "black", vjust = -0.5) +
  geom_signif(data = stat_test,
              aes(xmin = group, xmax = group_to, annotations = p_label, y_position = y.position),
              manual = TRUE,
              textsize = 6, size = 1.2, tip_length = 0.012, color = "black") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(ylim = y_axis_limits) +
  labs(x = NULL, y = y_axis_label, title = "") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22, face = "bold", color = "black"),
    axis.text    = element_text(size = 20, face = "bold", color = "black"),
    legend.position = "none",
    panel.grid.major = element_line(size = 0.8, color = "grey80"),
    panel.grid.minor = element_blank(),
    axis.line.y.left = element_line(size = 2, color = "black"),
    axis.line.x.bottom = element_line(size = 2, color = "black")
  )

#=========================#
# Print the Plot
#=========================#
print(p)
save_fig("FigureS1_C_Oxadiazon_pM", p, width = 7, height = 5)


























###################Radadr PLot#####################################



#########Revised 2 Radar Plot (Final, Fixed Font)##############################
#########FINAL USED FOR THE MANUSCRIPT##############################

#=========================#
# Load Required Libraries
#=========================#
library(dplyr)
library(fmsb)
library(tibble)

#=========================#
# Load and Prepare Data
#=========================#
data_file <- Pesticides_Northern_hospital_samples_Madja_LHC_1
data_file$tissuetype <- gsub("PC", "HBC", data_file$tissuetype)

# Define pesticide columns and rename mapping
pesticide_cols <- c(
  "Pendimethalin_Concnormcreatinine",
  "metsulfuronmethyl_Conc.normcreatinine",
  "dichlorophenoxyaceticacid_conc_normcreatininecolourimetric",
  "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine",
  "bromodiolon_calcresponse",
  "butachlor_calcresponse",
  "oxadiazon_calcresponse",
  "Cypermethrin (nM)",
  "Glyphosate (nM)_M",
  "Phosphoric acid(nM)_M",
  "AMPA (nM)_M"
)

pesticide_rename <- c(
  "Pendimethalin_Concnormcreatinine" = "Pendimethalin (nM)",
  "metsulfuronmethyl_Conc.normcreatinine" = "Metsulfuron-methyl (nM)",
  "dichlorophenoxyaceticacid_conc_normcreatininecolourimetric" = "2,4-Dichlorophenoxyacetic acid (nM)",
  "flocoumafenpeakscombined_calculatedconc_forconcaboveLOQ_normcreatinine" = "Flocoumafen (nM)",
  "bromodiolon_calcresponse" = "Bromadiolone (pM)",
  "butachlor_calcresponse" = "Butachlor (pM)",
  "oxadiazon_calcresponse" = "Oxadiazon (pM)",
  "Cypermethrin (nM)" = "Cypermethrin (nM)",
  "Glyphosate (nM)_M" = "Glyphosate (nM)",
  "Phosphoric acid(nM)_M" = "PPA (nM)",
  "AMPA (nM)_M" = "AMPA (nM)"
)

#=========================#
# Compute Median per Group
#=========================#
df_medians <- data_file %>%
  filter(tissuetype %in% c("HBC", "CLD", "HCC")) %>%
  select(tissuetype, all_of(pesticide_cols)) %>%
  group_by(tissuetype) %>%
  summarise(across(all_of(pesticide_cols), ~ median(.x, na.rm = TRUE))) %>%
  column_to_rownames("tissuetype") %>%
  t() %>% as.data.frame()

rownames(df_medians) <- pesticide_rename[rownames(df_medians)]

#=========================#
# Normalize Each Pesticide (Row-wise)
#=========================#
df_scaled <- as.data.frame(t(apply(df_medians, 1, function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) rep(0.5, length(x)) else (x - rng[1]) / (rng[2] - rng[1])
})))

# Transpose back so rows = groups
df_scaled <- as.data.frame(t(df_scaled))

# Add fmsb-required max/min rows
radar_data <- rbind(rep(1, ncol(df_scaled)), rep(0, ncol(df_scaled)), df_scaled)

#=========================#
# Plot Radar Chart
#=========================#

# Enlarge and bold axis tick label text (0, 0.25, ..., 1)
par(cex.axis = 3.5, font.axis = 2)  # ✅ Safe for all systems, no font warning

# Global font size (impacts axis labels too)
par(cex = 1.2, text.font = 2)  # ✅ Enlarges caxislabels and title/labels


# Define colors
colors_border <- c("blue", "green", "red")
colors_fill <- adjustcolor(colors_border, alpha.f = 0.25)

# Generate radar plot
radarchart(
  radar_data,
  axistype = 1,
  pcol = colors_border,
  pfcol = colors_fill,
  plwd = 1.8,
  plty = 1,
  cglcol = "black",
  cglty = 2,
  axislabcol = "black",
  caxislabels = seq(0, 1, 0.25),
  cglwd = 1.5,
  vlcex = 1.6,
  title = "Radar Plot of Median Pesticide Concentrations"
)

#=========================#
# Enlarged Bold Legend
#=========================#
legend("topright",
       legend = c("CLD", "HBC", "HCC"),
       bty = "n",
       pch = 20,
       col = colors_border,
       text.col = "black",
       text.font = 2,
       pt.cex = 5,
       cex = 2.5
)

#=========================#
# Exposure Level Scale Panel
#=========================#
legend("bottomleft",
       title = "Exposure Level (normalized)",
       legend = c("High (1.0)", "↑", "Medium (0.5)", "↑", "Low (0.0)"),
       text.col = "black",
       text.font = 2,
       cex = 1.1,
       bty = "n"
)


###############################################################################
###############################################################################
# SUPPLEMENTARY FIGURE S2: ADJUSTED ODDS RATIOS (FOREST PLOTS), HIGH VS LOW EXPOSURE
###############################################################################
###############################################################################

#####################################NO HBV and  HCV#######################################################
#==============================#
# TIGER-LC: PC vs HCC
# Model: PLSACAT (0 = Low, 1 = High)
#==============================#

# Load Required Libraries
library(tidyverse)
library(broom)

# Load & Prepare Data
df <- Copy_of_TIGERL_2025_Madja_IJC_Dax01 %>%
  filter(tissuetypeORG %in% c("PC", "HCC")) %>%
  mutate(
    CaseStatus     = ifelse(tissuetypeORG == "HCC", 1, 0),
    PLSACAT        = as.numeric(PLSACAT),
    Sex            = ifelse(gender_MS == "Male", 1, 0),
    Occupation_new = ifelse(`Occupation new` == 1, 1, 0),  # 1 = Agri
    Alcoholcat     = as.numeric(Alcoholcat),
    Age            = as.numeric(Agecat_DUPLICATE)
  ) %>%
  drop_na(CaseStatus, PLSACAT, Sex, Occupation_new, Alcoholcat, Age)

# Fit Logistic Regression Model
model <- glm(CaseStatus ~ PLSACAT + Sex + Occupation_new + Alcoholcat + Age,
             data = df, family = binomial())

# Create Odds Ratio Table
or_table <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      term == "(Intercept)"      ~ "Intercept",
      term == "PLSACAT"          ~ "PLS11 (Low / High)",
      term == "Sex"              ~ "Sex (Male / Female)",
      term == "Occupation_new"   ~ "Occupation (Agri / Non)",
      term == "Alcoholcat"       ~ "Alcohol (Ever / Never)",
      term == "Age"              ~ "Age (<55 / ≥55)",
      TRUE ~ term
    ),
    label = case_when(
      p.value < 0.0001 ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' < 0.0001')",
                                 estimate, conf.low, conf.high),
      TRUE ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' = %.4f')",
                     estimate, conf.low, conf.high, p.value)
    )
  ) %>%
  mutate(term = factor(term, levels = rev(unique(term))))

# Plot Forest Plot (Log Scale)
p_s2 <- ggplot(or_table, aes(x = estimate, y = term)) +
  geom_point(size = 3.5, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(x = 15, label = label),
            parse = TRUE, hjust = 0, size = 6.5, color = "black") +
  scale_x_log10(
    limits = c(0.01, 20),
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1", "10")
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = expression("Adjusted Odds Ratios (log scale) for PLS"[11]*" (Low / High): HCC vs PC"),
    x = "Odds Ratio (log scale)",
    y = "Variables"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    plot.margin  = margin(10, 300, 10, 10)
  )
print(p_s2)
save_fig("FigureS2_D_PLS11_HCC", p_s2, width = 16, height = 7)













#==============================#
# TIGER-LC: PC vs HCC
# Model: PLSBCAT (0 = Low, 1 = High)
#==============================#

# Load Required Libraries
library(tidyverse)
library(broom)

# Load & Prepare Data
df <- Copy_of_TIGERL_2025_Madja_IJC_Dax01 %>%
  filter(tissuetypeORG %in% c("PC", "HCC")) %>%
  mutate(
    CaseStatus     = ifelse(tissuetypeORG == "HCC", 1, 0),
    PLSBCAT        = as.numeric(PLSBCAT),
    Sex            = ifelse(gender_MS == "Male", 1, 0),
    Occupation_new = ifelse(`Occupation new` == 1, 1, 0),  # 1 = Agri
    Alcoholcat     = as.numeric(Alcoholcat),
    Age            = as.numeric(Agecat_DUPLICATE)
  ) %>%
  drop_na(CaseStatus, PLSBCAT, Sex, Occupation_new, Alcoholcat, Age)

# Fit Logistic Regression Model
model <- glm(CaseStatus ~ PLSBCAT + Sex + Occupation_new + Alcoholcat + Age,
             data = df, family = binomial())

# Create Odds Ratio Table
or_table <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      term == "(Intercept)"      ~ "Intercept",
      term == "PLSBCAT"          ~ "PLS8 (Low / High)",
      term == "Sex"              ~ "Sex (Male / Female)",
      term == "Occupation_new"   ~ "Occupation (Agri / Non)",
      term == "Alcoholcat"       ~ "Alcohol (Ever / Never)",
      term == "Age"              ~ "Age (<55 / ≥55)",
      TRUE ~ term
    ),
    label = case_when(
      p.value < 0.0001 ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' < 0.0001')",
                                 estimate, conf.low, conf.high),
      TRUE ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' = %.4f')",
                     estimate, conf.low, conf.high, p.value)
    )
  ) %>%
  mutate(term = factor(term, levels = rev(unique(term))))

# Plot Forest Plot (Log Scale)
p_s2 <- ggplot(or_table, aes(x = estimate, y = term)) +
  geom_point(size = 3.5, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(x = 15, label = label),
            parse = TRUE, hjust = 0, size = 6.5, color = "black") +
  scale_x_log10(
    limits = c(0.01, 20),
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1", "10")
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = expression("Adjusted Odds Ratios (log scale) for PLS"[8]*" (Low / High): HCC vs PC"),
    x = "Odds Ratio (log scale)",
    y = "Variables"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    plot.margin  = margin(10, 300, 10, 10)
  )
print(p_s2)
save_fig("FigureS2_E_PLS8_HCC", p_s2, width = 16, height = 7)








#==============================#
# TIGER-LC: PC vs HCC
# Model: CypermethrinCAT (0 = Low, 1 = High)
#==============================#

# Load Required Libraries
library(tidyverse)
library(broom)

# Load & Prepare Data
df <- Copy_of_TIGERL_2025_Madja_IJC_Dax01 %>%
  filter(tissuetypeORG %in% c("PC", "HCC")) %>%
  mutate(
    CaseStatus       = ifelse(tissuetypeORG == "HCC", 1, 0),
    CypermethrinCAT  = as.numeric(CypermethrinCAT),
    Sex              = ifelse(gender_MS == "Male", 1, 0),
    Occupation_new   = ifelse(`Occupation new` == 1, 1, 0),
    Alcoholcat       = as.numeric(Alcoholcat),
    Age              = as.numeric(Agecat_DUPLICATE)
  ) %>%
  drop_na(CaseStatus, CypermethrinCAT, Sex, Occupation_new, Alcoholcat, Age)

# Fit Logistic Regression Model
model <- glm(CaseStatus ~ CypermethrinCAT + Sex + Occupation_new + Alcoholcat + Age,
             data = df, family = binomial())

# Create Odds Ratio Table
or_table <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      term == "(Intercept)"       ~ "Intercept",
      term == "CypermethrinCAT"   ~ "Cypermethrin (Low / High)",
      term == "Sex"               ~ "Sex (Male / Female)",
      term == "Occupation_new"    ~ "Occupation (Agri / Non)",
      term == "Alcoholcat"        ~ "Alcohol (Ever / Never)",
      term == "Age"               ~ "Age (<55 / ≥55)",
      TRUE ~ term
    ),
    label = case_when(
      p.value < 0.0001 ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' < 0.0001')",
                                 estimate, conf.low, conf.high),
      TRUE ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' = %.4f')",
                     estimate, conf.low, conf.high, p.value)
    )
  ) %>%
  mutate(term = factor(term, levels = rev(unique(term))))

# Plot Forest Plot (Log Scale)
p_s2 <- ggplot(or_table, aes(x = estimate, y = term)) +
  geom_point(size = 3.5, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(x = 15, label = label),
            parse = TRUE, hjust = 0, size = 6.5, color = "black") +
  scale_x_log10(
    limits = c(0.01, 20),
    breaks = c(0.01, 0.1, 1, 10, 10),
    labels = c("0.01", "0.1", "1", "10", "10")
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Adjusted Odds Ratios (log scale) for Cypermethrin (Low / High): HCC vs PC",
    x = "Odds Ratio (log scale)",
    y = "Variables"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    plot.margin  = margin(10, 300, 10, 10)
  )
print(p_s2)
save_fig("FigureS2_F_Cypermethrin_HCC", p_s2, width = 16, height = 7)



















#==============================#
# TIGER-LC: PC vs CLD
# Model: CypermethrinCAT (0 = Low, 1 = High)
#==============================#

# Load Required Libraries
library(tidyverse)
library(broom)

# Load & Prepare Data
df <- Copy_of_TIGERL_2025_Madja_IJC_Dax01 %>%
  filter(tissuetypeORG %in% c("PC", "CLD")) %>%
  mutate(
    CaseStatus        = ifelse(tissuetypeORG == "CLD", 1, 0),
    CypermethrinCAT   = as.numeric(CypermethrinCAT),
    Sex               = ifelse(gender_MS == "Male", 0, 1),
    Occupationcat     = ifelse(Occupationcat == 1, 1, 0),         # 1 = Agri
    Alcoholcat        = as.numeric(Alcoholcat),
    Age               = as.numeric(Agecat_DUPLICATE)
  ) %>%
  drop_na(CaseStatus, CypermethrinCAT, Sex, Occupationcat, Alcoholcat, Age)

# Fit Logistic Regression Model
model <- glm(CaseStatus ~ CypermethrinCAT + Sex + Occupationcat + Alcoholcat + Age,
             data = df, family = binomial())

# Create OR Table
or_table <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      term == "(Intercept)"        ~ "Intercept",
      term == "CypermethrinCAT"    ~ "Cypermethrin (Low / High)",
      term == "Sex"                ~ "Sex (Male / Female)",
      term == "Occupationcat"      ~ "Occupation (Agri / Non)",
      term == "Alcoholcat"         ~ "Alcohol (Ever / Never)",
      term == "Age"                ~ "Age (<55 / ≥55)",
      TRUE ~ term
    ),
    label = case_when(
      p.value < 0.0001 ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' < 0.0001')",
                                 estimate, conf.low, conf.high),
      TRUE ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' = %.4f')",
                     estimate, conf.low, conf.high, p.value)
    )
  ) %>%
  mutate(term = factor(term, levels = rev(unique(term))))

# Plot Forest Plot (Log Scale)
p_s2 <- ggplot(or_table, aes(x = estimate, y = term)) +
  geom_point(size = 3.5, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(x = 15, label = label),
            parse = TRUE, hjust = 0, size = 6.5, color = "black") +
  scale_x_log10(
    limits = c(0.01, 20),
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1", "10")
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Adjusted Odds Ratios (log scale) for Cypermethrin (Low / High): CLD vs PC",
    x = "Odds Ratio (log scale)",
    y = "Variables"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    plot.margin  = margin(10, 300, 10, 10)
  )
print(p_s2)
save_fig("FigureS2_C_Cypermethrin_CLD", p_s2, width = 16, height = 7)














##################################### WITHOUT HBV/HCV FILTER #######################################################
#==============================#
# TIGER-LC: PC vs CLD
# Model: PLSACAT (0 = Low, 1 = High)
#==============================#

# Load Required Libraries
library(tidyverse)
library(broom)

# Load & Prepare Data
df <- Copy_of_TIGERL_2025_Madja_IJC_Dax01 %>%
  filter(tissuetypeORG %in% c("PC", "CLD")) %>%
  mutate(
    CaseStatus     = ifelse(tissuetypeORG == "CLD", 1, 0),
    PLSACAT        = as.numeric(PLSACAT),
    Sex            = ifelse(gender_MS == "Male", 0, 1),
    Occupationcat  = as.numeric(Occupationcat),  # New: Occupation (1 = Agri)
    Alcoholcat     = as.numeric(Alcoholcat),
    Age            = as.numeric(Agecat_DUPLICATE)
  ) %>%
  drop_na(CaseStatus, PLSACAT, Sex, Occupationcat, Alcoholcat, Age)

# Fit Logistic Regression Model
model <- glm(CaseStatus ~ PLSACAT + Sex + Occupationcat + Alcoholcat + Age,
             data = df, family = binomial())

# Create Odds Ratio Table
or_table <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      term == "(Intercept)"      ~ "Intercept",
      term == "PLSACAT"          ~ "PLS11 (Low / High)",
      term == "Sex"              ~ "Sex (Male / Female)",
      term == "Occupationcat"    ~ "Occupation (Agri / Non)",
      term == "Alcoholcat"       ~ "Alcohol (Ever / Never)",
      term == "Age"              ~ "Age (<55 / ≥55)",
      TRUE ~ term
    ),
    label = case_when(
      p.value < 0.0001 ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' < 0.0001')",
                                 estimate, conf.low, conf.high),
      TRUE ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' = %.4f')",
                     estimate, conf.low, conf.high, p.value)
    )
  ) %>%
  mutate(term = factor(term, levels = rev(unique(term))))

# Plot Forest Plot (Log Scale)
p_s2 <- ggplot(or_table, aes(x = estimate, y = term)) +
  geom_point(size = 3.5, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(x = 15, label = label),
            parse = TRUE, hjust = 0, size = 6.5, color = "black") +
  scale_x_log10(
    limits = c(0.01, 20),
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1", "10")
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = expression("Adjusted Odds Ratios (log scale) for PLS"[11]*" (Low / High): CLD vs PC"),
    x = "Odds Ratio (log scale)",
    y = "Variables"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    plot.margin  = margin(10, 300, 10, 10)
  )
print(p_s2)
save_fig("FigureS2_A_PLS11_CLD", p_s2, width = 16, height = 7)





















##################################### WITHOUT HBV/HCV FILTER #######################################################
#==============================#
# TIGER-LC: PC vs CLD
# Model: PLSBCAT (0 = Low, 1 = High)
#==============================#

# Load Required Libraries
library(tidyverse)
library(broom)

# Load & Prepare Data
df <- Copy_of_TIGERL_2025_Madja_IJC_Dax01 %>%
  filter(tissuetypeORG %in% c("PC", "CLD")) %>%
  mutate(
    CaseStatus     = ifelse(tissuetypeORG == "CLD", 1, 0),
    PLSBCAT        = as.numeric(PLSBCAT),
    Sex            = ifelse(gender_MS == "Male", 0, 1),
    Occupationcat  = as.numeric(Occupationcat),  # 1 = Agri
    Alcoholcat     = as.numeric(Alcoholcat),
    Age            = as.numeric(Agecat_DUPLICATE)
  ) %>%
  drop_na(CaseStatus, PLSBCAT, Sex, Occupationcat, Alcoholcat, Age)

# Fit Logistic Regression Model
model <- glm(CaseStatus ~ PLSBCAT + Sex + Occupationcat + Alcoholcat + Age,
             data = df, family = binomial())

# Create Odds Ratio Table
or_table <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      term == "(Intercept)"      ~ "Intercept",
      term == "PLSBCAT"          ~ "PLS8 (Low / High)",
      term == "Sex"              ~ "Sex (Male / Female)",
      term == "Occupationcat"    ~ "Occupation (Agri / Non)",
      term == "Alcoholcat"       ~ "Alcohol (Ever / Never)",
      term == "Age"              ~ "Age (<55 / ≥55)",
      TRUE ~ term
    ),
    label = case_when(
      p.value < 0.0001 ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' < 0.0001')",
                                 estimate, conf.low, conf.high),
      TRUE ~ sprintf("bold('OR = %.2f (%.2f–%.2f); ') * bold(italic('p')) * bold(' = %.4f')",
                     estimate, conf.low, conf.high, p.value)
    )
  ) %>%
  mutate(term = factor(term, levels = rev(unique(term))))

# Plot Forest Plot (Log Scale)
p_s2 <- ggplot(or_table, aes(x = estimate, y = term)) +
  geom_point(size = 3.5, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.2, size = 1.2, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_text(aes(x = 15, label = label),
            parse = TRUE, hjust = 0, size = 6.5, color = "black") +
  scale_x_log10(
    limits = c(0.01, 20),
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1", "10")
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = expression("Adjusted Odds Ratios (log scale) for PLS"[8]*" (Low / High): CLD vs PC"),
    x = "Odds Ratio (log scale)",
    y = "Variables"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title   = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    axis.text.x  = element_text(size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, face = "bold", color = "black"),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 18, face = "bold", color = "black"),
    plot.margin  = margin(10, 300, 10, 10)
  )
print(p_s2)
save_fig("FigureS2_B_PLS8_CLD", p_s2, width = 16, height = 7)


###############################################################################
###############################################################################
# SUPPLEMENTARY FIGURE S5: SHAP SUMMARY PLOTS
###############################################################################
###############################################################################

###########################################################################################
############################################################################################
############################FOr the manuscript##############################################
###########################################################################################



# ==========================#######################
# SHAP Summary Plot: PLS₁₁ + Covariates (CLD and HCC Separately)
# ==========================#######################

# 1. Load Required Libraries
library(xgboost)
library(SHAPforxgboost)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

# 2. Load and Prepare Data
df <- PILCRS_1_MASTER
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

# 3. Define Covariates and Label Names
covariates <- c("PLS_GroupA", "Gender_Cat", "Alcoholcat", "Ageorg", "Occupation new")
label_names <- c("PLS11", "Sex", "Alcohol", "Age", "Occupation")  # Plain text for plot labels

# 4. SHAP Prep Function
prepare_shap_plot <- function(shap_long, covars, labels) {
  shap_long$variable <- factor(shap_long$variable, levels = covars, labels = labels)
  mean_shap <- shap_long %>%
    group_by(variable) %>%
    summarize(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
    arrange(desc(mean_abs_shap))
  shap_long$variable <- factor(shap_long$variable, levels = mean_shap$variable)
  return(shap_long)
}

# 5. SHAP: CLD vs HBC
df_cld <- df %>%
  filter(tissuetypeORG %in% c("HBC (n=249)", "CLD (n=228)")) %>%
  filter(complete.cases(across(all_of(covariates)))) %>%
  mutate(RiskLabel = ifelse(tissuetypeORG == "HBC (n=249)", 0, 1))

X_cld <- as.matrix(df_cld[, covariates])
y_cld <- df_cld$RiskLabel

model_cld <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1),
                       data = xgb.DMatrix(X_cld, label = y_cld),
                       nrounds = 300, verbose = 0)

shap_values_cld <- shap.values(xgb_model = model_cld, X_train = X_cld)
shap_long_cld <- shap.prep(shap_contrib = shap_values_cld$shap_score, X_train = X_cld)
shap_long_cld <- prepare_shap_plot(shap_long_cld, covariates, label_names)

p_cld <- shap.plot.summary(shap_long_cld) +
  ggtitle(expression("SHAP Summary: CLD vs HBC (PLS"[11]*")")) +
  ylab("SHAP value (impact on model output)") +
  
  theme(
    plot.title    = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 11, face = "bold", color = "black")
  )

# 6. SHAP: HCC vs HBC
df_hcc <- df %>%
  filter(tissuetypeORG %in% c("HBC (n=249)", "HCC (n=116)")) %>%
  filter(complete.cases(across(all_of(covariates)))) %>%
  mutate(RiskLabel = ifelse(tissuetypeORG == "HBC (n=249)", 0, 1))

X_hcc <- as.matrix(df_hcc[, covariates])
y_hcc <- df_hcc$RiskLabel

model_hcc <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1),
                       data = xgb.DMatrix(X_hcc, label = y_hcc),
                       nrounds = 300, verbose = 0)

shap_values_hcc <- shap.values(xgb_model = model_hcc, X_train = X_hcc)
shap_long_hcc <- shap.prep(shap_contrib = shap_values_hcc$shap_score, X_train = X_hcc)
shap_long_hcc <- prepare_shap_plot(shap_long_hcc, covariates, label_names)

p_hcc <- shap.plot.summary(shap_long_hcc) +
  ggtitle(expression("SHAP Summary: HCC vs HBC (PLS"[11]*")")) +
  ylab("SHAP value (impact on model output)") +
  
  theme(
    plot.title    = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 11, face = "bold", color = "black")
  )

# 7. Combine Both SHAP Plots (display in R only)
print(p_cld + p_hcc)
save_fig("FigureS5_SHAP_PLS11", p_cld + p_hcc, width = 14, height = 6)





























# ==========================#######################
# SHAP Summary Plot: PLS₈ + Covariates (CLD and HCC Separately)
# ==========================#######################

# 1. Load Required Libraries
library(xgboost)
library(SHAPforxgboost)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

# 2. Load and Prepare Data
df <- PILCRS_1_MASTER
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

# 3. Define Covariates and Label Names
covariates <- c("PLS_GroupB", "Gender_Cat", "Alcoholcat", "Ageorg", "Occupation new")
label_names <- c("PLS8", "Sex", "Alcohol", "Age", "Occupation")  # Plain text for plot labels

# 4. SHAP Prep Function
prepare_shap_plot <- function(shap_long, covars, labels) {
  shap_long$variable <- factor(shap_long$variable, levels = covars, labels = labels)
  mean_shap <- shap_long %>%
    group_by(variable) %>%
    summarize(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
    arrange(desc(mean_abs_shap))
  shap_long$variable <- factor(shap_long$variable, levels = mean_shap$variable)
  return(shap_long)
}

# 5. SHAP: CLD vs HBC
df_cld <- df %>%
  filter(tissuetypeORG %in% c("HBC (n=249)", "CLD (n=228)")) %>%
  filter(complete.cases(across(all_of(covariates)))) %>%
  mutate(RiskLabel = ifelse(tissuetypeORG == "HBC (n=249)", 0, 1))

X_cld <- as.matrix(df_cld[, covariates])
y_cld <- df_cld$RiskLabel

model_cld <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1),
                       data = xgb.DMatrix(X_cld, label = y_cld),
                       nrounds = 300, verbose = 0)

shap_values_cld <- shap.values(xgb_model = model_cld, X_train = X_cld)
shap_long_cld <- shap.prep(shap_contrib = shap_values_cld$shap_score, X_train = X_cld)
shap_long_cld <- prepare_shap_plot(shap_long_cld, covariates, label_names)

p_cld <- shap.plot.summary(shap_long_cld) +
  ggtitle(expression("SHAP Summary: CLD vs HBC (PLS"[8]*")")) +
  ylab("SHAP value (impact on model output)") +
  theme(
    plot.title    = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 11, face = "bold", color = "black")
  )

# 6. SHAP: HCC vs HBC
df_hcc <- df %>%
  filter(tissuetypeORG %in% c("HBC (n=249)", "HCC (n=116)")) %>%
  filter(complete.cases(across(all_of(covariates)))) %>%
  mutate(RiskLabel = ifelse(tissuetypeORG == "HBC (n=249)", 0, 1))

X_hcc <- as.matrix(df_hcc[, covariates])
y_hcc <- df_hcc$RiskLabel

model_hcc <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1),
                       data = xgb.DMatrix(X_hcc, label = y_hcc),
                       nrounds = 300, verbose = 0)

shap_values_hcc <- shap.values(xgb_model = model_hcc, X_train = X_hcc)
shap_long_hcc <- shap.prep(shap_contrib = shap_values_hcc$shap_score, X_train = X_hcc)
shap_long_hcc <- prepare_shap_plot(shap_long_hcc, covariates, label_names)

p_hcc <- shap.plot.summary(shap_long_hcc) +
  ggtitle(expression("SHAP Summary: HCC vs HBC (PLS"[8]*")")) +
  ylab("SHAP value (impact on model output)") +
  theme(
    plot.title    = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 11, face = "bold", color = "black")
  )







# 7. Combine Both SHAP Plots (display in R only)
print(p_cld + p_hcc)
save_fig("FigureS5_SHAP_PLS8", p_cld + p_hcc, width = 14, height = 6)

































############################For the manuscript##############################

# ==========================#######################
# SHAP Summary Plot: Cypermethrin + Covariates (CLD and HCC Separately)
# ==========================#######################

# 1. Load Required Libraries
library(xgboost)
library(SHAPforxgboost)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)

# 2. Load and Prepare Data
df <- PILCRS_1_MASTER
df$tissuetypeORG <- recode(df$tissuetypeORG,
                           "PC" = "HBC (n=249)",
                           "CLD" = "CLD (n=228)",
                           "HCC" = "HCC (n=116)")
df$tissuetypeORG <- factor(df$tissuetypeORG,
                           levels = c("HBC (n=249)", "CLD (n=228)", "HCC (n=116)"))

# 3. Define Covariates and Label Names
covariates <- c("Cypermethrin_nM", "Gender_Cat", "Alcoholcat", "Ageorg", "Occupation new")
label_names <- c("Cypermethrin", "Sex", "Alcohol", "Age", "Occupation")

# 4. SHAP Prep Function
prepare_shap_plot <- function(shap_long, covars, labels) {
  shap_long$variable <- factor(shap_long$variable, levels = covars, labels = labels)
  mean_shap <- shap_long %>%
    group_by(variable) %>%
    summarize(mean_abs_shap = mean(abs(value)), .groups = "drop") %>%
    arrange(desc(mean_abs_shap))
  shap_long$variable <- factor(shap_long$variable, levels = mean_shap$variable)
  return(shap_long)
}

# 5. SHAP: CLD vs HBC
df_cld <- df %>%
  filter(tissuetypeORG %in% c("HBC (n=249)", "CLD (n=228)")) %>%
  filter(complete.cases(across(all_of(covariates)))) %>%
  mutate(RiskLabel = ifelse(tissuetypeORG == "HBC (n=249)", 0, 1))

X_cld <- as.matrix(df_cld[, covariates])
y_cld <- df_cld$RiskLabel

model_cld <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1),
                       data = xgb.DMatrix(X_cld, label = y_cld),
                       nrounds = 300, verbose = 0)

shap_values_cld <- shap.values(xgb_model = model_cld, X_train = X_cld)
shap_long_cld <- shap.prep(shap_contrib = shap_values_cld$shap_score, X_train = X_cld)
shap_long_cld <- prepare_shap_plot(shap_long_cld, covariates, label_names)

p_cld <- shap.plot.summary(shap_long_cld) +
  ggtitle(expression("SHAP Summary: CLD vs HBC (Cypermethrin)")) +
  ylab("SHAP value (impact on model output)") +
  theme(
    plot.title    = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 11, face = "bold", color = "black")
  )

# 6. SHAP: HCC vs HBC
df_hcc <- df %>%
  filter(tissuetypeORG %in% c("HBC (n=249)", "HCC (n=116)")) %>%
  filter(complete.cases(across(all_of(covariates)))) %>%
  mutate(RiskLabel = ifelse(tissuetypeORG == "HBC (n=249)", 0, 1))

X_hcc <- as.matrix(df_hcc[, covariates])
y_hcc <- df_hcc$RiskLabel

model_hcc <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, eta = 0.1),
                       data = xgb.DMatrix(X_hcc, label = y_hcc),
                       nrounds = 300, verbose = 0)

shap_values_hcc <- shap.values(xgb_model = model_hcc, X_train = X_hcc)
shap_long_hcc <- shap.prep(shap_contrib = shap_values_hcc$shap_score, X_train = X_hcc)
shap_long_hcc <- prepare_shap_plot(shap_long_hcc, covariates, label_names)

p_hcc <- shap.plot.summary(shap_long_hcc) +
  ggtitle(expression("SHAP Summary: HCC vs HBC (Cypermethrin)")) +
  ylab("SHAP value (impact on model output)") +
  theme(
    plot.title    = element_text(size = 14, face = "bold", color = "black", hjust = 0.5),
    axis.title.x  = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 11, face = "bold", color = "black")
  )

# 7. Combine Both SHAP Plots (display in R only)
print(p_cld + p_hcc)
save_fig("FigureS5_SHAP_Cypermethrin", p_cld + p_hcc, width = 14, height = 6)


###############################################################################
###############################################################################
# SUPPLEMENTARY FIGURE S6: ROC CURVES, HCC VS CLD
###############################################################################
###############################################################################

################################################################################
# Figure S6: ROC curves, HCC vs CLD (Reviewer #1, Comment 5)
# Uses the data already loaded in R (PILCRS_1_MASTER1)
################################################################################

library(pROC)
library(ggplot2)

# ---- 1. Data (already in environment) ----------------------------------------
df <- as.data.frame(PILCRS_1_MASTER1)
df <- df[df$tissuetypeORG %in% c("CLD", "HCC") & !is.na(df$PLS_GroupA), ]

# ---- 2. Variables -------------------------------------------------------------
df$HCC      <- ifelse(df$tissuetypeORG == "HCC", 1, 0)   # HCC = 1, CLD = 0
df$Sex      <- as.numeric(df$Gender_Cat)                  # 1 = male
df$Alcohol  <- as.numeric(df$Alcoholcat)                  # 1 = ever
df$Age      <- as.numeric(df$Ageorg)
df$Occ      <- as.numeric(df$`Occupation new`)            # 1 = agriculture
df$logPLS11 <- log2(as.numeric(df$PLS_GroupA) + 1)

afp <- as.numeric(df$afp)
iu  <- !is.na(df$afpUnit) & df$afpUnit == "IU/mL"
afp[iu] <- afp[iu] * 1.21                                 # IU/mL -> ng/mL
df$logAFP <- log10(afp)

vars <- c("HCC", "logPLS11", "logAFP", "Age", "Sex", "Alcohol", "Occ")
cc <- df[complete.cases(df[, vars]), ]
cat("Sample:", nrow(cc), "(HCC =", sum(cc$HCC), "; CLD =", sum(cc$HCC == 0), ")\n")

# ---- 3. Models and ROC --------------------------------------------------------
m_pls  <- glm(HCC ~ logPLS11 + Age + Sex + Alcohol + Occ,          family = binomial, data = cc)
m_afp  <- glm(HCC ~ logAFP   + Age + Sex + Alcohol + Occ,          family = binomial, data = cc)
m_both <- glm(HCC ~ logPLS11 + logAFP + Age + Sex + Alcohol + Occ, family = binomial, data = cc)

r_pls  <- roc(cc$HCC, fitted(m_pls),  quiet = TRUE, ci = TRUE)
r_afp  <- roc(cc$HCC, fitted(m_afp),  quiet = TRUE, ci = TRUE)
r_both <- roc(cc$HCC, fitted(m_both), quiet = TRUE, ci = TRUE)

print(data.frame(
  Model = c("PLS11 + covariates", "AFP + covariates", "PLS11 + AFP + covariates"),
  AUC   = sprintf("%.3f (%.3f-%.3f)",
                  c(r_pls$ci[2], r_afp$ci[2], r_both$ci[2]),
                  c(r_pls$ci[1], r_afp$ci[1], r_both$ci[1]),
                  c(r_pls$ci[3], r_afp$ci[3], r_both$ci[3]))))
cat("DeLong, PLS11+AFP vs PLS11: P =", signif(roc.test(r_both, r_pls)$p.value, 3), "\n")
cat("DeLong, PLS11+AFP vs AFP:   P =", signif(roc.test(r_both, r_afp)$p.value, 3), "\n")

# ---- 4. Figure ----------------------------------------------------------------
lab <- function(r, name) sprintf("%s (AUC %.3f)", name, auc(r))
roc_df <- rbind(
  data.frame(FPR = 1 - r_pls$specificities,  TPR = r_pls$sensitivities,  Model = lab(r_pls,  "PLS11 + covariates")),
  data.frame(FPR = 1 - r_afp$specificities,  TPR = r_afp$sensitivities,  Model = lab(r_afp,  "AFP + covariates")),
  data.frame(FPR = 1 - r_both$specificities, TPR = r_both$sensitivities, Model = lab(r_both, "PLS11 + AFP + covariates")))
roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))
roc_df <- roc_df[order(roc_df$Model, roc_df$FPR, roc_df$TPR), ]

p <- ggplot(roc_df, aes(FPR, TPR, color = Model)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_step(linewidth = 1.3, direction = "vh") +
  scale_color_manual(values = c("#33a02c", "#1f78b4", "#e31a1c")) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = sprintf("HCC vs CLD (n = %d; %d HCC, %d CLD)", nrow(cc), sum(cc$HCC), sum(cc$HCC == 0))) +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.62, 0.15), legend.title = element_blank(),
        axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))
print(p)
save_fig("FigureS6_ROC_HCC_vs_CLD", p, width = 7, height = 7)

# ---- 5. Save (optional) -------------------------------------------------------
# ggsave("FigureS6_HCC_vs_CLD_ROC.tiff", p, width = 7, height = 7, dpi = 600, compression = "lzw")