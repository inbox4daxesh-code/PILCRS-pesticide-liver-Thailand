###############################################################################
# PILCRS_Tables_11222026_42.R
# Builds all main and supplementary TABLES, in manuscript order.
# Manuscript: "Urinary Pesticide Profiles and Liver Disease Risk in Thailand:
#   A Machine Learning Risk Prediction Model" (Cancer Research Communications)
#
# HOW TO RUN
#   1. Click "Source" (Ctrl+Shift+S). For each dataset that is not already loaded,
#      a window opens; pick the matching Excel file (first sheet is read).
#   2. Each table is printed in the Console, drawn in the Plots pane, and saved to
#      the "output" folder as PDF and CSV. At the end, every table is also kept in
#      the R environment as a data frame (Table_1, Table_2, Table_3, Table_S3 ...)
#      and opened in the RStudio data viewer.
#   3. A summary at the end lists which steps ran and which (if any) failed.
###############################################################################

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
  "Copy_of_TIGERL_2025_Madja_IJC_Dax01",   # Copy_of_TIGERL_2025_Madja_IJC_Dax01.xlsx
  "PILCRS_1_MASTER1"                       # PILCRS_1_MASTER1.xlsx
))

## ---------------- Packages (installs any that are missing) ----------------
pkgs <- c("dplyr", "tidyr", "tibble", "gridExtra", "logistf", "pROC", "mice")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs)
}

## ---------------- Output folder ----------------
##OUT_DIR <- "output"
###dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------- Output folder ----------------
OUT_DIR <- "H:/Sep2023_Desktop Backup/TGLC_Chris_01112021/Pesticide Load Project/2025_Madja_Dax_Combine/Manuscript/Submissions/Cancer Research Communication/Major Revision_1/Data Analysis Final/Github/R/Tables"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)




## ---------------- Show tables in RStudio AND save to PDF; step runner ----------------
## ---------------- Show every plot in RStudio AND save it to PDF ----------------
# Plots are drawn in the RStudio Plots pane (use the arrows there to page back).
# Each finished page is recorded, and at the end of the step all pages are
# written to one PDF in OUT_DIR.
.cap <- new.env()
.cap_hooks <- list()
.cap_record <- function(...) {
  if (isTRUE(.cap$dirty)) { .cap$plots[[length(.cap$plots) + 1]] <- recordPlot(); .cap$dirty <- FALSE }
}
.cap_before_base <- function(...) if (isTRUE(par("page"))) .cap_record()   # only when a NEW page starts
.cap_mark <- function(...) .cap$dirty <- TRUE
start_capture <- function() {
  .cap$plots <- list(); .cap$dirty <- FALSE
  if (is.null(dev.list())) dev.new()
  dev.control("enable")
  setHook("before.plot.new", .cap_before_base); setHook("before.grid.newpage", .cap_record)
  setHook("plot.new", .cap_mark);               setHook("grid.newpage", .cap_mark)
}
stop_capture <- function() {
  for (h in c("before.plot.new", "before.grid.newpage", "plot.new", "grid.newpage")) {
    keep <- Filter(function(f) !(identical(f, .cap_before_base) || identical(f, .cap_record) || identical(f, .cap_mark)), getHook(h))
    setHook(h, NULL, "replace"); for (f in keep) setHook(h, f)
  }
}
finish_capture <- function(pdf_file, width, height) {
  if (!is.null(dev.list())) .cap_record()                       # last page
  stop_capture()
  if (!length(.cap$plots)) { cat("No plots were drawn for", pdf_file, "
"); return(invisible(NULL)) }
  screen <- dev.cur()
  cairo_pdf(file.path(OUT_DIR, pdf_file), width = width, height = height, onefile = TRUE)
  for (p in .cap$plots) suppressWarnings(replayPlot(p))
  dev.off(); if (screen %in% dev.list()) dev.set(screen)
  cat("Saved:", file.path(OUT_DIR, pdf_file), "(", length(.cap$plots), "page(s) )
")
}

## ---------------- Runner: one failure does not stop the rest ----------------
run_log <- data.frame(Step = character(0), Status = character(0), Message = character(0), stringsAsFactors = FALSE)
run_step <- function(name, expr, pdf_file = NULL, width = 7, height = 6) {
  cat("

#####################################################################
#", name,
      "
#####################################################################
")
  if (!is.null(pdf_file)) start_capture()
  msg <- tryCatch({ eval(expr, envir = globalenv()); "" }, error = function(e) conditionMessage(e))
  if (!is.null(pdf_file)) {
    if (nzchar(msg)) stop_capture() else
      msg <- tryCatch({ finish_capture(pdf_file, width, height); "" }, error = function(e) conditionMessage(e))
  }
  # close any PDF file device left open by a failed step (keeps the RStudio Plots pane)
  for (d in dev.list()) if (names(dev.list())[dev.list() == d] %in% c("pdf", "cairo_pdf")) dev.off(d)
  run_log[nrow(run_log) + 1, ] <<- list(name, if (nzchar(msg)) "FAILED" else "OK", msg)
  if (nzchar(msg)) cat("
*** FAILED:", name, "
***", msg, "
")
  invisible(NULL)
}

#############################################################################
# MAIN TABLES
#############################################################################


############################## 00 Setup analytic dataset (n = 593) ##############################
run_step("00 Setup analytic dataset", quote({
  
  ###############################################################################
  # 00_setup_analytic_data.R
  # Builds the analytic dataset (n = 593) used by all table scripts.
  # Requires these objects already loaded in the RStudio environment:
  #   Copy_of_TIGERL_2025_Madja_IJC_Dax01   (exposures, covariates, HBV/HCV)
  #   PILCRS_1_MASTER1                      (region, smoking, liver fluke)
  # Run this script first; the table scripts use the object `dat`.
  # Install once if needed: install.packages(c("dplyr", "gridExtra", "logistf", "pROC"))
  ###############################################################################
  
  suppressPackageStartupMessages({ library(dplyr); library(gridExtra); library(grid) })
  
  
  # Draws tables in the RStudio Plots pane and saves them as a multipage PDF (24 rows per page)
  save_tables_pdf <- function(tables, file, rows_per_page = 24) {
    start_capture(); on.exit(stop_capture(), add = TRUE)
    for (nm in names(tables)) {
      df <- as.data.frame(tables[[nm]], stringsAsFactors = FALSE)
      rn <- rownames(tables[[nm]])
      if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))) df <- cbind(Row = rn, df)
      pages <- split(seq_len(nrow(df)), ceiling(seq_len(nrow(df)) / rows_per_page))
      for (i in seq_along(pages)) {
        grid.newpage()
        grid.text(if (length(pages) > 1) paste0(nm, " (page ", i, " of ", length(pages), ")") else nm,
                  y = unit(0.96, "npc"), gp = gpar(fontsize = 13, fontface = "bold"))
        grid.draw(tableGrob(df[pages[[i]], , drop = FALSE], rows = NULL,
                            theme = ttheme_default(base_size = 10, padding = unit(c(3, 2), "mm"))))
      }
    }
    finish_capture(file, width = 11, height = 8.5)
    # Also save each table as a CSV (opens in Excel) next to the PDF
    stem <- sub("\\.pdf$", "", file)
    for (k in seq_along(tables)) {
      df <- as.data.frame(tables[[k]], stringsAsFactors = FALSE)
      rn <- rownames(tables[[k]])
      if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))) df <- cbind(Row = rn, df)
      csv <- if (length(tables) > 1) paste0(stem, "_part", k, ".csv") else paste0(stem, ".csv")
      write.csv(df, file.path(OUT_DIR, csv), row.names = FALSE)
      cat("Saved:", file.path(OUT_DIR, csv), "\n")
    }
    invisible(file.path(OUT_DIR, file))
  }
  
  # Column lookup that tolerates duplicate-name suffixes added on import
  pick <- function(df, name) {
    hit <- names(df)[names(df) == name]
    if (length(hit) == 0) hit <- grep(paste0("^", gsub("([.()])", "\\\\\\1", name), "(\\.\\.\\.[0-9]+)?$"), names(df), value = TRUE)
    if (length(hit) == 0) stop("Column not found: ", name)
    df[[hit[1]]]
  }
  
  src <- Copy_of_TIGERL_2025_Madja_IJC_Dax01
  dat <- data.frame(
    id    = pick(src, "TIGER-LC ID"),
    grp   = pick(src, "tissuetypeORG"),
    plsA  = pick(src, "PLS_GroupA"),
    pls11 = as.numeric(pick(src, "PLSFC_GroupAM")),   # PLS11 fold change (Table 2 quartiles)
    pls8  = as.numeric(pick(src, "PLSFC_GroupBM")),   # PLS8 fold change
    cyp   = as.numeric(pick(src, "Cypermethrin_nM")), # creatinine-normalized cypermethrin, NA = not measured
    age   = as.numeric(pick(src, "Ageorg")),
    sex   = ifelse(is.na(pick(src, "gender_MS")), NA, as.numeric(pick(src, "gender_MS") == "Male")),
    alc   = as.numeric(pick(src, "Alcoholcat")),      # 1 = ever, 0 = never
    occ   = as.numeric(pick(src, "Occupationcat")),   # 1 = agriculture
    hbv   = as.numeric(pick(src, "HBVCAT")),          # 1 = HBsAg positive
    hcv   = as.numeric(pick(src, "HCVCAT")),          # 1 = anti-HCV positive
    stringsAsFactors = FALSE
  ) %>%
    filter(grp %in% c("PC", "CLD", "HCC"), !is.na(plsA)) %>%
    mutate(grp = factor(grp, levels = c("PC", "CLD", "HCC")))
  
  m1 <- PILCRS_1_MASTER1
  extra <- data.frame(
    id     = pick(m1, "TIGER-LC ID"),
    region = pick(m1, "region_Cat"),   # 3 levels
    smoke  = pick(m1, "smoke_Cat"),    # never / past / current
    ov     = pick(m1, "ov_Cat"),       # self-reported liver fluke, 1 = yes
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(id)) %>% distinct(id, .keep_all = TRUE)
  
  dat <- left_join(dat, extra, by = "id")
  
  # Quartile cutpoints from Table 2 (upper bounds of Q1 to Q3)
  CUTS <- list(pls11 = c(0.16, 0.77, 2.22),
               pls8  = c(6.90e-8, 0.46, 1.78),
               cyp   = c(12.39, 27.78, 63.78))
  
  # missing_to_Q1 = TRUE reproduces the earlier coding (missing cypermethrin placed in Q1)
  make_q <- function(x, cuts, missing_to_Q1 = FALSE) {
    if (missing_to_Q1) x[is.na(x)] <- 0
    cut(x, c(-Inf, cuts, Inf), labels = c("Q1", "Q2", "Q3", "Q4"), right = TRUE)
  }
  
  cat("Analytic dataset:", nrow(dat), "participants\n")
  print(table(dat$grp))
  cat("Participants with measured cypermethrin by group:\n")
  print(table(dat$grp, measured = !is.na(dat$cyp)))
  
  
}))

############################## Table 1 ##############################
run_step("Table 1", quote({
  
  ###############################################################################
  # T01_Table1_counts.R
  # Recomputes Table 1 categorical counts, n (%), on the 593-participant cohort.
  # Run 00_setup_analytic_data.R first.
  ###############################################################################
  
  tab_n_pct <- function(var, label, levels = NULL, names_lv = NULL) {
    x <- dat[[var]]
    out <- sapply(c("PC", "CLD", "HCC"), function(g) {
      v <- x[dat$grp == g]; n <- length(v)
      lv <- if (is.null(levels)) sort(unique(v[!is.na(v)])) else levels
      c(sapply(lv, function(l) sprintf("%d (%.1f)", sum(v == l, na.rm = TRUE), 100 * sum(v == l, na.rm = TRUE) / n)),
        Missing = sprintf("%d (%.1f)", sum(is.na(v)), 100 * sum(is.na(v)) / n))
    })
    rn <- if (is.null(names_lv)) c(as.character(if (is.null(levels)) sort(unique(x[!is.na(x)])) else levels), "Missing") else c(names_lv, "Missing")
    rownames(out) <- paste(label, rn, sep = ": ")
    colnames(out) <- c("HBC", "CLD", "HCC")
    out
  }
  
  cat("Group sizes:\n"); print(table(dat$grp))
  cat("\nAge, mean (SD):\n")
  print(tapply(dat$age, dat$grp, function(a) sprintf("%.1f (%.1f); missing %d", mean(a, na.rm = TRUE), sd(a, na.rm = TRUE), sum(is.na(a)))))
  
  tbl1 <- rbind(
    tab_n_pct("sex",    "Sex",        c(1, 0), c("Male", "Female")),
    tab_n_pct("alc",    "Alcohol",    c(1, 0), c("Ever drinker", "Never")),
    tab_n_pct("occ",    "Occupation", c(1, 0), c("Agriculture", "Other")),
    tab_n_pct("hbv",    "HBsAg",      c(1, 0), c("Positive", "Negative")),
    tab_n_pct("hcv",    "Anti-HCV",   c(1, 0), c("Positive", "Negative")),
    tab_n_pct("region", "Region (code)"),
    tab_n_pct("smoke",  "Smoking (code)"),
    tab_n_pct("ov",     "Liver fluke", c(1, 0), c("Yes", "No"))
  )
  print(noquote(tbl1))
  
  cat("\nP values (chi-square; ANOVA for age):\n")
  for (v in c("sex", "alc", "occ", "hbv", "hcv", "region", "smoke", "ov"))
    cat(sprintf("  %-7s P = %.3g\n", v, suppressWarnings(chisq.test(table(dat$grp, dat[[v]]))$p.value)))
  cat(sprintf("  %-7s P = %.3g\n", "age", anova(lm(age ~ grp, data = dat))[1, "Pr(>F)"]))
  
  cat("\nMissing data, % of 593 (Methods text):\n")
  print(round(100 * colMeans(is.na(dat[, c("age", "sex", "alc", "occ", "hbv", "hcv", "region", "smoke", "ov")])), 1))
  
  pv <- sapply(c("sex", "alc", "occ", "hbv", "hcv", "region", "smoke", "ov"),
               function(v) signif(suppressWarnings(chisq.test(table(dat$grp, dat[[v]]))$p.value), 3))
  age_tab <- data.frame(Group = levels(dat$grp),
                        Age_mean_SD = tapply(dat$age, dat$grp, function(a) sprintf("%.1f (%.1f)", mean(a, na.rm = TRUE), sd(a, na.rm = TRUE))),
                        Missing = tapply(dat$age, dat$grp, function(a) sum(is.na(a))))
  save_tables_pdf(list(`Table 1: n (%) by group` = tbl1,
                       `Table 1: age` = age_tab,
                       `Table 1: P values (chi-square; ANOVA for age)` =
                         data.frame(Variable = c(names(pv), "age"),
                                    P = c(pv, signif(anova(lm(age ~ grp, data = dat))[1, "Pr(>F)"], 3)))),
                  "Table1_11222026_42.pdf")
  
  
}))

############################## Table 2 ##############################
run_step("Table 2", quote({
  
  ###############################################################################
  # T02_Table2_TableS6_quartile_ORs.R
  # Table 2 (unadjusted) and Supplementary Table S6 (adjusted for age, sex,
  # alcohol use, and agricultural occupation): quartile odds ratios for PLS11,
  # PLS8, and cypermethrin, each case group versus HBC.
  # Cypermethrin models use participants with a measured concentration.
  # Run 00_setup_analytic_data.R first.
  ###############################################################################
  
  quartile_table <- function(expo, case, covs = character(0), missing_to_Q1 = FALSE) {
    d <- dat[dat$grp %in% c("PC", case), ]
    d$y <- as.integer(d$grp == case)
    d$q <- make_q(d[[expo]], CUTS[[expo]], missing_to_Q1)
    d$qn <- as.numeric(d$q)
    d <- d[complete.cases(d[, c("y", "q", covs)]), ]
    rhs <- paste(c("q", covs), collapse = " + ")
    fit  <- glm(as.formula(paste("y ~", rhs)), data = d, family = binomial)
    fitT <- glm(as.formula(paste("y ~", paste(c("qn", covs), collapse = " + "))), data = d, family = binomial)
    co <- summary(fit)$coefficients; ci <- suppressMessages(confint.default(fit))
    hbc <- table(d$q[d$y == 0])
    res <- data.frame(
      Exposure = expo, Outcome = case, Quartile = levels(d$q),
      HBC_n_pct = sprintf("%d (%.1f)", as.integer(hbc), 100 * as.integer(hbc) / sum(hbc)),
      OR_95CI = c("1.00 (Reference)", sapply(c("qQ2", "qQ3", "qQ4"), function(k)
        sprintf("%.2f (%.2f-%.2f)", exp(co[k, 1]), exp(ci[k, 1]), exp(ci[k, 2])))),
      P = c("-", ifelse(co[c("qQ2", "qQ3", "qQ4"), 4] < 0.0001, "<0.0001", sprintf("%.4f", co[c("qQ2", "qQ3", "qQ4"), 4]))),
      stringsAsFactors = FALSE)
    attr(res, "p_trend") <- summary(fitT)$coefficients["qn", 4]
    attr(res, "n") <- c(N = nrow(d), cases = sum(d$y))
    res
  }
  
  run_set <- function(covs, title) {
    cat("\n==============================", title, "==============================\n")
    out <- list()
    for (expo in c("pls11", "pls8", "cyp")) for (case in c("CLD", "HCC")) {
      r <- quartile_table(expo, case, covs); n <- attr(r, "n")
      cat(sprintf("\n%s, %s vs HBC (N = %d; cases = %d)\n", expo, case, n["N"], n["cases"]))
      print(r[, -(1:2)], row.names = FALSE)
      cat(sprintf("P trend = %.4g\n", attr(r, "p_trend")))
      r$N_cases <- c(sprintf("%d (%d)", n["N"], n["cases"]), "", "", "")
      r$P_trend <- c(sprintf("%.4g", attr(r, "p_trend")), "", "", "")
      out[[paste(expo, case)]] <- r
    }
    res <- do.call(rbind, out); rownames(res) <- NULL; res
  }
  
  base_covs <- c("age", "sex", "alc", "occ")
  t2  <- run_set(character(0), "TABLE 2: UNADJUSTED")
  save_tables_pdf(list(`Table 2: quartile odds ratios, unadjusted` = t2), "Table2_11222026_42.pdf")
  
}))


############################## Table 3 ##############################
run_step("Table 3", quote({
  
  ###############################################################################
  # T07_remaining_tables.R
  # Table 3 and Supplementary Tables S3, S4, S5, S8, S9, S12, S13, S14, S15,
  # S16, S17. Uses objects already loaded in RStudio:
  #   Copy_of_TIGERL_2025_Madja_IJC_Dax01, PILCRS_1_MASTER1
  # Run 00_setup_analytic_data.R first (creates `dat` and `pick()`).
  ###############################################################################
  
  suppressPackageStartupMessages({ library(dplyr); library(pROC); library(splines) })
  set.seed(2025)
  
  src <- Copy_of_TIGERL_2025_Madja_IJC_Dax01
  extra2 <- data.frame(
    id    = pick(src, "TIGER-LC ID"),
    cypFC = suppressWarnings(as.numeric(pick(src, "CypermethrinFC_nM"))),
    gly   = suppressWarnings(as.numeric(pick(src, "Glyphosate (nM)"))),
    ampa  = suppressWarnings(as.numeric(pick(src, "AMPA (nM)"))),
    ppa   = suppressWarnings(as.numeric(pick(src, "Phosphoric acid(nM)"))),
    site  = pick(src, "Site...15"),
    fib4  = suppressWarnings(as.numeric(pick(src, "FIB-4"))),
    apri  = suppressWarnings(as.numeric(pick(src, "APRI"))),
    ast   = suppressWarnings(as.numeric(pick(src, "AST (SGOT) U/L"))),
    alt   = suppressWarnings(as.numeric(pick(src, "alt"))),
    afp   = suppressWarnings(as.numeric(pick(src, "afp"))),
    stringsAsFactors = FALSE
  ) %>% filter(!is.na(id)) %>% distinct(id, .keep_all = TRUE)
  dat <- left_join(dat[, setdiff(names(dat), setdiff(names(extra2), "id"))], extra2, by = "id")
  
  EXPO <- c(PLS11 = "pls11", PLS8 = "pls8", Cypermethrin = "cyp")
  BASE <- c("age", "sex", "alc", "occ")
  
  prep <- function(case, vars) {
    d <- dat[dat$grp %in% c("PC", case), ]
    d$y <- as.integer(d$grp == case)
    d[complete.cases(d[, c("y", vars)]), ]
  }
  cv_auc <- function(d, f, k = 10) {
    folds <- sample(rep(1:k, length.out = nrow(d))); p <- numeric(nrow(d)); fa <- c()
    for (i in 1:k) {
      m <- glm(f, data = d[folds != i, ], family = binomial)
      p[folds == i] <- predict(m, d[folds == i, ], type = "response")
      if (length(unique(d$y[folds == i])) == 2) fa <- c(fa, as.numeric(pROC::auc(d$y[folds == i], p[folds == i], quiet = TRUE)))
    }
    list(auc = as.numeric(pROC::auc(d$y, p, quiet = TRUE)), mean = mean(fa), sd = sd(fa), brier = mean((p - d$y)^2), pred = p)
  }
  fml <- function(vars) as.formula(paste("y ~", paste(vars, collapse = " + ")))
  
  ## ---------------- Table 3: cross-validated performance ----------------
  t3 <- do.call(rbind, lapply(names(EXPO), function(e) do.call(rbind, lapply(c("CLD","HCC"), function(o) {
    v <- c(EXPO[[e]], BASE); d <- prep(o, v); r <- cv_auc(d, fml(v))
    data.frame(Model = paste(e, "+ covariates"), Outcome = o, N = nrow(d), Cases = sum(d$y),
               AUC_CV = round(r$auc, 3), Mean_SD = sprintf("%.3f +/- %.3f", r$mean, r$sd), Brier = round(r$brier, 3))
  }))))
  cat("\n=== Table 3 ===\n"); print(t3, row.names = FALSE)
  save_tables_pdf(list(`Table 3: 10-fold cross-validated performance` = t3), "Table3_11222026_42.pdf")
  
}))


#############################################################################
# SUPPLEMENTARY TABLES
# Table S1 (LC-MS/MS assay parameters) is descriptive and has no code.
#############################################################################


############################## Supplementary Table S2 ##############################
run_step("Supplementary Table S2", quote({
  
  ###############################################################################
  # Supplementary Table S2: multiple imputation by chained equations (MICE) for
  # creatinine-normalized cypermethrin (44.5% missing).
  # 20 imputed datasets, 10 iterations; imputation predictors: age, sex, disease
  # group, occupation, HBV, HCV, and PLS11. Pooled with Rubin's rules.
  # Odds ratio: highest vs lowest cypermethrin quartile (Table 2 cutpoints),
  # case group vs HBC, unadjusted as in Table 2. Three analyses are shown:
  #   - Table 2 coding: participants without cypermethrin placed in Q1
  #   - Complete case:  participants with measured cypermethrin only
  #   - MICE (m = 20):  missing cypermethrin imputed
  # The random-number state is restored afterwards so later tables are unchanged.
  ###############################################################################
  
  suppressPackageStartupMessages(library(mice))
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  
  imp <- mice(dat[, c("cyp", "age", "sex", "grp", "occ", "hbv", "hcv", "pls11")],
              m = 20, maxit = 10, printFlag = FALSE, seed = 2025)
  
  q4_fit <- function(d, case, missing_to_Q1 = FALSE) {
    d <- d[d$grp %in% c("PC", case), ]
    d$y <- as.integer(d$grp == case)
    d$q <- make_q(d$cyp, CUTS$cyp, missing_to_Q1)
    glm(y ~ q, data = d, family = binomial)
  }
  fmt  <- function(est, lo, hi) sprintf("%.2f (%.2f-%.2f)", exp(est), exp(lo), exp(hi))
  fmtp <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  
  s2 <- do.call(rbind, lapply(c("CLD", "HCC"), function(o) {
    f1 <- q4_fit(dat, o, missing_to_Q1 = TRUE)
    f2 <- q4_fit(dat, o, missing_to_Q1 = FALSE)
    ci1 <- confint.default(f1)["qQ4", ]; ci2 <- confint.default(f2)["qQ4", ]
    fits <- lapply(seq_len(imp$m), function(i) q4_fit(mice::complete(imp, i), o))
    pl <- summary(pool(as.mira(fits)), conf.int = TRUE); r <- pl[pl$term == "qQ4", ]
    data.frame(
      Analysis = c("Table 2 coding (missing in Q1)", "Complete case (measured only)", "MICE (m = 20)"),
      Outcome  = o,
      OR_Q4_vs_Q1 = c(fmt(coef(f1)["qQ4"], ci1[1], ci1[2]), fmt(coef(f2)["qQ4"], ci2[1], ci2[2]),
                      fmt(r$estimate, log(exp(r[["2.5 %"]])), log(exp(r[["97.5 %"]])))),
      P = c(fmtp(summary(f1)$coefficients["qQ4", 4]), fmtp(summary(f2)$coefficients["qQ4", 4]), fmtp(r$p.value)),
      N = c(nobs(f1), nobs(f2), nobs(fits[[1]])),
      stringsAsFactors = FALSE)
  }))
  print(s2, row.names = FALSE)
  save_tables_pdf(list(`Supplementary Table S2: MICE for cypermethrin (Q4 vs Q1)` = s2), "TableS2_11222026_42.pdf")
  
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())
  
}))


############################## Supplementary Tables S3, S4, S5 ##############################
run_step("Supplementary Tables S3-S5", quote({
  # Uses helper functions defined in the Table 3 step.
  
  ## ---------------- Table S3: Cramer's V ----------------
  cramer <- function(a, b) { t <- table(a, b); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]; cs <- suppressWarnings(chisq.test(t, correct = FALSE))
  c(V = sqrt(cs$statistic / (sum(t) * (min(dim(t)) - 1))), chi2 = cs$statistic, df = cs$parameter, p = cs$p.value) }
  cat("\n=== Table S3 (check that `site` codes map to the three hospital networks) ===\n"); print(table(dat$site))
  s3 <- round(rbind(`Network vs region` = cramer(dat$site, dat$region),
                    `Region vs group`   = cramer(dat$region, dat$grp),
                    `Network vs group`  = cramer(dat$site, dat$grp)), 3)
  colnames(s3) <- c("Cramer_V", "Chi_square", "df", "P"); print(s3)
  
  ## ---------------- Table S4: Spearman among exposure metrics ----------------
  pairs <- list(c("pls11","pls8"), c("pls11","cyp"), c("pls8","cyp"), c("pls11","gly"), c("pls11","ampa"), c("pls11","ppa"),
                c("pls8","gly"), c("pls8","ampa"), c("pls8","ppa"), c("cyp","gly"), c("cyp","ampa"), c("cyp","ppa"))
  s4 <- do.call(rbind, lapply(pairs, function(p) { ct <- suppressWarnings(cor.test(dat[[p[1]]], dat[[p[2]]], method = "spearman"))
  data.frame(Pair = paste(p, collapse = " vs "), rho = round(ct$estimate, 3), p = signif(ct$p.value, 3)) }))
  cat("\n=== Table S4 ===\n"); print(s4, row.names = FALSE)
  
  ## ---------------- Table S5: sex-stratified CV AUC ----------------
  s5 <- do.call(rbind, lapply(c(1, 0), function(s) do.call(rbind, lapply(names(EXPO), function(e) do.call(rbind, lapply(c("CLD","HCC"), function(o) {
    v <- c(EXPO[[e]], "age", "alc", "occ"); d <- prep(o, c(v, "sex")); d <- d[d$sex == s, ]; r <- cv_auc(d, fml(v))
    data.frame(Sex = ifelse(s == 1, "Male", "Female"), Exposure = e, Outcome = o, N = nrow(d), Cases = sum(d$y), AUC_CV = round(r$auc, 3))
  }))))))
  cat("\n=== Table S5 ===\n"); print(s5, row.names = FALSE)
  save_tables_pdf(list(`Supplementary Table S3: Cramer's V` = s3), "TableS3_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S4: Spearman correlations among exposure metrics` = s4), "TableS4_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S5: sex-stratified cross-validated AUC` = s5), "TableS5_11222026_42.pdf")
  
}))


############################## Supplementary Table S6 ##############################
run_step("Supplementary Table S6", quote({
  
  # Quartile odds ratios adjusted for age, sex, alcohol use, and agricultural
  # occupation. Uses quartile_table() and run_set() defined in the Table 2 step.
  ts6 <- run_set(base_covs, "SUPPLEMENTARY TABLE S6: ADJUSTED (age, sex, alcohol, occupation)")
  save_tables_pdf(list(`Supplementary Table S6: quartile odds ratios adjusted for age, sex, alcohol, occupation` = ts6),
                  "TableS6_11222026_42.pdf")
  
}))


############################## Supplementary Table S7 ##############################
run_step("Supplementary Table S7", quote({
  
  ###############################################################################
  # T03_TableS7_HBV_HCV_adjusted.R
  # Highest versus lowest quartile odds ratios additionally adjusted for HBV and
  # HCV serostatus (age, sex, alcohol use, agricultural occupation, HBsAg,
  # anti-HCV), with Firth penalized-likelihood estimates (logistf) as a
  # sensitivity analysis because HBsAg positivity is rare among controls.
  # Run 00_setup_analytic_data.R first.
  ###############################################################################
  
  # Note: glm may warn "fitted probabilities numerically 0 or 1 occurred" because HBsAg
  # positivity is rare among controls; this is expected and is why Firth estimates are shown.
  suppressPackageStartupMessages(library(logistf))
  
  covs_viral <- c("age", "sex", "alc", "occ", "hbv", "hcv")
  
  viral_model <- function(expo, case, covs = covs_viral) {
    d <- dat[dat$grp %in% c("PC", case), ]
    d$y <- as.integer(d$grp == case)
    d$q <- make_q(d[[expo]], CUTS[[expo]])
    d <- d[complete.cases(d[, c("y", "q", covs)]), ]
    f <- as.formula(paste("y ~", paste(c("q", covs), collapse = " + ")))
    ml <- glm(f, data = d, family = binomial)
    co <- summary(ml)$coefficients["qQ4", ]; ci <- suppressMessages(confint.default(ml))["qQ4", ]
    fi <- logistf(f, data = d)
    k  <- which(names(coef(fi)) == "qQ4")
    data.frame(
      Exposure = expo, Outcome = case,
      N_cases  = sprintf("%d (%d)", nrow(d), sum(d$y)),
      OR_95CI  = sprintf("%.2f (%.2f-%.2f)", exp(co[1]), exp(ci[1]), exp(ci[2])),
      P        = sprintf("%.3f", co[4]),
      Firth_OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(coef(fi)[k]), exp(fi$ci.lower[k]), exp(fi$ci.upper[k])),
      Firth_P  = sprintf("%.3f", fi$prob[k]),
      stringsAsFactors = FALSE)
  }
  
  res <- do.call(rbind, lapply(c("pls11", "pls8", "cyp"), function(e)
    do.call(rbind, lapply(c("CLD", "HCC"), function(o) viral_model(e, o)))))
  cat("\nQ4 vs Q1, adjusted for age, sex, alcohol, occupation, HBV, and HCV\n")
  print(res, row.names = FALSE)
  
  cat("\nControls with known HBV / HCV status and positives:\n")
  print(with(dat[dat$grp == "PC", ], c(HBV_known = sum(!is.na(hbv)), HBV_pos = sum(hbv == 1, na.rm = TRUE),
                                       HCV_known = sum(!is.na(hcv)), HCV_pos = sum(hcv == 1, na.rm = TRUE))))
  
  save_tables_pdf(list(`Supplementary Table S7: Q4 vs Q1 adjusted for age, sex, alcohol, occupation, HBV, HCV (with Firth)` = res), "TableS7_11222026_42.pdf")
  
  
}))

############################## Supplementary Tables S8, S9 ##############################
run_step("Supplementary Tables S8-S9", quote({
  # Uses helper functions defined in the Table 3 step.
  
  ## ---------------- Table S8: + HBV and HCV ----------------
  s8 <- do.call(rbind, lapply(names(EXPO), function(e) do.call(rbind, lapply(c("CLD","HCC"), function(o) {
    v <- c(EXPO[[e]], BASE, "hbv", "hcv"); d <- prep(o, v); r <- cv_auc(d, fml(v))
    data.frame(Model = paste(e, "+ covariates + HBV/HCV"), Outcome = o, N = nrow(d), Cases = sum(d$y), AUC_CV = round(r$auc, 3), Brier = round(r$brier, 3))
  }))))
  cat("\n=== Table S8 ===\n"); print(s8, row.names = FALSE)
  
  ## ---------------- Table S9: within-group Spearman with liver markers ----------------
  s9 <- do.call(rbind, lapply(c(PLS11 = "pls11", PLS8 = "pls8", `Cypermethrin FC` = "cypFC"), function(e)
    do.call(rbind, lapply(c("CLD","HCC"), function(g) do.call(rbind, lapply(c(`FIB-4` = "fib4", APRI = "apri", AST = "ast", ALT = "alt", AFP = "afp"), function(b) {
      d <- dat[dat$grp == g & !is.na(dat[[e]]) & !is.na(dat[[b]]), ]
      ct <- suppressWarnings(cor.test(d[[e]], d[[b]], method = "spearman"))
      data.frame(Exposure = e, Biomarker = b, Group = g, N = nrow(d), rho = round(ct$estimate, 3), p = signif(ct$p.value, 3))
    }))))))
  cat("\n=== Table S9 ===\n"); print(s9, row.names = FALSE)
  save_tables_pdf(list(`Supplementary Table S8: models with HBV and HCV` = s8), "TableS8_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S9: within-group Spearman with liver markers` = s9), "TableS9_11222026_42.pdf")
  
}))


############################## Supplementary Table S10 ##############################
run_step("Supplementary Table S10", quote({
  
  ###############################################################################
  # T04_TableS10_Evalues.R
  # E-values (VanderWeele and Ding) for the adjusted Q4 vs Q1 odds ratios of
  # Supplementary Table S6, treating the OR as an approximation of the risk ratio.
  # Run 00_setup_analytic_data.R first.
  ###############################################################################
  
  evalue <- function(rr) ifelse(rr <= 1, 1, rr + sqrt(rr * (rr - 1)))
  
  covs <- c("age", "sex", "alc", "occ")
  ev <- do.call(rbind, lapply(c("pls11", "pls8", "cyp"), function(expo) {
    do.call(rbind, lapply(c("CLD", "HCC"), function(case) {
      d <- dat[dat$grp %in% c("PC", case), ]
      d$y <- as.integer(d$grp == case)
      d$q <- make_q(d[[expo]], CUTS[[expo]])
      d <- d[complete.cases(d[, c("y", "q", covs)]), ]
      fit <- glm(as.formula(paste("y ~", paste(c("q", covs), collapse = " + "))), data = d, family = binomial)
      or <- exp(coef(fit)["qQ4"]); lo <- exp(suppressMessages(confint.default(fit))["qQ4", 1])
      data.frame(Exposure = paste(expo, "Q4 vs Q1"), Outcome = case,
                 Adjusted_OR = round(or, 2), CI_lower = round(lo, 2),
                 E_value = round(evalue(or), 2), E_value_CI = round(evalue(lo), 2))
    }))
  }))
  print(ev, row.names = FALSE)
  
  save_tables_pdf(list(`Supplementary Table S10: E-values` = ev), "TableS10_11222026_42.pdf")
  
  
}))

############################## Supplementary Table S11 ##############################
run_step("Supplementary Table S11", quote({
  
  ###############################################################################
  # T05_TableS11_sensitivity.R
  # Supplementary Table S11: Q4 vs Q1 odds ratios with further adjustment for
  # region of residence, smoking status, and self-reported liver fluke history.
  # Base model: age, sex, alcohol use, agricultural occupation.
  # Cypermethrin models use participants with a measured concentration.
  # Run 00_setup_analytic_data.R first.
  ###############################################################################
  
  base <- c("age", "sex", "alc", "occ")
  specs <- list(
    "Base"                            = list(covs = base,                          sample = NULL),
    "+ Region"                        = list(covs = c(base, "region"),             sample = NULL),
    "+ Smoking"                       = list(covs = c(base, "smoke"),              sample = NULL),
    "+ Region + smoking"              = list(covs = c(base, "region", "smoke"),    sample = NULL),
    "Base (liver fluke-known sample)" = list(covs = base,                          sample = "ov"),
    "+ Liver fluke history"           = list(covs = c(base, "ov"),                 sample = NULL)
  )
  
  fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  
  s11 <- do.call(rbind, lapply(c("pls11", "pls8", "cyp"), function(expo) {
    do.call(rbind, lapply(c("CLD", "HCC"), function(case) {
      do.call(rbind, lapply(names(specs), function(nm) {
        sp <- specs[[nm]]
        d <- dat[dat$grp %in% c("PC", case), ]
        d$y <- as.integer(d$grp == case)
        d$q <- make_q(d[[expo]], CUTS[[expo]]); d$qn <- as.numeric(d$q)
        d$region <- factor(d$region); d$smoke <- factor(d$smoke)
        d <- d[complete.cases(d[, c("y", "q", sp$covs, sp$sample)]), ]
        fit  <- glm(as.formula(paste("y ~", paste(c("q", sp$covs), collapse = " + "))), data = d, family = binomial)
        fitT <- glm(as.formula(paste("y ~", paste(c("qn", sp$covs), collapse = " + "))), data = d, family = binomial)
        co <- summary(fit)$coefficients["qQ4", ]; ci <- suppressMessages(confint.default(fit))["qQ4", ]
        data.frame(Exposure = expo, Outcome = case, Model = nm,
                   N_cases = sprintf("%d (%d)", nrow(d), sum(d$y)),
                   OR_95CI = sprintf("%.2f (%.2f-%.2f)", exp(co[1]), exp(ci[1]), exp(ci[2])),
                   P = fmt_p(co[4]), P_trend = fmt_p(summary(fitT)$coefficients["qn", 4]))
      }))
    }))
  }))
  print(s11, row.names = FALSE)
  
  save_tables_pdf(list(`Supplementary Table S11: sensitivity analyses` = s11), "TableS11_11222026_42.pdf")
  
  
}))

############################## Supplementary Tables S12-S17 ##############################
run_step("Supplementary Tables S12-S17", quote({
  # Uses helper functions defined in the Table 3 step.
  
  ## ---------------- Table S12: restricted cubic splines ----------------
  s12 <- do.call(rbind, lapply(names(EXPO), function(e) do.call(rbind, lapply(c("CLD","HCC"), function(o) {
    d <- prep(o, c(EXPO[[e]], BASE)); d$x <- d[[EXPO[[e]]]]
    m1 <- glm(fml(c("x", BASE)), data = d, family = binomial)
    m2 <- glm(fml(c("ns(x, df = 4)", BASE)), data = d, family = binomial)
    data.frame(Exposure = e, Outcome = o, AUC_linear = round(as.numeric(pROC::auc(d$y, fitted(m1), quiet = TRUE)), 3),
               AUC_spline = round(as.numeric(pROC::auc(d$y, fitted(m2), quiet = TRUE)), 3),
               LR_p = signif(anova(m1, m2, test = "LRT")[2, "Pr(>Chi)"], 3))
  }))))
  cat("\n=== Table S12 ===\n"); print(s12, row.names = FALSE)
  
  ## ---------------- Table S13: PAF (Miettinen), above control median ----------------
  s13 <- do.call(rbind, lapply(names(EXPO), function(e) do.call(rbind, lapply(c("CLD","HCC"), function(o) {
    d <- prep(o, c(EXPO[[e]], BASE)); med <- median(d[[EXPO[[e]]]][d$y == 0]); d$hi <- as.integer(d[[EXPO[[e]]]] > med)
    or <- exp(coef(glm(fml(c("hi", BASE)), data = d, family = binomial))["hi"]); pc <- mean(d$hi[d$y == 1])
    data.frame(Exposure = e, Outcome = o, Cases_exposed_pct = round(100 * pc), Adjusted_OR = round(or, 2), PAF_pct = round(100 * pc * (or - 1) / or, 1))
  }))))
  cat("\n=== Table S13 ===\n"); print(s13, row.names = FALSE)
  
  ## ---------------- Tables S14 and S15: age- and alcohol-stratified ----------------
  strat <- function(labels, covs) do.call(rbind, lapply(names(labels), function(lv) do.call(rbind, lapply(names(EXPO), function(e) do.call(rbind, lapply(c("CLD","HCC"), function(o) {
    d <- prep(o, c(EXPO[[e]], BASE)); d <- d[eval(parse(text = labels[[lv]]), d), ]
    d$z <- as.numeric(scale(d[[EXPO[[e]]]])); m <- glm(fml(c("z", covs)), data = d, family = binomial)
    data.frame(Stratum = lv, Exposure = e, Outcome = o, N = nrow(d), Cases = sum(d$y),
               AUC = round(as.numeric(pROC::auc(d$y, fitted(m), quiet = TRUE)), 3), Exposure_OR_perSD = round(exp(coef(m)["z"]), 2))
  }))))))
  s14 <- strat(list(`<55 yr` = "age < 55", `>=55 yr` = "age >= 55"), c("sex","alc","occ"))
  cat("\n=== Table S14 (age) ===\n"); print(s14, row.names = FALSE)
  s15 <- strat(list(`Non-drinker` = "alc == 0", Drinker = "alc == 1"), c("age","sex","occ"))
  cat("\n=== Table S15 (alcohol) ===\n"); print(s15, row.names = FALSE)
  
  ## ---------------- Table S16: NRI and IDI (10-fold CV predictions) ----------------
  s16 <- do.call(rbind, lapply(names(EXPO), function(e) do.call(rbind, lapply(c("CLD","HCC"), function(o) {
    d <- prep(o, c(EXPO[[e]], BASE)); folds <- sample(rep(1:10, length.out = nrow(d)))
    pb <- pn <- numeric(nrow(d))
    for (i in 1:10) {
      tr <- d[folds != i, ]; te <- d[folds == i, ]
      pb[folds == i] <- predict(glm(fml(BASE), data = tr, family = binomial), te, type = "response")
      pn[folds == i] <- predict(glm(fml(c(EXPO[[e]], BASE)), data = tr, family = binomial), te, type = "response")
    }
    up <- pn > pb; ev <- d$y == 1
    nri <- (mean(up[ev]) - mean(!up[ev])) + (mean(!up[!ev]) - mean(up[!ev]))
    idi <- (mean(pn[ev]) - mean(pn[!ev])) - (mean(pb[ev]) - mean(pb[!ev]))
    data.frame(Exposure = e, Outcome = o, AUC_base = round(as.numeric(pROC::auc(d$y, pb, quiet = TRUE)), 3),
               AUC_plus = round(as.numeric(pROC::auc(d$y, pn, quiet = TRUE)), 3), NRI = round(nri, 3), IDI = round(idi, 3))
  }))))
  cat("\n=== Table S16 ===\n"); print(s16, row.names = FALSE)
  
  ## ---------------- Table S17: incremental AUC (same participants for all three steps) ----------------
  s17 <- do.call(rbind, lapply(c("CLD","HCC"), function(o) do.call(rbind, lapply(names(EXPO), function(e) {
    steps <- list(`Demographics only` = BASE, `+ HBV/HCV` = c(BASE, "hbv", "hcv"), `+ HBV/HCV + exposure` = c(BASE, "hbv", "hcv", EXPO[[e]]))
    out <- do.call(rbind, lapply(names(steps), function(s) { d <- prep(o, steps[[3]]); r <- cv_auc(d, fml(steps[[s]]))
    data.frame(Outcome = o, Exposure = e, Model = s, N = nrow(d), CV_AUC = round(r$auc, 3)) }))
    out$dAUC <- c(NA, round(diff(out$CV_AUC), 3)); out
  }))))
  cat("\n=== Table S17 ===\n"); print(s17, row.names = FALSE)
  
  
  save_tables_pdf(list(`Supplementary Table S12: restricted cubic splines` = s12), "TableS12_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S13: population attributable fraction` = s13), "TableS13_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S14: age-stratified` = s14), "TableS14_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S15: alcohol-stratified` = s15), "TableS15_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S16: NRI and IDI` = s16), "TableS16_11222026_42.pdf")
  save_tables_pdf(list(`Supplementary Table S17: incremental AUC` = s17), "TableS17_11222026_42.pdf")
  
}))



## ---------------- Keep every table in the R environment and open it ----------------
tab_or_null <- function(x) if (exists(x, envir = globalenv())) as.data.frame(get(x, envir = globalenv())) else NULL
with_rownames <- function(df) { if (is.null(df)) return(NULL); rn <- rownames(df)
  if (!identical(rn, as.character(seq_len(nrow(df))))) df <- cbind(Row = rn, df); rownames(df) <- NULL; df }
All_tables <- Filter(Negate(is.null), list(
  Table_1       = with_rownames(tab_or_null("tbl1")),
  Table_1_age   = with_rownames(tab_or_null("age_tab")),
  Table_1_P     = if (exists("pv")) data.frame(Variable = names(pv), P = unname(pv)) else NULL,
  Table_2       = tab_or_null("t2"),
  Table_3       = tab_or_null("t3"),
  Table_S2      = tab_or_null("s2"),
  Table_S3      = with_rownames(tab_or_null("s3")),
  Table_S4      = tab_or_null("s4"),
  Table_S5      = tab_or_null("s5"),
  Table_S6      = tab_or_null("ts6"),
  Table_S7      = tab_or_null("res"),
  Table_S8      = tab_or_null("s8"),
  Table_S9      = tab_or_null("s9"),
  Table_S10     = with_rownames(tab_or_null("ev")),
  Table_S11     = tab_or_null("s11"),
  Table_S12     = tab_or_null("s12"),
  Table_S13     = with_rownames(tab_or_null("s13")),
  Table_S14     = tab_or_null("s14"),
  Table_S15     = tab_or_null("s15"),
  Table_S16     = tab_or_null("s16"),
  Table_S17     = tab_or_null("s17")))
for (nm in names(All_tables)) {
  assign(nm, All_tables[[nm]], envir = globalenv())
  if (interactive()) View(All_tables[[nm]], title = nm)
}
cat("\nTables kept in the R environment:", paste(names(All_tables), collapse = ", "), "\n")

## ---------------- Summary ----------------
cat("\n\n========================= SUMMARY =========================\n")
print(run_log[, c("Step", "Status")], row.names = FALSE)
if (any(run_log$Status == "FAILED")) {
  cat("\nErrors:\n"); failed <- run_log[run_log$Status == "FAILED", ]
  for (i in seq_len(nrow(failed))) cat(" -", failed$Step[i], ":", failed$Message[i], "\n")
}
cat("\nPDF files are in:", OUT_DIR, "\n")
print(list.files(OUT_DIR, pattern = "\\.pdf$"))
