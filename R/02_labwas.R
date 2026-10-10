################################################################################
# 02_labwas.R
# Laboratory-wide association study (LabWAS) of STS copy number variants.
#
# Input : linear-regression summary statistics (inverse-normal-transformed
#         values; beta in SD units) for 552 laboratory measures per stratum
# Output: Figure 2  - Manhattan plot, male deletion hemizygotes
#         Figure S4 - Manhattan plot, male duplication carriers
#         Table S4  - complete LabWAS results, all strata
################################################################################

source("R/00_config.R")
source("R/functions.R")

# ---- 1. Laboratory categories (13) ------------------------------------------------
# Each measure is assigned to the first category whose prefix it matches.
lab_categories <- list(
  "Lipids"             = c("hdlc", "ldlc", "totchol", "trig", "lpa"),
  "CBC-RBC"            = c("rbc", "hemoglobin", "hct", "mcv", "mch", "mchc", "rdw"),
  "CBC-WBC"            = c("wbc", "neut_abs", "neut_fra", "lymph_abs", "lymph_fra",
                           "mono_abs", "mono_fra", "eos_abs", "eos_fra",
                           "baso_abs", "baso_fra"),
  "Platelets"          = c("platelet", "mpv"),
  "Kidney Function"    = c("creat_bsp", "bun_bsp", "egfr"),
  "Electrolytes"       = c("sodium_bsp", "potas_bsp", "chloride_bsp", "calcium_mg",
                           "calcium_mmol", "bicarbonate", "mg_bsp"),
  "Liver Function"     = c("alt", "ast", "albumin"),
  "Cardiac Biomarkers" = c("troponin", "troponini", "troponint", "bnp", "probnp",
                           "totalck", "ckmb_abs", "ckmb_fra"),
  "Metabolic/Diabetes" = c("glucose", "glucose_fs", "glucose_fasting", "a1c",
                           "c_peptide", "bohb", "gad"),
  "Iron Studies"       = c("iron", "ferritin", "tsat", "stfr", "total_iron"),
  "Inflammation"       = c("crp_l", "crp_dl", "esr", "rf", "ccp"),
  "Coagulation"        = c("pt", "inr"),
  "Other"              = c("uricacid", "psa", "hcv_rnavl", "hcv_bdnavl")
)

assign_lab_category <- function(lab) {
  for (cat_name in names(lab_categories)) {
    if (any(startsWith(lab, lab_categories[[cat_name]]))) return(cat_name)
  }
  "Other"
}

# Display names for labelled points in Figures 2 and S4
lab_display_names <- c(
  hdlc = "HDL-C", ldlc = "LDL-C", totchol = "Total Cholesterol", trig = "Triglycerides",
  calcium_mg = "Calcium", calcium_mmol = "Calcium", platelet = "Platelets",
  creat_bsp = "Creatinine", egfr = "eGFR", psa = "PSA", rdw = "RDW",
  lymph_fra = "Lymphocyte %", lymph_abs = "Lymphocyte (abs)", mono_fra = "Monocyte %",
  mono_abs = "Monocyte (abs)", eos_fra = "Eosinophil %", eos_abs = "Eosinophil (abs)",
  neut_fra = "Neutrophil %", neut_abs = "Neutrophil (abs)", baso_fra = "Basophil %",
  baso_abs = "Basophil (abs)", totalck = "Total CK", a1c = "HbA1c",
  sodium_bsp = "Sodium", potas_bsp = "Potassium", chloride_bsp = "Chloride",
  wbc = "WBC", rbc = "RBC", hemoglobin = "Haemoglobin", hct = "Haematocrit",
  mcv = "MCV", mch = "MCH", mchc = "MCHC", uricacid = "Uric Acid",
  glucose = "Glucose", albumin = "Albumin", alt = "ALT", ast = "AST"
)

# ---- 2. Import -------------------------------------------------------------------
labwas_raw <- read_strata(config$phewas_labwas_file, "labwas") %>%
  map(~ .x %>%
        mutate(phenotype = as.character(phenotype),
               across(c(beta, SE, p, n_total), as.numeric)) %>%
        filter(!is.na(beta), !is.na(SE), !is.na(p), !is.na(n_total)) %>%
        mutate(lab_category = map_chr(phenotype, assign_lab_category)))

# ---- 3. QC (>= 30 individuals) and multiple-testing correction ---------------------
labwas <- labwas_raw %>%
  map(~ .x %>%
        filter(n_total >= params$labwas_min_n) %>%
        add_multiple_testing())

cat("\nLabWAS (measures with >=", params$labwas_min_n, "individuals)\n")
report_strata(labwas, "measures")

# ---- 4. Figures --------------------------------------------------------------------
# Label each Bonferroni-significant laboratory test once (summary statistic
# with the lowest p-value, e.g. first / last / mean / median value).
labwas_labels <- function(df) {
  df %>%
    filter(sig_bonferroni) %>%
    mutate(base_lab = sub("(_fst|_lst|_mean|_median|_min|_max)$", "", phenotype)) %>%
    group_by(base_lab) %>%
    slice_min(p, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label_text = coalesce(lab_display_names[base_lab],
                                 tools::toTitleCase(gsub("_", " ", base_lab))))
}

labwas_figure <- function(df) {
  manhattan_plot(
    df, category_col = "lab_category", effect_col = "beta",
    label_data = labwas_labels(df),
    shape_labels = c(up = "Value Increased", down = "Value Decreased"),
    x_title = "Laboratory Category"
  )
}

save_figure(labwas_figure(labwas$male_deletion),    "Figure2_LabWAS_male_deletion",     10, 7)
save_figure(labwas_figure(labwas$male_duplication), "FigureS4_LabWAS_male_duplication", 10, 7)

# ---- 5. Table S4: complete results, all strata -------------------------------------
supp_table <- function(df) {
  df %>%
    arrange(p) %>%
    transmute(
      `Lab Test` = phenotype, `Lab Category` = lab_category,
      `Effect Direction` = ifelse(beta > 0, "Increased", "Decreased"),
      `Beta (SD)` = round(beta, 4), `Standard Error` = round(SE, 4),
      `Beta 95% CI lower` = round(beta - 1.96 * SE, 4),
      `Beta 95% CI upper` = round(beta + 1.96 * SE, 4),
      `P-value` = signif(p, 3), `FDR Q-value` = signif(fdr_q, 3),
      `Bonferroni Significant` = ifelse(sig_bonferroni, "Yes", "No"),
      `FDR Significant` = ifelse(sig_fdr, "Yes", "No"),
      `Sample Size` = n_total
    )
}

readme_s4 <- tibble(Information = c(
  "Table S4. Laboratory-wide association results for STS copy number variants.",
  "",
  "Carriers were matched 1:20 to controls on age, sex, and HARE-designated genetic ancestry.",
  "Linear regression on inverse-normal-transformed values (beta in SD units), adjusted for age,",
  "HARE-designated genetic ancestry, and the top 10 genetic principal components.",
  sprintf("Laboratory measures with fewer than %d individuals were excluded.", params$labwas_min_n),
  "Bonferroni threshold = 0.05 / number of measures tested in each stratum; FDR by Benjamini-Hochberg.",
  "Male strata are the primary analysis; female strata are exploratory (insufficient sample size).",
  "",
  "Bonferroni thresholds by stratum:",
  imap_chr(labwas, ~ sprintf("  %s: %d measures, p < %.2e",
                             group_labels[[.y]], nrow(.x), bonf_threshold(.x))),
  "",
  "Lab Test suffixes: _fst first, _lst last, _mean, _median, _min, _max value per person.",
  "Positive beta = higher value in carriers than in matched controls.",
  "Sample Size = carriers plus controls with at least one measurement."
))

write_xlsx(
  c(list(README = readme_s4),
    set_names(map(labwas, supp_table), gsub(" ", "_", group_labels[names(labwas)]))),
  out_path("TableS4_LabWAS_complete_results.xlsx")
)
cat("Saved: TableS4_LabWAS_complete_results.xlsx\n")

writeLines(capture.output(sessionInfo()), out_path("sessionInfo_labwas.txt"))
