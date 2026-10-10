################################################################################
# 01_phewas.R
# Phenome-wide association study (PheWAS) of STS copy number variants.
#
# Input : logistic-regression summary statistics for 1875 phecodes per stratum
#         (carriers vs 1:20 matched controls; see Methods)
# Output: Figure 1  - Manhattan plot, male deletion hemizygotes
#         Figure S3 - Manhattan plot, male duplication carriers
#         Table 1   - FDR-significant phecodes, male deletion hemizygotes
#         Table S3  - complete PheWAS results, all strata
################################################################################

source("R/00_config.R")
source("R/functions.R")

# ---- 1. Import -----------------------------------------------------------------
phewas_cols <- c("phecode", "phenotype", "category", "coef", "se_coef",
                 "lower95ci", "upper95ci", "p", "case_with_phen",
                 "case_without_phen", "control_with_phen", "control_without_phen")

phewas_raw <- read_strata(config$phewas_labwas_file, "phewas") %>%
  map(~ .x %>%
        select(any_of(c(phewas_cols, "group"))) %>%
        mutate(
          across(c(phecode, phenotype, category), as.character),
          p = clean_pvalue(p),
          across(c(coef, se_coef, lower95ci, upper95ci, case_with_phen,
                   case_without_phen, control_with_phen, control_without_phen),
                 as.numeric)
        ))

# ---- 2. QC: at least 5 carrier cases per phecode --------------------------------
# ---- 3. Multiple-testing correction within each stratum -------------------------
phewas <- phewas_raw %>%
  map(~ .x %>%
        filter(case_with_phen >= params$phewas_min_cases) %>%
        add_multiple_testing() %>%
        mutate(
          OR     = exp(coef),
          OR_lo  = exp(lower95ci),   # CIs are exported on the log-odds scale
          OR_hi  = exp(upper95ci),
          extreme_p = p < params$extreme_p
        ))

cat("\nPheWAS (phecodes with >=", params$phewas_min_cases, "carrier cases)\n")
report_strata(phewas, "phecodes")

# ---- 4. Figures -------------------------------------------------------------------
# Labels: all Bonferroni-significant phecodes, plus two FDR-significant phecodes
# discussed in the text (corneal oedema, atrial fibrillation).
phewas_labels <- function(df) {
  df %>%
    filter(sig_bonferroni |
             (sig_fdr & grepl("corneal edema|atrial fibrillation", phenotype,
                              ignore.case = TRUE))) %>%
    distinct(phecode, .keep_all = TRUE) %>%
    mutate(label_text = ifelse(extreme_p, paste0(phenotype, "*"), phenotype))
}

phewas_figure <- function(df) {
  lab <- phewas_labels(df)
  manhattan_plot(
    df, category_col = "category", effect_col = "coef", label_data = lab,
    shape_labels = c(up = "Risk Increasing (OR > 1)", down = "Risk Decreasing (OR < 1)"),
    x_title = "Disease Category",
    caption = if (any(lab$extreme_p)) expression("*" ~ p < 1 %*% 10^-15) else NULL
  )
}

save_figure(phewas_figure(phewas$male_deletion),    "Figure1_PheWAS_male_deletion",     10, 7)
save_figure(phewas_figure(phewas$male_duplication), "FigureS3_PheWAS_male_duplication", 10, 7)

# ---- 5. Table 1: FDR-significant phecodes, male deletion -------------------------
table1 <- phewas$male_deletion %>%
  filter(sig_fdr) %>%
  arrange(p) %>%
  transmute(
    Phecode      = phecode,
    Phenotype    = phenotype,
    Category     = category,
    Direction    = ifelse(coef > 0, "↑", "↓"),
    `OR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", OR, OR_lo, OR_hi),
    `P-value`    = fmt_p(p),
    Sig          = ifelse(sig_bonferroni, "a", "b"),
    `Cases (n)`  = case_with_phen,
    `Controls (n)` = control_with_phen
  )

ft <- flextable(table1) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 8, part = "all") %>%
  bold(part = "header") %>%
  align(j = c("Direction", "Sig"), align = "center", part = "all") %>%
  align(j = c("OR (95% CI)", "P-value", "Cases (n)", "Controls (n)"),
        align = "right", part = "all") %>%
  border_remove() %>%
  hline_top(border = fp_border(width = 1.5), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "header") %>%
  hline_bottom(border = fp_border(width = 1.5), part = "body") %>%
  add_footer_lines(c(
    "OR = odds ratio. Direction: ↑ increased risk (OR > 1); ↓ decreased risk (OR < 1).",
    sprintf("a Bonferroni-significant (p < %.2e); b FDR-significant (q < 0.05).",
            bonf_threshold(phewas$male_deletion)),
    "Cases (n) and Controls (n): number of individuals with each phenotype."
  )) %>%
  fontsize(size = 8, part = "footer") %>%
  autofit()

read_docx() %>%
  body_add_flextable(ft) %>%
  print(target = out_path("Table1_PheWAS_male_deletion.docx"))
write_csv(table1, out_path("Table1_PheWAS_male_deletion.csv"))
cat("Saved: Table1_PheWAS_male_deletion (.docx, .csv)\n")

# ---- 6. Table S3: complete results, all strata ------------------------------------
supp_table <- function(df) {
  df %>%
    arrange(p) %>%
    transmute(
      Phecode = phecode, Phenotype = phenotype, Category = category,
      `Effect Direction` = ifelse(coef > 0, "Risk Increasing", "Risk Decreasing"),
      `Beta (log OR)` = round(coef, 4),
      OR = round(OR, 3), `OR 95% CI lower` = round(OR_lo, 3), `OR 95% CI upper` = round(OR_hi, 3),
      `P-value` = signif(p, 3), `FDR Q-value` = signif(fdr_q, 3),
      `Bonferroni Significant` = ifelse(sig_bonferroni, "Yes", "No"),
      `FDR Significant` = ifelse(sig_fdr, "Yes", "No"),
      `Cases with Phenotype` = case_with_phen,
      `Cases without Phenotype` = case_without_phen,
      `Controls with Phenotype` = control_with_phen,
      `Controls without Phenotype` = control_without_phen
    )
}

readme_s3 <- tibble(Information = c(
  "Table S3. Complete phenome-wide association results for STS copy number variants.",
  "",
  "Carriers were matched 1:20 to controls on age, sex, and HARE-designated genetic ancestry.",
  "Logistic regression adjusted for age, HARE-designated genetic ancestry, and the top 10 genetic principal components.",
  sprintf("Phecodes with fewer than %d carrier cases were excluded.", params$phewas_min_cases),
  "Bonferroni threshold = 0.05 / number of phecodes tested in each stratum; FDR by Benjamini-Hochberg.",
  "Male strata are the primary analysis; female strata are exploratory (insufficient sample size).",
  "",
  "Bonferroni thresholds by stratum:",
  imap_chr(phewas, ~ sprintf("  %s: %d phecodes, p < %.2e",
                             group_labels[[.y]], nrow(.x), bonf_threshold(.x))),
  "",
  "Columns: Beta (log OR) = logistic regression coefficient; OR = odds ratio;",
  "P-value = nominal p-value; FDR Q-value = Benjamini-Hochberg adjusted p-value;",
  "Cases/Controls with/without Phenotype = counts of carriers and matched controls."
))

write_xlsx(
  c(list(README = readme_s3),
    set_names(map(phewas, supp_table), gsub(" ", "_", group_labels[names(phewas)]))),
  out_path("TableS3_PheWAS_complete_results.xlsx")
)
cat("Saved: TableS3_PheWAS_complete_results.xlsx\n")

writeLines(capture.output(sessionInfo()), out_path("sessionInfo_phewas.txt"))
