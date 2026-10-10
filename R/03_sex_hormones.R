################################################################################
# 03_sex_hormones.R
# Targeted analysis of testosterone and sex hormone-binding globulin (SHBG).
#
# Five pre-specified measures: total, free, and bioavailable testosterone,
# percent free testosterone, and SHBG (dihydrotestosterone and epitestosterone
# excluded for insufficient data). Because the measures are few and correlated,
# nominal p < 0.05 is the primary threshold; Bonferroni (0.05 / 5) is reported
# for transparency.
#
# Output: Figure 3 - forest plot, male deletion and duplication carriers
#         Sex hormone results table (all strata)
################################################################################

source("R/00_config.R")
source("R/functions.R")

hormone_labels <- c(
  total        = "Total Testosterone",
  bioavailable = "Bioavailable Testosterone",
  free         = "Free Testosterone",
  percent_free = "Percent Free Testosterone",
  shbg         = "Sex Hormone Binding Globulin (SHBG)"
)

# ---- 1. Import and QC ------------------------------------------------------------
hormones <- read_strata(config$hormone_file, "testosterone") %>%
  map(~ .x %>%
        mutate(phenotype = as.character(phenotype),
               across(c(beta, SE, p, n_total), as.numeric)) %>%
        filter(!is.na(beta), !is.na(SE), !is.na(p), !is.na(n_total),
               !phenotype %in% params$hormone_excluded,
               n_total >= params$hormone_min_n) %>%
        mutate(
          measure  = coalesce(hormone_labels[phenotype], phenotype),
          lower_ci = beta - 1.96 * SE,
          upper_ci = beta + 1.96 * SE
        ))

# Bonferroni uses the largest number of measures tested in any stratum (5)
n_tests  <- max(map_int(hormones, nrow))
bonf_hormone <- params$alpha / n_tests

hormones <- map(hormones, ~ mutate(.x,
  sig_nominal    = p < params$alpha,
  sig_bonferroni = p < bonf_hormone
))

cat(sprintf("\nSex hormones: %d measures; nominal p < 0.05; Bonferroni p < %.3f\n",
            n_tests, bonf_hormone))
iwalk(hormones, ~ cat(sprintf("%-20s nominal: %d | Bonferroni: %d\n", .y,
                             sum(.x$sig_nominal), sum(.x$sig_bonferroni))))

# ---- 2. Figure 3: forest plot, males -------------------------------------------------
male_hormones <- bind_rows(
  mutate(hormones$male_deletion,    cnv = "Deletion"),
  mutate(hormones$male_duplication, cnv = "Duplication")
) %>%
  mutate(
    significance = factor(case_when(sig_bonferroni ~ "Bonferroni",
                                    sig_nominal    ~ "Nominal",
                                    TRUE           ~ "Not Significant"),
                          levels = c("Bonferroni", "Nominal", "Not Significant")),
    cnv = factor(cnv, levels = c("Deletion", "Duplication"))
  )

figure3 <- ggplot(male_hormones, aes(x = beta, y = measure, color = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.3, linewidth = 0.8) +
  geom_point(size = 4) +
  facet_grid(. ~ cnv, scales = "free_x") +
  scale_color_manual(values = c("Bonferroni" = "#4472C4", "Nominal" = "#9DC3E6",
                                "Not Significant" = "gray60"), drop = FALSE) +
  labs(x = "Effect Size (Beta, SD units)", y = NULL, color = "Significance") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14)
  )

save_figure(figure3, "Figure3_sex_hormones_males", 12, 6)

# ---- 3. Results table (all strata) ----------------------------------------------------
hormone_table <- imap_dfr(hormones, ~ .x %>%
  arrange(p) %>%
  transmute(
    Stratum = group_labels[[.y]],
    Measure = measure,
    `Beta (SD)` = round(beta, 3),
    `95% CI` = sprintf("%.3f to %.3f", lower_ci, upper_ci),
    `P-value` = signif(p, 3),
    `Nominal (p < 0.05)` = ifelse(sig_nominal, "Yes", "No"),
    `Bonferroni (p < 0.01)` = ifelse(sig_bonferroni, "Yes", "No"),
    `Sample Size` = n_total
  ))

write_xlsx(
  list(
    README = tibble(Information = c(
      "Sex hormone association results for STS copy number variants.",
      "Linear regression on inverse-normal-transformed values (beta in SD units),",
      "adjusted for age, HARE-designated genetic ancestry, and the top 10 genetic principal components.",
      sprintf("Measures with fewer than %d individuals excluded; DHT and epitestosterone excluded (insufficient data).",
              params$hormone_min_n),
      sprintf("Primary threshold: nominal p < 0.05. Bonferroni p < %.2f (0.05 / %d) reported for transparency.",
              bonf_hormone, n_tests),
      "Male strata are the primary analysis; female strata are exploratory."
    )),
    Results = hormone_table
  ),
  out_path("Sex_hormone_results.xlsx")
)
cat("Saved: Sex_hormone_results.xlsx\n")

writeLines(capture.output(sessionInfo()), out_path("sessionInfo_sex_hormones.txt"))
