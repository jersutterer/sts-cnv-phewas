################################################################################
# make_mock_data.R
# Writes SYNTHETIC input files in the expected format so the pipeline can be
# run end to end without access to Million Veteran Program data.
# The values are random and have no scientific meaning.
#
# Usage (from the repository root):
#   Rscript tests/make_mock_data.R
#   Rscript R/01_phewas.R; Rscript R/02_labwas.R; Rscript R/03_sex_hormones.R
################################################################################

suppressPackageStartupMessages({ library(tidyverse); library(writexl) })
set.seed(2026)
dir.create("data", showWarnings = FALSE)

strata <- c("sts del females", "sts del males", "sts dup females", "sts dup males")

mock_phewas <- function(n = 1875) {
  cats <- c("circulatory system", "dermatologic", "genitourinary", "sense organs",
            "endocrine/metabolic", "hematopoietic", "digestive", "neoplasms",
            "mental disorders", "respiratory", "infectious diseases", "musculoskeletal")
  coef <- rnorm(n, 0, 0.3)
  coef[1:5] <- c(6.4, 5.1, 2.1, -0.75, 1.7)   # a few strong signals
  se <- runif(n, 0.15, 0.6)
  p <- 2 * pnorm(-abs(coef / se))
  tibble(
    phecode = sprintf("phe_%d", 100 + seq_len(n)), phecode_mod = phecode,
    phenotype = c("Ichthyosis congenita", "Congenital anomalies of skin",
                  "Corneal opacity", "Coronary atherosclerosis", "Atrial fibrillation",
                  sprintf("Phenotype %d", 6:n)),
    category = sample(cats, n, replace = TRUE),
    coef = coef, se_coef = se,
    lower95ci = coef - 1.96 * se, upper95ci = coef + 1.96 * se,
    p = ifelse(p < 2e-16, "<2e-16", as.character(p)),
    case_with_phen = rpois(n, 12), case_without_phen = 241 - case_with_phen,
    control_with_phen = rpois(n, 150), control_without_phen = 4820 - control_with_phen
  )
}

mock_labwas <- function(n = 552) {
  base <- c("hdlc", "trig", "calcium_mg", "egfr", "creat_bsp", "rdw", "platelet",
            "lymph_fra", "mono_fra", "psa", "sodium_bsp", "a1c", "wbc", "alt",
            "ferritin", "crp_l", "inr", "troponin", "uricacid", "mcv")
  suf <- c("_fst", "_lst", "_mean", "_median", "_min", "_max")
  labs <- head(c(outer(base, suf, paste0), sprintf("misc_%d", 1:500)), n)
  beta <- rnorm(n, 0, 0.08); beta[1:3] <- c(-0.64, 0.38, -0.62)
  se <- runif(n, 0.05, 0.12)
  tibble(phenotype = labs, beta = beta, SE = se,
         p = 2 * pnorm(-abs(beta / se)), n_total = sample(20:5000, n, replace = TRUE))
}

mock_hormones <- function() {
  ph <- c("total", "bioavailable", "free", "percent_free", "shbg", "dht", "epi_t")
  beta <- c(0.05, 0.2, 0.35, 0.3, -0.62, 0.1, NA)
  se <- c(0.12, 0.12, 0.11, 0.12, 0.14, 0.3, NA)
  tibble(phenotype = ph, beta = beta, SE = se,
         p = 2 * pnorm(-abs(beta / se)), n_total = c(400, 120, 300, 110, 150, 15, 5))
}

write_xlsx(c(set_names(map(strata, ~ mock_phewas()), paste(strata, "phewas")),
             set_names(map(strata, ~ mock_labwas()), paste(strata, "labwas"))),
           "data/sts del and dup labwas, t-was, phewas 29jun2023.xlsx")
write_xlsx(set_names(map(strata, ~ mock_hormones()), paste(strata, "testosterone")),
           "data/testosterone data.xlsx")
cat("Synthetic input files written to data/\n")
