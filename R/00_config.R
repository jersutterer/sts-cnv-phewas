################################################################################
# 00_config.R
# Shared settings for the STS copy number variant PheWAS, LabWAS, and
# sex hormone analyses. Edit the input paths below before running.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, tidyr, ggplot2, readr, stringr, purrr
  library(readxl)
  library(writexl)
  library(ggrepel)
  library(cowplot)
  library(flextable)
  library(officer)
})

# ---- Input files --------------------------------------------------------------
# Association summary statistics exported from the VA Informatics and Computing
# Infrastructure (VINCI). These files are not distributed with this repository;
# see data/README.md for the expected format.
config <- list(
  phewas_labwas_file = "data/sts del and dup labwas, t-was, phewas 29jun2023.xlsx",
  hormone_file       = "data/testosterone data.xlsx",
  output_dir         = "output"
)

# Sheet names for each sex x CNV stratum, in the order used throughout.
groups <- c("female_deletion", "male_deletion",
            "female_duplication", "male_duplication")

group_labels <- c(
  female_deletion    = "Female Deletion",
  male_deletion      = "Male Deletion",
  female_duplication = "Female Duplication",
  male_duplication   = "Male Duplication"
)

sheet_name <- function(group, analysis) {
  # e.g. ("male_deletion", "phewas") -> "sts del males phewas"
  cnv <- if (grepl("deletion", group)) "del" else "dup"
  sex <- if (grepl("^female", group)) "females" else "males"
  sprintf("sts %s %s %s", cnv, sex, analysis)
}

# ---- Analysis parameters (as reported in Methods and Table S1) ---------------
params <- list(
  phewas_min_cases     = 5,     # phecodes with fewer than 5 carrier cases excluded
  labwas_min_n         = 30,    # laboratory measures with fewer than 30 individuals excluded
  hormone_min_n        = 25,    # hormone measures with fewer than 25 individuals excluded
  hormone_excluded     = c("dht", "epi_t"),  # insufficient data
  fdr_q                = 0.05,
  alpha                = 0.05,
  extreme_p            = 1e-15  # flagged with an asterisk in Figure 1
)

dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- function(...) file.path(config$output_dir, ...)
