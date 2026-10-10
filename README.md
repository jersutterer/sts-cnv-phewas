# STS copy number variants: PheWAS, LabWAS, and sex hormone analyses

Analysis code accompanying:

> Sutterer JT, Gorman BR, Teerlink CC, et al. Phenome- and laboratory-wide association study of steroid sulfatase copy number variants in 650,427 participants of the Million Veteran Program reveals multi-system features of X-linked ichthyosis. *eBioMedicine* (submitted).

## What this repository contains

The code here starts from association **summary statistics** (odds ratios, betas, standard errors, p-values, and counts) and produces the multiple-testing correction, figures, and tables reported in the paper.

Cohort extraction, copy number variant calling, case-control matching, and the regression models were run inside the VA Informatics and Computing Infrastructure (VINCI) and are **not** included. That code is available from the corresponding authors on request, subject to VA approval.

| Script | Produces |
|---|---|
| `R/01_phewas.R` | Figure 1, Figure S3, Table 1, Table S3 |
| `R/02_labwas.R` | Figure 2, Figure S4, Table S4 |
| `R/03_sex_hormones.R` | Figure 3, sex hormone results table |
| `R/00_config.R` | Input paths and analysis parameters (QC thresholds, alpha) |
| `R/functions.R` | Shared helpers (multiple-testing correction, Manhattan plot) |
| `tests/make_mock_data.R` | Synthetic input files for testing the pipeline |

## Methods implemented

- **QC filters:** phecodes with at least 5 carrier cases; laboratory measures with at least 30 individuals; hormone measures with at least 25 individuals (dihydrotestosterone and epitestosterone excluded for insufficient data).
- **PheWAS and LabWAS:** the Bonferroni threshold is 0.05 divided by the number of tests passing QC, computed separately for each sex × CNV stratum. FDR uses the Benjamini–Hochberg method (q < 0.05).
- **Sex hormones:** the primary threshold is nominal p < 0.05. Bonferroni (0.05/5) is also reported.
- **Strata:** male strata are the primary analysis; female strata are reported as exploratory.

## Data

Individual-level Million Veteran Program data cannot be shared publicly. Qualified researchers can obtain them through the MVP data access process (https://www.mvp.va.gov/). The summary statistics needed to reproduce every figure and table are provided as Tables S3 and S4 of the paper. `data/README.md` describes the input format the scripts expect.

## Running

Requires R (4.3 or later was tested) and the packages tidyverse, readxl, writexl, ggrepel, cowplot, flextable, and officer.

```r
# from the repository root
source("R/01_phewas.R")
source("R/02_labwas.R")
source("R/03_sex_hormones.R")
```

Outputs are written to `output/`, together with `sessionInfo` files recording package versions.

To test the pipeline without MVP data, generate synthetic inputs first. The resulting values are random and meaningless.

```r
source("tests/make_mock_data.R")
```

## Citation and licence

Please cite the paper above (see `CITATION.cff`). The code is released under the MIT licence.
