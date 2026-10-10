# Input data (not distributed)

The scripts expect two Excel workbooks in this folder. File names are set in `R/00_config.R`.

## 1. `sts del and dup labwas, t-was, phewas 29jun2023.xlsx`

There is one sheet per stratum and analysis, named `sts {del|dup} {males|females} {phewas|labwas}`, for example `sts del males phewas`.

**PheWAS sheets** (one row per phecode; logistic regression):

| Column | Description |
|---|---|
| `phecode`, `phenotype`, `category` | Phecode identifier, description, and disease category |
| `coef`, `se_coef` | Log odds ratio and its standard error |
| `lower95ci`, `upper95ci` | 95% CI on the log-odds scale |
| `p` | p-value (values such as `<2e-16` are accepted) |
| `case_with_phen`, `case_without_phen` | Carriers with and without the phenotype |
| `control_with_phen`, `control_without_phen` | Matched controls with and without the phenotype |

**LabWAS sheets** (one row per laboratory measure; linear regression on inverse-normal-transformed values):

| Column | Description |
|---|---|
| `phenotype` | Laboratory measure, e.g. `hdlc_median`, `calcium_mg_mean` |
| `beta`, `SE` | Effect in SD units and its standard error |
| `p` | p-value |
| `n_total` | Carriers plus controls with at least one measurement |

## 2. `testosterone data.xlsx`

There is one sheet per stratum, named `sts {del|dup} {males|females} testosterone`. The columns are the same as the LabWAS sheets. `phenotype` takes the values `total`, `bioavailable`, `free`, `percent_free`, `shbg`, `dht`, and `epi_t`.

To create synthetic files in this format, run `tests/make_mock_data.R`.
