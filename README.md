# Early ICU IVME and Delirium After Cardiac Surgery

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22752940.svg)](https://doi.org/10.5281/zenodo.22752940)

Reproducibility code for **Early ICU Intravenous Opioid Dose and Subsequent CAM-ICU Positivity After Cardiac Surgery: A Retrospective 48-Hour Landmark Cohort Study**.

Authors: Yang Li and Yan Yan. Yan Yan is the corresponding author and the credentialed PhysioNet user who accessed MIMIC-IV for this analysis.

The archived `v1.0.0` release is available from Zenodo at <https://doi.org/10.5281/zenodo.22752940>.

## Public-release boundary

This repository contains **code and documentation only**. It contains no MIMIC-IV source data, database exports, cohort files, intermediate tables, patient-level derived data, model inputs, predictions, or bootstrap records. The scripts create restricted working files locally; those outputs remain subject to the PhysioNet Credentialed Health Data License 1.5.0 and must never be committed or redistributed.

Anyone reproducing the study must independently complete the required training, obtain credentialed access, and sign the applicable data use agreement for [MIMIC-IV version 3.1](https://physionet.org/content/mimiciv/3.1/).

## Repository structure

```text
mimic-ivme-delirium/
├── README.md
├── LICENSE
├── CITATION.cff
├── DATA_ACCESS_AND_PRIVACY.md
├── RELEASE_CHECKLIST.md
├── run_all.ps1
├── sql/
│   ├── 00_identification_counts.sql
│   ├── 01_cohort_exposure_outcome.sql
│   └── 02_covariates.sql
├── R/
│   ├── 00_check_environment.R
│   ├── 01_build_analysis_dataset.R
│   ├── 02_extract_covariates.R
│   ├── 03_overlap_weighting_bootstrap.R
│   ├── 04_spline_and_aipw.R
│   ├── 05_binary_sensitivity.R
│   ├── 06_threshold_analysis.R
│   ├── 07_figures_tables.R
│   └── 08_validate_dose_extraction.R
├── docs/
│   ├── analysis-methods.md
│   └── analytic-amendment-2026-09-05.md
└── environment/
    └── sessionInfo.txt
```

The PostgreSQL cohort query is deliberately kept as one ordered CTE pipeline so that the cohort, exposure, and outcome definitions are evaluated from the same landmark population without exporting intermediate tables. Clearly labelled sections cover the procedure dictionary and cohort selection, 0-48-hour IVME exposure, post-landmark CAM-ICU outcome, and final analysis cohort. `sql/02_covariates.sql` is executed by `R/02_extract_covariates.R` using a temporary in-session identifier table; no identifier table is written to the repository.

## Requirements

- MIMIC-IV version 3.1, `hosp` and `icu` modules, loaded in PostgreSQL under `mimiciv_hosp` and `mimiciv_icu`.
- R 4.6.0.
- PowerShell 7 or Windows PowerShell 5.1.

| Package | Version |
|---|---:|
| DBI | 1.3.0 |
| RPostgres | 1.4.10 |
| data.table | 1.18.6.1 |
| WeightIt | 2.0.0 |
| cobalt | 5.0.0 |
| ggplot2 | 4.0.3 |
| sandwich | 3.1.3 |
| Hmisc | 5.2.6 |
| glmnet | 5.0 |

See `environment/sessionInfo.txt` for the frozen session record.

## Reproduction

1. Obtain independent MIMIC-IV v3.1 access and load the `hosp` and `icu` modules into PostgreSQL.
2. Clone this repository into a private, access-controlled analysis environment.
3. Install R 4.6.0 and the package versions listed above.
4. Set database credentials only in the current PowerShell session:

```powershell
$env:MIMIC_DB_HOST = "localhost"
$env:MIMIC_DB_PORT = "5432"
$env:MIMIC_DB_NAME = "mimiciv"
$env:MIMIC_DB_USER = "postgres"
$env:PGPASSWORD = Read-Host "PostgreSQL password" -MaskInput
```

5. Run `./run_all.ps1` from the repository root.

All generated files are written below `outputs/`, which is excluded by `.gitignore`. Do not move patient-level files into a tracked directory.

The pipeline stops if the frozen benchmark counts do not match:

- 9,662 literature-mapped procedure-ICU pairs;
- 1,939 patients in the primary analysis cohort;
- 322 post-landmark CAM-ICU-positive events.

A mismatch should trigger investigation of database version, schema naming, SQL semantics, or package versions; benchmark values must not be edited merely to make a run pass.

## Analysis sequence

1. Environment and package validation.
2. Procedure-ICU identification, cohort construction, 0-48-hour opioid extraction, IVME conversion, and CAM-ICU outcome construction.
3. Pre-exposure comorbidity and laboratory covariate extraction.
4. Multinomial propensity scores, generalized overlap weighting, entropy calibration, and 500-replicate full-pipeline bootstrap.
5. Restricted cubic spline dose-response analysis using knots at the 10th, 50th, and 90th percentiles; covariate-standardized predictions.
6. Cross-fitted AIPW, binary sensitivity analyses, missing-outcome bounds and tipping-point analysis, E-values, and exploratory threshold analysis.
7. Table and figure generation plus independent dose-extraction audit.

The estimates are observational associations. The 10 mg and 25 mg group boundaries are analysis cut points, not validated clinical prescribing thresholds.

## Citation

Please cite `CITATION.cff`, the associated manuscript, and the four MIMIC-IV and PhysioNet references used in the manuscript:

1. Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic health record dataset. *Scientific Data*. 2023;10:1. <https://doi.org/10.1038/s41597-022-01899-x>.
2. Johnson A, Bulgarelli L, Pollard T, et al. MIMIC-IV version 3.1. PhysioNet; 2024. <https://doi.org/10.13026/kpb9-mt58>.
3. Pollard T, Moody BE, Lehman L, Gow BJ, Fernandes C, Xie C, et al. PhysioNet as a global platform for biomedical research. *Nature Health*. 2026;1:792-795. <https://doi.org/10.1038/s44360-026-00096-z>.
4. Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PC, Mark RG, et al. PhysioBank, PhysioToolkit, and PhysioNet. *Circulation*. 2000;101:e215-e220. <https://doi.org/10.1161/01.CIR.101.23.e215>.

## License

The authors' code and documentation are released under the MIT License. This license does not grant any right to redistribute MIMIC-IV data or patient-level derivatives.
