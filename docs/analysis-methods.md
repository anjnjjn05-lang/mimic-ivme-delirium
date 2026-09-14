# Analysis Methods and Code Map

## Design

Retrospective 48-hour landmark cohort study using MIMIC-IV version 3.1. Exposure is the recorded cumulative intravenous morphine-equivalent dose during ICU hours 0-48. The outcome is the first recorded positive CAM-ICU assessment after the landmark and before day 7, ICU discharge, hospital discharge, or death, whichever occurs first.

## SQL

| File | Purpose |
|---|---|
| `sql/00_identification_counts.sql` | Counts cardiac procedure-ICU pairs and unique patients using the frozen ICD-9/10 procedure dictionary. |
| `sql/01_cohort_exposure_outcome.sql` | Selects the first eligible postoperative ICU stay, applies adult and 48-hour eligibility rules, extracts and standardizes fentanyl/hydromorphone/morphine from `inputevents`, constructs IVME, and evaluates CAM-ICU item 228332. |
| `sql/02_covariates.sql` | Extracts prior comorbidity flags and the most recent laboratory values in the 24 hours before ICU admission for a temporary in-session cohort table. |

Key medication item IDs are fentanyl 221744 and 225942, hydromorphone 221833, and morphine 225154. Infusion doses are clipped to the intersection of the recorded infusion and ICU hours 0-48 before rate-by-duration conversion. Unsupported unit combinations are flagged and excluded from the final model cohort.

## R workflow

| File | Purpose |
|---|---|
| `R/00_check_environment.R` | Verifies all required packages. |
| `R/01_build_analysis_dataset.R` | Runs the SQL cohort pipeline, writes restricted working files below ignored `outputs/`, and checks the frozen 9,662/1,939/322 benchmarks. |
| `R/02_extract_covariates.R` | Loads eligible identifiers into a PostgreSQL temporary table and runs `sql/02_covariates.sql`. |
| `R/03_overlap_weighting_bootstrap.R` | Fits multinomial propensity scores, generalized overlap weights, entropy calibration, weighted outcome models, 500-replicate full-pipeline bootstrap, E-values, missing-outcome bounds, and tipping-point analyses. |
| `R/04_spline_and_aipw.R` | Fits the positive-dose restricted cubic spline and repeated cross-fitted binary AIPW analysis. |
| `R/05_binary_sensitivity.R` | Performs the >25 mg versus <=25 mg overlap-weighted sensitivity analysis and its bootstrap. |
| `R/06_threshold_analysis.R` | Performs exploratory threshold-location analyses without changing the prespecified four-group analysis. |
| `R/07_figures_tables.R` | Generates the cohort flow diagram, balance displays, reported figures, and tables. |
| `R/08_validate_dose_extraction.R` | Independently requeries `inputevents`, recalculates the three opioid components, and emits disclosure-reviewed aggregate validation summaries. |

## Exposure groups

- G1: 0 mg IVME;
- G2: >0 to 10 mg;
- G3: >10 to 25 mg;
- G4: >25 mg.

The boundaries are analysis cut points and must not be interpreted as clinical prescribing thresholds.

## Output handling

Every generated file is written under `outputs/`. Patient-level working files, model objects, propensity scores, predictions, and bootstrap replicates are restricted derivatives and must not be distributed. Aggregate tables and figures require manual disclosure review before use in a manuscript or supplement. No generated output belongs in the public repository.
