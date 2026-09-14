# Data Access and Privacy

## Absolute rule

This public repository contains code and documentation only. No data are distributed with it.

MIMIC-IV version 3.1 is a credentialed-access resource hosted by PhysioNet. Every person who runs this code must independently complete the required training, obtain credentialed access, sign the applicable data use agreement, and comply with the PhysioNet Credentialed Health Data License 1.5.0.

## Files that must never enter this repository

- MIMIC-IV source tables or exports;
- cohort CSV, TSV, Parquet, RDS, database, or spreadsheet files;
- intermediate extraction tables or cached query results;
- patient-level covariates, timestamps, medication events, outcomes, model inputs, predictions, or bootstrap records;
- identifier values, including values from `subject_id`, `hadm_id`, and `stay_id` columns;
- personal CITI reports or certificates, learner or record identifiers, signatures, or signed data use agreements;
- database passwords, connection strings, access tokens, `.env` files, screenshots, or terminal logs;
- figures or tables that have not undergone disclosure review.

Running the pipeline writes working files below `outputs/`. The directory is ignored by Git, but `.gitignore` is not a security boundary. Investigators must keep the working copy in an access-controlled environment and inspect `git status` before every commit.

## Sharing results

Only author-reviewed aggregate results may be submitted with the manuscript. Patient-level source data and patient-level derived data cannot be provided by the authors. Eligible researchers should obtain the source data directly from PhysioNet and rerun the scripts.

The MIT License covers only the authors' code and documentation. It does not alter or supersede the license or data use agreement governing MIMIC-IV.
