# Analytic Amendment Frozen on 2026-09-05

This document records a post-inspection analytic amendment. It is not a prospective registration and must not be described as one.

## Frozen changes

1. The cohort query exports all patients meeting the 48-hour landmark criteria before requiring an interpretable post-landmark CAM-ICU assessment.
2. The primary cohort excludes recorded CAM-ICU positivity during the first 48 hours, requires at least one interpretable CAM-ICU assessment after the landmark, and excludes patients with unsupported opioid-unit records.
3. The 0-48-hour IVME exposure includes recorded intravenous fentanyl, hydromorphone, and morphine. Infusions are clipped to the landmark window before rate-by-duration conversion.
4. The 10 mg and 25 mg cut points originated as empirical positive-dose quartiles in an earlier cohort, were fixed before the final cohort rerun, and are retained as analysis categories rather than clinical thresholds.
5. The primary four-group analysis uses multinomial generalized overlap weighting followed by entropy calibration. The nonparametric bootstrap refits the propensity model, calibration, and outcome model in each of 500 replicates.
6. Supporting analyses include raw-dose restricted cubic splines, binary exposure analysis, cross-fitted AIPW, missing-outcome bounds and tipping-point scenarios, E-values, upper-tail trimming, and exploratory threshold analysis.

## Reproducibility benchmarks

An authorized MIMIC-IV version 3.1 rerun should return 9,662 procedure-ICU pairs, 1,939 patients in the primary analysis cohort, and 322 post-landmark CAM-ICU-positive events.

This document contains no patient-level data, identifiers, database extracts, or credentialing information.
