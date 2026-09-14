# Public Release Checklist

- [x] Repository contains code and documentation only.
- [x] No CSV, TSV, Parquet, RDS, database, spreadsheet, image, PDF, archive, or patient-level output is tracked.
- [x] `outputs/`, `data/`, `derived/`, and `tmp/` are absent from the Git index.
- [x] No patient identifier values, passwords, tokens, connection strings, local absolute paths, personal emails, training identifiers, certificates, signatures, or signed DUA are present.
- [ ] A credentialed independent rerun reproduces 9,662 procedure-ICU pairs, 1,939 primary-cohort patients, and 322 events.
- [x] `CITATION.cff` validates and identifies Yang Li as given name Yang, family name Li, and names Yan Yan.
- [x] MIT License is present and explicitly does not cover MIMIC-IV data.
- [x] Git tag and release are both named `v1.0.0`.
- [x] Zenodo has supplied the version DOI: `10.5281/zenodo.22752940`.
- [x] Repository URL and archived DOI have been copied into the final manuscript.
