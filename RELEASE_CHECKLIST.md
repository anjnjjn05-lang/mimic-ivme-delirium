# Public Release Checklist

- [ ] Repository contains code and documentation only.
- [ ] No CSV, TSV, Parquet, RDS, database, spreadsheet, image, PDF, archive, or patient-level output is tracked.
- [ ] `outputs/`, `data/`, `derived/`, and `tmp/` are absent from the Git index.
- [ ] No patient identifier values, passwords, tokens, connection strings, local absolute paths, personal emails, training identifiers, certificates, signatures, or signed DUA are present.
- [ ] A credentialed independent rerun reproduces 9,662 procedure-ICU pairs, 1,939 primary-cohort patients, and 322 events.
- [ ] `CITATION.cff` validates and names Li Yang and Yan Yan.
- [ ] MIT License is present and explicitly does not cover MIMIC-IV data.
- [ ] Git tag and release are both named `v1.0.0`.
- [ ] Zenodo has supplied the version DOI.
- [ ] Repository URL and archived DOI have been copied into the final manuscript.
