#!/usr/bin/env Rscript
# Build exposure-preceding covariates for the full assessment-eligible cohort,
# then retain assessed patients for the primary analysis dataset.
suppressPackageStartupMessages({ library(DBI); library(RPostgres); library(data.table) })
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
args <- commandArgs(trailingOnly = TRUE)
root <- Sys.getenv('ANALYSIS_ROOT','..')
raw <- Sys.getenv('ANALYSIS_RAW',file.path(root,'outputs','restricted_patient_level'))
supp <- Sys.getenv('ANALYSIS_SUPP',file.path(root,'outputs','aggregate_review','supplementary'))
if (!length(args)) args <- c(file.path(raw,'assessment_eligible.csv'),
  file.path(raw,'baseline_covariates.csv'),file.path(supp,'covariate_availability.csv'),
  file.path(raw,'assessment_eligible_covariates.csv'))
if (length(args) != 4L) stop("Expected: <assessment_eligible.csv> <primary_covariates.csv> <audit.csv> <full_covariates.csv>")
needed <- c("MIMIC_DB_HOST", "MIMIC_DB_PORT", "MIMIC_DB_NAME", "MIMIC_DB_USER", "PGPASSWORD")
absent <- needed[!nzchar(Sys.getenv(needed))]
if (length(absent)) stop("Missing environment variables: ", paste(absent, collapse = ", "))

cohort <- fread(args[[1]])
id_cols <- c("subject_id", "hadm_id", "stay_id")
if (!all(id_cols %in% names(cohort))) stop("Final cohort lacks required identifiers.")
ids <- unique(cohort[, ..id_cols])
con <- dbConnect(RPostgres::Postgres(), host = Sys.getenv("MIMIC_DB_HOST"),
  port = as.integer(Sys.getenv("MIMIC_DB_PORT")), dbname = Sys.getenv("MIMIC_DB_NAME"),
  user = Sys.getenv("MIMIC_DB_USER"), password = Sys.getenv("PGPASSWORD"))
on.exit(dbDisconnect(con), add = TRUE)
dbWriteTable(con, "final_cohort_ids", as.data.frame(ids), temporary = TRUE, overwrite = TRUE)

sql <- "
WITH cohort AS (
  SELECT subject_id, hadm_id, stay_id FROM final_cohort_ids
),
prior_diagnoses AS (
  SELECT c.stay_id, d.icd_version, d.icd_code
  FROM cohort c
  JOIN mimiciv_hosp.admissions a_index ON a_index.hadm_id = c.hadm_id
  JOIN mimiciv_hosp.admissions a_prev ON a_prev.subject_id = c.subject_id AND a_prev.dischtime < a_index.admittime
  JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = a_prev.hadm_id
),
prior_flags AS (
  SELECT stay_id,
    bool_or((icd_version = 10 AND icd_code LIKE 'I50%') OR (icd_version = 9 AND icd_code LIKE '428%')) AS prior_chf,
    bool_or((icd_version = 10 AND (icd_code LIKE 'I21%' OR icd_code LIKE 'I22%' OR icd_code LIKE 'I252%')) OR (icd_version = 9 AND (icd_code LIKE '410%' OR icd_code LIKE '412%'))) AS prior_mi,
    bool_or((icd_version = 10 AND (icd_code LIKE 'I60%' OR icd_code LIKE 'I61%' OR icd_code LIKE 'I62%' OR icd_code LIKE 'I63%' OR icd_code LIKE 'I64%' OR icd_code LIKE 'I65%' OR icd_code LIKE 'I66%' OR icd_code LIKE 'I67%' OR icd_code LIKE 'I68%' OR icd_code LIKE 'I69%')) OR (icd_version = 9 AND (icd_code LIKE '430%' OR icd_code LIKE '431%' OR icd_code LIKE '432%' OR icd_code LIKE '433%' OR icd_code LIKE '434%' OR icd_code LIKE '435%' OR icd_code LIKE '436%' OR icd_code LIKE '437%' OR icd_code LIKE '438%'))) AS prior_cerebrovascular_disease,
    bool_or((icd_version = 10 AND (icd_code LIKE 'J40%' OR icd_code LIKE 'J41%' OR icd_code LIKE 'J42%' OR icd_code LIKE 'J43%' OR icd_code LIKE 'J44%' OR icd_code LIKE 'J45%' OR icd_code LIKE 'J46%' OR icd_code LIKE 'J47%')) OR (icd_version = 9 AND (icd_code LIKE '490%' OR icd_code LIKE '491%' OR icd_code LIKE '492%' OR icd_code LIKE '493%' OR icd_code LIKE '494%' OR icd_code LIKE '495%' OR icd_code LIKE '496%'))) AS prior_chronic_pulmonary_disease,
    bool_or((icd_version = 10 AND (icd_code LIKE 'E10%' OR icd_code LIKE 'E11%' OR icd_code LIKE 'E12%' OR icd_code LIKE 'E13%' OR icd_code LIKE 'E14%')) OR (icd_version = 9 AND icd_code LIKE '250%')) AS prior_diabetes,
    bool_or((icd_version = 10 AND (icd_code LIKE 'N18%' OR icd_code LIKE 'N19%')) OR (icd_version = 9 AND (icd_code LIKE '585%' OR icd_code LIKE '586%'))) AS prior_chronic_kidney_disease,
    bool_or((icd_version = 10 AND (icd_code LIKE 'K70%' OR icd_code LIKE 'K71%' OR icd_code LIKE 'K72%' OR icd_code LIKE 'K73%' OR icd_code LIKE 'K74%' OR icd_code LIKE 'K75%' OR icd_code LIKE 'K76%' OR icd_code LIKE 'K77%')) OR (icd_version = 9 AND icd_code LIKE '571%')) AS prior_liver_disease
  FROM prior_diagnoses GROUP BY stay_id
),
lab_rows AS (
  SELECT c.stay_id,
    CASE WHEN le.itemid = 50912 THEN 'creatinine' WHEN le.itemid = 51222 THEN 'hemoglobin'
         WHEN le.itemid = 50983 THEN 'sodium' WHEN le.itemid IN (51300,51301) THEN 'wbc'
         WHEN le.itemid IN (50809,50931) THEN 'glucose' WHEN le.itemid IN (50813,52442) THEN 'lactate'
         WHEN le.itemid = 50862 THEN 'albumin' END AS lab_name,
    le.valuenum,
    row_number() OVER (PARTITION BY c.stay_id, CASE WHEN le.itemid = 50912 THEN 'creatinine' WHEN le.itemid = 51222 THEN 'hemoglobin' WHEN le.itemid = 50983 THEN 'sodium' WHEN le.itemid IN (51300,51301) THEN 'wbc' WHEN le.itemid IN (50809,50931) THEN 'glucose' WHEN le.itemid IN (50813,52442) THEN 'lactate' WHEN le.itemid = 50862 THEN 'albumin' END ORDER BY le.charttime DESC, le.storetime DESC) AS rn
  FROM cohort c JOIN mimiciv_icu.icustays i ON i.stay_id = c.stay_id
  JOIN mimiciv_hosp.labevents le ON le.hadm_id = c.hadm_id
    AND le.charttime >= i.intime - interval '24 hour' AND le.charttime < i.intime
    AND le.itemid IN (50912,51222,50983,51300,51301,50809,50931,50813,52442,50862) AND le.valuenum IS NOT NULL
),
labs AS (
  SELECT stay_id,
    max(valuenum) FILTER (WHERE lab_name = 'creatinine' AND rn = 1) AS baseline_creatinine,
    max(valuenum) FILTER (WHERE lab_name = 'hemoglobin' AND rn = 1) AS baseline_hemoglobin,
    max(valuenum) FILTER (WHERE lab_name = 'sodium' AND rn = 1) AS baseline_sodium,
    max(valuenum) FILTER (WHERE lab_name = 'wbc' AND rn = 1) AS baseline_wbc,
    max(valuenum) FILTER (WHERE lab_name = 'glucose' AND rn = 1) AS baseline_glucose,
    max(valuenum) FILTER (WHERE lab_name = 'lactate' AND rn = 1) AS baseline_lactate,
    max(valuenum) FILTER (WHERE lab_name = 'albumin' AND rn = 1) AS baseline_albumin
  FROM lab_rows GROUP BY stay_id
)
SELECT c.stay_id,
  coalesce(p.prior_chf,false)::int AS prior_chf, coalesce(p.prior_mi,false)::int AS prior_mi,
  coalesce(p.prior_cerebrovascular_disease,false)::int AS prior_cerebrovascular_disease,
  coalesce(p.prior_chronic_pulmonary_disease,false)::int AS prior_chronic_pulmonary_disease,
  coalesce(p.prior_diabetes,false)::int AS prior_diabetes,
  coalesce(p.prior_chronic_kidney_disease,false)::int AS prior_chronic_kidney_disease,
  coalesce(p.prior_liver_disease,false)::int AS prior_liver_disease,
  l.baseline_creatinine, l.baseline_hemoglobin, l.baseline_sodium, l.baseline_wbc,
  l.baseline_glucose, l.baseline_lactate, l.baseline_albumin
FROM cohort c LEFT JOIN prior_flags p USING (stay_id) LEFT JOIN labs l USING (stay_id)
"
# The external SQL file is authoritative. The embedded string above is retained
# only to make historical comparisons straightforward and is not executed.
sql <- paste(
  readLines(file.path(root, "sql", "02_covariates.sql"), encoding = "UTF-8", warn = FALSE),
  collapse = "\n"
)
covariates <- as.data.table(dbGetQuery(con, sql))
model_data <- merge(cohort, covariates, by = "stay_id", all.x = TRUE, sort = FALSE)
if (nrow(model_data) != nrow(cohort) || anyDuplicated(model_data$stay_id)) stop("Covariate join changed the cohort.")
model_data[, race_group := fifelse(grepl("^WHITE", race), "White", fifelse(grepl("^BLACK", race), "Black", fifelse(grepl("^HISPANIC", race), "Hispanic", fifelse(grepl("^ASIAN", race), "Asian", "Other_or_unknown"))))]
model_data[, admission_group := fifelse(admission_type %chin% c("ELECTIVE", "SURGICAL SAME DAY ADMISSION"), "Elective", fifelse(grepl("URGENT|EW|EMER", admission_type), "Urgent_or_emergency", "Other"))]
if (!"valid_cam_assessment_n" %in% names(model_data)) stop("Assessment-eligible cohort lacks valid_cam_assessment_n.")
primary_data <- model_data[valid_cam_assessment_n > 0]
if (!nrow(primary_data)) stop("No assessed patients remained after covariate derivation.")
lab_cols <- grep("^baseline_", names(primary_data), value = TRUE)
core_cols <- c("age_at_admission", "gender", "race_group", "admission_group", "procedure_groups")
audit <- rbindlist(list(
  data.table(variable = core_cols, missing_n = sapply(core_cols, function(x) sum(is.na(primary_data[[x]]) | primary_data[[x]] == ""))),
  data.table(variable = lab_cols, missing_n = sapply(lab_cols, function(x) sum(is.na(primary_data[[x]]))))
))
audit[, total_n := nrow(primary_data)]
audit[, available_n := total_n - missing_n]
audit[, available_pct := round(100 * available_n / total_n, 1)]
fwrite(primary_data, args[[2]])
fwrite(audit, args[[3]])
fwrite(model_data, args[[4]])
cat("assessment_eligible_covariate_rows=", nrow(model_data), "\n", sep = "")
cat("primary_covariate_rows=", nrow(primary_data), "\n", sep = "")
cat("assessment_eligible_unique_stays=", uniqueN(model_data$stay_id), "\n", sep = "")
