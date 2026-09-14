suppressPackageStartupMessages({library(DBI);library(RPostgres);library(data.table)})
options(warn=1)
try(Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8'),silent=TRUE)
root <- Sys.getenv('ANALYSIS_ROOT','..')
raw <- Sys.getenv('ANALYSIS_RAW',file.path(root,'outputs','restricted_patient_level'))
supp <- Sys.getenv('ANALYSIS_SUPP',file.path(root,'outputs','aggregate_review','supplementary'))
needed <- c('MIMIC_DB_HOST','MIMIC_DB_PORT','MIMIC_DB_NAME','MIMIC_DB_USER','PGPASSWORD')
absent <- needed[!nzchar(Sys.getenv(needed))]
if (length(absent)) stop('Missing environment variables: ', paste(absent, collapse=', '))
dir.create(raw,recursive=TRUE,showWarnings=FALSE);dir.create(supp,recursive=TRUE,showWarnings=FALSE)
con <- dbConnect(Postgres(),host=Sys.getenv('MIMIC_DB_HOST'),port=as.integer(Sys.getenv('MIMIC_DB_PORT')),
 dbname=Sys.getenv('MIMIC_DB_NAME'),user=Sys.getenv('MIMIC_DB_USER'),password=Sys.getenv('PGPASSWORD'))
identification_sql <- paste(readLines(file.path(root,'sql','00_identification_counts.sql'),encoding='UTF-8',warn=FALSE),collapse='\n')
identification_counts <- as.data.table(dbGetQuery(con,identification_sql))
stopifnot(nrow(identification_counts)==1,
          identification_counts$procedure_icu_pairs >= identification_counts$unique_patients)
fwrite(identification_counts,file.path(supp,'identification_counts.csv'))
sql <- paste(readLines(file.path(root,'sql','01_cohort_exposure_outcome.sql'),encoding='UTF-8',warn=FALSE),collapse='\n')
d <- as.data.table(dbGetQuery(con,sql))
dbDisconnect(con)
stopifnot(!anyDuplicated(d$subject_id),!anyDuplicated(d$stay_id))
fwrite(d,file.path(raw,'landmark_all.csv'))
flow <- data.table(step=c('Adult first-per-person cardiac ICU; alive and in ICU at 48h; no recorded prior dementia',
 'No recorded positive CAM-ICU in first 48h','At least one valid post-landmark CAM-ICU','Complete interpretable opioid dose records'),
 n=c(nrow(d),nrow(d[early_cam_positive==0]),nrow(d[early_cam_positive==0 & valid_cam_assessment_n>0]),
 nrow(d[early_cam_positive==0 & valid_cam_assessment_n>0 & unsupported_event_n==0])))
fwrite(flow,file.path(supp,'cohort_flow.csv'))
eligible <- d[early_cam_positive==0 & unsupported_event_n==0]
fwrite(eligible,file.path(raw,'assessment_eligible.csv'))
final <- eligible[valid_cam_assessment_n>0]
stopifnot(!anyNA(final$incident_cam_positive))
expected_pairs <- as.integer(Sys.getenv('EXPECTED_PROCEDURE_ICU_PAIRS','9662'))
expected_n <- as.integer(Sys.getenv('EXPECTED_FINAL_N','1939'))
expected_events <- as.integer(Sys.getenv('EXPECTED_EVENTS','322'))
if (identification_counts$procedure_icu_pairs != expected_pairs) {
  stop('Procedure-ICU pair benchmark mismatch: expected ', expected_pairs,
       ', observed ', identification_counts$procedure_icu_pairs)
}
if (nrow(final) != expected_n || sum(final$incident_cam_positive) != expected_events) {
  stop('Primary cohort benchmark mismatch: expected n/events ', expected_n, '/', expected_events,
       ', observed ', nrow(final), '/', sum(final$incident_cam_positive))
}
fwrite(final,file.path(raw,'final_cohort.csv'))
print(flow)
cat('Unsupported-event patients:',sum(d$unsupported_event_n>0),'\n')
