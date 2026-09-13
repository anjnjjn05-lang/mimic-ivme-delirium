# Independent source-table validation; generates aggregate outputs only.
# This script does not alter the frozen analysis datasets or model results.
suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(data.table)
})

options(warn = 1)
try(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"), silent = TRUE)

root <- Sys.getenv("ANALYSIS_ROOT")
if (!nzchar(root)) stop("ANALYSIS_ROOT is required")
raw_dir <- Sys.getenv("ANALYSIS_RAW")
if (!nzchar(raw_dir)) raw_dir <- file.path(root, "outputs", "restricted_patient_level")
out_dir <- Sys.getenv("ANALYSIS_VALIDATION")
if (!nzchar(out_dir)) out_dir <- file.path(root, "outputs", "aggregate_review", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

needed <- c("MIMIC_DB_HOST", "MIMIC_DB_PORT", "MIMIC_DB_NAME", "MIMIC_DB_USER", "PGPASSWORD")
absent <- needed[!nzchar(Sys.getenv(needed))]
if (length(absent)) stop("Missing environment variables: ", paste(absent, collapse = ", "))

final <- fread(file.path(raw_dir, "final_cohort.csv"))
expected_n <- as.integer(Sys.getenv("EXPECTED_FINAL_N", "1939"))
stopifnot(nrow(final) == expected_n, !anyDuplicated(final$stay_id))
stay_values <- paste0("(", final$stay_id, ")", collapse = ",")

con <- dbConnect(
  Postgres(),
  host = Sys.getenv("MIMIC_DB_HOST"),
  port = as.integer(Sys.getenv("MIMIC_DB_PORT")),
  dbname = Sys.getenv("MIMIC_DB_NAME"),
  user = Sys.getenv("MIMIC_DB_USER"),
  password = Sys.getenv("PGPASSWORD")
)
on.exit(dbDisconnect(con), add = TRUE)

sql <- sprintf("\
WITH selected_stays(stay_id) AS (VALUES %s)
SELECT ie.stay_id, ie.orderid, ie.linkorderid, ie.itemid,
       di.label AS item_label,
       ie.starttime, ie.endtime,
       ie.amount, ie.amountuom, ie.rate, ie.rateuom,
       ie.patientweight, ie.statusdescription,
       CASE
         WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mcg/hour'
              AND ie.itemid IN (221744, 225942)
           THEN ie.rate * extract(epoch FROM (least(ie.endtime, i.intime + interval '48 hour') - greatest(ie.starttime, i.intime))) / 3600.0
         WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mcg/kg/hour'
              AND ie.itemid IN (221744, 225942) AND ie.patientweight > 0
           THEN ie.rate * ie.patientweight * extract(epoch FROM (least(ie.endtime, i.intime + interval '48 hour') - greatest(ie.starttime, i.intime))) / 3600.0
         WHEN ie.rate IS NULL AND ie.starttime >= i.intime AND ie.amount > 0
              AND ie.itemid IN (221744, 225942) AND ie.amountuom = 'mcg'
           THEN ie.amount
         WHEN ie.rate IS NULL AND ie.starttime >= i.intime AND ie.amount > 0
              AND ie.itemid IN (221744, 225942) AND ie.amountuom = 'mg'
           THEN ie.amount * 1000.0
         ELSE NULL
       END AS fentanyl_mcg,
       CASE
         WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mg/hour' AND ie.itemid = 221833
           THEN ie.rate * extract(epoch FROM (least(ie.endtime, i.intime + interval '48 hour') - greatest(ie.starttime, i.intime))) / 3600.0
         WHEN ie.rate > 0 AND ie.rateuom = 'mg/min' AND ie.itemid = 221833
           THEN ie.rate * extract(epoch FROM (least(ie.endtime, i.intime + interval '48 hour') - greatest(ie.starttime, i.intime))) / 60.0
         WHEN ie.rate IS NULL AND ie.starttime >= i.intime AND ie.amount > 0
              AND ie.itemid = 221833 AND ie.amountuom = 'mg'
           THEN ie.amount
         WHEN ie.rate IS NULL AND ie.starttime >= i.intime AND ie.amount > 0
              AND ie.itemid = 221833 AND ie.amountuom = 'mcg'
           THEN ie.amount / 1000.0
         ELSE NULL
       END AS hydromorphone_mg,
       CASE
         WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mg/hour' AND ie.itemid = 225154
           THEN ie.rate * extract(epoch FROM (least(ie.endtime, i.intime + interval '48 hour') - greatest(ie.starttime, i.intime))) / 3600.0
         WHEN ie.rate IS NULL AND ie.starttime >= i.intime AND ie.amount > 0
              AND ie.itemid = 225154 AND ie.amountuom = 'mg'
           THEN ie.amount
         ELSE NULL
       END AS morphine_mg
FROM selected_stays s
JOIN mimiciv_icu.icustays i USING (stay_id)
JOIN mimiciv_icu.inputevents ie USING (stay_id)
LEFT JOIN mimiciv_icu.d_items di ON di.itemid = ie.itemid
WHERE ie.itemid IN (221744, 225942, 221833, 225154)
  AND ie.starttime < i.intime + interval '48 hour'
  AND (ie.endtime > i.intime OR (ie.rate IS NULL AND ie.starttime >= i.intime))
  AND coalesce(ie.statusdescription, '') NOT IN ('Rewritten', 'Cancelled')
ORDER BY ie.stay_id, ie.starttime, ie.orderid", stay_values)

events <- as.data.table(dbGetQuery(con, sql))
events[, drug := fcase(
  itemid %in% c(221744L, 225942L), "fentanyl",
  itemid == 221833L, "hydromorphone",
  itemid == 225154L, "morphine",
  default = "other"
)]
events[, ivme_event := fifelse(!is.na(fentanyl_mcg), fentanyl_mcg * 0.1, 0) +
                         fifelse(!is.na(hydromorphone_mg), hydromorphone_mg * 5, 0) +
                         fifelse(!is.na(morphine_mg), morphine_mg, 0)]
events[, supported := !(is.na(fentanyl_mcg) & is.na(hydromorphone_mg) & is.na(morphine_mg))]

calc <- events[, .(
  recognized_events = .N,
  unsupported_events = sum(!supported),
  fentanyl_mcg_calc = sum(fentanyl_mcg, na.rm = TRUE),
  hydromorphone_mg_calc = sum(hydromorphone_mg, na.rm = TRUE),
  morphine_mg_calc = sum(morphine_mg, na.rm = TRUE),
  ivme_calc = sum(ivme_event, na.rm = TRUE)
), by = stay_id]
calc <- merge(final[, .(
  stay_id, recognized_opioid_event_n, unsupported_event_n,
  fentanyl_mcg_48h, hydromorphone_mg_48h, morphine_mg_48h, ivme_nci_48h
)], calc, by = "stay_id", all.x = TRUE)
for (v in c("recognized_events", "unsupported_events", "fentanyl_mcg_calc",
            "hydromorphone_mg_calc", "morphine_mg_calc", "ivme_calc")) {
  set(calc, which(is.na(calc[[v]])), v, 0)
}

validation <- data.table(
  check = c(
    "recognized event count",
    "unsupported event count",
    "fentanyl total rounded to 0.001",
    "hydromorphone total rounded to 0.001",
    "morphine total rounded to 0.001",
    "IVME total rounded to 0.001"
  ),
  mismatched_stays = c(
    sum(calc$recognized_opioid_event_n != calc$recognized_events),
    sum(calc$unsupported_event_n != calc$unsupported_events),
    sum(calc$fentanyl_mcg_48h != round(calc$fentanyl_mcg_calc, 3)),
    sum(calc$hydromorphone_mg_48h != round(calc$hydromorphone_mg_calc, 3)),
    sum(calc$morphine_mg_48h != round(calc$morphine_mg_calc, 3)),
    sum(calc$ivme_nci_48h != round(calc$ivme_calc, 3))
  )
)
fwrite(validation, file.path(out_dir, "raw_recalculation_validation.csv"))
if (any(validation$mismatched_stays != 0L)) {
  stop("Raw inputevents recalculation did not match the frozen cohort extract.")
}

positive_unrounded <- calc[ivme_calc > 0, ivme_calc]
positive <- round(positive_unrounded, 3)
quantile_types <- rbindlist(lapply(c("database_double", "rounded_0.001"), function(scale_name) {
  x <- if (scale_name == "database_double") positive_unrounded else positive
  rbindlist(lapply(1:9, function(tp) data.table(
    value_scale = scale_name,
    quantile_type = tp,
    p25 = unname(quantile(x, 0.25, type = tp)),
    p50 = unname(quantile(x, 0.50, type = tp)),
    p75 = unname(quantile(x, 0.75, type = tp))
  )))
}))
fwrite(quantile_types, file.path(out_dir, "empirical_quantiles_types_1_to_9.csv"))

raw_nominal_check <- rbindlist(lapply(c(10, 25), function(cutpoint) {
  x <- positive_unrounded[round(positive_unrounded, 3) == cutpoint]
  data.table(
    cutpoint_mg = cutpoint,
    n_rounding_to_cutpoint = length(x),
    minimum_unrounded_mg = min(x),
    maximum_unrounded_mg = max(x),
    maximum_absolute_deviation_mg = max(abs(x - cutpoint))
  )
}))
fwrite(raw_nominal_check, file.path(out_dir, "raw_double_precision_near_10_25.csv"))

cdf_crossings <- data.table(
  cutpoint_mg = c(10, 25),
  n_positive = length(positive),
  n_below = c(sum(positive < 10), sum(positive < 25)),
  n_equal = c(sum(positive == 10), sum(positive == 25)),
  n_at_or_below = c(sum(positive <= 10), sum(positive <= 25)),
  pct_below = 100 * c(mean(positive < 10), mean(positive < 25)),
  pct_equal = 100 * c(mean(positive == 10), mean(positive == 25)),
  pct_at_or_below = 100 * c(mean(positive <= 10), mean(positive <= 25))
)
fwrite(cdf_crossings, file.path(out_dir, "cdf_crossings_10_25.csv"))

calc[, drug_pattern := paste0(
  fifelse(round(fentanyl_mcg_calc, 6) > 0, "F", ""),
  fifelse(round(hydromorphone_mg_calc, 6) > 0, "H", ""),
  fifelse(round(morphine_mg_calc, 6) > 0, "M", "")
)]
calc[drug_pattern == "", drug_pattern := "none"]
exact <- calc[round(ivme_calc, 3) %in% c(10, 25)]
exact[, total_ivme_mg := round(ivme_calc, 3)]
exact_compositions <- exact[, .N, by = .(
  total_ivme_mg,
  fentanyl_mcg = round(fentanyl_mcg_calc, 3),
  hydromorphone_mg = round(hydromorphone_mg_calc, 3),
  morphine_mg = round(morphine_mg_calc, 3),
  drug_pattern,
  recognized_events
)][order(total_ivme_mg, -N, fentanyl_mcg, hydromorphone_mg, morphine_mg)]
fwrite(exact_compositions, file.path(out_dir, "exact_10_25_aggregate_compositions.csv"))

exact_patterns <- exact[, .N, by = .(total_ivme_mg, drug_pattern)][order(total_ivme_mg, -N)]
exact_patterns[, pct_within_cutpoint := 100 * N / sum(N), by = total_ivme_mg]
fwrite(exact_patterns, file.path(out_dir, "exact_10_25_drug_patterns.csv"))

bolus_grid <- events[is.na(rate) & supported & !is.na(amount) & amount > 0,
  .N, by = .(drug, itemid, item_label, amount, amountuom)][order(drug, -N)]
bolus_grid[, rank_within_drug := seq_len(.N), by = drug]
fwrite(bolus_grid[rank_within_drug <= 25], file.path(out_dir, "raw_bolus_amount_grid_top25.csv"))

infusion_grid <- events[!is.na(rate) & supported & rate > 0,
  .N, by = .(drug, itemid, item_label, rate, rateuom)][order(drug, -N)]
infusion_grid[, rank_within_drug := seq_len(.N), by = drug]
fwrite(infusion_grid[rank_within_drug <= 25], file.path(out_dir, "raw_infusion_rate_grid_top25.csv"))

event_ivme_grid <- events[supported == TRUE,
  .N, by = .(drug, ivme_event_mg = round(ivme_event, 3), administration = fifelse(is.na(rate), "bolus", "infusion"))][order(drug, administration, -N)]
event_ivme_grid[, rank_within_stratum := seq_len(.N), by = .(drug, administration)]
fwrite(event_ivme_grid[rank_within_stratum <= 25], file.path(out_dir, "standardized_event_ivme_grid_top25.csv"))

summary_out <- data.table(
  metric = c(
    "final cohort patients", "raw opioid inputevents", "supported raw events",
    "positive-dose patients", "exactly 10 mg patients", "exactly 25 mg patients",
    "all quantile types return P25=10", "all quantile types return P75=25"
  ),
  value = c(
    nrow(final), nrow(events), sum(events$supported), length(positive),
    sum(positive == 10), sum(positive == 25),
    all(quantile_types[value_scale == "rounded_0.001", p25] == 10),
    all(quantile_types[value_scale == "rounded_0.001", p75] == 25)
  )
)
fwrite(summary_out, file.path(out_dir, "raw_source_audit_summary.csv"))

print(validation)
print(quantile_types)
print(raw_nominal_check)
print(cdf_crossings)
print(exact_patterns)
cat("RAW SOURCE DOSE AUDIT COMPLETE\n")
