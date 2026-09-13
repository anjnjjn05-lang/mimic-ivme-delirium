
-- Revised cohort v3, 2026-09-05. Exports ALL landmark patients, including unassessed.
-- This amendment follows inspection of earlier results; it is not prospective preregistration.
-- Procedure dictionary: v1.0-Zhong-2026-09-04.
-- Dose dictionary: v1.0-NCI-parenteral-2026-09-04.
-- Estimand: association of 0-48 h recorded ICU IVME group with incident CAM-ICU positivity after the 48 h landmark.
-- This exports enrollment, exposure, outcome and basic baseline fields. Add covariates in a separate,
-- versioned script; do not overwrite this cohort definition.

-- SECTION 1: Cardiac procedure dictionary, first eligible postoperative ICU stay,
-- adult/48-hour landmark eligibility, and pre-landmark exclusions.
WITH procedure_dictionary AS (
    SELECT 'CABG'::text AS procedure_group, 9 AS icd_version, '3610'::text AS code_lo, '3619'::text AS code_hi
    UNION ALL SELECT 'valve', 9, '3599', '3599' UNION ALL SELECT 'valve', 9, '3500', '3504'
    UNION ALL SELECT 'valve', 9, '3510', '3514' UNION ALL SELECT 'valve', 9, '3520', '3528'
    UNION ALL SELECT 'septal_defect_repair', 9, '3550', '3555' UNION ALL SELECT 'septal_defect_repair', 9, '3560', '3560'
    UNION ALL SELECT 'septal_defect_repair', 9, '3570', '3570' UNION ALL SELECT 'septal_defect_repair', 9, '3571', '3571'
    UNION ALL SELECT 'septal_defect_repair', 9, '3598', '3598' UNION ALL SELECT 'aortic_replacement', 9, '3845', '3845'
),
icd10_prefix_dictionary AS (
    SELECT 'CABG'::text AS procedure_group, '0210'::text AS prefix
    UNION ALL SELECT 'CABG', '0211' UNION ALL SELECT 'CABG', '0212' UNION ALL SELECT 'CABG', '0213'
    UNION ALL SELECT 'valve', '02QF' UNION ALL SELECT 'valve', '02QG' UNION ALL SELECT 'valve', '02QH' UNION ALL SELECT 'valve', '02QJ'
    UNION ALL SELECT 'valve', '02RF' UNION ALL SELECT 'valve', '02RG' UNION ALL SELECT 'valve', '02RH' UNION ALL SELECT 'valve', '02RJ'
    UNION ALL SELECT 'aortic_replacement', '02RW' UNION ALL SELECT 'aortic_replacement', '02RX'
),
matched_procedures AS MATERIALIZED (
    SELECT p.hadm_id, p.chartdate, p.icd_version, p.icd_code, d.procedure_group
    FROM mimiciv_hosp.procedures_icd p JOIN procedure_dictionary d
      ON p.icd_version = d.icd_version AND p.icd_code BETWEEN d.code_lo AND d.code_hi
    UNION ALL
    SELECT p.hadm_id, p.chartdate, p.icd_version, p.icd_code, d.procedure_group
    FROM mimiciv_hosp.procedures_icd p JOIN icd10_prefix_dictionary d
      ON p.icd_version = 10 AND left(p.icd_code, 4) = d.prefix
),
candidate_pairs AS MATERIALIZED (
    SELECT DISTINCT i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime,
           a.admittime, a.dischtime, a.deathtime, a.admission_type, a.race, a.insurance
    FROM mimiciv_icu.icustays i
    JOIN mimiciv_hosp.admissions a ON a.hadm_id = i.hadm_id
    JOIN matched_procedures m ON m.hadm_id = i.hadm_id
                            AND m.chartdate BETWEEN i.intime::date - 1 AND i.intime::date
),
first_postop AS MATERIALIZED (
    SELECT * FROM (
        SELECT c.*, row_number() OVER (PARTITION BY c.subject_id ORDER BY c.intime, c.stay_id) AS postop_icu_rank
        FROM candidate_pairs c
    ) x WHERE postop_icu_rank = 1
),
prior_dementia AS MATERIALIZED (
    SELECT DISTINCT c.stay_id
    FROM first_postop c
    JOIN mimiciv_hosp.admissions a_prev ON a_prev.subject_id = c.subject_id AND a_prev.dischtime < c.admittime
    JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = a_prev.hadm_id
    WHERE (d.icd_version = 10 AND (d.icd_code LIKE 'F00%' OR d.icd_code LIKE 'F01%' OR d.icd_code LIKE 'F02%'
                                   OR d.icd_code LIKE 'F03%' OR d.icd_code LIKE 'G30%'))
       OR (d.icd_version = 9 AND d.icd_code LIKE '290%')
),
landmark AS MATERIALIZED (
    SELECT c.*, p.gender, p.anchor_year_group,
           EXISTS (SELECT 1 FROM mimiciv_hosp.admissions ap WHERE ap.subject_id=c.subject_id AND ap.dischtime<c.admittime)::int AS prior_admission_observed,
           (p.anchor_age + extract(year FROM c.admittime)::int - p.anchor_year)::int AS age_at_admission
    FROM first_postop c
    JOIN mimiciv_hosp.patients p ON p.subject_id = c.subject_id
    LEFT JOIN prior_dementia pd ON pd.stay_id = c.stay_id
    WHERE pd.stay_id IS NULL
      AND c.outtime >= c.intime + interval '48 hour'
      AND c.dischtime > c.intime + interval '48 hour'
      AND (c.deathtime IS NULL OR c.deathtime > c.intime + interval '48 hour')
      AND (p.anchor_age + extract(year FROM c.admittime)::int - p.anchor_year) >= 18
),
procedure_groups AS (
    SELECT l.stay_id, string_agg(DISTINCT m.procedure_group, '+' ORDER BY m.procedure_group) AS procedure_groups
    FROM landmark l JOIN matched_procedures m ON m.hadm_id = l.hadm_id
                                             AND m.chartdate BETWEEN l.intime::date - 1 AND l.intime::date
    GROUP BY l.stay_id
),
-- SECTION 2: Post-landmark CAM-ICU outcome and follow-up censoring.
followup_cam AS (
    SELECT l.stay_id,
           COUNT(*) FILTER (WHERE lower(btrim(coalesce(ce.value, ''))) IN ('positive', 'negative')) AS valid_cam_assessment_n,
           BOOL_OR(lower(btrim(coalesce(ce.value, ''))) = 'positive') AS incident_cam_positive
    FROM landmark l
    LEFT JOIN mimiciv_icu.chartevents ce ON ce.stay_id = l.stay_id AND ce.itemid = 228332
                                        AND ce.charttime >= l.intime + interval '48 hour'
                                        AND ce.charttime < least(l.intime + interval '7 day', l.outtime, l.dischtime, coalesce(l.deathtime,l.outtime))
    GROUP BY l.stay_id
),
early_cam AS (
    SELECT l.stay_id, BOOL_OR(lower(btrim(coalesce(ce.value, ''))) = 'positive') AS early_cam_positive,
           COUNT(*) FILTER (WHERE lower(btrim(coalesce(ce.value,'')))='negative') AS early_negative_n
    FROM landmark l
    LEFT JOIN mimiciv_icu.chartevents ce ON ce.stay_id = l.stay_id AND ce.itemid = 228332
                                        AND ce.charttime >= l.intime
                                        AND ce.charttime < l.intime + interval '48 hour'
    GROUP BY l.stay_id
),
-- SECTION 3: Recorded 0-48-hour intravenous opioid exposure, window clipping,
-- unit validation, rate-by-duration conversion, and IVME components.
standardized_events AS (
    SELECT l.stay_id, ie.itemid,
           CASE
               WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mcg/hour'
                    AND ie.itemid IN (221744, 225942)
                    THEN ie.rate * extract(epoch FROM (least(ie.endtime, l.intime + interval '48 hour') - greatest(ie.starttime, l.intime))) / 3600.0
               WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mcg/kg/hour'
                    AND ie.itemid IN (221744,225942) AND ie.patientweight > 0
                    THEN ie.rate * ie.patientweight * extract(epoch FROM (least(ie.endtime,l.intime+interval '48 hour')-greatest(ie.starttime,l.intime))) / 3600.0
               WHEN ie.rate IS NULL AND ie.starttime >= l.intime AND ie.amount > 0 AND ie.itemid IN (221744,225942) AND ie.amountuom = 'mcg' THEN ie.amount
               WHEN ie.rate IS NULL AND ie.starttime >= l.intime AND ie.amount > 0 AND ie.itemid IN (221744,225942) AND ie.amountuom = 'mg' THEN ie.amount * 1000.0
               ELSE NULL
           END AS fentanyl_mcg,
           CASE
               WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mg/hour' AND ie.itemid = 221833
                    THEN ie.rate * extract(epoch FROM (least(ie.endtime, l.intime + interval '48 hour') - greatest(ie.starttime, l.intime))) / 3600.0
               WHEN ie.rate > 0 AND ie.rateuom='mg/min' AND ie.itemid=221833
                    THEN ie.rate * extract(epoch FROM (least(ie.endtime,l.intime+interval '48 hour')-greatest(ie.starttime,l.intime))) / 60.0
               WHEN ie.rate IS NULL AND ie.starttime >= l.intime AND ie.amount > 0 AND ie.itemid = 221833 AND ie.amountuom = 'mg' THEN ie.amount
               WHEN ie.rate IS NULL AND ie.starttime >= l.intime AND ie.amount > 0 AND ie.itemid = 221833 AND ie.amountuom = 'mcg' THEN ie.amount/1000.0
               ELSE NULL
           END AS hydromorphone_mg,
           CASE
               WHEN ie.rate IS NOT NULL AND ie.rate > 0 AND ie.rateuom = 'mg/hour' AND ie.itemid = 225154
                    THEN ie.rate * extract(epoch FROM (least(ie.endtime, l.intime + interval '48 hour') - greatest(ie.starttime, l.intime))) / 3600.0
               WHEN ie.rate IS NULL AND ie.starttime >= l.intime AND ie.amount > 0 AND ie.itemid = 225154 AND ie.amountuom = 'mg' THEN ie.amount
               ELSE NULL
           END AS morphine_mg
    FROM landmark l
    JOIN mimiciv_icu.inputevents ie ON ie.stay_id = l.stay_id
                                  AND ie.itemid IN (221744, 225942, 221833, 225154)
                                  AND ie.starttime < l.intime + interval '48 hour'
                                  AND (ie.endtime > l.intime OR (ie.rate IS NULL AND ie.starttime >= l.intime))
                                  AND coalesce(ie.statusdescription,'') NOT IN ('Rewritten','Cancelled')
),
exposure AS (
    SELECT stay_id,
           COUNT(*) AS recognized_opioid_event_n,
           COUNT(*) FILTER (WHERE fentanyl_mcg IS NULL AND hydromorphone_mg IS NULL AND morphine_mg IS NULL) AS unsupported_event_n,
           coalesce(sum(fentanyl_mcg), 0) AS fentanyl_mcg_48h,
           coalesce(sum(hydromorphone_mg), 0) AS hydromorphone_mg_48h,
           coalesce(sum(morphine_mg), 0) AS morphine_mg_48h
    FROM standardized_events
    GROUP BY stay_id
),
-- SECTION 4: Merge cohort, procedure class, outcome, exposure, and basic
-- admission covariates into the patient-level result returned to R.
analysis_eligible AS (
    SELECT l.subject_id, l.hadm_id, l.stay_id, l.intime, l.admittime, l.outtime, l.dischtime, l.deathtime, l.gender, l.age_at_admission,
           l.anchor_year_group, l.prior_admission_observed, ec.early_negative_n,
           coalesce(ec.early_cam_positive,false)::int AS early_cam_positive,
           l.admission_type, l.race, l.insurance, pg.procedure_groups,
           c.valid_cam_assessment_n, c.incident_cam_positive,
           coalesce(e.recognized_opioid_event_n, 0) AS recognized_opioid_event_n,
           coalesce(e.unsupported_event_n, 0) AS unsupported_event_n,
           coalesce(e.fentanyl_mcg_48h, 0) AS fentanyl_mcg_48h,
           coalesce(e.hydromorphone_mg_48h, 0) AS hydromorphone_mg_48h,
           coalesce(e.morphine_mg_48h, 0) AS morphine_mg_48h,
           (coalesce(e.fentanyl_mcg_48h, 0) * 0.1 + coalesce(e.hydromorphone_mg_48h, 0) * 5.0 + coalesce(e.morphine_mg_48h, 0)) AS ivme_nci_48h
    FROM landmark l
    JOIN procedure_groups pg USING (stay_id)
    JOIN followup_cam c USING (stay_id)
    JOIN early_cam ec USING (stay_id)
    LEFT JOIN exposure e USING (stay_id)
),
final_cohort AS (
    -- The 10/25 mg boundaries originated as the empirical positive-dose P25/P75 values
    -- in the earlier v1 cohort (n=1,559), were fixed in v2, and were retained unchanged
    -- through the v3 cohort revision. In the final v3 cohort they still coincide with the
    -- empirical positive-dose P25/P75 because both percentile ranks fall within tied-dose
    -- blocks at 10 and 25 mg. The values are not dynamically recomputed in this SQL and
    -- should not be interpreted as clinically validated prescribing thresholds.
    SELECT a.*,
           CASE
               WHEN round(ivme_nci_48h::numeric, 3) = 0 THEN 'G1'
               WHEN round(ivme_nci_48h::numeric, 3) <= 10.0 THEN 'G2'
               WHEN round(ivme_nci_48h::numeric, 3) <= 25.0 THEN 'G3'
               ELSE 'G4'
           END AS ivme_group,
           CASE WHEN round(ivme_nci_48h::numeric, 3) > 25.0 THEN 1 ELSE 0 END AS high_ivme_nci_q4
    FROM analysis_eligible a
)
SELECT 'v3-2026-09-05'::text AS cohort_version,
       subject_id, hadm_id, stay_id, intime, admittime, outtime, dischtime, deathtime, gender, age_at_admission,
       anchor_year_group, prior_admission_observed, early_negative_n, early_cam_positive, unsupported_event_n,
       admission_type, race, insurance, procedure_groups,
       valid_cam_assessment_n, incident_cam_positive::int AS incident_cam_positive,
       recognized_opioid_event_n,
       round(fentanyl_mcg_48h::numeric, 3) AS fentanyl_mcg_48h,
       round(hydromorphone_mg_48h::numeric, 3) AS hydromorphone_mg_48h,
       round(morphine_mg_48h::numeric, 3) AS morphine_mg_48h,
       round(ivme_nci_48h::numeric, 3) AS ivme_nci_48h,
       ivme_group,
       high_ivme_nci_q4,
       -- Legacy column names retained for backward compatibility with frozen derived files.
       10.0::numeric AS nci_positive_user_p25_ivme,
       25.0::numeric AS nci_positive_user_p75_ivme
FROM final_cohort
ORDER BY stay_id;
