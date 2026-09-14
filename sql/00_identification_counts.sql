-- Counts preceding the v3 landmark cohort. Uses the frozen literature-reproduction dictionary.
WITH mp AS (
 SELECT hadm_id,chartdate FROM mimiciv_hosp.procedures_icd
 WHERE (icd_version=9 AND ((icd_code BETWEEN '3610' AND '3619') OR (icd_code BETWEEN '3500' AND '3504') OR
  (icd_code BETWEEN '3510' AND '3514') OR (icd_code BETWEEN '3520' AND '3528') OR
  icd_code IN ('3599','3560','3570','3571','3598','3845') OR icd_code BETWEEN '3550' AND '3555'))
 OR (icd_version=10 AND left(icd_code,4) IN ('0210','0211','0212','0213','02QF','02QG','02QH','02QJ','02RF','02RG','02RH','02RJ','02RW','02RX'))
), cp AS (
 SELECT DISTINCT i.subject_id,i.hadm_id,i.stay_id FROM mimiciv_icu.icustays i JOIN mp m ON m.hadm_id=i.hadm_id
 AND m.chartdate BETWEEN i.intime::date-1 AND i.intime::date
)
SELECT count(*) AS procedure_icu_pairs,count(DISTINCT subject_id) AS unique_patients FROM cp;
