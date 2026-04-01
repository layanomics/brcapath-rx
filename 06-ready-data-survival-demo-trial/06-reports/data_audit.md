# Data Audit

Generated: 2026-04-01T21:09:26+00:00

## Cohort Summary

   cohort  patients_total  os_complete_n  os_event_n primary_subtype_field  primary_subtype_non_missing_n os_time_unit                                                             notes
TCGA-BRCA            1076           1076         150      PAM50Call_RNAseq                            826         days Subtype is augmented from the UCSC Xena TCGA hub clinical matrix.
 METABRIC            2509           1981        1144       CLAUDIN_SUBTYPE                           1980       months     Subtype comes directly from cBioPortal patient clinical data.

## Survival Endpoint Summary

   cohort endpoint     time_column time_unit  complete_n  event_n  median_time_observed
TCGA-BRCA       OS    os_time_days      days        1076      150            860.000000
 METABRIC       OS  os_time_months    months        1981     1144            116.466667
 METABRIC      RFS rfs_time_months    months        2380      964            100.050000

## Subtype Distribution Summary

   cohort subtype_label  patients_n  pct_of_cohort subtype_field_used
TCGA-BRCA          LumA         415          38.57   PAM50Call_RNAseq
TCGA-BRCA       Missing         250          23.23   PAM50Call_RNAseq
TCGA-BRCA          LumB         189          17.57   PAM50Call_RNAseq
TCGA-BRCA         Basal         138          12.83   PAM50Call_RNAseq
TCGA-BRCA          Her2          62           5.76   PAM50Call_RNAseq
TCGA-BRCA        Normal          22           2.04   PAM50Call_RNAseq
 METABRIC          LumA         700          27.90    CLAUDIN_SUBTYPE
 METABRIC       Missing         529          21.08    CLAUDIN_SUBTYPE
 METABRIC          LumB         475          18.93    CLAUDIN_SUBTYPE
 METABRIC          Her2         224           8.93    CLAUDIN_SUBTYPE
 METABRIC   claudin-low         218           8.69    CLAUDIN_SUBTYPE
 METABRIC         Basal         209           8.33    CLAUDIN_SUBTYPE
 METABRIC        Normal         148           5.90    CLAUDIN_SUBTYPE
 METABRIC            NC           6           0.24    CLAUDIN_SUBTYPE

## Treatment Availability Summary

   cohort                                field             value  patients_n  pct_of_cohort
TCGA-BRCA                  prior_treatment_raw                No        1062          98.70
TCGA-BRCA                  prior_treatment_raw               Yes          12           1.12
TCGA-BRCA                  prior_treatment_raw      Not Reported           2           0.19
TCGA-BRCA history_of_neoadjuvant_treatment_raw                No        1063          98.79
TCGA-BRCA history_of_neoadjuvant_treatment_raw               Yes          12           1.12
TCGA-BRCA history_of_neoadjuvant_treatment_raw           Missing           1           0.09
TCGA-BRCA tcga_pharmaceutical_therapy_observed               yes         843          78.35
TCGA-BRCA tcga_pharmaceutical_therapy_observed                no         143          13.29
TCGA-BRCA tcga_pharmaceutical_therapy_observed      not reported          90           8.36
TCGA-BRCA      tcga_radiation_therapy_observed               yes         566          52.60
TCGA-BRCA      tcga_radiation_therapy_observed                no         427          39.68
TCGA-BRCA      tcga_radiation_therapy_observed      not reported          83           7.71
 METABRIC            chemotherapy_observed_raw                NO        1568          62.50
 METABRIC            chemotherapy_observed_raw           Missing         529          21.08
 METABRIC            chemotherapy_observed_raw               YES         412          16.42
 METABRIC         hormone_therapy_observed_raw               YES        1216          48.47
 METABRIC         hormone_therapy_observed_raw                NO         764          30.45
 METABRIC         hormone_therapy_observed_raw           Missing         529          21.08
 METABRIC           radio_therapy_observed_raw               YES        1173          46.75
 METABRIC           radio_therapy_observed_raw                NO         807          32.16
 METABRIC           radio_therapy_observed_raw           Missing         529          21.08
 METABRIC                   breast_surgery_raw        MASTECTOMY        1170          46.63
 METABRIC                   breast_surgery_raw BREAST CONSERVING         785          31.29
 METABRIC                   breast_surgery_raw           Missing         554          22.08

## Audit Notes

- TCGA patient table is collapsed to one row per patient using first non-null values.
- TCGA clinical duplicate-conflict counts after patient grouping check: age_at_earliest_diagnosis_in_years.diagnoses.xena_derived=0, ajcc_pathologic_stage.diagnoses=0, prior_treatment.diagnoses=0, treatment_type.treatments.diagnoses=0, treatment_or_therapy.treatments.diagnoses=0.
- TCGA auxiliary duplicate-conflict counts after tumor-row restriction: PAM50Call_RNAseq=0, breast_carcinoma_estrogen_receptor_status=0, breast_carcinoma_progesterone_receptor_status=0, lab_proc_her2_neu_immunohistochemistry_receptor_status=0, history_of_neoadjuvant_treatment=0.
- TCGA auxiliary subtype table is restricted to primary-tumor rows before patient collapse to avoid tumor-vs-normal PAM50 conflicts.
- METABRIC attribute count downloaded from cBioPortal: 36.
- METABRIC patient-level clinical long rows: 51528.
- METABRIC sample-level clinical long rows: 27868.
- METABRIC merged wide table rows: 2509.

## Explicit Caveats

- TCGA subtype information required a documented auxiliary UCSC Xena TCGA-hub table because the GDC-hub clinical table did not carry subtype labels.
- TCGA auxiliary subtype rows were restricted to primary-tumor samples before collapse to avoid tumor-vs-normal subtype conflicts.
- METABRIC treatment variables are observational binary descriptors and not randomized exposures.
