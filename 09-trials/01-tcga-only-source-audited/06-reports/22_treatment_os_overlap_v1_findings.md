# TCGA-BRCA Treatment-OS Overlap V1 Findings

This note records human-reviewed findings from the TCGA-BRCA treatment-OS overlap audit v1 workflow.

Important reminders:

- this document reviews overlap between saved treatment biotab records and the saved OS endpoint v1 cohort only
- this document does not normalize drug names or freeze treatment arms
- this document does not perform survival modeling, treatment-effect estimation, or causal analysis
- this document does not replace source biotabs, prior endpoint outputs, prior baseline model-input outputs, saved audit TSVs, or run logs
- this document remains TCGA-only and source-audited

## Reviewed Treatment-OS Overlap V1 Run ID

- Treatment-OS overlap latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_treatment_os_overlap_v1_latest.json`
- Reviewed treatment-OS overlap v1 run id: `20260414T172924Z`
- Source OS endpoint v1 run id: `20260414T162349Z`
- Source clinical biotab parse run id: `20260412T010932Z`
- Source baseline model-input v1 run id: `20260414T013108Z`
- Source minimal cohort v1 build id: `20260413T202134Z`
- Source clinical supplement/source run id: `20260412T000556Z`
- Review date: `2026-04-14`
- Reviewed by: `[fill in name or initials]`

## Treatment-OS Overlap Inputs Used

- Treatment-OS overlap latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_treatment_os_overlap_v1_latest.json`
- OS endpoint v1 latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_os_endpoint_v1_latest.json`
- Clinical biotabs latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_biotabs_latest.json`
- Baseline model-input v1 latest pointer: `01-data/audit/tcga-brca/model-input/tcga_brca_baseline_model_input_v1_latest.json`
- Processed OS endpoint TSV: `01-data/processed/tcga-brca/endpoint-prep/os_endpoint_v1_runs/20260414T162349Z/os_endpoint_v1.tsv`
- Processed clinical drug TSV: `01-data/processed/tcga-brca/clinical/biotab_runs/20260412T010932Z/clinical_drug.tsv`
- Processed clinical radiation TSV: `01-data/processed/tcga-brca/clinical/biotab_runs/20260412T010932Z/clinical_radiation.tsv`
- Processed clinical patient TSV: `01-data/processed/tcga-brca/clinical/biotab_runs/20260412T010932Z/clinical_patient.tsv`
- Audit overlap summary TSV: `01-data/audit/tcga-brca/treatment-prep/treatment_os_overlap_v1_runs/20260414T172924Z/treatment_os_overlap_summary.tsv`
- Audit overlap counts TSV: `01-data/audit/tcga-brca/treatment-prep/treatment_os_overlap_v1_runs/20260414T172924Z/treatment_os_overlap_counts.tsv`
- Audit drug therapy-type distribution TSV: `01-data/audit/tcga-brca/treatment-prep/treatment_os_overlap_v1_runs/20260414T172924Z/drug_therapy_type_distribution.tsv`
- Audit drug therapy-type patient summary TSV: `01-data/audit/tcga-brca/treatment-prep/treatment_os_overlap_v1_runs/20260414T172924Z/drug_therapy_type_patient_summary.tsv`
- Audit summary TSV: `01-data/audit/tcga-brca/treatment-prep/treatment_os_overlap_v1_runs/20260414T172924Z/treatment_os_overlap_audit_summary.tsv`
- Audit run log: `01-data/audit/tcga-brca/treatment-prep/treatment_os_overlap_v1_runs/20260414T172924Z/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
- `107_treatment_os_overlap_summary.tsv`
- `108_treatment_os_overlap_counts.tsv`
- `109_drug_therapy_type_distribution.tsv`
- `110_drug_therapy_type_patient_summary.tsv`
- `111_treatment_os_overlap_audit_summary.tsv`

## Audit Rules Used

- Unit of analysis: patient/case.
- Treatment scope: `clinical_drug.tsv` and `clinical_radiation.tsv` only.
- OS cohort source: saved OS endpoint v1 table only.
- Confounder snapshot source: `clinical_patient.tsv` only.
- Missing-like therapy-type normalization uses `value.strip().lower()` and treats the following as missing-like: `""`, `[discrepancy]`, `[not applicable]`, `[not available]`, `[not evaluated]`, `[unknown]`, `n/a`, `na`, `nan`, `none`, `null`.
- This workflow does not normalize `pharmaceutical_therapy_drug_name`, does not aggregate treatment arms, and does not perform modeling.

## Overlap Coverage Counts

- Final OS cohort size reviewed here: `1097` patients.
- OS event count in the reviewed cohort: `152`.
- OS censor count in the reviewed cohort: `945`.
- Patients with any drug row: `780` (`71.1%`).
- Patients with no drug rows: `317`.
- Patients with any radiation row: `528`.
- Patients with no radiation rows: `569`.
- Total clinical drug rows overlapping the OS cohort: `2406`.
- Total clinical radiation rows overlapping the OS cohort: `618`.
- Every OS cohort patient appears once in the overlap summary, and no duplicate overlap-summary barcodes were reported.

## Therapy-Type Structure

- Patients with a single non-missing therapy-type value: `423`.
- Patients with mixed non-missing therapy-type values: `355`.
- Patients with drug rows but only missing-like therapy-type values: `2`.
- Patients with an identifiable dominant therapy type: `729` total (`66.5%` of the full OS cohort; `93.5%` of treated patients).
- Dominant therapy groups observed:
- `Chemotherapy`: `526`
- `Hormone Therapy`: `199`
- `Chemotherapy|Hormone Therapy`: `2`
- `Ancillary`: `1`
- `Immunotherapy`: `1`
- Treated patients without a usable dominant therapy label:
- `mixed_no_dominant_therapy_type`: `49`
- `missing_therapy_type_despite_drug_rows`: `2`
- Patients with no drug rows at all: `317`

## Row-Count Structure

- Median drug-row count across all OS patients: `2.0`.
- Median drug-row count among patients with any drug row: `3.0`.
- Maximum drug-row count observed for one patient: `23`.
- Maximum radiation-row count observed for one patient: `5`.
- Drug-row count distribution shows a large structural no-drug bucket (`317` patients) plus a long upper tail of repeated treatment rows.

## Structural Cautions

- OS event rate differs sharply by drug-record presence: `59/780 (7.6%)` for `has_drug_row = yes` versus `93/317 (29.3%)` for `has_drug_row = no`.
- OS event rate also differs by radiation-record presence: `44/528 (8.3%)` for `has_radiation_row = yes` versus `108/569 (19.0%)` for `has_radiation_row = no`.
- These differences should be interpreted as structural overlap observations only. This review does not establish treatment exposure, treatment effect, confounding control, or causal direction.
- A substantial mixed-therapy subgroup remains (`355` patients), so any next-step grouping must keep explicit provenance for mixed versus single versus unknown treatment structure.

## Whether Treatment-OS Overlap V1 Is Ready For The Next Step

This section stays explicit that readiness means first-pass patient-level treatment aggregation only.

- Yes: the saved overlap audit is ready for patient-level treatment aggregation and coarse exploratory grouping.
- The saved audit summary reports `feasibility_interpretation == ready_for_coarse_exploratory_grouping`.
- The supporting rationale is direct: `71.1%` of the OS cohort has at least one drug row, and `729/780` treated patients (`93.5%`) already have an identifiable dominant therapy type.
- A `treatment_unknown` bucket is still required for the `317` patients with no drug rows.
- Mixed cases must stay explicitly flagged during aggregation, especially the `49` treated patients with mixed therapy-type values but no single dominant label.
- This workflow is not yet ready for normalized treatment-arm freezing or modeling.

## What Still Blocks Stronger Treatment Claims

- `pharmaceutical_therapy_drug_name` has not yet been audited or normalized, so arm-level treatment grouping is still blocked.
- `355` treated patients carry mixed therapy-type structure, and `49` of those do not resolve to a single dominant therapy type from the saved audit alone.
- The overlap audit is join-ready for the next aggregation step, but it is not evidence of balanced treatment groups, interpretable exposure windows, or causal comparability.

## Validation Checks Completed

- [x] Treatment-OS overlap latest pointer resolves to the reviewed run.
- [x] Reviewed run references the intended OS endpoint v1 run, clinical biotab parse run, baseline model-input v1 run, minimal cohort v1 build, and source clinical supplement run.
- [x] `treatment_os_overlap_summary.tsv`, `treatment_os_overlap_counts.tsv`, `drug_therapy_type_distribution.tsv`, `drug_therapy_type_patient_summary.tsv`, `treatment_os_overlap_audit_summary.tsv`, and `run_log.json` all exist and are non-empty.
- [x] `run_log.json` reports `status == completed`.
- [x] `run_log.json` reports `validation.passed == true`.
- [x] `run_log.json` reports `required_upstream_pointers_found == true`.
- [x] `run_log.json` reports `required_source_tables_found == true`.
- [x] `run_log.json` reports `os_endpoint_run_log_completed == true`.
- [x] `run_log.json` reports `clinical_biotabs_run_log_completed == true`.
- [x] `run_log.json` reports `baseline_model_input_run_log_completed == true`.
- [x] `run_log.json` reports `os_cohort_row_count_positive == true`.
- [x] `run_log.json` reports `drug_row_count_positive == true`.
- [x] `run_log.json` reports `radiation_row_count_positive == true`.
- [x] `run_log.json` reports `overlap_summary_row_count_matches_os_cohort == true`.
- [x] `run_log.json` reports `all_os_patients_in_overlap_summary == true`.
- [x] `run_log.json` reports `no_overlap_summary_duplicate_barcodes == true`.
- [x] `run_log.json` reports `patient_counts_reconcile == true`.
- [x] `run_log.json` reports `no_prior_run_overwrite == true`.
- [x] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [x] Review tables `107` through `111` were regenerated from the saved overlap outputs on disk only.
- [x] `111_treatment_os_overlap_audit_summary.tsv` reports `feasibility_interpretation == ready_for_coarse_exploratory_grouping`.
- [x] `111_treatment_os_overlap_audit_summary.tsv` reports `coarse_grouping_next_step == patient_level_treatment_aggregation`.
- [x] `111_treatment_os_overlap_audit_summary.tsv` reports `remaining_blocker == drug_name_not_normalized`.
- [x] This note remains an overlap-audit artifact only and does not present normalized treatment arms, treatment effects, or a modeling result.

## Current Interpretation

This note remains a treatment-OS overlap review artifact only. It answers whether the saved treatment biotab record structure overlaps the saved OS endpoint cohort well enough to justify the next audit step. The answer is yes for patient-level treatment aggregation: the overlap table is cohort-complete, treated-patient coverage is substantial, and most treated patients already have an identifiable dominant therapy-type label. The answer is still no for normalized treatment-arm freezing or modeling: drug names remain unnormalized, mixed therapy-type cases still need patient-level resolution, and the observed event-rate differences by record presence are structural rather than causal findings.
