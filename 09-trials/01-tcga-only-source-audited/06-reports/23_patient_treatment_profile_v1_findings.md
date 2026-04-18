# TCGA-BRCA Patient Treatment Profile V1 Findings

This note records human-reviewed findings from the TCGA-BRCA patient treatment profile v1 workflow.

Important reminders:

- this document reviews patient-level treatment profile construction only
- this document does not normalize drug names or freeze final treatment arms
- this document does not perform treatment-effect estimation, causal analysis, or modeling
- this document does not replace raw treatment tables, prior overlap outputs, prior OS endpoint outputs, saved audit TSVs, or run logs
- this document remains TCGA-only and source-audited

## Reviewed Patient Treatment Profile V1 Run ID

- Patient treatment profile latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_profile_v1_latest.json`
- Reviewed patient treatment profile v1 run id: `[fill in patient treatment profile v1 run id]`
- Source treatment-OS overlap v1 run id: `[fill in treatment-OS overlap v1 run id]`
- Source OS endpoint v1 run id: `[fill in OS endpoint v1 run id]`
- Source clinical biotab parse run id: `[fill in clinical biotab parse run id]`
- Source baseline model-input v1 run id: `[fill in baseline model-input v1 run id]`
- Source cohort v1 build id: `[fill in cohort v1 build id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Patient Treatment Profile V1 Inputs Used

- Patient treatment profile latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_profile_v1_latest.json`
- Treatment-OS overlap latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_treatment_os_overlap_v1_latest.json`
- OS endpoint v1 latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_os_endpoint_v1_latest.json`
- Clinical biotabs latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_biotabs_latest.json`
- Processed patient treatment profile TSV: `01-data/processed/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1.tsv`
- Audit patient treatment profile spec TSV: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_spec.tsv`
- Audit patient treatment profile summary TSV: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_summary.tsv`
- Audit patient treatment profile conflict TSV: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_conflict_audit.tsv`
- Audit run log: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
- `112_patient_treatment_profile_v1.tsv`
- `113_patient_treatment_profile_v1_conflict_audit.tsv`
- `114_patient_treatment_profile_v1_spec.tsv`
- `115_patient_treatment_profile_v1_summary.tsv`

## Construction Rules Used

Fill from:

- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_spec.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/114_patient_treatment_profile_v1_spec.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/115_patient_treatment_profile_v1_summary.tsv`

- `[fill in unit-of-analysis statement]`
- `[fill in therapy-type single vs mixed rule summary]`
- `[fill in dominant-type provisional rule summary]`
- `[fill in regimen-context scope summary]`
- `[fill in manual-review scope summary]`
- `[fill in no-drug-record explicit interpretation]`

## Coverage And Profile Counts

Fill from:

- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/112_patient_treatment_profile_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/115_patient_treatment_profile_v1_summary.tsv`

- `[fill in OS cohort size]`
- `[fill in patients with any drug row]`
- `[fill in no-drug-record burden]`
- `[fill in patients with any radiation row]`
- `[fill in patients with no drug or radiation rows]`
- `[fill in patients requiring manual review]`

## Therapy-Type Structure

Fill from:

- `01-data/processed/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/112_patient_treatment_profile_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/115_patient_treatment_profile_v1_summary.tsv`

- `[fill in single vs mixed therapy-type findings]`
- `[fill in missing-type-only finding]`
- `[fill in dominant-type identifiability finding]`
- `[fill in whether compound raw therapy-type labels remain visible]`

## Regimen Context And Timing Coverage

Fill from:

- `01-data/processed/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/<PATIENT_TREATMENT_PROFILE_V1_RUN_ID>/patient_treatment_profile_v1_conflict_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/112_patient_treatment_profile_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/113_patient_treatment_profile_v1_conflict_audit.tsv`

- `[fill in regimen-context coverage]`
- `[fill in single vs mixed regimen-context finding]`
- `[fill in timing coverage]`
- `[fill in any aggregated timing inversion note]`

## Whether The Repo Is Ready For The Next Step

This section must stay explicit that readiness means coarse exploratory grouping only, not treatment-arm freezing or modeling.

- `[fill in whether coarse exploratory treatment grouping now looks feasible]`
- `[fill in whether treatment-arm freeze review is still blocked]`
- `[fill in the main rationale for the current readiness state]`

## What Still Blocks Treatment-Arm Freeze Review

- `[fill in drug-name normalization blocker]`
- `[fill in remaining manual-review burden]`
- `[fill in any additional structural ambiguity note]`

## Validation Checks Completed

- [ ] Patient treatment profile latest pointer resolves to the reviewed run.
- [ ] Reviewed run references the intended treatment-OS overlap v1 run, OS endpoint v1 run, and clinical biotab parse run.
- [ ] `patient_treatment_profile_v1.tsv`, `patient_treatment_profile_v1_spec.tsv`, `patient_treatment_profile_v1_summary.tsv`, `patient_treatment_profile_v1_conflict_audit.tsv`, and `run_log.json` all exist and are non-empty.
- [ ] `run_log.json` reports `status == completed`.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `required_upstream_pointers_found == true`.
- [ ] `run_log.json` reports `required_source_tables_found == true`.
- [ ] `run_log.json` reports `treatment_os_overlap_run_log_completed == true`.
- [ ] `run_log.json` reports `os_endpoint_run_log_completed == true`.
- [ ] `run_log.json` reports `clinical_biotabs_run_log_completed == true`.
- [ ] `run_log.json` reports `profile_row_count_matches_os_cohort == true`.
- [ ] `run_log.json` reports `row_order_preserved_from_os_endpoint == true`.
- [ ] `run_log.json` reports `all_os_patients_in_profile == true`.
- [ ] `run_log.json` reports `no_profile_duplicate_barcodes == true`.
- [ ] `run_log.json` reports `summary_counts_reconcile_to_patient_table == true`.
- [ ] `run_log.json` reports `conflict_audit_row_count_matches_manual_review == true`.
- [ ] `run_log.json` reports `spec_row_count_matches_profile_columns == true`.
- [ ] `run_log.json` reports `overlap_continuity_profile_vs_overlap == true`.
- [ ] `run_log.json` reports `overlap_therapy_type_summary_reconciles == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] `115_patient_treatment_profile_v1_summary.tsv` reports `provisional_readiness_interpretation == ready_for_coarse_exploratory_grouping_next` or `not_ready`, as actually saved.
- [ ] `115_patient_treatment_profile_v1_summary.tsv` reports `treatment_arm_freeze_review_status == blocked`.
- [ ] `115_patient_treatment_profile_v1_summary.tsv` reports `treatment_recommendation_modeling_status == out_of_scope`.
- [ ] Review tables `112` through `115` were regenerated from saved outputs on disk only.
- [ ] This note remains a treatment-profile construction artifact only and does not present final treatment arms or treatment recommendations.

## Current Interpretation

This note remains a patient treatment profile v1 review artifact only. It should answer whether the saved TCGA-BRCA treatment tables can now be represented as one row per OS cohort patient while keeping mixed, ambiguous, no-drug-record, and missing-type cases explicit. It must not present the current layer as a frozen treatment-arm dataset, a treatment recommendation layer, or a causal treatment-analysis result.
