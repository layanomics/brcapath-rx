# TCGA-BRCA Patient Treatment Grouping V1 Findings

This note records human-reviewed findings from the TCGA-BRCA provisional patient-level treatment grouping v1 workflow.

Important reminders:

- this document reviews provisional patient-level treatment grouping only
- this document does not normalize drug names or freeze final treatment arms
- this document does not perform treatment-effect estimation, causal analysis, or modeling
- this document does not replace raw treatment tables, prior patient treatment profile outputs, prior OS endpoint outputs, saved audit TSVs, or run logs
- this document remains TCGA-only and source-audited

## Reviewed Patient Treatment Grouping V1 Run ID

- Patient treatment grouping latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_grouping_v1_latest.json`
- Reviewed patient treatment grouping v1 run id: `[fill in patient treatment grouping v1 run id]`
- Source patient treatment profile v1 run id: `[fill in patient treatment profile v1 run id]`
- Source treatment-OS overlap v1 run id: `[fill in treatment-OS overlap v1 run id]`
- Source OS endpoint v1 run id: `[fill in OS endpoint v1 run id]`
- Source baseline model-input v1 run id: `[fill in baseline model-input v1 run id]`
- Source cohort v1 build id: `[fill in cohort v1 build id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Patient Treatment Grouping V1 Inputs Used

- Patient treatment grouping latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_grouping_v1_latest.json`
- Patient treatment profile latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_profile_v1_latest.json`
- OS endpoint v1 latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_os_endpoint_v1_latest.json`
- Baseline model-input v1 latest pointer: `01-data/audit/tcga-brca/model-input/tcga_brca_baseline_model_input_v1_latest.json`
- Processed patient treatment grouping TSV: `01-data/processed/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1.tsv`
- Audit patient treatment grouping spec TSV: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_spec.tsv`
- Audit patient treatment grouping arm summary TSV: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_arm_summary.tsv`
- Audit patient treatment grouping conflict TSV: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_conflict_audit.tsv`
- Audit patient treatment grouping summary TSV: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_summary.tsv`
- Audit run log: `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
- `116_patient_treatment_grouping_v1.tsv`
- `117_patient_treatment_grouping_v1_arm_summary.tsv`
- `118_patient_treatment_grouping_v1_conflict_audit.tsv`
- `119_patient_treatment_grouping_v1_spec.tsv`
- `120_patient_treatment_grouping_v1_summary.tsv`

## Grouping Rules Used

Fill from:

- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_spec.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/119_patient_treatment_grouping_v1_spec.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/120_patient_treatment_grouping_v1_summary.tsv`

- `[fill in unit-of-analysis statement]`
- `[fill in exact approved grouping categories statement]`
- `[fill in no-drug-record rule summary]`
- `[fill in mixed-multi-type handling summary]`
- `[fill in single-type exact raw dominant-type mapping summary]`
- `[fill in ancillary-or-other fallback summary]`
- `[fill in manual-review scope summary]`

## Group Counts

Fill from:

- `01-data/processed/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_arm_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/116_patient_treatment_grouping_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/117_patient_treatment_grouping_v1_arm_summary.tsv`

- `[fill in cohort size]`
- `[fill in no_drug_record count]`
- `[fill in single_chemotherapy count]`
- `[fill in single_hormone_therapy count]`
- `[fill in single_targeted_therapy count]`
- `[fill in single_immunotherapy count]`
- `[fill in single_ancillary_or_other count]`
- `[fill in mixed_multi_type count]`
- `[fill in missing_type_only count]`

## Manual Review Burden

Fill from:

- `01-data/processed/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_conflict_audit.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/118_patient_treatment_grouping_v1_conflict_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/120_patient_treatment_grouping_v1_summary.tsv`

- `[fill in total manual-review count and fraction]`
- `[fill in whether the upstream 357-patient burden was preserved]`
- `[fill in mixed_multi_type manual-review burden]`
- `[fill in missing_type_only manual-review burden]`
- `[fill in compound raw-type visibility note]`

## Radiation Overlap And Coverage

Fill from:

- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_arm_summary.tsv`
- `01-data/audit/tcga-brca/treatment-prep/patient_treatment_grouping_v1_runs/<PATIENT_TREATMENT_GROUPING_V1_RUN_ID>/patient_treatment_grouping_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/117_patient_treatment_grouping_v1_arm_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/120_patient_treatment_grouping_v1_summary.tsv`

- `[fill in overall radiation overlap count and fraction]`
- `[fill in radiation overlap by group]`
- `[fill in regimen-context coverage by group]`
- `[fill in treatment-timing coverage by group]`

## Whether Grouped Descriptive Review Is Now Feasible

This section must stay explicit that readiness means grouped descriptive review only, not treatment-arm freeze or modeling.

- `[fill in whether grouped descriptive review is now feasible]`
- `[fill in saved readiness interpretation]`
- `[fill in the main rationale for the current readiness state]`

## What Still Blocks Treatment-Arm Freeze Review

- `[fill in drug-name normalization blocker]`
- `[fill in unresolved manual-review burden blocker]`
- `[fill in mixed-case handling blocker]`

## Validation Checks Completed

- [ ] Patient treatment grouping latest pointer resolves to the reviewed run.
- [ ] Reviewed run references the intended patient treatment profile v1 run, treatment-OS overlap v1 run, OS endpoint v1 run, baseline model-input v1 run, and cohort v1 build id.
- [ ] `patient_treatment_grouping_v1.tsv`, `patient_treatment_grouping_v1_spec.tsv`, `patient_treatment_grouping_v1_arm_summary.tsv`, `patient_treatment_grouping_v1_conflict_audit.tsv`, `patient_treatment_grouping_v1_summary.tsv`, and `run_log.json` all exist and are non-empty.
- [ ] `run_log.json` reports `status == completed`.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `required_upstream_pointers_found == true`.
- [ ] `run_log.json` reports `required_source_tables_found == true`.
- [ ] `run_log.json` reports `patient_treatment_profile_run_log_completed == true`.
- [ ] `run_log.json` reports `os_endpoint_run_log_completed == true`.
- [ ] `run_log.json` reports `baseline_model_input_run_log_completed == true`.
- [ ] `run_log.json` reports `grouping_row_count_matches_os_cohort == true`.
- [ ] `run_log.json` reports `row_order_preserved_from_os_endpoint == true`.
- [ ] `run_log.json` reports `all_os_patients_in_grouping == true`.
- [ ] `run_log.json` reports `no_grouping_duplicate_barcodes == true`.
- [ ] `run_log.json` reports `group_categories_restricted_to_approved_set == true`.
- [ ] `run_log.json` reports `manual_review_scope_matches_upstream_profile == true`.
- [ ] `run_log.json` reports `group_counts_reconcile_to_grouping_table == true`.
- [ ] `run_log.json` reports `arm_summary_counts_reconcile_to_grouping_table == true`.
- [ ] `run_log.json` reports `summary_counts_reconcile_to_grouping_table == true`.
- [ ] `run_log.json` reports `conflict_audit_row_count_matches_manual_review == true`.
- [ ] `run_log.json` reports `spec_row_count_matches_grouping_columns == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] `120_patient_treatment_grouping_v1_summary.tsv` reports `provisional_readiness_interpretation == ready_for_grouped_descriptive_review`.
- [ ] `120_patient_treatment_grouping_v1_summary.tsv` reports `treatment_arm_freeze_review_status == blocked`.
- [ ] `120_patient_treatment_grouping_v1_summary.tsv` reports `treatment_recommendation_modeling_status == out_of_scope`.
- [ ] Review tables `116` through `120` were regenerated from saved outputs on disk only.
- [ ] This note remains a provisional grouping artifact only and does not present final treatment arms or treatment recommendations.

## Current Interpretation

This note remains a provisional patient-level treatment grouping v1 review artifact only. It should answer whether the saved TCGA-BRCA patient treatment profile can now support coarse grouped descriptive review while keeping no-drug, mixed-multi-type, missing-type-only, and ancillary-or-other fallback cases explicit. It must not present the current layer as a frozen treatment-arm dataset, a treatment recommendation layer, or a treatment-effect analysis result.
