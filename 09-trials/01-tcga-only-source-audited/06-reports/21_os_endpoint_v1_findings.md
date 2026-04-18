# TCGA-BRCA OS Endpoint V1 Findings

This note records human-reviewed findings from the TCGA-BRCA OS endpoint v1 freeze workflow.

Important reminders:

- this document reviews OS endpoint v1 only
- this document does not claim final universal endpoint truth
- this document does not perform survival modeling
- this document does not add treatment detail
- this document does not replace raw files, XML-derived prep tables, prior audit outputs, prior cohort outputs, or prior model-input outputs
- this document does not include PFS, DFS, recurrence, treatment, or METABRIC work

## Reviewed OS Endpoint V1 Run ID

- OS endpoint v1 latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_os_endpoint_v1_latest.json`
- Reviewed OS endpoint v1 run id: `[fill in os endpoint v1 run id]`
- Source endpoint-target prep v1 run id: `[fill in endpoint-target prep v1 run id]`
- Source minimal cohort v1 build id: `[fill in cohort v1 build id]`
- Source baseline model-input v1 run id: `[fill in baseline model-input v1 run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## OS Endpoint V1 Inputs Used

- OS endpoint v1 latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_os_endpoint_v1_latest.json`
- Endpoint-target prep latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_endpoint_target_prep_v1_latest.json`
- Minimal cohort v1 latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_minimal_cohort_v1_latest.json`
- Baseline model-input v1 latest pointer: `01-data/audit/tcga-brca/model-input/tcga_brca_baseline_model_input_v1_latest.json`
- Processed OS endpoint TSV: `01-data/processed/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1.tsv`
- Audit OS rule spec TSV: `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_rule_spec.tsv`
- Audit OS conflict TSV: `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_conflict_audit.tsv`
- Audit OS summary TSV: `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_summary.tsv`
- Audit OS run log: `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `102_os_endpoint_v1.tsv`
  - `103_os_endpoint_v1_conflict_audit.tsv`
  - `104_os_endpoint_v1_rule_spec.tsv`
  - `105_os_endpoint_v1_manual_review_cases.tsv`
  - `106_os_endpoint_v1_summary.tsv`

## Rule Used

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_rule_spec.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/104_os_endpoint_v1_rule_spec.tsv`

- `[fill in event rule summary]`
- `[fill in time rule summary]`
- `[fill in conflict policy summary]`
- `[fill in inclusion-status policy summary]`

## Event / Censor Counts

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/102_os_endpoint_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/106_os_endpoint_v1_summary.tsv`

- `[fill in final row count]`
- `[fill in os_event = 1 count]`
- `[fill in os_event = 0 count]`
- `[fill in non-missing os_time_days count]`
- `[fill in included clean vs included conflict vs excluded counts]`

## Time-Source Distribution

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/102_os_endpoint_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/106_os_endpoint_v1_summary.tsv`

- `[fill in patients using days_to_death]`
- `[fill in patients using censored last-contact time]`
- `[fill in whether follow-up max last-contact often exceeds patient header]`
- `[fill in whether multiple follow-up versions materially contributed endpoint evidence]`

## Conflict Observations

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_conflict_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/103_os_endpoint_v1_conflict_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/105_os_endpoint_v1_manual_review_cases.tsv`

- `[fill in vital-status disagreement count and pattern]`
- `[fill in dead-without-death-time count and note]`
- `[fill in whether any no-usable-time or no-clear-status cases remain]`
- `[fill in multiple-source reconciliation note]`

## Manual-Review Burden

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/os_endpoint_v1_runs/<OS_ENDPOINT_V1_RUN_ID>/os_endpoint_v1_conflict_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/105_os_endpoint_v1_manual_review_cases.tsv`

- `[fill in patients requiring manual review]`
- `[fill in whether the manual-review set is mostly vital-status disagreement]`
- `[fill in whether provisional values were retained under Include + Flag]`

## Whether OS Endpoint V1 Is Ready To Join To X

This section must stay explicit that readiness means auditable first-pass join readiness only.

- `[fill in whether join keys are complete for all rows]`
- `[fill in whether included clean rows are sufficient for first-pass X + Y experiments]`
- `[fill in whether included conflict-flagged rows should be optionally filtered]`

## What Still Blocks Stronger Endpoint Claims

- `[fill in remaining manual reconciliation blocker]`
- `[fill in why this is still not stronger endpoint truth]`
- `[fill in unchanged treatment / design caveat]`

## Validation Checks Completed

- [ ] OS endpoint v1 latest pointer resolves to the reviewed run.
- [ ] Reviewed run references the intended endpoint-target prep v1 run, minimal cohort v1 build, and baseline model-input v1 run.
- [ ] `os_endpoint_v1.tsv`, `os_endpoint_v1_rule_spec.tsv`, `os_endpoint_v1_conflict_audit.tsv`, `os_endpoint_v1_summary.tsv`, and `run_log.json` all exist and are non-empty.
- [ ] `run_log.json` reports `status == completed`.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `required_source_tables_found == true`.
- [ ] `run_log.json` reports `os_row_count_matches_minimal_cohort == true`.
- [ ] `run_log.json` reports `os_row_count_matches_baseline_model_input == true`.
- [ ] `run_log.json` reports `patient_ids_unique == true`.
- [ ] `run_log.json` reports `manual_review_patients_covered_by_conflict_audit == true`.
- [ ] `run_log.json` reports `summary_counts_reconcile_to_os_table == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] `106_os_endpoint_v1_summary.tsv` reports `os_endpoint_v1_readiness_interpretation == join_ready_with_explicit_conflict_flags`.
- [ ] `106_os_endpoint_v1_summary.tsv` reports `stronger_endpoint_claims_status == blocked_pending_manual_reconciliation`.
- [ ] `106_os_endpoint_v1_summary.tsv` reports `treatment_feasibility_status == unchanged_not_addressed`.
- [ ] Review tables `102` through `106` were regenerated from OS endpoint outputs on disk only.
- [ ] This note remains OS endpoint v1 only and does not present final universal endpoint truth or a modeling result.

## Current Interpretation

This note remains an OS endpoint v1 review artifact only. It should answer whether the saved rule-based freeze creates a usable patient-level OS target with explicit provenance and conflict flags, how many patients remain conflicted or weakly supported, and whether the table looks join-ready for first-pass prognostic experiments with the current baseline X matrix. It must not present OS endpoint v1 as final endpoint truth, a treatment-ready dataset, or a finished survival-analysis result.
