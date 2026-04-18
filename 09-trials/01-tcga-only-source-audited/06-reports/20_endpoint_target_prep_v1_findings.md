# TCGA-BRCA Endpoint-Target Prep V1 Findings

This note records human-reviewed findings from the TCGA-BRCA clinical XML endpoint-target prep v1 workflow.

Important reminders:

- this document reviews endpoint-target preparation only
- this document does not freeze the final endpoint
- this document does not perform model training
- this document does not add treatment detail back into the patient-level matrix
- this document does not replace raw source files, parsed XML-derived TSVs, prior audit TSVs, prior cohort outputs, prior model-input outputs, or run logs
- this document covers only saved clinical XML patient-header and follow-up reconciliation outputs from the current TCGA-only source-audited restart phase
- this document does not include new raw downloads, METABRIC work, treatment grouping, survival modeling, or final endpoint selection

## Reviewed Endpoint-Target Prep V1 Run ID

- Endpoint-target prep latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_endpoint_target_prep_v1_latest.json`
- Reviewed endpoint-target prep v1 run id: `[fill in endpoint-target prep v1 run id]`
- XML source run used: `[fill in source run id]`
- Source endpoint crosswalk run id: `[fill in endpoint crosswalk run id]`
- Source ambiguity-resolution run id: `[fill in ambiguity-resolution run id]`
- Source minimal cohort v1 build id: `[fill in cohort v1 build id]`
- Source baseline model-input v1 run id: `[fill in baseline model-input v1 run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Endpoint-Target Prep Inputs Used

- Endpoint-target prep latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_endpoint_target_prep_v1_latest.json`
- Source supplements latest pointer: `01-data/audit/tcga-brca/source/tcga_brca_source_supplements_latest.json`
- Endpoint crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_endpoint_crosswalk_latest.json`
- Ambiguity-resolution latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_blueprint_ambiguity_resolution_latest.json`
- Minimal cohort v1 latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_minimal_cohort_v1_latest.json`
- Baseline model-input v1 latest pointer: `01-data/audit/tcga-brca/model-input/tcga_brca_baseline_model_input_v1_latest.json`
- Processed patient-header TSV: `01-data/processed/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/clinical_xml_patient_endpoint_fields.tsv`
- Processed follow-up long TSV: `01-data/processed/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/clinical_xml_followup_fields_long.tsv`
- Processed endpoint-target prep TSV: `01-data/processed/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1.tsv`
- Audit follow-up coverage TSV: `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/clinical_xml_followup_version_coverage.tsv`
- Audit overlap TSV: `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_overlap_audit.tsv`
- Audit summary TSV: `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_summary.tsv`
- Audit run log: `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `97_clinical_xml_patient_endpoint_fields.tsv`
  - `98_clinical_xml_followup_version_coverage.tsv`
  - `99_endpoint_target_prep_v1.tsv`
  - `100_endpoint_target_prep_v1_overlap_audit.tsv`
  - `101_endpoint_target_prep_v1_summary.tsv`

## Patient-Header Coverage

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/101_endpoint_target_prep_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/97_clinical_xml_patient_endpoint_fields.tsv`

- `[fill in parsed XML patient file count]`
- `[fill in patient-header vital_status coverage]`
- `[fill in patient-header days_to_last_followup coverage]`
- `[fill in patient-header days_to_death coverage]`
- `[fill in patient-header tumor-status coverage]`

## Follow-Up Coverage By Version

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/clinical_xml_followup_version_coverage.tsv`
- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/98_clinical_xml_followup_version_coverage.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/101_endpoint_target_prep_v1_summary.tsv`

- `[fill in patients with any XML follow-up]`
- `[fill in v1.5 patient and record coverage]`
- `[fill in v2.1 patient and record coverage]`
- `[fill in v4.0 patient and record coverage]`
- `[fill in older-version-only and overlapping-version notes]`

## New Coverage Gained Beyond Current Biotab v4.0

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/clinical_xml_followup_version_coverage.tsv`
- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/98_clinical_xml_followup_version_coverage.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/100_endpoint_target_prep_v1_overlap_audit.tsv`

- `[fill in patients with XML follow-up and no current grouped biotab v4 match]`
- `[fill in patients with any XML gain-over-v4 flag]`
- `[fill in whether the gain is mainly older-version coverage, later last-contact values, days_to_death evidence, or a mix]`

## Overlap / Disagreement Observations

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_overlap_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/100_endpoint_target_prep_v1_overlap_audit.tsv`

- `[fill in patient-header vs follow-up vital-status comparison note]`
- `[fill in patient-header vs follow-up last-contact comparison note]`
- `[fill in XML v4 vs current grouped biotab v4 agreement/disagreement note]`
- `[fill in any multi-version overlap note]`

## Days_To_Death Observations

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/99_endpoint_target_prep_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/100_endpoint_target_prep_v1_overlap_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/101_endpoint_target_prep_v1_summary.tsv`

- `[fill in count of patients with any non-missing XML days_to_death]`
- `[fill in whether days_to_death appears mainly in patient header, follow-up, or both]`
- `[fill in whether this materially improves later OS-style endpoint review]`

## Whether OS-Style Endpoint Freeze Is Now More Realistic

Fill from:

- `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/<ENDPOINT_TARGET_PREP_V1_RUN_ID>/endpoint_target_prep_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/101_endpoint_target_prep_v1_summary.tsv`

This section must stay explicit that the saved workflow is still endpoint-target preparation only.

- `[fill in readiness interpretation from saved summary]`
- `[fill in why endpoint freeze is closer but still blocked]`
- `[fill in what the patient-level prep table is now good enough for]`

## What Still Blocks Final Endpoint Freeze

Summarize only what remains unresolved after this endpoint-target prep layer. Keep this section at review level only.

- `[fill in unresolved precedence / reconciliation blocker]`
- `[fill in unresolved censoring / final endpoint-definition blocker]`
- `[fill in any remaining treatment or design blocker that was intentionally left unchanged]`

## Validation Checks Completed

- [ ] Endpoint-target prep latest pointer resolves to the reviewed run.
- [ ] Reviewed run references the intended XML source run, endpoint crosswalk run, ambiguity-resolution run, minimal cohort v1 build, and baseline model-input v1 run.
- [ ] Patient-header TSV, follow-up long TSV, follow-up coverage TSV, endpoint-target prep TSV, overlap audit TSV, summary TSV, and run log all exist and are non-empty.
- [ ] `run_log.json` reports `status == completed`.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `xml_files_found == true`.
- [ ] `run_log.json` reports `xml_patient_file_count_matches_minimal_cohort == true`.
- [ ] `run_log.json` reports `endpoint_prep_row_count_matches_minimal_cohort == true`.
- [ ] `run_log.json` reports `overlap_audit_row_count_matches_endpoint_prep == true`.
- [ ] `run_log.json` reports `row_bridge_fully_matched_to_baseline_feature_set_v1_audit_map == true`.
- [ ] `run_log.json` reports `followup_versions_restricted_to_supported_set == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] `101_endpoint_target_prep_v1_summary.tsv` reports `endpoint_freeze_status == blocked_pending_manual_reconciliation`.
- [ ] `101_endpoint_target_prep_v1_summary.tsv` reports `treatment_feasibility_status == unchanged_not_addressed`.
- [ ] Review tables `97` through `101` were regenerated from endpoint-target prep outputs on disk only.
- [ ] Mixed-version patient examples such as `TCGA-AO-A0J5` or `TCGA-AR-A0TT` were manually spot-checked from saved outputs.
- [ ] A v1.5-only patient such as `TCGA-BH-A0HF` was manually spot-checked from saved outputs.
- [ ] This note remains endpoint-target preparation only and does not freeze the endpoint or perform model training.

## Current Interpretation

This note remains an endpoint-target preparation review artifact only. It should answer whether the saved XML patient-header and follow-up parsing workflow materially improves endpoint-like coverage, whether the new patient-level endpoint-target prep table is good enough for manual endpoint-freeze review, and whether OS-style endpoint preparation now looks more realistic inside the current TCGA-only audited repo state. It must not present this layer as a final endpoint freeze, a treatment-ready dataset, or a modeling result.
