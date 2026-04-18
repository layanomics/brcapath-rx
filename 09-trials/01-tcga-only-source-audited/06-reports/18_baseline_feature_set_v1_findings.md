# TCGA-BRCA Baseline Feature-Set V1 Findings

This note records human-reviewed findings from the TCGA-BRCA baseline feature-set v1 workflow.

Important reminders:

- this document reviews a baseline feature-set preparation layer only
- this document covers first-pass baseline model-input preparation only
- this document does not perform model training
- this document does not freeze the final endpoint
- this document does not add treatment detail back into the patient-level table
- this document does not replace source files, parsed TSVs, prior audit TSVs, prior baseline-analysis outputs, prior baseline-profile outputs, feature-set outputs, or run logs
- this document does not include new raw parsing, XML parsing, SSF parsing, METABRIC work, treatment modeling, endpoint freeze, or outcome modeling

## Reviewed Baseline Feature-Set V1 Run ID

- Baseline feature-set latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_feature_set_v1_latest.json`
- Reviewed baseline feature-set v1 run id: `[fill in baseline feature-set v1 run id]`
- Source baseline-profile run id: `[fill in baseline-profile run id]`
- Source baseline-analysis-prep run id: `[fill in baseline-analysis-prep run id]`
- Source minimal cohort v1 build id: `[fill in cohort v1 build id]`
- Source ambiguity-resolution run id: `[fill in ambiguity-resolution run id]`
- Source clinical shortlist run id: `[fill in shortlist run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Baseline Feature-Set V1 Inputs Used

- Baseline feature-set latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_feature_set_v1_latest.json`
- Baseline-profile latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_profile_v1_latest.json`
- Baseline-analysis-prep latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_analysis_v1_latest.json`
- Processed feature-set TSV: `01-data/processed/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1.tsv`
- Audit feature-set spec TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_spec.tsv`
- Audit feature-set missingness TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_missingness.tsv`
- Audit feature-set summary TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_summary.tsv`
- Audit feature-set audit-map TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_audit_map.tsv`
- Audit run log: `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `86_baseline_feature_set_v1_preview.tsv`
  - `87_baseline_feature_set_v1_spec.tsv`
  - `88_baseline_feature_set_v1_missingness.tsv`
  - `89_baseline_feature_set_v1_decision_summary.tsv`
  - `90_baseline_feature_set_v1_summary.tsv`

## Final Row Count

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/90_baseline_feature_set_v1_summary.tsv`

Record only saved baseline feature-set statements.

- `[fill in final row count]`
- `[fill in patient/case unit statement]`
- `[fill in baseline-analysis carry-forward statement]`

## Final Included Feature List

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_spec.tsv`
- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/87_baseline_feature_set_v1_spec.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/90_baseline_feature_set_v1_summary.tsv`

Record the final included baseline clinical fields exactly as saved.

- `[fill in included feature count]`
- `[fill in auto-included feature summary]`
- `[fill in review-decided included feature summary]`
- `[fill in included feature names statement]`

## Review-Decided Field Outcomes

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_spec.tsv`
- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/89_baseline_feature_set_v1_decision_summary.tsv`

Be explicit about the saved include/exclude decision for each review-needed baseline clinical field.

- `[fill in ethnicity outcome]`
- `[fill in gender outcome]`
- `[fill in icd_10 outcome]`
- `[fill in icd_o_3_site outcome]`

## Missingness Observations

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_missingness.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/88_baseline_feature_set_v1_missingness.tsv`

Record only saved missingness observations. Do not impute or reinterpret values in this note.

- `[fill in highest-missingness retained feature note]`
- `[fill in lowest-missingness retained feature note]`
- `[fill in any notable review-decided missingness note]`

## What This Feature Set Is Ready For

Summarize only what this saved feature-set honestly supports now.

- `[fill in first-pass baseline model-input readiness]`
- `[fill in descriptive baseline review readiness]`
- `[fill in missingness review readiness]`

## What Is Still Deliberately Excluded

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_spec.tsv`
- `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/87_baseline_feature_set_v1_spec.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/90_baseline_feature_set_v1_summary.tsv`

Make the deliberate exclusions explicit.

- `[fill in audit/id exclusion note]`
- `[fill in endpoint exclusion note]`
- `[fill in follow-up evidence exclusion note]`
- `[fill in biospecimen evidence exclusion note]`
- `[fill in fixed review-field exclusion note]`

## Endpoint And Treatment Caveats

This section must remain explicit that the workflow is feature-set prep only.

- `[fill in endpoint-freeze blocker]`
- `[fill in treatment exclusion statement]`
- `[fill in feature-set-prep-only statement]`

## Validation Checks Completed

- [ ] Baseline feature-set latest pointer resolves to the reviewed run.
- [ ] Reviewed run references the intended baseline-profile and baseline-analysis-prep runs.
- [ ] Feature-set TSV, spec TSV, missingness TSV, summary TSV, audit-map TSV, and run log all exist and are non-empty.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] Final feature-set row count is greater than zero.
- [ ] Final feature-set row count equals the saved `baseline_analysis_v1` row count.
- [ ] Included spec rows reconcile exactly to the saved feature-set output columns.
- [ ] Audit-map row count equals the saved feature-set row count.
- [ ] Review tables `86` through `90` were regenerated from feature-set outputs on disk only.
- [ ] This note remains a baseline feature-set v1 review only and does not perform model training or final endpoint freeze.

## Current Interpretation

This note remains a baseline feature-set v1 review artifact only. It should answer which exact baseline clinical fields are included in the first feature-set v1 table, which review-needed baseline fields were included or excluded and why, and whether the saved table is ready for first-pass baseline model input. It must not present this layer as model training, a frozen endpoint dataset, a treatment-ready dataset, or a stronger clinical-effect analysis layer.
