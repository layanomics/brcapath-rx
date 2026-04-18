# TCGA-BRCA Baseline Model-Input V1 Findings

This note records human-reviewed findings from the TCGA-BRCA baseline model-input v1 workflow.

Important reminders:

- this document reviews a baseline model-input preparation layer only
- this document covers first-pass baseline model-input preparation only
- this document does not perform model training
- this document does not freeze the final endpoint
- this document does not add treatment detail back into the patient-level matrix
- this document does not replace source files, parsed TSVs, prior audit TSVs, prior baseline-analysis outputs, prior baseline-profile outputs, prior baseline feature-set outputs, or run logs
- this document does not include new raw parsing, XML parsing, SSF parsing, METABRIC work, treatment modeling, endpoint freeze, or outcome modeling

## Reviewed Baseline Model-Input V1 Run ID

- Baseline model-input latest pointer: `01-data/audit/tcga-brca/model-input/tcga_brca_baseline_model_input_v1_latest.json`
- Reviewed baseline model-input v1 run id: `[fill in baseline model-input v1 run id]`
- Source baseline feature-set v1 run id: `[fill in baseline feature-set v1 run id]`
- Source baseline-profile v1 run id: `[fill in baseline-profile v1 run id]`
- Source baseline-analysis-prep v1 run id: `[fill in baseline-analysis-prep v1 run id]`
- Source minimal cohort v1 build id: `[fill in cohort v1 build id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Baseline Model-Input V1 Inputs Used

- Baseline model-input latest pointer: `01-data/audit/tcga-brca/model-input/tcga_brca_baseline_model_input_v1_latest.json`
- Baseline feature-set latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_feature_set_v1_latest.json`
- Processed model-input matrix TSV: `01-data/processed/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1.tsv`
- Audit feature dictionary TSV: `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_feature_dictionary.tsv`
- Audit encoding spec TSV: `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_encoding_spec.tsv`
- Audit missingness actions TSV: `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_missingness_actions.tsv`
- Audit summary TSV: `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_summary.tsv`
- Audit run log: `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/run_log.json`
- Upstream feature-set audit map TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_feature_set_v1_runs/<BASELINE_FEATURE_SET_V1_RUN_ID>/baseline_feature_set_v1_audit_map.tsv`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `91_baseline_model_input_v1_feature_dictionary.tsv`
  - `92_baseline_model_input_v1_encoding_spec.tsv`
  - `93_baseline_model_input_v1_missingness_actions.tsv`
  - `94_baseline_model_input_v1_matrix_preview.tsv`
  - `95_baseline_model_input_v1_manual_review_fields.tsv`
  - `96_baseline_model_input_v1_summary.tsv`

## Final Matrix Shape

Fill from:

- `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/96_baseline_model_input_v1_summary.tsv`

Record only saved matrix statements.

- `[fill in final row count]`
- `[fill in final output column count]`
- `[fill in patient/case unit statement]`
- `[fill in row-order / audit-linkage statement]`

## Feature Typing Overview

Fill from:

- `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_feature_dictionary.tsv`
- `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/91_baseline_model_input_v1_feature_dictionary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/96_baseline_model_input_v1_summary.tsv`

- `[fill in count of binary source features]`
- `[fill in count of nominal categorical source features]`
- `[fill in count of numeric continuous source features]`
- `[fill in count of numeric count source features]`
- `[fill in any zero-count feature type note]`

## Encoding Plan

Fill from:

- `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_encoding_spec.tsv`
- `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/92_baseline_model_input_v1_encoding_spec.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/96_baseline_model_input_v1_summary.tsv`

- `[fill in numeric passthrough summary]`
- `[fill in binary map summary]`
- `[fill in exact-label one-hot summary]`
- `[fill in reversibility statement]`

## Missingness-Action Plan

Fill from:

- `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_missingness_actions.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/93_baseline_model_input_v1_missingness_actions.tsv`

- `[fill in keep-as-is fields note]`
- `[fill in allow-missing-category fields note]`
- `[fill in eligible-for-simple-imputation-later fields note]`
- `[fill in no-imputation-in-v1 statement]`

## Fields Still Needing Manual Review

Fill from:

- `01-data/audit/tcga-brca/model-input/baseline_v1_runs/<BASELINE_MODEL_INPUT_V1_RUN_ID>/baseline_model_input_v1_feature_dictionary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/95_baseline_model_input_v1_manual_review_fields.tsv`

- `[fill in AJCC review note]`
- `[fill in anatomic subdivision review note]`
- `[fill in histology review note]`
- `[fill in any additional reviewer interpretation]`

## What This Matrix Is Ready For

Summarize only what this saved matrix honestly supports now.

- `[fill in first-pass baseline modeling readiness]`
- `[fill in encoding review readiness]`
- `[fill in missingness review readiness]`
- `[fill in reversible audit-ready representation statement]`

## Endpoint And Treatment Caveats

This section must remain explicit that the workflow is model-input prep only.

- `[fill in endpoint-freeze blocker]`
- `[fill in treatment exclusion statement]`
- `[fill in model-input-prep-only statement]`

## Validation Checks Completed

- [ ] Baseline model-input latest pointer resolves to the reviewed run.
- [ ] Reviewed run references the intended baseline feature-set v1 run.
- [ ] Matrix TSV, feature dictionary TSV, encoding spec TSV, missingness actions TSV, summary TSV, and run log all exist and are non-empty.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `required_upstream_pointers_found == true`.
- [ ] `run_log.json` reports `required_source_tables_found == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] `run_log.json` reports `row_order_preserved_from_baseline_feature_set_v1 == true`.
- [ ] Final matrix row count is greater than zero.
- [ ] Final matrix row count equals the saved `baseline_feature_set_v1.tsv` row count.
- [ ] Feature dictionary row count equals the final matrix column count.
- [ ] Encoding spec output fields cover every final matrix column.
- [ ] Missingness action rows reconcile exactly to the 25 input source features.
- [ ] `96_baseline_model_input_v1_summary.tsv` reports `endpoint_freeze_status == blocked`.
- [ ] `96_baseline_model_input_v1_summary.tsv` reports `treatment_inclusion_status == excluded`.
- [ ] Review tables `91` through `96` were regenerated from model-input outputs on disk only.
- [ ] This note remains a baseline model-input v1 review only and does not perform model training or final endpoint freeze.

## Current Interpretation

This note remains a baseline model-input v1 review artifact only. It should answer whether the saved 25-feature baseline clinical set can be turned into a clean first-pass model-input matrix, which fields are straightforward versus still flagged for manual review, and whether the saved matrix is ready for first-pass baseline-only modeling experiments. It must not present this layer as model training, a frozen endpoint dataset, a treatment-ready dataset, or a stronger clinical-effect analysis layer.
