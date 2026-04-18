# TCGA-BRCA Baseline Profile V1 Findings

This note records human-reviewed findings from the TCGA-BRCA baseline-profile v1 workflow.

Important reminders:

- this document reviews a descriptive baseline-profile layer only
- this document covers candidate-field selection only
- this document does not perform model training
- this document does not freeze the final feature set
- this document does not freeze the final endpoint
- this document does not add treatment detail back into the patient-level table
- this document does not replace source files, parsed TSVs, prior audit TSVs, prior baseline-analysis outputs, profile outputs, or run logs
- this document does not include new raw parsing, XML parsing, SSF parsing, METABRIC work, treatment modeling, endpoint freeze, or outcome modeling

## Reviewed Baseline-Profile Run ID

- Baseline-profile latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_profile_v1_latest.json`
- Reviewed baseline-profile run id: `[fill in baseline-profile run id]`
- Source baseline-analysis-prep run id: `[fill in baseline-analysis-prep run id]`
- Source minimal cohort v1 build id: `[fill in cohort v1 build id]`
- Source ambiguity-resolution run id: `[fill in ambiguity-resolution run id]`
- Source clinical shortlist run id: `[fill in shortlist run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Baseline-Profile Inputs Used

- Baseline-analysis-prep latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_analysis_v1_latest.json`
- Minimal cohort v1 latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_minimal_cohort_v1_latest.json`
- Ambiguity-resolution latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_blueprint_ambiguity_resolution_latest.json`
- Clinical shortlist latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_shortlist_latest.json`
- Processed field summary TSV: `01-data/processed/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_field_summary.tsv`
- Processed missingness-ranked TSV: `01-data/processed/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_missingness_ranked.tsv`
- Processed value-summary TSV: `01-data/processed/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_value_summary.tsv`
- Audit candidate-fields TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_candidate_fields.tsv`
- Audit excluded-fields TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_excluded_fields.tsv`
- Audit summary TSV: `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_summary.tsv`
- Audit run log: `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `80_baseline_profile_v1_field_summary.tsv`
  - `81_baseline_profile_v1_missingness_ranked.tsv`
  - `82_baseline_profile_v1_value_summary.tsv`
  - `83_baseline_profile_v1_candidate_fields.tsv`
  - `84_baseline_profile_v1_excluded_fields.tsv`
  - `85_baseline_profile_v1_summary.tsv`

## Final Row Count

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/85_baseline_profile_v1_summary.tsv`

Record only saved baseline-profile statements.

- `[fill in final row count]`
- `[fill in patient/case unit statement]`
- `[fill in baseline-analysis source-table statement]`

## Retained Field Categories

Fill from:

- `01-data/processed/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_field_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/80_baseline_profile_v1_field_summary.tsv`
- `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_summary.tsv`

Record the retained categories exactly as saved.

- `[fill in audit/id field summary]`
- `[fill in baseline clinical field summary]`
- `[fill in follow-up join-evidence field summary]`
- `[fill in endpoint candidate field summary]`
- `[fill in biospecimen sample-anchor evidence summary]`
- `[fill in derived prep-flag summary]`

## Highest-Value Baseline Clinical Fields

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_candidate_fields.tsv`
- `01-data/processed/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_value_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/82_baseline_profile_v1_value_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/83_baseline_profile_v1_candidate_fields.tsv`

Be explicit that these remain first-pass candidate fields only.

- `[fill in strongest direct-candidate field group]`
- `[fill in strongest receptor / biomarker candidate group]`
- `[fill in strongest stage / pathology candidate group]`

## Missingness Observations

Fill from:

- `01-data/processed/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_missingness_ranked.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/81_baseline_profile_v1_missingness_ranked.tsv`

Record only saved missingness observations. Do not impute or reinterpret values in this note.

- `[fill in highest-missingness retained field note]`
- `[fill in baseline-clinical completeness note]`
- `[fill in near-constant review-needed note]`

## Weak Or Excluded Fields

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_candidate_fields.tsv`
- `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_excluded_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/83_baseline_profile_v1_candidate_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/84_baseline_profile_v1_excluded_fields.tsv`

Separate review-needed baseline fields from fields excluded for now.

- `[fill in excluded audit/id field note]`
- `[fill in excluded derived-flag note]`
- `[fill in any sparse or low-information note]`

## Endpoint And Biospecimen Fields Kept Out Of First Modeling Prep

Fill from:

- `01-data/audit/tcga-brca/analysis-prep/baseline_profile_v1_runs/<BASELINE_PROFILE_V1_RUN_ID>/baseline_profile_v1_candidate_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/83_baseline_profile_v1_candidate_fields.tsv`

Be explicit that these remain review-only and must stay out of first-pass baseline modeling prep.

- `[fill in endpoint-candidate review note]`
- `[fill in follow-up join-evidence review note]`
- `[fill in biospecimen sample-anchor review note]`

## What The Table Is Ready For Now

Summarize only what this baseline-profile layer honestly supports now.

- `[fill in descriptive reporting readiness]`
- `[fill in missingness review readiness]`
- `[fill in first-pass baseline feature-selection readiness]`

## What Still Blocks Endpoint Freeze Or Treatment Modeling

Summarize only what remains unresolved after this baseline-profile layer. This section must stay at review level only.

- `[fill in endpoint-freeze blocker]`
- `[fill in treatment-modeling blocker]`
- `[fill in any remaining ambiguity blocker]`

## Validation Checks Completed

- [ ] Baseline-profile latest pointer resolves to the reviewed run.
- [ ] Reviewed run references the intended baseline-analysis-prep, cohort v1, ambiguity-resolution, and clinical shortlist runs.
- [ ] Field summary TSV, missingness-ranked TSV, value-summary TSV, candidate-fields TSV, excluded-fields TSV, summary TSV, and run log all exist and are non-empty.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] Final profile row count is greater than zero.
- [ ] Final profile row count equals the saved `baseline_analysis_v1` row count.
- [ ] Candidate buckets reconcile exactly to the profiled retained fields.
- [ ] Endpoint freeze remains blocked in the saved summary.
- [ ] Treatment inclusion remains excluded in the saved summary.
- [ ] Review tables `80` through `85` were regenerated from baseline-profile outputs on disk only.
- [ ] This note remains a profiling and candidate-field-selection review only and does not perform model training or final feature freeze.

## Current Interpretation

This note remains a baseline-profile v1 review artifact only. It should answer which retained baseline fields look strong enough for first-pass baseline feature selection, which fields remain weak or provisional, and whether the current `baseline_analysis_v1` table is ready for descriptive reporting and first-pass baseline feature selection. It must not present this layer as model training, a frozen endpoint dataset, a treatment-ready dataset, or a final feature freeze.
