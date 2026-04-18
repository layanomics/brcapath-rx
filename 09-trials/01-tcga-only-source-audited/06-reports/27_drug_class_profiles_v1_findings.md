# TCGA-BRCA Drug-Class Profile Descriptive Review And Arm-Freeze Candidate Audit V1 Findings

This note records the first descriptive review of TCGA-BRCA patient-level drug-class profiles
and the first arm-freeze candidate audit derived from the provisional drug-name normalization layer.

Important reminders:

- this document is descriptive review and arm-freeze candidate audit only
- this document does not normalize drug names
- this document does not redefine treatment groups
- this document does not freeze treatment arms
- this document does not silently exclude mixed or unknown patients
- this document does not perform causal analysis, treatment-effect estimation, or modeling
- this document does not produce treatment recommendations
- bucket definitions come only from the saved patient drug-class profile, not from fresh raw re-derivation
- raw treatment tables remain unchanged and are not replaced by this review layer

## Reviewed Run

- Drug-class profile latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_drug_class_profile_v1_latest.json`
- Reviewed drug-class profile v1 run id: `[fill in drug_class_profile_v1_run_id]`
- Source drug-name normalization v1 run id: `[fill in drug_name_normalization_v1_run_id]`
- Source patient treatment grouping v1 run id: `[fill in patient_treatment_grouping_v1_run_id]`
- Source OS endpoint v1 run id: `[fill in os_endpoint_v1_run_id]`
- Source baseline model-input v1 run id: `[fill in baseline_model_input_v1_run_id]`
- Source baseline analysis v1 run id: `[fill in baseline_analysis_v1_run_id]`
- Source cohort v1 build id: `[fill in cohort_v1_build_id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Inputs Used

- Drug-class profile latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_drug_class_profile_v1_latest.json`
- Drug-name normalization latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_drug_name_normalization_v1_latest.json`
- Patient treatment grouping latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_grouping_v1_latest.json`
- Baseline analysis latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_analysis_v1_latest.json`
- Processed drug-class profile TSV: `01-data/processed/tcga-brca/treatment-prep/drug_name_normalization_v1_runs/<DRUG_NAME_NORMALIZATION_V1_RUN_ID>/patient_drug_class_profile_v1.tsv`
- Audit outputs in `01-data/audit/tcga-brca/treatment-prep/drug_class_profile_v1_runs/<DRUG_CLASS_PROFILE_V1_RUN_ID>/`:
  - `drug_class_profile_v1_group_summary.tsv`
  - `drug_class_profile_v1_confounder_coverage.tsv`
  - `drug_class_profile_v1_value_composition.tsv`
  - `drug_class_profile_v1_arm_freeze_candidate_audit.tsv`
  - `drug_class_profile_v1_summary.tsv`
  - `run_log.json`
- Notebook-regenerated review tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `131_drug_class_profile_v1_group_summary.tsv`
  - `132_drug_class_profile_v1_confounder_coverage.tsv`
  - `133_drug_class_profile_v1_value_composition.tsv`
  - `134_drug_class_profile_v1_arm_freeze_candidate_audit.tsv`
  - `135_drug_class_profile_v1_summary.tsv`

## Patient-Level Drug-Class Profile Buckets

Fill from:

- `drug_class_profile_v1_group_summary.tsv`
- `131_drug_class_profile_v1_group_summary.tsv`
- `135_drug_class_profile_v1_summary.tsv`

Bucket definitions (from saved patient_drug_class_profile_v1.tsv only — no raw re-derivation):

- `single_known_class_clean`: single_or_multi_drug_class == single_class AND both manual-review flags == no
- `single_known_class_plus_unknown`: single_or_multi_drug_class == single_class AND at least one review flag == yes
- `multi_known_class`: single_or_multi_drug_class == multi_class
- `unknown_only`: single_or_multi_drug_class == unknown_only AND no missing-like flag
- `missing_like_only`: single_or_multi_drug_class == unknown_only AND contains_missing_like_raw_drug_name flag set

Patient counts by bucket:

- Total treated patients reviewed: `[fill in total — should be 780]`
- `single_known_class_clean`: `[fill in patient_count and patient_fraction]`
- `single_known_class_plus_unknown`: `[fill in patient_count and patient_fraction]`
- `multi_known_class`: `[fill in patient_count and patient_fraction]`
- `unknown_only`: `[fill in patient_count and patient_fraction]`
- `missing_like_only`: `[fill in patient_count and patient_fraction]`
- Bucket count sum check: `[fill in sum — should equal 780]`

## Patient Counts By Provisional Drug Class

Fill from:

- `135_drug_class_profile_v1_summary.tsv` (classes section)

For any-class membership (class appears in patient's provisional_drug_class_values_json):

- `alkylating_agent`: `[fill in count]`
- `anthracycline`: `[fill in count]`
- `taxane`: `[fill in count]`
- `platinum`: `[fill in count]`
- `antimetabolite`: `[fill in count]`
- `endocrine_serm`: `[fill in count]`
- `endocrine_aromatase_inhibitor`: `[fill in count]`
- `endocrine_other`: `[fill in count]`
- `her2_targeted`: `[fill in count]`
- `immunotherapy`: `[fill in count]`
- `ancillary_supportive`: `[fill in count]`
- `other_cytotoxic`: `[fill in count]`

For dominant-class assignment (dominant_provisional_drug_class_if_any):

- `endocrine_serm`: `[fill in dominant count]`
- `endocrine_aromatase_inhibitor`: `[fill in dominant count]`
- Other dominant classes with N >= 50: `[fill in note]`
- Major classes for audit (N >= 20 any-membership): `[fill in list]`
- Classes with dominant N >= 50: `[fill in count and names]`

## OS Event Distribution By Bucket

Fill from:

- `drug_class_profile_v1_group_summary.tsv`
- `131_drug_class_profile_v1_group_summary.tsv`

- `single_known_class_clean`: `[fill in os_event_count and os_event_fraction]`
- `single_known_class_plus_unknown`: `[fill in os_event_count and os_event_fraction]`
- `multi_known_class`: `[fill in os_event_count and os_event_fraction]`
- `unknown_only`: `[fill in os_event_count and os_event_fraction]`
- `missing_like_only`: `[fill in os_event_count and os_event_fraction]`
- Main descriptive observation about OS distribution across buckets: `[fill in short note]`
- Note on adequacy of OS event counts for any future arm-freeze analysis: `[fill in note]`

## Radiation, Regimen Context, And Timing Coverage By Bucket

Fill from:

- `drug_class_profile_v1_group_summary.tsv`
- `131_drug_class_profile_v1_group_summary.tsv`

- `single_known_class_clean`: radiation=`[fill]`, regimen=`[fill]`, timing=`[fill]`
- `single_known_class_plus_unknown`: radiation=`[fill]`, regimen=`[fill]`, timing=`[fill]`
- `multi_known_class`: radiation=`[fill]`, regimen=`[fill]`, timing=`[fill]`
- `unknown_only`: radiation=`[fill]`, regimen=`[fill]`, timing=`[fill]`
- `missing_like_only`: radiation=`[fill]`, regimen=`[fill]`, timing=`[fill]`
- Buckets with the weakest treatment-context detail: `[fill in note]`

## Baseline Confounder Coverage By Bucket And Major Class

Fill from:

- `drug_class_profile_v1_confounder_coverage.tsv`
- `132_drug_class_profile_v1_confounder_coverage.tsv`
- `135_drug_class_profile_v1_summary.tsv`

Descriptive non-missing coverage threshold used: `0.85`

For single_known_class_clean bucket (primary candidate subset):

- `age_at_diagnosis` coverage: `[fill in non_missing_count/row_count and fraction]`
- `er_status_by_ihc` coverage: `[fill in non_missing_count/row_count and fraction]`
- `pr_status_by_ihc` coverage: `[fill in non_missing_count/row_count and fraction]`
- `her2_status_by_ihc` coverage: `[fill in non_missing_count/row_count and fraction]`
- `ajcc_pathologic_tumor_stage` coverage: `[fill in non_missing_count/row_count and fraction]`
- `histological_type` coverage: `[fill in non_missing_count/row_count and fraction]`

For all_treated (overall):

- `age_at_diagnosis` coverage: `[fill in]`
- `er_status_by_ihc` coverage: `[fill in]`
- `pr_status_by_ihc` coverage: `[fill in]`
- `her2_status_by_ihc` coverage: `[fill in]`
- `ajcc_pathologic_tumor_stage` coverage: `[fill in]`
- `histological_type` coverage: `[fill in]`

For major provisional drug classes (any-class membership):

- `endocrine_serm` (any-class): ER/PR/HER2/stage coverage: `[fill in note]`
- `endocrine_aromatase_inhibitor` (any-class): ER/PR/HER2/stage coverage: `[fill in note]`
- Other major classes noted: `[fill in note]`

Key coverage observations:

- Buckets meeting all 4 coverage-ready fields >= 85%: `[fill in note]`
- Notable coverage weakness still relevant to next step: `[fill in note]`

## Composition Summaries By Bucket

Fill from:

- `drug_class_profile_v1_value_composition.tsv`
- `133_drug_class_profile_v1_value_composition.tsv`

For `single_known_class_clean` bucket:

- ER status distribution: `[fill in top values]`
- PR status distribution: `[fill in top values]`
- HER2 status distribution: `[fill in top values]`
- AJCC tumor stage distribution: `[fill in top values]`
- Histological type distribution: `[fill in top values]`
- Age at diagnosis: min=`[fill]`, median=`[fill]`, max=`[fill]`

For `multi_known_class` bucket (descriptive note only):

- Predominant ER/PR/HER2 profile: `[fill in short note]`
- Age range: `[fill in]`

## Arm-Freeze Candidate Audit

Fill from:

- `drug_class_profile_v1_arm_freeze_candidate_audit.tsv`
- `134_drug_class_profile_v1_arm_freeze_candidate_audit.tsv`

Thresholds used for candidate assessment:

- Minimum patient count for consideration: `50`
- Manual-review fraction tolerance: `0.20`
- Confounder non-missing coverage tolerance: `0.85`

Candidate statuses:

- `promising`: N >= 50, manual_review_fraction <= 0.20, confounder_coverage_ok == yes
- `needs_review`: N >= 50 but fails manual-review or coverage threshold
- `not_recommended`: N < 50 or fails all thresholds

Promising candidates (`arm_freeze_candidate_status == promising`):

- `[fill in candidate_label 1]`: n=`[fill]`, os_events=`[fill]`, manual_review_count=`[fill]`, confounder_coverage_ok=`[fill]`
  - Why: `[fill in why_candidate_or_not]`
  - Recommended: `[fill in recommended_next_handling]`
- `[fill in candidate_label 2]`: n=`[fill]`, os_events=`[fill]`, manual_review_count=`[fill]`, confounder_coverage_ok=`[fill]`
  - Why: `[fill in why_candidate_or_not]`
  - Recommended: `[fill in recommended_next_handling]`
- `[fill in additional promising candidates if any]`

Needs-review candidates (`arm_freeze_candidate_status == needs_review`):

- `[fill in candidate_label]`: n=`[fill]`, why=`[fill]`, what blocks: `[fill in]`

Total candidates assessed: `[fill in]`
Total promising: `[fill in]`
Total needs_review: `[fill in]`
Total not_recommended: `[fill in]`

## What Still Blocks Treatment-Arm Freeze Review

- Multi-class patient majority: `[fill in confirmation that 553/780 are multi_class and not yet resolved]`
- Remaining manual-review burden in the full cohort: `[fill in total manual_review_count and fraction]`
- Unknown-only and missing-like-only patients remain unresolved: `[fill in counts]`
- Whether single_known_class_clean subset alone is sufficient for a useful freeze: `[fill in note]`
- Whether dominant-class subsets have adequate OS event counts for analysis: `[fill in note]`
- Final blocker statement: `[fill in concise sentence — e.g., treatment-arm freeze review requires resolving the dominant-class distribution within the clean bucket and confirming OS event adequacy before proceeding]`

## Validation Checks Completed

- [ ] Latest pointer resolves to the reviewed run.
- [ ] `run_log.json` reports `status == completed`.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] Required upstream pointers were found.
- [ ] Required source tables were found.
- [ ] Drug-name normalization, patient treatment grouping, and baseline analysis run logs all report completion.
- [ ] Drug-class profile row count is positive.
- [ ] Grouping and baseline analysis row counts are positive.
- [ ] All drug-class-profile barcodes resolved in grouping and baseline tables.
- [ ] Bucket counts sum to the total treated patient count.
- [ ] Group summary has exactly 5 rows (one per bucket).
- [ ] Confounder coverage row count is positive.
- [ ] Value composition row count is positive.
- [ ] Arm-freeze candidate audit row count is positive.
- [ ] Summary readiness interpretation is one of the allowed descriptive statuses.
- [ ] No silent overwrite of existing run directory.
- [ ] Latest pointer updated after all outputs succeeded.
- [ ] Notebook review tables `131` through `135` were regenerated from saved outputs on disk only.
- [ ] This note remains a descriptive review and arm-freeze candidate audit only and does not present frozen treatment arms, perform modeling, or produce treatment recommendations.

## Current Interpretation

This note should answer one question only: whether the current TCGA-BRCA drug-class profile layer
provides a strong enough descriptive basis to justify beginning treatment-arm freeze candidate discussion.
It must stay explicit that completing this review does not itself freeze treatment arms, does not resolve
the multi-class patient majority, does not complete the manual-review burden, and does not authorize
treatment recommendation modeling.

The promising candidates identified here — if confirmed during review — would indicate that a subset of
patients with a single clean drug-class assignment and adequate confounder coverage exists and could
support a future arm-freeze decision. However, any actual freeze must be preceded by explicit arm-freeze
review, OS event adequacy assessment, and resolution of remaining manual-review flags.
