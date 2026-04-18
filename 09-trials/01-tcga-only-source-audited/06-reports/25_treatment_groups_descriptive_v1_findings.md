# TCGA-BRCA Treatment Groups Descriptive V1 Findings

This note records the first grouped descriptive review of the provisional TCGA-BRCA patient treatment groups.

Important reminders:

- this document is descriptive review only
- this document does not normalize drug names
- this document does not redefine treatment groups
- this document does not silently exclude mixed patients
- this document does not freeze treatment arms
- this document does not perform causal analysis, treatment-effect estimation, or modeling
- this document does not produce treatment recommendations
- raw treatment tables remain unchanged and are not replaced by this review layer

## Reviewed Run

- Treatment groups descriptive latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_treatment_groups_descriptive_v1_latest.json`
- Reviewed treatment groups descriptive v1 run id: `[fill in treatment groups descriptive v1 run id]`
- Source patient treatment grouping v1 run id: `[fill in patient treatment grouping v1 run id]`
- Source patient treatment profile v1 run id: `[fill in patient treatment profile v1 run id]`
- Source treatment-OS overlap v1 run id: `[fill in treatment-OS overlap v1 run id]`
- Source OS endpoint v1 run id: `[fill in OS endpoint v1 run id]`
- Source baseline model-input v1 run id: `[fill in baseline model-input v1 run id]`
- Source baseline feature-set v1 run id: `[fill in baseline feature-set v1 run id]`
- Source baseline analysis v1 run id: `[fill in baseline analysis v1 run id]`
- Source cohort v1 build id: `[fill in cohort v1 build id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Inputs Used

- Treatment groups descriptive latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_treatment_groups_descriptive_v1_latest.json`
- Patient treatment grouping latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_grouping_v1_latest.json`
- OS endpoint latest pointer: `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_os_endpoint_v1_latest.json`
- Baseline model-input latest pointer: `01-data/audit/tcga-brca/model-input/tcga_brca_baseline_model_input_v1_latest.json`
- Baseline feature-set latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_feature_set_v1_latest.json`
- Baseline analysis latest pointer: `01-data/audit/tcga-brca/analysis-prep/tcga_brca_baseline_analysis_v1_latest.json`
- Audit run directory: `01-data/audit/tcga-brca/treatment-prep/treatment_groups_descriptive_v1_runs/<TREATMENT_GROUPS_DESCRIPTIVE_V1_RUN_ID>/`
- Audit outputs:
- `treatment_groups_descriptive_v1_group_summary.tsv`
- `treatment_groups_descriptive_v1_confounder_coverage.tsv`
- `treatment_groups_descriptive_v1_value_composition.tsv`
- `treatment_groups_descriptive_v1_manual_review_summary.tsv`
- `treatment_groups_descriptive_v1_summary.tsv`
- `run_log.json`
- Notebook-regenerated review tables in `09-trials/01-tcga-only-source-audited/05-results/`:
- `121_treatment_groups_descriptive_v1_group_summary.tsv`
- `122_treatment_groups_descriptive_v1_confounder_coverage.tsv`
- `123_treatment_groups_descriptive_v1_value_composition.tsv`
- `124_treatment_groups_descriptive_v1_manual_review_summary.tsv`
- `125_treatment_groups_descriptive_v1_summary.tsv`

## Group Counts

Fill from:

- `treatment_groups_descriptive_v1_group_summary.tsv`
- `121_treatment_groups_descriptive_v1_group_summary.tsv`

- Total OS-linked cohort size: `[fill in total cohort size]`
- `no_drug_record`: `[fill in patient_count and patient_fraction]`
- `single_chemotherapy`: `[fill in patient_count and patient_fraction]`
- `single_hormone_therapy`: `[fill in patient_count and patient_fraction]`
- `single_targeted_therapy`: `[fill in patient_count and patient_fraction]`
- `single_immunotherapy`: `[fill in patient_count and patient_fraction]`
- `single_ancillary_or_other`: `[fill in patient_count and patient_fraction]`
- `mixed_multi_type`: `[fill in patient_count and patient_fraction]`
- `missing_type_only`: `[fill in patient_count and patient_fraction]`
- Which provisional groups currently have no assigned patients: `[fill in zero-count groups]`

## OS Event Distribution By Group

Fill from:

- `treatment_groups_descriptive_v1_group_summary.tsv`
- `121_treatment_groups_descriptive_v1_group_summary.tsv`

- `no_drug_record`: `[fill in os_event_count and os_event_fraction]`
- `single_chemotherapy`: `[fill in os_event_count and os_event_fraction]`
- `single_hormone_therapy`: `[fill in os_event_count and os_event_fraction]`
- `single_targeted_therapy`: `[fill in os_event_count and os_event_fraction]`
- `single_immunotherapy`: `[fill in os_event_count and os_event_fraction]`
- `single_ancillary_or_other`: `[fill in os_event_count and os_event_fraction]`
- `mixed_multi_type`: `[fill in os_event_count and os_event_fraction]`
- `missing_type_only`: `[fill in os_event_count and os_event_fraction]`
- Main descriptive contrast to note: `[fill in short interpretation]`

## Manual-Review Burden By Group

Fill from:

- `treatment_groups_descriptive_v1_group_summary.tsv`
- `treatment_groups_descriptive_v1_manual_review_summary.tsv`
- `121_treatment_groups_descriptive_v1_group_summary.tsv`
- `124_treatment_groups_descriptive_v1_manual_review_summary.tsv`

- Overall manual-review count and fraction: `[fill in total count and fraction]`
- Whether the current 357-patient burden was preserved: `[fill in yes/no and note]`
- `no_drug_record`: `[fill in manual_review_count and fraction]`
- `single_chemotherapy`: `[fill in manual_review_count and fraction]`
- `single_hormone_therapy`: `[fill in manual_review_count and fraction]`
- `single_targeted_therapy`: `[fill in manual_review_count and fraction]`
- `single_immunotherapy`: `[fill in manual_review_count and fraction]`
- `single_ancillary_or_other`: `[fill in manual_review_count and fraction]`
- `mixed_multi_type`: `[fill in manual_review_count and fraction]`
- `missing_type_only`: `[fill in manual_review_count and fraction]`
- Conflict subtype counts carried forward from grouping v1 conflict audit:
- `[fill in mixed_multi_type_with_provisional_dominant count]`
- `[fill in mixed_multi_type_tied_dominant count]`
- `[fill in missing_type_only count]`

## Radiation, Regimen Context, And Timing Coverage

Fill from:

- `treatment_groups_descriptive_v1_group_summary.tsv`
- `121_treatment_groups_descriptive_v1_group_summary.tsv`

- Radiation overlap by group: `[fill in short table or bullets]`
- Regimen-context coverage by group: `[fill in short table or bullets]`
- Treatment-timing coverage by group: `[fill in short table or bullets]`
- Which groups have the weakest treatment-context detail: `[fill in note]`

## Baseline Confounder Coverage By Group

Fill from:

- `treatment_groups_descriptive_v1_confounder_coverage.tsv`
- `122_treatment_groups_descriptive_v1_confounder_coverage.tsv`
- `treatment_groups_descriptive_v1_summary.tsv`
- `125_treatment_groups_descriptive_v1_summary.tsv`

- Chosen descriptive non-missing coverage threshold: `[fill in threshold]`
- `age_at_diagnosis` coverage by group: `[fill in note]`
- `er_status_by_ihc` coverage by group: `[fill in note]`
- `pr_status_by_ihc` coverage by group: `[fill in note]`
- `her2_status_by_ihc` coverage by group: `[fill in note]`
- `ajcc_pathologic_tumor_stage` coverage by group: `[fill in note]`
- `histological_type` coverage by group: `[fill in note]`
- Groups meeting adequate ER/PR/HER2/stage coverage together: `[fill in count and group names]`
- Any notable coverage weakness that still matters for the next step: `[fill in note]`

## Composition Summaries By Group

Fill from:

- `treatment_groups_descriptive_v1_value_composition.tsv`
- `123_treatment_groups_descriptive_v1_value_composition.tsv`

- ER composition highlights by group: `[fill in note]`
- PR composition highlights by group: `[fill in note]`
- HER2 composition highlights by group: `[fill in note]`
- AJCC pathologic tumor stage composition highlights by group: `[fill in note]`
- Histological type composition highlights by group: `[fill in note]`
- Age-at-diagnosis summaries by group: `[fill in note]`

## Which Groups Look Strongest For The Next Step

- Groups with patient_count >= 50: `[fill in count and group names]`
- Groups with manual_review_fraction within the chosen descriptive tolerance: `[fill in count and group names]`
- Large groups meeting both manual-review and receptor/stage coverage tolerances: `[fill in count and group names]`
- Groups that look descriptively strongest if treatment-arm freeze review were considered next: `[fill in note]`

## What Still Blocks Treatment-Arm Freeze Review

- Drug names remain unnormalized: `[fill in confirmation]`
- Remaining manual-review burden: `[fill in confirmation]`
- Mixed-group handling remains provisional: `[fill in confirmation]`
- Any weaker confounder coverage that still matters: `[fill in confirmation]`
- Final blocker statement: `[fill in concise sentence]`

## Validation Checks Completed

- [ ] Latest pointer resolves to the reviewed run.
- [ ] `run_log.json` reports `status == completed`.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] Required upstream pointers were found.
- [ ] Required source tables were found.
- [ ] Grouping, OS endpoint, baseline model-input, baseline feature-set, and baseline analysis run logs all report completion.
- [ ] Grouping row count matches the OS cohort row count.
- [ ] Row order is preserved from the OS endpoint table.
- [ ] All OS cohort patients are accounted for once via grouping linkage.
- [ ] Grouping lineage matches the OS endpoint and baseline inputs.
- [ ] Group counts reconcile to the patient treatment grouping table.
- [ ] Confounder coverage row counts reconcile to group counts.
- [ ] Manual-review group totals reconcile to the group summary.
- [ ] Manual-review conflict subtype counts sum to the manual-review totals.
- [ ] Summary counts reconcile to the group summary.
- [ ] Value composition includes an age summary row for every non-empty provisional group.
- [ ] The summary readiness interpretation is one of the allowed descriptive statuses.
- [ ] Treatment recommendation modeling remains out of scope.
- [ ] Notebook review tables `121` through `125` were regenerated from saved outputs on disk only.

## Current Interpretation

This note should answer one question only: whether the current provisional treatment groups are descriptively usable enough to justify the next step, treatment-arm freeze review. It must stay explicit that grouped descriptive review completion does not itself freeze treatment arms, does not normalize raw drug names, does not resolve the remaining manual-review burden, and does not authorize treatment recommendation modeling.
