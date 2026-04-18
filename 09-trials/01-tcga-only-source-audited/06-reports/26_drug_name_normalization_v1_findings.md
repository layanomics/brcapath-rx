# TCGA-BRCA Drug Name Normalization V1 Findings

This note records the first provisional TCGA-BRCA drug-name normalization and patient-level drug-class summary review.

Important reminders:

- this document is a data-organization and normalization-audit artifact only
- this document does not modify raw treatment rows
- this document does not split compound or regimen-like raw drug-name values
- this document does not freeze final treatment arms
- this document does not perform causal analysis, treatment-effect estimation, or modeling
- this document does not produce treatment recommendations
- this document remains TCGA-only and source-audited

## Reviewed Run

- Drug-name normalization latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_drug_name_normalization_v1_latest.json`
- Reviewed drug-name normalization v1 run id: `[fill in drug-name normalization v1 run id]`
- Source patient treatment profile v1 run id: `[fill in patient treatment profile v1 run id]`
- Source patient treatment grouping v1 run id: `[fill in patient treatment grouping v1 run id]`
- Source treatment-OS overlap v1 run id: `[fill in treatment-OS overlap v1 run id]`
- Source OS endpoint v1 run id: `[fill in OS endpoint v1 run id]`
- Source clinical biotab parse run id: `[fill in clinical biotab parse run id]`
- Source clinical source run id: `[fill in clinical source run id]`
- Source baseline model-input v1 run id: `[fill in baseline model-input v1 run id]`
- Source cohort v1 build id: `[fill in cohort v1 build id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Inputs Used

- Drug-name normalization latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_drug_name_normalization_v1_latest.json`
- Patient treatment profile latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_profile_v1_latest.json`
- Patient treatment grouping latest pointer: `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_grouping_v1_latest.json`
- Clinical biotabs latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_biotabs_latest.json`
- Processed patient drug-class profile TSV: `01-data/processed/tcga-brca/treatment-prep/drug_name_normalization_v1_runs/<DRUG_NAME_NORMALIZATION_V1_RUN_ID>/patient_drug_class_profile_v1.tsv`
- Audit outputs:
- `drug_name_inventory_v1.tsv`
- `drug_name_normalization_map_v1.tsv`
- `drug_name_class_map_v1.tsv`
- `drug_name_normalization_v1_summary.tsv`
- `run_log.json`
- Notebook-regenerated review tables in `09-trials/01-tcga-only-source-audited/05-results/`:
- `126_drug_name_inventory_v1.tsv`
- `127_drug_name_normalization_map_v1.tsv`
- `128_drug_name_class_map_v1.tsv`
- `129_patient_drug_class_profile_v1.tsv`
- `130_drug_name_normalization_v1_summary.tsv`

## Distinct Raw Drug-Name Inventory

Fill from:

- `drug_name_inventory_v1.tsv`
- `127_drug_name_normalization_map_v1.tsv`
- `130_drug_name_normalization_v1_summary.tsv`

- Total distinct raw drug names: `[fill in count]`
- Total mapped confidently: `[fill in count]`
- Total mapped but review: `[fill in count]`
- Total unmapped: `[fill in count]`
- Top confident brand/generic or formatting collapses observed: `[fill in note]`
- Top regimen-like, compound, placebo, or other preserved-review labels: `[fill in note]`
- Top still-unmapped raw names: `[fill in note]`

## Provisional Normalized Names And Drug Classes

Fill from:

- `drug_name_normalization_map_v1.tsv`
- `drug_name_class_map_v1.tsv`
- `130_drug_name_normalization_v1_summary.tsv`

- Total distinct provisional normalized names: `[fill in count]`
- Total distinct provisional drug classes observed: `[fill in count]`
- Major provisional classes observed: `[fill in note]`
- Which normalized names remain `unknown_or_review_needed`: `[fill in note]`
- Whether the class vocabulary looks adequate for the next review step: `[fill in note]`

## Patient-Level Drug-Class Burden

Fill from:

- `patient_drug_class_profile_v1.tsv`
- `130_drug_name_normalization_v1_summary.tsv`

- Total treated patients reviewed: `[fill in count]`
- Patients with any normalized name: `[fill in count]`
- Patients with any unknown/review-needed mapping: `[fill in count and fraction]`
- Patients with single known drug class: `[fill in count and fraction]`
- Patients with multiple known drug classes: `[fill in count and fraction]`
- Patients with unknown-only drug class profile: `[fill in count and fraction]`
- Patients with missing-like-only raw drug names: `[fill in count]`
- Main patient-level review flags still present: `[fill in note]`

## Whether This Reduces The Treatment-Arm Freeze Blocker

This section must stay explicit that readiness here is still provisional review readiness only, not arm finalization.

- Saved readiness interpretation: `[fill in value from 130 summary]`
- Whether the repo is now ready for drug-class review: `[fill in yes/no with note]`
- Whether the repo is now ready for treatment-arm freeze review next: `[fill in yes/no with note]`
- Main reason for the current readiness state: `[fill in concise note]`

## What Still Remains Unresolved

- Remaining unmapped raw drug names: `[fill in note]`
- Remaining review-needed normalized names: `[fill in note]`
- Remaining `unknown_or_review_needed` class burden: `[fill in note]`
- Remaining mixed-patient burden after this layer: `[fill in note]`
- Final blocker statement: `[fill in concise sentence]`

## Validation Checks Completed

- [ ] Latest pointer resolves to the reviewed run.
- [ ] `run_log.json` reports `status == completed`.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] Required upstream pointers were found.
- [ ] Required source tables were found.
- [ ] Patient treatment profile, patient treatment grouping, and clinical biotabs run logs all report completion.
- [ ] `clinical_drug.tsv` row count is positive.
- [ ] Treated-patient count is positive.
- [ ] Every distinct raw drug name appears once in `drug_name_inventory_v1.tsv`.
- [ ] Every distinct raw drug name appears once in `drug_name_normalization_map_v1.tsv`.
- [ ] Inventory and normalization-map raw-name coverage match exactly.
- [ ] Every provisional normalized name appears once in `drug_name_class_map_v1.tsv`.
- [ ] Class-map normalized-name coverage matches the normalization map exactly.
- [ ] Treated-patient set matches the patient set in `clinical_drug.tsv`.
- [ ] `patient_drug_class_profile_v1.tsv` row count matches the treated-patient set.
- [ ] No duplicate patient barcodes exist in `patient_drug_class_profile_v1.tsv`.
- [ ] Patient-profile row order is preserved from the treated subset of `patient_treatment_profile_v1.tsv`.
- [ ] Summary counts reconcile to the inventory, class map, and patient profile tables.
- [ ] Normalization-confidence values are restricted to the allowed set.
- [ ] Mapping-status values are restricted to the allowed set.
- [ ] Class values are restricted to the allowed provisional class vocabulary.
- [ ] Summary readiness interpretation is one of the allowed values.
- [ ] Review tables `126` through `130` were regenerated from saved outputs on disk only.
- [ ] This note remains a provisional normalization/mapping artifact only and does not present frozen treatment arms or treatment recommendations.

## Current Interpretation

This note should answer one question only: whether the current TCGA-BRCA treatment-prep layer now has an auditable provisional drug-name normalization and patient-level drug-class summary strong enough to support the next review step. It must stay explicit that this layer does not replace raw treatment tables, does not silently remove ambiguity, does not finalize treatment arms, and does not authorize treatment recommendation modeling.
