# TCGA-BRCA Indexed Clinical vs Local Treatment Findings V1

This note is a focused comparison of the current local TCGA-BRCA treatment layers against the official indexed GDC clinical `cases` treatment surface.

It is not a broad source audit.
It does not build treatment arms.
It does not do modeling.
It forces a decision.

## Reviewed Run

- Indexed clinical vs local treatment run id: `20260416T004139Z`
- Review date: `2026-04-16`
- GDC release used for the live indexed pull: `Data Release 45.0 - December 04, 2025`
- GDC tag: `8.3.1`

## Artifacts Created

- Script:
  - `09-trials/01-tcga-only-source-audited/03-scripts/26_compare_tcga_brca_indexed_clinical_vs_local_treatment_v1.py`
- Notebook:
  - `09-trials/01-tcga-only-source-audited/02-notebooks/26_review_tcga_brca_indexed_clinical_vs_local_treatment_v1.ipynb`
- Latest pointer:
  - `01-data/audit/tcga-brca/source/tcga_brca_indexed_clinical_vs_local_treatment_v1_latest.json`
- Run log:
  - `01-data/audit/tcga-brca/source/indexed_clinical_vs_local_treatment_v1_runs/20260416T004139Z/run_log.json`
- Primary run outputs:
  - `01-data/audit/tcga-brca/source/indexed_clinical_vs_local_treatment_v1_runs/20260416T004139Z/indexed_clinical_treatment_patient_level_v1.tsv`
  - `01-data/audit/tcga-brca/source/indexed_clinical_vs_local_treatment_v1_runs/20260416T004139Z/indexed_vs_local_treatment_overlap_v1.tsv`
  - `01-data/audit/tcga-brca/source/indexed_clinical_vs_local_treatment_v1_runs/20260416T004139Z/indexed_vs_local_treatment_gap_resolution_v1.tsv`
  - `01-data/audit/tcga-brca/source/indexed_clinical_vs_local_treatment_v1_runs/20260416T004139Z/indexed_vs_local_treatment_summary_v1.tsv`
- Review-layer exports written by the notebook:
  - `09-trials/01-tcga-only-source-audited/05-results/136_indexed_clinical_treatment_patient_level_v1.tsv`
  - `09-trials/01-tcga-only-source-audited/05-results/137_indexed_vs_local_treatment_overlap_v1.tsv`
  - `09-trials/01-tcga-only-source-audited/05-results/138_indexed_vs_local_treatment_gap_resolution_v1.tsv`
  - `09-trials/01-tcga-only-source-audited/05-results/139_indexed_vs_local_treatment_summary_v1.tsv`

## Direct Answers

### 1. How many patients with no local drug row gain any treatment evidence from indexed clinical?

- `317 of 317`
- That is `100.0%` of the current no-local-drug patients.
- Within those `317` patients:
  - `247` gain potentially useful new evidence because indexed clinical adds broad druglike and/or radiation coverage that is missing locally.
  - `70` gain only coarse recovery, all driven by surgery-only indexed rows.

### 2. How many patients with no local drug-or-radiation rows gain any treatment evidence from indexed clinical?

- `278 of 278`
- That is `100.0%` of the current no-local-drug-or-radiation patients.
- Within those `278` patients:
  - `208` gain potentially useful new evidence.
  - `70` gain only coarse recovery, again driven by surgery-only indexed rows.

### 3. Does indexed clinical add new treatment categories, new patient coverage, or only coarse duplicates?

- It adds new patient coverage.
- It also adds new indexed treatment categories not represented in the current local drug layer, especially:
  - `Surgery, NOS`
  - `Pharmaceutical Therapy, NOS`
  - `Radiation Therapy, NOS`
  - `Radiation, External Beam`
  - smaller counts of `Bisphosphonate Therapy`, `Immunotherapy (Including Vaccines)`, and related radiation subtypes
- But the added detail is mostly broad category-level treatment typing, not rich regimen detail.

Important detail:

- Among the `278` patients with no local drug-or-radiation rows, indexed clinical adds:
  - `208` with broad druglike plus radiation recovery
  - `70` with surgery-only recovery
  - only `1` with any therapeutic-agent detail
  - `0` with regimen-or-line detail
  - `6` with any timing detail

So the value is real, but it is mainly broad coverage recovery, not deep treatment annotation.

### 4. Does indexed clinical reduce the current practical blocker enough to matter?

- Yes.
- The workflow decision rule was:
  - `material_improvement_possible_now` if useful recovery reaches at least `20%` of the no-local-drug group or the no-local-drug-or-radiation group.
- Actual useful recovery was:
  - `247 / 317 = 77.9%` for no-local-drug
  - `208 / 278 = 74.8%` for no-local-drug-or-radiation

That is well above the threshold.

## What The Overlap Actually Looks Like

- Any treatment evidence on the current `1097`-patient local profile cohort:
  - `819` both local and indexed
  - `278` indexed only
  - `0` local only
  - `0` neither
- Druglike coverage:
  - `780` both
  - `247` indexed only
  - `0` local only
  - `70` neither
- Radiation coverage:
  - `528` both
  - `495` indexed only
  - `0` local only
  - `74` neither

Interpretation:

- Indexed clinical is not just repeating the current local drug/radiation coverage.
- It recovers a large block of patients who currently look untreated in the local biotab-derived layers.
- The recovery is still coarse, but it is large enough to change the practical next step.

## Forced Conclusion

`material_improvement_possible_now`

## Recommended Next Step

Integrate indexed patient-level treatment flags and raw indexed treatment-type sets as the next official supplemental treatment layer.
