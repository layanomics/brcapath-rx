# TCGA-BRCA Indexed Clinical Treatment Layer Findings V1

This note documents the first official integration of indexed GDC clinical treatment evidence into the current TCGA-BRCA patient-level treatment layer.

It is integration only.
It does not build treatment arms.
It does not do modeling.
It improves broad treatment availability and status flags.
It does not treat indexed treatment types as full regimen reconstruction.

## Reviewed Run

- Patient treatment master run id: `20260416T014627Z`
- Review date: `2026-04-16`
- Anchored local profile run id: `20260414T200843Z`
- Indexed clinical comparison run id: `20260416T004139Z`
- Source-audit lineage run id: `20260415T233526Z`
- GDC release carried forward from the indexed comparison layer: `Data Release 45.0 - December 04, 2025`
- GDC tag: `8.3.1`

## Artifacts Created

- Script:
  - `09-trials/01-tcga-only-source-audited/03-scripts/27_integrate_tcga_brca_indexed_clinical_treatment_layer_v1.py`
- Notebook:
  - `09-trials/01-tcga-only-source-audited/02-notebooks/27_review_tcga_brca_indexed_clinical_treatment_layer_v1.ipynb`
- Latest pointer:
  - `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_master_v1_latest.json`
- Run log:
  - `01-data/audit/tcga-brca/treatment-prep/patient_treatment_master_v1_runs/20260416T014627Z/run_log.json`
- Primary run outputs:
  - `01-data/processed/tcga-brca/treatment-prep/patient_treatment_master_v1_runs/20260416T014627Z/patient_treatment_master_v1.tsv`
  - `01-data/audit/tcga-brca/treatment-prep/patient_treatment_master_v1_runs/20260416T014627Z/patient_treatment_master_v1_source_audit.tsv`
  - `01-data/audit/tcga-brca/treatment-prep/patient_treatment_master_v1_runs/20260416T014627Z/patient_treatment_master_v1_summary.tsv`
  - `01-data/audit/tcga-brca/treatment-prep/patient_treatment_master_v1_runs/20260416T014627Z/patient_treatment_master_v1_missingness.tsv`
- Review-layer exports written by the notebook:
  - `09-trials/01-tcga-only-source-audited/05-results/140_patient_treatment_master_v1.tsv`
  - `09-trials/01-tcga-only-source-audited/05-results/141_patient_treatment_master_v1_source_audit.tsv`
  - `09-trials/01-tcga-only-source-audited/05-results/142_patient_treatment_master_v1_summary.tsv`
  - `09-trials/01-tcga-only-source-audited/05-results/143_patient_treatment_master_v1_missingness.tsv`

## Direct Answers

### 1. How much patient-level treatment coverage improved

- `any_treatment`: `819 -> 1097` for a gain of `278`
- `druglike`: `780 -> 1027` for a gain of `247`
- `radiation`: `528 -> 1023` for a gain of `495`
- `non_surgical`: `819 -> 1027` for a gain of `208`

Interpretation:

- The integrated layer now closes the broad any-treatment gap across the fixed `1097`-patient profile cohort.
- The biggest incremental gain is radiation coverage.
- These gains are broad source-audited treatment presence flags, not detailed regimen recovery.

### 2. How many patients are now supported by local only, indexed only, or both

- `any_treatment`: `0 local_only / 278 indexed_only / 819 both / 0 neither`
- `druglike`: `0 local_only / 247 indexed_only / 780 both / 70 neither`
- `radiation`: `0 local_only / 495 indexed_only / 528 both / 74 neither`
- `non_surgical`: `0 local_only / 208 indexed_only / 819 both / 70 neither`

Additional integration detail:

- `1` indexed-only case remains outside the anchored local profile cohort and was reported but not integrated.
- `1096` patients now have indexed surgery evidence, but surgery remains an indexed-only broad presence signal in this version.

### 3. What still remains weak

- `70` patients are still surgery-only recoveries without broad non-surgical support.
- `70` patients still lack any broad druglike evidence after integration.
- `74` patients still lack any broad radiation evidence after integration.
- Indexed regimen-line detail remains sparse:
  - `3` patients have any indexed regimen-line value in the integrated cohort.
  - `0` locally missing patients gain regimen-line detail from indexed clinical.
- Indexed incremental detail remains limited:
  - `1` locally missing patient gains therapeutic-agent detail.
  - `6` locally missing patients gain any timing detail.

Interpretation:

- Indexed clinical materially improves coverage, but its added detail remains mostly broad treatment presence and type information.
- It should not be interpreted as full treatment-arm or regimen reconstruction.

### 4. Whether the dataset is now strong enough for the next stage of treatment organization

`ready_for_next_treatment_organization_step_with_broad_source_audited_flags_only`

Interpretation:

- Yes for the next treatment-organization step, if that step uses broad, provenance-aware treatment flags.
- No for final treatment-arm construction, causal analysis, treatment recommendation modeling, or claims of detailed regimen completeness.

## Weaknesses

- Surgery recovery is broad and useful for presence auditing, but not sufficient to describe actionable non-surgical management.
- Indexed treatment-type sets improve status coverage, yet remain coarser than a regimen-ready layer.
- Local-only evidence does not rescue additional patients in the master overlap counts because indexed clinical fully covers the local broad any-treatment cohort.
- The integration remains intentionally conservative and keeps prior source-specific outputs untouched.

## Final Recommendation

Keep `patient_treatment_master_v1` as the next official supplemental treatment layer for TCGA-BRCA patient-level work.

Use it as a broad, source-audited treatment master table for the next organization step.
Do not collapse it into treatment arms in this phase.
Do not present indexed treatment types as regimen reconstruction.
Preserve the local and indexed provenance fields in all downstream treatment organization work.
