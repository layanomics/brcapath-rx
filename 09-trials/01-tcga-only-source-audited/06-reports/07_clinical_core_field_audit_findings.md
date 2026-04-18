# TCGA-BRCA Clinical Core Field Audit Findings

This note records human-reviewed findings from the TCGA-BRCA clinical core field-audit workflow.

Important reminders:

- this document remains part of source audit
- this document does not freeze the cohort
- this document does not freeze the endpoint
- this document does not replace source files, parsed TSVs, schema TSVs, audit TSVs, or run logs
- this document covers only the 4 core parsed clinical biotab tables audited in this workflow
- this document does not include XML parsing, cohort construction, treatment recommendation logic, or modeling

## Reviewed Audit Run ID

- Clinical core field-audit latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_core_field_audit_latest.json`
- Reviewed audit run id: `[fill in audit run id]`
- Upstream clinical biotab parse run id: `[fill in parse run id]`
- Source supplement run id: `[fill in source run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Source Tables Audited

- Upstream clinical biotab latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_biotabs_latest.json`
- Clinical core audit run log: `01-data/audit/tcga-brca/variables/clinical_core_field_audit_runs/<AUDIT_RUN_ID>/run_log.json`
- Combined field audit TSV: `01-data/audit/tcga-brca/variables/clinical_core_field_audit_runs/<AUDIT_RUN_ID>/clinical_core_field_audit.tsv`
- Field-group summary TSV: `01-data/audit/tcga-brca/variables/clinical_core_field_audit_runs/<AUDIT_RUN_ID>/clinical_core_field_group_summary.tsv`
- Missingness summary TSV: `01-data/audit/tcga-brca/variables/clinical_core_field_audit_runs/<AUDIT_RUN_ID>/clinical_core_missingness_summary.tsv`
- Per-table field inventories:
  - `clinical_patient_field_inventory.tsv`
  - `clinical_drug_field_inventory.tsv`
  - `clinical_radiation_field_inventory.tsv`
  - `clinical_follow_up_v4_0_field_inventory.tsv`
- Core source tables reviewed:
  - `clinical_patient.tsv`
  - `clinical_drug.tsv`
  - `clinical_radiation.tsv`
  - `clinical_follow_up_v4_0.tsv`

## Field Counts by Table

Fill from:

- `01-data/audit/tcga-brca/variables/clinical_core_field_audit_runs/<AUDIT_RUN_ID>/clinical_core_field_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/18_clinical_core_field_counts_by_table.tsv`

Suggested review table:

| table_name | field_count | row_count | notes |
| --- | ---: | ---: | --- |
| `[fill in]` | `[fill in]` | `[fill in]` | `[fill in]` |

## Major Field Groups Observed

Fill from:

- `01-data/audit/tcga-brca/variables/clinical_core_field_audit_runs/<AUDIT_RUN_ID>/clinical_core_field_group_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/20_clinical_core_field_group_distribution.tsv`

Record only audit-layer grouping based on exact field names and source-table context.

- demographics:
  - `[fill in observed field-name patterns]`
- diagnosis / pathology:
  - `[fill in observed field-name patterns]`
- stage:
  - `[fill in observed field-name patterns]`
- receptor / biomarker:
  - `[fill in observed field-name patterns]`
- treatment / drug:
  - `[fill in observed field-name patterns]`
- radiation:
  - `[fill in observed field-name patterns]`
- follow-up / outcome-like:
  - `[fill in observed field-name patterns]`
- identifier / admin:
  - `[fill in observed field-name patterns]`
- other / unclear:
  - `[fill in observed field-name patterns]`

## Likely High-Value Fields

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/21_clinical_core_high_completeness_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/23_clinical_core_likely_high_value_fields.tsv`

List exact field names only unless a source-backed interpretation has already been validated.

- `[fill in exact field names]`
- `[fill in exact field names]`
- `[fill in exact field names]`

## Likely Weak or Sparse Fields

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/22_clinical_core_high_missingness_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/24_clinical_core_likely_weak_fields.tsv`

Record only source-table observations such as repeated missing-like tokens, near-empty fields, or very low distinct-value structure.

- `[fill in exact field names or direct audit observations]`
- `[fill in exact field names or direct audit observations]`
- `[fill in exact field names or direct audit observations]`

## Treatment-Field Observations

Use only exact field names and saved audit evidence from `clinical_drug` and treatment-related fields in the other audited tables.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Follow-Up-Field Observations

Use only exact field names and saved audit evidence from `clinical_follow_up_v4_0` and follow-up / outcome-like fields in the other audited tables.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Unresolved Source Questions

- Which core clinical fields appear important enough to cross-check later against source supplements or XML without freezing a cohort now?
- Are any fields grouped as `other / unclear` actually high-priority source-review targets despite ambiguous names?
- Which treatment and follow-up fields need later provenance review before they can be considered for cohort or endpoint decisions?
- Are there fields with strong names but weak completeness that should be explicitly tracked as likely unusable?

## Validation Checks Completed

- [ ] Clinical core field-audit latest pointer resolves to the reviewed audit run.
- [ ] The reviewed audit run references the intended upstream clinical biotab parse run.
- [ ] All 4 per-table field inventory TSVs open and their row counts match upstream core-table column counts.
- [ ] Combined field audit TSV row count matches the sum of the 4 audited table column counts.
- [ ] Missingness summary TSV row count matches the combined field audit TSV row count.
- [ ] Review tables `18` through `24` were regenerated from saved audit outputs on disk only.
- [ ] No raw downloads, parsed clinical TSVs, XML parsing, cohort freeze, endpoint freeze, or modeling were performed in this note.

## Current Interpretation

This note remains a source-audit artifact only. It should summarize which exact source fields exist in the 4 core clinical tables, how those fields cluster at a light audit level, and how usable they appear based on direct missingness-like and table-structure evidence. It must not freeze the cohort, freeze the endpoint, or promote any field into final modeling use without later provenance review.
