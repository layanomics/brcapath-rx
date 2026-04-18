# TCGA-BRCA Clinical Biotab Findings

This note is the review template for the first TCGA-BRCA Clinical Supplement biotab parsing workflow.

Important reminders:

- this document remains part of source audit
- this document does not freeze the cohort
- this document does not freeze the endpoint
- this document does not replace source files, parsed manifests, schema inventories, or run logs
- this document covers Clinical Supplement biotab text tables only
- this document does not include Clinical XML, OMF XML, or any final cohort construction

## Reviewed Parsing Run ID

- Clinical biotab latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_biotabs_latest.json`
- Reviewed parsing run id: `[fill in parse run id]`
- Source supplement run id: `[fill in source run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Source Files Parsed

- Source supplement latest pointer: `01-data/audit/tcga-brca/source/tcga_brca_source_supplements_latest.json`
- Clinical biotab manifest: `01-data/audit/tcga-brca/variables/clinical_biotab_runs/<PARSE_RUN_ID>/clinical_biotab_table_manifest.tsv`
- Clinical biotab run log: `01-data/audit/tcga-brca/variables/clinical_biotab_runs/<PARSE_RUN_ID>/run_log.json`
- Discovered clinical biotab source files from the reviewed source run:
  - `nationwidechildrens.org_clinical_patient_brca.txt`
  - `nationwidechildrens.org_clinical_drug_brca.txt`
  - `nationwidechildrens.org_clinical_radiation_brca.txt`
  - `nationwidechildrens.org_clinical_follow_up_v1.5_brca.txt`
  - `nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt`
  - `nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt`
  - `nationwidechildrens.org_clinical_nte_brca.txt`
  - `nationwidechildrens.org_clinical_follow_up_v4.0_nte_brca.txt`
  - `nationwidechildrens.org_clinical_omf_v4.0_brca.txt`

## Parsed Tables Generated

- Processed run directory: `01-data/processed/tcga-brca/clinical/biotab_runs/<PARSE_RUN_ID>/`
- Audit run directory: `01-data/audit/tcga-brca/variables/clinical_biotab_runs/<PARSE_RUN_ID>/`
- Parsed primary TSV outputs:
  - `clinical_patient.tsv`
  - `clinical_drug.tsv`
  - `clinical_radiation.tsv`
  - `clinical_follow_up_v1_5.tsv`
  - `clinical_follow_up_v2_1.tsv`
  - `clinical_follow_up_v4_0.tsv`
  - `clinical_nte.tsv`
  - `clinical_follow_up_v4_0_nte.tsv`
  - `clinical_omf_v4_0.tsv`
- Schema TSV outputs:
  - `schema_clinical_patient.tsv`
  - `schema_clinical_drug.tsv`
  - `schema_clinical_radiation.tsv`
  - `schema_clinical_follow_up_v1_5.tsv`
  - `schema_clinical_follow_up_v2_1.tsv`
  - `schema_clinical_follow_up_v4_0.tsv`
  - `schema_clinical_nte.tsv`
  - `schema_clinical_follow_up_v4_0_nte.tsv`
  - `schema_clinical_omf_v4_0.tsv`

## Row Counts

Fill from:

- `01-data/audit/tcga-brca/variables/clinical_biotab_runs/<PARSE_RUN_ID>/clinical_biotab_table_manifest.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/15_clinical_biotab_table_summary.tsv`

Suggested review table:

| table_name | row_count | column_count | source_filename | notes |
| --- | ---: | ---: | --- | --- |
| `[fill in]` | `[fill in]` | `[fill in]` | `[fill in]` | `[fill in]` |

## Major Field Groups Observed

Fill from column-name review only. Do not overinterpret field meaning at this stage.

- demographics:
  - `[fill in observed column-name groups]`
- stage and diagnosis:
  - `[fill in observed column-name groups]`
- receptor and biomarker fields:
  - `[fill in observed column-name groups]`
- drug and treatment fields:
  - `[fill in observed column-name groups]`
- radiation fields:
  - `[fill in observed column-name groups]`
- follow-up and outcome-like fields:
  - `[fill in observed column-name groups]`
- new tumor event / NTE fields:
  - `[fill in observed column-name groups]`

## Likely High-Value Fields for Later Audit

List field names only unless a source-backed interpretation has already been validated.

- `[fill in exact field names from parsed tables]`
- `[fill in exact field names from parsed tables]`
- `[fill in exact field names from parsed tables]`

## Missingness / Structure Observations

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/17_clinical_biotab_missingness_summary.tsv`

Record only direct source-table observations such as:

- repeated `[Not Available]` / `[Not Applicable]` / `[Unknown]` patterns
- differences across follow-up table versions
- table-level row-count differences
- fields that appear present in names but are sparsely populated
- any structural quirks carried in the biotab preamble or schema rows

## Unresolved Parsing Questions

- Are any important clinical variables represented only in Clinical XML and not in the project-level biotab tables?
- How should `clinical_nte` and `clinical_follow_up_v4.0_nte` be cross-checked later against XML without collapsing them now?
- What specific value does `clinical_omf_v4.0` add relative to patient and follow-up biotab tables?
- Are any fields duplicated across patient, follow-up, NTE, and OMF tables in ways that will need later provenance review?

## Validation Checks Completed

- [ ] Latest clinical biotab pointer resolves to the reviewed parsing run.
- [ ] Clinical biotab manifest opens and lists one row per parsed table.
- [ ] Parsed TSV row counts match the manifest row counts.
- [ ] Schema TSV files preserve line-2 aliases and line-3 CDE metadata.
- [ ] Review tables were regenerated from parsed outputs on disk only.
- [ ] No XML parsing, cohort freeze, endpoint freeze, or modeling was performed in this note.

## Current Interpretation

The parsed clinical biotab layer is already informative enough to prioritize a first structured clinical field audit. The main high-value tables are `clinical_patient`, `clinical_drug`, `clinical_radiation`, and `clinical_follow_up_v4_0`. These appear to capture the core patient-level, treatment-level, radiation-level, and follow-up-level source structure needed for the next audit phase. Older follow-up versions, NTE tables, and OMF tables remain important supporting layers but should not drive the first-pass audit alone. The next recommended step is a field-level audit of the Tier 1 tables, with missingness review built in, before moving on to biospecimen biotab parsing.
