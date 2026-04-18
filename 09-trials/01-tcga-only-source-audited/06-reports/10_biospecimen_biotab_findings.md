# TCGA-BRCA Biospecimen Biotab Findings

This note is the review template for the first TCGA-BRCA core Biospecimen Supplement biotab parsing workflow.

Important reminders:

- this document remains part of source audit
- this document does not freeze the cohort
- this document does not freeze the endpoint
- this document does not define a final patient-, sample-, or slide-level analysis unit
- this document does not replace source files, parsed TSVs, source inventories, schema inventories, or run logs
- this document covers only core Biospecimen Supplement biotab text tables
- this document does not include Biospecimen XML, SSF XML, METABRIC work, cohort construction, or modeling

## Reviewed Parsing Run ID

- Biospecimen biotab latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_biotabs_latest.json`
- Reviewed parsing run id: `[fill in parse run id]`
- Source supplement run id: `[fill in source run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Source Files Parsed

- Source supplement latest pointer: `01-data/audit/tcga-brca/source/tcga_brca_source_supplements_latest.json`
- Biospecimen biotab source inventory: `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/biospecimen_biotab_source_inventory.tsv`
- Biospecimen biotab manifest: `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/biospecimen_biotab_table_manifest.tsv`
- Biospecimen biotab run log: `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/run_log.json`
- Parsed core biospecimen source files from the reviewed source run:
  - `nationwidechildrens.org_biospecimen_sample_brca.txt`
  - `nationwidechildrens.org_biospecimen_portion_brca.txt`
  - `nationwidechildrens.org_biospecimen_analyte_brca.txt`
  - `nationwidechildrens.org_biospecimen_slide_brca.txt`
  - `nationwidechildrens.org_biospecimen_aliquot_brca.txt`
  - `nationwidechildrens.org_biospecimen_protocol_brca.txt`
  - `nationwidechildrens.org_biospecimen_shipment_portion_brca.txt`
  - `nationwidechildrens.org_biospecimen_diagnostic_slides_brca.txt`
- Discovered but deferred in this v1 workflow:
  - `nationwidechildrens.org_ssf_tumor_samples_brca.txt`
  - `nationwidechildrens.org_ssf_normal_controls_brca.txt`

## Parsed Tables Generated

- Processed run directory: `01-data/processed/tcga-brca/biospecimen/biotab_runs/<PARSE_RUN_ID>/`
- Audit run directory: `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/`
- Parsed primary TSV outputs:
  - `biospecimen_sample.tsv`
  - `biospecimen_portion.tsv`
  - `biospecimen_analyte.tsv`
  - `biospecimen_slide.tsv`
  - `biospecimen_aliquot.tsv`
  - `biospecimen_protocol.tsv`
  - `biospecimen_shipment_portion.tsv`
  - `biospecimen_diagnostic_slides.tsv`
- Schema TSV outputs:
  - `schema_biospecimen_sample.tsv`
  - `schema_biospecimen_portion.tsv`
  - `schema_biospecimen_analyte.tsv`
  - `schema_biospecimen_slide.tsv`
  - `schema_biospecimen_aliquot.tsv`
  - `schema_biospecimen_protocol.tsv`
  - `schema_biospecimen_shipment_portion.tsv`
  - `schema_biospecimen_diagnostic_slides.tsv`

## Row Counts

Fill from:

- `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/biospecimen_biotab_table_manifest.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/37_biospecimen_biotab_table_summary.tsv`

Suggested review table:

| table_name | row_count | column_count | source_filename | notes |
| --- | ---: | ---: | --- | --- |
| `[fill in]` | `[fill in]` | `[fill in]` | `[fill in]` | `[fill in]` |

## Major Identifier Groups Observed

Fill from column-name review only. Do not overinterpret identifier meaning at this stage.

- patient-level identifier columns:
  - `[fill in observed column-name groups]`
- sample-level identifier columns:
  - `[fill in observed column-name groups]`
- portion-level identifier columns:
  - `[fill in observed column-name groups]`
- analyte-level identifier columns:
  - `[fill in observed column-name groups]`
- slide-level identifier columns:
  - `[fill in observed column-name groups]`
- aliquot-level identifier columns:
  - `[fill in observed column-name groups]`
- protocol or shipment-related linkage columns:
  - `[fill in observed column-name groups]`

## Likely Linkage Path Observed

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/40_biospecimen_biotab_linkage_fields.tsv`

Keep this strictly at field-name level and source-audit level.

- likely backbone path:
  - `[fill in field-name-based path such as patient -> sample -> portion / analyte / slide / aliquot]`
- supporting linkage side paths:
  - `[fill in analyte -> protocol or shipment-linked observations if present]`
- note on patient layer:
  - `No standalone biospecimen patient biotab table was parsed in this workflow; patient-level linkage evidence comes from fields such as exact patient barcode or patient UUID column names where present.`

## Likely High-Value Linkage Fields for Later Audit

List field names only unless a source-backed interpretation has already been validated.

- `[fill in exact field names from parsed tables]`
- `[fill in exact field names from parsed tables]`
- `[fill in exact field names from parsed tables]`

## Missingness / Structure Observations

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/39_biospecimen_biotab_missingness_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/41_biospecimen_biotab_preamble_summary.tsv`

Record only direct source-table observations such as:

- repeated `[Not Available]` / `[Not Applicable]` / `[Unknown]` patterns
- row-count differences across sample, portion, analyte, slide, and aliquot tables
- columns that appear identifier-like by name but are incompletely populated
- line-2 `CDE_ID` coverage patterns preserved in schema TSVs
- any table-level preamble or structural quirks documented by the parser

## Unresolved Parsing Questions

- Are any important biospecimen linkage details represented only in Biospecimen XML or SSF XML rather than the core project-level biotab tables?
- How should `biospecimen_shipment_portion` be interpreted relative to portion- and aliquot-level linkage without harmonizing it prematurely?
- What additional linkage value, if any, should be reviewed later in the deferred `ssf_*.txt` text tables?
- Which exact fields will need the next manual identifier-hierarchy audit before any later cohort construction?

## Validation Checks Completed

- [ ] Biospecimen biotab latest pointer resolves to the reviewed parsing run.
- [ ] Biospecimen source inventory opens and lists both parsed core files and deferred `ssf_*.txt` text files.
- [ ] Biospecimen biotab manifest opens and lists one row per parsed core table.
- [ ] Parsed TSV row counts match the manifest row counts.
- [ ] Schema TSV files preserve line-2 `CDE_ID` metadata while primary parsed TSVs preserve line-1 headers exactly.
- [ ] Review tables `37` through `41` were regenerated from parsed outputs on disk only.
- [ ] No Biospecimen XML, SSF XML, cohort freeze, endpoint freeze, identifier harmonization, or modeling was performed in this note.

## Current Interpretation

This note remains a source-audit artifact only. The parsed core biospecimen biotab layer should be used to review which identifier-like fields exist by exact name, how the project-level sample, portion, analyte, slide, aliquot, protocol, shipment, and diagnostic-slide tables appear to relate, and where obvious missingness or structural constraints might affect later linkage audit. It must not define a final harmonized specimen hierarchy, freeze a cohort, or convert these source tables into a modeling-ready analysis dataset.
