# TCGA-BRCA Biospecimen Identifier Crosswalk Findings

This note records human-reviewed findings from the TCGA-BRCA biospecimen identifier crosswalk audit workflow.

Important reminders:

- this document remains part of source audit
- this document does not freeze the cohort
- this document does not freeze the endpoint
- this document does not define a final patient-, sample-, portion-, analyte-, slide-, or aliquot-level cohort
- this document does not replace source files, parsed TSVs, schema TSVs, prior biospecimen review tables, crosswalk TSVs, or run logs
- this document covers only the parsed core biospecimen biotab layer plus approved side-linkage tables
- this document does not include XML parsing, SSF parsing, METABRIC work, cohort construction, endpoint freeze, identifier harmonization, or modeling

## Reviewed Biospecimen Identifier Crosswalk Run ID

- Biospecimen identifier crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_identifier_crosswalk_latest.json`
- Reviewed crosswalk run id: `[fill in crosswalk run id]`
- Source biospecimen parse run id: `[fill in parse run id]`
- Source supplement run id: `[fill in source run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Source Parsing Run Used

- Biospecimen biotab latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_biotabs_latest.json`
- Biospecimen biotab source inventory: `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/biospecimen_biotab_source_inventory.tsv`
- Biospecimen biotab manifest: `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/biospecimen_biotab_table_manifest.tsv`
- Biospecimen biotab run log: `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/<PARSE_RUN_ID>/run_log.json`
- Crosswalk run log: `01-data/audit/tcga-brca/variables/biospecimen_identifier_crosswalk_runs/<CROSSWALK_RUN_ID>/run_log.json`
- Crosswalk inventory TSV: `01-data/audit/tcga-brca/variables/biospecimen_identifier_crosswalk_runs/<CROSSWALK_RUN_ID>/biospecimen_identifier_inventory.tsv`
- Crosswalk TSV: `01-data/audit/tcga-brca/variables/biospecimen_identifier_crosswalk_runs/<CROSSWALK_RUN_ID>/biospecimen_identifier_crosswalk.tsv`
- Crosswalk summary TSV: `01-data/audit/tcga-brca/variables/biospecimen_identifier_crosswalk_runs/<CROSSWALK_RUN_ID>/biospecimen_identifier_crosswalk_summary.tsv`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `42_biospecimen_identifier_inventory.tsv`
  - `43_biospecimen_identifier_crosswalk_by_family.tsv`
  - `44_biospecimen_identifier_primary_candidates.tsv`
  - `45_biospecimen_identifier_overlapping_candidates.tsv`
  - `46_biospecimen_identifier_ambiguous_candidates.tsv`
  - `47_biospecimen_identifier_likely_linkage_chain.tsv`

## Identifier Candidates Found

Fill from:

- `01-data/audit/tcga-brca/variables/biospecimen_identifier_crosswalk_runs/<CROSSWALK_RUN_ID>/biospecimen_identifier_inventory.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/42_biospecimen_identifier_inventory.tsv`

Record exact table names, field names, candidate-source labels, completeness observations, and direct row-level evidence only.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Linkage Families Observed

Fill from:

- `01-data/audit/tcga-brca/variables/biospecimen_identifier_crosswalk_runs/<CROSSWALK_RUN_ID>/biospecimen_identifier_crosswalk_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/43_biospecimen_identifier_crosswalk_by_family.tsv`

Do not claim field equivalence beyond the saved crosswalk role, rule, table context, and reversible pattern checks.

- patient_identifier_like:
  - `[fill in exact observed fields]`
- sample_identifier_like:
  - `[fill in exact observed fields]`
- portion_identifier_like:
  - `[fill in exact observed fields]`
- analyte_identifier_like:
  - `[fill in exact observed fields]`
- slide_identifier_like:
  - `[fill in exact observed fields]`
- aliquot_identifier_like:
  - `[fill in exact observed fields]`
- multi_level_or_unclear_identifier_like:
  - `[fill in exact observed fields if present]`

## Strongest Likely Linkage Keys

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/44_biospecimen_identifier_primary_candidates.tsv`

Keep this at crosswalk-audit level only. This section must not define a final harmonized specimen hierarchy.

- patient:
  - `[fill in strongest barcode and/or uuid candidates]`
- sample:
  - `[fill in strongest barcode and/or uuid candidates]`
- portion:
  - `[fill in strongest barcode and/or uuid candidates]`
- analyte:
  - `[fill in strongest barcode and/or uuid candidates]`
- slide:
  - `[fill in strongest barcode and/or uuid candidates]`
- aliquot:
  - `[fill in strongest barcode and/or uuid candidates]`

## Overlapping or Duplicated Linkage Fields

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/45_biospecimen_identifier_overlapping_candidates.tsv`
- `01-data/audit/tcga-brca/variables/biospecimen_identifier_crosswalk_runs/<CROSSWALK_RUN_ID>/biospecimen_identifier_crosswalk.tsv`

Describe repeated identifier-like fields across tables without claiming identity by default.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Side-Table Linkage Observations

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/43_biospecimen_identifier_crosswalk_by_family.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/47_biospecimen_identifier_likely_linkage_chain.tsv`

Focus on `biospecimen_protocol`, `biospecimen_shipment_portion`, and `biospecimen_diagnostic_slides` only when the saved crosswalk marked them as side-link candidates.

- `[fill in]`
- `[fill in]`

## Ambiguous Linkage Candidates Needing Manual Review

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/46_biospecimen_identifier_ambiguous_candidates.tsv`

Be explicit when a field remained ambiguous because it was generic, multi-level, or only weakly supported by reversible checks.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Likely Central Linkage Chain

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/47_biospecimen_identifier_likely_linkage_chain.tsv`

Keep this section at linkage-preparation level only.

- likely central backbone:
  - `[fill in likely patient -> sample -> portion -> analyte -> slide / aliquot path]`
- likely supporting side links:
  - `[fill in analyte -> protocol, diagnostic slide, or shipment observations if justified by saved evidence]`
- unresolved chain details:
  - `[fill in any steps that remain ambiguous or partial]`

## Validation Checks Completed

- [ ] Biospecimen identifier crosswalk latest pointer resolves to the reviewed crosswalk run.
- [ ] Reviewed crosswalk run references the intended biospecimen parse run.
- [ ] Upstream biospecimen parse run log is marked completed.
- [ ] Required central biospecimen tables were found in the manifest.
- [ ] Inventory, crosswalk, and summary TSV row counts are all greater than zero.
- [ ] Summary family counts sum to the crosswalk TSV row count.
- [ ] Review tables `42` through `47` were regenerated from crosswalk outputs on disk only.
- [ ] Repeated identifier fields across tables were not treated as semantically identical by default.
- [ ] Multi-level or generic identifier-like fields remained ambiguous unless stronger saved evidence justified another role.
- [ ] No raw downloads, XML parsing, SSF parsing, cohort construction, endpoint freeze, METABRIC work, identifier harmonization, or modeling were performed in this note.

## Current Interpretation

This note remains a source-audit and linkage-preparation artifact only. It should summarize which identifier-like fields exist across the parsed biospecimen core tables, how those fields appear to support the patient -> sample -> portion -> analyte -> slide / aliquot linkage chain, where strong barcode- and UUID-based candidates appear to repeat across tables, and which fields still require manual review before any later cohort construction. It must not define a final cohort, freeze a final endpoint, or convert the audited biospecimen source tables into a harmonized analysis dataset.
