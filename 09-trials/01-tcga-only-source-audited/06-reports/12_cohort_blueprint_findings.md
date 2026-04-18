# TCGA-BRCA Cohort Blueprint Findings

This note records human-reviewed findings from the TCGA-BRCA cohort blueprint workflow.

Important reminders:

- this document remains a blueprint and study-design artifact only
- this document does not create a final cohort table
- this document does not freeze the final endpoint
- this document does not harmonize source fields into model-ready variables
- this document does not replace source files, parsed TSVs, prior audit TSVs, blueprint TSVs, or run logs
- this document covers only the saved clinical shortlist, endpoint crosswalk, and biospecimen identifier crosswalk layers plus their blueprint synthesis
- this document does not include new raw parsing, XML parsing, SSF parsing, METABRIC work, cohort construction, or modeling

## Reviewed Blueprint Run ID

- Cohort blueprint latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_cohort_blueprint_latest.json`
- Reviewed blueprint run id: `[fill in blueprint run id]`
- Source clinical shortlist run id: `[fill in shortlist run id]`
- Source clinical core audit run id: `[fill in core audit run id]`
- Source endpoint crosswalk run id: `[fill in endpoint crosswalk run id]`
- Source biospecimen identifier crosswalk run id: `[fill in biospecimen crosswalk run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Blueprint Inputs Used

- Clinical shortlist latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_shortlist_latest.json`
- Endpoint crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_endpoint_crosswalk_latest.json`
- Biospecimen identifier crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_identifier_crosswalk_latest.json`
- Cohort blueprint run log: `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/run_log.json`
- Required fields TSV: `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_required_fields.tsv`
- Optional fields TSV: `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_optional_fields.tsv`
- Deferred fields TSV: `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_deferred_fields.tsv`
- Join path TSV: `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_join_path.tsv`
- Ambiguities TSV: `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_ambiguities.tsv`
- Summary TSV: `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_summary.tsv`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `48_cohort_blueprint_required_fields.tsv`
  - `49_cohort_blueprint_optional_fields.tsv`
  - `50_cohort_blueprint_deferred_fields.tsv`
  - `51_cohort_blueprint_join_path.tsv`
  - `52_cohort_blueprint_ambiguities.tsv`
  - `53_cohort_blueprint_summary.tsv`

## Proposed Unit of Analysis

Fill from:

- `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/53_cohort_blueprint_summary.tsv`

Record only saved blueprint statements.

- `[fill in patient/case analysis-unit statement]`
- `[fill in sample anchor statement]`
- `[fill in any explicit not-final-endpoint note]`

## Required Fields for First Baseline

Fill from:

- `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_required_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/48_cohort_blueprint_required_fields.tsv`

Record exact table names, field names, proposed roles, and rationale summaries only.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Optional Fields for First Baseline

Fill from:

- `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_optional_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/49_cohort_blueprint_optional_fields.tsv`

Keep treatment as treatment-proxy only and endpoint as candidate/join-preparation only.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Deferred Fields

Fill from:

- `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_deferred_fields.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/50_cohort_blueprint_deferred_fields.tsv`

Be explicit when a field is deferred because it belongs to later biospecimen child layers or unresolved endpoint/timing logic.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Proposed Join Path

Fill from:

- `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_join_path.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/51_cohort_blueprint_join_path.tsv`

Describe the saved join steps only. Do not convert them into a final merged cohort here.

- `[fill in case -> clinical path]`
- `[fill in case -> biospecimen sample anchor path]`
- `[fill in later-expansion specimen child-layer path]`

## Strongest Current Evidence

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/51_cohort_blueprint_join_path.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/53_cohort_blueprint_summary.tsv`

Record only blueprint-supported statements about strong current evidence.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Major Unresolved Ambiguities

Fill from:

- `01-data/audit/tcga-brca/cohort/cohort_blueprint_runs/<BLUEPRINT_RUN_ID>/cohort_blueprint_ambiguities.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/52_cohort_blueprint_ambiguities.tsv`

Be explicit when ambiguity remains because of overlapping endpoint fields, sparse timing fields, patient barcode vs UUID handling, treatment-table multiplicity, or missing direct child-layer keys.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Implications for Later Cohort Construction

Summarize only what this blueprint changes for the next cohort-construction step. This section must still stay at planning level and must not create the final cohort.

- `[fill in concise implication]`
- `[fill in concise implication]`
- `[fill in concise implication]`

## Validation Checks Completed

- [ ] Cohort blueprint latest pointer resolves to the reviewed blueprint run.
- [ ] Reviewed blueprint run references the intended shortlist, endpoint crosswalk, and biospecimen crosswalk runs.
- [ ] All six blueprint TSV outputs are present and have row counts greater than zero.
- [ ] Required, optional, and deferred field assignments are disjoint by source layer, table, and field.
- [ ] Join-path rows use only `strong_current_evidence`, `workable_but_needs_review`, or `deferred`.
- [ ] Review tables `48` through `53` were regenerated from blueprint outputs on disk only.
- [ ] The proposed analysis unit remains patient/case only.
- [ ] The required biospecimen anchor remains sample only.
- [ ] Treatment fields remain proxy-only and endpoint fields remain candidate/join-preparation only.
- [ ] No final cohort table, endpoint freeze, harmonized variable table, XML/SSF parsing, METABRIC work, or modeling were performed in this note.

## Current Interpretation

This note remains a cohort-construction blueprint only. It should summarize how the first TCGA-only baseline study is currently best framed from the saved audited evidence: patient/case as the analysis unit, sample as the required biospecimen anchor, clinical baseline fields as the first required covariate layer, treatment represented only as source-audited proxy fields, and endpoint represented only as candidate/join-preparation fields. It must not present a final merged cohort, freeze a final endpoint, or convert the blueprint into a modeling-ready dataset.
