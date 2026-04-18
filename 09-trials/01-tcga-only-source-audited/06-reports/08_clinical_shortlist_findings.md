# TCGA-BRCA Clinical Shortlist Findings

This note records human-reviewed findings from the TCGA-BRCA clinical shortlist audit workflow.

Important reminders:

- this document remains part of source audit
- this document does not freeze the cohort
- this document does not freeze the endpoint
- this document does not recommend treatment-effect claims or treatment-response modeling
- this document does not replace source files, parsed TSVs, prior audit TSVs, shortlist TSVs, or run logs
- this document covers only the shortlist layer built from the 4 core clinical tables
- this document does not include XML parsing, biospecimen parsing, cohort construction, or modeling

## Reviewed Shortlist Run ID

- Clinical shortlist latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_shortlist_latest.json`
- Reviewed shortlist run id: `[fill in shortlist run id]`
- Source clinical core audit run id: `[fill in core audit run id]`
- Upstream clinical biotab parse run id: `[fill in parse run id]`
- Source supplement run id: `[fill in source run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Source Audit Run Used

- Clinical core audit latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_core_field_audit_latest.json`
- Shortlist run log: `01-data/audit/tcga-brca/variables/clinical_shortlist_runs/<SHORTLIST_RUN_ID>/run_log.json`
- Shortlist TSV: `01-data/audit/tcga-brca/variables/clinical_shortlist_runs/<SHORTLIST_RUN_ID>/clinical_shortlist.tsv`
- Shortlist summary TSV: `01-data/audit/tcga-brca/variables/clinical_shortlist_runs/<SHORTLIST_RUN_ID>/clinical_shortlist_summary.tsv`
- Shortlist by-table TSV: `01-data/audit/tcga-brca/variables/clinical_shortlist_runs/<SHORTLIST_RUN_ID>/clinical_shortlist_by_table.tsv`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `25_clinical_shortlist_bucket_counts.tsv`
  - `26_clinical_shortlist_bucket_counts_by_table.tsv`
  - `27_clinical_shortlist_usable_baseline_fields.tsv`
  - `28_clinical_shortlist_usable_treatment_proxy_fields.tsv`
  - `29_clinical_shortlist_usable_endpoint_candidate_fields.tsv`
  - `30_clinical_shortlist_weak_or_unusable_fields.tsv`
  - `31_clinical_shortlist_unclear_manual_review_fields.tsv`

## Usable Baseline Fields

Fill from:

- `01-data/audit/tcga-brca/variables/clinical_shortlist_runs/<SHORTLIST_RUN_ID>/clinical_shortlist.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/27_clinical_shortlist_usable_baseline_fields.tsv`

List exact field names and direct shortlist observations only.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Usable Treatment-Proxy Fields

Fill from:

- `01-data/audit/tcga-brca/variables/clinical_shortlist_runs/<SHORTLIST_RUN_ID>/clinical_shortlist.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/28_clinical_shortlist_usable_treatment_proxy_fields.tsv`

Record only shortlist-level observations. Do not claim treatment-effect validity.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Usable Endpoint-Candidate Fields

Fill from:

- `01-data/audit/tcga-brca/variables/clinical_shortlist_runs/<SHORTLIST_RUN_ID>/clinical_shortlist.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/29_clinical_shortlist_usable_endpoint_candidate_fields.tsv`

Record only endpoint-candidate status from source audit. Do not freeze the endpoint here.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Weak / Unusable Fields

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/30_clinical_shortlist_weak_or_unusable_fields.tsv`

Record direct shortlist reasons such as all-missing, extreme missingness, or near-constant structure.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Unclear Fields Needing Manual Review

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/31_clinical_shortlist_unclear_manual_review_fields.tsv`

Record exact field names and why ambiguity remains at the shortlist layer.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Major Next-Step Recommendation

Summarize the next source-audit move only. Keep this out of cohort freeze, endpoint freeze, and modeling.

- `[fill in concise next-step recommendation]`

## Validation Checks Completed

- [ ] Clinical shortlist latest pointer resolves to the reviewed shortlist run.
- [ ] Reviewed shortlist run references the intended clinical core audit run.
- [ ] Shortlist TSV row count matches the combined core field-audit row count.
- [ ] Every field has exactly one shortlist bucket and one manual review priority.
- [ ] Shortlist summary TSV bucket counts sum to the total shortlist row count.
- [ ] Shortlist by-table TSV bucket counts sum to each table's audited field count.
- [ ] Review result tables `25` through `31` were regenerated from shortlist outputs on disk only.
- [ ] No raw downloads, parsed clinical TSVs, XML parsing, biospecimen parsing, cohort freeze, endpoint freeze, or modeling were performed in this note.

## Current Interpretation

This note remains an audit and prioritization artifact only. It should summarize which exact source fields look most promising for baseline review, treatment-proxy review, and endpoint-candidate review, while explicitly preserving ambiguity where source names or completeness are not strong enough to support firmer claims. It must not define the final cohort, freeze the endpoint, or convert shortlist buckets into final standardized variables.
