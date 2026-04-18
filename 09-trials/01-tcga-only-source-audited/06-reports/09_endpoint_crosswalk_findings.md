# TCGA-BRCA Endpoint Crosswalk Findings

This note records human-reviewed findings from the TCGA-BRCA endpoint crosswalk audit workflow.

Important reminders:

- this document remains part of source audit
- this document does not freeze the endpoint
- this document does not freeze the cohort
- this document does not define a patient-level cohort
- this document does not replace source files, parsed TSVs, prior audit TSVs, shortlist TSVs, crosswalk TSVs, or run logs
- this document covers only endpoint-candidate reconciliation preparation across `clinical_patient` and `clinical_follow_up_v4_0`
- this document does not include XML parsing, biospecimen parsing, METABRIC work, modeling, or treatment-effect claims

## Reviewed Endpoint Crosswalk Run ID

- Endpoint crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_endpoint_crosswalk_latest.json`
- Reviewed endpoint crosswalk run id: `[fill in crosswalk run id]`
- Source shortlist run id: `[fill in shortlist run id]`
- Source core audit run id: `[fill in core audit run id]`
- Upstream clinical biotab parse run id: `[fill in parse run id]`
- Source supplement run id: `[fill in source run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Source Audit Runs Used

- Clinical shortlist latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_shortlist_latest.json`
- Clinical core audit latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_core_field_audit_latest.json`
- Endpoint crosswalk run log: `01-data/audit/tcga-brca/variables/endpoint_crosswalk_runs/<CROSSWALK_RUN_ID>/run_log.json`
- Endpoint candidate inventory TSV: `01-data/audit/tcga-brca/variables/endpoint_crosswalk_runs/<CROSSWALK_RUN_ID>/endpoint_candidate_inventory.tsv`
- Endpoint crosswalk TSV: `01-data/audit/tcga-brca/variables/endpoint_crosswalk_runs/<CROSSWALK_RUN_ID>/endpoint_crosswalk.tsv`
- Endpoint crosswalk summary TSV: `01-data/audit/tcga-brca/variables/endpoint_crosswalk_runs/<CROSSWALK_RUN_ID>/endpoint_crosswalk_summary.tsv`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `32_endpoint_candidate_inventory.tsv`
  - `33_endpoint_crosswalk_by_family.tsv`
  - `34_endpoint_crosswalk_primary_candidates.tsv`
  - `35_endpoint_crosswalk_overlapping_candidates.tsv`
  - `36_endpoint_crosswalk_ambiguous_candidates.tsv`

## Endpoint-Candidate Fields Found

Fill from:

- `01-data/audit/tcga-brca/variables/endpoint_crosswalk_runs/<CROSSWALK_RUN_ID>/endpoint_candidate_inventory.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/32_endpoint_candidate_inventory.tsv`

Record exact field names, table names, candidate source labels, and direct completeness observations only.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Signal Families Observed

Fill from:

- `01-data/audit/tcga-brca/variables/endpoint_crosswalk_runs/<CROSSWALK_RUN_ID>/endpoint_crosswalk_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/33_endpoint_crosswalk_by_family.tsv`

Do not claim semantic identity beyond what field names and table context support.

- survival_status_like:
  - `[fill in exact observed fields]`
- last_contact_like:
  - `[fill in exact observed fields]`
- death_time_like:
  - `[fill in exact observed fields]`
- progression_or_tumor_status_like:
  - `[fill in exact observed fields]`
- new_tumor_event_like:
  - `[fill in exact observed fields]`
- followup_loss_like:
  - `[fill in exact observed fields]`
- other_endpoint_like:
  - `[fill in exact observed fields if present]`

## Likely Strongest Candidate Fields

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/34_endpoint_crosswalk_primary_candidates.tsv`

Record only crosswalk-layer observations. This is not endpoint freeze.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Overlapping / Conflicting Candidates

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/35_endpoint_crosswalk_overlapping_candidates.tsv`
- `01-data/audit/tcga-brca/variables/endpoint_crosswalk_runs/<CROSSWALK_RUN_ID>/endpoint_crosswalk.tsv`

Describe where the same field name or closely aligned source signal appears across patient and follow-up tables. Do not claim equivalence beyond the saved crosswalk rule and source naming context.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Ambiguous Candidates Needing Manual Review

Fill from:

- `09-trials/01-tcga-only-source-audited/05-results/36_endpoint_crosswalk_ambiguous_candidates.tsv`

Be explicit when sparse `death_days_to`, sparse `new_tumor_event_dx_indicator`, or all-missing progression-time fields were carried forward only as reconciliation targets.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Implications for Later Endpoint Freeze

Summarize only what this crosswalk changes for later endpoint review sequencing. This section must stay at audit/reconciliation level and must not define the final endpoint event, censoring rule, or cohort.

- `[fill in concise implication]`
- `[fill in concise implication]`

## Validation Checks Completed

- [ ] Endpoint crosswalk latest pointer resolves to the reviewed crosswalk run.
- [ ] Reviewed crosswalk run references the intended shortlist run and core audit run.
- [ ] Upstream shortlist and core audit run logs are marked completed.
- [ ] Parsed `clinical_patient.tsv` and `clinical_follow_up_v4_0.tsv` headers were validated against selected crosswalk fields.
- [ ] Endpoint candidate inventory row count is greater than zero and includes at least one `usable_endpoint_candidate`.
- [ ] Endpoint crosswalk TSV row count is greater than zero.
- [ ] Endpoint crosswalk summary TSV family counts sum to the crosswalk TSV row count.
- [ ] Review tables `32` through `36` were regenerated from crosswalk outputs on disk only.
- [ ] Duplicated fields across `clinical_patient` and `clinical_follow_up_v4_0` were not treated as semantically identical by default.
- [ ] `death_days_to` and progression-time fields, when present, were carried forward only as reconciliation-preparation fields rather than endpoint freeze decisions.
- [ ] No raw downloads, XML parsing, biospecimen parsing, METABRIC work, cohort freeze, endpoint freeze, or modeling were performed in this note.

## Current Interpretation

This note remains a source-audit and reconciliation-preparation artifact only. It should summarize which endpoint-like source fields are present across `clinical_patient` and `clinical_follow_up_v4_0`, where those fields overlap by exact name or apparent signal family, how complete they look from saved audit outputs, and which fields still require manual reconciliation before any later endpoint freeze. It must not define the final endpoint, censoring rule, time origin, or final analytic cohort.
