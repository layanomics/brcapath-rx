# TCGA-BRCA Blueprint Ambiguity Resolution Findings

This note records human-reviewed findings from the TCGA-BRCA blueprint ambiguity-resolution workflow.

Important reminders:

- this document remains part of source audit and ambiguity resolution only
- this document does not create a final cohort table
- this document does not freeze the final endpoint
- this document does not replace source files, parsed TSVs, prior audit TSVs, blueprint TSVs, ambiguity-resolution TSVs, or run logs
- this document covers only saved ambiguity triage from the current cohort blueprint plus the current clinical shortlist, endpoint crosswalk, and biospecimen identifier crosswalk layers
- this document does not include raw parsing, XML parsing, SSF parsing, METABRIC work, cohort execution, endpoint freeze, or modeling

## Reviewed Ambiguity-Resolution Run ID

- Ambiguity-resolution latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_blueprint_ambiguity_resolution_latest.json`
- Reviewed ambiguity-resolution run id: `[fill in ambiguity-resolution run id]`
- Source cohort blueprint run id: `[fill in blueprint run id]`
- Source clinical shortlist run id: `[fill in shortlist run id]`
- Source endpoint crosswalk run id: `[fill in endpoint crosswalk run id]`
- Source biospecimen identifier crosswalk run id: `[fill in biospecimen crosswalk run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Ambiguity-Resolution Inputs Used

- Cohort blueprint latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_cohort_blueprint_latest.json`
- Clinical shortlist latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_shortlist_latest.json`
- Endpoint crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_endpoint_crosswalk_latest.json`
- Biospecimen identifier crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_identifier_crosswalk_latest.json`
- Ambiguity-resolution run log: `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/run_log.json`
- Inventory TSV: `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_inventory.tsv`
- Priority TSV: `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_priority.tsv`
- Actions TSV: `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_actions.tsv`
- Summary TSV: `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_summary.tsv`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `54_blueprint_ambiguity_resolution_inventory.tsv`
  - `55_blueprint_blocking_ambiguities.tsv`
  - `56_blueprint_provisional_rule_candidates.tsv`
  - `57_blueprint_manual_review_ambiguities.tsv`
  - `58_blueprint_later_xml_validation_candidates.tsv`
  - `59_blueprint_ambiguity_resolution_actions.tsv`
  - `60_blueprint_ambiguity_resolution_summary.tsv`

## Number of Blocking Ambiguities

Fill from:

- `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/60_blueprint_ambiguity_resolution_summary.tsv`

Record the current `blocking_now`, `review_before_build`, and `non_blocking_for_minimal_build` counts exactly as saved.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Provisional Rules That Appear Safe Enough for a Minimal Build

Fill from:

- `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_inventory.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/56_blueprint_provisional_rule_candidates.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/59_blueprint_ambiguity_resolution_actions.tsv`

Record only saved provisional-rule statements and their evidence basis. Do not convert them into a final cohort execution rule here.

- `[fill in provisional patient/case join rule]`
- `[fill in provisional endpoint coexistence rule]`
- `[fill in any other saved provisional rule]`

## Ambiguities That Still Require Manual Review

Fill from:

- `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_inventory.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/57_blueprint_manual_review_ambiguities.tsv`

Be explicit about what evidence is still missing and why the ambiguity remains unresolved.

- `[fill in]`
- `[fill in]`
- `[fill in]`

## Ambiguities Deferred From the Minimal Build

Fill from:

- `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_inventory.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/58_blueprint_later_xml_validation_candidates.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/59_blueprint_ambiguity_resolution_actions.tsv`

Separate explicit minimal-build deferrals from later XML-validation candidates.

- `[fill in sparse timing deferral]`
- `[fill in treatment-detail deferral]`
- `[fill in later XML-validation candidate]`

## Recommended Actions Before Minimal Cohort Construction

Fill from:

- `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_priority.tsv`
- `01-data/audit/tcga-brca/cohort/ambiguity_resolution_runs/<AMBIGUITY_RESOLUTION_RUN_ID>/ambiguity_resolution_actions.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/59_blueprint_ambiguity_resolution_actions.tsv`

List the grouped actions in saved priority order only.

- `[fill in priority action 1]`
- `[fill in priority action 2]`
- `[fill in priority action 3]`

## What Still Blocks Final Cohort Freeze

Summarize only what remains unresolved after this ambiguity-resolution layer. Keep this section at planning level only and do not present a final cohort or endpoint freeze.

- `[fill in case-count / identifier blocker if still present]`
- `[fill in endpoint-freeze blocker if still present]`
- `[fill in any later child-layer blocker if still relevant]`

## Validation Checks Completed

- [ ] Ambiguity-resolution latest pointer resolves to the reviewed run.
- [ ] Reviewed ambiguity-resolution run references the intended cohort blueprint, clinical shortlist, endpoint crosswalk, and biospecimen crosswalk runs.
- [ ] Inventory, priority, actions, and summary TSV row counts are all greater than zero.
- [ ] Inventory row count equals the current blueprint ambiguity count.
- [ ] Every ambiguity maps to exactly one linked action id.
- [ ] Linked action lists cover each ambiguity id exactly once.
- [ ] Summary counts reconcile to the saved inventory/action tables.
- [ ] Review tables `54` through `60` were regenerated from ambiguity-resolution outputs on disk only.
- [ ] No final cohort table, endpoint freeze, XML/SSF parsing, METABRIC work, or modeling were performed in this note.

## Current Interpretation

This note remains an ambiguity-resolution and resolution-planning artifact only. It should summarize what can be provisionally resolved now, what still requires manual review, what can be deferred from a lean minimal build, and whether the saved evidence suggests that a minimal dry-run cohort build is justified after the documented pre-build actions are addressed. It must not create a final merged cohort, freeze a final endpoint, or convert the audited evidence into a modeling-ready dataset.
