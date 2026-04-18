# TCGA-BRCA Minimal Dry-Run Cohort Findings

This note records human-reviewed findings from the TCGA-BRCA minimal dry-run cohort workflow.

Important reminders:

- this document remains a minimal dry-run cohort review artifact only
- this document does not create a final frozen cohort
- this document does not freeze the final endpoint
- this document does not replace source files, parsed TSVs, prior audit TSVs, blueprint TSVs, ambiguity-resolution TSVs, dry-run TSVs, or run logs
- this document covers only the saved pointer-driven dry-run cohort build produced from currently audited TCGA-only evidence
- this document does not include new raw parsing, XML parsing, SSF parsing, METABRIC work, treatment aggregation, final cohort harmonization, or modeling

## Reviewed Dry-Run Build ID

- Minimal dry-run latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_minimal_dry_run_cohort_latest.json`
- Reviewed dry-run build id: `[fill in dry-run build id]`
- Source cohort blueprint run id: `[fill in blueprint run id]`
- Source ambiguity-resolution run id: `[fill in ambiguity-resolution run id]`
- Source clinical shortlist run id: `[fill in shortlist run id]`
- Source endpoint crosswalk run id: `[fill in endpoint crosswalk run id]`
- Source biospecimen identifier crosswalk run id: `[fill in biospecimen crosswalk run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Dry-Run Inputs Used

- Cohort blueprint latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_cohort_blueprint_latest.json`
- Ambiguity-resolution latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_blueprint_ambiguity_resolution_latest.json`
- Clinical shortlist latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_shortlist_latest.json`
- Endpoint crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_endpoint_crosswalk_latest.json`
- Biospecimen identifier crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_identifier_crosswalk_latest.json`
- Minimal dry-run run log: `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/run_log.json`
- Minimal dry-run cohort TSV: `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_cohort.tsv`
- Minimal dry-run join audit TSV: `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_join_audit.tsv`
- Minimal dry-run row counts TSV: `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_row_counts.tsv`
- Minimal dry-run exclusions TSV: `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_exclusions.tsv`
- Minimal dry-run summary TSV: `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_summary.tsv`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `61_minimal_dry_run_cohort_preview.tsv`
  - `62_minimal_dry_run_join_audit.tsv`
  - `63_minimal_dry_run_row_counts.tsv`
  - `64_minimal_dry_run_exclusions.tsv`
  - `65_minimal_dry_run_ambiguity_flags.tsv`
  - `66_minimal_dry_run_summary.tsv`

## Provisional Patient Universe Used

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/66_minimal_dry_run_summary.tsv`

Record only saved dry-run statements.

- `[fill in patient/case unit-of-analysis statement]`
- `[fill in clinical_patient provisional-universe statement]`
- `[fill in source-row-order / provisional row-id statement]`

## Required Fields Included

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_cohort.tsv`
- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/61_minimal_dry_run_cohort_preview.tsv`

Record the minimal baseline field set exactly as included in the dry run. Keep both patient identifiers explicit and do not rewrite them into one canonical id.

- `[fill in baseline field set summary]`
- `[fill in bcr_patient_barcode note]`
- `[fill in bcr_patient_uuid note]`

## Endpoint Candidate Fields Attached

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_cohort.tsv`
- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_join_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/62_minimal_dry_run_join_audit.tsv`

Be explicit that these remain source-specific candidate columns only and are not a frozen endpoint.

- `[fill in patient-layer candidate columns]`
- `[fill in follow-up-layer candidate columns]`
- `[fill in side-by-side endpoint note]`

## Biospecimen Sample Anchor Attached

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_cohort.tsv`
- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_join_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/62_minimal_dry_run_join_audit.tsv`

Record only the saved sample-anchor evidence behavior. Do not extend this section into child-layer expansion.

- `[fill in biospecimen_sample UUID-anchor statement]`
- `[fill in grouped sample barcode evidence statement]`
- `[fill in no portion/analyte/slide/aliquot expansion statement]`

## Join Success Summary

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_join_audit.tsv`
- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_row_counts.tsv`
- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/62_minimal_dry_run_join_audit.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/63_minimal_dry_run_row_counts.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/66_minimal_dry_run_summary.tsv`

Record the saved counts exactly and keep mismatches visible.

- `[fill in final dry-run row count]`
- `[fill in matched follow-up patient and row counts]`
- `[fill in matched biospecimen sample patient and row counts]`

## Exclusions Made

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_exclusions.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/64_minimal_dry_run_exclusions.tsv`

Separate explicit dry-run exclusions from unresolved signals that remain only as audit notes.

- `[fill in sparse timing field exclusions]`
- `[fill in treatment-detail exclusions]`
- `[fill in excluded child-layer biospecimen tables]`

## Ambiguities Still Carried Forward

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_cohort.tsv`
- `01-data/audit/tcga-brca/cohort/minimal_dry_run_runs/<DRY_RUN_BUILD_ID>/minimal_dry_run_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/65_minimal_dry_run_ambiguity_flags.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/66_minimal_dry_run_summary.tsv`

Be explicit that this minimal dry run carries unresolved ambiguity signals forward instead of resolving them.

- `[fill in retained parallel patient-id note]`
- `[fill in 1097 / 1098 / 1101 mismatch note]`
- `[fill in endpoint overlap / no-freeze note]`

## What Still Blocks Final Cohort Freeze

Summarize only what remains unresolved after this minimal dry-run build. This section must stay at review/planning level only and must not convert the dry run into a final frozen cohort.

- `[fill in unresolved patient-count / identifier blocker]`
- `[fill in unresolved endpoint-freeze blocker]`
- `[fill in unresolved child-layer or treatment-design blocker]`

## Validation Checks Completed

- [ ] Minimal dry-run latest pointer resolves to the reviewed run.
- [ ] Reviewed dry-run run references the intended blueprint, ambiguity-resolution, shortlist, endpoint crosswalk, and biospecimen crosswalk runs.
- [ ] Dry-run cohort, join audit, row counts, exclusions, summary, and run log files all exist and are non-empty.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] Final dry-run cohort row count is greater than zero.
- [ ] Final dry-run cohort row count equals the `clinical_patient` row count.
- [ ] Join-audit counts reconcile with saved row-count and summary tables.
- [ ] The `1097 / 1098 / 1101` mismatch remains explicitly visible in saved outputs.
- [ ] Deferred sparse timing fields, treatment-detail exclusions, and excluded child-layer biospecimen tables appear in the exclusions output.
- [ ] Review tables `61` through `66` were regenerated from dry-run outputs on disk only.
- [ ] This note remains a minimal dry run review only and does not freeze the final cohort or endpoint.

## Current Interpretation

This note remains a minimal dry-run cohort review artifact only. It should answer whether the currently audited evidence can support a traceable patient-level dry-run join test without pretending unresolved issues are solved. It must not present the dry run as a final cohort freeze, a final endpoint definition, a treatment aggregation layer, a child-layer biospecimen expansion, or a modeling-ready dataset.
