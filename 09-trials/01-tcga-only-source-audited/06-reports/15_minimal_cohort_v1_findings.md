# TCGA-BRCA Minimal Cohort V1 Findings

This note records human-reviewed findings from the TCGA-BRCA minimal cohort v1 workflow.

Important reminders:

- this document reviews a provisional minimal cohort v1 only
- this document does not create a final frozen cohort
- this document does not freeze the final endpoint
- this document does not replace source files, parsed TSVs, prior audit TSVs, blueprint TSVs, ambiguity-resolution TSVs, dry-run TSVs, cohort v1 TSVs, or run logs
- this document covers only the saved pointer-driven minimal cohort v1 build produced from currently audited TCGA-only evidence
- this document does not include new raw parsing, XML parsing, SSF parsing, METABRIC work, treatment aggregation, final cohort harmonization, endpoint freeze, or modeling

## Reviewed V1 Build ID

- Minimal cohort v1 latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_minimal_cohort_v1_latest.json`
- Reviewed cohort v1 build id: `[fill in cohort v1 build id]`
- Source minimal dry-run build id: `[fill in dry-run build id]`
- Source cohort blueprint run id: `[fill in blueprint run id]`
- Source ambiguity-resolution run id: `[fill in ambiguity-resolution run id]`
- Source clinical shortlist run id: `[fill in shortlist run id]`
- Source endpoint crosswalk run id: `[fill in endpoint crosswalk run id]`
- Source biospecimen identifier crosswalk run id: `[fill in biospecimen crosswalk run id]`
- Review date: `[fill in review date]`
- Reviewed by: `[fill in name or initials]`

## Cohort V1 Inputs Used

- Minimal dry-run latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_minimal_dry_run_cohort_latest.json`
- Cohort blueprint latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_cohort_blueprint_latest.json`
- Ambiguity-resolution latest pointer: `01-data/audit/tcga-brca/cohort/tcga_brca_blueprint_ambiguity_resolution_latest.json`
- Clinical shortlist latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_clinical_shortlist_latest.json`
- Endpoint crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_endpoint_crosswalk_latest.json`
- Biospecimen identifier crosswalk latest pointer: `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_identifier_crosswalk_latest.json`
- Processed cohort TSV: `01-data/processed/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1.tsv`
- Audit spec TSV: `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_spec.tsv`
- Audit join audit TSV: `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_join_audit.tsv`
- Audit exclusions TSV: `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_exclusions.tsv`
- Audit summary TSV: `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_summary.tsv`
- Audit run log: `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/run_log.json`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `67_minimal_cohort_v1_preview.tsv`
  - `68_minimal_cohort_v1_spec.tsv`
  - `69_minimal_cohort_v1_join_audit.tsv`
  - `70_minimal_cohort_v1_exclusions.tsv`
  - `71_minimal_cohort_v1_ambiguity_flags.tsv`
  - `72_minimal_cohort_v1_summary.tsv`

## Final V1 Row Count

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/72_minimal_cohort_v1_summary.tsv`

Record only saved cohort v1 statements.

- `[fill in final row count]`
- `[fill in patient/case unit statement]`
- `[fill in provisional clinical_patient universe statement]`

## Baseline Fields Included

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_spec.tsv`
- `01-data/processed/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/68_minimal_cohort_v1_spec.tsv`

Record the saved baseline field set exactly as included in cohort v1.

- `[fill in baseline field set summary]`
- `[fill in bcr_patient_barcode note]`
- `[fill in bcr_patient_uuid note]`

## Endpoint Candidate Fields Included

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_spec.tsv`
- `01-data/processed/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/68_minimal_cohort_v1_spec.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/69_minimal_cohort_v1_join_audit.tsv`

Be explicit that these remain source-specific candidate columns only and are not a frozen endpoint.

- `[fill in patient-layer candidate columns]`
- `[fill in follow-up-layer candidate columns]`
- `[fill in side-by-side endpoint note]`

## Biospecimen Sample Anchor Fields Included

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_spec.tsv`
- `01-data/processed/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/68_minimal_cohort_v1_spec.tsv`

Record only the saved sample-anchor evidence behavior. Do not extend this section into child-layer expansion.

- `[fill in biospecimen_sample UUID-anchor statement]`
- `[fill in grouped sample barcode evidence statement]`
- `[fill in no portion/analyte/slide/aliquot expansion statement]`

## Exclusions

Fill from:

- `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_exclusions.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/70_minimal_cohort_v1_exclusions.tsv`

Separate explicit v1 exclusions from unresolved signals that remain audit notes only.

- `[fill in sparse timing field exclusions]`
- `[fill in treatment exclusions]`
- `[fill in biospecimen sample-helper exclusions]`
- `[fill in excluded child-layer biospecimen tables]`

## Ambiguities Still Carried Forward

Fill from:

- `01-data/processed/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1.tsv`
- `01-data/audit/tcga-brca/cohort/minimal_cohort_v1_runs/<COHORT_V1_BUILD_ID>/minimal_cohort_v1_summary.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/71_minimal_cohort_v1_ambiguity_flags.tsv`
- `09-trials/01-tcga-only-source-audited/05-results/72_minimal_cohort_v1_summary.tsv`

Be explicit that this v1 build carries unresolved ambiguity signals forward instead of resolving them.

- `[fill in retained parallel patient-id note]`
- `[fill in 1097 / 1098 / 1101 mismatch note]`
- `[fill in endpoint overlap / no-freeze note]`

## What V1 Is Good Enough For

Summarize only what this provisional build honestly supports now.

- `[fill in baseline-analysis preparation use]`
- `[fill in identifier/join audit use]`
- `[fill in explicit non-final interpretation]`

## What Still Blocks Final Cohort Freeze

Summarize only what remains unresolved after this cohort v1 build. This section must stay at review level only and must not convert v1 into a final frozen cohort.

- `[fill in unresolved patient-count / identifier blocker]`
- `[fill in unresolved endpoint-freeze blocker]`
- `[fill in unresolved child-layer or treatment-design blocker]`

## Validation Checks Completed

- [ ] Minimal cohort v1 latest pointer resolves to the reviewed run.
- [ ] Reviewed cohort v1 run references the intended dry-run, blueprint, ambiguity-resolution, shortlist, endpoint crosswalk, and biospecimen crosswalk runs.
- [ ] Cohort TSV, spec TSV, join audit TSV, exclusions TSV, summary TSV, and run log all exist and are non-empty.
- [ ] `run_log.json` reports `validation.passed == true`.
- [ ] `run_log.json` reports `no_prior_run_overwrite == true`.
- [ ] `run_log.json` reports `latest_pointer_written_after_success_only == true`.
- [ ] Final cohort v1 row count is greater than zero.
- [ ] Final cohort v1 row count equals the `clinical_patient` row count.
- [ ] Cohort v1 contract still matches `1097`, `619`, `716`, `1097`, `2293`, and `1097 / 1098 / 1101`.
- [ ] Deferred sparse timing fields, treatment exclusions, biospecimen sample-helper exclusions, and excluded child-layer biospecimen tables appear in the saved outputs.
- [ ] Review tables `67` through `72` were regenerated from cohort v1 outputs on disk only.
- [ ] This note remains a provisional cohort v1 review only and does not freeze the final cohort or endpoint.

## Current Interpretation

This note remains a provisional minimal cohort v1 review artifact only. It should answer whether the currently audited evidence can support a real, traceable patient-level cohort table for downstream baseline-analysis preparation without pretending unresolved issues are solved. It must not present cohort v1 as a final cohort freeze, a final endpoint definition, a treatment aggregation layer, a child-layer biospecimen expansion, or a modeling-ready dataset.
