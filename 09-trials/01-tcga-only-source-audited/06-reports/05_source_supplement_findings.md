# TCGA-BRCA Source Supplement Findings

This note records the human-reviewed findings from the Clinical Supplement and Biospecimen Supplement source-acquisition step.

Important reminders:

- this document is still part of source audit
- this document does not freeze the cohort
- this document does not freeze the endpoint
- this document does not replace source files, manifests, or run logs

## Reviewed Run ID

- Reviewed run id: `20260412T000556Z`
- Review date: `2026-04-12`
- Reviewed by: `[fill in name or initials]`

## Source Files Generated

- Latest supplement manifest pointer: `01-data/audit/tcga-brca/source/tcga_brca_source_supplements_latest.json`
- Supplement run log: `01-data/audit/tcga-brca/source/supplements/runs/20260412T000556Z/run_log.json`
- Combined supplement metadata TSV: `01-data/audit/tcga-brca/source/supplements/runs/20260412T000556Z/tcga_brca_source_supplements_metadata.tsv`
- Clinical Supplement manifest: `01-data/raw/tcga-brca/gdc/manifests/clinical/20260412T000556Z/clinical_supplement_manifest.tsv`
- Biospecimen Supplement manifest: `01-data/raw/tcga-brca/gdc/manifests/biospecimen/20260412T000556Z/biospecimen_supplement_manifest.tsv`
- Clinical raw download folder, if used: `01-data/raw/tcga-brca/gdc/downloads/clinical/20260412T000556Z`
- Biospecimen raw download folder, if used: `01-data/raw/tcga-brca/gdc/downloads/biospecimen/20260412T000556Z`
- Review result tables in `09-trials/01-tcga-only-source-audited/05-results/`:
  - `08_source_supplement_counts.tsv`
  - `09_source_supplement_case_coverage.tsv`
  - `10_source_supplement_data_format_counts.tsv`
  - `11_source_supplement_access_state_counts.tsv`
  - `12_source_supplement_case_count_distribution.tsv`
  - `13_source_supplement_filename_patterns.tsv`
  - `14_source_supplement_download_status.tsv`

## Supplement Counts

- Clinical Supplement file count: `1183`
- Biospecimen Supplement file count: `2205`
- Notes on format mix across source classes:
  - Clinical Supplement includes:
    - `1097` BCR XML
    - `77` BCR OMF XML
    - `9` BCR Biotab
  - Biospecimen Supplement includes:
    - `1098` BCR XML
    - `1097` BCR SSF XML
    - `10` BCR Biotab

## Case Coverage Observations

- Clinical Supplement unique case submitter coverage: `1098`
- Biospecimen Supplement unique case submitter coverage: `1098`
- Notes on one-file-per-case versus project-level source tables:
  - Most XML-like source files are effectively one-file-per-case.
  - Clinical distribution shows:
    - `1174` files with `case_count = 1`
    - `9` files with `case_count = 1098`
  - Biospecimen distribution shows:
    - `2195` files with `case_count = 1`
    - `10` files with `case_count = 1098`
  - This strongly suggests a mixed source structure:
    - per-case XML / XML-like files
    - project-level tabular biotab files spanning the full cohort

## Identifier Observations

- Observed case identifier patterns:
  - Per-case files are consistently named with `TCGA-...` case identifiers embedded in the filename.
  - Project-level biotab files do not represent a single case and instead summarize cohort-wide content.

- Observed filename patterns:
  - Clinical:
    - `nationwidechildrens.org_clinical` (`1097`)
    - `nationwidechildrens.org_omf` (`77`)
    - project-level biotabs such as:
      - `nationwidechildrens.org_clinical_patient_brca.txt`
      - `nationwidechildrens.org_clinical_drug_brca.txt`
      - `nationwidechildrens.org_clinical_radiation_brca.txt`
      - multiple follow-up tables
  - Biospecimen:
    - `nationwidechildrens.org_biospecimen` (`1098`)
    - `nationwidechildrens.org_ssf` (`1097`)
    - project-level biotabs such as:
      - `nationwidechildrens.org_biospecimen_sample_brca.txt`
      - `nationwidechildrens.org_biospecimen_portion_brca.txt`
      - `nationwidechildrens.org_biospecimen_analyte_brca.txt`
      - `nationwidechildrens.org_biospecimen_slide_brca.txt`
      - `nationwidechildrens.org_biospecimen_aliquot_brca.txt`

- Observed distinctions between XML, SSF XML, OMF XML, and Biotab files:
  - `BCR XML` appears to be the main per-case source format.
  - `BCR OMF XML` appears as an additional clinical source subtype with much smaller count (`77`).
  - `BCR SSF XML` appears as a major biospecimen-related source subtype with near case-level coverage (`1097`).
  - `BCR Biotab` files are project-level tabular summaries and appear to be strong early parsing candidates because they are easier to inspect and align with source audit.

## Likely Next Parsing Targets

- Candidate 1: `Clinical Biotab tables`
  Rationale: Lowest-friction starting point for source audit of patient, drug, radiation, and follow-up fields. These are easier to parse than per-case XML and directly relevant to cohort, treatment, and endpoint feasibility review.

- Candidate 2: `Biospecimen Biotab tables`
  Rationale: Best early target for understanding identifier hierarchy across patient, sample, portion, analyte, slide, and aliquot layers before any deeper linkage work.

- Candidate 3: `Per-case BCR XML / SSF XML validation layer`
  Rationale: Should be used after biotab-first parsing to resolve ambiguities, confirm provenance, and inspect fields not captured clearly enough in the project-level tables.

## Unresolved Source Questions

- Which variables are present only in per-case XML versus also represented in the project-level biotabs?
- What is the exact added value of `BCR OMF XML` for clinical review and `BCR SSF XML` for biospecimen review?
- Should the first parsing workflow begin with project-level biotabs only, or should it include a small XML validation subset from the start?

## Validation Checks Completed

- [x] Latest supplement pointer resolves to the reviewed run.
- [x] Combined supplement metadata TSV opens and contains both source classes.
- [x] Clinical and Biospecimen manifest row counts match metadata row counts.
- [x] All reviewed files remain restricted to `TCGA-BRCA`.
- [x] Review tables were regenerated from disk only.
- [x] No cohort freeze or endpoint freeze decisions were made in this note.

## Current Interpretation

The supplement acquisition step confirms that TCGA-BRCA source supplements are broad, cleanly retrievable, and cohort-complete at the case-submitter level. The source structure is mixed: many per-case XML-style files plus a smaller number of project-level biotab tables that span the full cohort. For the next phase, the most defensible parsing strategy is biotab-first for rapid source audit, followed by selective XML/SSF/OMF review for provenance confirmation and ambiguity resolution.
