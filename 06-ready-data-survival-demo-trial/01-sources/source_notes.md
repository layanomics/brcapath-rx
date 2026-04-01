# Source Notes

Generated: 2026-04-01T21:08:25+00:00

## Scope

- Fast Stage 1 prognostic survival demo only.
- No treatment-effect estimation.
- No treatment recommendation claim.
- No image or WSI data.
- GDC XML is not the main path in this folder.

## Primary Sources Used Tonight

1. TCGA-BRCA from UCSC Xena GDC hub clinical and survival tables.
2. METABRIC from cBioPortal clinical endpoints.

## Auxiliary Source Used

- UCSC Xena TCGA hub `TCGA.BRCA.sampleMap/BRCA_clinicalMatrix`.
- Reason: the GDC-hub TCGA clinical table does not expose PAM50/receptor fields needed for subtype summaries.
- Restriction: this auxiliary table is used only to augment TCGA subtype and receptor annotations.

## Downloaded Files

### tcga_brca_xena_gdc_clinical

- Cohort: TCGA-BRCA
- Platform: UCSC Xena GDC hub
- Role: primary
- Purpose: TCGA ready-made clinical descriptors and observational treatment availability
- URL: https://gdc.xenahubs.net/download/TCGA-BRCA.clinical.tsv.gz
- Local file: `02-raw-downloads/tcga_brca_xena_gdc_clinical.tsv.gz`
- Retrieved at: 2026-04-01T21:08:03+00:00
- Records: 1255
- SHA256: `a510e8013897ac75dde1e16f6fa400bf39629377c7b50c6766bfa347ca888783`
- Notes: Primary TCGA clinical table selected for tonight's Stage 1 non-XML workflow.

### tcga_brca_xena_gdc_survival

- Cohort: TCGA-BRCA
- Platform: UCSC Xena GDC hub
- Role: primary
- Purpose: TCGA overall survival endpoint
- URL: https://gdc.xenahubs.net/download/TCGA-BRCA.survival.tsv.gz
- Local file: `02-raw-downloads/tcga_brca_xena_gdc_survival.tsv.gz`
- Retrieved at: 2026-04-01T21:08:05+00:00
- Records: 1232
- SHA256: `63b073096c9a0a0ae7b332acc14ae20e7e53a544a57648a993b6cbacd5f13046`
- Notes: Primary TCGA survival table selected for tonight's Stage 1 workflow.

### tcga_brca_xena_tcgahub_clinical_matrix

- Cohort: TCGA-BRCA
- Platform: UCSC Xena TCGA hub
- Role: auxiliary
- Purpose: TCGA subtype and receptor augmentation only
- URL: https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix
- Local file: `02-raw-downloads/tcga_brca_xena_tcgahub_clinical_matrix.tsv`
- Retrieved at: 2026-04-01T21:08:09+00:00
- Records: 1247
- SHA256: `39eb3be0fb86e6a577bd2cc01502a7fa5a271e1e1cba294e9dc644ad99580d7f`
- Notes: Used only because PAM50 and receptor fields are absent from the GDC-hub TCGA clinical table.

### metabric_cbio_clinical_attributes

- Cohort: METABRIC
- Platform: cBioPortal REST API
- Role: primary
- Purpose: Clinical field metadata
- URL: https://www.cbioportal.org/api/studies/brca_metabric/clinical-attributes
- Local file: `02-raw-downloads/metabric_cbio_clinical_attributes.json`
- Retrieved at: 2026-04-01T21:08:10+00:00
- Records: 36
- SHA256: `4f0b622daff788340cd621d1d135a239d7dde609b85421feb3c3cee7b7142741`
- Notes: Clinical attribute definitions for METABRIC.

### metabric_cbio_clinical_data_sample

- Cohort: METABRIC
- Platform: cBioPortal REST API
- Role: primary
- Purpose: Sample-level clinical values including stage and receptor fields
- URL: https://www.cbioportal.org/api/studies/brca_metabric/clinical-data
- Local file: `02-raw-downloads/metabric_cbio_clinical_data_sample.json`
- Retrieved at: 2026-04-01T21:08:16+00:00
- Records: 27868
- SHA256: `9886391740cde1f60b3607c049fba79dbdc39079e59a8f54dab424106ab38c6f`
- Notes: Sample-level clinical long table.

### metabric_cbio_clinical_data_patient

- Cohort: METABRIC
- Platform: cBioPortal REST API
- Role: primary
- Purpose: Patient-level clinical values including survival, subtype, and treatment fields
- URL: https://www.cbioportal.org/api/studies/brca_metabric/clinical-data?clinicalDataType=PATIENT
- Local file: `02-raw-downloads/metabric_cbio_clinical_data_patient.json`
- Retrieved at: 2026-04-01T21:08:25+00:00
- Records: 51528
- SHA256: `05e8afd056ecc9cb5049e158c0493d303983371dbee12723fb5d88d57f82addb`
- Notes: Patient-level clinical long table.

## Meeting-Safe Framing

- Treatment fields in these public cohorts are observational descriptors only.
- Treatment fields are retained for availability summaries and context, not for causal claims.
- All prognostic modeling in this trial is cohort-specific and non-causal.
