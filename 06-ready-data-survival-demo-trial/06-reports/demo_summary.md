# Demo Summary

Generated: 2026-04-01T21:09:26+00:00

## Scope

- This package is a Stage 1 prognostic survival demo.
- It uses public ready-made non-image data only.
- It does not estimate treatment effects or recommend treatments.

## Cohorts

- TCGA-BRCA: 1076 patients with 1076 analyzable OS records.
- METABRIC: 2509 patients with 1981 analyzable OS records.

## Subtype Choices

- TCGA-BRCA subtype field used for plots: `PAM50Call_RNAseq` from UCSC Xena TCGA hub clinical matrix.
- METABRIC subtype field used for plots: `CLAUDIN_SUBTYPE` from cBioPortal patient clinical data.

## Kaplan-Meier Plots

- TCGA-BRCA plotted subtype groups: LumA, LumB, Basal, Her2, Normal
- METABRIC plotted subtype groups: LumA, LumB, Her2, claudin-low, Basal, Normal

## Baseline Cox Models

- TCGA-BRCA: complete cases=805, concordance_index=0.729, formula=`age_at_diagnosis_years + C(ajcc_stage_simple) + C(subtype_label)`
- METABRIC: complete cases=1980, concordance_index=0.652, formula=`age_at_diagnosis_years + npi + C(subtype_label)`

## Meeting-Safe Interpretation

- These models are prognostic association models within each cohort.
- Treatment variables are retained only as observational context and availability summaries.
- No causal treatment recommendation claim is made anywhere in this package.
