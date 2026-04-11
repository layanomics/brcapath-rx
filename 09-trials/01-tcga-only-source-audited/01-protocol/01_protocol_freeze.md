# TCGA-Only Protocol Freeze

## Study type
Prognostic modeling study using TCGA-BRCA only in the initial trial.

## Target population
Patients in TCGA-BRCA eligible for a source-audited cohort build.

## Intended prediction task
To be finalized after source audit.

## Unit of analysis
Choose one and do not mix:
- patient / case
- sample
- slide

## Candidate outcomes
- overall survival
- disease-specific survival
- progression-free interval
- recurrence-related endpoint if truly available and auditable

## Time zero
To be finalized explicitly before modeling.

## Candidate predictor domains
- demographics
- diagnosis / stage / grade
- receptor status
- treatment exposure fields
- biospecimen linkage fields
- omics-derived labels only if provenance is explicit

## Inclusion principles
- must belong to TCGA-BRCA
- must have auditable identifiers
- must have source-traceable variables used in the trial

## Exclusion principles
- duplicate or ambiguous identifiers
- post-baseline leakage variables if prediction is baseline-oriented
- variables with unclear provenance

## Validation gates before modeling
- cohort counts reconciled
- identifier hierarchy reconciled
- source files logged
- variable dictionary started
- endpoint definition frozen
- missingness reviewed
