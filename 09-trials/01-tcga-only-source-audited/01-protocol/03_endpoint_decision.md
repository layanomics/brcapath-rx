# Endpoint Decision Worksheet

## Candidate endpoints
- Overall Survival (OS)
- Disease-Specific Survival (DSS)
- Progression-Free Interval (PFI)
- Recurrence-related endpoint only if source-traceable and sufficiently complete

## Current endpoint position
No endpoint is selected yet.
Endpoint selection must follow TCGA source audit.

## Required checks for each endpoint
- exact source field(s)
- exact event definition
- exact censoring rule
- exact time origin
- completeness / missingness
- biological and clinical plausibility
- suitability for TCGA-only baseline modeling

## Preferred decision order
1. OS
2. DSS
3. PFI
4. recurrence-related endpoint only if auditable

## Decision rule
Only choose an endpoint after source audit confirms:
1. traceable provenance
2. adequate completeness
3. interpretable timing
4. no major ambiguity
