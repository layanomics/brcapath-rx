# TCGA-BRCA Near-Term Action Plan

This plan translates the broad inventory into near-term actions.
It does not freeze the final project scope.

## Current position
The project remains a broad TCGA-BRCA inventory and prioritization effort.
The goal now is to choose the most defensible first deep-dive steps for a TCGA-only baseline workflow.

## Ranked near-term priorities

### Priority 1 — Clinical Supplement source audit
Why this is first:
- directly relevant to cohort definition
- directly relevant to treatment-field review
- directly relevant to endpoint feasibility review
- relatively low technical burden
- high source-audit value

What this step should answer:
- what clinical variables truly exist in source form
- what treatment-related fields are present
- what outcome-related fields are present
- how complete and interpretable these fields are

### Priority 2 — Biospecimen Supplement source audit
Why this is second:
- needed for identifier hierarchy validation
- needed for case/sample/specimen linkage understanding
- important before deeper multi-layer integration
- relatively low technical burden
- high source-audit value

What this step should answer:
- how case, sample, portion, analyte, and slide identifiers relate
- whether linkage structure is clean enough for later integration work
- what biospecimen-derived constraints may affect cohort construction

### Priority 3 — Expression-layer review
Why this is third:
- strong candidate for later biological integration
- especially relevant if PAM50 remains important
- more technically involved than source supplements
- should follow source audit rather than replace it

What this step should answer:
- how many expression quantification files are available
- which workflows are present
- whether expression is the most sensible next layer after source audit

## What remains intentionally deferred
- final endpoint choice
- final unit of analysis
- cohort freeze
- survival modeling
- causal treatment-effect claims
- merged TCGA+METABRIC design
- heavy raw omics download classes
- whole-project narrowing to a single data layer

## Working interpretation
The broad inventory shows that TCGA-BRCA has enough depth for several possible directions.
The best next move is not to chase the largest data classes first.
The best next move is to begin with high-value source-audit layers, then decide which biological layer deserves the next deep review.

## Immediate next implementation target
The next concrete work item should be:
1. Clinical Supplement source acquisition workflow
2. Biospecimen Supplement source acquisition workflow
3. review of expression-layer readiness after source audit starts
