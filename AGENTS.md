# AGENTS

BRCAPath-Rx is currently in a TCGA-only, source-audited restart phase.

- Work inventory-first: scripts generate data, notebooks review data, markdown records decisions.
- Reuse existing config and payload files before adding new settings.
- Keep raw source outputs immutable. New retrievals should create new run folders instead of overwriting prior raw outputs.
- Do not silently invent file names, columns, endpoint semantics, or scientific facts.
- Keep changes easy for a human researcher to validate step by step.
- No modeling, cohort construction, endpoint selection, or TCGA+METABRIC merging during this phase.
