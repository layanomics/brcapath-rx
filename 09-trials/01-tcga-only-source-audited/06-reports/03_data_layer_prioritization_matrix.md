# TCGA-BRCA Data Layer Prioritization Matrix

This matrix is for prioritization only.
It does not freeze final scope or final download order.

| Data layer | Broad availability in inventory | Likely value for TCGA-only baseline study | Access burden | Technical burden | Source-audit value | Suggested near-term priority | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Clinical Supplement | High | High | Lower | Low | High | High | Strong source-audit candidate |
| Biospecimen Supplement | High | High | Lower | Low | High | High | Important for identifier hierarchy and linkage |
| Gene Expression Quantification | Present | High | Moderate | Moderate | Moderate | Medium-High | Strong later candidate if PAM50 remains important |
| Slide Image | Present | Medium-High | Lower | Moderate | Lower | Medium | Valuable, but not required before source audit |
| Annotated Somatic Mutation | Very High | Medium | Higher | Moderate | Lower | Medium-Low | Rich but not first source-audit target |
| Raw Simple Somatic Mutation | Very High | Medium | Higher | Higher | Lower | Low | Too heavy for early phase |
| Copy Number Variation | Very High | Medium | Moderate | Moderate | Lower | Medium-Low | Useful later, not first audit layer |
| DNA Methylation | High | Medium | Moderate | Moderate | Lower | Low-Medium | Later-layer candidate |
| Proteome Profiling | Lower | Medium | Moderate | Moderate | Lower | Low-Medium | Interesting but smaller and later-stage |
| Sequencing Reads | Very High | Low for first phase | High | Very High | Lower | Low | Too heavy for first phase |

## Current interpretation
- The project remains broad.
- Source-audit-first layers are currently the most defensible near-term priority.
- Expression remains a strong candidate for the next layer after source audit.
- Large raw omics classes should not be first just because they are numerous.

## Open decisions
- Should expression be prioritized immediately after source audit?
- Should slide linkage be reviewed before expression?
- Which layers are sufficiently open and interpretable for early-phase review?
