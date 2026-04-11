# Broad TCGA-BRCA Inventory Interpretation

## Purpose
This note summarizes the full TCGA-BRCA inventory at a broad level.

Important:
- this is still the whole-inventory phase
- this does not freeze the final project scope
- this does not choose the final endpoint
- this does not restrict the project to only clinical and biospecimen data

## Broad inventory summary
The TCGA-BRCA inventory confirms that the project has a large multi-omic and multi-modal data universe rather than a narrow clinical-only dataset.

## Data categories observed
- Simple Nucleotide Variation: 21132
- Copy Number Variation: 14346
- Sequencing Reads: 9282
- Structural Variation: 5772
- Biospecimen: 5317
- Transcriptome Profiling: 4876
- DNA Methylation: 3714
- Somatic Structural Variation: 3128
- Clinical: 2288
- Proteome Profiling: 919

## Data types observed
Top observed file types include:
- Annotated Somatic Mutation: 10134
- Aligned Reads: 9282
- Raw Simple Somatic Mutation: 6751
- Transcript Fusion: 4924
- Structural Rearrangement: 3976
- Gene Level Copy Number: 3314
- Copy Number Segment: 3256
- Slide Image: 3112
- Masked Intensities: 2476
- Raw Intensities: 2263

## Experimental strategies observed
- WXS: 17049
- Genotyping Array: 14329
- WGS: 12383
- RNA-Seq: 11079
- <missing>: 4493
- Methylation Array: 3714
- miRNA-Seq: 3621
- Tissue Slide: 1979
- Diagnostic Slide: 1133
- Reverse Phase Protein Array: 919

## Access overview
- controlled: 42843
- open: 27931

## Sample type overview
- Primary Tumor: 58244 linked files
- Blood Derived Normal: 33721 linked files
- Solid Tissue Normal: 4525 linked files
- Metastatic: 309 linked files

## Workflow overview
Top observed workflow types include:
- DNAcopy: 4458
- BWA with Mark Duplicates and BQSR: 4382
- SeSAMe Methylation Beta Estimation: 3714
- VarScan2 Annotation: 2986
- Arriba: 2462
- STAR - Counts: 2462
- STAR-Fusion: 2462
- BCGSC miRNA Profiling: 2414
- Birdseed: 2263
- ASCAT2: 2168

## Broad interpretation
- TCGA-BRCA clearly contains enough breadth for multiple possible project directions.
- The inventory should remain broad at this stage.
- Clinical Supplement and Biospecimen Supplement remain strong first source-audit candidates, but they do not define the whole project scope.
- Transcriptome, slide image, mutation, copy-number, methylation, and proteome layers all remain in scope for later evidence-based prioritization.
- The main question now is not whether TCGA-BRCA has enough data, but which data layers should be examined first in depth for a defensible TCGA-only baseline study.

## What remains open
- final unit of analysis
- final endpoint
- exact role of PAM50
- which data layer should be examined after source audit
- whether first deep-dive should prioritize source supplements, expression, slide linkage, or another layer

## Current position
At this phase, the project remains an all-available-data inventory and prioritization exercise, not a narrow precommitment to one source class.
