# TCGA Source Acquisition Plan

## Official source priority
1. GDC Data Portal / GDC API
2. TCGA clinical supplement files
3. TCGA biospecimen supplement files
4. GDC-generated manifest files
5. GDC Data Transfer Tool logs

## First download scope
- TCGA-BRCA clinical supplement files
- TCGA-BRCA biospecimen supplement files

## Rationale
- supplement files may contain fields not fully represented in indexed API clinical data
- biospecimen XML may contain linkage details needed for sample/case hierarchy
- raw source files must be preserved unchanged before parsing

## Download policy
- keep all raw files immutable
- save manifest files separately
- save download logs separately
- do not mix parsed outputs with raw source files
