# TCGA-BRCA Inventory Findings

This template records human-reviewed findings from the TCGA-BRCA inventory step.

Important reminders:
- this document is for inventory and audit only
- this document does not freeze a cohort
- this document does not choose an endpoint
- interpretation should be based on saved source outputs and review tables

## Source Files Generated

- Latest manifest: `[fill in path]`
- Run log: `[fill in path]`
- Case inventory TSV: `[fill in path]`
- File inventory TSV: `[fill in path]`
- Raw response snapshot folder: `[fill in path]`
- Review summary tables in `09-trials/01-tcga-only-source-audited/05-results/`: `[fill in filenames used]`

## Case Count

- Case count from saved inventory: `[fill in count]`
- Unique case submitter count: `[fill in count]`
- Notes on any count discrepancies or reruns: `[fill in notes]`

## File Count

- File count from saved inventory: `[fill in count]`
- Notes on pagination completeness: `[fill in notes]`
- Notes on access mix or release-state mix: `[fill in notes]`

## Top Available Data Categories

| Rank | data_category | file_count | Notes |
| --- | --- | --- | --- |
| 1 | `[fill in]` | `[fill in]` | `[fill in]` |
| 2 | `[fill in]` | `[fill in]` | `[fill in]` |
| 3 | `[fill in]` | `[fill in]` | `[fill in]` |

## Top Available Data Types

| Rank | data_type | file_count | Notes |
| --- | --- | --- | --- |
| 1 | `[fill in]` | `[fill in]` | `[fill in]` |
| 2 | `[fill in]` | `[fill in]` | `[fill in]` |
| 3 | `[fill in]` | `[fill in]` | `[fill in]` |

## Availability Observations

- `[fill in observation about broad inventory coverage]`
- `[fill in observation about controlled vs open access]`
- `[fill in observation about sample linkage or workflow coverage]`

## Likely First Download Candidates

- Candidate 1: `[fill in source class or data grouping]`
  Rationale: `[fill in]`
- Candidate 2: `[fill in source class or data grouping]`
  Rationale: `[fill in]`
- Candidate 3: `[fill in source class or data grouping]`
  Rationale: `[fill in]`

## Unresolved Questions

- `[fill in question about source completeness]`
- `[fill in question about file classes worth downloading first]`
- `[fill in question about any identifier or access constraints]`

## Validation Checks Completed

- [ ] Latest manifest points to the reviewed run.
- [ ] Case inventory TSV opens and contains expected columns.
- [ ] File inventory TSV opens and contains expected columns.
- [ ] File inventory row count matches the run log pagination total.
- [ ] Review summary tables were regenerated from disk.
- [ ] No cohort rules or endpoint choices were made in this inventory note.
