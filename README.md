# BRCAPath-Rx

Multimodal deep learning system for personalized breast cancer
treatment decision support. Integrates clinical, genomic, and
treatment history data from TCGA-BRCA and METABRIC to rank
treatment options by predicted benefit for individual patients.

## Setup

    conda env create -f environment.yml
    conda activate brcapath-rx

## Project Structure

| Directory          | Contents                         | Git tracked? |
|--------------------|----------------------------------|--------------|
| 01-data/raw/       | Raw data from GDC and cBioPortal | No           |
| 01-data/processed/ | Cleaned and merged datasets      | No           |
| 01-data/audit/     | Discovery and audit reports      | No           |
| 02-scripts/        | All executable scripts           | Yes          |
| 03-notebooks/      | EDA and reporting notebooks      | Yes          |
| 04-results/        | Figures and tables               | Yes          |
| 05-models/         | Model architecture code          | Yes          |

## Data Sources

| Dataset   | Source                          | Access         |
|-----------|---------------------------------|----------------|
| TCGA-BRCA | GDC — portal.gdc.cancer.gov     | Open (Level 3) |
| METABRIC  | cBioPortal — brca_metabric      | Open           |

## Phase 1: Data Acquisition

Step 1 — Run discovery (produces inventory, downloads nothing):

    python 02-scripts/discovery/00_discover_sources.py

Step 2 — Review the report:

    01-data/audit/source_discovery_report.md

Step 3 — Download scripts written after reviewing the report.

## Git Policy

Raw data is never committed. Always re-downloaded from canonical
sources using scripts in 02-scripts/download/.

Note: 01-data/audit/ is excluded from Git (generated outputs).
After reviewing the discovery report, copy any report you want
to keep permanently into 02-scripts/discovery/ for Git tracking.
