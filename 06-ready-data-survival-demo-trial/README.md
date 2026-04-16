# Ready-Data Survival Demo Trial

This folder contains a fast, meeting-safe Stage 1 prognostic survival demo for BRCAPath-Rx.

Scope:
- Prognostic survival/risk modeling only
- No treatment-effect estimation
- No treatment recommendation claims
- No image or WSI data
- Ready-made public non-image data only

Primary sources used here:
- UCSC Xena GDC hub TCGA-BRCA clinical and survival tables
- cBioPortal METABRIC clinical endpoints

Documented auxiliary source:
- UCSC Xena TCGA hub `TCGA.BRCA.sampleMap/BRCA_clinicalMatrix`
- Used only to recover TCGA subtype and receptor fields that are absent from the GDC-hub TCGA clinical table

Run order:

```powershell
& "C:\Users\layan\AppData\Local\Programs\Python\Python311\python.exe" `
  "d:\Projects\brcapath-rx\06-ready-data-survival-demo-trial\01_download_ready_data.py"

& "C:\Users\layan\AppData\Local\Programs\Python\Python311\python.exe" `
  "d:\Projects\brcapath-rx\06-ready-data-survival-demo-trial\02_build_stage1_survival_demo.py"
```

Main outputs:
- `01-sources/dataset_manifest.tsv`
- `01-sources/source_notes.md`
- `03-processed/*.tsv`
- `04-models/*.tsv`
- `05-figures/*.png`
- `06-reports/*.md`
- `07-demo-assets/*.csv`
- `07-demo-assets/demo_metrics.json`
