# TCGA-BRCA Treatment Source Audit Findings

This note records a combined offline + online audit of whether the local BRCAPath-Rx repo already covers the meaningful GDC treatment-related source layers for `TCGA-BRCA`.

Important reminders:

- this is an audit only
- this does not modify raw source files
- this does not build treatment arms
- this does not perform modeling
- this does not assume missingness just because a layer has not been used yet
- this does not assume XML/BCR solves everything

## Reviewed Audit Run

- Treatment source audit latest pointer: `01-data/audit/tcga-brca/source/tcga_brca_treatment_source_audit_latest.json`
- Reviewed treatment source audit run id: `20260415T233526Z`
- Review date: `2026-04-16`
- Reviewed by: `Codex`

## Audit Artifacts Created

- Audit run log: `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/run_log.json`
- Official GDC source inventory TSV: `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/official_gdc_source_inventory.tsv`
- Local treatment source inventory TSV: `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/local_treatment_source_inventory.tsv`
- Gap analysis TSV: `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/treatment_source_gap_analysis.tsv`
- Audit summary TSV: `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/treatment_source_audit_summary.tsv`
- Official references TSV: `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/treatment_source_official_references.tsv`
- Live GDC query snapshots: `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/`

## What Was Inspected Locally

- Source inventory pointer and full run-specific inventory:
  - `01-data/audit/tcga-brca/source/tcga_brca_inventory_latest.json`
  - `01-data/audit/tcga-brca/source/runs/20260411T224824Z/run_log.json`
  - `01-data/audit/tcga-brca/source/runs/20260411T224824Z/tcga_brca_file_inventory.tsv`
- Source supplement pointer, metadata, manifests, and download logs:
  - `01-data/audit/tcga-brca/source/tcga_brca_source_supplements_latest.json`
  - `01-data/audit/tcga-brca/source/supplements/runs/20260412T000556Z/tcga_brca_source_supplements_metadata.tsv`
  - `01-data/raw/tcga-brca/gdc/manifests/clinical/20260412T000556Z/clinical_supplement_manifest.tsv`
  - `01-data/raw/tcga-brca/gdc/manifests/biospecimen/20260412T000556Z/biospecimen_supplement_manifest.tsv`
  - `01-data/raw/tcga-brca/gdc/downloads/clinical/20260412T000556Z/`
  - `01-data/raw/tcga-brca/gdc/downloads/biospecimen/20260412T000556Z/`
- Parsed clinical and biospecimen biotab outputs:
  - `01-data/audit/tcga-brca/variables/tcga_brca_clinical_biotabs_latest.json`
  - `01-data/audit/tcga-brca/variables/clinical_biotab_runs/20260412T010932Z/clinical_biotab_table_manifest.tsv`
  - `01-data/audit/tcga-brca/variables/tcga_brca_biospecimen_biotabs_latest.json`
  - `01-data/audit/tcga-brca/variables/biospecimen_biotab_runs/20260413T070359Z/biospecimen_biotab_source_inventory.tsv`
- XML-derived endpoint outputs and treatment-prep lineage:
  - `01-data/audit/tcga-brca/endpoint-prep/tcga_brca_endpoint_target_prep_v1_latest.json`
  - `01-data/audit/tcga-brca/endpoint-prep/xml_followup_v1_runs/20260414T025234Z/run_log.json`
  - `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_profile_v1_latest.json`
  - `01-data/audit/tcga-brca/treatment-prep/patient_treatment_profile_v1_runs/20260414T200843Z/run_log.json`
  - `01-data/audit/tcga-brca/treatment-prep/tcga_brca_patient_treatment_grouping_v1_latest.json`
- Protocol note:
  - `09-trials/01-tcga-only-source-audited/01-protocol/04_source_acquisition_plan.md`

Note:
- For full Clinical-category counts I relied on the run-specific inventory referenced by `tcga_brca_inventory_latest.json`. The run log confirms the full `70774`-row inventory for `20260411T224824Z`.

## What Was Checked Online

- Official docs:
  - GDC Clinical Data [R1]
  - GDC Clinical Supplement [R2]
  - GDC Biospecimen Data [R3]
  - GDC Repository / Additional Data Download [R4]
- Official live dictionary endpoints:
  - `clinical_supplement` [R5]
  - `biospecimen_supplement` [R6]
  - `treatment` [R7]
  - `follow_up` [R8]
- Official live API evidence:
  - GDC Files API [R9]
  - GDC Cases API [R10]
  - GDC API status [R11]

## 1. Official GDC Source Inventory

| Official source layer | TCGA-BRCA availability | What it contains in general | Treatment relevance | Relative richness | Key refs |
| --- | --- | --- | --- | --- | --- |
| Indexed GDC clinical surface (`cases` API / portal Clinical TSV-JSON) | Available, `1098` cases | Model-aligned clinical entities, including `diagnoses`, `diagnoses.treatments`, `follow_ups`, `exposures`, and `demographic` | High | Coarser than raw TCGA supplement files; useful official indexed treatment surface | `R1,R4,R7,R8,R10` |
| Clinical Supplement / `BCR XML` | Available, `1097` files | Clinical metadata information in a raw per-case XML subtype | High | Richer raw source than indexed clinical; harder to parse | `R2,R5,R9` |
| Clinical Supplement / `BCR OMF XML` | Available, `77` files | Clinical metadata information in a smaller XML subtype | Possible | Rich raw source; harder to parse | `R2,R5,R9` |
| Clinical Supplement / `BCR Biotab` | Available, `9` files | Clinical metadata information in tabular form | High | Coarser and easier than XML | `R2,R5,R9` |
| Biospecimen Supplement / `BCR XML` | Available, `1098` files | Biospecimen metadata information in raw XML form | Low, indirect | Richer provenance/linkage source than biotab | `R3,R6,R9` |
| Biospecimen Supplement / `BCR SSF XML` | Available, `1097` files | Biospecimen metadata information in SSF XML form | Low, indirect | Richer provenance/linkage source than biotab | `R3,R6,R9` |
| Biospecimen Supplement / `BCR Biotab` | Available, `10` files | Biospecimen metadata information in tabular form | Low, indirect | Coarser and easier than biospecimen XML | `R3,R6,R9` |
| Pathology Report / `PDF` | Available, `1105` files | Unstructured clinical pathology reports | Possible but low for structured extraction | Unstructured narrative layer | `R4,R9` |

Broad-to-narrow discovery note:
- Live category discovery showed `Clinical` contains `Clinical Supplement` and `Pathology Report`, while `Biospecimen` contains `Biospecimen Supplement` and `Slide Image`.
- `Slide Image` is real and official, but not a meaningful treatment source layer for this audit.

Indexed-clinical note:
- The live `cases` audit snapshot confirms that `diagnoses.treatments` is populated in TCGA-BRCA, not just defined in the dictionary. Top `treatment_type` buckets in the saved live query include `surgery, nos (1096)`, `radiation therapy, nos (589)`, `chemotherapy (583)`, `hormone therapy (521)`, and `pharmaceutical therapy, nos (367)`.

## 2. Local Repo Inventory

| Local source layer | Path / scope | What appears to be there | Download status | Parsed status | Treatment value |
| --- | --- | --- | --- | --- | --- |
| Indexed GDC clinical surface | No saved local clinical TSV/JSON export identified | Repo has source inventory/query artifacts, but no dedicated indexed-clinical export artifact | Not downloaded locally | Not applicable locally | High if obtained |
| Clinical Supplement / `BCR XML` | `01-data/raw/tcga-brca/gdc/downloads/clinical/20260412T000556Z/` | `1097` clinical BCR XML files in supplement metadata | Downloaded completely | Only endpoint/follow-up exploitation identified so far | High |
| Clinical Supplement / `BCR OMF XML` | `01-data/raw/tcga-brca/gdc/downloads/clinical/20260412T000556Z/` | `77` clinical BCR OMF XML files in supplement metadata | Downloaded completely | Not yet parsed/exploited for treatment | Possible |
| Clinical Supplement / `BCR Biotab` | Clinical manifest plus `clinical_biotab_table_manifest.tsv` | `9` clinical biotabs downloaded and parsed | Downloaded completely | Already parsed | High |
| Biospecimen Supplement / `BCR XML` | `01-data/raw/tcga-brca/gdc/downloads/biospecimen/20260412T000556Z/` | `1098` biospecimen BCR XML files | Downloaded completely | Not parsed | Low, indirect |
| Biospecimen Supplement / `BCR SSF XML` | `01-data/raw/tcga-brca/gdc/downloads/biospecimen/20260412T000556Z/` | `1097` biospecimen BCR SSF XML files | Downloaded completely | Not parsed | Low, indirect |
| Biospecimen Supplement / `BCR Biotab` | Biospecimen manifest plus source inventory / table manifest | `10` biospecimen biotabs downloaded; `8` parsed, `2` deferred | Downloaded completely | Partially parsed | Low, indirect |
| Pathology Report / `PDF` | Full file inventory says available; raw GDC tree has no PDFs | `1105` pathology reports listed in full file inventory, `0` local PDFs found | Not downloaded locally | Not parsed | Low for structured extraction |

Important local specifics:

- Clinical supplement downloads reconcile locally:
  - `1183` clinical manifest rows
  - `1183` validated clinical downloads
- Biospecimen supplement downloads reconcile locally:
  - `2205` biospecimen manifest rows
  - `2205` validated biospecimen downloads
- Clinical biotab parsing is complete for the `9` clinical biotab files and includes:
  - `clinical_drug`, `clinical_radiation`, `clinical_patient`, `clinical_follow_up_v1_5`, `clinical_follow_up_v2_1`, `clinical_follow_up_v4_0`, `clinical_follow_up_v4_0_nte`, `clinical_nte`, `clinical_omf_v4_0`
- Biospecimen biotab parsing is partial by design:
  - `8` tables parsed
  - `2` deferred SSF text tables: `nationwidechildrens.org_ssf_normal_controls_brca.txt`, `nationwidechildrens.org_ssf_tumor_samples_brca.txt`
- The only identified saved clinical-XML-derived outputs are endpoint/follow-up outputs under:
  - `01-data/processed/tcga-brca/endpoint-prep/xml_followup_v1_runs/20260414T025234Z/`
- The endpoint XML workflow explicitly restricts its XML scope to `BCR XML` and excludes `BCR OMF XML`.
- I found no saved dedicated clinical XML treatment parse output and no saved local indexed-clinical TSV/JSON export.

## 3. Gap Analysis

| Official source layer | Local status | Parsed status | Likely value | Classification | Notes |
| --- | --- | --- | --- | --- | --- |
| Indexed GDC clinical surface | No local export identified | No local exploitation | High | Not downloaded but available in GDC | This is the strongest remaining official treatment surface to check next |
| Clinical Supplement / `BCR XML` | Downloaded completely | Only endpoint/follow-up exploitation identified | High | Already downloaded but not yet fully parsed / not yet fully exploited | Download gap is not the issue here |
| Clinical Supplement / `BCR OMF XML` | Downloaded completely | No dedicated exploitation identified | Possible to moderate | Already downloaded but not yet parsed / not yet exploited | Raw OMF XML still stands as an untested layer |
| Clinical Supplement / `BCR Biotab` | Downloaded completely | Parsed | High | Already downloaded and already parsed | This is the main structured treatment layer currently in use |
| Biospecimen Supplement / `BCR XML` | Downloaded completely | Not parsed | Low, indirect | Already downloaded but not yet parsed / not yet exploited | Low priority for treatment itself |
| Biospecimen Supplement / `BCR SSF XML` | Downloaded completely | Not parsed | Low, indirect | Already downloaded but not yet parsed / not yet exploited | Low priority for treatment itself |
| Biospecimen Supplement / `BCR Biotab` | Downloaded completely | `8/10` parsed | Low, indirect | Already downloaded but not yet fully parsed / not yet fully exploited | Not the main treatment-data gap |
| Pathology Report / `PDF` | Available in GDC but not downloaded locally | Not parsed | Low for structured extraction | Not primary / probably not worth chasing first | Could matter only if you want manual narrative review |

Patient-level missingness interpretation:
- Current treatment-prep missingness is not a supplement download failure.
- `patient_treatment_profile_v1` reports:
  - `780` patients with any drug row
  - `528` patients with any radiation row
  - `278` patients with no drug or radiation rows
- Because indexed clinical and most clinical XML/OMF XML treatment surfaces are still unchecked locally, those patient gaps should not yet be labeled source-level absent across all meaningful GDC layers.

## 4. Interpretation

- You already have the full meaningful supplement download coverage:
  - complete `Clinical Supplement`
  - complete `Biospecimen Supplement`
- You already parsed the most usable structured clinical treatment tables:
  - `clinical_drug`
  - `clinical_radiation`
  - clinical patient/follow-up supporting tables
- The main remaining gap is not another supplement download. The main remaining gap is that two official treatment-adjacent surfaces are still underused:
  - the indexed GDC clinical surface
  - the already downloaded raw clinical XML layers, especially `BCR XML` and `BCR OMF XML`
- Biospecimen XML layers are real and downloaded, but they are mainly provenance/linkage layers rather than primary treatment sources.
- Pathology reports are officially available but are probably not worth chasing first for structured treatment decision-support.
- The repo currently shows a mixed limitation:
  - download coverage looks strong for the meaningful supplement layers
  - source exploitation is still incomplete
  - some patient-level treatment sparsity may still be genuine source-level sparsity, but it is too early to claim that yet

## 5. Final Answers

- What I inspected locally:
  - full run-specific TCGA-BRCA file inventory
  - supplement pointers, metadata, manifests, download logs, raw download directories
  - parsed clinical biotab outputs
  - parsed biospecimen biotab outputs
  - XML endpoint/follow-up outputs
  - treatment-prep lineage outputs and run logs
- What I checked online:
  - official GDC docs for Clinical Data, Clinical Supplement, Biospecimen Data, and Repository Additional Data Download
  - official live GDC data dictionary endpoints for `clinical_supplement`, `biospecimen_supplement`, `treatment`, and `follow_up`
  - live GDC Files API, Cases API, and API status
- What surprised me most:
  - You appear to have already downloaded the complete meaningful supplement file layers, but you still do not have a saved local indexed-clinical export, and the already-downloaded clinical XML layers remain only lightly exploited.
- Whether you likely downloaded the meaningful treatment-related GDC layers already:
  - Yes, for the meaningful structured supplement file layers, very likely yes.
  - The main official treatment-relevant layer still missing locally is the indexed GDC clinical export/API surface, which is not a supplement download.
  - Pathology PDFs are also missing locally, but they are not the main structured-treatment gap.
- Best next step:
  - First compare the currently parsed clinical treatment biotab coverage against the indexed GDC clinical treatment surface.
  - Then, only if material treatment evidence still looks missing, inspect whether downloaded `Clinical Supplement / BCR XML` or `Clinical Supplement / BCR OMF XML` add treatment fields not already represented in the parsed biotabs or indexed clinical surface.
  - Do not treat the current `278` no-drug-or-radiation patients as a download failure, and do not call them source-level absent yet.

## References

- `R1` GDC Clinical Data: <https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Data/>
- `R2` GDC Clinical Supplement: <https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Supplement/>
- `R3` GDC Biospecimen Data: <https://docs.gdc.cancer.gov/Encyclopedia/pages/Biospecimen_Data/>
- `R4` GDC Repository Users Guide: <https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Repository/>
- `R5` GDC Data Dictionary `clinical_supplement`: <https://api.gdc.cancer.gov/v0/submission/_dictionary/clinical_supplement>
- `R6` GDC Data Dictionary `biospecimen_supplement`: <https://api.gdc.cancer.gov/v0/submission/_dictionary/biospecimen_supplement>
- `R7` GDC Data Dictionary `treatment`: <https://api.gdc.cancer.gov/v0/submission/_dictionary/treatment>
- `R8` GDC Data Dictionary `follow_up`: <https://api.gdc.cancer.gov/v0/submission/_dictionary/follow_up>
- `R9` GDC Files API: <https://api.gdc.cancer.gov/files>
- `R10` GDC Cases API: <https://api.gdc.cancer.gov/cases>
- `R11` GDC API Status: <https://api.gdc.cancer.gov/status>

The exact filtered API request/response snapshots used for this audit are saved under:
- `01-data/audit/tcga-brca/source/treatment_source_audit_runs/20260415T233526Z/`
