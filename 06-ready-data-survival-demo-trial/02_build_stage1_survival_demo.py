from __future__ import annotations

import ast
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter


ROOT = Path(__file__).resolve().parent
RAW_DIR = ROOT / "02-raw-downloads"
PROCESSED_DIR = ROOT / "03-processed"
MODELS_DIR = ROOT / "04-models"
FIGURES_DIR = ROOT / "05-figures"
REPORTS_DIR = ROOT / "06-reports"
DEMO_DIR = ROOT / "07-demo-assets"

TCGA_GDC_CLINICAL = RAW_DIR / "tcga_brca_xena_gdc_clinical.tsv.gz"
TCGA_GDC_SURVIVAL = RAW_DIR / "tcga_brca_xena_gdc_survival.tsv.gz"
TCGA_XENA_AUX = RAW_DIR / "tcga_brca_xena_tcgahub_clinical_matrix.tsv"
METABRIC_ATTR = RAW_DIR / "metabric_cbio_clinical_attributes.json"
METABRIC_SAMPLE = RAW_DIR / "metabric_cbio_clinical_data_sample.json"
METABRIC_PATIENT = RAW_DIR / "metabric_cbio_clinical_data_patient.json"

MIN_GROUP_FOR_PLOT = 20


@dataclass
class ModelArtifacts:
    summary: pd.DataFrame
    metrics: dict[str, Any]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def ensure_dirs() -> None:
    for path in [PROCESSED_DIR, MODELS_DIR, FIGURES_DIR, REPORTS_DIR, DEMO_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def require_inputs() -> None:
    required = [
        TCGA_GDC_CLINICAL,
        TCGA_GDC_SURVIVAL,
        TCGA_XENA_AUX,
        METABRIC_ATTR,
        METABRIC_SAMPLE,
        METABRIC_PATIENT,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing downloaded inputs. Run 01_download_ready_data.py first.\n"
            + "\n".join(missing)
        )


def first_non_null(series: pd.Series) -> Any:
    values = [value for value in series.tolist() if pd.notna(value)]
    return values[0] if values else np.nan


def unique_non_null_count(series: pd.Series) -> int:
    return int(pd.Series(series.dropna().astype(str).unique()).shape[0])


def simplify_stage(raw_value: Any) -> Any:
    if pd.isna(raw_value):
        return np.nan
    text = str(raw_value).strip()
    if not text:
        return np.nan
    text = text.replace("Stage ", "").replace("STAGE ", "").upper()
    if text.startswith("IV"):
        return "IV"
    if text.startswith("III"):
        return "III"
    if text.startswith("II"):
        return "II"
    if text.startswith("I"):
        return "I"
    if text.startswith("X"):
        return "X"
    return text


def parse_python_list_literal(raw_value: Any) -> list[str]:
    if pd.isna(raw_value):
        return []
    text = str(raw_value).strip()
    if not text:
        return []
    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return [text]
    if isinstance(parsed, list):
        return [str(item) for item in parsed]
    return [str(parsed)]


def derive_tcga_treatment_flags(type_value: Any, therapy_value: Any) -> dict[str, Any]:
    treatment_types = parse_python_list_literal(type_value)
    treatment_flags = parse_python_list_literal(therapy_value)
    observed = {
        "tcga_pharmaceutical_therapy_observed": np.nan,
        "tcga_radiation_therapy_observed": np.nan,
    }
    for treatment_type, treatment_flag in zip(treatment_types, treatment_flags):
        normalized_type = treatment_type.lower()
        normalized_flag = treatment_flag.lower()
        if "pharmaceutical" in normalized_type:
            observed["tcga_pharmaceutical_therapy_observed"] = normalized_flag
        if "radiation" in normalized_type:
            observed["tcga_radiation_therapy_observed"] = normalized_flag
    return observed


def parse_status_before_colon(raw_value: Any) -> Any:
    if pd.isna(raw_value):
        return np.nan
    text = str(raw_value).strip()
    if not text:
        return np.nan
    prefix = text.split(":", 1)[0].strip()
    try:
        return int(prefix)
    except ValueError:
        return np.nan


def to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def write_tsv(frame: pd.DataFrame, destination: Path) -> None:
    frame.to_csv(destination, sep="\t", index=False)


def frame_to_text(frame: pd.DataFrame) -> str:
    return frame.to_string(index=False)


def build_tcga_table() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    clinical = pd.read_csv(TCGA_GDC_CLINICAL, sep="\t")
    survival = pd.read_csv(TCGA_GDC_SURVIVAL, sep="\t")
    aux = pd.read_csv(TCGA_XENA_AUX, sep="\t")

    clinical["patient_id"] = clinical["submitter_id"].fillna(
        clinical["sample"].astype(str).str.slice(0, 12)
    )
    survival["patient_id"] = survival["_PATIENT"]
    aux["patient_id"] = aux["_PATIENT"].fillna(aux["bcr_patient_barcode"])

    aux = aux[
        aux.get("sample_type", pd.Series("", index=aux.index)).fillna("").eq("Primary Tumor")
        | aux.get("bcr_sample_barcode", pd.Series("", index=aux.index))
        .fillna("")
        .str.contains("-01", regex=False)
    ].copy()

    clinical_conflict_columns = [
        "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived",
        "ajcc_pathologic_stage.diagnoses",
        "prior_treatment.diagnoses",
        "treatment_type.treatments.diagnoses",
        "treatment_or_therapy.treatments.diagnoses",
    ]
    aux_conflict_columns = [
        "PAM50Call_RNAseq",
        "breast_carcinoma_estrogen_receptor_status",
        "breast_carcinoma_progesterone_receptor_status",
        "lab_proc_her2_neu_immunohistochemistry_receptor_status",
        "history_of_neoadjuvant_treatment",
    ]

    tcga_clinical_conflicts = []
    for column in clinical_conflict_columns:
        conflicts = int(
            clinical.groupby("patient_id")[column].agg(unique_non_null_count).gt(1).sum()
        )
        tcga_clinical_conflicts.append(f"{column}={conflicts}")

    tcga_aux_conflicts = []
    for column in aux_conflict_columns:
        conflicts = int(aux.groupby("patient_id")[column].agg(unique_non_null_count).gt(1).sum())
        tcga_aux_conflicts.append(f"{column}={conflicts}")

    clinical_patient = clinical.groupby("patient_id", as_index=False).agg(
        {
            "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived": first_non_null,
            "ajcc_pathologic_stage.diagnoses": first_non_null,
            "prior_treatment.diagnoses": first_non_null,
            "treatment_type.treatments.diagnoses": first_non_null,
            "treatment_or_therapy.treatments.diagnoses": first_non_null,
        }
    )
    survival_patient = survival.groupby("patient_id", as_index=False).agg(
        {"OS.time": first_non_null, "OS": first_non_null}
    )
    aux_patient = aux.groupby("patient_id", as_index=False).agg(
        {
            "PAM50Call_RNAseq": first_non_null,
            "breast_carcinoma_estrogen_receptor_status": first_non_null,
            "breast_carcinoma_progesterone_receptor_status": first_non_null,
            "lab_proc_her2_neu_immunohistochemistry_receptor_status": first_non_null,
            "history_of_neoadjuvant_treatment": first_non_null,
        }
    )

    tcga = (
        survival_patient.merge(clinical_patient, on="patient_id", how="left")
        .merge(aux_patient, on="patient_id", how="left")
        .rename(
            columns={
                "OS.time": "os_time_days",
                "OS": "os_event",
                "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived": "age_at_diagnosis_years",
                "ajcc_pathologic_stage.diagnoses": "ajcc_pathologic_stage_raw",
                "PAM50Call_RNAseq": "pam50_subtype_raw",
                "breast_carcinoma_estrogen_receptor_status": "er_status_raw",
                "breast_carcinoma_progesterone_receptor_status": "pr_status_raw",
                "lab_proc_her2_neu_immunohistochemistry_receptor_status": "her2_ihc_status_raw",
                "history_of_neoadjuvant_treatment": "history_of_neoadjuvant_treatment_raw",
                "prior_treatment.diagnoses": "prior_treatment_raw",
                "treatment_type.treatments.diagnoses": "gdc_treatment_type_raw",
                "treatment_or_therapy.treatments.diagnoses": "gdc_treatment_or_therapy_raw",
            }
        )
    )

    tcga["cohort"] = "TCGA-BRCA"
    tcga["os_time_days"] = to_numeric(tcga["os_time_days"])
    tcga["os_event"] = to_numeric(tcga["os_event"]).fillna(0).astype(int)
    tcga["age_at_diagnosis_years"] = to_numeric(tcga["age_at_diagnosis_years"])
    tcga["ajcc_stage_simple"] = tcga["ajcc_pathologic_stage_raw"].map(simplify_stage)
    tcga["subtype_label"] = tcga["pam50_subtype_raw"]
    tcga["subtype_field_used"] = np.where(tcga["pam50_subtype_raw"].notna(), "PAM50Call_RNAseq", None)
    tcga["subtype_source_dataset"] = np.where(
        tcga["pam50_subtype_raw"].notna(),
        "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
        None,
    )

    treatment_flags = tcga.apply(
        lambda row: derive_tcga_treatment_flags(
            row["gdc_treatment_type_raw"], row["gdc_treatment_or_therapy_raw"]
        ),
        axis=1,
        result_type="expand",
    )
    tcga = pd.concat([tcga, treatment_flags], axis=1)
    tcga["analysis_ready_os"] = tcga["os_time_days"].notna() & tcga["os_event"].isin([0, 1])

    inventory_rows = [
        {
            "cohort": "TCGA-BRCA",
            "source_dataset": "TCGA-BRCA.survival.tsv",
            "standardized_domain": "survival",
            "output_column": "os_time_days",
            "source_column": "OS.time",
            "field_role": "chosen_primary",
            "cohort_specific": True,
            "unit_or_encoding": "days",
            "non_missing_n": int(tcga["os_time_days"].notna().sum()),
            "notes": "Primary overall survival duration from UCSC Xena GDC hub.",
        },
        {
            "cohort": "TCGA-BRCA",
            "source_dataset": "TCGA-BRCA.survival.tsv",
            "standardized_domain": "survival",
            "output_column": "os_event",
            "source_column": "OS",
            "field_role": "chosen_primary",
            "cohort_specific": True,
            "unit_or_encoding": "0=censored, 1=event",
            "non_missing_n": int(tcga["os_event"].notna().sum()),
            "notes": "Primary overall survival event indicator from UCSC Xena GDC hub.",
        },
        {
            "cohort": "TCGA-BRCA",
            "source_dataset": "TCGA-BRCA.clinical.tsv",
            "standardized_domain": "age",
            "output_column": "age_at_diagnosis_years",
            "source_column": "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived",
            "field_role": "chosen_primary",
            "cohort_specific": False,
            "unit_or_encoding": "years",
            "non_missing_n": int(tcga["age_at_diagnosis_years"].notna().sum()),
            "notes": "Xena-derived age in years from GDC-hub clinical table.",
        },
        {
            "cohort": "TCGA-BRCA",
            "source_dataset": "TCGA-BRCA.clinical.tsv",
            "standardized_domain": "stage",
            "output_column": "ajcc_stage_simple",
            "source_column": "ajcc_pathologic_stage.diagnoses",
            "field_role": "chosen_primary",
            "cohort_specific": True,
            "unit_or_encoding": "simplified AJCC stage group",
            "non_missing_n": int(tcga["ajcc_stage_simple"].notna().sum()),
            "notes": "Derived by collapsing raw pathologic stage labels to I/II/III/IV/X.",
        },
        {
            "cohort": "TCGA-BRCA",
            "source_dataset": "TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
            "standardized_domain": "subtype",
            "output_column": "subtype_label",
            "source_column": "PAM50Call_RNAseq",
            "field_role": "chosen_auxiliary",
            "cohort_specific": True,
            "unit_or_encoding": "PAM50 label",
            "non_missing_n": int(tcga["subtype_label"].notna().sum()),
            "notes": "Auxiliary Xena TCGA-hub subtype field because GDC-hub clinical table does not expose subtype labels.",
        },
        {
            "cohort": "TCGA-BRCA",
            "source_dataset": "TCGA-BRCA.clinical.tsv",
            "standardized_domain": "treatment",
            "output_column": "tcga_pharmaceutical_therapy_observed",
            "source_column": "treatment_or_therapy.treatments.diagnoses + treatment_type.treatments.diagnoses",
            "field_role": "descriptive_only",
            "cohort_specific": True,
            "unit_or_encoding": "yes/no/not reported",
            "non_missing_n": int(tcga["tcga_pharmaceutical_therapy_observed"].notna().sum()),
            "notes": "Derived observational descriptor only. Not used for causal or treatment-benefit claims.",
        },
        {
            "cohort": "TCGA-BRCA",
            "source_dataset": "TCGA-BRCA.clinical.tsv",
            "standardized_domain": "treatment",
            "output_column": "tcga_radiation_therapy_observed",
            "source_column": "treatment_or_therapy.treatments.diagnoses + treatment_type.treatments.diagnoses",
            "field_role": "descriptive_only",
            "cohort_specific": True,
            "unit_or_encoding": "yes/no/not reported",
            "non_missing_n": int(tcga["tcga_radiation_therapy_observed"].notna().sum()),
            "notes": "Derived observational descriptor only. Not used for causal or treatment-benefit claims.",
        },
    ]

    audit_notes = [
        "TCGA patient table is collapsed to one row per patient using first non-null values.",
        f"TCGA clinical duplicate-conflict counts after patient grouping check: {', '.join(tcga_clinical_conflicts)}.",
        f"TCGA auxiliary duplicate-conflict counts after tumor-row restriction: {', '.join(tcga_aux_conflicts)}.",
        "TCGA auxiliary subtype table is restricted to primary-tumor rows before patient collapse to avoid tumor-vs-normal PAM50 conflicts.",
    ]
    return tcga, pd.DataFrame(inventory_rows), audit_notes


def load_json(path: Path) -> list[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_metabric_table() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    attrs = pd.DataFrame(load_json(METABRIC_ATTR))
    sample_long = pd.DataFrame(load_json(METABRIC_SAMPLE))
    patient_long = pd.DataFrame(load_json(METABRIC_PATIENT))

    sample_wide = sample_long.pivot_table(
        index="patientId",
        columns="clinicalAttributeId",
        values="value",
        aggfunc="first",
    ).reset_index()
    sample_wide.columns.name = None
    patient_wide = patient_long.pivot_table(
        index="patientId",
        columns="clinicalAttributeId",
        values="value",
        aggfunc="first",
    ).reset_index()
    patient_wide.columns.name = None

    metabric = patient_wide.merge(sample_wide, on="patientId", how="left", suffixes=("", "_sample"))
    metabric = metabric.rename(
        columns={
            "patientId": "patient_id",
            "AGE_AT_DIAGNOSIS": "age_at_diagnosis_years",
            "OS_MONTHS": "os_time_months",
            "OS_STATUS": "os_status_raw",
            "RFS_MONTHS": "rfs_time_months",
            "RFS_STATUS": "rfs_status_raw",
            "CLAUDIN_SUBTYPE": "claudin_subtype_raw",
            "THREEGENE": "threegene_subtype_raw",
            "INTCLUST": "intclust_raw",
            "ER_IHC": "er_status_ihc_raw",
            "HER2_SNP6": "her2_snp6_raw",
            "CHEMOTHERAPY": "chemotherapy_observed_raw",
            "HORMONE_THERAPY": "hormone_therapy_observed_raw",
            "RADIO_THERAPY": "radio_therapy_observed_raw",
            "BREAST_SURGERY": "breast_surgery_raw",
            "INFERRED_MENOPAUSAL_STATE": "menopausal_state_raw",
            "NPI": "npi_raw",
            "TUMOR_STAGE": "tumor_stage_raw",
            "HISTOLOGICAL_SUBTYPE": "histological_subtype_raw",
        }
    )

    metabric["cohort"] = "METABRIC"
    metabric["age_at_diagnosis_years"] = to_numeric(metabric["age_at_diagnosis_years"])
    metabric["os_time_months"] = to_numeric(metabric["os_time_months"])
    metabric["os_event"] = metabric["os_status_raw"].map(parse_status_before_colon)
    metabric["rfs_time_months"] = to_numeric(metabric["rfs_time_months"])
    metabric["rfs_event"] = metabric["rfs_status_raw"].map(parse_status_before_colon)
    metabric["npi"] = to_numeric(metabric["npi_raw"])
    metabric["subtype_label"] = metabric["claudin_subtype_raw"]
    metabric["subtype_field_used"] = np.where(
        metabric["claudin_subtype_raw"].notna(),
        "CLAUDIN_SUBTYPE",
        None,
    )
    metabric["subtype_source_dataset"] = np.where(
        metabric["claudin_subtype_raw"].notna(),
        "cBioPortal patient-level clinical data",
        None,
    )
    metabric["analysis_ready_os"] = metabric["os_time_months"].notna() & metabric["os_event"].isin([0, 1])
    metabric["analysis_ready_rfs"] = metabric["rfs_time_months"].notna() & metabric["rfs_event"].isin([0, 1])

    inventory_rows = [
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "survival",
            "output_column": "os_time_months",
            "source_column": "OS_MONTHS",
            "field_role": "chosen_primary",
            "cohort_specific": False,
            "unit_or_encoding": "months",
            "non_missing_n": int(metabric["os_time_months"].notna().sum()),
            "notes": "Primary METABRIC overall survival duration.",
        },
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "survival",
            "output_column": "os_event",
            "source_column": "OS_STATUS",
            "field_role": "chosen_primary",
            "cohort_specific": False,
            "unit_or_encoding": "0=living, 1=deceased",
            "non_missing_n": int(metabric["os_event"].notna().sum()),
            "notes": "Primary METABRIC overall survival event indicator.",
        },
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "survival",
            "output_column": "rfs_time_months",
            "source_column": "RFS_MONTHS",
            "field_role": "available_secondary",
            "cohort_specific": True,
            "unit_or_encoding": "months",
            "non_missing_n": int(metabric["rfs_time_months"].notna().sum()),
            "notes": "Retained for endpoint inventory only. The baseline model uses OS.",
        },
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "subtype",
            "output_column": "subtype_label",
            "source_column": "CLAUDIN_SUBTYPE",
            "field_role": "chosen_primary",
            "cohort_specific": True,
            "unit_or_encoding": "subtype label",
            "non_missing_n": int(metabric["subtype_label"].notna().sum()),
            "notes": "Chosen METABRIC subtype field for distribution summary and KM plots.",
        },
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "subtype",
            "output_column": "threegene_subtype_raw",
            "source_column": "THREEGENE",
            "field_role": "available_secondary",
            "cohort_specific": True,
            "unit_or_encoding": "subtype label",
            "non_missing_n": int(metabric["threegene_subtype_raw"].notna().sum()),
            "notes": "Retained as a cohort-specific secondary subtype field.",
        },
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "treatment",
            "output_column": "chemotherapy_observed_raw",
            "source_column": "CHEMOTHERAPY",
            "field_role": "descriptive_only",
            "cohort_specific": True,
            "unit_or_encoding": "YES/NO",
            "non_missing_n": int(metabric["chemotherapy_observed_raw"].notna().sum()),
            "notes": "Observational treatment descriptor only. Not used for treatment-effect claims.",
        },
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "treatment",
            "output_column": "hormone_therapy_observed_raw",
            "source_column": "HORMONE_THERAPY",
            "field_role": "descriptive_only",
            "cohort_specific": True,
            "unit_or_encoding": "YES/NO",
            "non_missing_n": int(metabric["hormone_therapy_observed_raw"].notna().sum()),
            "notes": "Observational treatment descriptor only. Not used for treatment-effect claims.",
        },
        {
            "cohort": "METABRIC",
            "source_dataset": "cBioPortal patient clinical-data",
            "standardized_domain": "treatment",
            "output_column": "radio_therapy_observed_raw",
            "source_column": "RADIO_THERAPY",
            "field_role": "descriptive_only",
            "cohort_specific": True,
            "unit_or_encoding": "YES/NO",
            "non_missing_n": int(metabric["radio_therapy_observed_raw"].notna().sum()),
            "notes": "Observational treatment descriptor only. Not used for treatment-effect claims.",
        },
    ]

    audit_notes = [
        f"METABRIC attribute count downloaded from cBioPortal: {attrs.shape[0]}.",
        f"METABRIC patient-level clinical long rows: {patient_long.shape[0]}.",
        f"METABRIC sample-level clinical long rows: {sample_long.shape[0]}.",
        f"METABRIC merged wide table rows: {metabric.shape[0]}.",
    ]
    return metabric, pd.DataFrame(inventory_rows), audit_notes


def fit_cox_model(
    cohort_name: str,
    analysis_df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    covariate_columns: list[str],
    formula: str,
) -> tuple[pd.DataFrame, ModelArtifacts]:
    model_frame = analysis_df[[duration_col, event_col, *covariate_columns]].copy()
    model_frame = model_frame.dropna(axis=0).copy()
    if model_frame.empty:
        empty_summary = pd.DataFrame(
            [
                {
                    "cohort": cohort_name,
                    "term": "NO_MODEL",
                    "exp(coef)": np.nan,
                    "coef": np.nan,
                    "p": np.nan,
                    "coef lower 95%": np.nan,
                    "coef upper 95%": np.nan,
                    "notes": "No complete cases available for the requested baseline Cox model.",
                }
            ]
        )
        empty_metrics = {
            "cohort": cohort_name,
            "duration_column": duration_col,
            "event_column": event_col,
            "complete_cases": 0,
            "concordance_index": np.nan,
            "covariates": covariate_columns,
            "formula": formula,
        }
        return model_frame, ModelArtifacts(summary=empty_summary, metrics=empty_metrics)

    cph = CoxPHFitter()
    cph.fit(model_frame, duration_col=duration_col, event_col=event_col, formula=formula)
    summary = cph.summary.reset_index().rename(columns={"covariate": "term"})
    summary.insert(0, "cohort", cohort_name)
    summary["notes"] = "Baseline prognostic Cox model. Treatment fields excluded from the model."
    metrics = {
        "cohort": cohort_name,
        "duration_column": duration_col,
        "event_column": event_col,
        "complete_cases": int(model_frame.shape[0]),
        "concordance_index": float(cph.concordance_index_),
        "covariates": covariate_columns,
        "formula": formula,
    }
    risk_score = cph.predict_log_partial_hazard(model_frame).rename("baseline_risk_score")
    model_frame = model_frame.join(risk_score)
    return model_frame, ModelArtifacts(summary=summary, metrics=metrics)


def attach_risk_scores(
    analysis_df: pd.DataFrame,
    patient_id_col: str,
    model_frame: pd.DataFrame,
) -> pd.DataFrame:
    enriched = analysis_df.copy()
    enriched["baseline_risk_score"] = np.nan
    enriched["cox_complete_case"] = False
    if "baseline_risk_score" in model_frame.columns:
        joinable = model_frame[["baseline_risk_score"]].copy()
        joinable[patient_id_col] = enriched.loc[joinable.index, patient_id_col].values
        joinable = joinable.drop_duplicates(subset=[patient_id_col])
        enriched = enriched.merge(joinable, on=patient_id_col, how="left", suffixes=("", "_model"))
        if "baseline_risk_score_model" in enriched.columns:
            enriched["baseline_risk_score"] = enriched["baseline_risk_score_model"]
            enriched = enriched.drop(columns=["baseline_risk_score_model"])
        enriched["cox_complete_case"] = enriched["baseline_risk_score"].notna()
    return enriched


def build_cohort_summary(tcga: pd.DataFrame, metabric: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "cohort": "TCGA-BRCA",
            "patients_total": int(tcga.shape[0]),
            "os_complete_n": int(tcga["analysis_ready_os"].sum()),
            "os_event_n": int(tcga.loc[tcga["analysis_ready_os"], "os_event"].sum()),
            "primary_subtype_field": "PAM50Call_RNAseq",
            "primary_subtype_non_missing_n": int(tcga["subtype_label"].notna().sum()),
            "os_time_unit": "days",
            "notes": "Subtype is augmented from the UCSC Xena TCGA hub clinical matrix.",
        },
        {
            "cohort": "METABRIC",
            "patients_total": int(metabric.shape[0]),
            "os_complete_n": int(metabric["analysis_ready_os"].sum()),
            "os_event_n": int(metabric.loc[metabric["analysis_ready_os"], "os_event"].sum()),
            "primary_subtype_field": "CLAUDIN_SUBTYPE",
            "primary_subtype_non_missing_n": int(metabric["subtype_label"].notna().sum()),
            "os_time_unit": "months",
            "notes": "Subtype comes directly from cBioPortal patient clinical data.",
        },
    ]
    return pd.DataFrame(rows)


def build_survival_summary(tcga: pd.DataFrame, metabric: pd.DataFrame) -> pd.DataFrame:
    rows = [
        {
            "cohort": "TCGA-BRCA",
            "endpoint": "OS",
            "time_column": "os_time_days",
            "time_unit": "days",
            "complete_n": int(tcga["analysis_ready_os"].sum()),
            "event_n": int(tcga.loc[tcga["analysis_ready_os"], "os_event"].sum()),
            "median_time_observed": float(tcga.loc[tcga["analysis_ready_os"], "os_time_days"].median()),
        },
        {
            "cohort": "METABRIC",
            "endpoint": "OS",
            "time_column": "os_time_months",
            "time_unit": "months",
            "complete_n": int(metabric["analysis_ready_os"].sum()),
            "event_n": int(metabric.loc[metabric["analysis_ready_os"], "os_event"].sum()),
            "median_time_observed": float(
                metabric.loc[metabric["analysis_ready_os"], "os_time_months"].median()
            ),
        },
        {
            "cohort": "METABRIC",
            "endpoint": "RFS",
            "time_column": "rfs_time_months",
            "time_unit": "months",
            "complete_n": int(metabric["analysis_ready_rfs"].sum()),
            "event_n": int(metabric.loc[metabric["analysis_ready_rfs"], "rfs_event"].sum()),
            "median_time_observed": float(
                metabric.loc[metabric["analysis_ready_rfs"], "rfs_time_months"].median()
            ),
        },
    ]
    return pd.DataFrame(rows)


def build_subtype_summary(cohort_df: pd.DataFrame, cohort_name: str) -> pd.DataFrame:
    grouped = (
        cohort_df["subtype_label"]
        .fillna("Missing")
        .value_counts(dropna=False)
        .rename_axis("subtype_label")
        .reset_index(name="patients_n")
    )
    grouped.insert(0, "cohort", cohort_name)
    grouped["pct_of_cohort"] = (grouped["patients_n"] / cohort_df.shape[0] * 100).round(2)
    grouped["subtype_field_used"] = cohort_df["subtype_field_used"].dropna().iloc[0]
    return grouped


def build_treatment_summary(tcga: pd.DataFrame, metabric: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort_name, frame, columns in [
        (
            "TCGA-BRCA",
            tcga,
            [
                "prior_treatment_raw",
                "history_of_neoadjuvant_treatment_raw",
                "tcga_pharmaceutical_therapy_observed",
                "tcga_radiation_therapy_observed",
            ],
        ),
        (
            "METABRIC",
            metabric,
            [
                "chemotherapy_observed_raw",
                "hormone_therapy_observed_raw",
                "radio_therapy_observed_raw",
                "breast_surgery_raw",
            ],
        ),
    ]:
        for column in columns:
            value_counts = frame[column].fillna("Missing").astype(str).value_counts()
            for value, count in value_counts.items():
                rows.append(
                    {
                        "cohort": cohort_name,
                        "field": column,
                        "value": value,
                        "patients_n": int(count),
                        "pct_of_cohort": round(count / frame.shape[0] * 100, 2),
                    }
                )
    return pd.DataFrame(rows)


def plot_km_by_subtype(
    cohort_df: pd.DataFrame,
    cohort_name: str,
    duration_col: str,
    event_col: str,
    time_unit: str,
    destination: Path,
) -> list[str]:
    plot_df = cohort_df[[duration_col, event_col, "subtype_label"]].dropna().copy()
    group_sizes = plot_df["subtype_label"].value_counts()
    keep_groups = group_sizes[group_sizes >= MIN_GROUP_FOR_PLOT].index.tolist()
    plot_df = plot_df[plot_df["subtype_label"].isin(keep_groups)]
    if plot_df.empty:
        return []

    plt.figure(figsize=(10, 7))
    kmf = KaplanMeierFitter()
    for subtype in sorted(keep_groups):
        subtype_df = plot_df[plot_df["subtype_label"] == subtype]
        kmf.fit(subtype_df[duration_col], subtype_df[event_col], label=f"{subtype} (n={len(subtype_df)})")
        kmf.plot(ci_show=False)

    plt.title(f"{cohort_name} overall survival by subtype")
    plt.xlabel(f"Time ({time_unit})")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.savefig(destination, dpi=200)
    plt.close()
    return keep_groups


def plot_subtype_distribution(cohort_df: pd.DataFrame, cohort_name: str, destination: Path) -> None:
    plot_df = cohort_df["subtype_label"].fillna("Missing").value_counts().reset_index()
    plot_df.columns = ["subtype_label", "patients_n"]

    plt.figure(figsize=(10, 6))
    sns.barplot(data=plot_df, x="patients_n", y="subtype_label", color="#4C78A8")
    plt.title(f"{cohort_name} subtype distribution")
    plt.xlabel("Patients")
    plt.ylabel("Subtype")
    plt.tight_layout()
    plt.savefig(destination, dpi=200)
    plt.close()


def plot_cox_coefficients(summary_df: pd.DataFrame, cohort_name: str, destination: Path) -> None:
    plot_df = summary_df.copy()
    plot_df = plot_df[plot_df["term"] != "NO_MODEL"]
    if plot_df.empty:
        return

    plot_df = plot_df.sort_values("exp(coef)")
    y_positions = np.arange(plot_df.shape[0])
    lower = plot_df["exp(coef)"] - np.exp(plot_df["coef lower 95%"])
    upper = np.exp(plot_df["coef upper 95%"]) - plot_df["exp(coef)"]

    plt.figure(figsize=(10, max(4, 0.5 * plot_df.shape[0] + 2)))
    plt.errorbar(
        plot_df["exp(coef)"],
        y_positions,
        xerr=[lower, upper],
        fmt="o",
        color="#E45756",
        ecolor="#999999",
        capsize=3,
    )
    plt.axvline(1.0, color="black", linestyle="--", linewidth=1)
    plt.yticks(y_positions, plot_df["term"])
    plt.xlabel("Hazard ratio")
    plt.title(f"{cohort_name} baseline Cox model coefficients")
    plt.tight_layout()
    plt.savefig(destination, dpi=200)
    plt.close()


def build_demo_summary(
    cohort_summary: pd.DataFrame,
    model_metrics: list[dict[str, Any]],
    km_groups: dict[str, list[str]],
) -> None:
    tcga_row = cohort_summary.loc[cohort_summary["cohort"] == "TCGA-BRCA"].iloc[0]
    meta_row = cohort_summary.loc[cohort_summary["cohort"] == "METABRIC"].iloc[0]

    lines = [
        "# Demo Summary",
        "",
        f"Generated: {utc_now_iso()}",
        "",
        "## Scope",
        "",
        "- This package is a Stage 1 prognostic survival demo.",
        "- It uses public ready-made non-image data only.",
        "- It does not estimate treatment effects or recommend treatments.",
        "",
        "## Cohorts",
        "",
        f"- TCGA-BRCA: {int(tcga_row['patients_total'])} patients with {int(tcga_row['os_complete_n'])} analyzable OS records.",
        f"- METABRIC: {int(meta_row['patients_total'])} patients with {int(meta_row['os_complete_n'])} analyzable OS records.",
        "",
        "## Subtype Choices",
        "",
        "- TCGA-BRCA subtype field used for plots: `PAM50Call_RNAseq` from UCSC Xena TCGA hub clinical matrix.",
        "- METABRIC subtype field used for plots: `CLAUDIN_SUBTYPE` from cBioPortal patient clinical data.",
        "",
        "## Kaplan-Meier Plots",
        "",
        f"- TCGA-BRCA plotted subtype groups: {', '.join(km_groups.get('TCGA-BRCA', [])) or 'none'}",
        f"- METABRIC plotted subtype groups: {', '.join(km_groups.get('METABRIC', [])) or 'none'}",
        "",
        "## Baseline Cox Models",
        "",
    ]
    for metric in model_metrics:
        c_index = metric["concordance_index"]
        c_index_text = "NA" if pd.isna(c_index) else f"{c_index:.3f}"
        lines.append(
            f"- {metric['cohort']}: complete cases={metric['complete_cases']}, concordance_index={c_index_text}, formula=`{metric['formula']}`"
        )

    lines.extend(
        [
            "",
            "## Meeting-Safe Interpretation",
            "",
            "- These models are prognostic association models within each cohort.",
            "- Treatment variables are retained only as observational context and availability summaries.",
            "- No causal treatment recommendation claim is made anywhere in this package.",
            "",
        ]
    )
    (REPORTS_DIR / "demo_summary.md").write_text("\n".join(lines), encoding="utf-8")


def build_data_audit(
    cohort_summary: pd.DataFrame,
    survival_summary: pd.DataFrame,
    subtype_summary: pd.DataFrame,
    treatment_summary: pd.DataFrame,
    audit_notes: list[str],
) -> None:
    lines = [
        "# Data Audit",
        "",
        f"Generated: {utc_now_iso()}",
        "",
        "## Cohort Summary",
        "",
        frame_to_text(cohort_summary),
        "",
        "## Survival Endpoint Summary",
        "",
        frame_to_text(survival_summary),
        "",
        "## Subtype Distribution Summary",
        "",
        frame_to_text(subtype_summary),
        "",
        "## Treatment Availability Summary",
        "",
        frame_to_text(treatment_summary),
        "",
        "## Audit Notes",
        "",
    ]
    lines.extend(f"- {note}" for note in audit_notes)
    lines.extend(
        [
            "",
            "## Explicit Caveats",
            "",
            "- TCGA subtype information required a documented auxiliary UCSC Xena TCGA-hub table because the GDC-hub clinical table did not carry subtype labels.",
            "- TCGA auxiliary subtype rows were restricted to primary-tumor samples before collapse to avoid tumor-vs-normal subtype conflicts.",
            "- METABRIC treatment variables are observational binary descriptors and not randomized exposures.",
            "",
        ]
    )
    (REPORTS_DIR / "data_audit.md").write_text("\n".join(lines), encoding="utf-8")


def build_safe_claims_file() -> None:
    lines = [
        "# Safe Claims For Meeting",
        "",
        "## What You Can Say",
        "",
        "- This is a Stage 1 prognostic survival demo built from public ready-made BRCA cohort tables.",
        "- The package shows cohort-specific survival endpoint preparation, subtype summaries, Kaplan-Meier curves, and baseline interpretable Cox models.",
        "- Treatment fields are present only as observational and incomplete clinical descriptors.",
        "- Treatment variables are summarized for data availability and context only.",
        "- No treatment-effect estimation is being claimed.",
        "- No treatment recommendation is being made.",
        "",
        "## What You Should Not Say",
        "",
        "- Do not say the model improves survival.",
        "- Do not say the model predicts treatment benefit.",
        "- Do not say a patient should receive a specific treatment based on these outputs.",
        "- Do not describe this as causal inference or treatment decision support at this stage.",
        "",
        "## Safe One-Paragraph Version",
        "",
        "This is a prognostic Stage 1 demo using public BRCA cohort tables from UCSC Xena and cBioPortal. "
        "It demonstrates a clean workflow for downloading, auditing, harmonizing, and modeling survival endpoints with subtype-aware summaries. "
        "Treatment-related fields are observational and incomplete, so they are included only for descriptive context. "
        "No causal treatment recommendation claim is being made.",
        "",
    ]
    (REPORTS_DIR / "safe_claims_for_meeting.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    require_inputs()
    sns.set_theme(style="whitegrid")

    tcga, tcga_inventory, tcga_audit_notes = build_tcga_table()
    metabric, metabric_inventory, metabric_audit_notes = build_metabric_table()

    tcga_model_frame, tcga_model = fit_cox_model(
        cohort_name="TCGA-BRCA",
        analysis_df=tcga,
        duration_col="os_time_days",
        event_col="os_event",
        covariate_columns=["age_at_diagnosis_years", "ajcc_stage_simple", "subtype_label"],
        formula="age_at_diagnosis_years + C(ajcc_stage_simple) + C(subtype_label)",
    )
    metabric_model_frame, metabric_model = fit_cox_model(
        cohort_name="METABRIC",
        analysis_df=metabric,
        duration_col="os_time_months",
        event_col="os_event",
        covariate_columns=["age_at_diagnosis_years", "npi", "subtype_label"],
        formula="age_at_diagnosis_years + npi + C(subtype_label)",
    )

    tcga = attach_risk_scores(tcga.reset_index(drop=True), "patient_id", tcga_model_frame)
    metabric = attach_risk_scores(metabric.reset_index(drop=True), "patient_id", metabric_model_frame)

    cohort_summary = build_cohort_summary(tcga, metabric)
    survival_summary = build_survival_summary(tcga, metabric)
    subtype_summary = pd.concat(
        [build_subtype_summary(tcga, "TCGA-BRCA"), build_subtype_summary(metabric, "METABRIC")],
        ignore_index=True,
    )
    treatment_summary = build_treatment_summary(tcga, metabric)
    inventory = pd.concat([tcga_inventory, metabric_inventory], ignore_index=True)

    write_tsv(cohort_summary, PROCESSED_DIR / "cohort_summary.tsv")
    write_tsv(survival_summary, PROCESSED_DIR / "survival_endpoint_summary.tsv")
    write_tsv(subtype_summary, PROCESSED_DIR / "subtype_distribution_summary.tsv")
    write_tsv(treatment_summary, PROCESSED_DIR / "treatment_availability_summary.tsv")
    write_tsv(inventory, PROCESSED_DIR / "harmonized_cohort_inventory.tsv")
    write_tsv(tcga, PROCESSED_DIR / "tcga_brca_analysis_table.tsv")
    write_tsv(metabric, PROCESSED_DIR / "metabric_analysis_table.tsv")
    write_tsv(tcga_model.summary, MODELS_DIR / "tcga_brca_baseline_cox_summary.tsv")
    write_tsv(metabric_model.summary, MODELS_DIR / "metabric_baseline_cox_summary.tsv")
    write_tsv(
        pd.DataFrame([tcga_model.metrics, metabric_model.metrics]),
        MODELS_DIR / "baseline_model_metrics.tsv",
    )

    tcga.to_csv(DEMO_DIR / "demo_table_tcga_brca.csv", index=False)
    metabric.to_csv(DEMO_DIR / "demo_table_metabric.csv", index=False)

    km_groups = {
        "TCGA-BRCA": plot_km_by_subtype(
            tcga,
            "TCGA-BRCA",
            "os_time_days",
            "os_event",
            "days",
            FIGURES_DIR / "tcga_brca_km_by_subtype.png",
        ),
        "METABRIC": plot_km_by_subtype(
            metabric,
            "METABRIC",
            "os_time_months",
            "os_event",
            "months",
            FIGURES_DIR / "metabric_km_by_subtype.png",
        ),
    }
    plot_subtype_distribution(tcga, "TCGA-BRCA", FIGURES_DIR / "tcga_brca_subtype_distribution.png")
    plot_subtype_distribution(metabric, "METABRIC", FIGURES_DIR / "metabric_subtype_distribution.png")
    plot_cox_coefficients(
        tcga_model.summary, "TCGA-BRCA", FIGURES_DIR / "tcga_brca_cox_coefficients.png"
    )
    plot_cox_coefficients(
        metabric_model.summary, "METABRIC", FIGURES_DIR / "metabric_cox_coefficients.png"
    )

    build_demo_summary(
        cohort_summary=cohort_summary,
        model_metrics=[tcga_model.metrics, metabric_model.metrics],
        km_groups=km_groups,
    )
    build_data_audit(
        cohort_summary=cohort_summary,
        survival_summary=survival_summary,
        subtype_summary=subtype_summary,
        treatment_summary=treatment_summary,
        audit_notes=tcga_audit_notes + metabric_audit_notes,
    )
    build_safe_claims_file()

    metrics = {
        "generated_at_utc": utc_now_iso(),
        "scope": "Stage 1 prognostic survival demo",
        "tcga_brca": {
            "patients_total": int(tcga.shape[0]),
            "os_complete_n": int(tcga["analysis_ready_os"].sum()),
            "os_event_n": int(tcga.loc[tcga["analysis_ready_os"], "os_event"].sum()),
            "os_time_unit": "days",
            "subtype_field_used": "PAM50Call_RNAseq",
            "subtype_non_missing_n": int(tcga["subtype_label"].notna().sum()),
            "km_groups_plotted": km_groups["TCGA-BRCA"],
            "cox_complete_cases": tcga_model.metrics["complete_cases"],
            "cox_concordance_index": tcga_model.metrics["concordance_index"],
        },
        "metabric": {
            "patients_total": int(metabric.shape[0]),
            "os_complete_n": int(metabric["analysis_ready_os"].sum()),
            "os_event_n": int(metabric.loc[metabric["analysis_ready_os"], "os_event"].sum()),
            "rfs_complete_n": int(metabric["analysis_ready_rfs"].sum()),
            "rfs_event_n": int(metabric.loc[metabric["analysis_ready_rfs"], "rfs_event"].sum()),
            "os_time_unit": "months",
            "subtype_field_used": "CLAUDIN_SUBTYPE",
            "subtype_non_missing_n": int(metabric["subtype_label"].notna().sum()),
            "km_groups_plotted": km_groups["METABRIC"],
            "cox_complete_cases": metabric_model.metrics["complete_cases"],
            "cox_concordance_index": metabric_model.metrics["concordance_index"],
        },
        "meeting_safe_language": {
            "prognostic_demo_only": True,
            "treatment_effect_estimation_claim": False,
            "treatment_recommendation_claim": False,
        },
    }
    with (DEMO_DIR / "demo_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, indent=2)

    print("Stage 1 demo package built successfully.")
    print(f"  - {PROCESSED_DIR}")
    print(f"  - {MODELS_DIR}")
    print(f"  - {FIGURES_DIR}")
    print(f"  - {REPORTS_DIR}")
    print(f"  - {DEMO_DIR}")


if __name__ == "__main__":
    main()
