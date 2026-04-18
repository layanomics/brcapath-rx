#!/usr/bin/env python
"""Compare indexed GDC clinical treatment coverage against local TCGA-BRCA treatment layers."""

from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

PROJECT_ID = "TCGA-BRCA"
GDC_API_BASE_URL = "https://api.gdc.cancer.gov"
REQUEST_TIMEOUT_SECONDS = 120
USER_AGENT = "brcapath-rx-indexed-clinical-vs-local-treatment/1.0"
INDEXED_CASE_PAGE_SIZE = 200

LOCAL_MISSING_LIKE_TOKENS = {
    "",
    "[discrepancy]",
    "[not applicable]",
    "[not available]",
    "[not evaluated]",
    "[unknown]",
    "n/a",
    "na",
    "nan",
    "none",
    "null",
}
INDEXED_MISSING_LIKE_TOKENS = LOCAL_MISSING_LIKE_TOKENS.union(
    {
        "not allowed to collect",
        "not applicable",
        "not reported",
        "unknown",
        "unspecified",
    }
)
INDEXED_DRUGLIKE_TREATMENT_TYPES = {
    "ancillary treatment",
    "bisphosphonate therapy",
    "chemotherapy",
    "hormone therapy",
    "immunotherapy (including vaccines)",
    "pharmaceutical therapy, nos",
    "targeted molecular therapy",
}
INDEXED_RADIATION_TREATMENT_TYPE_TOKENS = (
    "brachytherapy",
    "gamma knife",
    "radioisotope",
    "radiation",
    "stereotactic",
)

PATIENT_LEVEL_FIELDNAMES = [
    "indexed_clinical_vs_local_treatment_v1_run_id",
    "treatment_source_audit_run_id",
    "patient_treatment_profile_v1_run_id",
    "patient_treatment_grouping_v1_run_id",
    "clinical_biotab_parse_run_id",
    "cohort_join_status",
    "cohort_scope",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "indexed_case_id",
    "local_clinical_any_drug_row",
    "local_clinical_drug_row_count",
    "local_clinical_any_radiation_row",
    "local_clinical_radiation_row_count",
    "local_clinical_drug_therapy_type_values_json",
    "local_profile_has_any_drug_row",
    "local_profile_drug_row_count",
    "local_profile_has_any_radiation_row",
    "local_profile_radiation_row_count",
    "local_profile_drug_therapy_type_values_json",
    "local_profile_regimen_context_values_json",
    "local_profile_has_any_treatment_timing",
    "local_profile_status",
    "local_profile_requires_manual_review",
    "local_grouping_treatment_group_v1",
    "local_grouping_treatment_group_v1_rule",
    "local_grouping_requires_manual_review",
    "local_has_any_treatment_evidence",
    "indexed_diagnosis_count",
    "indexed_treatment_row_count",
    "indexed_follow_up_row_count",
    "indexed_treatment_type_values_json",
    "indexed_treatment_type_distinct_count",
    "indexed_treatment_intent_type_values_json",
    "indexed_therapeutic_agents_values_json",
    "indexed_regimen_line_values_json",
    "indexed_followup_timepoint_category_values_json",
    "indexed_has_any_treatment_row",
    "indexed_has_any_druglike_treatment",
    "indexed_has_any_radiation_treatment",
    "indexed_has_any_surgery_treatment",
    "indexed_has_any_non_surgical_treatment",
    "indexed_has_any_therapeutic_agent_value",
    "indexed_has_any_regimen_line_value",
    "indexed_has_any_treatment_timing",
    "indexed_has_prior_treatment_flag",
    "indexed_has_followup_post_initial_treatment",
    "indexed_has_followup_explicit_treatment_field",
    "indexed_adds_missing_druglike_coverage",
    "indexed_adds_missing_radiation_coverage",
    "indexed_adds_agent_detail_to_locally_missing_patient",
    "indexed_adds_regimen_line_detail_to_locally_missing_patient",
    "indexed_adds_timing_detail_to_locally_missing_patient",
    "indexed_surgery_only_flag",
    "indexed_context_only_flag",
    "indexed_vs_local_any_treatment_class",
    "indexed_incremental_value_class",
]
OVERLAP_FIELDNAMES = [
    "indexed_clinical_vs_local_treatment_v1_run_id",
    "comparison_scope",
    "dimension",
    "local_status",
    "indexed_status",
    "overlap_class",
    "patient_count",
    "patient_fraction_of_local_profile",
    "notes",
]
GAP_RESOLUTION_FIELDNAMES = [
    "indexed_clinical_vs_local_treatment_v1_run_id",
    "gap_group",
    "denominator_patient_count",
    "gap_metric",
    "patient_count",
    "patient_fraction_of_gap_group",
    "notes",
]
SUMMARY_FIELDNAMES = [
    "indexed_clinical_vs_local_treatment_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]


class IndexedClinicalVsLocalTreatmentV1Error(RuntimeError):
    """Raised when the indexed-clinical comparison workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    repo_root: Path
    trial_config: Path
    audit_root: Path
    source_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    clinical_biotabs_latest_pointer: Path
    patient_treatment_profile_latest_pointer: Path
    patient_treatment_grouping_latest_pointer: Path
    treatment_source_audit_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    clinical_biotabs_pointer: dict[str, Any]
    clinical_biotabs_run_log: dict[str, Any]
    clinical_drug_rows: list[dict[str, str]]
    clinical_radiation_rows: list[dict[str, str]]
    clinical_patient_rows: list[dict[str, str]]
    patient_treatment_profile_pointer: dict[str, Any]
    patient_treatment_profile_run_log: dict[str, Any]
    patient_treatment_profile_rows: list[dict[str, str]]
    patient_treatment_grouping_pointer: dict[str, Any]
    patient_treatment_grouping_run_log: dict[str, Any]
    patient_treatment_grouping_rows: list[dict[str, str]]
    treatment_source_audit_pointer: dict[str, Any]
    treatment_source_audit_run_log: dict[str, Any]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class IndexedEvidenceOutputs:
    gdc_status_payload: dict[str, Any]
    cases_combined_payload: dict[str, Any]
    cases_request_manifest_path: str
    cases_combined_response_path: str
    page_request_paths: list[str]
    page_response_paths: list[str]
    fields: list[str]
    filters: dict[str, Any]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise IndexedClinicalVsLocalTreatmentV1Error(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise IndexedClinicalVsLocalTreatmentV1Error(f"Unable to create an import spec for: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def normalize_barcode(value: str) -> str:
    return value.strip().upper()


def normalize_token(value: Any) -> str:
    return str(value or "").strip().lower()


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def format_fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.0000"
    return f"{numerator / denominator:.4f}"


def json_list(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=True)


def ordered_unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def is_local_missing_like(value: Any) -> bool:
    return normalize_token(value) in LOCAL_MISSING_LIKE_TOKENS


def is_indexed_missing_like(value: Any) -> bool:
    return normalize_token(value) in INDEXED_MISSING_LIKE_TOKENS


def collect_ordered_non_missing_values(
    raw_values: list[Any],
    *,
    missing_like_fn,
) -> list[str]:
    values: list[str] = []
    for raw_value in raw_values:
        stripped = str(raw_value or "").strip()
        if missing_like_fn(stripped):
            continue
        values.append(stripped)
    return ordered_unique(values)


def build_unique_lookup_by_barcode(
    rows: list[dict[str, str]],
    *,
    barcode_field: str,
    label: str,
) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        barcode = normalize_barcode(row.get(barcode_field, ""))
        if not barcode:
            raise IndexedClinicalVsLocalTreatmentV1Error(f"{label} contains an empty barcode value.")
        if barcode in lookup:
            raise IndexedClinicalVsLocalTreatmentV1Error(f"{label} contains a duplicate barcode: {barcode}")
        lookup[barcode] = row
    return lookup


def group_rows_by_barcode(
    rows: list[dict[str, str]],
    *,
    barcode_field: str,
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        barcode = normalize_barcode(row.get(barcode_field, ""))
        if barcode:
            grouped[barcode].append(row)
    return grouped


def parse_json_list(raw_value: str) -> list[str]:
    stripped = str(raw_value).strip()
    if not stripped:
        return []
    parsed = json.loads(stripped)
    if not isinstance(parsed, list):
        raise IndexedClinicalVsLocalTreatmentV1Error(f"Expected JSON list value, got: {raw_value!r}")
    return [str(value) for value in parsed]


def require_columns(rows: list[dict[str, str]], required_columns: set[str], label: str) -> None:
    if not rows:
        raise IndexedClinicalVsLocalTreatmentV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise IndexedClinicalVsLocalTreatmentV1Error(f"{label} is missing required columns: {sorted(missing)}")


def require_completed_pointer_and_run_log(
    pointer_path: Path,
    label: str,
    helper_module: Any,
    repo_root: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    if not pointer_path.exists():
        raise IndexedClinicalVsLocalTreatmentV1Error(f"Required {label} pointer not found: {pointer_path}")
    pointer = helper_module.load_json(pointer_path)
    run_log_relative = str(pointer.get("run_log_json") or "").strip()
    if not run_log_relative:
        raise IndexedClinicalVsLocalTreatmentV1Error(f"{label} pointer is missing run_log_json: {pointer_path}")
    run_log_path = helper_module.resolve_existing_path(repo_root, run_log_relative, f"{label} run log")
    run_log = helper_module.load_json(run_log_path)
    if str(run_log.get("status") or "") != "completed":
        raise IndexedClinicalVsLocalTreatmentV1Error(f"{label} run log is not completed: {run_log_path}")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise IndexedClinicalVsLocalTreatmentV1Error(
            f"{label} run log does not report validation.passed == true: {run_log_path}"
        )
    return pointer, run_log


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    if not trial_config.exists():
        raise IndexedClinicalVsLocalTreatmentV1Error(f"Required trial config not found: {trial_config}")
    trial_config_data = helper_module.load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    source_root = audit_root / "tcga-brca" / "source"
    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        audit_root=audit_root,
        source_root=source_root,
        audit_runs_root=source_root / "indexed_clinical_vs_local_treatment_v1_runs",
        latest_pointer=source_root / "tcga_brca_indexed_clinical_vs_local_treatment_v1_latest.json",
        clinical_biotabs_latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_biotabs_latest.json",
        patient_treatment_profile_latest_pointer=(
            audit_root / "tcga-brca" / "treatment-prep" / "tcga_brca_patient_treatment_profile_v1_latest.json"
        ),
        patient_treatment_grouping_latest_pointer=(
            audit_root / "tcga-brca" / "treatment-prep" / "tcga_brca_patient_treatment_grouping_v1_latest.json"
        ),
        treatment_source_audit_latest_pointer=source_root / "tcga_brca_treatment_source_audit_latest.json",
    )


def load_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    clinical_biotabs_pointer, clinical_biotabs_run_log = require_completed_pointer_and_run_log(
        paths.clinical_biotabs_latest_pointer,
        "clinical biotabs latest",
        helper_module,
        paths.repo_root,
    )
    patient_treatment_profile_pointer, patient_treatment_profile_run_log = require_completed_pointer_and_run_log(
        paths.patient_treatment_profile_latest_pointer,
        "patient treatment profile latest",
        helper_module,
        paths.repo_root,
    )
    patient_treatment_grouping_pointer, patient_treatment_grouping_run_log = require_completed_pointer_and_run_log(
        paths.patient_treatment_grouping_latest_pointer,
        "patient treatment grouping latest",
        helper_module,
        paths.repo_root,
    )
    treatment_source_audit_pointer, treatment_source_audit_run_log = require_completed_pointer_and_run_log(
        paths.treatment_source_audit_latest_pointer,
        "treatment source audit latest",
        helper_module,
        paths.repo_root,
    )

    helper_module.require_keys(
        clinical_biotabs_pointer,
        {"parse_run_id", "processed_run_directory", "run_log_json"},
        "clinical biotabs latest pointer",
        paths.clinical_biotabs_latest_pointer,
    )
    helper_module.require_keys(
        patient_treatment_profile_pointer,
        {
            "patient_treatment_profile_v1_run_id",
            "patient_treatment_profile_v1_tsv",
            "run_log_json",
        },
        "patient treatment profile latest pointer",
        paths.patient_treatment_profile_latest_pointer,
    )
    helper_module.require_keys(
        patient_treatment_grouping_pointer,
        {
            "patient_treatment_grouping_v1_run_id",
            "patient_treatment_grouping_v1_tsv",
            "run_log_json",
        },
        "patient treatment grouping latest pointer",
        paths.patient_treatment_grouping_latest_pointer,
    )
    helper_module.require_keys(
        treatment_source_audit_pointer,
        {
            "treatment_source_audit_run_id",
            "run_log_json",
            "gdc_data_release",
            "gdc_tag",
        },
        "treatment source audit latest pointer",
        paths.treatment_source_audit_latest_pointer,
    )

    clinical_processed_dir = Path(str(clinical_biotabs_pointer["processed_run_directory"]))
    clinical_drug_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (clinical_processed_dir / "clinical_drug.tsv").as_posix(),
        "clinical_drug.tsv",
    )
    clinical_radiation_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (clinical_processed_dir / "clinical_radiation.tsv").as_posix(),
        "clinical_radiation.tsv",
    )
    clinical_patient_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (clinical_processed_dir / "clinical_patient.tsv").as_posix(),
        "clinical_patient.tsv",
    )
    patient_treatment_profile_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(patient_treatment_profile_pointer["patient_treatment_profile_v1_tsv"]),
        "patient_treatment_profile_v1.tsv",
    )
    patient_treatment_grouping_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(patient_treatment_grouping_pointer["patient_treatment_grouping_v1_tsv"]),
        "patient_treatment_grouping_v1.tsv",
    )

    clinical_drug_rows = helper_module.read_tsv_dict_rows(clinical_drug_tsv)
    clinical_radiation_rows = helper_module.read_tsv_dict_rows(clinical_radiation_tsv)
    clinical_patient_rows = helper_module.read_tsv_dict_rows(clinical_patient_tsv)
    patient_treatment_profile_rows = helper_module.read_tsv_dict_rows(patient_treatment_profile_tsv)
    patient_treatment_grouping_rows = helper_module.read_tsv_dict_rows(patient_treatment_grouping_tsv)

    require_columns(
        clinical_drug_rows,
        {"bcr_patient_barcode", "bcr_patient_uuid", "pharmaceutical_therapy_type"},
        "clinical_drug.tsv",
    )
    require_columns(
        clinical_radiation_rows,
        {"bcr_patient_barcode", "bcr_patient_uuid"},
        "clinical_radiation.tsv",
    )
    require_columns(
        clinical_patient_rows,
        {"bcr_patient_barcode", "bcr_patient_uuid"},
        "clinical_patient.tsv",
    )
    require_columns(
        patient_treatment_profile_rows,
        {
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "has_any_drug_row",
            "drug_row_count",
            "has_any_radiation_row",
            "radiation_row_count",
            "drug_therapy_type_values_json",
            "regimen_context_values_json",
            "has_any_treatment_timing",
            "treatment_profile_requires_manual_review",
            "treatment_profile_status",
        },
        "patient_treatment_profile_v1.tsv",
    )
    require_columns(
        patient_treatment_grouping_rows,
        {
            "bcr_patient_barcode",
            "treatment_group_v1",
            "treatment_group_v1_rule",
            "treatment_group_v1_requires_manual_review",
        },
        "patient_treatment_grouping_v1.tsv",
    )

    return WorkflowInputs(
        clinical_biotabs_pointer=clinical_biotabs_pointer,
        clinical_biotabs_run_log=clinical_biotabs_run_log,
        clinical_drug_rows=clinical_drug_rows,
        clinical_radiation_rows=clinical_radiation_rows,
        clinical_patient_rows=clinical_patient_rows,
        patient_treatment_profile_pointer=patient_treatment_profile_pointer,
        patient_treatment_profile_run_log=patient_treatment_profile_run_log,
        patient_treatment_profile_rows=patient_treatment_profile_rows,
        patient_treatment_grouping_pointer=patient_treatment_grouping_pointer,
        patient_treatment_grouping_run_log=patient_treatment_grouping_run_log,
        patient_treatment_grouping_rows=patient_treatment_grouping_rows,
        treatment_source_audit_pointer=treatment_source_audit_pointer,
        treatment_source_audit_run_log=treatment_source_audit_run_log,
        input_paths={
            "clinical_drug_tsv": clinical_drug_tsv,
            "clinical_radiation_tsv": clinical_radiation_tsv,
            "clinical_patient_tsv": clinical_patient_tsv,
            "patient_treatment_profile_v1_tsv": patient_treatment_profile_tsv,
            "patient_treatment_grouping_v1_tsv": patient_treatment_grouping_tsv,
        },
    )


def save_json_snapshot(helper_module: Any, repo_root: Path, path: Path, payload: Any) -> str:
    helper_module.write_json(path, payload)
    return repo_relative(path, repo_root)


def gdc_get_json(
    session: requests.Session,
    url: str,
    *,
    params: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    response = session.get(url, params=params, timeout=REQUEST_TIMEOUT_SECONDS)
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        preview = response.text[:500]
        raise IndexedClinicalVsLocalTreatmentV1Error(
            f"GDC request failed for {url} with HTTP {response.status_code}. Body preview: {preview}"
        ) from exc
    try:
        payload = response.json()
    except json.JSONDecodeError as exc:
        raise IndexedClinicalVsLocalTreatmentV1Error(f"Unable to decode JSON response from {url}") from exc
    return payload, {"final_url": response.url, "params": params or {}}


def fetch_indexed_cases(
    audit_run_dir: Path,
    paths: WorkflowPaths,
    helper_module: Any,
) -> IndexedEvidenceOutputs:
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})

    status_payload, status_request_meta = gdc_get_json(session, f"{GDC_API_BASE_URL}/status")
    save_json_snapshot(helper_module, paths.repo_root, audit_run_dir / "gdc_status_request.json", status_request_meta)
    save_json_snapshot(helper_module, paths.repo_root, audit_run_dir / "gdc_status_response.json", status_payload)

    fields = [
        "submitter_id",
        "case_id",
        "diagnoses.submitter_id",
        "diagnoses.prior_treatment",
        "diagnoses.treatments.submitter_id",
        "diagnoses.treatments.treatment_type",
        "diagnoses.treatments.treatment_or_therapy",
        "diagnoses.treatments.treatment_intent_type",
        "diagnoses.treatments.therapeutic_agents",
        "diagnoses.treatments.regimen_or_line_of_therapy",
        "diagnoses.treatments.days_to_treatment_start",
        "diagnoses.treatments.days_to_treatment_end",
        "diagnoses.treatments.reason_treatment_not_given",
        "follow_ups.submitter_id",
        "follow_ups.timepoint_category",
        "follow_ups.treatment_emergent_adverse_event",
        "follow_ups.hormone_replacement_therapy_type",
    ]
    filters = {"op": "in", "content": {"field": "project.project_id", "value": [PROJECT_ID]}}
    request_manifest = {
        "endpoint": f"{GDC_API_BASE_URL}/cases",
        "project_id": PROJECT_ID,
        "filters": filters,
        "fields": fields,
        "sort": "submitter_id:asc",
        "page_size": INDEXED_CASE_PAGE_SIZE,
    }
    request_manifest_path = save_json_snapshot(
        helper_module,
        paths.repo_root,
        audit_run_dir / "indexed_cases_request_manifest.json",
        request_manifest,
    )

    all_hits: list[dict[str, Any]] = []
    page_request_paths: list[str] = []
    page_response_paths: list[str] = []
    warnings: list[Any] = []
    seen_barcodes: set[str] = set()
    pagination_total: int | None = None
    page_index = 0
    offset = 0

    while True:
        params = {
            "filters": json.dumps(filters),
            "fields": ",".join(fields),
            "sort": "submitter_id:asc",
            "size": INDEXED_CASE_PAGE_SIZE,
            "from": offset,
        }
        page_payload, page_request_meta = gdc_get_json(session, f"{GDC_API_BASE_URL}/cases", params=params)
        page_name = f"indexed_cases_page_{page_index + 1:04d}"
        page_request_paths.append(
            save_json_snapshot(helper_module, paths.repo_root, audit_run_dir / f"{page_name}_request.json", page_request_meta)
        )
        page_response_paths.append(
            save_json_snapshot(helper_module, paths.repo_root, audit_run_dir / f"{page_name}_response.json", page_payload)
        )

        hits = page_payload.get("data", {}).get("hits", [])
        if not isinstance(hits, list):
            raise IndexedClinicalVsLocalTreatmentV1Error("Indexed cases payload is missing a list of hits.")
        pagination = page_payload.get("data", {}).get("pagination", {})
        if pagination_total is None:
            pagination_total = int(pagination.get("total", 0))
            if pagination_total <= 0:
                raise IndexedClinicalVsLocalTreatmentV1Error("Indexed cases pagination.total must be positive.")
        for hit in hits:
            barcode = normalize_barcode(hit.get("submitter_id", ""))
            if not barcode:
                raise IndexedClinicalVsLocalTreatmentV1Error("Indexed cases payload contains an empty submitter_id.")
            if barcode in seen_barcodes:
                raise IndexedClinicalVsLocalTreatmentV1Error(f"Indexed cases payload contains a duplicate submitter_id: {barcode}")
            seen_barcodes.add(barcode)
            all_hits.append(hit)
        warnings.extend(page_payload.get("warnings", []))

        if len(all_hits) >= pagination_total:
            break
        if not hits:
            raise IndexedClinicalVsLocalTreatmentV1Error(
                "Indexed cases pagination ended early before reaching pagination.total."
            )
        offset += len(hits)
        page_index += 1

    combined_payload = {
        "data": {
            "pagination": {
                "count": len(all_hits),
                "from": 0,
                "page_size": INDEXED_CASE_PAGE_SIZE,
                "sort": "submitter_id:asc",
                "total": pagination_total,
            },
            "hits": all_hits,
        },
        "warnings": warnings,
    }
    combined_response_path = save_json_snapshot(
        helper_module,
        paths.repo_root,
        audit_run_dir / "indexed_cases_combined_response.json",
        combined_payload,
    )
    return IndexedEvidenceOutputs(
        gdc_status_payload=status_payload,
        cases_combined_payload=combined_payload,
        cases_request_manifest_path=request_manifest_path,
        cases_combined_response_path=combined_response_path,
        page_request_paths=page_request_paths,
        page_response_paths=page_response_paths,
        fields=fields,
        filters=filters,
    )


def indexed_treatment_row_is_explicit(treatment_row: dict[str, Any]) -> bool:
    if treatment_row.get("days_to_treatment_start") is not None:
        return True
    if treatment_row.get("days_to_treatment_end") is not None:
        return True
    for field_name in [
        "treatment_type",
        "treatment_or_therapy",
        "treatment_intent_type",
        "therapeutic_agents",
        "regimen_or_line_of_therapy",
        "reason_treatment_not_given",
    ]:
        if not is_indexed_missing_like(treatment_row.get(field_name)):
            return True
    return False


def is_indexed_druglike_treatment_type(raw_type: str) -> bool:
    return normalize_token(raw_type) in INDEXED_DRUGLIKE_TREATMENT_TYPES


def is_indexed_radiation_treatment_type(raw_type: str) -> bool:
    normalized = normalize_token(raw_type)
    return any(token in normalized for token in INDEXED_RADIATION_TREATMENT_TYPE_TOKENS)


def is_indexed_surgery_treatment_type(raw_type: str) -> bool:
    return "surgery" in normalize_token(raw_type)


def summarize_indexed_case(case_hit: dict[str, Any]) -> dict[str, Any]:
    case_id = str(case_hit.get("case_id") or "").strip()
    diagnoses = case_hit.get("diagnoses") or []
    follow_ups = case_hit.get("follow_ups") or []
    if not isinstance(diagnoses, list):
        raise IndexedClinicalVsLocalTreatmentV1Error("Indexed case diagnoses must be a list.")
    if not isinstance(follow_ups, list):
        raise IndexedClinicalVsLocalTreatmentV1Error("Indexed case follow_ups must be a list.")

    explicit_treatment_rows: list[dict[str, Any]] = []
    treatment_type_raw_values: list[Any] = []
    treatment_intent_raw_values: list[Any] = []
    therapeutic_agents_raw_values: list[Any] = []
    regimen_line_raw_values: list[Any] = []
    prior_treatment_values: list[str] = []

    for diagnosis in diagnoses:
        if not isinstance(diagnosis, dict):
            continue
        prior_treatment_values.append(str(diagnosis.get("prior_treatment") or ""))
        treatment_rows = diagnosis.get("treatments") or []
        if not isinstance(treatment_rows, list):
            raise IndexedClinicalVsLocalTreatmentV1Error("Indexed case diagnosis.treatments must be a list.")
        for treatment_row in treatment_rows:
            if not isinstance(treatment_row, dict):
                continue
            if not indexed_treatment_row_is_explicit(treatment_row):
                continue
            explicit_treatment_rows.append(treatment_row)
            treatment_type_raw_values.append(treatment_row.get("treatment_type"))
            treatment_intent_raw_values.append(treatment_row.get("treatment_intent_type"))
            therapeutic_agents_raw_values.append(treatment_row.get("therapeutic_agents"))
            regimen_line_raw_values.append(treatment_row.get("regimen_or_line_of_therapy"))

    followup_timepoint_categories = collect_ordered_non_missing_values(
        [follow_up.get("timepoint_category") for follow_up in follow_ups if isinstance(follow_up, dict)],
        missing_like_fn=is_indexed_missing_like,
    )
    has_followup_explicit_treatment_field = False
    for follow_up in follow_ups:
        if not isinstance(follow_up, dict):
            continue
        if not is_indexed_missing_like(follow_up.get("treatment_emergent_adverse_event")):
            has_followup_explicit_treatment_field = True
            break
        if not is_indexed_missing_like(follow_up.get("hormone_replacement_therapy_type")):
            has_followup_explicit_treatment_field = True
            break

    treatment_type_values = collect_ordered_non_missing_values(
        treatment_type_raw_values,
        missing_like_fn=is_indexed_missing_like,
    )
    treatment_intent_values = collect_ordered_non_missing_values(
        treatment_intent_raw_values,
        missing_like_fn=is_indexed_missing_like,
    )
    therapeutic_agent_values = collect_ordered_non_missing_values(
        therapeutic_agents_raw_values,
        missing_like_fn=is_indexed_missing_like,
    )
    regimen_line_values = collect_ordered_non_missing_values(
        regimen_line_raw_values,
        missing_like_fn=is_indexed_missing_like,
    )
    normalized_treatment_types = [normalize_token(value) for value in treatment_type_values]

    has_any_druglike_treatment = bool(therapeutic_agent_values) or any(
        is_indexed_druglike_treatment_type(value) for value in treatment_type_values
    )
    has_any_radiation_treatment = any(is_indexed_radiation_treatment_type(value) for value in treatment_type_values)
    has_any_surgery_treatment = any(is_indexed_surgery_treatment_type(value) for value in treatment_type_values)
    has_any_non_surgical_treatment = bool(
        has_any_druglike_treatment
        or has_any_radiation_treatment
        or any(
            normalized_value and not is_indexed_surgery_treatment_type(normalized_value)
            for normalized_value in normalized_treatment_types
        )
    )
    has_any_treatment_timing = any(
        treatment_row.get("days_to_treatment_start") is not None
        or treatment_row.get("days_to_treatment_end") is not None
        for treatment_row in explicit_treatment_rows
    )
    has_prior_treatment_flag = any(normalize_token(value) == "yes" for value in prior_treatment_values)
    has_followup_post_initial_treatment = any(
        normalize_token(value) == "post initial treatment" for value in followup_timepoint_categories
    )

    return {
        "indexed_case_id": case_id,
        "indexed_diagnosis_count": len(diagnoses),
        "indexed_treatment_row_count": len(explicit_treatment_rows),
        "indexed_follow_up_row_count": len(follow_ups),
        "indexed_treatment_type_values": treatment_type_values,
        "indexed_treatment_intent_type_values": treatment_intent_values,
        "indexed_therapeutic_agents_values": therapeutic_agent_values,
        "indexed_regimen_line_values": regimen_line_values,
        "indexed_followup_timepoint_category_values": followup_timepoint_categories,
        "indexed_has_any_treatment_row": bool(explicit_treatment_rows),
        "indexed_has_any_druglike_treatment": has_any_druglike_treatment,
        "indexed_has_any_radiation_treatment": has_any_radiation_treatment,
        "indexed_has_any_surgery_treatment": has_any_surgery_treatment,
        "indexed_has_any_non_surgical_treatment": has_any_non_surgical_treatment,
        "indexed_has_any_therapeutic_agent_value": bool(therapeutic_agent_values),
        "indexed_has_any_regimen_line_value": bool(regimen_line_values),
        "indexed_has_any_treatment_timing": has_any_treatment_timing,
        "indexed_has_prior_treatment_flag": has_prior_treatment_flag,
        "indexed_has_followup_post_initial_treatment": has_followup_post_initial_treatment,
        "indexed_has_followup_explicit_treatment_field": has_followup_explicit_treatment_field,
    }


def build_indexed_lookup(indexed_hits: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    indexed_lookup: dict[str, dict[str, Any]] = {}
    for case_hit in indexed_hits:
        barcode = normalize_barcode(case_hit.get("submitter_id", ""))
        if not barcode:
            raise IndexedClinicalVsLocalTreatmentV1Error("Indexed case hit is missing submitter_id.")
        if barcode in indexed_lookup:
            raise IndexedClinicalVsLocalTreatmentV1Error(f"Duplicate indexed case submitter_id detected: {barcode}")
        indexed_lookup[barcode] = summarize_indexed_case(case_hit)
    return indexed_lookup


def build_patient_level_rows(
    run_id: str,
    inputs: WorkflowInputs,
    indexed_lookup: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, int], dict[str, list[str]], dict[str, Any]]:
    clinical_patient_lookup = build_unique_lookup_by_barcode(
        inputs.clinical_patient_rows,
        barcode_field="bcr_patient_barcode",
        label="clinical_patient.tsv",
    )
    profile_lookup = build_unique_lookup_by_barcode(
        inputs.patient_treatment_profile_rows,
        barcode_field="bcr_patient_barcode",
        label="patient_treatment_profile_v1.tsv",
    )
    grouping_lookup = build_unique_lookup_by_barcode(
        inputs.patient_treatment_grouping_rows,
        barcode_field="bcr_patient_barcode",
        label="patient_treatment_grouping_v1.tsv",
    )
    drug_rows_by_barcode = group_rows_by_barcode(inputs.clinical_drug_rows, barcode_field="bcr_patient_barcode")
    radiation_rows_by_barcode = group_rows_by_barcode(
        inputs.clinical_radiation_rows,
        barcode_field="bcr_patient_barcode",
    )

    profile_order = [normalize_barcode(row["bcr_patient_barcode"]) for row in inputs.patient_treatment_profile_rows]
    grouped_order = [normalize_barcode(row["bcr_patient_barcode"]) for row in inputs.patient_treatment_grouping_rows]
    if profile_order != grouped_order:
        raise IndexedClinicalVsLocalTreatmentV1Error(
            "Patient treatment grouping row order must match the patient treatment profile row order."
        )

    indexed_only_barcodes = sorted(set(indexed_lookup).difference(profile_lookup))
    patient_order = profile_order + indexed_only_barcodes

    patient_rows: list[dict[str, Any]] = []
    counts = Counter()
    category_lists: dict[str, list[str]] = {
        "indexed_treatment_types": [],
        "local_profile_drug_therapy_types": [],
    }

    for barcode in patient_order:
        profile_row = profile_lookup.get(barcode)
        grouping_row = grouping_lookup.get(barcode)
        clinical_patient_row = clinical_patient_lookup.get(barcode)
        indexed_row = indexed_lookup.get(barcode)
        local_drug_rows = drug_rows_by_barcode.get(barcode, [])
        local_radiation_rows = radiation_rows_by_barcode.get(barcode, [])

        if profile_row is not None and clinical_patient_row is None:
            raise IndexedClinicalVsLocalTreatmentV1Error(
                f"Profile barcode {barcode} is missing from clinical_patient.tsv."
            )
        if profile_row is not None and grouping_row is None:
            raise IndexedClinicalVsLocalTreatmentV1Error(
                f"Profile barcode {barcode} is missing from patient_treatment_grouping_v1.tsv."
            )
        if profile_row is not None and indexed_row is None:
            raise IndexedClinicalVsLocalTreatmentV1Error(
                f"Profile barcode {barcode} is missing from the indexed clinical cases pull."
            )

        local_clinical_any_drug_row = bool(local_drug_rows)
        local_clinical_any_radiation_row = bool(local_radiation_rows)
        local_profile_has_any_drug_row = profile_row is not None and profile_row["has_any_drug_row"] == "yes"
        local_profile_has_any_radiation_row = profile_row is not None and profile_row["has_any_radiation_row"] == "yes"
        if profile_row is not None:
            if local_clinical_any_drug_row != local_profile_has_any_drug_row:
                raise IndexedClinicalVsLocalTreatmentV1Error(
                    f"Local drug-row flag mismatch between clinical_drug.tsv and patient profile for {barcode}."
                )
            if len(local_drug_rows) != int(profile_row["drug_row_count"]):
                raise IndexedClinicalVsLocalTreatmentV1Error(
                    f"Local drug-row count mismatch between clinical_drug.tsv and patient profile for {barcode}."
                )
            if local_clinical_any_radiation_row != local_profile_has_any_radiation_row:
                raise IndexedClinicalVsLocalTreatmentV1Error(
                    f"Local radiation-row flag mismatch between clinical_radiation.tsv and patient profile for {barcode}."
                )
            if len(local_radiation_rows) != int(profile_row["radiation_row_count"]):
                raise IndexedClinicalVsLocalTreatmentV1Error(
                    f"Local radiation-row count mismatch between clinical_radiation.tsv and patient profile for {barcode}."
                )

        local_profile_drug_therapy_type_values = parse_json_list(profile_row["drug_therapy_type_values_json"]) if profile_row is not None else []
        local_clinical_drug_therapy_type_values = [str(row.get("pharmaceutical_therapy_type", "")) for row in local_drug_rows]
        if profile_row is not None and local_profile_drug_therapy_type_values != local_clinical_drug_therapy_type_values:
            raise IndexedClinicalVsLocalTreatmentV1Error(
                f"Local drug therapy type values mismatch between clinical_drug.tsv and patient profile for {barcode}."
            )

        local_has_any_treatment_evidence = local_profile_has_any_drug_row or local_profile_has_any_radiation_row
        indexed_has_any_treatment_row = bool(indexed_row and indexed_row["indexed_has_any_treatment_row"])
        indexed_has_any_druglike_treatment = bool(indexed_row and indexed_row["indexed_has_any_druglike_treatment"])
        indexed_has_any_radiation_treatment = bool(indexed_row and indexed_row["indexed_has_any_radiation_treatment"])
        indexed_has_any_therapeutic_agent_value = bool(indexed_row and indexed_row["indexed_has_any_therapeutic_agent_value"])
        indexed_has_any_regimen_line_value = bool(indexed_row and indexed_row["indexed_has_any_regimen_line_value"])
        indexed_has_any_treatment_timing = bool(indexed_row and indexed_row["indexed_has_any_treatment_timing"])
        indexed_has_any_surgery_treatment = bool(indexed_row and indexed_row["indexed_has_any_surgery_treatment"])
        indexed_has_any_non_surgical_treatment = bool(indexed_row and indexed_row["indexed_has_any_non_surgical_treatment"])
        indexed_has_prior_treatment_flag = bool(indexed_row and indexed_row["indexed_has_prior_treatment_flag"])
        indexed_has_followup_post_initial_treatment = bool(indexed_row and indexed_row["indexed_has_followup_post_initial_treatment"])
        indexed_has_followup_explicit_treatment_field = bool(
            indexed_row and indexed_row["indexed_has_followup_explicit_treatment_field"]
        )

        indexed_adds_missing_druglike_coverage = indexed_has_any_druglike_treatment and not local_profile_has_any_drug_row
        indexed_adds_missing_radiation_coverage = indexed_has_any_radiation_treatment and not local_profile_has_any_radiation_row
        indexed_adds_agent_detail_to_locally_missing_patient = not local_has_any_treatment_evidence and indexed_has_any_therapeutic_agent_value
        indexed_adds_regimen_line_detail_to_locally_missing_patient = not local_has_any_treatment_evidence and indexed_has_any_regimen_line_value
        indexed_adds_timing_detail_to_locally_missing_patient = not local_has_any_treatment_evidence and indexed_has_any_treatment_timing
        indexed_surgery_only_flag = indexed_has_any_treatment_row and indexed_has_any_surgery_treatment and not indexed_has_any_non_surgical_treatment
        indexed_context_only_flag = (
            not indexed_has_any_treatment_row
            and (
                indexed_has_prior_treatment_flag
                or indexed_has_followup_post_initial_treatment
                or indexed_has_followup_explicit_treatment_field
            )
        )
        indexed_adds_useful_new_evidence = bool(
            indexed_adds_missing_druglike_coverage
            or indexed_adds_missing_radiation_coverage
            or indexed_adds_agent_detail_to_locally_missing_patient
            or indexed_adds_regimen_line_detail_to_locally_missing_patient
            or indexed_adds_timing_detail_to_locally_missing_patient
        )

        if local_has_any_treatment_evidence and indexed_has_any_treatment_row:
            indexed_vs_local_any_treatment_class = "both"
        elif local_has_any_treatment_evidence and not indexed_has_any_treatment_row:
            indexed_vs_local_any_treatment_class = "local_only"
        elif not local_has_any_treatment_evidence and indexed_has_any_treatment_row:
            indexed_vs_local_any_treatment_class = "indexed_only"
        else:
            indexed_vs_local_any_treatment_class = "neither"

        if indexed_adds_useful_new_evidence:
            indexed_incremental_value_class = "potentially_useful_new_evidence"
        elif (
            indexed_has_any_treatment_row
            or indexed_context_only_flag
            or indexed_has_prior_treatment_flag
            or indexed_has_followup_post_initial_treatment
            or indexed_has_followup_explicit_treatment_field
        ):
            indexed_incremental_value_class = "coarse_only"
        else:
            indexed_incremental_value_class = "none"

        if profile_row is not None and indexed_row is not None:
            cohort_join_status = "local_profile_and_indexed_case"
            cohort_scope = "local_profile"
        elif profile_row is None and indexed_row is not None:
            cohort_join_status = "indexed_case_only_not_in_local_profile"
            cohort_scope = "indexed_only"
        else:
            cohort_join_status = "local_profile_only_missing_indexed_case"
            cohort_scope = "local_profile"

        row = {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "treatment_source_audit_run_id": str(inputs.treatment_source_audit_pointer["treatment_source_audit_run_id"]),
            "patient_treatment_profile_v1_run_id": str(inputs.patient_treatment_profile_pointer["patient_treatment_profile_v1_run_id"]),
            "patient_treatment_grouping_v1_run_id": str(inputs.patient_treatment_grouping_pointer["patient_treatment_grouping_v1_run_id"]),
            "clinical_biotab_parse_run_id": str(inputs.clinical_biotabs_pointer["parse_run_id"]),
            "cohort_join_status": cohort_join_status,
            "cohort_scope": cohort_scope,
            "bcr_patient_barcode": barcode,
            "bcr_patient_uuid": str((clinical_patient_row or {}).get("bcr_patient_uuid", "")),
            "indexed_case_id": str((indexed_row or {}).get("indexed_case_id", "")),
            "local_clinical_any_drug_row": yes_no(local_clinical_any_drug_row),
            "local_clinical_drug_row_count": str(len(local_drug_rows)),
            "local_clinical_any_radiation_row": yes_no(local_clinical_any_radiation_row),
            "local_clinical_radiation_row_count": str(len(local_radiation_rows)),
            "local_clinical_drug_therapy_type_values_json": json_list(local_clinical_drug_therapy_type_values),
            "local_profile_has_any_drug_row": str((profile_row or {}).get("has_any_drug_row", "no")),
            "local_profile_drug_row_count": str((profile_row or {}).get("drug_row_count", "0")),
            "local_profile_has_any_radiation_row": str((profile_row or {}).get("has_any_radiation_row", "no")),
            "local_profile_radiation_row_count": str((profile_row or {}).get("radiation_row_count", "0")),
            "local_profile_drug_therapy_type_values_json": str((profile_row or {}).get("drug_therapy_type_values_json", "[]")),
            "local_profile_regimen_context_values_json": str((profile_row or {}).get("regimen_context_values_json", "[]")),
            "local_profile_has_any_treatment_timing": str((profile_row or {}).get("has_any_treatment_timing", "no")),
            "local_profile_status": str((profile_row or {}).get("treatment_profile_status", "")),
            "local_profile_requires_manual_review": str((profile_row or {}).get("treatment_profile_requires_manual_review", "no")),
            "local_grouping_treatment_group_v1": str((grouping_row or {}).get("treatment_group_v1", "")),
            "local_grouping_treatment_group_v1_rule": str((grouping_row or {}).get("treatment_group_v1_rule", "")),
            "local_grouping_requires_manual_review": str((grouping_row or {}).get("treatment_group_v1_requires_manual_review", "no")),
            "local_has_any_treatment_evidence": yes_no(local_has_any_treatment_evidence),
            "indexed_diagnosis_count": str((indexed_row or {}).get("indexed_diagnosis_count", 0)),
            "indexed_treatment_row_count": str((indexed_row or {}).get("indexed_treatment_row_count", 0)),
            "indexed_follow_up_row_count": str((indexed_row or {}).get("indexed_follow_up_row_count", 0)),
            "indexed_treatment_type_values_json": json_list(list((indexed_row or {}).get("indexed_treatment_type_values", []))),
            "indexed_treatment_type_distinct_count": str(len((indexed_row or {}).get("indexed_treatment_type_values", []))),
            "indexed_treatment_intent_type_values_json": json_list(list((indexed_row or {}).get("indexed_treatment_intent_type_values", []))),
            "indexed_therapeutic_agents_values_json": json_list(list((indexed_row or {}).get("indexed_therapeutic_agents_values", []))),
            "indexed_regimen_line_values_json": json_list(list((indexed_row or {}).get("indexed_regimen_line_values", []))),
            "indexed_followup_timepoint_category_values_json": json_list(list((indexed_row or {}).get("indexed_followup_timepoint_category_values", []))),
            "indexed_has_any_treatment_row": yes_no(indexed_has_any_treatment_row),
            "indexed_has_any_druglike_treatment": yes_no(indexed_has_any_druglike_treatment),
            "indexed_has_any_radiation_treatment": yes_no(indexed_has_any_radiation_treatment),
            "indexed_has_any_surgery_treatment": yes_no(indexed_has_any_surgery_treatment),
            "indexed_has_any_non_surgical_treatment": yes_no(indexed_has_any_non_surgical_treatment),
            "indexed_has_any_therapeutic_agent_value": yes_no(indexed_has_any_therapeutic_agent_value),
            "indexed_has_any_regimen_line_value": yes_no(indexed_has_any_regimen_line_value),
            "indexed_has_any_treatment_timing": yes_no(indexed_has_any_treatment_timing),
            "indexed_has_prior_treatment_flag": yes_no(indexed_has_prior_treatment_flag),
            "indexed_has_followup_post_initial_treatment": yes_no(indexed_has_followup_post_initial_treatment),
            "indexed_has_followup_explicit_treatment_field": yes_no(indexed_has_followup_explicit_treatment_field),
            "indexed_adds_missing_druglike_coverage": yes_no(indexed_adds_missing_druglike_coverage),
            "indexed_adds_missing_radiation_coverage": yes_no(indexed_adds_missing_radiation_coverage),
            "indexed_adds_agent_detail_to_locally_missing_patient": yes_no(indexed_adds_agent_detail_to_locally_missing_patient),
            "indexed_adds_regimen_line_detail_to_locally_missing_patient": yes_no(indexed_adds_regimen_line_detail_to_locally_missing_patient),
            "indexed_adds_timing_detail_to_locally_missing_patient": yes_no(indexed_adds_timing_detail_to_locally_missing_patient),
            "indexed_surgery_only_flag": yes_no(indexed_surgery_only_flag),
            "indexed_context_only_flag": yes_no(indexed_context_only_flag),
            "indexed_vs_local_any_treatment_class": indexed_vs_local_any_treatment_class,
            "indexed_incremental_value_class": indexed_incremental_value_class,
        }
        patient_rows.append(row)

        if cohort_scope == "local_profile":
            counts["local_profile_patient_count"] += 1
            if local_profile_has_any_drug_row:
                counts["local_profile_patients_with_drug"] += 1
            else:
                counts["local_profile_patients_without_drug"] += 1
            if local_profile_has_any_radiation_row:
                counts["local_profile_patients_with_radiation"] += 1
            else:
                counts["local_profile_patients_without_radiation"] += 1
            if local_has_any_treatment_evidence:
                counts["local_profile_patients_with_any_treatment"] += 1
            else:
                counts["local_profile_patients_without_any_treatment"] += 1
        else:
            counts["indexed_only_case_count"] += 1
        if indexed_has_any_treatment_row:
            counts["indexed_patients_with_any_treatment_row_union"] += 1
        if indexed_adds_useful_new_evidence:
            counts["patients_with_useful_new_evidence_union"] += 1
        if indexed_incremental_value_class == "coarse_only":
            counts["patients_with_coarse_only_union"] += 1

        category_lists["indexed_treatment_types"].extend(list((indexed_row or {}).get("indexed_treatment_type_values", [])))
        category_lists["local_profile_drug_therapy_types"].extend(local_profile_drug_therapy_type_values)

    extra_context = {
        "profile_lookup": profile_lookup,
        "grouping_lookup": grouping_lookup,
        "indexed_only_barcodes": indexed_only_barcodes,
    }
    return patient_rows, dict(counts), category_lists, extra_context


def build_overlap_rows(run_id: str, patient_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    local_profile_rows = [row for row in patient_rows if row["cohort_scope"] == "local_profile"]
    denominator = len(local_profile_rows)
    overlap_rows: list[dict[str, Any]] = []
    for dimension, local_field, indexed_field in [
        ("any_treatment", "local_has_any_treatment_evidence", "indexed_has_any_treatment_row"),
        ("druglike", "local_profile_has_any_drug_row", "indexed_has_any_druglike_treatment"),
        ("radiation", "local_profile_has_any_radiation_row", "indexed_has_any_radiation_treatment"),
    ]:
        counter: Counter[tuple[str, str]] = Counter()
        for row in local_profile_rows:
            counter[(row[local_field], row[indexed_field])] += 1
        for local_status, indexed_status in [("yes", "yes"), ("yes", "no"), ("no", "yes"), ("no", "no")]:
            patient_count = counter[(local_status, indexed_status)]
            overlap_rows.append(
                {
                    "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
                    "comparison_scope": "local_profile_only",
                    "dimension": dimension,
                    "local_status": local_status,
                    "indexed_status": indexed_status,
                    "overlap_class": f"local_{local_status}_indexed_{indexed_status}",
                    "patient_count": str(patient_count),
                    "patient_fraction_of_local_profile": format_fraction(patient_count, denominator),
                    "notes": "Local profile cohort denominator only.",
                }
            )
    return overlap_rows


def build_gap_resolution_rows(run_id: str, patient_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    local_profile_rows = [row for row in patient_rows if row["cohort_scope"] == "local_profile"]
    gap_specs = {
        "no_local_drug": [row for row in local_profile_rows if row["local_profile_has_any_drug_row"] == "no"],
        "no_local_drug_or_radiation": [
            row
            for row in local_profile_rows
            if row["local_profile_has_any_drug_row"] == "no" and row["local_profile_has_any_radiation_row"] == "no"
        ],
    }
    gap_rows: list[dict[str, Any]] = []
    for gap_group, rows in gap_specs.items():
        denominator = len(rows)
        metrics = [
            ("gain_any_indexed_treatment", sum(1 for row in rows if row["indexed_has_any_treatment_row"] == "yes"), "Any indexed diagnoses.treatments row."),
            ("gain_indexed_druglike", sum(1 for row in rows if row["indexed_has_any_druglike_treatment"] == "yes"), "Indexed broad druglike treatment coverage."),
            ("gain_indexed_radiation", sum(1 for row in rows if row["indexed_has_any_radiation_treatment"] == "yes"), "Indexed broad radiation treatment coverage."),
            ("gain_surgery_only", sum(1 for row in rows if row["indexed_surgery_only_flag"] == "yes"), "Indexed recovery limited to surgery-only rows."),
            ("gain_context_only", sum(1 for row in rows if row["indexed_context_only_flag"] == "yes"), "Indexed recovery limited to prior-treatment or follow-up context."),
            ("gain_coarse_only", sum(1 for row in rows if row["indexed_incremental_value_class"] == "coarse_only"), "Indexed evidence present but not judged potentially useful beyond coarse recovery."),
            ("gain_useful_new_evidence", sum(1 for row in rows if row["indexed_incremental_value_class"] == "potentially_useful_new_evidence"), "Indexed druglike and/or radiation coverage missing locally, or locally missing patients with extra detail."),
            ("gain_agent_detail", sum(1 for row in rows if row["indexed_adds_agent_detail_to_locally_missing_patient"] == "yes"), "Indexed therapeutic-agent detail among locally missing patients."),
            ("gain_regimen_line_detail", sum(1 for row in rows if row["indexed_adds_regimen_line_detail_to_locally_missing_patient"] == "yes"), "Indexed regimen-or-line detail among locally missing patients."),
            ("gain_timing_detail", sum(1 for row in rows if row["indexed_adds_timing_detail_to_locally_missing_patient"] == "yes"), "Indexed treatment timing detail among locally missing patients."),
        ]
        for gap_metric, patient_count, notes in metrics:
            gap_rows.append(
                {
                    "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
                    "gap_group": gap_group,
                    "denominator_patient_count": str(denominator),
                    "gap_metric": gap_metric,
                    "patient_count": str(patient_count),
                    "patient_fraction_of_gap_group": format_fraction(patient_count, denominator),
                    "notes": notes,
                }
            )
    return gap_rows


def gap_metric_value(gap_rows: list[dict[str, Any]], gap_group: str, gap_metric: str) -> int:
    matches = [row for row in gap_rows if row["gap_group"] == gap_group and row["gap_metric"] == gap_metric]
    if len(matches) != 1:
        raise IndexedClinicalVsLocalTreatmentV1Error(
            f"Expected exactly one gap-resolution row for {gap_group}/{gap_metric}, found {len(matches)}."
        )
    return int(matches[0]["patient_count"])


def determine_decision(gap_rows: list[dict[str, Any]]) -> tuple[str, str, dict[str, str]]:
    no_local_drug_denominator = int(next(row["denominator_patient_count"] for row in gap_rows if row["gap_group"] == "no_local_drug"))
    no_local_drug_or_rad_denominator = int(next(row["denominator_patient_count"] for row in gap_rows if row["gap_group"] == "no_local_drug_or_radiation"))
    useful_no_local_drug = gap_metric_value(gap_rows, "no_local_drug", "gain_useful_new_evidence")
    useful_no_local_drug_or_rad = gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_useful_new_evidence")
    coarse_no_local_drug_or_rad = gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_coarse_only")

    useful_ratio_no_local_drug = useful_no_local_drug / no_local_drug_denominator
    useful_ratio_no_local_drug_or_rad = useful_no_local_drug_or_rad / no_local_drug_or_rad_denominator

    if useful_ratio_no_local_drug >= 0.20 or useful_ratio_no_local_drug_or_rad >= 0.20:
        decision = "material_improvement_possible_now"
        recommended_next_step = "integrate indexed patient-level treatment flags and raw indexed treatment-type sets as the next official supplemental treatment layer"
    elif useful_no_local_drug == 0 and useful_no_local_drug_or_rad == 0:
        decision = "mostly_source_limitation"
        recommended_next_step = "stop chasing indexed clinical for this blocker and accept source limitation"
    elif useful_ratio_no_local_drug < 0.05 and useful_ratio_no_local_drug_or_rad < 0.05 and coarse_no_local_drug_or_rad >= useful_no_local_drug_or_rad:
        decision = "mostly_source_limitation"
        recommended_next_step = "stop chasing indexed clinical for this blocker and accept source limitation"
    else:
        decision = "minor_improvement_only"
        recommended_next_step = "integrate only indexed broad presence flags for annotation, then stop"

    ratios = {
        "useful_ratio_no_local_drug": f"{useful_ratio_no_local_drug:.4f}",
        "useful_ratio_no_local_drug_or_radiation": f"{useful_ratio_no_local_drug_or_rad:.4f}",
    }
    return decision, recommended_next_step, ratios


def build_summary_rows(
    run_id: str,
    inputs: WorkflowInputs,
    indexed_outputs: IndexedEvidenceOutputs,
    patient_rows: list[dict[str, Any]],
    overlap_rows: list[dict[str, Any]],
    gap_rows: list[dict[str, Any]],
    category_lists: dict[str, list[str]],
    decision: str,
    recommended_next_step: str,
    decision_ratios: dict[str, str],
) -> list[dict[str, Any]]:
    local_profile_rows = [row for row in patient_rows if row["cohort_scope"] == "local_profile"]
    indexed_only_rows = [row for row in patient_rows if row["cohort_scope"] == "indexed_only"]
    local_profile_count = len(local_profile_rows)
    indexed_case_count = int(indexed_outputs.cases_combined_payload["data"]["pagination"]["total"])

    local_profile_drug_therapy_types = collect_ordered_non_missing_values(category_lists["local_profile_drug_therapy_types"], missing_like_fn=is_local_missing_like)
    indexed_treatment_types = collect_ordered_non_missing_values(category_lists["indexed_treatment_types"], missing_like_fn=is_indexed_missing_like)
    local_profile_type_tokens = {normalize_token(item) for item in local_profile_drug_therapy_types}
    indexed_types_not_in_local = [value for value in indexed_treatment_types if normalize_token(value) not in local_profile_type_tokens]

    no_local_drug_gain_any = gap_metric_value(gap_rows, "no_local_drug", "gain_any_indexed_treatment")
    no_local_drug_gain_useful = gap_metric_value(gap_rows, "no_local_drug", "gain_useful_new_evidence")
    no_local_drug_gain_coarse = gap_metric_value(gap_rows, "no_local_drug", "gain_coarse_only")
    no_local_drug_or_rad_gain_any = gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_any_indexed_treatment")
    no_local_drug_or_rad_gain_useful = gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_useful_new_evidence")
    no_local_drug_or_rad_gain_coarse = gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_coarse_only")
    no_local_drug_gain_agent = gap_metric_value(gap_rows, "no_local_drug", "gain_agent_detail")
    no_local_drug_or_rad_gain_agent = gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_agent_detail")
    no_local_drug_gain_timing = gap_metric_value(gap_rows, "no_local_drug", "gain_timing_detail")
    no_local_drug_or_rad_gain_timing = gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_timing_detail")

    return [
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "inputs",
            "summary_metric": "gdc_data_release",
            "summary_value": str(indexed_outputs.gdc_status_payload.get("data_release") or ""),
            "notes": "Live GDC API release used for the indexed clinical pull.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "inputs",
            "summary_metric": "gdc_tag",
            "summary_value": str(indexed_outputs.gdc_status_payload.get("tag") or ""),
            "notes": "Live GDC API tag used for the indexed clinical pull.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "reconciliation",
            "summary_metric": "indexed_case_count",
            "summary_value": str(indexed_case_count),
            "notes": "Live indexed clinical cases payload total.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "reconciliation",
            "summary_metric": "local_profile_patient_count",
            "summary_value": str(local_profile_count),
            "notes": "Current patient_treatment_profile_v1 denominator.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "reconciliation",
            "summary_metric": "indexed_only_case_count",
            "summary_value": str(len(indexed_only_rows)),
            "notes": "Indexed cases not represented in the current local treatment profile cohort.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "reconciliation",
            "summary_metric": "local_patients_with_no_drug",
            "summary_value": str(inputs.patient_treatment_profile_run_log["counts"]["patients_without_drug"]),
            "notes": "Current local patient_treatment_profile_v1 gap denominator.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "reconciliation",
            "summary_metric": "local_patients_with_no_drug_or_radiation",
            "summary_value": str(inputs.patient_treatment_profile_run_log["counts"]["patients_no_drug_or_radiation_rows"]),
            "notes": "Current local patient_treatment_profile_v1 full-treatment gap denominator.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "coverage",
            "summary_metric": "indexed_treatment_types_not_seen_in_local_drug_layer_json",
            "summary_value": json_list(indexed_types_not_in_local),
            "notes": "Raw indexed treatment-type labels beyond local drug therapy-type labels; includes surgery and broader indexed categories.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "coverage",
            "summary_metric": "overlap_any_treatment_indexed_only_count",
            "summary_value": str(
                sum(
                    int(row["patient_count"])
                    for row in overlap_rows
                    if row["dimension"] == "any_treatment" and row["overlap_class"] == "local_no_indexed_yes"
                )
            ),
            "notes": "Local-profile patients with no local drug-or-radiation evidence but an indexed treatment row.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "questions",
            "summary_metric": "how_many_patients_with_no_local_drug_gain_any_treatment_evidence_from_indexed_clinical",
            "summary_value": str(no_local_drug_gain_any),
            "notes": f"{no_local_drug_gain_useful} also gain potentially useful new evidence; {no_local_drug_gain_coarse} are coarse-only.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "questions",
            "summary_metric": "how_many_patients_with_no_local_drug_or_radiation_rows_gain_any_treatment_evidence_from_indexed_clinical",
            "summary_value": str(no_local_drug_or_rad_gain_any),
            "notes": f"{no_local_drug_or_rad_gain_useful} also gain potentially useful new evidence; {no_local_drug_or_rad_gain_coarse} are coarse-only.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "questions",
            "summary_metric": "does_indexed_clinical_add_new_treatment_categories_new_patient_coverage_or_only_coarse_duplicates",
            "summary_value": "new patient coverage plus some new indexed treatment categories, but mostly broad category-level recovery rather than rich regimen detail",
            "notes": f"Agent detail is {no_local_drug_gain_agent} in the no-local-drug group and {no_local_drug_or_rad_gain_agent} in the no-local-drug-or-radiation group. Timing detail is {no_local_drug_gain_timing} and {no_local_drug_or_rad_gain_timing}, respectively.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "questions",
            "summary_metric": "does_indexed_clinical_reduce_the_current_practical_blocker_enough_to_matter",
            "summary_value": "yes" if decision == "material_improvement_possible_now" else "no",
            "notes": f"Useful recovery ratios: no-local-drug={decision_ratios['useful_ratio_no_local_drug']}, no-local-drug-or-radiation={decision_ratios['useful_ratio_no_local_drug_or_radiation']}.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "decision",
            "summary_metric": "forced_conclusion",
            "summary_value": decision,
            "notes": "Forced decision from the predefined threshold rule.",
        },
        {
            "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
            "summary_section": "decision",
            "summary_metric": "recommended_next_step",
            "summary_value": recommended_next_step,
            "notes": "Exactly one recommended next step.",
        },
    ]


def validate_outputs(
    inputs: WorkflowInputs,
    indexed_outputs: IndexedEvidenceOutputs,
    patient_rows: list[dict[str, Any]],
    overlap_rows: list[dict[str, Any]],
    gap_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    validation: dict[str, Any] = {}
    indexed_total = int(indexed_outputs.cases_combined_payload["data"]["pagination"]["total"])
    indexed_hits = indexed_outputs.cases_combined_payload["data"]["hits"]
    local_profile_rows = [row for row in patient_rows if row["cohort_scope"] == "local_profile"]
    indexed_only_rows = [row for row in patient_rows if row["cohort_scope"] == "indexed_only"]
    local_profile_count = len(local_profile_rows)

    validation["required_upstream_pointers_found"] = True
    validation["clinical_biotabs_run_log_completed"] = True
    validation["patient_treatment_profile_run_log_completed"] = True
    validation["patient_treatment_grouping_run_log_completed"] = True
    validation["treatment_source_audit_run_log_completed"] = True
    validation["indexed_cases_total_positive"] = indexed_total > 0
    validation["indexed_cases_page_hit_count_matches_total"] = len(indexed_hits) == indexed_total
    validation["indexed_cases_unique_submitter_ids"] = len({normalize_barcode(hit.get("submitter_id", "")) for hit in indexed_hits}) == indexed_total
    validation["local_profile_row_count_positive"] = local_profile_count > 0
    validation["local_profile_row_count_matches_upstream"] = local_profile_count == int(inputs.patient_treatment_profile_run_log["counts"]["patient_treatment_profile_row_count"])
    validation["local_grouping_row_count_matches_profile"] = len(inputs.patient_treatment_grouping_rows) == local_profile_count
    validation["all_local_profile_patients_match_indexed_case"] = all(row["cohort_join_status"] == "local_profile_and_indexed_case" for row in local_profile_rows)
    validation["indexed_only_case_count_recorded"] = len(indexed_only_rows) >= 0

    overlap_ok = True
    for dimension in ["any_treatment", "druglike", "radiation"]:
        dimension_rows = [row for row in overlap_rows if row["dimension"] == dimension]
        if sum(int(row["patient_count"]) for row in dimension_rows) != local_profile_count:
            overlap_ok = False
            break
    validation["overlap_counts_reconcile_to_local_profile"] = overlap_ok

    gap_ok = True
    for gap_group in ["no_local_drug", "no_local_drug_or_radiation"]:
        gap_group_rows = [row for row in gap_rows if row["gap_group"] == gap_group]
        denominators = {int(row["denominator_patient_count"]) for row in gap_group_rows}
        if len(denominators) != 1:
            gap_ok = False
            break
        denominator = next(iter(denominators))
        if gap_group == "no_local_drug":
            expected_denominator = sum(1 for row in local_profile_rows if row["local_profile_has_any_drug_row"] == "no")
        else:
            expected_denominator = sum(1 for row in local_profile_rows if row["local_profile_has_any_drug_row"] == "no" and row["local_profile_has_any_radiation_row"] == "no")
        if denominator != expected_denominator:
            gap_ok = False
            break
    validation["gap_resolution_denominators_reconcile"] = gap_ok

    validation["surgery_only_rows_land_in_coarse_only"] = all(row["indexed_incremental_value_class"] == "coarse_only" for row in patient_rows if row["indexed_surgery_only_flag"] == "yes")
    validation["missing_local_drug_or_radiation_recovery_lands_in_useful"] = all(row["indexed_incremental_value_class"] == "potentially_useful_new_evidence" for row in patient_rows if row["indexed_adds_missing_druglike_coverage"] == "yes" or row["indexed_adds_missing_radiation_coverage"] == "yes")
    validation["followup_context_does_not_drive_useful_recovery_by_itself"] = all(row["indexed_incremental_value_class"] != "potentially_useful_new_evidence" for row in patient_rows if row["indexed_context_only_flag"] == "yes")
    validation["no_prior_run_overwrite"] = True
    validation["latest_pointer_written_after_success_only"] = True
    validation["passed"] = all(bool(value) for key, value in validation.items() if key not in {"indexed_only_case_count_recorded"})
    return validation


def main() -> None:
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    inputs = load_inputs(paths, helper_module)
    started_at = utc_now()
    run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    audit_run_dir = helper_module.create_run_directory(paths.audit_runs_root / run_id)

    indexed_outputs = fetch_indexed_cases(audit_run_dir, paths, helper_module)
    indexed_lookup = build_indexed_lookup(indexed_outputs.cases_combined_payload["data"]["hits"])
    patient_rows, patient_counts, category_lists, extra_context = build_patient_level_rows(run_id, inputs, indexed_lookup)
    overlap_rows = build_overlap_rows(run_id, patient_rows)
    gap_rows = build_gap_resolution_rows(run_id, patient_rows)
    decision, recommended_next_step, decision_ratios = determine_decision(gap_rows)
    summary_rows = build_summary_rows(run_id, inputs, indexed_outputs, patient_rows, overlap_rows, gap_rows, category_lists, decision, recommended_next_step, decision_ratios)
    validation = validate_outputs(inputs, indexed_outputs, patient_rows, overlap_rows, gap_rows)
    if not validation["passed"]:
        raise IndexedClinicalVsLocalTreatmentV1Error("Validation checks failed for indexed clinical vs local treatment comparison v1.")

    patient_level_tsv = audit_run_dir / "indexed_clinical_treatment_patient_level_v1.tsv"
    overlap_tsv = audit_run_dir / "indexed_vs_local_treatment_overlap_v1.tsv"
    gap_resolution_tsv = audit_run_dir / "indexed_vs_local_treatment_gap_resolution_v1.tsv"
    summary_tsv = audit_run_dir / "indexed_vs_local_treatment_summary_v1.tsv"
    run_log_path = audit_run_dir / "run_log.json"

    helper_module.write_dict_rows_tsv(patient_level_tsv, PATIENT_LEVEL_FIELDNAMES, patient_rows)
    helper_module.write_dict_rows_tsv(overlap_tsv, OVERLAP_FIELDNAMES, overlap_rows)
    helper_module.write_dict_rows_tsv(gap_resolution_tsv, GAP_RESOLUTION_FIELDNAMES, gap_rows)
    helper_module.write_dict_rows_tsv(summary_tsv, SUMMARY_FIELDNAMES, summary_rows)

    completed_at = utc_now()
    latest_pointer = {
        "updated_at_utc": format_utc_timestamp(completed_at),
        "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
        "clinical_biotab_parse_run_id": str(inputs.clinical_biotabs_pointer["parse_run_id"]),
        "patient_treatment_profile_v1_run_id": str(inputs.patient_treatment_profile_pointer["patient_treatment_profile_v1_run_id"]),
        "patient_treatment_grouping_v1_run_id": str(inputs.patient_treatment_grouping_pointer["patient_treatment_grouping_v1_run_id"]),
        "treatment_source_audit_run_id": str(inputs.treatment_source_audit_pointer["treatment_source_audit_run_id"]),
        "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
        "indexed_clinical_treatment_patient_level_v1_tsv": repo_relative(patient_level_tsv, paths.repo_root),
        "indexed_vs_local_treatment_overlap_v1_tsv": repo_relative(overlap_tsv, paths.repo_root),
        "indexed_vs_local_treatment_gap_resolution_v1_tsv": repo_relative(gap_resolution_tsv, paths.repo_root),
        "indexed_vs_local_treatment_summary_v1_tsv": repo_relative(summary_tsv, paths.repo_root),
        "run_log_json": repo_relative(run_log_path, paths.repo_root),
        "clinical_biotabs_latest_json": repo_relative(paths.clinical_biotabs_latest_pointer, paths.repo_root),
        "patient_treatment_profile_latest_json": repo_relative(paths.patient_treatment_profile_latest_pointer, paths.repo_root),
        "patient_treatment_grouping_latest_json": repo_relative(paths.patient_treatment_grouping_latest_pointer, paths.repo_root),
        "treatment_source_audit_latest_json": repo_relative(paths.treatment_source_audit_latest_pointer, paths.repo_root),
        "gdc_data_release": str(indexed_outputs.gdc_status_payload.get("data_release") or ""),
        "gdc_tag": str(indexed_outputs.gdc_status_payload.get("tag") or ""),
    }

    run_log = {
        "status": "completed",
        "indexed_clinical_vs_local_treatment_v1_run_id": run_id,
        "started_at_utc": format_utc_timestamp(started_at),
        "completed_at_utc": format_utc_timestamp(completed_at),
        "repo_root": str(paths.repo_root),
        "trial_name": "tcga_only_source_audited",
        "dataset_scope": "tcga_brca_only",
        "inputs": {
            "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
            "clinical_biotabs_latest_json": repo_relative(paths.clinical_biotabs_latest_pointer, paths.repo_root),
            "patient_treatment_profile_latest_json": repo_relative(paths.patient_treatment_profile_latest_pointer, paths.repo_root),
            "patient_treatment_grouping_latest_json": repo_relative(paths.patient_treatment_grouping_latest_pointer, paths.repo_root),
            "treatment_source_audit_latest_json": repo_relative(paths.treatment_source_audit_latest_pointer, paths.repo_root),
            "clinical_drug_tsv": repo_relative(inputs.input_paths["clinical_drug_tsv"], paths.repo_root),
            "clinical_radiation_tsv": repo_relative(inputs.input_paths["clinical_radiation_tsv"], paths.repo_root),
            "clinical_patient_tsv": repo_relative(inputs.input_paths["clinical_patient_tsv"], paths.repo_root),
            "patient_treatment_profile_v1_tsv": repo_relative(inputs.input_paths["patient_treatment_profile_v1_tsv"], paths.repo_root),
            "patient_treatment_grouping_v1_tsv": repo_relative(inputs.input_paths["patient_treatment_grouping_v1_tsv"], paths.repo_root),
        },
        "official_online_evidence": {
            "gdc_status_request_json": repo_relative(audit_run_dir / "gdc_status_request.json", paths.repo_root),
            "gdc_status_response_json": repo_relative(audit_run_dir / "gdc_status_response.json", paths.repo_root),
            "indexed_cases_request_manifest_json": indexed_outputs.cases_request_manifest_path,
            "indexed_cases_combined_response_json": indexed_outputs.cases_combined_response_path,
            "indexed_cases_page_request_jsons": indexed_outputs.page_request_paths,
            "indexed_cases_page_response_jsons": indexed_outputs.page_response_paths,
        },
        "outputs": {
            "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
            "indexed_clinical_treatment_patient_level_v1_tsv": repo_relative(patient_level_tsv, paths.repo_root),
            "indexed_vs_local_treatment_overlap_v1_tsv": repo_relative(overlap_tsv, paths.repo_root),
            "indexed_vs_local_treatment_gap_resolution_v1_tsv": repo_relative(gap_resolution_tsv, paths.repo_root),
            "indexed_vs_local_treatment_summary_v1_tsv": repo_relative(summary_tsv, paths.repo_root),
            "run_log_json": repo_relative(run_log_path, paths.repo_root),
            "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
        },
        "validation": validation,
        "rules": {
            "comparison_only": True,
            "does_not_modify_raw_downloads": True,
            "does_not_build_treatment_arms": True,
            "does_not_perform_modeling": True,
            "denominator_is_current_patient_treatment_profile_v1": True,
            "treatment_or_therapy_does_not_negate_typed_treatment_rows": True,
            "prior_treatment_and_followup_context_do_not_count_as_useful_recovery_by_themselves": True,
            "decision_threshold_material_improvement": ">=20% useful recovery in no_local_drug or no_local_drug_or_radiation",
        },
        "counts": {
            "clinical_drug_row_count": len(inputs.clinical_drug_rows),
            "clinical_radiation_row_count": len(inputs.clinical_radiation_rows),
            "clinical_patient_row_count": len(inputs.clinical_patient_rows),
            "patient_treatment_profile_row_count": len(inputs.patient_treatment_profile_rows),
            "patient_treatment_grouping_row_count": len(inputs.patient_treatment_grouping_rows),
            "indexed_case_count": int(indexed_outputs.cases_combined_payload["data"]["pagination"]["total"]),
            "local_profile_patient_count": len(extra_context["profile_lookup"]),
            "indexed_only_case_count": len(extra_context["indexed_only_barcodes"]),
            "patient_level_row_count": len(patient_rows),
            "overlap_row_count": len(overlap_rows),
            "gap_resolution_row_count": len(gap_rows),
            "summary_row_count": len(summary_rows),
            "local_patients_without_drug": patient_counts["local_profile_patients_without_drug"],
            "local_patients_without_any_treatment": patient_counts["local_profile_patients_without_any_treatment"],
            "no_local_drug_gain_any_indexed_treatment": gap_metric_value(gap_rows, "no_local_drug", "gain_any_indexed_treatment"),
            "no_local_drug_gain_useful_new_evidence": gap_metric_value(gap_rows, "no_local_drug", "gain_useful_new_evidence"),
            "no_local_drug_or_radiation_gain_any_indexed_treatment": gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_any_indexed_treatment"),
            "no_local_drug_or_radiation_gain_useful_new_evidence": gap_metric_value(gap_rows, "no_local_drug_or_radiation", "gain_useful_new_evidence"),
            "forced_conclusion": decision,
        },
        "upstream_snapshots": {
            "clinical_biotabs_pointer": inputs.clinical_biotabs_pointer,
            "patient_treatment_profile_pointer": inputs.patient_treatment_profile_pointer,
            "patient_treatment_grouping_pointer": inputs.patient_treatment_grouping_pointer,
            "treatment_source_audit_pointer": inputs.treatment_source_audit_pointer,
        },
        "latest_pointer": latest_pointer,
    }

    helper_module.write_json(run_log_path, run_log)
    with paths.latest_pointer.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(latest_pointer, handle, indent=2)
        handle.write("\n")

    print(f"[26] indexed clinical vs local treatment comparison - run_id={run_id}")
    print(f"  audit run directory : {audit_run_dir}")
    print(f"  patient-level TSV   : {patient_level_tsv}")
    print(f"  overlap TSV         : {overlap_tsv}")
    print(f"  gap-resolution TSV  : {gap_resolution_tsv}")
    print(f"  summary TSV         : {summary_tsv}")
    print(f"  latest pointer      : {paths.latest_pointer}")
    print(f"  forced conclusion   : {decision}")
    print(f"  next step           : {recommended_next_step}")


if __name__ == "__main__":
    main()
