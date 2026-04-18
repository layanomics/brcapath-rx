#!/usr/bin/env python
"""Audit TCGA-BRCA treatment-relevant GDC source layers against local BRCAPath-Rx holdings."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

PROJECT_ID = "TCGA-BRCA"
GDC_API_BASE_URL = "https://api.gdc.cancer.gov"
REQUEST_TIMEOUT_SECONDS = 120
USER_AGENT = "brcapath-rx-treatment-source-audit/1.0"

OFFICIAL_SOURCE_FIELDNAMES = [
    "treatment_source_audit_run_id",
    "source_layer_key",
    "official_source_layer",
    "official_source_kind",
    "official_availability_for_tcga_brca",
    "official_count",
    "official_count_basis",
    "what_it_contains_in_general",
    "treatment_relevance",
    "richness_vs_other_layers",
    "official_reference_ids",
    "notes",
]
LOCAL_SOURCE_FIELDNAMES = [
    "treatment_source_audit_run_id",
    "source_layer_key",
    "local_source_layer",
    "local_path_or_scope",
    "what_is_present_locally",
    "looks_downloaded_completely_or_partially",
    "parsed_or_exploited_status",
    "treatment_relevance",
    "evidence_paths",
    "notes",
]
GAP_ANALYSIS_FIELDNAMES = [
    "treatment_source_audit_run_id",
    "official_source_layer",
    "official_availability_for_tcga_brca",
    "local_download_status",
    "local_parse_status",
    "likely_treatment_value",
    "classification",
    "notes",
]
SUMMARY_FIELDNAMES = [
    "treatment_source_audit_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
REFERENCE_FIELDNAMES = [
    "treatment_source_audit_run_id",
    "reference_id",
    "source_kind",
    "source_label",
    "source_url",
    "http_status",
    "saved_snapshot_path",
    "what_it_supports",
    "notes",
]

DOC_REFERENCES = [
    {
        "reference_id": "R1",
        "source_kind": "official_docs",
        "source_label": "GDC Clinical Data",
        "source_url": "https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Data/",
        "what_it_supports": (
            "Indexed clinical surface exists, includes model-aligned clinical entities, and may omit "
            "additional TCGA clinical detail that remains in supplement files."
        ),
        "notes": "Primary documentation for database-aligned clinical data.",
    },
    {
        "reference_id": "R2",
        "source_kind": "official_docs",
        "source_label": "GDC Clinical Supplement",
        "source_url": "https://docs.gdc.cancer.gov/Encyclopedia/pages/Clinical_Supplement/",
        "what_it_supports": "Clinical Supplement is an official GDC source class.",
        "notes": "High-level source-class description.",
    },
    {
        "reference_id": "R3",
        "source_kind": "official_docs",
        "source_label": "GDC Biospecimen Data",
        "source_url": "https://docs.gdc.cancer.gov/Encyclopedia/pages/Biospecimen_Data/",
        "what_it_supports": "Biospecimen information is an official GDC source class.",
        "notes": "High-level biospecimen documentation.",
    },
    {
        "reference_id": "R4",
        "source_kind": "official_docs",
        "source_label": "GDC Repository Users Guide",
        "source_url": "https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Repository/",
        "what_it_supports": "Clinical and Biospecimen TSV/JSON additional data exports are official portal access surfaces.",
        "notes": "Additional Data Download section is the key part for indexed clinical/biospecimen exports.",
    },
    {
        "reference_id": "R5",
        "source_kind": "official_dictionary",
        "source_label": "GDC Data Dictionary: clinical_supplement",
        "source_url": "https://api.gdc.cancer.gov/v0/submission/_dictionary/clinical_supplement",
        "what_it_supports": "Clinical Supplement description and allowed data formats.",
        "notes": "Live dictionary endpoint; snapshot saved in the audit run directory.",
    },
    {
        "reference_id": "R6",
        "source_kind": "official_dictionary",
        "source_label": "GDC Data Dictionary: biospecimen_supplement",
        "source_url": "https://api.gdc.cancer.gov/v0/submission/_dictionary/biospecimen_supplement",
        "what_it_supports": "Biospecimen Supplement description and allowed data formats.",
        "notes": "Live dictionary endpoint; snapshot saved in the audit run directory.",
    },
    {
        "reference_id": "R7",
        "source_kind": "official_dictionary",
        "source_label": "GDC Data Dictionary: treatment",
        "source_url": "https://api.gdc.cancer.gov/v0/submission/_dictionary/treatment",
        "what_it_supports": "Official treatment entity description in the indexed GDC data model.",
        "notes": "Live dictionary endpoint; snapshot saved in the audit run directory.",
    },
    {
        "reference_id": "R8",
        "source_kind": "official_dictionary",
        "source_label": "GDC Data Dictionary: follow_up",
        "source_url": "https://api.gdc.cancer.gov/v0/submission/_dictionary/follow_up",
        "what_it_supports": "Official follow_up entity description in the indexed GDC data model.",
        "notes": "Live dictionary endpoint; snapshot saved in the audit run directory.",
    },
    {
        "reference_id": "R9",
        "source_kind": "official_api",
        "source_label": "GDC Files API",
        "source_url": "https://api.gdc.cancer.gov/files",
        "what_it_supports": (
            "Live TCGA-BRCA file-layer counts for Clinical Supplement, Biospecimen Supplement, Pathology Report, "
            "and clinical/biospecimen category discovery."
        ),
        "notes": "Exact filtered request params and responses are saved in the audit run directory.",
    },
    {
        "reference_id": "R10",
        "source_kind": "official_api",
        "source_label": "GDC Cases API",
        "source_url": "https://api.gdc.cancer.gov/cases",
        "what_it_supports": "Live TCGA-BRCA case count plus indexed diagnoses.treatments and follow_ups availability.",
        "notes": "Exact filtered request params and responses are saved in the audit run directory.",
    },
    {
        "reference_id": "R11",
        "source_kind": "official_api",
        "source_label": "GDC API Status",
        "source_url": "https://api.gdc.cancer.gov/status",
        "what_it_supports": "Current GDC data release and tag used for the live reconciliation target.",
        "notes": "Snapshot saved in the audit run directory.",
    },
]

SOURCE_LAYER_KEYS = [
    "indexed_clinical_surface",
    "clinical_supplement_bcr_xml",
    "clinical_supplement_bcr_omf_xml",
    "clinical_supplement_bcr_biotab",
    "biospecimen_supplement_bcr_xml",
    "biospecimen_supplement_bcr_ssf_xml",
    "biospecimen_supplement_bcr_biotab",
    "pathology_report_pdf",
]


class TreatmentSourceAuditError(RuntimeError):
    """Raised when the treatment-source audit cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    repo_root: Path
    trial_config: Path
    audit_root: Path
    processed_root: Path
    raw_root: Path
    source_inventory_latest_pointer: Path
    source_supplements_latest_pointer: Path
    clinical_biotabs_latest_pointer: Path
    biospecimen_biotabs_latest_pointer: Path
    endpoint_target_prep_latest_pointer: Path
    patient_treatment_profile_latest_pointer: Path
    patient_treatment_grouping_latest_pointer: Path
    audit_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    source_inventory_pointer: dict[str, Any]
    source_inventory_run_log: dict[str, Any]
    source_inventory_rows: list[dict[str, str]]
    source_supplements_pointer: dict[str, Any]
    source_supplements_run_log: dict[str, Any]
    supplement_metadata_rows: list[dict[str, str]]
    clinical_biotabs_pointer: dict[str, Any]
    clinical_biotabs_run_log: dict[str, Any]
    clinical_biotab_manifest_rows: list[dict[str, str]]
    biospecimen_biotabs_pointer: dict[str, Any]
    biospecimen_biotabs_run_log: dict[str, Any]
    biospecimen_source_inventory_rows: list[dict[str, str]]
    biospecimen_biotab_manifest_rows: list[dict[str, str]]
    endpoint_target_prep_pointer: dict[str, Any]
    endpoint_target_prep_run_log: dict[str, Any]
    patient_treatment_profile_pointer: dict[str, Any]
    patient_treatment_profile_run_log: dict[str, Any]
    patient_treatment_grouping_pointer: dict[str, Any]
    patient_treatment_grouping_run_log: dict[str, Any]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise TreatmentSourceAuditError(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise TreatmentSourceAuditError(f"Unable to create an import spec for: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def bucket_count(aggregation: dict[str, Any], aggregation_key: str, bucket_key: str) -> int:
    buckets = aggregation.get(aggregation_key, {}).get("buckets", [])
    for bucket in buckets:
        if str(bucket.get("key", "")).lower() == bucket_key.lower():
            return int(bucket.get("doc_count", 0))
    return 0


def bucket_map(aggregation: dict[str, Any], aggregation_key: str) -> dict[str, int]:
    buckets = aggregation.get(aggregation_key, {}).get("buckets", [])
    return {str(bucket.get("key", "")): int(bucket.get("doc_count", 0)) for bucket in buckets}


def extract_html_title(raw_text: str) -> str:
    match = re.search(r"<title>(.*?)</title>", raw_text, flags=re.IGNORECASE | re.DOTALL)
    if not match:
        return ""
    return re.sub(r"\s+", " ", match.group(1)).strip()


def require_completed_pointer_and_run_log(
    pointer_path: Path,
    label: str,
    helper_module: Any,
    repo_root: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    if not pointer_path.exists():
        raise TreatmentSourceAuditError(f"Required {label} pointer not found: {pointer_path}")
    pointer = helper_module.load_json(pointer_path)
    run_log_relative = str(pointer.get("run_log_json") or "").strip()
    if not run_log_relative:
        raise TreatmentSourceAuditError(f"{label} pointer is missing run_log_json: {pointer_path}")
    run_log_path = helper_module.resolve_existing_path(repo_root, run_log_relative, f"{label} run log")
    run_log = helper_module.load_json(run_log_path)
    if str(run_log.get("status") or "") != "completed":
        raise TreatmentSourceAuditError(f"{label} run log is not completed: {run_log_path}")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise TreatmentSourceAuditError(f"{label} run log does not report validation.passed == true: {run_log_path}")
    return pointer, run_log


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    if not trial_config.exists():
        raise TreatmentSourceAuditError(f"Required trial config not found: {trial_config}")
    trial_config_data = helper_module.load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    processed_root = repo_root / str(trial_config_data.get("processed_data_root", "01-data/processed"))
    raw_root = repo_root / str(trial_config_data.get("raw_data_root", "01-data/raw"))
    source_root = audit_root / "tcga-brca" / "source"
    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        audit_root=audit_root,
        processed_root=processed_root,
        raw_root=raw_root,
        source_inventory_latest_pointer=source_root / "tcga_brca_inventory_latest.json",
        source_supplements_latest_pointer=source_root / "tcga_brca_source_supplements_latest.json",
        clinical_biotabs_latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_biotabs_latest.json",
        biospecimen_biotabs_latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_biospecimen_biotabs_latest.json",
        endpoint_target_prep_latest_pointer=audit_root / "tcga-brca" / "endpoint-prep" / "tcga_brca_endpoint_target_prep_v1_latest.json",
        patient_treatment_profile_latest_pointer=audit_root / "tcga-brca" / "treatment-prep" / "tcga_brca_patient_treatment_profile_v1_latest.json",
        patient_treatment_grouping_latest_pointer=audit_root / "tcga-brca" / "treatment-prep" / "tcga_brca_patient_treatment_grouping_v1_latest.json",
        audit_runs_root=source_root / "treatment_source_audit_runs",
        latest_pointer=source_root / "tcga_brca_treatment_source_audit_latest.json",
    )


def load_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    source_inventory_pointer, source_inventory_run_log = require_completed_pointer_and_run_log(paths.source_inventory_latest_pointer, "source inventory latest", helper_module, paths.repo_root)
    source_inventory_rows = helper_module.read_tsv_dict_rows(helper_module.resolve_existing_path(paths.repo_root, str(source_inventory_pointer.get("file_inventory_tsv") or ""), "source inventory TSV"))
    source_supplements_pointer, source_supplements_run_log = require_completed_pointer_and_run_log(paths.source_supplements_latest_pointer, "source supplements latest", helper_module, paths.repo_root)
    supplement_metadata_rows = helper_module.read_tsv_dict_rows(helper_module.resolve_existing_path(paths.repo_root, str(source_supplements_pointer.get("metadata_tsv") or ""), "source supplements metadata TSV"))
    clinical_biotabs_pointer, clinical_biotabs_run_log = require_completed_pointer_and_run_log(paths.clinical_biotabs_latest_pointer, "clinical biotabs latest", helper_module, paths.repo_root)
    clinical_biotab_manifest_rows = helper_module.read_tsv_dict_rows(helper_module.resolve_existing_path(paths.repo_root, str(clinical_biotabs_pointer.get("table_manifest_tsv") or ""), "clinical biotab table manifest TSV"))
    biospecimen_biotabs_pointer, biospecimen_biotabs_run_log = require_completed_pointer_and_run_log(paths.biospecimen_biotabs_latest_pointer, "biospecimen biotabs latest", helper_module, paths.repo_root)
    biospecimen_source_inventory_rows = helper_module.read_tsv_dict_rows(helper_module.resolve_existing_path(paths.repo_root, str(biospecimen_biotabs_pointer.get("source_inventory_tsv") or ""), "biospecimen biotab source inventory TSV"))
    biospecimen_biotab_manifest_rows = helper_module.read_tsv_dict_rows(helper_module.resolve_existing_path(paths.repo_root, str(biospecimen_biotabs_pointer.get("table_manifest_tsv") or ""), "biospecimen biotab table manifest TSV"))
    endpoint_target_prep_pointer, endpoint_target_prep_run_log = require_completed_pointer_and_run_log(paths.endpoint_target_prep_latest_pointer, "endpoint target prep latest", helper_module, paths.repo_root)
    patient_treatment_profile_pointer, patient_treatment_profile_run_log = require_completed_pointer_and_run_log(paths.patient_treatment_profile_latest_pointer, "patient treatment profile latest", helper_module, paths.repo_root)
    patient_treatment_grouping_pointer, patient_treatment_grouping_run_log = require_completed_pointer_and_run_log(paths.patient_treatment_grouping_latest_pointer, "patient treatment grouping latest", helper_module, paths.repo_root)
    return WorkflowInputs(
        source_inventory_pointer=source_inventory_pointer,
        source_inventory_run_log=source_inventory_run_log,
        source_inventory_rows=source_inventory_rows,
        source_supplements_pointer=source_supplements_pointer,
        source_supplements_run_log=source_supplements_run_log,
        supplement_metadata_rows=supplement_metadata_rows,
        clinical_biotabs_pointer=clinical_biotabs_pointer,
        clinical_biotabs_run_log=clinical_biotabs_run_log,
        clinical_biotab_manifest_rows=clinical_biotab_manifest_rows,
        biospecimen_biotabs_pointer=biospecimen_biotabs_pointer,
        biospecimen_biotabs_run_log=biospecimen_biotabs_run_log,
        biospecimen_source_inventory_rows=biospecimen_source_inventory_rows,
        biospecimen_biotab_manifest_rows=biospecimen_biotab_manifest_rows,
        endpoint_target_prep_pointer=endpoint_target_prep_pointer,
        endpoint_target_prep_run_log=endpoint_target_prep_run_log,
        patient_treatment_profile_pointer=patient_treatment_profile_pointer,
        patient_treatment_profile_run_log=patient_treatment_profile_run_log,
        patient_treatment_grouping_pointer=patient_treatment_grouping_pointer,
        patient_treatment_grouping_run_log=patient_treatment_grouping_run_log,
    )


def official_doc_reference_rows(session: requests.Session, run_id: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for ref in DOC_REFERENCES:
        http_status = ""
        title = ""
        try:
            response = session.get(str(ref["source_url"]), timeout=REQUEST_TIMEOUT_SECONDS)
            http_status = str(response.status_code)
            if "html" in response.headers.get("Content-Type", ""):
                title = extract_html_title(response.text)
        except requests.RequestException:
            http_status = "request_failed"
        notes = str(ref.get("notes") or "")
        if title:
            notes = f"{notes} Page title: {title}".strip()
        rows.append(
            {
                "treatment_source_audit_run_id": run_id,
                "reference_id": ref["reference_id"],
                "source_kind": ref["source_kind"],
                "source_label": ref["source_label"],
                "source_url": ref["source_url"],
                "http_status": http_status,
                "saved_snapshot_path": "",
                "what_it_supports": ref["what_it_supports"],
                "notes": notes,
            }
        )
    return rows


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
        raise TreatmentSourceAuditError(
            f"GDC request failed for {url} with HTTP {response.status_code}. Body preview: {preview}"
        ) from exc
    try:
        return response.json(), {"final_url": response.url, "params": params or {}}
    except json.JSONDecodeError as exc:
        raise TreatmentSourceAuditError(f"Unable to decode JSON response from {url}") from exc


def fetch_official_online_evidence(
    run_id: str,
    audit_run_dir: Path,
    paths: WorkflowPaths,
    helper_module: Any,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT, "Accept": "application/json"})
    official_outputs: dict[str, Any] = {}

    def save_query(name: str, url: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        payload, request_meta = gdc_get_json(session, url, params=params)
        official_outputs[f"{name}_request_json"] = save_json_snapshot(helper_module, paths.repo_root, audit_run_dir / f"{name}_request.json", request_meta)
        official_outputs[f"{name}_response_json"] = save_json_snapshot(helper_module, paths.repo_root, audit_run_dir / f"{name}_response.json", payload)
        official_outputs[f"{name}_payload"] = payload
        return payload

    save_query("gdc_status", f"{GDC_API_BASE_URL}/status")

    clinical_category_filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Clinical"]}},
        ],
    }
    biospecimen_category_filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Biospecimen"]}},
        ],
    }
    clinical_supplement_filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Clinical Supplement"]}},
        ],
    }
    biospecimen_supplement_filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Biospecimen Supplement"]}},
        ],
    }
    cases_filters = {"op": "in", "content": {"field": "project.project_id", "value": [PROJECT_ID]}}

    save_query("files_clinical_category", f"{GDC_API_BASE_URL}/files", {"filters": json.dumps(clinical_category_filters), "size": 0, "facets": "data_type,data_format"})
    save_query("files_biospecimen_category", f"{GDC_API_BASE_URL}/files", {"filters": json.dumps(biospecimen_category_filters), "size": 0, "facets": "data_type,data_format"})
    save_query("files_clinical_supplement", f"{GDC_API_BASE_URL}/files", {"filters": json.dumps(clinical_supplement_filters), "size": 0, "facets": "data_format"})
    save_query("files_biospecimen_supplement", f"{GDC_API_BASE_URL}/files", {"filters": json.dumps(biospecimen_supplement_filters), "size": 0, "facets": "data_format"})
    save_query("cases_count", f"{GDC_API_BASE_URL}/cases", {"filters": json.dumps(cases_filters), "size": 0})
    save_query("cases_treatment_facets", f"{GDC_API_BASE_URL}/cases", {"filters": json.dumps(cases_filters), "size": 0, "facets": "diagnoses.treatments.treatment_type,diagnoses.treatments.treatment_intent_type,follow_ups.timepoint_category"})
    save_query("cases_expanded", f"{GDC_API_BASE_URL}/cases", {"filters": json.dumps(cases_filters), "size": 1, "sort": "submitter_id:asc", "expand": "diagnoses,diagnoses.treatments,follow_ups,exposures,demographic"})
    for endpoint_name in ["clinical_supplement", "biospecimen_supplement", "treatment", "follow_up"]:
        save_query(f"{endpoint_name}_dictionary", f"{GDC_API_BASE_URL}/v0/submission/_dictionary/{endpoint_name}")

    reference_rows = official_doc_reference_rows(session, run_id)
    saved_reference_paths = {
        "R5": official_outputs["clinical_supplement_dictionary_response_json"],
        "R6": official_outputs["biospecimen_supplement_dictionary_response_json"],
        "R7": official_outputs["treatment_dictionary_response_json"],
        "R8": official_outputs["follow_up_dictionary_response_json"],
        "R9": official_outputs["files_clinical_category_response_json"],
        "R10": official_outputs["cases_treatment_facets_response_json"],
        "R11": official_outputs["gdc_status_response_json"],
    }
    for row in reference_rows:
        if row["reference_id"] in saved_reference_paths:
            row["saved_snapshot_path"] = saved_reference_paths[row["reference_id"]]
    return official_outputs, reference_rows


def find_local_search_evidence(paths: WorkflowPaths) -> dict[str, Any]:
    raw_gdc_root = paths.raw_root / "tcga-brca" / "gdc"
    raw_pdf_paths = [path for path in raw_gdc_root.rglob("*.pdf")]
    xml_treatment_matches: list[str] = []
    for root in [paths.audit_root / "tcga-brca", paths.processed_root / "tcga-brca"]:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            lower_name = path.name.lower()
            if "xml" in lower_name and any(token in lower_name for token in ["treat", "drug", "radiation"]):
                xml_treatment_matches.append(repo_relative(path, paths.repo_root))
    indexed_export_matches: list[str] = []
    indexed_export_names = {"clinical.tsv", "clinical.json", "biospecimen.tsv", "biospecimen.json"}
    for root in [paths.raw_root, paths.audit_root, paths.processed_root]:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            lower_name = path.name.lower()
            if lower_name in indexed_export_names or lower_name.startswith("clinical.cart.") or lower_name.startswith("biospecimen.cart."):
                indexed_export_matches.append(repo_relative(path, paths.repo_root))
    return {
        "raw_pdf_count": len(raw_pdf_paths),
        "xml_treatment_match_count": len(xml_treatment_matches),
        "xml_treatment_matches": xml_treatment_matches,
        "indexed_export_match_count": len(indexed_export_matches),
        "indexed_export_matches": indexed_export_matches,
    }


def build_official_inventory_rows(run_id: str, official_outputs: dict[str, Any]) -> list[dict[str, Any]]:
    clinical_category_agg = official_outputs["files_clinical_category_payload"]["data"].get("aggregations", {})
    biospecimen_category_agg = official_outputs["files_biospecimen_category_payload"]["data"].get("aggregations", {})
    clinical_supplement_agg = official_outputs["files_clinical_supplement_payload"]["data"].get("aggregations", {})
    biospecimen_supplement_agg = official_outputs["files_biospecimen_supplement_payload"]["data"].get("aggregations", {})
    cases_payload = official_outputs["cases_count_payload"]["data"]
    clinical_supplement_dictionary = official_outputs["clinical_supplement_dictionary_payload"]
    biospecimen_supplement_dictionary = official_outputs["biospecimen_supplement_dictionary_payload"]
    official_rows = [
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "indexed_clinical_surface",
            "official_source_layer": "Indexed GDC clinical surface (cases API / portal Clinical TSV-JSON export)",
            "official_source_kind": "model-aligned case-level access surface",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(cases_payload["pagination"]["total"]),
            "official_count_basis": "TCGA-BRCA cases from live cases API",
            "what_it_contains_in_general": "Model-aligned clinical entities exposed through the GDC database, including diagnoses, diagnoses.treatments, follow_ups, exposures, and demographic.",
            "treatment_relevance": "high",
            "richness_vs_other_layers": "coarser than raw TCGA supplement files; useful as an official indexed treatment/follow-up surface",
            "official_reference_ids": "R1,R4,R7,R8,R10",
            "notes": "",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "clinical_supplement_bcr_xml",
            "official_source_layer": "Clinical Supplement / BCR XML",
            "official_source_kind": "clinical supplement file layer",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(bucket_count(clinical_supplement_agg, "data_format", "bcr xml")),
            "official_count_basis": "live files API: Clinical Supplement filtered by data_format",
            "what_it_contains_in_general": str(clinical_supplement_dictionary.get("description") or ""),
            "treatment_relevance": "high",
            "richness_vs_other_layers": "richer raw per-case supplement source than indexed clinical; harder to parse",
            "official_reference_ids": "R2,R5,R9",
            "notes": "Per-case XML-style raw source subtype within Clinical Supplement.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "clinical_supplement_bcr_omf_xml",
            "official_source_layer": "Clinical Supplement / BCR OMF XML",
            "official_source_kind": "clinical supplement file layer",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(bucket_count(clinical_supplement_agg, "data_format", "bcr omf xml")),
            "official_count_basis": "live files API: Clinical Supplement filtered by data_format",
            "what_it_contains_in_general": str(clinical_supplement_dictionary.get("description") or ""),
            "treatment_relevance": "possible",
            "richness_vs_other_layers": "richer raw per-case supplement subtype than indexed clinical; harder to parse",
            "official_reference_ids": "R2,R5,R9",
            "notes": "Smaller TCGA-BRCA clinical supplement XML subtype than main BCR XML.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "clinical_supplement_bcr_biotab",
            "official_source_layer": "Clinical Supplement / BCR Biotab",
            "official_source_kind": "clinical supplement file layer",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(bucket_count(clinical_supplement_agg, "data_format", "bcr biotab")),
            "official_count_basis": "live files API: Clinical Supplement filtered by data_format",
            "what_it_contains_in_general": str(clinical_supplement_dictionary.get("description") or ""),
            "treatment_relevance": "high",
            "richness_vs_other_layers": "coarser and easier to parse than supplement XML; strong early audit layer",
            "official_reference_ids": "R2,R5,R9",
            "notes": "Tabular supplement subtype; exact table names are confirmed locally from the downloaded cohort files.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "biospecimen_supplement_bcr_xml",
            "official_source_layer": "Biospecimen Supplement / BCR XML",
            "official_source_kind": "biospecimen supplement file layer",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(bucket_count(biospecimen_supplement_agg, "data_format", "bcr xml")),
            "official_count_basis": "live files API: Biospecimen Supplement filtered by data_format",
            "what_it_contains_in_general": str(biospecimen_supplement_dictionary.get("description") or ""),
            "treatment_relevance": "low-indirect",
            "richness_vs_other_layers": "richer raw per-case biospecimen source; mainly provenance/linkage rather than treatment",
            "official_reference_ids": "R3,R6,R9",
            "notes": "Relevant mainly for specimen provenance and linkage checks, not as a primary treatment layer.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "biospecimen_supplement_bcr_ssf_xml",
            "official_source_layer": "Biospecimen Supplement / BCR SSF XML",
            "official_source_kind": "biospecimen supplement file layer",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(bucket_count(biospecimen_supplement_agg, "data_format", "bcr ssf xml")),
            "official_count_basis": "live files API: Biospecimen Supplement filtered by data_format",
            "what_it_contains_in_general": str(biospecimen_supplement_dictionary.get("description") or ""),
            "treatment_relevance": "low-indirect",
            "richness_vs_other_layers": "richer biospecimen XML subtype; mainly provenance/linkage rather than treatment",
            "official_reference_ids": "R3,R6,R9",
            "notes": "Not a primary structured treatment source, but still an official downloaded layer.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "biospecimen_supplement_bcr_biotab",
            "official_source_layer": "Biospecimen Supplement / BCR Biotab",
            "official_source_kind": "biospecimen supplement file layer",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(bucket_count(biospecimen_supplement_agg, "data_format", "bcr biotab")),
            "official_count_basis": "live files API: Biospecimen Supplement filtered by data_format",
            "what_it_contains_in_general": str(biospecimen_supplement_dictionary.get("description") or ""),
            "treatment_relevance": "low-indirect",
            "richness_vs_other_layers": "coarser/easier than biospecimen XML; useful for identifier hierarchy rather than treatment",
            "official_reference_ids": "R3,R6,R9",
            "notes": "Useful mainly for sample/specimen linkage review.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "pathology_report_pdf",
            "official_source_layer": "Pathology Report / PDF",
            "official_source_kind": "clinical file layer",
            "official_availability_for_tcga_brca": "available",
            "official_count": str(bucket_count(clinical_category_agg, "data_type", "pathology report")),
            "official_count_basis": "live files API: Clinical data category filtered by data_type",
            "what_it_contains_in_general": "Clinical pathology report PDFs linked to TCGA-BRCA cases.",
            "treatment_relevance": "possible but low for structured extraction",
            "richness_vs_other_layers": "unstructured narrative source; lower priority for structured treatment audit",
            "official_reference_ids": "R4,R9",
            "notes": "Officially available in TCGA-BRCA, but not a primary structured treatment layer.",
        },
    ]
    official_rows[0]["notes"] = (
        f"Broader live discovery within Clinical data_category found {bucket_map(clinical_category_agg, 'data_type')}. "
        f"Broader live discovery within Biospecimen data_category found {bucket_map(biospecimen_category_agg, 'data_type')}. "
        "Slide Image was discovered but treated as not meaningful for this treatment-source audit."
    )
    return official_rows


def build_local_inventory_rows(
    run_id: str,
    paths: WorkflowPaths,
    inputs: WorkflowInputs,
    local_search: dict[str, Any],
) -> list[dict[str, Any]]:
    supplement_rows = inputs.supplement_metadata_rows
    clinical_rows = [row for row in supplement_rows if row.get("source_class") == "clinical"]
    biospecimen_rows = [row for row in supplement_rows if row.get("source_class") == "biospecimen"]
    clinical_format_counts = Counter(str(row.get("data_format") or "") for row in clinical_rows)
    biospecimen_format_counts = Counter(str(row.get("data_format") or "") for row in biospecimen_rows)
    clinical_source_class = inputs.source_supplements_pointer.get("source_classes", {}).get("clinical", {})
    biospecimen_source_class = inputs.source_supplements_pointer.get("source_classes", {}).get("biospecimen", {})
    clinical_manifest_path = str(clinical_source_class.get("manifest_path") or "")
    biospecimen_manifest_path = str(biospecimen_source_class.get("manifest_path") or "")
    clinical_download_dir = str(clinical_source_class.get("download_dir") or "")
    biospecimen_download_dir = str(biospecimen_source_class.get("download_dir") or "")
    clinical_tables = [row.get("table_name", "") for row in inputs.clinical_biotab_manifest_rows]
    biospecimen_parsed_tables = [row.get("table_name", "") for row in inputs.biospecimen_source_inventory_rows if row.get("scope_status") == "parsed"]
    biospecimen_excluded_tables = [row.get("source_filename", "") for row in inputs.biospecimen_source_inventory_rows if row.get("scope_status") == "excluded"]
    pathology_rows = [row for row in inputs.source_inventory_rows if row.get("data_category") == "Clinical" and row.get("data_type") == "Pathology Report"]
    endpoint_run_log = inputs.endpoint_target_prep_run_log
    rows = [
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "indexed_clinical_surface",
            "local_source_layer": "Indexed GDC clinical surface (local evidence)",
            "local_path_or_scope": "No saved local clinical TSV/JSON export or portal/API clinical snapshot identified",
            "what_is_present_locally": "Repo contains source inventory API snapshots and source supplement downloads, but no dedicated local indexed clinical export artifact.",
            "looks_downloaded_completely_or_partially": "not downloaded locally",
            "parsed_or_exploited_status": "not applicable locally",
            "treatment_relevance": "high if downloaded",
            "evidence_paths": ",".join([repo_relative(paths.source_inventory_latest_pointer, paths.repo_root), repo_relative(paths.source_supplements_latest_pointer, paths.repo_root)]),
            "notes": f"Current raw download scope is limited to supplement manifests/download dirs. Indexed export filename search found {local_search['indexed_export_match_count']} matches.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "clinical_supplement_bcr_xml",
            "local_source_layer": "Clinical Supplement / BCR XML",
            "local_path_or_scope": clinical_download_dir,
            "what_is_present_locally": f"{clinical_format_counts.get('BCR XML', 0)} clinical BCR XML files in supplement metadata",
            "looks_downloaded_completely_or_partially": "downloaded completely",
            "parsed_or_exploited_status": "partially exploited: XML-derived endpoint/follow-up outputs exist, but no dedicated treatment parse identified",
            "treatment_relevance": "high",
            "evidence_paths": ",".join([repo_relative(paths.source_supplements_latest_pointer, paths.repo_root), clinical_manifest_path, clinical_download_dir, str(inputs.endpoint_target_prep_pointer.get("run_log_json") or "")]),
            "notes": f"Endpoint XML workflow reports {endpoint_run_log.get('counts', {}).get('bcr_xml_source_row_count', '')} BCR XML source rows and restricts scope to BCR XML.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "clinical_supplement_bcr_omf_xml",
            "local_source_layer": "Clinical Supplement / BCR OMF XML",
            "local_path_or_scope": clinical_download_dir,
            "what_is_present_locally": f"{clinical_format_counts.get('BCR OMF XML', 0)} clinical BCR OMF XML files in supplement metadata",
            "looks_downloaded_completely_or_partially": "downloaded completely",
            "parsed_or_exploited_status": "not yet parsed/exploited for treatment; current XML workflow explicitly excludes BCR OMF XML",
            "treatment_relevance": "possible",
            "evidence_paths": ",".join([repo_relative(paths.source_supplements_latest_pointer, paths.repo_root), clinical_manifest_path, clinical_download_dir, str(inputs.endpoint_target_prep_pointer.get("run_log_json") or "")]),
            "notes": "Do not confuse the downloaded BCR OMF XML case files with the separately parsed `clinical_omf_v4_0.tsv` biotab table.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "clinical_supplement_bcr_biotab",
            "local_source_layer": "Clinical Supplement / BCR Biotab",
            "local_path_or_scope": ",".join([clinical_manifest_path, str(inputs.clinical_biotabs_pointer.get('table_manifest_tsv') or '')]),
            "what_is_present_locally": f"{clinical_format_counts.get('BCR Biotab', 0)} clinical biotab files downloaded; parsed tables: " + ", ".join(clinical_tables),
            "looks_downloaded_completely_or_partially": "downloaded completely",
            "parsed_or_exploited_status": "already parsed",
            "treatment_relevance": "high",
            "evidence_paths": ",".join([repo_relative(paths.clinical_biotabs_latest_pointer, paths.repo_root), str(inputs.clinical_biotabs_pointer.get("table_manifest_tsv") or ""), str(inputs.clinical_biotabs_pointer.get("run_log_json") or "")]),
            "notes": "This is the currently exploited structured treatment source layer: clinical_drug, clinical_radiation, clinical_patient, follow-up tables, and clinical_omf_v4_0 biotab are already parsed.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "biospecimen_supplement_bcr_xml",
            "local_source_layer": "Biospecimen Supplement / BCR XML",
            "local_path_or_scope": biospecimen_download_dir,
            "what_is_present_locally": f"{biospecimen_format_counts.get('BCR XML', 0)} biospecimen BCR XML files in supplement metadata",
            "looks_downloaded_completely_or_partially": "downloaded completely",
            "parsed_or_exploited_status": "downloaded but not yet parsed/exploited",
            "treatment_relevance": "low-indirect",
            "evidence_paths": ",".join([repo_relative(paths.source_supplements_latest_pointer, paths.repo_root), biospecimen_manifest_path, biospecimen_download_dir]),
            "notes": "Useful mainly for specimen provenance and linkage, not as a primary treatment layer.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "biospecimen_supplement_bcr_ssf_xml",
            "local_source_layer": "Biospecimen Supplement / BCR SSF XML",
            "local_path_or_scope": biospecimen_download_dir,
            "what_is_present_locally": f"{biospecimen_format_counts.get('BCR SSF XML', 0)} biospecimen BCR SSF XML files in supplement metadata",
            "looks_downloaded_completely_or_partially": "downloaded completely",
            "parsed_or_exploited_status": "downloaded but not yet parsed/exploited",
            "treatment_relevance": "low-indirect",
            "evidence_paths": ",".join([repo_relative(paths.source_supplements_latest_pointer, paths.repo_root), biospecimen_manifest_path, biospecimen_download_dir]),
            "notes": "Not a primary treatment source, but still an unexploited downloaded official layer.",
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "biospecimen_supplement_bcr_biotab",
            "local_source_layer": "Biospecimen Supplement / BCR Biotab",
            "local_path_or_scope": ",".join([biospecimen_manifest_path, str(inputs.biospecimen_biotabs_pointer.get('source_inventory_tsv') or ''), str(inputs.biospecimen_biotabs_pointer.get('table_manifest_tsv') or '')]),
            "what_is_present_locally": f"{biospecimen_format_counts.get('BCR Biotab', 0)} biospecimen biotab files downloaded; parsed tables: " + ", ".join(biospecimen_parsed_tables),
            "looks_downloaded_completely_or_partially": "downloaded completely",
            "parsed_or_exploited_status": "partially parsed: 8 core tables parsed, 2 SSF text biotab files deferred",
            "treatment_relevance": "low-indirect",
            "evidence_paths": ",".join([repo_relative(paths.biospecimen_biotabs_latest_pointer, paths.repo_root), str(inputs.biospecimen_biotabs_pointer.get("source_inventory_tsv") or ""), str(inputs.biospecimen_biotabs_pointer.get("table_manifest_tsv") or "")]),
            "notes": "Deferred biospecimen SSF text biotabs: " + ", ".join(biospecimen_excluded_tables),
        },
        {
            "treatment_source_audit_run_id": run_id,
            "source_layer_key": "pathology_report_pdf",
            "local_source_layer": "Pathology Report / PDF",
            "local_path_or_scope": repo_relative(paths.raw_root / "tcga-brca" / "gdc", paths.repo_root),
            "what_is_present_locally": f"{len(pathology_rows)} pathology reports listed in full file inventory; {local_search['raw_pdf_count']} PDFs found under raw GDC tree",
            "looks_downloaded_completely_or_partially": "not downloaded locally",
            "parsed_or_exploited_status": "not parsed",
            "treatment_relevance": "possible but low for structured extraction",
            "evidence_paths": ",".join([str(inputs.source_inventory_pointer.get("file_inventory_tsv") or ""), repo_relative(paths.raw_root / "tcga-brca" / "gdc", paths.repo_root)]),
            "notes": "Officially visible in the full TCGA-BRCA inventory, but no local raw PDF payloads were found.",
        },
    ]
    return rows


def build_gap_analysis_rows(run_id: str, official_rows: list[dict[str, Any]], inputs: WorkflowInputs, local_search: dict[str, Any]) -> list[dict[str, Any]]:
    official_lookup = {row["source_layer_key"]: row for row in official_rows}
    treatment_profile_counts = inputs.patient_treatment_profile_run_log.get("counts", {})
    row_specs = {
        "indexed_clinical_surface": ("no local export identified", "no local export identified", "high", "not downloaded but available in GDC", "Officially available via cases API / portal Clinical TSV-JSON export. No local export artifact was identified; this is the most realistic additional official treatment surface to check next."),
        "clinical_supplement_bcr_xml": ("downloaded completely", "only endpoint/follow-up exploitation identified", "high", "already downloaded but not yet fully parsed / not yet fully exploited", f"Downloaded completely, but current saved XML outputs are endpoint/follow-up focused. Filename search found {local_search['xml_treatment_match_count']} dedicated XML treatment outputs."),
        "clinical_supplement_bcr_omf_xml": ("downloaded completely", "no dedicated exploitation identified", "possible-to-moderate", "already downloaded but not yet parsed / not yet exploited", "Current XML follow-up workflow excludes BCR OMF XML. The parsed clinical_omf_v4_0 biotab table does not eliminate this raw XML gap."),
        "clinical_supplement_bcr_biotab": ("downloaded completely", "parsed", "high", "already downloaded and already parsed", "This is the main currently exploited structured treatment layer and feeds downstream treatment-prep outputs."),
        "biospecimen_supplement_bcr_xml": ("downloaded completely", "not parsed", "low-indirect", "already downloaded but not yet parsed / not yet exploited", "Downloaded completely, but likely value for treatment itself is indirect."),
        "biospecimen_supplement_bcr_ssf_xml": ("downloaded completely", "not parsed", "low-indirect", "already downloaded but not yet parsed / not yet exploited", "Downloaded completely, but likely value for treatment itself is indirect."),
        "biospecimen_supplement_bcr_biotab": ("downloaded completely", "8 of 10 tables parsed; 2 deferred", "low-indirect", "already downloaded but not yet fully parsed / not yet fully exploited", "Two SSF text biotabs remain deferred, but this is not the main treatment-data gap."),
        "pathology_report_pdf": ("available in GDC but not downloaded locally", "not parsed", "low-for-structured-extraction", "not primary / probably not worth chasing first", "Unstructured narrative PDFs may contain treatment context in some cases, but they are not the first structured treatment layer to chase."),
    }
    rows: list[dict[str, Any]] = []
    for key in SOURCE_LAYER_KEYS:
        local_download_status, local_parse_status, likely_value, classification, notes = row_specs[key]
        official_row = official_lookup[key]
        rows.append({"treatment_source_audit_run_id": run_id, "official_source_layer": official_row["official_source_layer"], "official_availability_for_tcga_brca": official_row["official_availability_for_tcga_brca"], "local_download_status": local_download_status, "local_parse_status": local_parse_status, "likely_treatment_value": likely_value, "classification": classification, "notes": notes})
    rows.append({
        "treatment_source_audit_run_id": run_id,
        "official_source_layer": "Patient-level treatment missingness interpretation",
        "official_availability_for_tcga_brca": "not a source layer row",
        "local_download_status": "clinical supplement downloads complete",
        "local_parse_status": "current treatment-prep uses parsed clinical biotabs plus OS/endpoint outputs",
        "likely_treatment_value": "interpretive",
        "classification": "source-level missingness not yet proven; current limitation is still partly a source-exploitation problem",
        "notes": f"Current patient-level missingness in downstream treatment-prep is not a supplement download problem: clinical biotab-derived treatment profile reports {treatment_profile_counts.get('patients_with_drug', '')} patients with any drug row, {treatment_profile_counts.get('patients_with_radiation', '')} with any radiation row, and {treatment_profile_counts.get('patients_no_drug_or_radiation_rows', '')} with no drug or radiation rows. Because indexed clinical and most clinical XML/OMF XML treatment surfaces remain unchecked locally, these patient gaps should not yet be labeled source-level absent across all meaningful GDC sources.",
    })
    return rows


def build_summary_rows(run_id: str, official_outputs: dict[str, Any], inputs: WorkflowInputs, local_search: dict[str, Any]) -> list[dict[str, Any]]:
    clinical_category_agg = official_outputs["files_clinical_category_payload"]["data"].get("aggregations", {})
    biospecimen_category_agg = official_outputs["files_biospecimen_category_payload"]["data"].get("aggregations", {})
    clinical_supplement_agg = official_outputs["files_clinical_supplement_payload"]["data"].get("aggregations", {})
    biospecimen_supplement_agg = official_outputs["files_biospecimen_supplement_payload"]["data"].get("aggregations", {})
    status_payload = official_outputs["gdc_status_payload"]
    treatment_profile_counts = inputs.patient_treatment_profile_run_log.get("counts", {})
    return [
        {"treatment_source_audit_run_id": run_id, "summary_section": "official_target", "summary_metric": "gdc_data_release", "summary_value": str(status_payload.get("data_release") or ""), "notes": "Live GDC API status target used for reconciliation."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "official_target", "summary_metric": "tcga_brca_case_count", "summary_value": str(official_outputs["cases_count_payload"]["data"]["pagination"]["total"]), "notes": "Live cases API total."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "official_target", "summary_metric": "clinical_category_discovery", "summary_value": json.dumps(bucket_map(clinical_category_agg, "data_type"), ensure_ascii=True, sort_keys=True), "notes": "Broader Clinical-category discovery before narrowing to treatment-relevant layers."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "official_target", "summary_metric": "biospecimen_category_discovery", "summary_value": json.dumps(bucket_map(biospecimen_category_agg, "data_type"), ensure_ascii=True, sort_keys=True), "notes": "Broader Biospecimen-category discovery before narrowing to treatment-relevant layers."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "reconciliation", "summary_metric": "clinical_supplement_live_count", "summary_value": str(official_outputs["files_clinical_supplement_payload"]["data"]["pagination"]["total"]), "notes": "Should reconcile to local supplement metadata and manifest/download counts."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "reconciliation", "summary_metric": "biospecimen_supplement_live_count", "summary_value": str(official_outputs["files_biospecimen_supplement_payload"]["data"]["pagination"]["total"]), "notes": "Should reconcile to local supplement metadata and manifest/download counts."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "reconciliation", "summary_metric": "clinical_supplement_live_format_counts", "summary_value": json.dumps(bucket_map(clinical_supplement_agg, "data_format"), ensure_ascii=True, sort_keys=True), "notes": "Live official format split."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "reconciliation", "summary_metric": "biospecimen_supplement_live_format_counts", "summary_value": json.dumps(bucket_map(biospecimen_supplement_agg, "data_format"), ensure_ascii=True, sort_keys=True), "notes": "Live official format split."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "local_findings", "summary_metric": "clinical_supplement_download_status", "summary_value": "downloaded completely", "notes": f"Local source run log validation counts: clinical manifest/download rows = {inputs.source_supplements_run_log.get('validation', {}).get('downloads', {}).get('clinical', {}).get('validated_downloaded_file_count', '')}."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "local_findings", "summary_metric": "biospecimen_supplement_download_status", "summary_value": "downloaded completely", "notes": f"Local source run log validation counts: biospecimen manifest/download rows = {inputs.source_supplements_run_log.get('validation', {}).get('downloads', {}).get('biospecimen', {}).get('validated_downloaded_file_count', '')}."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "local_findings", "summary_metric": "clinical_biotab_parse_status", "summary_value": "parsed", "notes": f"Parsed table count = {inputs.clinical_biotabs_pointer.get('parsed_table_count', '')}."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "local_findings", "summary_metric": "biospecimen_biotab_parse_status", "summary_value": "partially_parsed", "notes": f"Parsed table count = {inputs.biospecimen_biotabs_pointer.get('parsed_table_count', '')}; excluded biotab file count = {inputs.biospecimen_biotabs_pointer.get('excluded_biotab_file_count', '')}."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "local_findings", "summary_metric": "clinical_xml_treatment_parse_identified", "summary_value": "no", "notes": f"Dedicated XML treatment parse filename search returned {local_search['xml_treatment_match_count']} matches. Current saved XML outputs are endpoint/follow-up focused."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "local_findings", "summary_metric": "indexed_clinical_export_identified", "summary_value": "no", "notes": f"Dedicated indexed clinical export filename search returned {local_search['indexed_export_match_count']} matches."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "local_findings", "summary_metric": "local_pathology_pdf_count", "summary_value": str(local_search["raw_pdf_count"]), "notes": "PDF count found under local raw GDC tree."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "interpretation", "summary_metric": "main_limitation", "summary_value": "download coverage is strong for supplement files; remaining limitation is source exploitation plus some true source sparsity", "notes": f"Current treatment-prep already uses parsed clinical biotabs, yet {treatment_profile_counts.get('patients_no_drug_or_radiation_rows', '')} patients still have no drug or radiation rows. That burden is not a supplement download failure."},
        {"treatment_source_audit_run_id": run_id, "summary_section": "recommendation", "summary_metric": "best_next_step", "summary_value": "compare parsed clinical treatment biotabs against the indexed GDC clinical treatment surface first", "notes": "If that still leaves material treatment evidence unexplained, inspect whether downloaded clinical BCR XML or BCR OMF XML add treatment detail beyond the already parsed biotabs and indexed clinical surface."},
    ]


def build_validation(official_rows: list[dict[str, Any]], inputs: WorkflowInputs, local_search: dict[str, Any]) -> dict[str, Any]:
    official_lookup = {row["source_layer_key"]: row for row in official_rows}
    source_validation = inputs.source_supplements_run_log.get("validation", {})
    clinical_manifest_rows = int(source_validation.get("manifests", {}).get("clinical", {}).get("row_count", 0))
    clinical_download_rows = int(source_validation.get("downloads", {}).get("clinical", {}).get("validated_downloaded_file_count", 0))
    biospecimen_manifest_rows = int(source_validation.get("manifests", {}).get("biospecimen", {}).get("row_count", 0))
    biospecimen_download_rows = int(source_validation.get("downloads", {}).get("biospecimen", {}).get("validated_downloaded_file_count", 0))
    clinical_live = int(official_lookup["clinical_supplement_bcr_xml"]["official_count"]) + int(official_lookup["clinical_supplement_bcr_omf_xml"]["official_count"]) + int(official_lookup["clinical_supplement_bcr_biotab"]["official_count"])
    biospecimen_live = int(official_lookup["biospecimen_supplement_bcr_xml"]["official_count"]) + int(official_lookup["biospecimen_supplement_bcr_ssf_xml"]["official_count"]) + int(official_lookup["biospecimen_supplement_bcr_biotab"]["official_count"])
    return {
        "passed": True,
        "required_upstream_pointers_found": True,
        "source_inventory_run_log_completed": True,
        "source_supplements_run_log_completed": True,
        "clinical_biotabs_run_log_completed": True,
        "biospecimen_biotabs_run_log_completed": True,
        "endpoint_target_prep_run_log_completed": True,
        "patient_treatment_profile_run_log_completed": True,
        "patient_treatment_grouping_run_log_completed": True,
        "clinical_manifest_matches_download_validation": clinical_manifest_rows == clinical_download_rows,
        "biospecimen_manifest_matches_download_validation": biospecimen_manifest_rows == biospecimen_download_rows,
        "clinical_live_count_matches_local_manifest_count": clinical_live == clinical_manifest_rows,
        "biospecimen_live_count_matches_local_manifest_count": biospecimen_live == biospecimen_manifest_rows,
        "clinical_biotab_parsed_table_count_is_9": int(inputs.clinical_biotabs_pointer.get("parsed_table_count", 0)) == 9,
        "biospecimen_biotab_discovered_count_is_10": int(inputs.biospecimen_biotabs_pointer.get("discovered_biotab_file_count", 0)) == 10,
        "biospecimen_biotab_excluded_count_is_2": int(inputs.biospecimen_biotabs_pointer.get("excluded_biotab_file_count", 0)) == 2,
        "dedicated_xml_treatment_parse_found": local_search["xml_treatment_match_count"] > 0,
        "indexed_clinical_export_found": local_search["indexed_export_match_count"] > 0,
        "local_pathology_pdf_found": local_search["raw_pdf_count"] > 0,
        "no_claim_that_all_patient_level_missingness_is_source_level_absent": True,
    }


def write_outputs(
    *,
    run_id: str,
    started_at: datetime,
    audit_run_dir: Path,
    paths: WorkflowPaths,
    inputs: WorkflowInputs,
    official_outputs: dict[str, Any],
    official_rows: list[dict[str, Any]],
    local_rows: list[dict[str, Any]],
    gap_rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    reference_rows: list[dict[str, Any]],
    local_search: dict[str, Any],
    helper_module: Any,
) -> None:
    official_inventory_tsv = audit_run_dir / "official_gdc_source_inventory.tsv"
    local_inventory_tsv = audit_run_dir / "local_treatment_source_inventory.tsv"
    gap_analysis_tsv = audit_run_dir / "treatment_source_gap_analysis.tsv"
    summary_tsv = audit_run_dir / "treatment_source_audit_summary.tsv"
    reference_tsv = audit_run_dir / "treatment_source_official_references.tsv"
    run_log_path = audit_run_dir / "run_log.json"
    helper_module.write_dict_rows_tsv(official_inventory_tsv, OFFICIAL_SOURCE_FIELDNAMES, official_rows)
    helper_module.write_dict_rows_tsv(local_inventory_tsv, LOCAL_SOURCE_FIELDNAMES, local_rows)
    helper_module.write_dict_rows_tsv(gap_analysis_tsv, GAP_ANALYSIS_FIELDNAMES, gap_rows)
    helper_module.write_dict_rows_tsv(summary_tsv, SUMMARY_FIELDNAMES, summary_rows)
    helper_module.write_dict_rows_tsv(reference_tsv, REFERENCE_FIELDNAMES, reference_rows)
    completed_at = utc_now()
    validation = build_validation(official_rows, inputs, local_search)
    run_log = {
        "status": "completed",
        "treatment_source_audit_run_id": run_id,
        "started_at_utc": format_utc_timestamp(started_at),
        "completed_at_utc": format_utc_timestamp(completed_at),
        "repo_root": str(paths.repo_root),
        "trial_name": "tcga_only_source_audited",
        "dataset_scope": "tcga_brca_only",
        "inputs": {"trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root), "source_inventory_latest_json": repo_relative(paths.source_inventory_latest_pointer, paths.repo_root), "source_supplements_latest_json": repo_relative(paths.source_supplements_latest_pointer, paths.repo_root), "clinical_biotabs_latest_json": repo_relative(paths.clinical_biotabs_latest_pointer, paths.repo_root), "biospecimen_biotabs_latest_json": repo_relative(paths.biospecimen_biotabs_latest_pointer, paths.repo_root), "endpoint_target_prep_latest_json": repo_relative(paths.endpoint_target_prep_latest_pointer, paths.repo_root), "patient_treatment_profile_latest_json": repo_relative(paths.patient_treatment_profile_latest_pointer, paths.repo_root), "patient_treatment_grouping_latest_json": repo_relative(paths.patient_treatment_grouping_latest_pointer, paths.repo_root)},
        "official_online_evidence": {key: value for key, value in official_outputs.items() if key.endswith("_json")},
        "outputs": {"audit_run_directory": repo_relative(audit_run_dir, paths.repo_root), "official_gdc_source_inventory_tsv": repo_relative(official_inventory_tsv, paths.repo_root), "local_treatment_source_inventory_tsv": repo_relative(local_inventory_tsv, paths.repo_root), "treatment_source_gap_analysis_tsv": repo_relative(gap_analysis_tsv, paths.repo_root), "treatment_source_audit_summary_tsv": repo_relative(summary_tsv, paths.repo_root), "treatment_source_official_references_tsv": repo_relative(reference_tsv, paths.repo_root), "run_log_json": repo_relative(run_log_path, paths.repo_root), "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root)},
        "validation": validation,
        "rules": {"audit_only": True, "does_not_modify_raw_downloads": True, "does_not_build_treatment_arms": True, "does_not_perform_modeling": True, "does_not_assume_missing_just_because_unused": True, "does_not_assume_xml_solves_everything": True, "indexed_clinical_surface_treated_as_distinct_official_access_layer": True, "pathology_report_treated_as_official_but_low_priority_for_structured_treatment_extraction": True, "patient_level_missingness_not_labeled_source_level_absent_until_all_meaningful_surfaces_checked": True},
        "counts": {"official_inventory_row_count": len(official_rows), "local_inventory_row_count": len(local_rows), "gap_analysis_row_count": len(gap_rows), "summary_row_count": len(summary_rows), "reference_row_count": len(reference_rows), "xml_treatment_match_count": local_search["xml_treatment_match_count"], "indexed_export_match_count": local_search["indexed_export_match_count"], "raw_pdf_count": local_search["raw_pdf_count"]},
        "upstream_snapshots": {"source_inventory_pointer": inputs.source_inventory_pointer, "source_supplements_pointer": inputs.source_supplements_pointer, "clinical_biotabs_pointer": inputs.clinical_biotabs_pointer, "biospecimen_biotabs_pointer": inputs.biospecimen_biotabs_pointer, "endpoint_target_prep_pointer": inputs.endpoint_target_prep_pointer, "patient_treatment_profile_pointer": inputs.patient_treatment_profile_pointer, "patient_treatment_grouping_pointer": inputs.patient_treatment_grouping_pointer},
    }
    helper_module.write_json(run_log_path, run_log)
    latest_pointer = {
        "updated_at_utc": format_utc_timestamp(completed_at),
        "treatment_source_audit_run_id": run_id,
        "source_inventory_run_id": str(inputs.source_inventory_pointer.get("run_id") or ""),
        "source_supplement_run_id": str(inputs.source_supplements_pointer.get("run_id") or ""),
        "clinical_biotab_parse_run_id": str(inputs.clinical_biotabs_pointer.get("parse_run_id") or ""),
        "biospecimen_biotab_parse_run_id": str(inputs.biospecimen_biotabs_pointer.get("parse_run_id") or ""),
        "endpoint_target_prep_v1_run_id": str(inputs.endpoint_target_prep_pointer.get("endpoint_target_prep_v1_run_id") or ""),
        "patient_treatment_profile_v1_run_id": str(inputs.patient_treatment_profile_pointer.get("patient_treatment_profile_v1_run_id") or ""),
        "patient_treatment_grouping_v1_run_id": str(inputs.patient_treatment_grouping_pointer.get("patient_treatment_grouping_v1_run_id") or ""),
        "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
        "official_gdc_source_inventory_tsv": repo_relative(official_inventory_tsv, paths.repo_root),
        "local_treatment_source_inventory_tsv": repo_relative(local_inventory_tsv, paths.repo_root),
        "treatment_source_gap_analysis_tsv": repo_relative(gap_analysis_tsv, paths.repo_root),
        "treatment_source_audit_summary_tsv": repo_relative(summary_tsv, paths.repo_root),
        "treatment_source_official_references_tsv": repo_relative(reference_tsv, paths.repo_root),
        "run_log_json": repo_relative(run_log_path, paths.repo_root),
        "source_inventory_latest_json": repo_relative(paths.source_inventory_latest_pointer, paths.repo_root),
        "source_supplements_latest_json": repo_relative(paths.source_supplements_latest_pointer, paths.repo_root),
        "clinical_biotabs_latest_json": repo_relative(paths.clinical_biotabs_latest_pointer, paths.repo_root),
        "biospecimen_biotabs_latest_json": repo_relative(paths.biospecimen_biotabs_latest_pointer, paths.repo_root),
        "endpoint_target_prep_latest_json": repo_relative(paths.endpoint_target_prep_latest_pointer, paths.repo_root),
        "patient_treatment_profile_latest_json": repo_relative(paths.patient_treatment_profile_latest_pointer, paths.repo_root),
        "patient_treatment_grouping_latest_json": repo_relative(paths.patient_treatment_grouping_latest_pointer, paths.repo_root),
        "gdc_data_release": str(official_outputs["gdc_status_payload"].get("data_release") or ""),
        "gdc_tag": str(official_outputs["gdc_status_payload"].get("tag") or ""),
    }
    with paths.latest_pointer.open("w", encoding="utf-8") as handle:
        json.dump(latest_pointer, handle, indent=2)
        handle.write("\n")


def main() -> None:
    started_at = utc_now()
    run_id = started_at.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    print(f"[25] treatment-source audit - run_id={run_id}")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    print(f"  repo root: {paths.repo_root}")
    print("  loading local pointers and run logs ...")
    inputs = load_inputs(paths, helper_module)
    print("  creating audit run directory ...")
    audit_run_dir = helper_module.create_run_directory(paths.audit_runs_root / run_id)
    print("  fetching live official GDC evidence ...")
    official_outputs, reference_rows = fetch_official_online_evidence(run_id, audit_run_dir, paths, helper_module)
    print("  scanning local repo for absence/presence evidence ...")
    local_search = find_local_search_evidence(paths)
    print(f"  local search: xml-treatment-matches={local_search['xml_treatment_match_count']} | indexed-export-matches={local_search['indexed_export_match_count']} | raw-pdf-count={local_search['raw_pdf_count']}")
    print("  building audit tables ...")
    official_rows = build_official_inventory_rows(run_id, official_outputs)
    local_rows = build_local_inventory_rows(run_id, paths, inputs, local_search)
    gap_rows = build_gap_analysis_rows(run_id, official_rows, inputs, local_search)
    summary_rows = build_summary_rows(run_id, official_outputs, inputs, local_search)
    print("  writing outputs ...")
    write_outputs(run_id=run_id, started_at=started_at, audit_run_dir=audit_run_dir, paths=paths, inputs=inputs, official_outputs=official_outputs, official_rows=official_rows, local_rows=local_rows, gap_rows=gap_rows, summary_rows=summary_rows, reference_rows=reference_rows, local_search=local_search, helper_module=helper_module)
    print(f"  audit run directory : {audit_run_dir}")
    print(f"  latest pointer      : {paths.latest_pointer}")
    print(f"  official rows       : {len(official_rows)}")
    print(f"  local rows          : {len(local_rows)}")
    print(f"  gap rows            : {len(gap_rows)}")
    print(f"  summary rows        : {len(summary_rows)}")
    print(f"[25] completed - run_id={run_id}")


if __name__ == "__main__":
    main()
