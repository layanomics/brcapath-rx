#!/usr/bin/env python
"""Parse TCGA-BRCA clinical XML follow-up layers and prepare endpoint-target prep v1."""

from __future__ import annotations

import importlib.util
import json
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

PATIENT_HEADER_SOURCE_LABEL = "patient_header"
SUPPORTED_FOLLOWUP_VERSIONS = ("1.5", "2.1", "4.0")
FOLLOWUP_SOURCE_LABEL_BY_VERSION = {
    "1.5": "followup_v1_5",
    "2.1": "followup_v2_1",
    "4.0": "followup_v4_0",
}
FOLLOWUP_VERSION_SUFFIX_BY_VERSION = {
    "1.5": "v1_5",
    "2.1": "v2_1",
    "4.0": "v4_0",
}
PATIENT_HEADER_FIELD_NAMES = [
    "vital_status",
    "days_to_last_followup",
    "days_to_last_known_alive",
    "days_to_death",
    "person_neoplasm_cancer_status",
    "new_tumor_event_after_initial_treatment",
    "days_to_new_tumor_event_after_initial_treatment",
]
FOLLOWUP_FIELD_NAMES = [
    "vital_status",
    "days_to_last_followup",
    "days_to_last_known_alive",
    "days_to_death",
    "person_neoplasm_cancer_status",
    "new_tumor_event_after_initial_treatment",
    "days_to_new_tumor_event_after_initial_treatment",
    "new_neoplasm_event_type",
    "new_neoplasm_event_occurrence_anatomic_site",
    "new_neoplasm_occurrence_anatomic_site_text",
    "lost_follow_up",
    "followup_case_report_form_submission_reason",
]
FOLLOWUP_CANDIDATE_FIELD_NAMES = [
    "vital_status",
    "days_to_last_followup",
    "days_to_last_known_alive",
    "days_to_death",
    "person_neoplasm_cancer_status",
    "new_tumor_event_after_initial_treatment",
]
PATIENT_HEADER_EXCLUDED_SUBTREES = {
    "follow_ups",
    "drugs",
    "radiations",
}
MISSING_LIKE_NORMALIZATION = "value.strip().lower()"
MISSING_LIKE_TOKENS = {
    "",
    "[not available]",
    "[not applicable]",
    "[unknown]",
    "[not evaluated]",
    "[discrepancy]",
    "na",
    "n/a",
    "null",
    "none",
    "nan",
}
SUMMARY_FIELDNAMES = [
    "endpoint_target_prep_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
COVERAGE_FIELDNAMES = [
    "endpoint_target_prep_v1_run_id",
    "coverage_section",
    "coverage_metric",
    "dimension_name",
    "dimension_value",
    "coverage_value",
    "notes",
]
PATIENT_ENDPOINT_FIELD_FIELDNAMES = [
    "endpoint_target_prep_v1_run_id",
    "source_run_id",
    "source_file_id",
    "source_filename",
    "source_path",
    "admin_file_uuid",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "source_label",
    *PATIENT_HEADER_FIELD_NAMES,
]
FOLLOWUP_LONG_FIELDNAMES = [
    "endpoint_target_prep_v1_run_id",
    "source_run_id",
    "source_file_id",
    "source_filename",
    "source_path",
    "admin_file_uuid",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "followup_version",
    "followup_source_label",
    "followup_sequence",
    "bcr_followup_barcode",
    "bcr_followup_uuid",
    "form_completion_day",
    "form_completion_month",
    "form_completion_year",
    "form_completion_iso_date",
    "xml_field_name",
    "xml_field_occurrence_index",
    "xml_field_value",
    "is_missing_like",
]
CURRENT_BIOTAB_V4_FIELD_MAP = {
    "bcr_followup_barcode": "clinical_follow_up_v4_0__bcr_followup_barcode_json",
    "bcr_followup_uuid": "clinical_follow_up_v4_0__bcr_followup_uuid_json",
    "lost_follow_up": "clinical_follow_up_v4_0__followup_lost_to_json",
    "days_to_last_followup": "clinical_follow_up_v4_0__last_contact_days_to_json",
    "new_tumor_event_after_initial_treatment": "clinical_follow_up_v4_0__new_tumor_event_dx_indicator_json",
    "person_neoplasm_cancer_status": "clinical_follow_up_v4_0__tumor_status_json",
    "vital_status": "clinical_follow_up_v4_0__vital_status_json",
}
ENDPOINT_FREEZE_STATUS = "blocked_pending_manual_reconciliation"
TREATMENT_FEASIBILITY_STATUS = "unchanged_not_addressed"
OS_STYLE_READINESS_INTERPRETATION = "closer_ready_for_manual_endpoint_freeze_review"


class EndpointTargetPrepV1Error(RuntimeError):
    """Raised when the endpoint-target prep workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the endpoint-target prep workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    source_latest_pointer: Path
    endpoint_crosswalk_latest_pointer: Path
    ambiguity_resolution_latest_pointer: Path
    minimal_cohort_latest_pointer: Path
    baseline_model_input_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved endpoint-target prep inputs loaded from saved audit layers."""

    source_latest_pointer: dict[str, Any]
    source_run_log: dict[str, Any]
    endpoint_crosswalk_latest_pointer: dict[str, Any]
    endpoint_crosswalk_run_log: dict[str, Any]
    ambiguity_resolution_latest_pointer: dict[str, Any]
    ambiguity_resolution_run_log: dict[str, Any]
    minimal_cohort_latest_pointer: dict[str, Any]
    minimal_cohort_run_log: dict[str, Any]
    baseline_model_input_latest_pointer: dict[str, Any]
    baseline_model_input_run_log: dict[str, Any]
    source_metadata_rows: list[dict[str, str]]
    minimal_cohort_rows: list[dict[str, str]]
    baseline_feature_set_audit_map_rows: list[dict[str, str]]
    baseline_model_input_row_count: int
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise EndpointTargetPrepV1Error(f"Required helper script not found: {script_path}")

    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise EndpointTargetPrepV1Error(f"Unable to create an import spec for: {script_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def bool_to_str(value: bool) -> str:
    return "true" if value else "false"


def ordered_unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def json_list(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=True)


def normalize_missing_like(value: str) -> str:
    return value.strip().lower()


def is_missing_like(value: str) -> bool:
    return normalize_missing_like(value) in MISSING_LIKE_TOKENS


def local_name(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", 1)[1]
    return tag


def namespace_uri(tag: str) -> str:
    if tag.startswith("{") and "}" in tag:
        return tag[1:].split("}", 1)[0]
    return ""


def parse_int_or_none(value: str) -> int | None:
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return int(stripped)
    except ValueError:
        return None


def parse_json_list_cell(value: str, label: str) -> list[str]:
    stripped = value.strip()
    if not stripped:
        return []
    try:
        parsed = json.loads(stripped)
    except json.JSONDecodeError as exc:
        raise EndpointTargetPrepV1Error(f"Unable to parse JSON list for {label}: {value!r}") from exc
    if not isinstance(parsed, list):
        raise EndpointTargetPrepV1Error(f"Expected a JSON list for {label}: {value!r}")
    return [str(item).strip() for item in parsed if str(item).strip()]


def normalize_barcode(value: str) -> str:
    return value.strip().upper()


def normalize_uuid(value: str) -> str:
    return value.strip().upper()


def extract_attr_by_local_name(element: ET.Element, attribute_local_name: str) -> str:
    for attribute_name, attribute_value in element.attrib.items():
        if local_name(attribute_name) == attribute_local_name:
            return str(attribute_value).strip()
    return ""


def is_nil_element(element: ET.Element) -> bool:
    return extract_attr_by_local_name(element, "nil").lower() == "true"


def element_text(element: ET.Element | None) -> str:
    if element is None:
        return ""
    if is_nil_element(element):
        return ""
    return (element.text or "").strip()


def iter_descendants_excluding(
    element: ET.Element,
    excluded_local_names: set[str] | None = None,
) -> Any:
    excluded = excluded_local_names or set()
    for child in list(element):
        yield child
        if local_name(child.tag) in excluded:
            continue
        yield from iter_descendants_excluding(child, excluded)


def find_descendant_elements(
    element: ET.Element,
    field_name: str,
    *,
    excluded_local_names: set[str] | None = None,
) -> list[ET.Element]:
    return [
        candidate
        for candidate in iter_descendants_excluding(element, excluded_local_names)
        if local_name(candidate.tag) == field_name
    ]


def first_descendant_text(
    element: ET.Element,
    field_name: str,
    *,
    excluded_local_names: set[str] | None = None,
) -> str:
    matches = find_descendant_elements(element, field_name, excluded_local_names=excluded_local_names)
    if not matches:
        return ""
    return element_text(matches[0])


def all_descendant_texts(
    element: ET.Element,
    field_name: str,
    *,
    excluded_local_names: set[str] | None = None,
) -> list[str]:
    return [
        element_text(candidate)
        for candidate in find_descendant_elements(
            element,
            field_name,
            excluded_local_names=excluded_local_names,
        )
    ]


def derive_iso_date(day_value: str, month_value: str, year_value: str) -> str:
    day_int = parse_int_or_none(day_value)
    month_int = parse_int_or_none(month_value)
    year_int = parse_int_or_none(year_value)
    if day_int is None or month_int is None or year_int is None:
        return ""
    try:
        return date(year_int, month_int, day_int).isoformat()
    except ValueError:
        return ""


def distinct_nonmissing(values: list[str]) -> list[str]:
    return ordered_unique([value for value in values if not is_missing_like(value)])


def canonical_value_set(values: list[str]) -> set[str]:
    return {value.strip() for value in values if not is_missing_like(value)}


def max_int_from_values(values: list[str]) -> int | None:
    parsed = [parse_int_or_none(value) for value in values]
    numeric_values = [value for value in parsed if value is not None]
    if not numeric_values:
        return None
    return max(numeric_values)


def parse_bool_string(value: str) -> bool:
    return value.strip().lower() == "true"


def build_unique_lookup(
    rows: list[dict[str, str]],
    *,
    key_name: str,
    key_builder: Any,
) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        key = key_builder(row)
        if not key:
            raise EndpointTargetPrepV1Error(f"{key_name} produced an empty key while building a lookup.")
        if key in lookup:
            raise EndpointTargetPrepV1Error(f"Duplicate {key_name} detected: {key}")
        lookup[key] = row
    return lookup


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_config = (
        repo_root
        / "09-trials"
        / "01-tcga-only-source-audited"
        / "04-config"
        / "trial_config.yaml"
    )
    if not trial_config.exists():
        raise EndpointTargetPrepV1Error(f"Required trial config not found: {trial_config}")

    trial_config_data = helper_module.load_yaml(trial_config)
    processed_root = repo_root / str(trial_config_data.get("processed_data_root", "01-data/processed"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    endpoint_prep_root = audit_root / "tcga-brca" / "endpoint-prep"

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        processed_runs_root=processed_root / "tcga-brca" / "endpoint-prep" / "xml_followup_v1_runs",
        audit_runs_root=endpoint_prep_root / "xml_followup_v1_runs",
        latest_pointer=endpoint_prep_root / "tcga_brca_endpoint_target_prep_v1_latest.json",
        source_latest_pointer=audit_root / "tcga-brca" / "source" / "tcga_brca_source_supplements_latest.json",
        endpoint_crosswalk_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_endpoint_crosswalk_latest.json"
        ),
        ambiguity_resolution_latest_pointer=(
            audit_root / "tcga-brca" / "cohort" / "tcga_brca_blueprint_ambiguity_resolution_latest.json"
        ),
        minimal_cohort_latest_pointer=(
            audit_root / "tcga-brca" / "cohort" / "tcga_brca_minimal_cohort_v1_latest.json"
        ),
        baseline_model_input_latest_pointer=(
            audit_root / "tcga-brca" / "model-input" / "tcga_brca_baseline_model_input_v1_latest.json"
        ),
    )


def require_pointer(
    pointer_path: Path,
    *,
    label: str,
    required_keys: set[str],
    helper_module: Any,
) -> dict[str, Any]:
    if not pointer_path.exists():
        raise EndpointTargetPrepV1Error(f"Required {label} latest pointer not found: {pointer_path}")
    payload = helper_module.load_json(pointer_path)
    helper_module.require_keys(payload, required_keys, label, pointer_path)
    return payload


def require_completed_run_log(
    *,
    repo_root: Path,
    run_log_relative_path: str,
    label: str,
    helper_module: Any,
) -> tuple[Path, dict[str, Any]]:
    run_log_path = helper_module.resolve_existing_path(repo_root, run_log_relative_path, label)
    run_log = helper_module.load_json(run_log_path)
    if run_log.get("status") != "completed":
        raise EndpointTargetPrepV1Error(f"{label} is not completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise EndpointTargetPrepV1Error(f"{label} does not report validation.passed == true.")
    return run_log_path, run_log


def require_columns(rows: list[dict[str, str]], required_columns: set[str], label: str) -> None:
    if not rows:
        raise EndpointTargetPrepV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise EndpointTargetPrepV1Error(f"{label} is missing required columns: {sorted(missing)}")


def select_bcr_xml_source_rows(source_metadata_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for row in source_metadata_rows:
        if str(row.get("source_class") or "") != "clinical":
            continue
        if str(row.get("data_format") or "") != "BCR XML":
            continue
        selected.append(row)

    selected.sort(
        key=lambda row: (
            parse_json_list_cell(str(row.get("case_submitter_ids") or ""), "case_submitter_ids")[0],
            str(row.get("file_name") or ""),
        )
    )
    return selected


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    source_latest_pointer = require_pointer(
        paths.source_latest_pointer,
        label="source supplements",
        required_keys={"run_id", "metadata_tsv", "run_log_json", "source_classes"},
        helper_module=helper_module,
    )
    clinical_source_class = source_latest_pointer.get("source_classes", {}).get("clinical")
    if not isinstance(clinical_source_class, dict):
        raise EndpointTargetPrepV1Error("The source supplements latest pointer has no clinical source class.")
    if not bool(clinical_source_class.get("download_completed", False)):
        raise EndpointTargetPrepV1Error("The clinical source class does not report download_completed == true.")

    endpoint_crosswalk_latest_pointer = require_pointer(
        paths.endpoint_crosswalk_latest_pointer,
        label="endpoint crosswalk",
        required_keys={
            "crosswalk_run_id",
            "endpoint_crosswalk_tsv",
            "endpoint_crosswalk_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )
    ambiguity_resolution_latest_pointer = require_pointer(
        paths.ambiguity_resolution_latest_pointer,
        label="blueprint ambiguity resolution",
        required_keys={
            "ambiguity_resolution_run_id",
            "ambiguity_resolution_actions_tsv",
            "ambiguity_resolution_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )
    minimal_cohort_latest_pointer = require_pointer(
        paths.minimal_cohort_latest_pointer,
        label="minimal cohort v1",
        required_keys={
            "cohort_v1_build_id",
            "minimal_cohort_v1_tsv",
            "minimal_cohort_v1_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )
    baseline_model_input_latest_pointer = require_pointer(
        paths.baseline_model_input_latest_pointer,
        label="baseline model-input v1",
        required_keys={
            "baseline_model_input_v1_run_id",
            "baseline_model_input_v1_tsv",
            "baseline_model_input_v1_summary_tsv",
            "baseline_feature_set_v1_audit_map_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )

    source_run_log_path, source_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(source_latest_pointer["run_log_json"]),
        label="source supplements run log",
        helper_module=helper_module,
    )
    endpoint_crosswalk_run_log_path, endpoint_crosswalk_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(endpoint_crosswalk_latest_pointer["run_log_json"]),
        label="endpoint crosswalk run log",
        helper_module=helper_module,
    )
    ambiguity_resolution_run_log_path, ambiguity_resolution_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(ambiguity_resolution_latest_pointer["run_log_json"]),
        label="ambiguity resolution run log",
        helper_module=helper_module,
    )
    minimal_cohort_run_log_path, minimal_cohort_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(minimal_cohort_latest_pointer["run_log_json"]),
        label="minimal cohort v1 run log",
        helper_module=helper_module,
    )
    baseline_model_input_run_log_path, baseline_model_input_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(baseline_model_input_latest_pointer["run_log_json"]),
        label="baseline model-input v1 run log",
        helper_module=helper_module,
    )

    input_paths = {
        "source_metadata_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(source_latest_pointer["metadata_tsv"]),
            "source supplements metadata TSV",
        ),
        "clinical_download_dir": helper_module.resolve_existing_path(
            paths.repo_root,
            str(clinical_source_class["download_dir"]),
            "clinical download directory",
        ),
        "source_run_log_json": source_run_log_path,
        "endpoint_crosswalk_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_crosswalk_tsv"]),
            "endpoint crosswalk TSV",
        ),
        "endpoint_crosswalk_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_crosswalk_summary_tsv"]),
            "endpoint crosswalk summary TSV",
        ),
        "endpoint_crosswalk_run_log_json": endpoint_crosswalk_run_log_path,
        "ambiguity_resolution_actions_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["ambiguity_resolution_actions_tsv"]),
            "ambiguity resolution actions TSV",
        ),
        "ambiguity_resolution_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["ambiguity_resolution_summary_tsv"]),
            "ambiguity resolution summary TSV",
        ),
        "ambiguity_resolution_run_log_json": ambiguity_resolution_run_log_path,
        "minimal_cohort_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["minimal_cohort_v1_tsv"]),
            "minimal cohort v1 TSV",
        ),
        "minimal_cohort_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["minimal_cohort_v1_summary_tsv"]),
            "minimal cohort v1 summary TSV",
        ),
        "minimal_cohort_v1_run_log_json": minimal_cohort_run_log_path,
        "baseline_model_input_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_model_input_latest_pointer["baseline_model_input_v1_tsv"]),
            "baseline model-input v1 TSV",
        ),
        "baseline_model_input_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_model_input_latest_pointer["baseline_model_input_v1_summary_tsv"]),
            "baseline model-input v1 summary TSV",
        ),
        "baseline_feature_set_v1_audit_map_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_model_input_latest_pointer["baseline_feature_set_v1_audit_map_tsv"]),
            "baseline feature-set v1 audit map TSV",
        ),
        "baseline_model_input_v1_run_log_json": baseline_model_input_run_log_path,
    }

    source_metadata_rows = helper_module.read_tsv_dict_rows(input_paths["source_metadata_tsv"])
    minimal_cohort_rows = helper_module.read_tsv_dict_rows(input_paths["minimal_cohort_v1_tsv"])
    baseline_feature_set_audit_map_rows = helper_module.read_tsv_dict_rows(
        input_paths["baseline_feature_set_v1_audit_map_tsv"]
    )
    baseline_model_input_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_model_input_v1_tsv"])

    require_columns(
        source_metadata_rows,
        {"source_class", "file_id", "file_name", "data_format", "case_ids", "case_submitter_ids"},
        "source supplements metadata TSV",
    )
    require_columns(
        minimal_cohort_rows,
        {
            "cohort_v1_build_id",
            "provisional_patient_row_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "followup_match_row_count",
            *CURRENT_BIOTAB_V4_FIELD_MAP.values(),
            "clinical_patient__vital_status",
            "clinical_patient__last_contact_days_to",
            "clinical_patient__tumor_status",
            "clinical_patient__new_tumor_event_dx_indicator",
        },
        "minimal_cohort_v1.tsv",
    )
    require_columns(
        baseline_feature_set_audit_map_rows,
        {
            "feature_set_v1_row_index",
            "baseline_analysis_v1_row_id",
            "provisional_patient_row_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
        },
        "baseline_feature_set_v1_audit_map.tsv",
    )

    return WorkflowInputs(
        source_latest_pointer=source_latest_pointer,
        source_run_log=source_run_log,
        endpoint_crosswalk_latest_pointer=endpoint_crosswalk_latest_pointer,
        endpoint_crosswalk_run_log=endpoint_crosswalk_run_log,
        ambiguity_resolution_latest_pointer=ambiguity_resolution_latest_pointer,
        ambiguity_resolution_run_log=ambiguity_resolution_run_log,
        minimal_cohort_latest_pointer=minimal_cohort_latest_pointer,
        minimal_cohort_run_log=minimal_cohort_run_log,
        baseline_model_input_latest_pointer=baseline_model_input_latest_pointer,
        baseline_model_input_run_log=baseline_model_input_run_log,
        source_metadata_rows=source_metadata_rows,
        minimal_cohort_rows=minimal_cohort_rows,
        baseline_feature_set_audit_map_rows=baseline_feature_set_audit_map_rows,
        baseline_model_input_row_count=len(baseline_model_input_rows),
        input_paths=input_paths,
    )


def validate_upstream_state(workflow_inputs: WorkflowInputs) -> None:
    bcr_xml_source_rows = select_bcr_xml_source_rows(workflow_inputs.source_metadata_rows)
    if not bcr_xml_source_rows:
        raise EndpointTargetPrepV1Error("No BCR XML clinical source rows were found in the source metadata TSV.")
    if len(workflow_inputs.minimal_cohort_rows) != len(workflow_inputs.baseline_feature_set_audit_map_rows):
        raise EndpointTargetPrepV1Error(
            "Minimal cohort v1 rows and baseline feature-set audit map rows do not match in count."
        )
    if len(workflow_inputs.minimal_cohort_rows) != workflow_inputs.baseline_model_input_row_count:
        raise EndpointTargetPrepV1Error(
            "Minimal cohort v1 rows and baseline model-input v1 rows do not match in count."
        )


def infer_followup_version(followup_element: ET.Element) -> str:
    version = extract_attr_by_local_name(followup_element, "version")
    if version:
        return version
    namespace = namespace_uri(followup_element.tag)
    if namespace:
        candidate = namespace.rstrip("/").split("/")[-1]
        if candidate:
            return candidate
    return ""


def parse_clinical_xml_patient(
    *,
    endpoint_target_prep_v1_run_id: str,
    source_run_id: str,
    repo_root: Path,
    source_row: dict[str, str],
    clinical_download_dir: Path,
) -> tuple[dict[str, str], list[dict[str, Any]], set[str]]:
    file_id = str(source_row.get("file_id") or "").strip()
    file_name = str(source_row.get("file_name") or "").strip()
    if not file_id or not file_name:
        raise EndpointTargetPrepV1Error(f"Source metadata row is missing file_id or file_name: {source_row}")

    xml_path = clinical_download_dir / file_id / file_name
    if not xml_path.exists():
        raise EndpointTargetPrepV1Error(f"Expected XML file was not found on disk: {xml_path}")

    metadata_barcodes = parse_json_list_cell(
        str(source_row.get("case_submitter_ids") or ""),
        f"case_submitter_ids for {file_name}",
    )
    metadata_case_ids = parse_json_list_cell(
        str(source_row.get("case_ids") or ""),
        f"case_ids for {file_name}",
    )
    metadata_barcode = normalize_barcode(metadata_barcodes[0]) if metadata_barcodes else ""
    metadata_patient_uuid = normalize_uuid(metadata_case_ids[0]) if metadata_case_ids else ""

    root = ET.parse(xml_path).getroot()
    patient_nodes = [child for child in list(root) if local_name(child.tag) == "patient"]
    if len(patient_nodes) != 1:
        raise EndpointTargetPrepV1Error(
            f"Expected exactly one patient element in XML file {xml_path}, found {len(patient_nodes)}."
        )
    patient_node = patient_nodes[0]

    admin_file_uuid = normalize_uuid(first_descendant_text(root, "file_uuid"))
    xml_barcode = normalize_barcode(
        first_descendant_text(
            patient_node,
            "bcr_patient_barcode",
            excluded_local_names=PATIENT_HEADER_EXCLUDED_SUBTREES,
        )
    )
    xml_patient_uuid = normalize_uuid(
        first_descendant_text(
            patient_node,
            "bcr_patient_uuid",
            excluded_local_names=PATIENT_HEADER_EXCLUDED_SUBTREES,
        )
    )
    patient_barcode = xml_barcode or metadata_barcode
    patient_uuid = xml_patient_uuid or metadata_patient_uuid

    if metadata_barcode and patient_barcode and metadata_barcode != patient_barcode:
        raise EndpointTargetPrepV1Error(
            f"Metadata/XML barcode mismatch for {xml_path}: {metadata_barcode} != {patient_barcode}"
        )
    if metadata_patient_uuid and patient_uuid and metadata_patient_uuid != patient_uuid:
        raise EndpointTargetPrepV1Error(
            f"Metadata/XML UUID mismatch for {xml_path}: {metadata_patient_uuid} != {patient_uuid}"
        )
    if not patient_barcode:
        raise EndpointTargetPrepV1Error(f"Unable to resolve bcr_patient_barcode for {xml_path}")

    patient_header_row = {
        "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
        "source_run_id": source_run_id,
        "source_file_id": file_id,
        "source_filename": file_name,
        "source_path": repo_relative(xml_path, repo_root),
        "admin_file_uuid": admin_file_uuid,
        "bcr_patient_barcode": patient_barcode,
        "bcr_patient_uuid": patient_uuid,
        "source_label": PATIENT_HEADER_SOURCE_LABEL,
    }
    for field_name in PATIENT_HEADER_FIELD_NAMES:
        patient_header_row[field_name] = first_descendant_text(
            patient_node,
            field_name,
            excluded_local_names=PATIENT_HEADER_EXCLUDED_SUBTREES,
        )

    followup_records: list[dict[str, Any]] = []
    unsupported_versions: set[str] = set()
    followup_container_nodes = [child for child in list(patient_node) if local_name(child.tag) == "follow_ups"]
    for container_node in followup_container_nodes:
        followup_nodes = [candidate for candidate in container_node.iter() if local_name(candidate.tag) == "follow_up"]
        for followup_node in followup_nodes:
            version = infer_followup_version(followup_node)
            if version not in SUPPORTED_FOLLOWUP_VERSIONS:
                unsupported_versions.add(version or "[missing]")
                continue

            followup_sequence = extract_attr_by_local_name(followup_node, "sequence")
            bcr_followup_barcode = normalize_barcode(first_descendant_text(followup_node, "bcr_followup_barcode"))
            bcr_followup_uuid = normalize_uuid(first_descendant_text(followup_node, "bcr_followup_uuid"))
            form_completion_day = first_descendant_text(followup_node, "day_of_form_completion")
            form_completion_month = first_descendant_text(followup_node, "month_of_form_completion")
            form_completion_year = first_descendant_text(followup_node, "year_of_form_completion")
            form_completion_iso_date = derive_iso_date(
                form_completion_day,
                form_completion_month,
                form_completion_year,
            )

            field_values: dict[str, list[str]] = {}
            for field_name in FOLLOWUP_FIELD_NAMES:
                values = all_descendant_texts(followup_node, field_name)
                field_values[field_name] = values if values else [""]

            followup_records.append(
                {
                    "followup_version": version,
                    "followup_source_label": FOLLOWUP_SOURCE_LABEL_BY_VERSION[version],
                    "followup_sequence": followup_sequence,
                    "bcr_followup_barcode": bcr_followup_barcode,
                    "bcr_followup_uuid": bcr_followup_uuid,
                    "form_completion_day": form_completion_day,
                    "form_completion_month": form_completion_month,
                    "form_completion_year": form_completion_year,
                    "form_completion_iso_date": form_completion_iso_date,
                    "field_values": field_values,
                }
            )

    return patient_header_row, followup_records, unsupported_versions


def build_followup_long_rows(
    *,
    patient_header_row: dict[str, str],
    followup_records: list[dict[str, Any]],
) -> list[dict[str, str]]:
    long_rows: list[dict[str, str]] = []
    for record in followup_records:
        for field_name in FOLLOWUP_FIELD_NAMES:
            values = list(record["field_values"].get(field_name, [""]))
            if not values:
                values = [""]
            for occurrence_index, value in enumerate(values, start=1):
                long_rows.append(
                    {
                        "endpoint_target_prep_v1_run_id": patient_header_row["endpoint_target_prep_v1_run_id"],
                        "source_run_id": patient_header_row["source_run_id"],
                        "source_file_id": patient_header_row["source_file_id"],
                        "source_filename": patient_header_row["source_filename"],
                        "source_path": patient_header_row["source_path"],
                        "admin_file_uuid": patient_header_row["admin_file_uuid"],
                        "bcr_patient_barcode": patient_header_row["bcr_patient_barcode"],
                        "bcr_patient_uuid": patient_header_row["bcr_patient_uuid"],
                        "followup_version": str(record["followup_version"]),
                        "followup_source_label": str(record["followup_source_label"]),
                        "followup_sequence": str(record["followup_sequence"]),
                        "bcr_followup_barcode": str(record["bcr_followup_barcode"]),
                        "bcr_followup_uuid": str(record["bcr_followup_uuid"]),
                        "form_completion_day": str(record["form_completion_day"]),
                        "form_completion_month": str(record["form_completion_month"]),
                        "form_completion_year": str(record["form_completion_year"]),
                        "form_completion_iso_date": str(record["form_completion_iso_date"]),
                        "xml_field_name": field_name,
                        "xml_field_occurrence_index": str(occurrence_index),
                        "xml_field_value": value,
                        "is_missing_like": bool_to_str(is_missing_like(value)),
                    }
                )
    return long_rows


def candidate_objects_for_field(
    followup_records: list[dict[str, Any]],
    *,
    version: str,
    field_name: str,
) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for record in followup_records:
        if str(record["followup_version"]) != version:
            continue
        values = list(record["field_values"].get(field_name, [""]))
        for value in values:
            if is_missing_like(value):
                continue
            candidates.append(
                {
                    "followup_sequence": str(record["followup_sequence"]),
                    "bcr_followup_barcode": str(record["bcr_followup_barcode"]),
                    "bcr_followup_uuid": str(record["bcr_followup_uuid"]),
                    "form_completion_iso_date": str(record["form_completion_iso_date"]),
                    "value": value,
                }
            )
    return candidates


def nonmissing_values_from_followups(
    followup_records: list[dict[str, Any]],
    *,
    field_name: str,
    version: str | None = None,
) -> list[str]:
    values: list[str] = []
    for record in followup_records:
        if version is not None and str(record["followup_version"]) != version:
            continue
        values.extend(list(record["field_values"].get(field_name, [])))
    return distinct_nonmissing(values)


def nonmissing_values_any_xml(
    patient_header_row: dict[str, str],
    followup_records: list[dict[str, Any]],
    *,
    field_name: str,
) -> list[str]:
    values: list[str] = [str(patient_header_row.get(field_name) or "")]
    values.extend(nonmissing_values_from_followups(followup_records, field_name=field_name))
    return distinct_nonmissing(values)


def build_status_flags(
    *,
    has_any_xml_followup: bool,
    has_only_older_xml_followup: bool,
    has_xml_gain_over_biotab_v4: bool,
    has_nonmissing_xml_days_to_death: bool,
    has_nonmissing_xml_last_contact_like: bool,
    has_os_style_ingredients_any_xml: bool,
    xml_followup_max_last_contact_exceeds_header: bool,
    xml_header_vs_followup_vital_status_disagreement: bool,
    xml_header_vs_followup_last_contact_disagreement: bool,
    xml_header_vs_followup_days_to_death_disagreement: bool,
    xml_followup_multi_version_overlap: bool,
    xml_followup_multi_vital_status_values: bool,
    xml_v4_vs_current_biotab_v4_vital_status_disagreement: bool,
    xml_v4_vs_current_biotab_v4_last_contact_disagreement: bool,
) -> list[str]:
    flags: list[str] = []
    if has_any_xml_followup:
        flags.append("has_any_xml_followup")
    if has_only_older_xml_followup:
        flags.append("has_only_older_xml_followup")
    if has_xml_gain_over_biotab_v4:
        flags.append("has_xml_gain_over_current_biotab_v4")
    if has_nonmissing_xml_days_to_death:
        flags.append("has_nonmissing_xml_days_to_death")
    if has_nonmissing_xml_last_contact_like:
        flags.append("has_nonmissing_xml_last_contact_like")
    if has_os_style_ingredients_any_xml:
        flags.append("has_os_style_ingredients_any_xml")
    if xml_followup_max_last_contact_exceeds_header:
        flags.append("xml_followup_max_last_contact_exceeds_header")
    if xml_header_vs_followup_vital_status_disagreement:
        flags.append("xml_header_vs_followup_vital_status_disagreement")
    if xml_header_vs_followup_last_contact_disagreement:
        flags.append("xml_header_vs_followup_last_contact_disagreement")
    if xml_header_vs_followup_days_to_death_disagreement:
        flags.append("xml_header_vs_followup_days_to_death_disagreement")
    if xml_followup_multi_version_overlap:
        flags.append("xml_followup_multi_version_overlap")
    if xml_followup_multi_vital_status_values:
        flags.append("xml_followup_multi_vital_status_values")
    if xml_v4_vs_current_biotab_v4_vital_status_disagreement:
        flags.append("xml_v4_vs_current_biotab_v4_vital_status_disagreement")
    if xml_v4_vs_current_biotab_v4_last_contact_disagreement:
        flags.append("xml_v4_vs_current_biotab_v4_last_contact_disagreement")
    if not has_os_style_ingredients_any_xml:
        flags.append("os_style_ingredients_still_incomplete")
    return flags


def build_endpoint_prep_fieldnames() -> list[str]:
    fieldnames = [
        "endpoint_target_prep_v1_run_id",
        "cohort_v1_build_id",
        "baseline_model_input_v1_run_id",
        "bcr_patient_barcode",
        "bcr_patient_uuid",
        "provisional_patient_row_id",
        "baseline_analysis_v1_row_id",
        "feature_set_v1_row_index",
        "patient_header_source_label",
        "patient_header_vital_status",
        "patient_header_days_to_last_followup",
        "patient_header_days_to_last_known_alive",
        "patient_header_days_to_death",
        "patient_header_person_neoplasm_cancer_status",
        "patient_header_new_tumor_event_after_initial_treatment",
        "patient_header_days_to_new_tumor_event_after_initial_treatment",
        "xml_followup_record_count_total",
        "xml_followup_record_count_v1_5",
        "xml_followup_record_count_v2_1",
        "xml_followup_record_count_v4_0",
        "xml_followup_versions_present_json",
    ]
    for version in SUPPORTED_FOLLOWUP_VERSIONS:
        version_suffix = FOLLOWUP_VERSION_SUFFIX_BY_VERSION[version]
        for field_name in FOLLOWUP_CANDIDATE_FIELD_NAMES:
            fieldnames.append(f"followup_{version_suffix}_{field_name}_candidates_json")
    fieldnames.extend(
        [
            "has_any_xml_followup",
            "has_xml_followup_v1_5",
            "has_xml_followup_v2_1",
            "has_xml_followup_v4_0",
            "has_only_older_xml_followup",
            "has_overlapping_xml_followup_versions",
            "has_nonmissing_xml_days_to_death",
            "has_nonmissing_xml_last_contact_like",
            "has_os_style_ingredients_any_xml",
            "has_xml_gain_over_biotab_v4",
            "xml_followup_multi_vital_status_values",
            "xml_header_vs_followup_vital_status_disagreement",
            "xml_header_vs_followup_last_contact_disagreement",
            "xml_header_vs_followup_days_to_death_disagreement",
            "xml_followup_max_days_to_last_followup",
            "xml_followup_max_days_to_last_known_alive",
            "xml_followup_max_days_to_death",
            "xml_followup_max_last_contact_exceeds_header",
            "xml_followup_max_last_contact_equals_header",
            "xml_followup_max_last_contact_less_than_header",
            "provisional_endpoint_prep_status_flags_json",
        ]
    )
    return fieldnames


def build_overlap_audit_fieldnames() -> list[str]:
    return [
        "endpoint_target_prep_v1_run_id",
        "cohort_v1_build_id",
        "bcr_patient_barcode",
        "bcr_patient_uuid",
        "provisional_patient_row_id",
        "baseline_analysis_v1_row_id",
        "feature_set_v1_row_index",
        "current_biotab_v4_followup_match_row_count",
        "current_biotab_v4_has_match",
        "current_biotab_v4_bcr_followup_barcode_json",
        "current_biotab_v4_bcr_followup_uuid_json",
        "current_biotab_v4_lost_follow_up_json",
        "current_biotab_v4_days_to_last_followup_json",
        "current_biotab_v4_new_tumor_event_after_initial_treatment_json",
        "current_biotab_v4_person_neoplasm_cancer_status_json",
        "current_biotab_v4_vital_status_json",
        "xml_patient_header_vital_status",
        "xml_patient_header_days_to_last_followup",
        "xml_patient_header_days_to_last_known_alive",
        "xml_patient_header_days_to_death",
        "xml_patient_header_person_neoplasm_cancer_status",
        "xml_patient_header_new_tumor_event_after_initial_treatment",
        "xml_followup_versions_present_json",
        "xml_followup_record_count_total",
        "xml_followup_record_count_v1_5",
        "xml_followup_record_count_v2_1",
        "xml_followup_record_count_v4_0",
        "xml_followup_vital_status_values_json",
        "xml_followup_days_to_last_followup_values_json",
        "xml_followup_days_to_last_known_alive_values_json",
        "xml_followup_days_to_death_values_json",
        "xml_followup_person_neoplasm_cancer_status_values_json",
        "xml_followup_new_tumor_event_after_initial_treatment_values_json",
        "xml_v4_followup_vital_status_values_json",
        "xml_v4_followup_days_to_last_followup_values_json",
        "xml_v4_followup_person_neoplasm_cancer_status_values_json",
        "xml_v4_followup_new_tumor_event_after_initial_treatment_values_json",
        "xml_followup_max_days_to_last_followup",
        "xml_followup_max_days_to_last_known_alive",
        "xml_followup_max_days_to_death",
        "xml_followup_max_last_contact_exceeds_header",
        "xml_followup_max_last_contact_equals_header",
        "xml_followup_max_last_contact_less_than_header",
        "xml_header_vs_followup_vital_status_disagreement",
        "xml_header_vs_followup_last_contact_disagreement",
        "xml_header_vs_followup_days_to_death_disagreement",
        "xml_v4_vs_current_biotab_v4_vital_status_disagreement",
        "xml_v4_vs_current_biotab_v4_last_contact_disagreement",
        "xml_v4_vs_current_biotab_v4_person_neoplasm_cancer_status_disagreement",
        "xml_v4_vs_current_biotab_v4_new_tumor_event_after_initial_treatment_disagreement",
        "has_any_xml_followup",
        "has_only_older_xml_followup",
        "has_overlapping_xml_followup_versions",
        "has_nonmissing_xml_days_to_death",
        "has_nonmissing_xml_last_contact_like",
        "has_os_style_ingredients_any_xml",
        "has_xml_gain_over_biotab_v4",
        "overlap_audit_flags_json",
    ]


def build_patient_level_outputs(
    *,
    endpoint_target_prep_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
    parsed_patients_by_barcode: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    minimal_cohort_lookup = build_unique_lookup(
        workflow_inputs.minimal_cohort_rows,
        key_name="minimal cohort provisional_patient_row_id",
        key_builder=lambda row: str(row["provisional_patient_row_id"]).strip(),
    )
    bridge_lookup = build_unique_lookup(
        workflow_inputs.baseline_feature_set_audit_map_rows,
        key_name="baseline feature-set audit map provisional_patient_row_id",
        key_builder=lambda row: str(row["provisional_patient_row_id"]).strip(),
    )

    endpoint_prep_rows: list[dict[str, str]] = []
    overlap_audit_rows: list[dict[str, str]] = []

    version_patient_counter: Counter[str] = Counter()
    version_record_counter: Counter[str] = Counter()
    combo_counter: Counter[str] = Counter()
    unsupported_current_biotab_death_layer = True
    xml_followup_patients = 0
    patients_with_current_biotab_v4 = 0
    patients_with_xml_followup_and_no_current_biotab_v4 = 0
    patients_with_xml_gain_over_biotab_v4 = 0
    patients_with_later_followup_than_header = 0
    patients_with_any_xml_days_to_death = 0
    patients_with_any_xml_last_contact_like = 0
    patients_with_os_style_ingredients = 0

    endpoint_prep_fieldnames = build_endpoint_prep_fieldnames()
    overlap_audit_fieldnames = build_overlap_audit_fieldnames()

    for provisional_patient_row_id in sorted(
        minimal_cohort_lookup.keys(),
        key=lambda value: int(value),
    ):
        minimal_row = minimal_cohort_lookup[provisional_patient_row_id]
        bridge_row = bridge_lookup.get(provisional_patient_row_id)
        if bridge_row is None:
            raise EndpointTargetPrepV1Error(
                f"Missing baseline feature-set audit-map row for provisional_patient_row_id={provisional_patient_row_id}"
            )

        minimal_barcode = normalize_barcode(str(minimal_row["bcr_patient_barcode"]))
        minimal_uuid = normalize_uuid(str(minimal_row["bcr_patient_uuid"]))
        if normalize_barcode(str(bridge_row["bcr_patient_barcode"])) != minimal_barcode:
            raise EndpointTargetPrepV1Error(
                f"Barcode mismatch between minimal cohort and audit map for provisional_patient_row_id={provisional_patient_row_id}"
            )
        if normalize_uuid(str(bridge_row["bcr_patient_uuid"])) != minimal_uuid:
            raise EndpointTargetPrepV1Error(
                f"UUID mismatch between minimal cohort and audit map for provisional_patient_row_id={provisional_patient_row_id}"
            )

        parsed_patient = parsed_patients_by_barcode.get(minimal_barcode)
        if parsed_patient is None:
            raise EndpointTargetPrepV1Error(
                f"Minimal cohort patient {minimal_barcode} was not found in the parsed XML patient set."
            )

        patient_header_row = parsed_patient["patient_header_row"]
        followup_records = list(parsed_patient["followup_records"])
        version_counts = Counter(str(record["followup_version"]) for record in followup_records)
        version_list = [version for version in SUPPORTED_FOLLOWUP_VERSIONS if version_counts.get(version, 0) > 0]
        combo_label = "|".join(version_list) if version_list else "none"

        if followup_records:
            xml_followup_patients += 1
        for version in SUPPORTED_FOLLOWUP_VERSIONS:
            if version_counts.get(version, 0) > 0:
                version_patient_counter[version] += 1
            version_record_counter[version] += version_counts.get(version, 0)
        combo_counter[combo_label] += 1

        current_biotab_v4_match_row_count = parse_int_or_none(str(minimal_row["followup_match_row_count"])) or 0
        current_biotab_v4_has_match = current_biotab_v4_match_row_count > 0
        if current_biotab_v4_has_match:
            patients_with_current_biotab_v4 += 1

        current_biotab_values = {
            field_name: parse_json_list_cell(str(minimal_row[column_name] or ""), column_name)
            for field_name, column_name in CURRENT_BIOTAB_V4_FIELD_MAP.items()
        }
        current_biotab_vital_values = distinct_nonmissing(current_biotab_values["vital_status"])
        current_biotab_last_contact_values = distinct_nonmissing(current_biotab_values["days_to_last_followup"])
        current_biotab_tumor_status_values = distinct_nonmissing(
            current_biotab_values["person_neoplasm_cancer_status"]
        )
        current_biotab_new_tumor_values = distinct_nonmissing(
            current_biotab_values["new_tumor_event_after_initial_treatment"]
        )
        current_biotab_max_last_contact = max_int_from_values(current_biotab_last_contact_values)

        followup_vital_values = nonmissing_values_from_followups(
            followup_records,
            field_name="vital_status",
        )
        followup_days_to_last_followup_values = nonmissing_values_from_followups(
            followup_records,
            field_name="days_to_last_followup",
        )
        followup_days_to_last_known_alive_values = nonmissing_values_from_followups(
            followup_records,
            field_name="days_to_last_known_alive",
        )
        followup_days_to_death_values = nonmissing_values_from_followups(
            followup_records,
            field_name="days_to_death",
        )
        followup_tumor_status_values = nonmissing_values_from_followups(
            followup_records,
            field_name="person_neoplasm_cancer_status",
        )
        followup_new_tumor_values = nonmissing_values_from_followups(
            followup_records,
            field_name="new_tumor_event_after_initial_treatment",
        )

        xml_v4_vital_values = nonmissing_values_from_followups(
            followup_records,
            field_name="vital_status",
            version="4.0",
        )
        xml_v4_last_contact_values = nonmissing_values_from_followups(
            followup_records,
            field_name="days_to_last_followup",
            version="4.0",
        )
        xml_v4_tumor_status_values = nonmissing_values_from_followups(
            followup_records,
            field_name="person_neoplasm_cancer_status",
            version="4.0",
        )
        xml_v4_new_tumor_values = nonmissing_values_from_followups(
            followup_records,
            field_name="new_tumor_event_after_initial_treatment",
            version="4.0",
        )

        patient_header_vital_status = str(patient_header_row.get("vital_status") or "")
        patient_header_days_to_last_followup = str(patient_header_row.get("days_to_last_followup") or "")
        patient_header_days_to_last_known_alive = str(
            patient_header_row.get("days_to_last_known_alive") or ""
        )
        patient_header_days_to_death = str(patient_header_row.get("days_to_death") or "")
        patient_header_tumor_status = str(patient_header_row.get("person_neoplasm_cancer_status") or "")
        patient_header_new_tumor = str(
            patient_header_row.get("new_tumor_event_after_initial_treatment") or ""
        )

        all_xml_days_to_death_values = nonmissing_values_any_xml(
            patient_header_row,
            followup_records,
            field_name="days_to_death",
        )
        all_xml_last_contact_like_values = distinct_nonmissing(
            [
                patient_header_days_to_last_followup,
                patient_header_days_to_last_known_alive,
                *followup_days_to_last_followup_values,
                *followup_days_to_last_known_alive_values,
            ]
        )
        all_xml_vital_status_values = distinct_nonmissing(
            [patient_header_vital_status, *followup_vital_values]
        )

        has_any_xml_followup = bool(followup_records)
        has_xml_followup_v1_5 = version_counts.get("1.5", 0) > 0
        has_xml_followup_v2_1 = version_counts.get("2.1", 0) > 0
        has_xml_followup_v4_0 = version_counts.get("4.0", 0) > 0
        has_only_older_xml_followup = has_any_xml_followup and not has_xml_followup_v4_0
        has_overlapping_xml_followup_versions = len(version_list) > 1
        has_nonmissing_xml_days_to_death = bool(all_xml_days_to_death_values)
        has_nonmissing_xml_last_contact_like = bool(all_xml_last_contact_like_values)
        has_os_style_ingredients_any_xml = bool(all_xml_vital_status_values) and (
            has_nonmissing_xml_days_to_death or has_nonmissing_xml_last_contact_like
        )

        if has_nonmissing_xml_days_to_death:
            patients_with_any_xml_days_to_death += 1
        if has_nonmissing_xml_last_contact_like:
            patients_with_any_xml_last_contact_like += 1
        if has_os_style_ingredients_any_xml:
            patients_with_os_style_ingredients += 1

        xml_followup_max_days_to_last_followup = max_int_from_values(followup_days_to_last_followup_values)
        xml_followup_max_days_to_last_known_alive = max_int_from_values(
            followup_days_to_last_known_alive_values
        )
        xml_followup_max_days_to_death = max_int_from_values(followup_days_to_death_values)
        header_days_to_last_followup_int = parse_int_or_none(patient_header_days_to_last_followup)

        xml_followup_max_last_contact_exceeds_header = (
            xml_followup_max_days_to_last_followup is not None
            and header_days_to_last_followup_int is not None
            and xml_followup_max_days_to_last_followup > header_days_to_last_followup_int
        )
        xml_followup_max_last_contact_equals_header = (
            xml_followup_max_days_to_last_followup is not None
            and header_days_to_last_followup_int is not None
            and xml_followup_max_days_to_last_followup == header_days_to_last_followup_int
        )
        xml_followup_max_last_contact_less_than_header = (
            xml_followup_max_days_to_last_followup is not None
            and header_days_to_last_followup_int is not None
            and xml_followup_max_days_to_last_followup < header_days_to_last_followup_int
        )
        if xml_followup_max_last_contact_exceeds_header:
            patients_with_later_followup_than_header += 1

        xml_followup_multi_vital_status_values = len(canonical_value_set(followup_vital_values)) > 1
        xml_header_vs_followup_vital_status_disagreement = (
            bool(canonical_value_set([patient_header_vital_status]))
            and bool(canonical_value_set(followup_vital_values))
            and canonical_value_set([patient_header_vital_status]) != canonical_value_set(followup_vital_values)
        )
        xml_header_vs_followup_last_contact_disagreement = (
            bool(canonical_value_set([patient_header_days_to_last_followup]))
            and bool(canonical_value_set(followup_days_to_last_followup_values))
            and canonical_value_set([patient_header_days_to_last_followup])
            != canonical_value_set(followup_days_to_last_followup_values)
        )
        xml_header_vs_followup_days_to_death_disagreement = (
            bool(canonical_value_set([patient_header_days_to_death]))
            and bool(canonical_value_set(followup_days_to_death_values))
            and canonical_value_set([patient_header_days_to_death])
            != canonical_value_set(followup_days_to_death_values)
        )

        xml_v4_vs_current_biotab_v4_vital_status_disagreement = (
            bool(canonical_value_set(xml_v4_vital_values))
            and bool(canonical_value_set(current_biotab_vital_values))
            and canonical_value_set(xml_v4_vital_values) != canonical_value_set(current_biotab_vital_values)
        )
        xml_v4_vs_current_biotab_v4_last_contact_disagreement = (
            bool(canonical_value_set(xml_v4_last_contact_values))
            and bool(canonical_value_set(current_biotab_last_contact_values))
            and canonical_value_set(xml_v4_last_contact_values)
            != canonical_value_set(current_biotab_last_contact_values)
        )
        xml_v4_vs_current_biotab_v4_tumor_status_disagreement = (
            bool(canonical_value_set(xml_v4_tumor_status_values))
            and bool(canonical_value_set(current_biotab_tumor_status_values))
            and canonical_value_set(xml_v4_tumor_status_values)
            != canonical_value_set(current_biotab_tumor_status_values)
        )
        xml_v4_vs_current_biotab_v4_new_tumor_event_disagreement = (
            bool(canonical_value_set(xml_v4_new_tumor_values))
            and bool(canonical_value_set(current_biotab_new_tumor_values))
            and canonical_value_set(xml_v4_new_tumor_values)
            != canonical_value_set(current_biotab_new_tumor_values)
        )

        has_xml_gain_over_biotab_v4 = has_any_xml_followup and (
            (not current_biotab_v4_has_match)
            or (
                xml_followup_max_days_to_last_followup is not None
                and (
                    current_biotab_max_last_contact is None
                    or xml_followup_max_days_to_last_followup > current_biotab_max_last_contact
                )
            )
            or has_nonmissing_xml_days_to_death
        )
        if has_any_xml_followup and not current_biotab_v4_has_match:
            patients_with_xml_followup_and_no_current_biotab_v4 += 1
        if has_xml_gain_over_biotab_v4:
            patients_with_xml_gain_over_biotab_v4 += 1

        status_flags = build_status_flags(
            has_any_xml_followup=has_any_xml_followup,
            has_only_older_xml_followup=has_only_older_xml_followup,
            has_xml_gain_over_biotab_v4=has_xml_gain_over_biotab_v4,
            has_nonmissing_xml_days_to_death=has_nonmissing_xml_days_to_death,
            has_nonmissing_xml_last_contact_like=has_nonmissing_xml_last_contact_like,
            has_os_style_ingredients_any_xml=has_os_style_ingredients_any_xml,
            xml_followup_max_last_contact_exceeds_header=xml_followup_max_last_contact_exceeds_header,
            xml_header_vs_followup_vital_status_disagreement=xml_header_vs_followup_vital_status_disagreement,
            xml_header_vs_followup_last_contact_disagreement=xml_header_vs_followup_last_contact_disagreement,
            xml_header_vs_followup_days_to_death_disagreement=xml_header_vs_followup_days_to_death_disagreement,
            xml_followup_multi_version_overlap=has_overlapping_xml_followup_versions,
            xml_followup_multi_vital_status_values=xml_followup_multi_vital_status_values,
            xml_v4_vs_current_biotab_v4_vital_status_disagreement=(
                xml_v4_vs_current_biotab_v4_vital_status_disagreement
            ),
            xml_v4_vs_current_biotab_v4_last_contact_disagreement=(
                xml_v4_vs_current_biotab_v4_last_contact_disagreement
            ),
        )

        endpoint_prep_row = {fieldname: "" for fieldname in endpoint_prep_fieldnames}
        endpoint_prep_row.update(
            {
                "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
                "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
                "baseline_model_input_v1_run_id": str(
                    workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
                ),
                "bcr_patient_barcode": minimal_barcode,
                "bcr_patient_uuid": minimal_uuid,
                "provisional_patient_row_id": provisional_patient_row_id,
                "baseline_analysis_v1_row_id": str(bridge_row["baseline_analysis_v1_row_id"]),
                "feature_set_v1_row_index": str(bridge_row["feature_set_v1_row_index"]),
                "patient_header_source_label": PATIENT_HEADER_SOURCE_LABEL,
                "patient_header_vital_status": patient_header_vital_status,
                "patient_header_days_to_last_followup": patient_header_days_to_last_followup,
                "patient_header_days_to_last_known_alive": patient_header_days_to_last_known_alive,
                "patient_header_days_to_death": patient_header_days_to_death,
                "patient_header_person_neoplasm_cancer_status": patient_header_tumor_status,
                "patient_header_new_tumor_event_after_initial_treatment": patient_header_new_tumor,
                "patient_header_days_to_new_tumor_event_after_initial_treatment": str(
                    patient_header_row.get("days_to_new_tumor_event_after_initial_treatment") or ""
                ),
                "xml_followup_record_count_total": str(len(followup_records)),
                "xml_followup_record_count_v1_5": str(version_counts.get("1.5", 0)),
                "xml_followup_record_count_v2_1": str(version_counts.get("2.1", 0)),
                "xml_followup_record_count_v4_0": str(version_counts.get("4.0", 0)),
                "xml_followup_versions_present_json": json_list(version_list),
                "has_any_xml_followup": bool_to_str(has_any_xml_followup),
                "has_xml_followup_v1_5": bool_to_str(has_xml_followup_v1_5),
                "has_xml_followup_v2_1": bool_to_str(has_xml_followup_v2_1),
                "has_xml_followup_v4_0": bool_to_str(has_xml_followup_v4_0),
                "has_only_older_xml_followup": bool_to_str(has_only_older_xml_followup),
                "has_overlapping_xml_followup_versions": bool_to_str(has_overlapping_xml_followup_versions),
                "has_nonmissing_xml_days_to_death": bool_to_str(has_nonmissing_xml_days_to_death),
                "has_nonmissing_xml_last_contact_like": bool_to_str(has_nonmissing_xml_last_contact_like),
                "has_os_style_ingredients_any_xml": bool_to_str(has_os_style_ingredients_any_xml),
                "has_xml_gain_over_biotab_v4": bool_to_str(has_xml_gain_over_biotab_v4),
                "xml_followup_multi_vital_status_values": bool_to_str(xml_followup_multi_vital_status_values),
                "xml_header_vs_followup_vital_status_disagreement": bool_to_str(
                    xml_header_vs_followup_vital_status_disagreement
                ),
                "xml_header_vs_followup_last_contact_disagreement": bool_to_str(
                    xml_header_vs_followup_last_contact_disagreement
                ),
                "xml_header_vs_followup_days_to_death_disagreement": bool_to_str(
                    xml_header_vs_followup_days_to_death_disagreement
                ),
                "xml_followup_max_days_to_last_followup": (
                    str(xml_followup_max_days_to_last_followup)
                    if xml_followup_max_days_to_last_followup is not None
                    else ""
                ),
                "xml_followup_max_days_to_last_known_alive": (
                    str(xml_followup_max_days_to_last_known_alive)
                    if xml_followup_max_days_to_last_known_alive is not None
                    else ""
                ),
                "xml_followup_max_days_to_death": (
                    str(xml_followup_max_days_to_death)
                    if xml_followup_max_days_to_death is not None
                    else ""
                ),
                "xml_followup_max_last_contact_exceeds_header": bool_to_str(
                    xml_followup_max_last_contact_exceeds_header
                ),
                "xml_followup_max_last_contact_equals_header": bool_to_str(
                    xml_followup_max_last_contact_equals_header
                ),
                "xml_followup_max_last_contact_less_than_header": bool_to_str(
                    xml_followup_max_last_contact_less_than_header
                ),
                "provisional_endpoint_prep_status_flags_json": json_list(status_flags),
            }
        )

        for version in SUPPORTED_FOLLOWUP_VERSIONS:
            version_suffix = FOLLOWUP_VERSION_SUFFIX_BY_VERSION[version]
            for field_name in FOLLOWUP_CANDIDATE_FIELD_NAMES:
                endpoint_prep_row[f"followup_{version_suffix}_{field_name}_candidates_json"] = json.dumps(
                    candidate_objects_for_field(
                        followup_records,
                        version=version,
                        field_name=field_name,
                    ),
                    ensure_ascii=True,
                )

        overlap_audit_row = {fieldname: "" for fieldname in overlap_audit_fieldnames}
        overlap_audit_row.update(
            {
                "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
                "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
                "bcr_patient_barcode": minimal_barcode,
                "bcr_patient_uuid": minimal_uuid,
                "provisional_patient_row_id": provisional_patient_row_id,
                "baseline_analysis_v1_row_id": str(bridge_row["baseline_analysis_v1_row_id"]),
                "feature_set_v1_row_index": str(bridge_row["feature_set_v1_row_index"]),
                "current_biotab_v4_followup_match_row_count": str(current_biotab_v4_match_row_count),
                "current_biotab_v4_has_match": bool_to_str(current_biotab_v4_has_match),
                "current_biotab_v4_bcr_followup_barcode_json": str(
                    minimal_row[CURRENT_BIOTAB_V4_FIELD_MAP["bcr_followup_barcode"]]
                ),
                "current_biotab_v4_bcr_followup_uuid_json": str(
                    minimal_row[CURRENT_BIOTAB_V4_FIELD_MAP["bcr_followup_uuid"]]
                ),
                "current_biotab_v4_lost_follow_up_json": str(
                    minimal_row[CURRENT_BIOTAB_V4_FIELD_MAP["lost_follow_up"]]
                ),
                "current_biotab_v4_days_to_last_followup_json": str(
                    minimal_row[CURRENT_BIOTAB_V4_FIELD_MAP["days_to_last_followup"]]
                ),
                "current_biotab_v4_new_tumor_event_after_initial_treatment_json": str(
                    minimal_row[CURRENT_BIOTAB_V4_FIELD_MAP["new_tumor_event_after_initial_treatment"]]
                ),
                "current_biotab_v4_person_neoplasm_cancer_status_json": str(
                    minimal_row[CURRENT_BIOTAB_V4_FIELD_MAP["person_neoplasm_cancer_status"]]
                ),
                "current_biotab_v4_vital_status_json": str(
                    minimal_row[CURRENT_BIOTAB_V4_FIELD_MAP["vital_status"]]
                ),
                "xml_patient_header_vital_status": patient_header_vital_status,
                "xml_patient_header_days_to_last_followup": patient_header_days_to_last_followup,
                "xml_patient_header_days_to_last_known_alive": patient_header_days_to_last_known_alive,
                "xml_patient_header_days_to_death": patient_header_days_to_death,
                "xml_patient_header_person_neoplasm_cancer_status": patient_header_tumor_status,
                "xml_patient_header_new_tumor_event_after_initial_treatment": patient_header_new_tumor,
                "xml_followup_versions_present_json": json_list(version_list),
                "xml_followup_record_count_total": str(len(followup_records)),
                "xml_followup_record_count_v1_5": str(version_counts.get("1.5", 0)),
                "xml_followup_record_count_v2_1": str(version_counts.get("2.1", 0)),
                "xml_followup_record_count_v4_0": str(version_counts.get("4.0", 0)),
                "xml_followup_vital_status_values_json": json_list(followup_vital_values),
                "xml_followup_days_to_last_followup_values_json": json_list(
                    followup_days_to_last_followup_values
                ),
                "xml_followup_days_to_last_known_alive_values_json": json_list(
                    followup_days_to_last_known_alive_values
                ),
                "xml_followup_days_to_death_values_json": json_list(followup_days_to_death_values),
                "xml_followup_person_neoplasm_cancer_status_values_json": json_list(
                    followup_tumor_status_values
                ),
                "xml_followup_new_tumor_event_after_initial_treatment_values_json": json_list(
                    followup_new_tumor_values
                ),
                "xml_v4_followup_vital_status_values_json": json_list(xml_v4_vital_values),
                "xml_v4_followup_days_to_last_followup_values_json": json_list(xml_v4_last_contact_values),
                "xml_v4_followup_person_neoplasm_cancer_status_values_json": json_list(
                    xml_v4_tumor_status_values
                ),
                "xml_v4_followup_new_tumor_event_after_initial_treatment_values_json": json_list(
                    xml_v4_new_tumor_values
                ),
                "xml_followup_max_days_to_last_followup": (
                    str(xml_followup_max_days_to_last_followup)
                    if xml_followup_max_days_to_last_followup is not None
                    else ""
                ),
                "xml_followup_max_days_to_last_known_alive": (
                    str(xml_followup_max_days_to_last_known_alive)
                    if xml_followup_max_days_to_last_known_alive is not None
                    else ""
                ),
                "xml_followup_max_days_to_death": (
                    str(xml_followup_max_days_to_death)
                    if xml_followup_max_days_to_death is not None
                    else ""
                ),
                "xml_followup_max_last_contact_exceeds_header": bool_to_str(
                    xml_followup_max_last_contact_exceeds_header
                ),
                "xml_followup_max_last_contact_equals_header": bool_to_str(
                    xml_followup_max_last_contact_equals_header
                ),
                "xml_followup_max_last_contact_less_than_header": bool_to_str(
                    xml_followup_max_last_contact_less_than_header
                ),
                "xml_header_vs_followup_vital_status_disagreement": bool_to_str(
                    xml_header_vs_followup_vital_status_disagreement
                ),
                "xml_header_vs_followup_last_contact_disagreement": bool_to_str(
                    xml_header_vs_followup_last_contact_disagreement
                ),
                "xml_header_vs_followup_days_to_death_disagreement": bool_to_str(
                    xml_header_vs_followup_days_to_death_disagreement
                ),
                "xml_v4_vs_current_biotab_v4_vital_status_disagreement": bool_to_str(
                    xml_v4_vs_current_biotab_v4_vital_status_disagreement
                ),
                "xml_v4_vs_current_biotab_v4_last_contact_disagreement": bool_to_str(
                    xml_v4_vs_current_biotab_v4_last_contact_disagreement
                ),
                "xml_v4_vs_current_biotab_v4_person_neoplasm_cancer_status_disagreement": bool_to_str(
                    xml_v4_vs_current_biotab_v4_tumor_status_disagreement
                ),
                "xml_v4_vs_current_biotab_v4_new_tumor_event_after_initial_treatment_disagreement": bool_to_str(
                    xml_v4_vs_current_biotab_v4_new_tumor_event_disagreement
                ),
                "has_any_xml_followup": bool_to_str(has_any_xml_followup),
                "has_only_older_xml_followup": bool_to_str(has_only_older_xml_followup),
                "has_overlapping_xml_followup_versions": bool_to_str(has_overlapping_xml_followup_versions),
                "has_nonmissing_xml_days_to_death": bool_to_str(has_nonmissing_xml_days_to_death),
                "has_nonmissing_xml_last_contact_like": bool_to_str(has_nonmissing_xml_last_contact_like),
                "has_os_style_ingredients_any_xml": bool_to_str(has_os_style_ingredients_any_xml),
                "has_xml_gain_over_biotab_v4": bool_to_str(has_xml_gain_over_biotab_v4),
                "overlap_audit_flags_json": json_list(status_flags),
            }
        )

        endpoint_prep_rows.append(endpoint_prep_row)
        overlap_audit_rows.append(overlap_audit_row)

    metrics = {
        "version_patient_counter": version_patient_counter,
        "version_record_counter": version_record_counter,
        "combo_counter": combo_counter,
        "xml_followup_patients": xml_followup_patients,
        "patients_with_current_biotab_v4": patients_with_current_biotab_v4,
        "patients_with_xml_followup_and_no_current_biotab_v4": (
            patients_with_xml_followup_and_no_current_biotab_v4
        ),
        "patients_with_xml_gain_over_biotab_v4": patients_with_xml_gain_over_biotab_v4,
        "patients_with_later_followup_than_header": patients_with_later_followup_than_header,
        "patients_with_any_xml_days_to_death": patients_with_any_xml_days_to_death,
        "patients_with_any_xml_last_contact_like": patients_with_any_xml_last_contact_like,
        "patients_with_os_style_ingredients": patients_with_os_style_ingredients,
        "unsupported_current_biotab_death_layer": unsupported_current_biotab_death_layer,
    }
    return endpoint_prep_rows, overlap_audit_rows, metrics


def build_followup_version_coverage_rows(
    *,
    endpoint_target_prep_v1_run_id: str,
    endpoint_prep_rows: list[dict[str, str]],
    overlap_audit_rows: list[dict[str, str]],
    patient_header_rows: list[dict[str, str]],
    patient_metrics: dict[str, Any],
) -> list[dict[str, str]]:
    coverage_rows: list[dict[str, str]] = []

    def add_row(
        coverage_section: str,
        coverage_metric: str,
        coverage_value: Any,
        *,
        dimension_name: str = "",
        dimension_value: str = "",
        notes: str = "",
    ) -> None:
        coverage_rows.append(
            {
                "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
                "coverage_section": coverage_section,
                "coverage_metric": coverage_metric,
                "dimension_name": dimension_name,
                "dimension_value": dimension_value,
                "coverage_value": str(coverage_value),
                "notes": notes,
            }
        )

    add_row(
        "xml_patient_files",
        "parsed_xml_patient_files",
        len(patient_header_rows),
        dimension_name="scope",
        dimension_value="all",
    )

    for version in SUPPORTED_FOLLOWUP_VERSIONS:
        add_row(
            "followup_version_patient_counts",
            "patients_with_followup_version",
            patient_metrics["version_patient_counter"].get(version, 0),
            dimension_name="followup_version",
            dimension_value=version,
        )
        add_row(
            "followup_version_record_counts",
            "followup_record_count",
            patient_metrics["version_record_counter"].get(version, 0),
            dimension_name="followup_version",
            dimension_value=version,
        )

    for combo_label, combo_count in sorted(patient_metrics["combo_counter"].items()):
        add_row(
            "followup_version_combo_patient_counts",
            "patients_with_followup_version_combo",
            combo_count,
            dimension_name="followup_version_combo",
            dimension_value=combo_label,
        )

    add_row(
        "older_version_only_counts",
        "patients_with_v1_5_only",
        sum(parse_bool_string(row["has_xml_followup_v1_5"]) and not (
            parse_bool_string(row["has_xml_followup_v2_1"]) or parse_bool_string(row["has_xml_followup_v4_0"])
        ) for row in endpoint_prep_rows),
        dimension_name="subset",
        dimension_value="v1_5_only",
    )
    add_row(
        "older_version_only_counts",
        "patients_with_v2_1_only",
        sum(parse_bool_string(row["has_xml_followup_v2_1"]) and not (
            parse_bool_string(row["has_xml_followup_v1_5"]) or parse_bool_string(row["has_xml_followup_v4_0"])
        ) for row in endpoint_prep_rows),
        dimension_name="subset",
        dimension_value="v2_1_only",
    )
    add_row(
        "older_version_only_counts",
        "patients_with_v4_0_only",
        sum(parse_bool_string(row["has_xml_followup_v4_0"]) and not (
            parse_bool_string(row["has_xml_followup_v1_5"]) or parse_bool_string(row["has_xml_followup_v2_1"])
        ) for row in endpoint_prep_rows),
        dimension_name="subset",
        dimension_value="v4_0_only",
    )
    add_row(
        "older_version_only_counts",
        "patients_with_only_older_xml_followup",
        sum(parse_bool_string(row["has_only_older_xml_followup"]) for row in endpoint_prep_rows),
        dimension_name="subset",
        dimension_value="older_only",
    )
    add_row(
        "older_version_only_counts",
        "patients_with_overlapping_versions",
        sum(parse_bool_string(row["has_overlapping_xml_followup_versions"]) for row in endpoint_prep_rows),
        dimension_name="subset",
        dimension_value="overlap",
    )

    add_row(
        "gain_vs_current_biotab_v4",
        "patients_with_current_biotab_v4_match",
        patient_metrics["patients_with_current_biotab_v4"],
        dimension_name="comparison_scope",
        dimension_value="current_biotab_v4_lineage",
    )
    add_row(
        "gain_vs_current_biotab_v4",
        "patients_without_current_biotab_v4_match",
        len(overlap_audit_rows) - patient_metrics["patients_with_current_biotab_v4"],
        dimension_name="comparison_scope",
        dimension_value="current_biotab_v4_lineage",
    )
    add_row(
        "gain_vs_current_biotab_v4",
        "patients_with_xml_followup_and_no_current_biotab_v4",
        patient_metrics["patients_with_xml_followup_and_no_current_biotab_v4"],
        dimension_name="comparison_scope",
        dimension_value="coverage_gain",
    )
    add_row(
        "gain_vs_current_biotab_v4",
        "patients_with_xml_gain_over_biotab_v4",
        patient_metrics["patients_with_xml_gain_over_biotab_v4"],
        dimension_name="comparison_scope",
        dimension_value="coverage_or_signal_gain",
    )
    add_row(
        "gain_vs_current_biotab_v4",
        "patients_with_xml_v4_vs_current_biotab_v4_last_contact_disagreement",
        sum(
            parse_bool_string(row["xml_v4_vs_current_biotab_v4_last_contact_disagreement"])
            for row in overlap_audit_rows
        ),
        dimension_name="comparison_scope",
        dimension_value="parsing_consistency_check",
    )
    add_row(
        "gain_vs_current_biotab_v4",
        "patients_with_xml_v4_vs_current_biotab_v4_vital_status_disagreement",
        sum(
            parse_bool_string(row["xml_v4_vs_current_biotab_v4_vital_status_disagreement"])
            for row in overlap_audit_rows
        ),
        dimension_name="comparison_scope",
        dimension_value="parsing_consistency_check",
    )

    return coverage_rows


def build_summary_rows(
    *,
    endpoint_target_prep_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
    patient_header_rows: list[dict[str, str]],
    followup_long_rows: list[dict[str, str]],
    endpoint_prep_rows: list[dict[str, str]],
    overlap_audit_rows: list[dict[str, str]],
    patient_metrics: dict[str, Any],
) -> list[dict[str, str]]:
    summary_rows: list[dict[str, str]] = []

    def add_row(summary_section: str, summary_metric: str, summary_value: Any, notes: str = "") -> None:
        summary_rows.append(
            {
                "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
                "summary_section": summary_section,
                "summary_metric": summary_metric,
                "summary_value": str(summary_value),
                "notes": notes,
            }
        )

    patient_header_lookup = build_unique_lookup(
        patient_header_rows,
        key_name="patient header barcode",
        key_builder=lambda row: str(row["bcr_patient_barcode"]).strip(),
    )

    add_row("input_runs", "source_run_id", workflow_inputs.source_latest_pointer["run_id"])
    add_row(
        "input_runs",
        "endpoint_crosswalk_run_id",
        workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"],
    )
    add_row(
        "input_runs",
        "ambiguity_resolution_run_id",
        workflow_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"],
    )
    add_row("input_runs", "cohort_v1_build_id", workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"])
    add_row(
        "input_runs",
        "baseline_model_input_v1_run_id",
        workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"],
    )

    add_row("row_counts", "total_xml_patient_files_parsed", len(patient_header_rows))
    add_row("row_counts", "clinical_xml_patient_endpoint_fields_row_count", len(patient_header_rows))
    add_row("row_counts", "clinical_xml_followup_fields_long_row_count", len(followup_long_rows))
    add_row("row_counts", "endpoint_target_prep_v1_row_count", len(endpoint_prep_rows))
    add_row("row_counts", "endpoint_target_prep_v1_overlap_audit_row_count", len(overlap_audit_rows))

    for field_name in PATIENT_HEADER_FIELD_NAMES:
        add_row(
            "patient_header_coverage",
            f"patient_header_nonmissing_{field_name}_count",
            sum(not is_missing_like(str(row[field_name])) for row in patient_header_rows),
        )

    add_row("followup_coverage", "patients_with_any_xml_followup", patient_metrics["xml_followup_patients"])
    for version in SUPPORTED_FOLLOWUP_VERSIONS:
        version_suffix = FOLLOWUP_VERSION_SUFFIX_BY_VERSION[version]
        add_row(
            "followup_coverage",
            f"patients_with_followup_{version_suffix}",
            patient_metrics["version_patient_counter"].get(version, 0),
        )
        add_row(
            "followup_coverage",
            f"followup_record_count_{version_suffix}",
            patient_metrics["version_record_counter"].get(version, 0),
        )
    add_row(
        "followup_coverage",
        "patients_with_only_older_xml_followup",
        sum(parse_bool_string(row["has_only_older_xml_followup"]) for row in endpoint_prep_rows),
    )
    add_row(
        "followup_coverage",
        "patients_with_overlapping_xml_followup_versions",
        sum(parse_bool_string(row["has_overlapping_xml_followup_versions"]) for row in endpoint_prep_rows),
    )
    add_row(
        "followup_coverage",
        "patients_with_no_xml_followup",
        len(endpoint_prep_rows) - patient_metrics["xml_followup_patients"],
    )

    add_row(
        "gain_vs_current_biotab_v4",
        "patients_gaining_followup_beyond_current_biotab_v4",
        patient_metrics["patients_with_xml_followup_and_no_current_biotab_v4"],
        notes="Defined as any XML follow-up present with no current grouped biotab v4 match in minimal_cohort_v1.",
    )
    add_row(
        "gain_vs_current_biotab_v4",
        "patients_with_xml_gain_over_biotab_v4",
        patient_metrics["patients_with_xml_gain_over_biotab_v4"],
        notes="Union flag covering older-version-only coverage, later XML last-contact evidence, or XML days_to_death evidence.",
    )
    add_row(
        "days_to_death",
        "patients_with_nonmissing_days_to_death_any_xml",
        patient_metrics["patients_with_any_xml_days_to_death"],
        notes="Counts patient-header or follow-up XML days_to_death evidence.",
    )
    add_row(
        "last_contact_review",
        "patients_with_later_followup_than_patient_header_days_to_last_followup",
        patient_metrics["patients_with_later_followup_than_header"],
        notes="Compares follow-up max days_to_last_followup against patient-header days_to_last_followup only.",
    )
    add_row(
        "os_style_endpoint_prep",
        "patients_with_nonmissing_last_contact_like_any_xml",
        patient_metrics["patients_with_any_xml_last_contact_like"],
        notes="days_to_last_followup is primary comparable field; days_to_last_known_alive is carried as auxiliary evidence only.",
    )
    add_row(
        "os_style_endpoint_prep",
        "patients_with_usable_os_style_ingredients_any_xml",
        patient_metrics["patients_with_os_style_ingredients"],
        notes="Defined as any XML vital_status plus any XML days_to_death or last-contact-like evidence.",
    )
    add_row(
        "os_style_endpoint_prep",
        "os_style_endpoint_prep_readiness_interpretation",
        OS_STYLE_READINESS_INTERPRETATION,
        notes="Endpoint freeze remains blocked; this workflow prepares a patient-level reconciliation layer only.",
    )
    add_row("status", "endpoint_freeze_status", ENDPOINT_FREEZE_STATUS)
    add_row("status", "treatment_feasibility_status", TREATMENT_FEASIBILITY_STATUS)

    sample_barcodes = ["TCGA-AO-A0J5", "TCGA-AR-A0TT", "TCGA-BH-A0HF"]
    for sample_barcode in sample_barcodes:
        patient_row = patient_header_lookup.get(sample_barcode)
        if patient_row is not None:
            add_row(
                "manual_validation_examples",
                f"sample_patient_present_{sample_barcode}",
                "true",
                notes=patient_row["source_path"],
            )

    return summary_rows


def validate_generated_outputs(
    *,
    workflow_inputs: WorkflowInputs,
    bcr_xml_source_rows: list[dict[str, str]],
    patient_header_rows: list[dict[str, str]],
    followup_long_rows: list[dict[str, str]],
    endpoint_prep_rows: list[dict[str, str]],
    overlap_audit_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    unsupported_followup_versions: set[str],
) -> dict[str, Any]:
    summary_lookup = {row["summary_metric"]: row for row in summary_rows}
    output_rows_positive = bool(
        patient_header_rows
        and followup_long_rows
        and endpoint_prep_rows
        and overlap_audit_rows
        and summary_rows
    )
    xml_files_found = len(bcr_xml_source_rows) > 0
    xml_patient_file_count_matches_minimal_cohort = len(bcr_xml_source_rows) == len(
        workflow_inputs.minimal_cohort_rows
    )
    patient_header_row_count_matches_xml_file_count = len(patient_header_rows) == len(bcr_xml_source_rows)
    endpoint_prep_row_count_matches_minimal_cohort = len(endpoint_prep_rows) == len(
        workflow_inputs.minimal_cohort_rows
    )
    overlap_audit_row_count_matches_endpoint_prep = len(overlap_audit_rows) == len(endpoint_prep_rows)
    row_bridge_fully_matched_to_baseline_feature_set_v1_audit_map = len(endpoint_prep_rows) == len(
        workflow_inputs.baseline_feature_set_audit_map_rows
    )
    followup_versions_restricted_to_supported_set = not unsupported_followup_versions
    summary_endpoint_freeze_blocked = (
        summary_lookup["endpoint_freeze_status"]["summary_value"] == ENDPOINT_FREEZE_STATUS
    )
    summary_treatment_unchanged = (
        summary_lookup["treatment_feasibility_status"]["summary_value"] == TREATMENT_FEASIBILITY_STATUS
    )

    return {
        "passed": (
            xml_files_found
            and xml_patient_file_count_matches_minimal_cohort
            and patient_header_row_count_matches_xml_file_count
            and output_rows_positive
            and endpoint_prep_row_count_matches_minimal_cohort
            and overlap_audit_row_count_matches_endpoint_prep
            and row_bridge_fully_matched_to_baseline_feature_set_v1_audit_map
            and followup_versions_restricted_to_supported_set
            and summary_endpoint_freeze_blocked
            and summary_treatment_unchanged
        ),
        "required_upstream_pointers_found": True,
        "source_latest_pointer_found": True,
        "source_run_log_completed": True,
        "source_validation_passed": True,
        "endpoint_crosswalk_latest_pointer_found": True,
        "endpoint_crosswalk_run_log_completed": True,
        "endpoint_crosswalk_validation_passed": True,
        "ambiguity_resolution_latest_pointer_found": True,
        "ambiguity_resolution_run_log_completed": True,
        "ambiguity_resolution_validation_passed": True,
        "minimal_cohort_v1_latest_pointer_found": True,
        "minimal_cohort_v1_run_log_completed": True,
        "minimal_cohort_v1_validation_passed": True,
        "baseline_model_input_v1_latest_pointer_found": True,
        "baseline_model_input_v1_run_log_completed": True,
        "baseline_model_input_v1_validation_passed": True,
        "xml_files_found": xml_files_found,
        "xml_patient_file_count_matches_minimal_cohort": xml_patient_file_count_matches_minimal_cohort,
        "patient_header_row_count_matches_xml_file_count": patient_header_row_count_matches_xml_file_count,
        "output_rows_positive": output_rows_positive,
        "endpoint_prep_row_count_matches_minimal_cohort": endpoint_prep_row_count_matches_minimal_cohort,
        "overlap_audit_row_count_matches_endpoint_prep": overlap_audit_row_count_matches_endpoint_prep,
        "row_bridge_fully_matched_to_baseline_feature_set_v1_audit_map": (
            row_bridge_fully_matched_to_baseline_feature_set_v1_audit_map
        ),
        "followup_versions_restricted_to_supported_set": followup_versions_restricted_to_supported_set,
        "unsupported_followup_versions_json": json_list(sorted(unsupported_followup_versions)),
        "summary_endpoint_freeze_blocked": summary_endpoint_freeze_blocked,
        "summary_treatment_unchanged": summary_treatment_unchanged,
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
    }


def build_latest_pointer_payload(
    *,
    endpoint_target_prep_v1_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
        "source_run_id": str(workflow_inputs.source_latest_pointer["run_id"]),
        "endpoint_crosswalk_run_id": str(workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]),
        "ambiguity_resolution_run_id": str(
            workflow_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
        "baseline_model_input_v1_run_id": str(
            workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
        ),
        "baseline_analysis_v1_run_id": str(
            workflow_inputs.baseline_model_input_run_log["baseline_analysis_v1_run_id"]
        ),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "clinical_xml_patient_endpoint_fields_tsv": repo_relative(
            output_paths["clinical_xml_patient_endpoint_fields_tsv"],
            paths.repo_root,
        ),
        "clinical_xml_followup_fields_long_tsv": repo_relative(
            output_paths["clinical_xml_followup_fields_long_tsv"],
            paths.repo_root,
        ),
        "clinical_xml_followup_version_coverage_tsv": repo_relative(
            output_paths["clinical_xml_followup_version_coverage_tsv"],
            paths.repo_root,
        ),
        "endpoint_target_prep_v1_tsv": repo_relative(
            output_paths["endpoint_target_prep_v1_tsv"],
            paths.repo_root,
        ),
        "endpoint_target_prep_v1_overlap_audit_tsv": repo_relative(
            output_paths["endpoint_target_prep_v1_overlap_audit_tsv"],
            paths.repo_root,
        ),
        "endpoint_target_prep_v1_summary_tsv": repo_relative(
            output_paths["endpoint_target_prep_v1_summary_tsv"],
            paths.repo_root,
        ),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "source_supplements_latest_json": repo_relative(paths.source_latest_pointer, paths.repo_root),
        "endpoint_crosswalk_latest_json": repo_relative(
            paths.endpoint_crosswalk_latest_pointer,
            paths.repo_root,
        ),
        "ambiguity_resolution_latest_json": repo_relative(
            paths.ambiguity_resolution_latest_pointer,
            paths.repo_root,
        ),
        "minimal_cohort_v1_latest_json": repo_relative(paths.minimal_cohort_latest_pointer, paths.repo_root),
        "baseline_model_input_v1_latest_json": repo_relative(
            paths.baseline_model_input_latest_pointer,
            paths.repo_root,
        ),
        "baseline_feature_set_v1_audit_map_tsv": str(
            workflow_inputs.baseline_model_input_latest_pointer["baseline_feature_set_v1_audit_map_tsv"]
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def build_upstream_snapshots(workflow_inputs: WorkflowInputs) -> dict[str, Any]:
    return {
        "source_supplements_latest_pointer": workflow_inputs.source_latest_pointer,
        "source_run_log_status": {
            "status": workflow_inputs.source_run_log.get("status"),
            "run_id": workflow_inputs.source_run_log.get("run_id"),
            "validation": workflow_inputs.source_run_log.get("validation"),
        },
        "endpoint_crosswalk_latest_pointer": workflow_inputs.endpoint_crosswalk_latest_pointer,
        "endpoint_crosswalk_run_log_status": {
            "status": workflow_inputs.endpoint_crosswalk_run_log.get("status"),
            "crosswalk_run_id": workflow_inputs.endpoint_crosswalk_run_log.get("crosswalk_run_id"),
            "validation": workflow_inputs.endpoint_crosswalk_run_log.get("validation"),
        },
        "ambiguity_resolution_latest_pointer": workflow_inputs.ambiguity_resolution_latest_pointer,
        "ambiguity_resolution_run_log_status": {
            "status": workflow_inputs.ambiguity_resolution_run_log.get("status"),
            "ambiguity_resolution_run_id": workflow_inputs.ambiguity_resolution_run_log.get(
                "ambiguity_resolution_run_id"
            ),
            "validation": workflow_inputs.ambiguity_resolution_run_log.get("validation"),
        },
        "minimal_cohort_v1_latest_pointer": workflow_inputs.minimal_cohort_latest_pointer,
        "minimal_cohort_v1_run_log_status": {
            "status": workflow_inputs.minimal_cohort_run_log.get("status"),
            "cohort_v1_build_id": workflow_inputs.minimal_cohort_run_log.get("cohort_v1_build_id"),
            "validation": workflow_inputs.minimal_cohort_run_log.get("validation"),
        },
        "baseline_model_input_v1_latest_pointer": workflow_inputs.baseline_model_input_latest_pointer,
        "baseline_model_input_v1_run_log_status": {
            "status": workflow_inputs.baseline_model_input_run_log.get("status"),
            "baseline_model_input_v1_run_id": workflow_inputs.baseline_model_input_run_log.get(
                "baseline_model_input_v1_run_id"
            ),
            "validation": workflow_inputs.baseline_model_input_run_log.get("validation"),
        },
    }


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    endpoint_target_prep_v1_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    processed_run_dir = paths.processed_runs_root / endpoint_target_prep_v1_run_id
    audit_run_dir = paths.audit_runs_root / endpoint_target_prep_v1_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)
        validate_upstream_state(workflow_inputs)

        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        patient_fields_path = processed_run_dir / "clinical_xml_patient_endpoint_fields.tsv"
        followup_long_path = processed_run_dir / "clinical_xml_followup_fields_long.tsv"
        coverage_path = audit_run_dir / "clinical_xml_followup_version_coverage.tsv"
        endpoint_prep_path = processed_run_dir / "endpoint_target_prep_v1.tsv"
        overlap_audit_path = audit_run_dir / "endpoint_target_prep_v1_overlap_audit.tsv"
        summary_path = audit_run_dir / "endpoint_target_prep_v1_summary.tsv"

        clinical_download_dir = workflow_inputs.input_paths["clinical_download_dir"]
        bcr_xml_source_rows = select_bcr_xml_source_rows(workflow_inputs.source_metadata_rows)
        parsed_patients_by_barcode: dict[str, dict[str, Any]] = {}
        patient_header_rows: list[dict[str, str]] = []
        followup_long_rows: list[dict[str, str]] = []
        unsupported_followup_versions: set[str] = set()

        for source_row in bcr_xml_source_rows:
            patient_header_row, followup_records, patient_unsupported_versions = parse_clinical_xml_patient(
                endpoint_target_prep_v1_run_id=endpoint_target_prep_v1_run_id,
                source_run_id=str(workflow_inputs.source_latest_pointer["run_id"]),
                repo_root=paths.repo_root,
                source_row=source_row,
                clinical_download_dir=clinical_download_dir,
            )
            barcode = patient_header_row["bcr_patient_barcode"]
            if barcode in parsed_patients_by_barcode:
                raise EndpointTargetPrepV1Error(f"Duplicate parsed XML patient barcode detected: {barcode}")
            parsed_patients_by_barcode[barcode] = {
                "patient_header_row": patient_header_row,
                "followup_records": followup_records,
            }
            patient_header_rows.append(patient_header_row)
            followup_long_rows.extend(
                build_followup_long_rows(
                    patient_header_row=patient_header_row,
                    followup_records=followup_records,
                )
            )
            unsupported_followup_versions.update(patient_unsupported_versions)

        endpoint_prep_rows, overlap_audit_rows, patient_metrics = build_patient_level_outputs(
            endpoint_target_prep_v1_run_id=endpoint_target_prep_v1_run_id,
            workflow_inputs=workflow_inputs,
            parsed_patients_by_barcode=parsed_patients_by_barcode,
        )
        coverage_rows = build_followup_version_coverage_rows(
            endpoint_target_prep_v1_run_id=endpoint_target_prep_v1_run_id,
            endpoint_prep_rows=endpoint_prep_rows,
            overlap_audit_rows=overlap_audit_rows,
            patient_header_rows=patient_header_rows,
            patient_metrics=patient_metrics,
        )
        summary_rows = build_summary_rows(
            endpoint_target_prep_v1_run_id=endpoint_target_prep_v1_run_id,
            workflow_inputs=workflow_inputs,
            patient_header_rows=patient_header_rows,
            followup_long_rows=followup_long_rows,
            endpoint_prep_rows=endpoint_prep_rows,
            overlap_audit_rows=overlap_audit_rows,
            patient_metrics=patient_metrics,
        )
        validation_payload = validate_generated_outputs(
            workflow_inputs=workflow_inputs,
            bcr_xml_source_rows=bcr_xml_source_rows,
            patient_header_rows=patient_header_rows,
            followup_long_rows=followup_long_rows,
            endpoint_prep_rows=endpoint_prep_rows,
            overlap_audit_rows=overlap_audit_rows,
            summary_rows=summary_rows,
            unsupported_followup_versions=unsupported_followup_versions,
        )

        helper_module.write_dict_rows_tsv(
            patient_fields_path,
            PATIENT_ENDPOINT_FIELD_FIELDNAMES,
            patient_header_rows,
        )
        helper_module.write_dict_rows_tsv(
            followup_long_path,
            FOLLOWUP_LONG_FIELDNAMES,
            followup_long_rows,
        )
        helper_module.write_dict_rows_tsv(
            coverage_path,
            COVERAGE_FIELDNAMES,
            coverage_rows,
        )
        helper_module.write_dict_rows_tsv(
            endpoint_prep_path,
            build_endpoint_prep_fieldnames(),
            endpoint_prep_rows,
        )
        helper_module.write_dict_rows_tsv(
            overlap_audit_path,
            build_overlap_audit_fieldnames(),
            overlap_audit_rows,
        )
        helper_module.write_dict_rows_tsv(
            summary_path,
            SUMMARY_FIELDNAMES,
            summary_rows,
        )

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "clinical_xml_patient_endpoint_fields_tsv": patient_fields_path,
            "clinical_xml_followup_fields_long_tsv": followup_long_path,
            "clinical_xml_followup_version_coverage_tsv": coverage_path,
            "endpoint_target_prep_v1_tsv": endpoint_prep_path,
            "endpoint_target_prep_v1_overlap_audit_tsv": overlap_audit_path,
            "endpoint_target_prep_v1_summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            endpoint_target_prep_v1_run_id=endpoint_target_prep_v1_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
            "source_run_id": str(workflow_inputs.source_latest_pointer["run_id"]),
            "endpoint_crosswalk_run_id": str(
                workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
            ),
            "ambiguity_resolution_run_id": str(
                workflow_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
            ),
            "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
            "baseline_model_input_v1_run_id": str(
                workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
            ),
            "baseline_analysis_v1_run_id": str(
                workflow_inputs.baseline_model_input_run_log["baseline_analysis_v1_run_id"]
            ),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "source_supplements_latest_json": repo_relative(paths.source_latest_pointer, paths.repo_root),
                "endpoint_crosswalk_latest_json": repo_relative(
                    paths.endpoint_crosswalk_latest_pointer,
                    paths.repo_root,
                ),
                "ambiguity_resolution_latest_json": repo_relative(
                    paths.ambiguity_resolution_latest_pointer,
                    paths.repo_root,
                ),
                "minimal_cohort_v1_latest_json": repo_relative(
                    paths.minimal_cohort_latest_pointer,
                    paths.repo_root,
                ),
                "baseline_model_input_v1_latest_json": repo_relative(
                    paths.baseline_model_input_latest_pointer,
                    paths.repo_root,
                ),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "clinical_xml_patient_endpoint_fields_tsv": repo_relative(
                    patient_fields_path,
                    paths.repo_root,
                ),
                "clinical_xml_followup_fields_long_tsv": repo_relative(
                    followup_long_path,
                    paths.repo_root,
                ),
                "clinical_xml_followup_version_coverage_tsv": repo_relative(
                    coverage_path,
                    paths.repo_root,
                ),
                "endpoint_target_prep_v1_tsv": repo_relative(endpoint_prep_path, paths.repo_root),
                "endpoint_target_prep_v1_overlap_audit_tsv": repo_relative(
                    overlap_audit_path,
                    paths.repo_root,
                ),
                "endpoint_target_prep_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation_payload,
            "rules": {
                "unit_of_analysis": "patient/case",
                "anchor_universe": "minimal_cohort_v1",
                "patient_header_source_label": PATIENT_HEADER_SOURCE_LABEL,
                "followup_source_labels_json": json_list(
                    [FOLLOWUP_SOURCE_LABEL_BY_VERSION[version] for version in SUPPORTED_FOLLOWUP_VERSIONS]
                ),
                "supported_followup_versions_json": json_list(list(SUPPORTED_FOLLOWUP_VERSIONS)),
                "xml_data_format_in_scope": "BCR XML",
                "xml_data_format_excluded": "BCR OMF XML",
                "raw_input_files_immutable": True,
                "parsed_outputs_immutable": True,
                "patient_header_vs_followup_precedence_decision_deferred": True,
                "patient_header_vs_followup_values_kept_side_by_side": True,
                "days_to_last_followup_primary_comparable_field": True,
                "days_to_last_known_alive_auxiliary_last_contact_like_field": True,
                "no_final_endpoint_column_created": True,
                "no_endpoint_freeze": True,
                "no_modeling": True,
                "no_treatment_group_construction": True,
                "no_treatment_reintegration_into_model_matrix": True,
                "no_metabric": True,
                "current_biotab_v4_comparison_derived_from_minimal_cohort_v1_lineage": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "missing_like_tokens_json": json_list(sorted(MISSING_LIKE_TOKENS)),
                "uuid_normalization_for_join_and_audit": "uppercase_preserve_hyphenation",
            },
            "counts": {
                "source_metadata_row_count": len(workflow_inputs.source_metadata_rows),
                "bcr_xml_source_row_count": len(bcr_xml_source_rows),
                "minimal_cohort_row_count": len(workflow_inputs.minimal_cohort_rows),
                "baseline_feature_set_audit_map_row_count": len(
                    workflow_inputs.baseline_feature_set_audit_map_rows
                ),
                "baseline_model_input_row_count": workflow_inputs.baseline_model_input_row_count,
                "clinical_xml_patient_endpoint_fields_row_count": len(patient_header_rows),
                "clinical_xml_followup_fields_long_row_count": len(followup_long_rows),
                "clinical_xml_followup_version_coverage_row_count": len(coverage_rows),
                "endpoint_target_prep_v1_row_count": len(endpoint_prep_rows),
                "endpoint_target_prep_v1_overlap_audit_row_count": len(overlap_audit_rows),
                "endpoint_target_prep_v1_summary_row_count": len(summary_rows),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": build_upstream_snapshots(workflow_inputs),
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "endpoint_target_prep_v1_run_id": endpoint_target_prep_v1_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_clinical_xml_endpoint_target_prep_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA clinical XML endpoint-target prep v1 workflow complete.")
    print(f"Endpoint-target prep v1 run ID: {run_log['endpoint_target_prep_v1_run_id']}")
    print(f"Source run ID: {run_log['source_run_id']}")
    print(f"Cohort v1 build ID: {run_log['cohort_v1_build_id']}")
    print(f"Baseline model-input v1 run ID: {run_log['baseline_model_input_v1_run_id']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(
        "Patient endpoint fields TSV: "
        f"{run_log['outputs']['clinical_xml_patient_endpoint_fields_tsv']}"
    )
    print(
        "Follow-up long TSV: "
        f"{run_log['outputs']['clinical_xml_followup_fields_long_tsv']}"
    )
    print(
        "Endpoint prep TSV: "
        f"{run_log['outputs']['endpoint_target_prep_v1_tsv']}"
    )
    print(
        "Overlap audit TSV: "
        f"{run_log['outputs']['endpoint_target_prep_v1_overlap_audit_tsv']}"
    )
    print(f"Summary TSV: {run_log['outputs']['endpoint_target_prep_v1_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
