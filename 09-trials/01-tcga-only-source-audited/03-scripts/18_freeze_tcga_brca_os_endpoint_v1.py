#!/usr/bin/env python
"""Freeze an auditable TCGA-BRCA overall-survival endpoint v1 from saved prep outputs."""

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

SUPPORTED_FOLLOWUP_VERSIONS = ("1.5", "2.1", "4.0")
FOLLOWUP_VERSION_SUFFIX_BY_VERSION = {
    "1.5": "v1_5",
    "2.1": "v2_1",
    "4.0": "v4_0",
}
FOLLOWUP_SOURCE_LABEL_BY_VERSION = {
    "1.5": "followup_v1_5",
    "2.1": "followup_v2_1",
    "4.0": "followup_v4_0",
}
FOLLOWUP_VERSION_BY_SOURCE_LABEL = {
    source_label: version for version, source_label in FOLLOWUP_SOURCE_LABEL_BY_VERSION.items()
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
RULE_SPEC_FIELDNAMES = [
    "os_endpoint_v1_run_id",
    "rule_name",
    "rule_order",
    "rule_description",
    "fields_used",
    "rationale",
    "notes",
]
CONFLICT_AUDIT_FIELDNAMES = [
    "os_endpoint_v1_run_id",
    "cohort_v1_build_id",
    "baseline_model_input_v1_run_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "provisional_patient_row_id",
    "baseline_analysis_v1_row_id",
    "feature_set_v1_row_index",
    "conflict_type",
    "source_evidence_summary",
    "provisional_os_event",
    "provisional_os_time_days",
    "os_event_source",
    "os_time_source",
    "os_last_contact_source_max",
    "os_endpoint_inclusion_status",
    "recommended_handling",
    "manual_review_priority",
    "os_requires_manual_review",
]
SUMMARY_FIELDNAMES = [
    "os_endpoint_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
OS_ENDPOINT_FIELDNAMES = [
    "os_endpoint_v1_run_id",
    "cohort_v1_build_id",
    "baseline_model_input_v1_run_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "provisional_patient_row_id",
    "baseline_analysis_v1_row_id",
    "feature_set_v1_row_index",
    "os_event",
    "os_time_days",
    "os_event_source",
    "os_time_source",
    "os_last_contact_source_max",
    "os_conflict_vital_status",
    "os_time_missing",
    "os_conflict_dead_without_death_time",
    "os_conflict_time_without_clear_status",
    "os_conflict_event_time_inconsistent",
    "os_requires_manual_review",
    "os_endpoint_inclusion_status",
    "os_header_vital_status",
    "os_followup_vital_status_values_json",
    "os_header_days_to_death",
    "os_followup_max_days_to_death",
    "os_header_days_to_last_followup",
    "os_followup_max_days_to_last_followup",
    "os_followup_max_days_to_last_known_alive",
    "os_followup_versions_present_json",
    "os_followup_versions_contributing_time_json",
    "os_multiple_followup_versions_contributed",
    "os_followup_max_last_contact_exceeds_header",
    "os_rule_status_flags_json",
]
ENDPOINT_READINESS_INTERPRETATION = "join_ready_with_explicit_conflict_flags"
STRONGER_ENDPOINT_CLAIMS_STATUS = "blocked_pending_manual_reconciliation"
TREATMENT_FEASIBILITY_STATUS = "unchanged_not_addressed"
INCLUSION_STATUS_CLEAN = "included_rule_based_clean"
INCLUSION_STATUS_CONFLICT = "included_rule_based_conflict_flagged"
INCLUSION_STATUS_EXCLUDED_NO_TIME = "excluded_no_usable_time"
INCLUSION_STATUS_EXCLUDED_NO_STATUS = "excluded_no_clear_status"
HIGH_PRIORITY_CONFLICT_TYPES = {
    "vital_status_disagreement",
    "dead_without_death_time",
    "time_without_clear_status",
    "no_usable_time",
}
LOW_PRIORITY_CONFLICT_TYPE = "multiple_source_reconciliation"


class OSEndpointV1Error(RuntimeError):
    """Raised when the OS endpoint v1 workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the OS endpoint v1 workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    endpoint_target_prep_latest_pointer: Path
    minimal_cohort_latest_pointer: Path
    baseline_model_input_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    endpoint_target_prep_latest_pointer: dict[str, Any]
    endpoint_target_prep_run_log: dict[str, Any]
    minimal_cohort_latest_pointer: dict[str, Any]
    minimal_cohort_run_log: dict[str, Any]
    baseline_model_input_latest_pointer: dict[str, Any]
    baseline_model_input_run_log: dict[str, Any]
    patient_fields_rows: list[dict[str, str]]
    followup_long_rows: list[dict[str, str]]
    endpoint_prep_rows: list[dict[str, str]]
    overlap_audit_rows: list[dict[str, str]]
    endpoint_summary_rows: list[dict[str, str]]
    minimal_cohort_rows: list[dict[str, str]]
    baseline_model_input_rows: list[dict[str, str]]
    baseline_feature_set_audit_map_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise OSEndpointV1Error(f"Required helper script not found: {script_path}")

    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise OSEndpointV1Error(f"Unable to create an import spec for: {script_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


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


def parse_bool_string(value: str) -> bool:
    return value.strip().lower() == "true"


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def metric_slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.strip().lower())
    return slug.strip("_") or "value"


def normalize_missing_like(value: str) -> str:
    return value.strip().lower()


def is_missing_like(value: str) -> bool:
    return normalize_missing_like(value) in MISSING_LIKE_TOKENS


def normalize_barcode(value: str) -> str:
    return value.strip().upper()


def normalize_uuid(value: str) -> str:
    return value.strip().upper()


def parse_int_or_none(value: str) -> int | None:
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return int(stripped)
    except ValueError:
        return None


def parse_json_list_cell(value: str, label: str) -> list[Any]:
    stripped = value.strip()
    if not stripped:
        return []
    try:
        parsed = json.loads(stripped)
    except json.JSONDecodeError as exc:
        raise OSEndpointV1Error(f"Unable to parse JSON list for {label}: {value!r}") from exc
    if not isinstance(parsed, list):
        raise OSEndpointV1Error(f"Expected a JSON list for {label}: {value!r}")
    return parsed


def parse_candidate_items(value: str, label: str) -> list[dict[str, str]]:
    parsed = parse_json_list_cell(value, label)
    normalized: list[dict[str, str]] = []
    for index, item in enumerate(parsed, start=1):
        if isinstance(item, dict):
            normalized_item = {
                str(key): "" if item_value is None else str(item_value).strip()
                for key, item_value in item.items()
            }
        else:
            normalized_item = {"value": str(item).strip()}
        normalized_item.setdefault("value", "")
        if is_missing_like(normalized_item["value"]):
            continue
        normalized_item["candidate_index"] = str(index)
        normalized.append(normalized_item)
    return normalized


def candidate_values(candidate_items: list[dict[str, str]]) -> list[str]:
    return [item["value"] for item in candidate_items if item["value"]]


def distinct_candidate_values(candidate_items: list[dict[str, str]]) -> list[str]:
    return ordered_unique(candidate_values(candidate_items))


def distinct_values(values: list[str]) -> list[str]:
    return ordered_unique([value for value in values if not is_missing_like(value)])


def max_int(values: list[int]) -> int | None:
    if not values:
        return None
    return max(values)


def numeric_values_from_candidate_items(candidate_items: list[dict[str, str]]) -> list[int]:
    numeric_values: list[int] = []
    for item in candidate_items:
        parsed = parse_int_or_none(item["value"])
        if parsed is not None:
            numeric_values.append(parsed)
    return numeric_values


def followup_version_from_evidence_label(label: str) -> str | None:
    for version, source_label in FOLLOWUP_SOURCE_LABEL_BY_VERSION.items():
        if label.startswith(f"{source_label}_"):
            return version
    return None


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
            raise OSEndpointV1Error(f"{key_name} produced an empty key while building a lookup.")
        if key in lookup:
            raise OSEndpointV1Error(f"Duplicate {key_name} detected: {key}")
        lookup[key] = row
    return lookup


def require_pointer(
    pointer_path: Path,
    *,
    label: str,
    required_keys: set[str],
    helper_module: Any,
) -> dict[str, Any]:
    if not pointer_path.exists():
        raise OSEndpointV1Error(f"Required {label} latest pointer not found: {pointer_path}")
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
        raise OSEndpointV1Error(f"{label} is not completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise OSEndpointV1Error(f"{label} does not report validation.passed == true.")
    return run_log_path, run_log


def require_columns(rows: list[dict[str, str]], required_columns: set[str], label: str) -> None:
    if not rows:
        raise OSEndpointV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise OSEndpointV1Error(f"{label} is missing required columns: {sorted(missing)}")


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
        raise OSEndpointV1Error(f"Required trial config not found: {trial_config}")

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
        processed_runs_root=processed_root / "tcga-brca" / "endpoint-prep" / "os_endpoint_v1_runs",
        audit_runs_root=endpoint_prep_root / "os_endpoint_v1_runs",
        latest_pointer=endpoint_prep_root / "tcga_brca_os_endpoint_v1_latest.json",
        endpoint_target_prep_latest_pointer=endpoint_prep_root / "tcga_brca_endpoint_target_prep_v1_latest.json",
        minimal_cohort_latest_pointer=audit_root / "tcga-brca" / "cohort" / "tcga_brca_minimal_cohort_v1_latest.json",
        baseline_model_input_latest_pointer=(
            audit_root / "tcga-brca" / "model-input" / "tcga_brca_baseline_model_input_v1_latest.json"
        ),
    )


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    endpoint_target_prep_latest_pointer = require_pointer(
        paths.endpoint_target_prep_latest_pointer,
        label="endpoint-target prep v1",
        required_keys={
            "endpoint_target_prep_v1_run_id",
            "cohort_v1_build_id",
            "baseline_model_input_v1_run_id",
            "clinical_xml_patient_endpoint_fields_tsv",
            "clinical_xml_followup_fields_long_tsv",
            "endpoint_target_prep_v1_tsv",
            "endpoint_target_prep_v1_overlap_audit_tsv",
            "endpoint_target_prep_v1_summary_tsv",
            "run_log_json",
            "source_run_id",
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

    _, endpoint_target_prep_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(endpoint_target_prep_latest_pointer["run_log_json"]),
        label="endpoint-target prep v1 run log",
        helper_module=helper_module,
    )
    _, minimal_cohort_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(minimal_cohort_latest_pointer["run_log_json"]),
        label="minimal cohort v1 run log",
        helper_module=helper_module,
    )
    _, baseline_model_input_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(baseline_model_input_latest_pointer["run_log_json"]),
        label="baseline model-input v1 run log",
        helper_module=helper_module,
    )

    input_paths = {
        "clinical_xml_patient_endpoint_fields_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(endpoint_target_prep_latest_pointer["clinical_xml_patient_endpoint_fields_tsv"]),
            "clinical_xml_patient_endpoint_fields.tsv",
        ),
        "clinical_xml_followup_fields_long_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(endpoint_target_prep_latest_pointer["clinical_xml_followup_fields_long_tsv"]),
            "clinical_xml_followup_fields_long.tsv",
        ),
        "endpoint_target_prep_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(endpoint_target_prep_latest_pointer["endpoint_target_prep_v1_tsv"]),
            "endpoint_target_prep_v1.tsv",
        ),
        "endpoint_target_prep_v1_overlap_audit_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(endpoint_target_prep_latest_pointer["endpoint_target_prep_v1_overlap_audit_tsv"]),
            "endpoint_target_prep_v1_overlap_audit.tsv",
        ),
        "endpoint_target_prep_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(endpoint_target_prep_latest_pointer["endpoint_target_prep_v1_summary_tsv"]),
            "endpoint_target_prep_v1_summary.tsv",
        ),
        "minimal_cohort_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["minimal_cohort_v1_tsv"]),
            "minimal_cohort_v1.tsv",
        ),
        "baseline_model_input_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_model_input_latest_pointer["baseline_model_input_v1_tsv"]),
            "baseline_model_input_v1.tsv",
        ),
        "baseline_feature_set_v1_audit_map_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_model_input_latest_pointer["baseline_feature_set_v1_audit_map_tsv"]),
            "baseline_feature_set_v1_audit_map.tsv",
        ),
    }

    patient_fields_rows = helper_module.read_tsv_dict_rows(input_paths["clinical_xml_patient_endpoint_fields_tsv"])
    followup_long_rows = helper_module.read_tsv_dict_rows(input_paths["clinical_xml_followup_fields_long_tsv"])
    endpoint_prep_rows = helper_module.read_tsv_dict_rows(input_paths["endpoint_target_prep_v1_tsv"])
    overlap_audit_rows = helper_module.read_tsv_dict_rows(input_paths["endpoint_target_prep_v1_overlap_audit_tsv"])
    endpoint_summary_rows = helper_module.read_tsv_dict_rows(input_paths["endpoint_target_prep_v1_summary_tsv"])
    minimal_cohort_rows = helper_module.read_tsv_dict_rows(input_paths["minimal_cohort_v1_tsv"])
    baseline_model_input_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_model_input_v1_tsv"])
    baseline_feature_set_audit_map_rows = helper_module.read_tsv_dict_rows(
        input_paths["baseline_feature_set_v1_audit_map_tsv"]
    )

    require_columns(
        patient_fields_rows,
        {
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "vital_status",
            "days_to_last_followup",
            "days_to_last_known_alive",
            "days_to_death",
        },
        "clinical_xml_patient_endpoint_fields.tsv",
    )
    require_columns(
        followup_long_rows,
        {
            "bcr_patient_barcode",
            "followup_version",
            "xml_field_name",
            "xml_field_value",
        },
        "clinical_xml_followup_fields_long.tsv",
    )
    require_columns(
        endpoint_prep_rows,
        {
            "cohort_v1_build_id",
            "baseline_model_input_v1_run_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
            "patient_header_vital_status",
            "patient_header_days_to_last_followup",
            "patient_header_days_to_last_known_alive",
            "patient_header_days_to_death",
            "xml_followup_versions_present_json",
            "xml_followup_max_days_to_last_followup",
            "xml_followup_max_days_to_last_known_alive",
            "xml_followup_max_days_to_death",
            "xml_followup_max_last_contact_exceeds_header",
        },
        "endpoint_target_prep_v1.tsv",
    )
    require_columns(
        overlap_audit_rows,
        {
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
        },
        "endpoint_target_prep_v1_overlap_audit.tsv",
    )
    require_columns(
        minimal_cohort_rows,
        {
            "cohort_v1_build_id",
            "provisional_patient_row_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
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
        endpoint_target_prep_latest_pointer=endpoint_target_prep_latest_pointer,
        endpoint_target_prep_run_log=endpoint_target_prep_run_log,
        minimal_cohort_latest_pointer=minimal_cohort_latest_pointer,
        minimal_cohort_run_log=minimal_cohort_run_log,
        baseline_model_input_latest_pointer=baseline_model_input_latest_pointer,
        baseline_model_input_run_log=baseline_model_input_run_log,
        patient_fields_rows=patient_fields_rows,
        followup_long_rows=followup_long_rows,
        endpoint_prep_rows=endpoint_prep_rows,
        overlap_audit_rows=overlap_audit_rows,
        endpoint_summary_rows=endpoint_summary_rows,
        minimal_cohort_rows=minimal_cohort_rows,
        baseline_model_input_rows=baseline_model_input_rows,
        baseline_feature_set_audit_map_rows=baseline_feature_set_audit_map_rows,
        input_paths=input_paths,
    )


def validate_upstream_state(workflow_inputs: WorkflowInputs) -> None:
    endpoint_prep_pointer = workflow_inputs.endpoint_target_prep_latest_pointer
    minimal_pointer = workflow_inputs.minimal_cohort_latest_pointer
    baseline_pointer = workflow_inputs.baseline_model_input_latest_pointer

    if str(endpoint_prep_pointer["cohort_v1_build_id"]) != str(minimal_pointer["cohort_v1_build_id"]):
        raise OSEndpointV1Error(
            "endpoint-target prep v1 and minimal cohort v1 latest pointers reference different cohort_v1_build_id values."
        )
    if str(endpoint_prep_pointer["baseline_model_input_v1_run_id"]) != str(
        baseline_pointer["baseline_model_input_v1_run_id"]
    ):
        raise OSEndpointV1Error(
            "endpoint-target prep v1 and baseline model-input v1 latest pointers reference different baseline_model_input_v1_run_id values."
        )
    if len(workflow_inputs.endpoint_prep_rows) != len(workflow_inputs.minimal_cohort_rows):
        raise OSEndpointV1Error("endpoint_target_prep_v1.tsv row count does not match minimal_cohort_v1.tsv.")
    if len(workflow_inputs.baseline_model_input_rows) != len(workflow_inputs.minimal_cohort_rows):
        raise OSEndpointV1Error("baseline_model_input_v1.tsv row count does not match minimal_cohort_v1.tsv.")
    if len(workflow_inputs.patient_fields_rows) != len(workflow_inputs.minimal_cohort_rows):
        raise OSEndpointV1Error(
            "clinical_xml_patient_endpoint_fields.tsv row count does not match minimal_cohort_v1.tsv."
        )
    if len(workflow_inputs.overlap_audit_rows) != len(workflow_inputs.endpoint_prep_rows):
        raise OSEndpointV1Error(
            "endpoint_target_prep_v1_overlap_audit.tsv row count does not match endpoint_target_prep_v1.tsv."
        )


def get_followup_candidate_items(
    endpoint_row: dict[str, str],
    *,
    version: str,
    field_name: str,
) -> list[dict[str, str]]:
    version_suffix = FOLLOWUP_VERSION_SUFFIX_BY_VERSION[version]
    column_name = f"followup_{version_suffix}_{field_name}_candidates_json"
    return parse_candidate_items(str(endpoint_row[column_name] or ""), column_name)


def build_rule_spec_rows(os_endpoint_v1_run_id: str) -> list[dict[str, str]]:
    rows = [
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "anchor_universe",
            "rule_order": "1",
            "rule_description": "Freeze one OS row per patient in the current minimal cohort / baseline lineage.",
            "fields_used": json_list(
                [
                    "endpoint_target_prep_v1.tsv",
                    "minimal_cohort_v1.tsv",
                    "baseline_feature_set_v1_audit_map.tsv",
                ]
            ),
            "rationale": "The OS target must stay join-compatible with the current X-side cohort and model-input lineage.",
            "notes": "No treatment fields, recurrence fields, modeling fields, or METABRIC work are introduced here.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "event_dead_if_any_dead_signal",
            "rule_order": "2",
            "rule_description": "Set os_event=1 if any patient-header or follow-up vital-status value is Dead, or any days_to_death value is non-missing.",
            "fields_used": json_list(
                [
                    "patient_header_vital_status",
                    "followup_v1_5_vital_status_candidates_json",
                    "followup_v2_1_vital_status_candidates_json",
                    "followup_v4_0_vital_status_candidates_json",
                    "patient_header_days_to_death",
                    "followup_v1_5_days_to_death_candidates_json",
                    "followup_v2_1_days_to_death_candidates_json",
                    "followup_v4_0_days_to_death_candidates_json",
                ]
            ),
            "rationale": "The v1 OS freeze uses a conservative death rule and does not silently discard older XML layers.",
            "notes": "Death-time evidence itself is treated as death evidence even if Alive appears elsewhere.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "event_alive_only_if_consistent",
            "rule_order": "3",
            "rule_description": "Set os_event=0 only when all non-missing vital-status evidence is Alive and no days_to_death evidence exists.",
            "fields_used": json_list(
                [
                    "patient_header_vital_status",
                    "followup_v1_5_vital_status_candidates_json",
                    "followup_v2_1_vital_status_candidates_json",
                    "followup_v4_0_vital_status_candidates_json",
                    "patient_header_days_to_death",
                    "followup_v1_5_days_to_death_candidates_json",
                    "followup_v2_1_days_to_death_candidates_json",
                    "followup_v4_0_days_to_death_candidates_json",
                ]
            ),
            "rationale": "Alive censoring should only be assigned when the available status evidence is internally consistent.",
            "notes": "Patients with no clear status remain flagged instead of being force-censored.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "event_blank_if_no_clear_status",
            "rule_order": "4",
            "rule_description": "Leave os_event blank if neither the death rule nor the Alive-only rule can be satisfied.",
            "fields_used": json_list(["os_event", "os_conflict_time_without_clear_status", "os_endpoint_inclusion_status"]),
            "rationale": "The workflow must keep unresolved cases explicit rather than guessing a censor/event state.",
            "notes": "These rows remain in the patient-level table but are marked excluded_no_clear_status.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "time_use_max_days_to_death",
            "rule_order": "5",
            "rule_description": "If any days_to_death value exists, set os_time_days to the maximum available days_to_death across header and follow-up versions.",
            "fields_used": json_list(
                [
                    "patient_header_days_to_death",
                    "followup_v1_5_days_to_death_candidates_json",
                    "followup_v2_1_days_to_death_candidates_json",
                    "followup_v4_0_days_to_death_candidates_json",
                ]
            ),
            "rationale": "The freeze rule keeps all XML layers in play and records exact death-time provenance instead of imposing silent precedence.",
            "notes": "Tied maxima keep multiple source labels in os_time_source.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "time_else_use_max_last_contact_like",
            "rule_order": "6",
            "rule_description": "If no days_to_death exists, set os_time_days to the maximum last-contact-like time across header and follow-up days_to_last_followup and days_to_last_known_alive.",
            "fields_used": json_list(
                [
                    "patient_header_days_to_last_followup",
                    "patient_header_days_to_last_known_alive",
                    "followup_v1_5_days_to_last_followup_candidates_json",
                    "followup_v2_1_days_to_last_followup_candidates_json",
                    "followup_v4_0_days_to_last_followup_candidates_json",
                    "followup_v1_5_days_to_last_known_alive_candidates_json",
                    "followup_v2_1_days_to_last_known_alive_candidates_json",
                    "followup_v4_0_days_to_last_known_alive_candidates_json",
                ]
            ),
            "rationale": "The v1 OS freeze uses a conservative censoring fallback and explicitly preserves the max last-contact provenance.",
            "notes": "days_to_last_followup remains primary; days_to_last_known_alive acts as auxiliary evidence only if it wins the maximum rule.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "conflict_vital_status",
            "rule_order": "7",
            "rule_description": "Flag os_conflict_vital_status=yes when Alive and Dead both appear across patient-header and follow-up evidence.",
            "fields_used": json_list(
                [
                    "os_header_vital_status",
                    "os_followup_vital_status_values_json",
                    "os_conflict_vital_status",
                    "os_requires_manual_review",
                ]
            ),
            "rationale": "The workflow must keep cross-source status disagreement explicit for manual endpoint reconciliation.",
            "notes": "Conflict rows retain provisional values under the Include + Flag policy.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "conflict_dead_without_death_time",
            "rule_order": "8",
            "rule_description": "Flag os_conflict_dead_without_death_time=yes when os_event=1 but the chosen os_time_days fell back to last-contact-like evidence because no days_to_death exists.",
            "fields_used": json_list(
                [
                    "os_event",
                    "os_time_days",
                    "os_time_source",
                    "os_conflict_dead_without_death_time",
                    "os_conflict_event_time_inconsistent",
                ]
            ),
            "rationale": "This is a weaker event/time pairing and should remain explicit even though a provisional OS value is retained.",
            "notes": "The current saved data is expected to have one such case.",
        },
        {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "rule_name": "inclusion_status_mapping",
            "rule_order": "9",
            "rule_description": "Assign inclusion status from the resolved event/time state and conflict flags.",
            "fields_used": json_list(
                [
                    "os_time_missing",
                    "os_conflict_time_without_clear_status",
                    "os_requires_manual_review",
                    "os_endpoint_inclusion_status",
                ]
            ),
            "rationale": "Downstream users need an explicit, auditable switch for filtering clean versus flagged OS rows.",
            "notes": "Allowed values are included_rule_based_clean, included_rule_based_conflict_flagged, excluded_no_usable_time, and excluded_no_clear_status.",
        },
    ]
    return rows


def build_os_outputs(
    *,
    os_endpoint_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    minimal_lookup = build_unique_lookup(
        workflow_inputs.minimal_cohort_rows,
        key_name="minimal cohort provisional_patient_row_id",
        key_builder=lambda row: str(row["provisional_patient_row_id"]).strip(),
    )
    endpoint_lookup = build_unique_lookup(
        workflow_inputs.endpoint_prep_rows,
        key_name="endpoint prep provisional_patient_row_id",
        key_builder=lambda row: str(row["provisional_patient_row_id"]).strip(),
    )
    audit_map_lookup = build_unique_lookup(
        workflow_inputs.baseline_feature_set_audit_map_rows,
        key_name="baseline audit-map provisional_patient_row_id",
        key_builder=lambda row: str(row["provisional_patient_row_id"]).strip(),
    )

    os_rows: list[dict[str, str]] = []
    conflict_rows: list[dict[str, str]] = []

    summary_counters: Counter[str] = Counter()
    time_source_counter: Counter[str] = Counter()
    time_type_counter: Counter[str] = Counter()
    event_source_counter: Counter[str] = Counter()
    inclusion_counter: Counter[str] = Counter()
    conflict_type_patients: dict[str, set[str]] = {
        conflict_type: set()
        for conflict_type in [*HIGH_PRIORITY_CONFLICT_TYPES, LOW_PRIORITY_CONFLICT_TYPE]
    }
    manual_review_patients: set[str] = set()
    patient_conflict_audit_coverage: dict[str, set[str]] = {}

    for provisional_patient_row_id in sorted(minimal_lookup.keys(), key=lambda value: int(value)):
        minimal_row = minimal_lookup[provisional_patient_row_id]
        endpoint_row = endpoint_lookup.get(provisional_patient_row_id)
        bridge_row = audit_map_lookup.get(provisional_patient_row_id)

        if endpoint_row is None:
            raise OSEndpointV1Error(
                f"Missing endpoint-target prep row for provisional_patient_row_id={provisional_patient_row_id}."
            )
        if bridge_row is None:
            raise OSEndpointV1Error(
                f"Missing baseline feature-set audit-map row for provisional_patient_row_id={provisional_patient_row_id}."
            )

        minimal_barcode = normalize_barcode(str(minimal_row["bcr_patient_barcode"]))
        minimal_uuid = normalize_uuid(str(minimal_row["bcr_patient_uuid"]))
        if normalize_barcode(str(endpoint_row["bcr_patient_barcode"])) != minimal_barcode:
            raise OSEndpointV1Error(
                f"Barcode mismatch between minimal cohort and endpoint prep for provisional_patient_row_id={provisional_patient_row_id}."
            )
        if normalize_uuid(str(endpoint_row["bcr_patient_uuid"])) != minimal_uuid:
            raise OSEndpointV1Error(
                f"UUID mismatch between minimal cohort and endpoint prep for provisional_patient_row_id={provisional_patient_row_id}."
            )
        if normalize_barcode(str(bridge_row["bcr_patient_barcode"])) != minimal_barcode:
            raise OSEndpointV1Error(
                f"Barcode mismatch between minimal cohort and baseline audit map for provisional_patient_row_id={provisional_patient_row_id}."
            )
        if normalize_uuid(str(bridge_row["bcr_patient_uuid"])) != minimal_uuid:
            raise OSEndpointV1Error(
                f"UUID mismatch between minimal cohort and baseline audit map for provisional_patient_row_id={provisional_patient_row_id}."
            )
        if str(endpoint_row["baseline_analysis_v1_row_id"]).strip() != str(bridge_row["baseline_analysis_v1_row_id"]).strip():
            raise OSEndpointV1Error(
                f"baseline_analysis_v1_row_id mismatch for provisional_patient_row_id={provisional_patient_row_id}."
            )
        if str(endpoint_row["feature_set_v1_row_index"]).strip() != str(bridge_row["feature_set_v1_row_index"]).strip():
            raise OSEndpointV1Error(
                f"feature_set_v1_row_index mismatch for provisional_patient_row_id={provisional_patient_row_id}."
            )

        header_status = str(endpoint_row["patient_header_vital_status"] or "").strip()
        header_days_to_death_raw = str(endpoint_row["patient_header_days_to_death"] or "").strip()
        header_days_to_last_followup_raw = str(endpoint_row["patient_header_days_to_last_followup"] or "").strip()
        header_days_to_last_known_alive_raw = str(endpoint_row["patient_header_days_to_last_known_alive"] or "").strip()

        header_days_to_death = parse_int_or_none(header_days_to_death_raw)
        header_days_to_last_followup = parse_int_or_none(header_days_to_last_followup_raw)
        header_days_to_last_known_alive = parse_int_or_none(header_days_to_last_known_alive_raw)

        followup_versions_present = [
            version
            for version in parse_json_list_cell(
                str(endpoint_row["xml_followup_versions_present_json"] or ""),
                "xml_followup_versions_present_json",
            )
            if str(version) in SUPPORTED_FOLLOWUP_VERSIONS
        ]
        followup_versions_present = ordered_unique([str(version) for version in followup_versions_present])

        followup_status_values_by_version: dict[str, list[str]] = {}
        followup_days_to_death_by_version: dict[str, list[int]] = {}
        followup_days_to_last_followup_by_version: dict[str, list[int]] = {}
        followup_days_to_last_known_alive_by_version: dict[str, list[int]] = {}

        followup_vital_values_overall: list[str] = []
        followup_alive_source_labels: list[str] = []
        followup_dead_source_labels: list[str] = []
        followup_days_to_death_source_labels: list[str] = []

        for version in SUPPORTED_FOLLOWUP_VERSIONS:
            status_items = get_followup_candidate_items(endpoint_row, version=version, field_name="vital_status")
            days_to_death_items = get_followup_candidate_items(endpoint_row, version=version, field_name="days_to_death")
            days_to_last_followup_items = get_followup_candidate_items(
                endpoint_row,
                version=version,
                field_name="days_to_last_followup",
            )
            days_to_last_known_alive_items = get_followup_candidate_items(
                endpoint_row,
                version=version,
                field_name="days_to_last_known_alive",
            )

            source_label = FOLLOWUP_SOURCE_LABEL_BY_VERSION[version]
            status_values = distinct_candidate_values(status_items)
            followup_status_values_by_version[version] = status_values
            followup_vital_values_overall.extend(status_values)

            if "Alive" in status_values:
                followup_alive_source_labels.append(f"{source_label}_vital_status")
            if "Dead" in status_values:
                followup_dead_source_labels.append(f"{source_label}_vital_status")

            death_values = numeric_values_from_candidate_items(days_to_death_items)
            followup_days_to_death_by_version[version] = death_values
            if death_values:
                followup_days_to_death_source_labels.append(f"{source_label}_days_to_death")

            followup_days_to_last_followup_by_version[version] = numeric_values_from_candidate_items(
                days_to_last_followup_items
            )
            followup_days_to_last_known_alive_by_version[version] = numeric_values_from_candidate_items(
                days_to_last_known_alive_items
            )

        followup_vital_values = distinct_values(followup_vital_values_overall)
        followup_versions_with_days_to_death = [
            version for version in SUPPORTED_FOLLOWUP_VERSIONS if followup_days_to_death_by_version[version]
        ]
        status_values_all = distinct_values([header_status, *followup_vital_values])
        has_header_dead = header_status == "Dead"
        has_header_alive = header_status == "Alive"
        has_any_days_to_death = header_days_to_death is not None or bool(followup_versions_with_days_to_death)
        any_dead_signal = has_header_dead or bool(followup_dead_source_labels) or has_any_days_to_death
        all_statuses_alive_only = bool(status_values_all) and set(status_values_all) == {"Alive"} and not has_any_days_to_death

        if any_dead_signal:
            os_event = "1"
        elif all_statuses_alive_only:
            os_event = "0"
        else:
            os_event = ""

        event_source_labels: list[str] = []
        if os_event == "1":
            if has_header_dead:
                event_source_labels.append("patient_header_vital_status")
            event_source_labels.extend(followup_dead_source_labels)
            if header_days_to_death is not None:
                event_source_labels.append("patient_header_days_to_death")
            event_source_labels.extend(followup_days_to_death_source_labels)
        elif os_event == "0":
            if has_header_alive:
                event_source_labels.append("patient_header_vital_status")
            event_source_labels.extend(followup_alive_source_labels)
        else:
            if header_status:
                event_source_labels.append("patient_header_vital_status")
            event_source_labels.extend(
                [
                    f"{FOLLOWUP_SOURCE_LABEL_BY_VERSION[version]}_vital_status"
                    for version in SUPPORTED_FOLLOWUP_VERSIONS
                    if followup_status_values_by_version[version]
                ]
            )
        event_source_labels = ordered_unique(event_source_labels)

        death_time_candidates: list[int] = []
        if header_days_to_death is not None:
            death_time_candidates.append(header_days_to_death)
        for version in SUPPORTED_FOLLOWUP_VERSIONS:
            death_time_candidates.extend(followup_days_to_death_by_version[version])
        max_days_to_death = max_int(death_time_candidates)

        last_contact_candidates: list[int] = []
        if header_days_to_last_followup is not None:
            last_contact_candidates.append(header_days_to_last_followup)
        if header_days_to_last_known_alive is not None:
            last_contact_candidates.append(header_days_to_last_known_alive)
        for version in SUPPORTED_FOLLOWUP_VERSIONS:
            last_contact_candidates.extend(followup_days_to_last_followup_by_version[version])
            last_contact_candidates.extend(followup_days_to_last_known_alive_by_version[version])
        max_last_contact = max_int(last_contact_candidates)

        last_contact_source_labels: list[str] = []
        if max_last_contact is not None:
            if header_days_to_last_followup == max_last_contact:
                last_contact_source_labels.append("patient_header_days_to_last_followup")
            if header_days_to_last_known_alive == max_last_contact:
                last_contact_source_labels.append("patient_header_days_to_last_known_alive")
            for version in SUPPORTED_FOLLOWUP_VERSIONS:
                source_label = FOLLOWUP_SOURCE_LABEL_BY_VERSION[version]
                if any(value == max_last_contact for value in followup_days_to_last_followup_by_version[version]):
                    last_contact_source_labels.append(f"{source_label}_days_to_last_followup")
                if any(value == max_last_contact for value in followup_days_to_last_known_alive_by_version[version]):
                    last_contact_source_labels.append(f"{source_label}_days_to_last_known_alive")
        last_contact_source_labels = ordered_unique(last_contact_source_labels)

        if max_days_to_death is not None:
            os_time_days = str(max_days_to_death)
            time_source_labels: list[str] = []
            if header_days_to_death == max_days_to_death:
                time_source_labels.append("patient_header_days_to_death")
            for version in SUPPORTED_FOLLOWUP_VERSIONS:
                source_label = FOLLOWUP_SOURCE_LABEL_BY_VERSION[version]
                if any(value == max_days_to_death for value in followup_days_to_death_by_version[version]):
                    time_source_labels.append(f"{source_label}_days_to_death")
            time_type = "days_to_death"
        elif max_last_contact is not None:
            os_time_days = str(max_last_contact)
            time_source_labels = list(last_contact_source_labels)
            time_type = "last_contact"
        else:
            os_time_days = ""
            time_source_labels = []
            time_type = "missing"
        time_source_labels = ordered_unique(time_source_labels)

        followup_versions_contributing_time = ordered_unique(
            [
                version
                for version in [followup_version_from_evidence_label(label) for label in time_source_labels]
                if version is not None
            ]
        )
        followup_versions_contributing_endpoint_evidence = ordered_unique(
            [
                version
                for version in [
                    followup_version_from_evidence_label(label)
                    for label in ordered_unique([*event_source_labels, *time_source_labels, *last_contact_source_labels])
                ]
                if version is not None
            ]
        )

        os_conflict_vital_status = "Alive" in status_values_all and "Dead" in status_values_all
        os_time_missing = os_time_days == ""
        os_conflict_time_without_clear_status = os_event == "" and not os_time_missing
        os_conflict_dead_without_death_time = os_event == "1" and max_days_to_death is None and not os_time_missing
        os_conflict_event_time_inconsistent = (
            (os_event == "1" and max_days_to_death is None and not os_time_missing)
            or (os_event == "0" and max_days_to_death is not None)
        )
        os_requires_manual_review = (
            os_conflict_vital_status
            or os_time_missing
            or os_conflict_dead_without_death_time
            or os_conflict_time_without_clear_status
            or os_conflict_event_time_inconsistent
        )

        if os_time_missing:
            os_endpoint_inclusion_status = INCLUSION_STATUS_EXCLUDED_NO_TIME
        elif os_event == "":
            os_endpoint_inclusion_status = INCLUSION_STATUS_EXCLUDED_NO_STATUS
        elif os_requires_manual_review:
            os_endpoint_inclusion_status = INCLUSION_STATUS_CONFLICT
        else:
            os_endpoint_inclusion_status = INCLUSION_STATUS_CLEAN

        os_multiple_followup_versions_contributed = len(followup_versions_contributing_endpoint_evidence) > 1
        os_followup_max_last_contact_exceeds_header = parse_bool_string(
            str(endpoint_row["xml_followup_max_last_contact_exceeds_header"] or "")
        )
        multiple_source_reconciliation = (
            not os_requires_manual_review
            and (
                len(time_source_labels) > 1
                or len(last_contact_source_labels) > 1
                or os_multiple_followup_versions_contributed
            )
        )

        rule_status_flags = ordered_unique(
            [
                label
                for label, is_true in [
                    ("os_conflict_vital_status", os_conflict_vital_status),
                    ("os_time_missing", os_time_missing),
                    ("os_conflict_dead_without_death_time", os_conflict_dead_without_death_time),
                    ("os_conflict_time_without_clear_status", os_conflict_time_without_clear_status),
                    ("os_conflict_event_time_inconsistent", os_conflict_event_time_inconsistent),
                    ("os_requires_manual_review", os_requires_manual_review),
                    ("os_multiple_followup_versions_contributed", os_multiple_followup_versions_contributed),
                    ("os_followup_max_last_contact_exceeds_header", os_followup_max_last_contact_exceeds_header),
                    ("multiple_source_reconciliation", multiple_source_reconciliation),
                ]
                if is_true
            ]
        )

        os_row = {
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
            "baseline_model_input_v1_run_id": str(
                workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
            ),
            "bcr_patient_barcode": minimal_barcode,
            "bcr_patient_uuid": minimal_uuid,
            "provisional_patient_row_id": provisional_patient_row_id,
            "baseline_analysis_v1_row_id": str(endpoint_row["baseline_analysis_v1_row_id"]).strip(),
            "feature_set_v1_row_index": str(endpoint_row["feature_set_v1_row_index"]).strip(),
            "os_event": os_event,
            "os_time_days": os_time_days,
            "os_event_source": "|".join(event_source_labels),
            "os_time_source": "|".join(time_source_labels),
            "os_last_contact_source_max": "|".join(last_contact_source_labels),
            "os_conflict_vital_status": yes_no(os_conflict_vital_status),
            "os_time_missing": yes_no(os_time_missing),
            "os_conflict_dead_without_death_time": yes_no(os_conflict_dead_without_death_time),
            "os_conflict_time_without_clear_status": yes_no(os_conflict_time_without_clear_status),
            "os_conflict_event_time_inconsistent": yes_no(os_conflict_event_time_inconsistent),
            "os_requires_manual_review": yes_no(os_requires_manual_review),
            "os_endpoint_inclusion_status": os_endpoint_inclusion_status,
            "os_header_vital_status": header_status,
            "os_followup_vital_status_values_json": json_list(followup_vital_values),
            "os_header_days_to_death": header_days_to_death_raw,
            "os_followup_max_days_to_death": str(endpoint_row["xml_followup_max_days_to_death"] or "").strip(),
            "os_header_days_to_last_followup": header_days_to_last_followup_raw,
            "os_followup_max_days_to_last_followup": str(
                endpoint_row["xml_followup_max_days_to_last_followup"] or ""
            ).strip(),
            "os_followup_max_days_to_last_known_alive": str(
                endpoint_row["xml_followup_max_days_to_last_known_alive"] or ""
            ).strip(),
            "os_followup_versions_present_json": json_list(followup_versions_present),
            "os_followup_versions_contributing_time_json": json_list(followup_versions_contributing_time),
            "os_multiple_followup_versions_contributed": yes_no(os_multiple_followup_versions_contributed),
            "os_followup_max_last_contact_exceeds_header": yes_no(os_followup_max_last_contact_exceeds_header),
            "os_rule_status_flags_json": json_list(rule_status_flags),
        }
        os_rows.append(os_row)

        summary_counters["os_rows"] += 1
        if os_event == "1":
            summary_counters["os_event_1"] += 1
        elif os_event == "0":
            summary_counters["os_event_0"] += 1
        else:
            summary_counters["os_event_missing"] += 1

        if os_time_missing:
            summary_counters["os_time_missing"] += 1
        else:
            summary_counters["os_time_nonmissing"] += 1

        if time_type == "days_to_death":
            summary_counters["time_using_days_to_death"] += 1
        elif time_type == "last_contact":
            summary_counters["time_using_last_contact"] += 1
        else:
            summary_counters["time_using_missing"] += 1

        if os_conflict_vital_status:
            summary_counters["os_conflict_vital_status"] += 1
            conflict_type_patients["vital_status_disagreement"].add(minimal_barcode)
        if os_conflict_dead_without_death_time:
            summary_counters["os_conflict_dead_without_death_time"] += 1
            conflict_type_patients["dead_without_death_time"].add(minimal_barcode)
        if os_conflict_time_without_clear_status:
            summary_counters["os_conflict_time_without_clear_status"] += 1
            conflict_type_patients["time_without_clear_status"].add(minimal_barcode)
        if os_time_missing:
            summary_counters["os_time_missing_patients"] += 1
            conflict_type_patients["no_usable_time"].add(minimal_barcode)
        if os_conflict_event_time_inconsistent:
            summary_counters["os_conflict_event_time_inconsistent"] += 1
        if os_requires_manual_review:
            summary_counters["os_requires_manual_review"] += 1
            manual_review_patients.add(minimal_barcode)
        if multiple_source_reconciliation:
            summary_counters["multiple_source_reconciliation"] += 1
            conflict_type_patients[LOW_PRIORITY_CONFLICT_TYPE].add(minimal_barcode)

        summary_counters[f"inclusion_status::{os_endpoint_inclusion_status}"] += 1
        summary_counters["followup_max_last_contact_exceeds_header"] += int(os_followup_max_last_contact_exceeds_header)
        summary_counters["multiple_followup_versions_contributed"] += int(os_multiple_followup_versions_contributed)
        summary_counters["patients_with_any_followup_version_time_source"] += int(bool(followup_versions_contributing_time))

        time_source_key = "|".join(time_source_labels) if time_source_labels else "<missing>"
        time_source_counter[time_source_key] += 1
        time_type_counter[time_type] += 1
        event_source_key = "|".join(event_source_labels) if event_source_labels else "<missing>"
        event_source_counter[event_source_key] += 1
        inclusion_counter[os_endpoint_inclusion_status] += 1

        evidence_summary = json.dumps(
            {
                "header_vital_status": header_status,
                "followup_vital_status_values": followup_vital_values,
                "header_days_to_death": header_days_to_death_raw,
                "followup_max_days_to_death": str(endpoint_row["xml_followup_max_days_to_death"] or "").strip(),
                "header_days_to_last_followup": header_days_to_last_followup_raw,
                "header_days_to_last_known_alive": header_days_to_last_known_alive_raw,
                "followup_max_days_to_last_followup": str(
                    endpoint_row["xml_followup_max_days_to_last_followup"] or ""
                ).strip(),
                "followup_max_days_to_last_known_alive": str(
                    endpoint_row["xml_followup_max_days_to_last_known_alive"] or ""
                ).strip(),
                "followup_versions_present": followup_versions_present,
                "os_event_source": os_row["os_event_source"],
                "os_time_source": os_row["os_time_source"],
                "os_last_contact_source_max": os_row["os_last_contact_source_max"],
            },
            ensure_ascii=True,
            sort_keys=True,
        )

        patient_conflict_types: set[str] = set()

        def add_conflict_row(conflict_type: str, recommended_handling: str, manual_review_priority: str) -> None:
            conflict_rows.append(
                {
                    "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
                    "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
                    "baseline_model_input_v1_run_id": str(
                        workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
                    ),
                    "bcr_patient_barcode": minimal_barcode,
                    "bcr_patient_uuid": minimal_uuid,
                    "provisional_patient_row_id": provisional_patient_row_id,
                    "baseline_analysis_v1_row_id": os_row["baseline_analysis_v1_row_id"],
                    "feature_set_v1_row_index": os_row["feature_set_v1_row_index"],
                    "conflict_type": conflict_type,
                    "source_evidence_summary": evidence_summary,
                    "provisional_os_event": os_event,
                    "provisional_os_time_days": os_time_days,
                    "os_event_source": os_row["os_event_source"],
                    "os_time_source": os_row["os_time_source"],
                    "os_last_contact_source_max": os_row["os_last_contact_source_max"],
                    "os_endpoint_inclusion_status": os_endpoint_inclusion_status,
                    "recommended_handling": recommended_handling,
                    "manual_review_priority": manual_review_priority,
                    "os_requires_manual_review": os_row["os_requires_manual_review"],
                }
            )
            patient_conflict_types.add(conflict_type)

        if os_conflict_vital_status:
            add_conflict_row(
                "vital_status_disagreement",
                "Retain provisional OS values, keep conflict flags explicit, and manually reconcile Alive vs Dead evidence.",
                "high",
            )
        if os_conflict_dead_without_death_time:
            add_conflict_row(
                "dead_without_death_time",
                "Retain provisional event=1 with last-contact fallback time, and manually verify whether a death-time source exists.",
                "high",
            )
        if os_conflict_time_without_clear_status:
            add_conflict_row(
                "time_without_clear_status",
                "Keep the patient row but exclude it from the usable OS subset until status evidence is reconciled.",
                "high",
            )
        if os_time_missing:
            add_conflict_row(
                "no_usable_time",
                "Keep the patient row but exclude it from the usable OS subset until a usable OS time is found.",
                "high",
            )
        if multiple_source_reconciliation:
            add_conflict_row(
                LOW_PRIORITY_CONFLICT_TYPE,
                "Retain the row as rule-based clean; this audit row documents concordant multi-source reconciliation only.",
                "low",
            )

        patient_conflict_audit_coverage[minimal_barcode] = patient_conflict_types

    metrics = {
        "summary_counters": summary_counters,
        "time_source_counter": time_source_counter,
        "time_type_counter": time_type_counter,
        "event_source_counter": event_source_counter,
        "inclusion_counter": inclusion_counter,
        "conflict_type_patients": conflict_type_patients,
        "manual_review_patients": manual_review_patients,
        "patient_conflict_audit_coverage": patient_conflict_audit_coverage,
    }
    return os_rows, conflict_rows, metrics


def build_summary_rows(
    *,
    os_endpoint_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
    os_rows: list[dict[str, str]],
    conflict_rows: list[dict[str, str]],
    rule_spec_rows: list[dict[str, str]],
    metrics: dict[str, Any],
) -> list[dict[str, str]]:
    summary_rows: list[dict[str, str]] = []
    summary_counters: Counter[str] = metrics["summary_counters"]
    time_source_counter: Counter[str] = metrics["time_source_counter"]
    inclusion_counter: Counter[str] = metrics["inclusion_counter"]
    conflict_type_patients: dict[str, set[str]] = metrics["conflict_type_patients"]

    def add_row(summary_section: str, summary_metric: str, summary_value: Any, notes: str = "") -> None:
        summary_rows.append(
            {
                "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
                "summary_section": summary_section,
                "summary_metric": summary_metric,
                "summary_value": str(summary_value),
                "notes": notes,
            }
        )

    add_row(
        "input_runs",
        "endpoint_target_prep_v1_run_id",
        workflow_inputs.endpoint_target_prep_latest_pointer["endpoint_target_prep_v1_run_id"],
    )
    add_row("input_runs", "cohort_v1_build_id", workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"])
    add_row(
        "input_runs",
        "baseline_model_input_v1_run_id",
        workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"],
    )
    add_row("input_runs", "source_run_id", workflow_inputs.endpoint_target_prep_latest_pointer["source_run_id"])

    add_row("row_counts", "os_endpoint_v1_row_count", len(os_rows))
    add_row("row_counts", "os_endpoint_v1_conflict_audit_row_count", len(conflict_rows))
    add_row("row_counts", "os_endpoint_v1_rule_spec_row_count", len(rule_spec_rows))
    add_row("row_counts", "clinical_xml_patient_endpoint_fields_row_count", len(workflow_inputs.patient_fields_rows))
    add_row("row_counts", "clinical_xml_followup_fields_long_row_count", len(workflow_inputs.followup_long_rows))

    add_row("event_counts", "patients_with_os_event_1", summary_counters["os_event_1"])
    add_row("event_counts", "patients_with_os_event_0", summary_counters["os_event_0"])
    add_row("event_counts", "patients_with_missing_os_event", summary_counters["os_event_missing"])
    add_row("time_counts", "patients_with_nonmissing_os_time_days", summary_counters["os_time_nonmissing"])
    add_row("time_counts", "patients_with_missing_os_time_days", summary_counters["os_time_missing"])
    add_row("time_counts", "patients_using_days_to_death", summary_counters["time_using_days_to_death"])
    add_row(
        "time_counts",
        "patients_using_censored_last_contact_time",
        summary_counters["time_using_last_contact"],
    )

    for time_source, count in sorted(time_source_counter.items()):
        add_row(
            "time_source_distribution",
            f"patients_with_time_source__{metric_slug(time_source)}",
            count,
            notes=time_source,
        )

    add_row(
        "conflict_counts",
        "patients_with_vital_status_conflict",
        len(conflict_type_patients["vital_status_disagreement"]),
    )
    add_row(
        "conflict_counts",
        "patients_with_dead_without_death_time",
        len(conflict_type_patients["dead_without_death_time"]),
        notes="Event=1 retained, but os_time_days fell back to last-contact-like evidence because no days_to_death was available.",
    )
    add_row(
        "conflict_counts",
        "patients_with_time_without_clear_status",
        len(conflict_type_patients["time_without_clear_status"]),
    )
    add_row(
        "conflict_counts",
        "patients_with_no_usable_time",
        len(conflict_type_patients["no_usable_time"]),
    )
    add_row(
        "conflict_counts",
        "patients_with_event_time_inconsistent",
        summary_counters["os_conflict_event_time_inconsistent"],
    )
    add_row(
        "conflict_counts",
        "patients_requiring_manual_review",
        summary_counters["os_requires_manual_review"],
    )
    add_row(
        "conflict_counts",
        "patients_with_multiple_source_reconciliation",
        len(conflict_type_patients[LOW_PRIORITY_CONFLICT_TYPE]),
        notes="Low-priority concordant multi-source reconciliation cases kept for audit visibility.",
    )

    for inclusion_status in [
        INCLUSION_STATUS_CLEAN,
        INCLUSION_STATUS_CONFLICT,
        INCLUSION_STATUS_EXCLUDED_NO_TIME,
        INCLUSION_STATUS_EXCLUDED_NO_STATUS,
    ]:
        add_row("inclusion_status", inclusion_status, inclusion_counter[inclusion_status])

    add_row(
        "provenance",
        "patients_with_followup_max_last_contact_exceeding_header",
        summary_counters["followup_max_last_contact_exceeds_header"],
    )
    add_row(
        "provenance",
        "patients_with_multiple_followup_versions_contributing_endpoint_evidence",
        summary_counters["multiple_followup_versions_contributed"],
    )
    add_row(
        "provenance",
        "patients_with_any_followup_version_contributing_chosen_time",
        summary_counters["patients_with_any_followup_version_time_source"],
    )

    add_row("status", "os_endpoint_v1_readiness_interpretation", ENDPOINT_READINESS_INTERPRETATION)
    add_row("status", "stronger_endpoint_claims_status", STRONGER_ENDPOINT_CLAIMS_STATUS)
    add_row("status", "treatment_feasibility_status", TREATMENT_FEASIBILITY_STATUS)
    add_row(
        "status",
        "workflow_scope",
        "os_endpoint_v1_only",
        notes="This workflow does not construct recurrence, DFS, PFS, treatment, or METABRIC-derived endpoints.",
    )

    return summary_rows


def validate_generated_outputs(
    *,
    workflow_inputs: WorkflowInputs,
    os_rows: list[dict[str, str]],
    conflict_rows: list[dict[str, str]],
    rule_spec_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    metrics: dict[str, Any],
) -> dict[str, Any]:
    summary_lookup = {
        str(row["summary_metric"]): str(row["summary_value"])
        for row in summary_rows
    }
    manual_review_patients = metrics["manual_review_patients"]
    conflict_audit_patients = {row["bcr_patient_barcode"] for row in conflict_rows}
    unique_barcodes = {row["bcr_patient_barcode"] for row in os_rows}
    unique_row_ids = {row["provisional_patient_row_id"] for row in os_rows}
    unique_feature_set_indices = {row["feature_set_v1_row_index"] for row in os_rows}
    allowed_inclusion_statuses = {
        INCLUSION_STATUS_CLEAN,
        INCLUSION_STATUS_CONFLICT,
        INCLUSION_STATUS_EXCLUDED_NO_TIME,
        INCLUSION_STATUS_EXCLUDED_NO_STATUS,
    }
    followup_versions_valid = True
    for row in os_rows:
        versions = parse_json_list_cell(
            str(row["os_followup_versions_present_json"] or ""),
            "os_followup_versions_present_json",
        )
        if any(str(version) not in SUPPORTED_FOLLOWUP_VERSIONS for version in versions):
            followup_versions_valid = False
            break

    validation = {
        "passed": True,
        "required_upstream_pointers_found": True,
        "endpoint_target_prep_latest_pointer_found": True,
        "endpoint_target_prep_run_log_completed": True,
        "endpoint_target_prep_validation_passed": True,
        "minimal_cohort_v1_latest_pointer_found": True,
        "minimal_cohort_v1_run_log_completed": True,
        "minimal_cohort_v1_validation_passed": True,
        "baseline_model_input_v1_latest_pointer_found": True,
        "baseline_model_input_v1_run_log_completed": True,
        "baseline_model_input_v1_validation_passed": True,
        "required_source_tables_found": all(path.exists() for path in workflow_inputs.input_paths.values()),
        "output_rows_positive": len(os_rows) > 0,
        "os_row_count_matches_minimal_cohort": len(os_rows) == len(workflow_inputs.minimal_cohort_rows),
        "os_row_count_matches_baseline_model_input": len(os_rows) == len(workflow_inputs.baseline_model_input_rows),
        "patient_ids_unique": len(unique_barcodes) == len(os_rows) and len(unique_row_ids) == len(os_rows),
        "feature_set_row_indices_unique": len(unique_feature_set_indices) == len(os_rows),
        "manual_review_patients_covered_by_conflict_audit": manual_review_patients.issubset(conflict_audit_patients),
        "rule_spec_row_count_positive": len(rule_spec_rows) > 0,
        "conflict_audit_row_count_positive": len(conflict_rows) > 0,
        "summary_row_count_positive": len(summary_rows) > 0,
        "summary_counts_reconcile_to_os_table": (
            summary_lookup.get("os_endpoint_v1_row_count") == str(len(os_rows))
            and summary_lookup.get("patients_with_os_event_1")
            == str(sum(row["os_event"] == "1" for row in os_rows))
            and summary_lookup.get("patients_with_os_event_0")
            == str(sum(row["os_event"] == "0" for row in os_rows))
            and summary_lookup.get("patients_with_nonmissing_os_time_days")
            == str(sum(bool(row["os_time_days"]) for row in os_rows))
            and summary_lookup.get("patients_requiring_manual_review")
            == str(sum(row["os_requires_manual_review"] == "yes" for row in os_rows))
        ),
        "os_inclusion_status_values_valid": all(
            row["os_endpoint_inclusion_status"] in allowed_inclusion_statuses for row in os_rows
        ),
        "followup_versions_restricted_to_supported_set": followup_versions_valid,
        "summary_readiness_matches_expected": (
            summary_lookup.get("os_endpoint_v1_readiness_interpretation") == ENDPOINT_READINESS_INTERPRETATION
        ),
        "summary_treatment_unchanged": (
            summary_lookup.get("treatment_feasibility_status") == TREATMENT_FEASIBILITY_STATUS
        ),
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
    }
    validation["passed"] = all(bool(value) for key, value in validation.items() if key != "passed")
    return validation


def build_latest_pointer_payload(
    *,
    os_endpoint_v1_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
        "endpoint_target_prep_v1_run_id": str(
            workflow_inputs.endpoint_target_prep_latest_pointer["endpoint_target_prep_v1_run_id"]
        ),
        "source_run_id": str(workflow_inputs.endpoint_target_prep_latest_pointer["source_run_id"]),
        "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
        "baseline_model_input_v1_run_id": str(
            workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
        ),
        "baseline_analysis_v1_run_id": str(workflow_inputs.baseline_model_input_run_log["baseline_analysis_v1_run_id"]),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "os_endpoint_v1_tsv": repo_relative(output_paths["os_endpoint_v1_tsv"], paths.repo_root),
        "os_endpoint_v1_rule_spec_tsv": repo_relative(output_paths["os_endpoint_v1_rule_spec_tsv"], paths.repo_root),
        "os_endpoint_v1_conflict_audit_tsv": repo_relative(
            output_paths["os_endpoint_v1_conflict_audit_tsv"],
            paths.repo_root,
        ),
        "os_endpoint_v1_summary_tsv": repo_relative(output_paths["os_endpoint_v1_summary_tsv"], paths.repo_root),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "endpoint_target_prep_v1_latest_json": repo_relative(paths.endpoint_target_prep_latest_pointer, paths.repo_root),
        "minimal_cohort_v1_latest_json": repo_relative(paths.minimal_cohort_latest_pointer, paths.repo_root),
        "baseline_model_input_v1_latest_json": repo_relative(
            paths.baseline_model_input_latest_pointer,
            paths.repo_root,
        ),
    }


def build_upstream_snapshots(workflow_inputs: WorkflowInputs) -> dict[str, Any]:
    return {
        "endpoint_target_prep_v1_latest_pointer": workflow_inputs.endpoint_target_prep_latest_pointer,
        "endpoint_target_prep_v1_run_log_status": {
            "status": workflow_inputs.endpoint_target_prep_run_log.get("status"),
            "endpoint_target_prep_v1_run_id": workflow_inputs.endpoint_target_prep_run_log.get(
                "endpoint_target_prep_v1_run_id"
            ),
            "validation": workflow_inputs.endpoint_target_prep_run_log.get("validation"),
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


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    os_endpoint_v1_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    processed_run_dir = paths.processed_runs_root / os_endpoint_v1_run_id
    audit_run_dir = paths.audit_runs_root / os_endpoint_v1_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)
        validate_upstream_state(workflow_inputs)

        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        os_endpoint_path = processed_run_dir / "os_endpoint_v1.tsv"
        rule_spec_path = audit_run_dir / "os_endpoint_v1_rule_spec.tsv"
        conflict_audit_path = audit_run_dir / "os_endpoint_v1_conflict_audit.tsv"
        summary_path = audit_run_dir / "os_endpoint_v1_summary.tsv"

        rule_spec_rows = build_rule_spec_rows(os_endpoint_v1_run_id)
        os_rows, conflict_rows, metrics = build_os_outputs(
            os_endpoint_v1_run_id=os_endpoint_v1_run_id,
            workflow_inputs=workflow_inputs,
        )
        summary_rows = build_summary_rows(
            os_endpoint_v1_run_id=os_endpoint_v1_run_id,
            workflow_inputs=workflow_inputs,
            os_rows=os_rows,
            conflict_rows=conflict_rows,
            rule_spec_rows=rule_spec_rows,
            metrics=metrics,
        )
        validation_payload = validate_generated_outputs(
            workflow_inputs=workflow_inputs,
            os_rows=os_rows,
            conflict_rows=conflict_rows,
            rule_spec_rows=rule_spec_rows,
            summary_rows=summary_rows,
            metrics=metrics,
        )
        if not validation_payload["passed"]:
            raise OSEndpointV1Error("OS endpoint v1 validation did not pass.")

        helper_module.write_dict_rows_tsv(os_endpoint_path, OS_ENDPOINT_FIELDNAMES, os_rows)
        helper_module.write_dict_rows_tsv(rule_spec_path, RULE_SPEC_FIELDNAMES, rule_spec_rows)
        helper_module.write_dict_rows_tsv(conflict_audit_path, CONFLICT_AUDIT_FIELDNAMES, conflict_rows)
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "os_endpoint_v1_tsv": os_endpoint_path,
            "os_endpoint_v1_rule_spec_tsv": rule_spec_path,
            "os_endpoint_v1_conflict_audit_tsv": conflict_audit_path,
            "os_endpoint_v1_summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            os_endpoint_v1_run_id=os_endpoint_v1_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "endpoint_target_prep_v1_run_id": str(
                workflow_inputs.endpoint_target_prep_latest_pointer["endpoint_target_prep_v1_run_id"]
            ),
            "source_run_id": str(workflow_inputs.endpoint_target_prep_latest_pointer["source_run_id"]),
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
                "endpoint_target_prep_v1_latest_json": repo_relative(
                    paths.endpoint_target_prep_latest_pointer,
                    paths.repo_root,
                ),
                "minimal_cohort_v1_latest_json": repo_relative(paths.minimal_cohort_latest_pointer, paths.repo_root),
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
                "os_endpoint_v1_tsv": repo_relative(os_endpoint_path, paths.repo_root),
                "os_endpoint_v1_rule_spec_tsv": repo_relative(rule_spec_path, paths.repo_root),
                "os_endpoint_v1_conflict_audit_tsv": repo_relative(conflict_audit_path, paths.repo_root),
                "os_endpoint_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation_payload,
            "rules": {
                "unit_of_analysis": "patient/case",
                "endpoint_scope": "overall_survival_v1_only",
                "raw_input_files_immutable": True,
                "parsed_endpoint_prep_outputs_immutable": True,
                "event_rule_dead_if_any_dead_or_days_to_death": True,
                "event_rule_alive_only_if_all_vital_status_alive_and_no_days_to_death": True,
                "time_rule_max_days_to_death_then_max_last_contact_like": True,
                "last_contact_like_fields_json": json_list(
                    [
                        "patient_header_days_to_last_followup",
                        "patient_header_days_to_last_known_alive",
                        "followup_v1_5_days_to_last_followup",
                        "followup_v2_1_days_to_last_followup",
                        "followup_v4_0_days_to_last_followup",
                        "followup_v1_5_days_to_last_known_alive",
                        "followup_v2_1_days_to_last_known_alive",
                        "followup_v4_0_days_to_last_known_alive",
                    ]
                ),
                "days_to_last_followup_primary": True,
                "days_to_last_known_alive_auxiliary": True,
                "include_plus_flag_conflict_policy": True,
                "manual_review_does_not_drop_rule_resolved_rows": True,
                "no_progression_or_recurrence_fields": True,
                "no_treatment_fields": True,
                "no_modeling": True,
                "no_metabric": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "missing_like_tokens_json": json_list(sorted(MISSING_LIKE_TOKENS)),
                "supported_followup_versions_json": json_list(list(SUPPORTED_FOLLOWUP_VERSIONS)),
            },
            "counts": {
                "clinical_xml_patient_endpoint_fields_row_count": len(workflow_inputs.patient_fields_rows),
                "clinical_xml_followup_fields_long_row_count": len(workflow_inputs.followup_long_rows),
                "endpoint_target_prep_v1_row_count": len(workflow_inputs.endpoint_prep_rows),
                "minimal_cohort_v1_row_count": len(workflow_inputs.minimal_cohort_rows),
                "baseline_model_input_v1_row_count": len(workflow_inputs.baseline_model_input_rows),
                "baseline_feature_set_audit_map_row_count": len(workflow_inputs.baseline_feature_set_audit_map_rows),
                "os_endpoint_v1_row_count": len(os_rows),
                "os_endpoint_v1_rule_spec_row_count": len(rule_spec_rows),
                "os_endpoint_v1_conflict_audit_row_count": len(conflict_rows),
                "os_endpoint_v1_summary_row_count": len(summary_rows),
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
            "os_endpoint_v1_run_id": os_endpoint_v1_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_os_endpoint_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA OS endpoint v1 workflow complete.")
    print(f"OS endpoint v1 run ID: {run_log['os_endpoint_v1_run_id']}")
    print(f"Endpoint-target prep v1 run ID: {run_log['endpoint_target_prep_v1_run_id']}")
    print(f"Cohort v1 build ID: {run_log['cohort_v1_build_id']}")
    print(f"Baseline model-input v1 run ID: {run_log['baseline_model_input_v1_run_id']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"OS endpoint TSV: {run_log['outputs']['os_endpoint_v1_tsv']}")
    print(f"Rule spec TSV: {run_log['outputs']['os_endpoint_v1_rule_spec_tsv']}")
    print(f"Conflict audit TSV: {run_log['outputs']['os_endpoint_v1_conflict_audit_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['os_endpoint_v1_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
