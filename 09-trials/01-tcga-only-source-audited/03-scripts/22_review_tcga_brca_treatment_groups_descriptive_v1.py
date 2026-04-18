#!/usr/bin/env python
"""Review provisional TCGA-BRCA treatment groups with descriptive summaries only."""

from __future__ import annotations

import importlib.util
import json
import statistics
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

APPROVED_TREATMENT_GROUP_ORDER = [
    "no_drug_record",
    "single_chemotherapy",
    "single_hormone_therapy",
    "single_targeted_therapy",
    "single_immunotherapy",
    "single_ancillary_or_other",
    "mixed_multi_type",
    "missing_type_only",
]
REQUIRED_CONFOUNDER_FIELDS = [
    "age_at_diagnosis",
    "er_status_by_ihc",
    "pr_status_by_ihc",
    "her2_status_by_ihc",
    "ajcc_pathologic_tumor_stage",
    "histological_type",
]
COVERAGE_READY_FIELDS = [
    "er_status_by_ihc",
    "pr_status_by_ihc",
    "her2_status_by_ihc",
    "ajcc_pathologic_tumor_stage",
]
CATEGORICAL_COMPOSITION_FIELDS = [
    "er_status_by_ihc",
    "pr_status_by_ihc",
    "her2_status_by_ihc",
    "ajcc_pathologic_tumor_stage",
    "histological_type",
]
GROUP_SUMMARY_FIELDNAMES = [
    "treatment_groups_descriptive_v1_run_id",
    "group_order",
    "treatment_group_v1",
    "patient_count",
    "patient_fraction",
    "os_event_count",
    "os_event_fraction",
    "manual_review_count",
    "manual_review_fraction",
    "radiation_overlap_count",
    "radiation_overlap_fraction",
    "regimen_context_count",
    "regimen_context_fraction",
    "timing_coverage_count",
    "timing_coverage_fraction",
    "notes",
]
CONFOUNDER_COVERAGE_FIELDNAMES = [
    "treatment_groups_descriptive_v1_run_id",
    "treatment_group_v1",
    "field_name",
    "row_count",
    "non_missing_count",
    "missing_like_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "notes",
]
VALUE_COMPOSITION_FIELDNAMES = [
    "treatment_groups_descriptive_v1_run_id",
    "treatment_group_v1",
    "field_name",
    "value_summary_type",
    "group_patient_count",
    "non_missing_count",
    "missing_like_count",
    "distinct_non_missing_count",
    "value_rank",
    "value_label",
    "value_count",
    "value_fraction",
    "numeric_min",
    "numeric_median",
    "numeric_max",
    "notes",
]
MANUAL_REVIEW_SUMMARY_FIELDNAMES = [
    "treatment_groups_descriptive_v1_run_id",
    "treatment_group_v1",
    "summary_level",
    "group_patient_count",
    "manual_review_count",
    "manual_review_fraction",
    "conflict_type",
    "conflict_type_count",
    "conflict_type_fraction_within_manual_review",
    "notes",
]
SUMMARY_FIELDNAMES = [
    "treatment_groups_descriptive_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
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
MIN_GROUP_SIZE_FOR_FREEZE_REVIEW = 50
MANUAL_REVIEW_FRACTION_TOLERANCE = 0.20
CONFOUNDER_NON_MISSING_TOLERANCE = 0.85
MAX_TOP_VALUES = 10
READINESS_NOT_READY = "not_ready"
READINESS_DESCRIPTIVE_COMPLETED = "ready_for_grouped_descriptive_review_completed"
READINESS_TREATMENT_ARM_FREEZE_REVIEW_NEXT = "ready_for_treatment_arm_freeze_review_next"
TREATMENT_RECOMMENDATION_MODELING_STATUS = "out_of_scope"


class TreatmentGroupsDescriptiveV1Error(RuntimeError):
    """Raised when the grouped descriptive review workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the grouped descriptive review workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    grouping_latest_pointer: Path
    os_endpoint_latest_pointer: Path
    baseline_model_input_latest_pointer: Path
    baseline_feature_set_latest_pointer: Path
    baseline_analysis_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    grouping_latest_pointer: dict[str, Any]
    grouping_run_log: dict[str, Any]
    os_endpoint_latest_pointer: dict[str, Any]
    os_endpoint_run_log: dict[str, Any]
    baseline_model_input_latest_pointer: dict[str, Any]
    baseline_model_input_run_log: dict[str, Any]
    baseline_feature_set_latest_pointer: dict[str, Any]
    baseline_feature_set_run_log: dict[str, Any]
    baseline_analysis_latest_pointer: dict[str, Any]
    baseline_analysis_run_log: dict[str, Any]
    grouping_rows: list[dict[str, str]]
    grouping_conflict_rows: list[dict[str, str]]
    os_endpoint_rows: list[dict[str, str]]
    baseline_feature_set_rows: list[dict[str, str]]
    baseline_feature_set_audit_map_rows: list[dict[str, str]]
    baseline_analysis_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class GroupMetrics:
    """Per-group and cohort-level descriptive metrics."""

    total_cohort_size: int
    group_counts: dict[str, int]
    os_event_counts: dict[str, int]
    manual_review_counts: dict[str, int]
    radiation_overlap_counts: dict[str, int]
    regimen_context_counts: dict[str, int]
    timing_coverage_counts: dict[str, int]
    present_groups: list[str]
    large_groups: list[str]
    groups_manual_review_within_tolerance: list[str]
    groups_all_ready_fields_within_tolerance: list[str]
    coverage_non_missing_counts: dict[tuple[str, str], int]
    coverage_missing_like_counts: dict[tuple[str, str], int]
    coverage_distinct_non_missing_counts: dict[tuple[str, str], int]
    coverage_non_missing_fractions: dict[tuple[str, str], float]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise TreatmentGroupsDescriptiveV1Error(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise TreatmentGroupsDescriptiveV1Error(f"Unable to create an import spec for: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def format_fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.0000"
    return f"{numerator / denominator:.4f}"


def format_float(value: float) -> str:
    return f"{value:.4f}"


def format_numeric(value: float | int | None) -> str:
    if value is None:
        return ""
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.1f}".rstrip("0").rstrip(".")


def normalize_barcode(value: str) -> str:
    return value.strip().upper()


def normalize_uuid(value: str) -> str:
    return value.strip().upper()


def normalize_missing_like(value: Any) -> str:
    return str(value or "").strip().lower()


def is_missing_like(value: Any) -> bool:
    return normalize_missing_like(value) in MISSING_LIKE_TOKENS


def parse_int(value: Any, label: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise TreatmentGroupsDescriptiveV1Error(
            f"Expected integer-like value for {label}: {value!r}"
        ) from exc


def parse_float_or_none(value: Any) -> float | None:
    stripped = str(value or "").strip()
    if not stripped or is_missing_like(stripped):
        return None
    try:
        return float(stripped)
    except ValueError:
        return None


def require_columns(rows: list[dict[str, str]], required_columns: set[str], label: str) -> None:
    if not rows:
        raise TreatmentGroupsDescriptiveV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise TreatmentGroupsDescriptiveV1Error(f"{label} is missing required columns: {sorted(missing)}")


def require_completed_run_log(
    *,
    repo_root: Path,
    run_log_relative_path: str,
    label: str,
    helper_module: Any,
) -> dict[str, Any]:
    run_log_path = helper_module.resolve_existing_path(repo_root, run_log_relative_path, label)
    run_log = helper_module.load_json(run_log_path)
    if run_log.get("status") != "completed":
        raise TreatmentGroupsDescriptiveV1Error(f"{label} does not report status == completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise TreatmentGroupsDescriptiveV1Error(f"{label} does not report validation.passed == true.")
    return run_log


def build_unique_lookup(
    rows: list[dict[str, str]],
    *,
    key_field: str,
    label: str,
    normalize_key: bool = False,
) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        raw_key = str(row.get(key_field, ""))
        key = normalize_barcode(raw_key) if normalize_key else raw_key.strip()
        if not key:
            raise TreatmentGroupsDescriptiveV1Error(f"Blank {key_field} encountered in {label}.")
        if key in lookup:
            raise TreatmentGroupsDescriptiveV1Error(
                f"Duplicate {key_field} '{key}' encountered in {label}."
            )
        lookup[key] = row
    return lookup


def has_any_regimen_context(row: dict[str, str]) -> bool:
    return str(row.get("regimen_context_values_json", "")).strip() not in {"", "[]"}


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    if not trial_config.exists():
        raise TreatmentGroupsDescriptiveV1Error(f"Required trial config not found: {trial_config}")

    trial_config_data = helper_module.load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    treatment_prep_root = audit_root / "tcga-brca" / "treatment-prep"
    analysis_prep_root = audit_root / "tcga-brca" / "analysis-prep"
    model_input_root = audit_root / "tcga-brca" / "model-input"

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        audit_runs_root=treatment_prep_root / "treatment_groups_descriptive_v1_runs",
        latest_pointer=treatment_prep_root / "tcga_brca_treatment_groups_descriptive_v1_latest.json",
        grouping_latest_pointer=treatment_prep_root / "tcga_brca_patient_treatment_grouping_v1_latest.json",
        os_endpoint_latest_pointer=audit_root
        / "tcga-brca"
        / "endpoint-prep"
        / "tcga_brca_os_endpoint_v1_latest.json",
        baseline_model_input_latest_pointer=model_input_root / "tcga_brca_baseline_model_input_v1_latest.json",
        baseline_feature_set_latest_pointer=analysis_prep_root
        / "tcga_brca_baseline_feature_set_v1_latest.json",
        baseline_analysis_latest_pointer=analysis_prep_root / "tcga_brca_baseline_analysis_v1_latest.json",
    )


def sort_baseline_analysis_rows(baseline_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    sorted_rows = sorted(
        baseline_rows,
        key=lambda row: parse_int(row["baseline_analysis_v1_row_id"], "baseline_analysis_v1_row_id"),
    )
    observed_ids = [
        parse_int(row["baseline_analysis_v1_row_id"], "baseline_analysis_v1_row_id")
        for row in sorted_rows
    ]
    if len(set(observed_ids)) != len(observed_ids):
        raise TreatmentGroupsDescriptiveV1Error(
            "baseline_analysis_v1.tsv does not preserve one row per baseline_analysis_v1_row_id."
        )
    if observed_ids != list(range(1, len(sorted_rows) + 1)):
        raise TreatmentGroupsDescriptiveV1Error(
            "baseline_analysis_v1.tsv does not preserve sequential baseline_analysis_v1_row_id values."
        )
    return sorted_rows


def validate_feature_set_lineage(
    baseline_analysis_rows: list[dict[str, str]],
    baseline_feature_set_rows: list[dict[str, str]],
    audit_map_rows: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    sorted_analysis_rows = sort_baseline_analysis_rows(baseline_analysis_rows)
    if len(sorted_analysis_rows) != len(baseline_feature_set_rows):
        raise TreatmentGroupsDescriptiveV1Error(
            "baseline_feature_set_v1.tsv row count does not match baseline_analysis_v1.tsv row count."
        )
    if len(sorted_analysis_rows) != len(audit_map_rows):
        raise TreatmentGroupsDescriptiveV1Error(
            "baseline_feature_set_v1_audit_map.tsv row count does not match baseline_analysis_v1.tsv row count."
        )

    confounder_lookup: dict[str, dict[str, str]] = {}
    for expected_index, (analysis_row, feature_row, audit_row) in enumerate(
        zip(sorted_analysis_rows, baseline_feature_set_rows, audit_map_rows),
        start=1,
    ):
        observed_index = parse_int(audit_row["feature_set_v1_row_index"], "feature_set_v1_row_index")
        if observed_index != expected_index:
            raise TreatmentGroupsDescriptiveV1Error(
                "baseline_feature_set_v1_audit_map.tsv does not preserve sequential "
                "feature_set_v1_row_index values."
            )
        for id_field in [
            "baseline_analysis_v1_row_id",
            "provisional_patient_row_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
        ]:
            left_value = str(audit_row[id_field]).strip()
            right_value = str(analysis_row[id_field]).strip()
            if normalize_barcode(left_value) != normalize_barcode(right_value) and id_field == "bcr_patient_barcode":
                raise TreatmentGroupsDescriptiveV1Error(
                    f"Feature-set audit-map barcode mismatch at row index {expected_index}."
                )
            if normalize_uuid(left_value) != normalize_uuid(right_value) and id_field == "bcr_patient_uuid":
                raise TreatmentGroupsDescriptiveV1Error(
                    f"Feature-set audit-map UUID mismatch at row index {expected_index}."
                )
            if id_field not in {"bcr_patient_barcode", "bcr_patient_uuid"} and left_value != right_value:
                raise TreatmentGroupsDescriptiveV1Error(
                    f"Feature-set audit-map identifier mismatch for {id_field} at row index {expected_index}."
                )
        for field_name in REQUIRED_CONFOUNDER_FIELDS:
            if str(feature_row.get(field_name, "")) != str(analysis_row.get(field_name, "")):
                raise TreatmentGroupsDescriptiveV1Error(
                    f"Feature-set value mismatch for {field_name} at feature_set_v1_row_index {expected_index}."
                )

        confounder_lookup[str(expected_index)] = {
            "feature_set_v1_row_index": str(expected_index),
            "baseline_analysis_v1_row_id": str(analysis_row["baseline_analysis_v1_row_id"]),
            "provisional_patient_row_id": str(analysis_row["provisional_patient_row_id"]),
            "bcr_patient_barcode": str(analysis_row["bcr_patient_barcode"]),
            "bcr_patient_uuid": str(analysis_row["bcr_patient_uuid"]),
            **{field_name: str(analysis_row.get(field_name, "")) for field_name in REQUIRED_CONFOUNDER_FIELDS},
        }

    return confounder_lookup


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    required_pointers = [
        (paths.grouping_latest_pointer, "patient treatment grouping v1 latest pointer"),
        (paths.os_endpoint_latest_pointer, "OS endpoint v1 latest pointer"),
        (paths.baseline_model_input_latest_pointer, "baseline model-input v1 latest pointer"),
        (paths.baseline_feature_set_latest_pointer, "baseline feature-set v1 latest pointer"),
        (paths.baseline_analysis_latest_pointer, "baseline analysis v1 latest pointer"),
    ]
    for pointer_path, label in required_pointers:
        if not pointer_path.exists():
            raise TreatmentGroupsDescriptiveV1Error(f"Required {label} not found: {pointer_path}")

    grouping_pointer = helper_module.load_json(paths.grouping_latest_pointer)
    helper_module.require_keys(
        grouping_pointer,
        {
            "patient_treatment_grouping_v1_run_id",
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "patient_treatment_grouping_v1_tsv",
            "patient_treatment_grouping_v1_conflict_audit_tsv",
            "run_log_json",
        },
        "patient treatment grouping v1 latest pointer",
        paths.grouping_latest_pointer,
    )
    os_pointer = helper_module.load_json(paths.os_endpoint_latest_pointer)
    helper_module.require_keys(
        os_pointer,
        {
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "baseline_analysis_v1_run_id",
            "cohort_v1_build_id",
            "os_endpoint_v1_tsv",
            "run_log_json",
        },
        "OS endpoint v1 latest pointer",
        paths.os_endpoint_latest_pointer,
    )
    baseline_model_input_pointer = helper_module.load_json(paths.baseline_model_input_latest_pointer)
    helper_module.require_keys(
        baseline_model_input_pointer,
        {
            "baseline_model_input_v1_run_id",
            "baseline_feature_set_v1_run_id",
            "baseline_feature_set_v1_latest_json",
            "baseline_feature_set_v1_audit_map_tsv",
            "cohort_v1_build_id",
            "run_log_json",
        },
        "baseline model-input v1 latest pointer",
        paths.baseline_model_input_latest_pointer,
    )
    baseline_feature_set_pointer = helper_module.load_json(paths.baseline_feature_set_latest_pointer)
    helper_module.require_keys(
        baseline_feature_set_pointer,
        {
            "baseline_feature_set_v1_run_id",
            "baseline_analysis_v1_run_id",
            "baseline_feature_set_v1_tsv",
            "baseline_feature_set_v1_audit_map_tsv",
            "cohort_v1_build_id",
            "run_log_json",
        },
        "baseline feature-set v1 latest pointer",
        paths.baseline_feature_set_latest_pointer,
    )
    baseline_analysis_pointer = helper_module.load_json(paths.baseline_analysis_latest_pointer)
    helper_module.require_keys(
        baseline_analysis_pointer,
        {
            "baseline_analysis_v1_run_id",
            "baseline_analysis_v1_tsv",
            "cohort_v1_build_id",
            "run_log_json",
        },
        "baseline analysis v1 latest pointer",
        paths.baseline_analysis_latest_pointer,
    )

    grouping_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(grouping_pointer["run_log_json"]),
        label="patient treatment grouping v1 run log",
        helper_module=helper_module,
    )
    os_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(os_pointer["run_log_json"]),
        label="OS endpoint v1 run log",
        helper_module=helper_module,
    )
    baseline_model_input_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(baseline_model_input_pointer["run_log_json"]),
        label="baseline model-input v1 run log",
        helper_module=helper_module,
    )
    baseline_feature_set_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(baseline_feature_set_pointer["run_log_json"]),
        label="baseline feature-set v1 run log",
        helper_module=helper_module,
    )
    baseline_analysis_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(baseline_analysis_pointer["run_log_json"]),
        label="baseline analysis v1 run log",
        helper_module=helper_module,
    )

    if str(grouping_pointer["os_endpoint_v1_run_id"]) != str(os_pointer["os_endpoint_v1_run_id"]):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between treatment grouping and OS endpoint pointers for os_endpoint_v1_run_id."
        )
    if str(grouping_pointer["baseline_model_input_v1_run_id"]) != str(
        os_pointer["baseline_model_input_v1_run_id"]
    ):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between treatment grouping and OS endpoint pointers for baseline_model_input_v1_run_id."
        )
    if str(grouping_pointer["baseline_model_input_v1_run_id"]) != str(
        baseline_model_input_pointer["baseline_model_input_v1_run_id"]
    ):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between treatment grouping and baseline model-input pointers for baseline_model_input_v1_run_id."
        )
    if str(grouping_pointer["cohort_v1_build_id"]) != str(os_pointer["cohort_v1_build_id"]):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between treatment grouping and OS endpoint pointers for cohort_v1_build_id."
        )
    if str(grouping_pointer["cohort_v1_build_id"]) != str(
        baseline_model_input_pointer["cohort_v1_build_id"]
    ):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between treatment grouping and baseline model-input pointers for cohort_v1_build_id."
        )
    if str(grouping_pointer["cohort_v1_build_id"]) != str(
        baseline_feature_set_pointer["cohort_v1_build_id"]
    ):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between treatment grouping and baseline feature-set pointers for cohort_v1_build_id."
        )
    if str(grouping_pointer["cohort_v1_build_id"]) != str(baseline_analysis_pointer["cohort_v1_build_id"]):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between treatment grouping and baseline analysis pointers for cohort_v1_build_id."
        )
    if str(os_pointer["baseline_analysis_v1_run_id"]) != str(
        baseline_analysis_pointer["baseline_analysis_v1_run_id"]
    ):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between OS endpoint and baseline analysis pointers for baseline_analysis_v1_run_id."
        )
    if str(baseline_model_input_pointer["baseline_feature_set_v1_run_id"]) != str(
        baseline_feature_set_pointer["baseline_feature_set_v1_run_id"]
    ):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between baseline model-input and baseline feature-set pointers for baseline_feature_set_v1_run_id."
        )
    if str(baseline_feature_set_pointer["baseline_analysis_v1_run_id"]) != str(
        baseline_analysis_pointer["baseline_analysis_v1_run_id"]
    ):
        raise TreatmentGroupsDescriptiveV1Error(
            "Mismatch between baseline feature-set and baseline analysis pointers for baseline_analysis_v1_run_id."
        )

    grouping_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(grouping_pointer["patient_treatment_grouping_v1_tsv"]),
        "patient_treatment_grouping_v1.tsv",
    )
    grouping_conflict_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(grouping_pointer["patient_treatment_grouping_v1_conflict_audit_tsv"]),
        "patient_treatment_grouping_v1_conflict_audit.tsv",
    )
    os_endpoint_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(os_pointer["os_endpoint_v1_tsv"]),
        "os_endpoint_v1.tsv",
    )
    baseline_feature_set_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(baseline_feature_set_pointer["baseline_feature_set_v1_tsv"]),
        "baseline_feature_set_v1.tsv",
    )
    baseline_feature_set_audit_map_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(baseline_feature_set_pointer["baseline_feature_set_v1_audit_map_tsv"]),
        "baseline_feature_set_v1_audit_map.tsv",
    )
    baseline_analysis_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(baseline_analysis_pointer["baseline_analysis_v1_tsv"]),
        "baseline_analysis_v1.tsv",
    )

    grouping_rows = helper_module.read_tsv_dict_rows(grouping_tsv)
    grouping_conflict_rows = helper_module.read_tsv_dict_rows(grouping_conflict_tsv)
    os_endpoint_rows = helper_module.read_tsv_dict_rows(os_endpoint_tsv)
    baseline_feature_set_rows = helper_module.read_tsv_dict_rows(baseline_feature_set_tsv)
    baseline_feature_set_audit_map_rows = helper_module.read_tsv_dict_rows(baseline_feature_set_audit_map_tsv)
    baseline_analysis_rows = helper_module.read_tsv_dict_rows(baseline_analysis_tsv)

    require_columns(
        grouping_rows,
        {
            "patient_treatment_grouping_v1_run_id",
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
            "os_event",
            "os_time_days",
            "treatment_group_v1",
            "treatment_group_v1_requires_manual_review",
            "has_any_radiation_row",
            "regimen_context_values_json",
            "has_any_treatment_timing",
        },
        "patient_treatment_grouping_v1.tsv",
    )
    require_columns(
        grouping_conflict_rows,
        {
            "treatment_group_v1",
            "conflict_type",
            "treatment_group_v1_requires_manual_review",
        },
        "patient_treatment_grouping_v1_conflict_audit.tsv",
    )
    require_columns(
        os_endpoint_rows,
        {
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
            "os_event",
            "os_time_days",
        },
        "os_endpoint_v1.tsv",
    )
    require_columns(
        baseline_feature_set_rows,
        set(REQUIRED_CONFOUNDER_FIELDS),
        "baseline_feature_set_v1.tsv",
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
    require_columns(
        baseline_analysis_rows,
        {
            "baseline_analysis_v1_row_id",
            "provisional_patient_row_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            *REQUIRED_CONFOUNDER_FIELDS,
        },
        "baseline_analysis_v1.tsv",
    )

    validate_feature_set_lineage(
        baseline_analysis_rows,
        baseline_feature_set_rows,
        baseline_feature_set_audit_map_rows,
    )

    return WorkflowInputs(
        grouping_latest_pointer=grouping_pointer,
        grouping_run_log=grouping_run_log,
        os_endpoint_latest_pointer=os_pointer,
        os_endpoint_run_log=os_run_log,
        baseline_model_input_latest_pointer=baseline_model_input_pointer,
        baseline_model_input_run_log=baseline_model_input_run_log,
        baseline_feature_set_latest_pointer=baseline_feature_set_pointer,
        baseline_feature_set_run_log=baseline_feature_set_run_log,
        baseline_analysis_latest_pointer=baseline_analysis_pointer,
        baseline_analysis_run_log=baseline_analysis_run_log,
        grouping_rows=grouping_rows,
        grouping_conflict_rows=grouping_conflict_rows,
        os_endpoint_rows=os_endpoint_rows,
        baseline_feature_set_rows=baseline_feature_set_rows,
        baseline_feature_set_audit_map_rows=baseline_feature_set_audit_map_rows,
        baseline_analysis_rows=baseline_analysis_rows,
        input_paths={
            "patient_treatment_grouping_v1_tsv": grouping_tsv,
            "patient_treatment_grouping_v1_conflict_audit_tsv": grouping_conflict_tsv,
            "os_endpoint_v1_tsv": os_endpoint_tsv,
            "baseline_feature_set_v1_tsv": baseline_feature_set_tsv,
            "baseline_feature_set_v1_audit_map_tsv": baseline_feature_set_audit_map_tsv,
            "baseline_analysis_v1_tsv": baseline_analysis_tsv,
        },
    )


def build_joined_rows(
    workflow_inputs: WorkflowInputs,
) -> tuple[list[dict[str, str]], bool, bool, bool, bool]:
    os_lookup = build_unique_lookup(
        workflow_inputs.os_endpoint_rows,
        key_field="bcr_patient_barcode",
        label="os_endpoint_v1.tsv",
        normalize_key=True,
    )
    baseline_analysis_lookup = build_unique_lookup(
        workflow_inputs.baseline_analysis_rows,
        key_field="baseline_analysis_v1_row_id",
        label="baseline_analysis_v1.tsv",
    )
    confounder_lookup = validate_feature_set_lineage(
        workflow_inputs.baseline_analysis_rows,
        workflow_inputs.baseline_feature_set_rows,
        workflow_inputs.baseline_feature_set_audit_map_rows,
    )

    joined_rows: list[dict[str, str]] = []
    grouping_barcodes: list[str] = []
    os_barcodes: list[str] = [normalize_barcode(row["bcr_patient_barcode"]) for row in workflow_inputs.os_endpoint_rows]
    grouping_matches_os = True
    grouping_matches_feature_map = True

    for grouping_row in workflow_inputs.grouping_rows:
        barcode = normalize_barcode(grouping_row["bcr_patient_barcode"])
        grouping_barcodes.append(barcode)
        os_row = os_lookup.get(barcode)
        if os_row is None:
            raise TreatmentGroupsDescriptiveV1Error(
                f"Grouping row barcode not found in os_endpoint_v1.tsv: {grouping_row['bcr_patient_barcode']}"
            )
        for field_name in [
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
            "os_event",
            "os_time_days",
        ]:
            left_value = str(grouping_row[field_name]).strip()
            right_value = str(os_row[field_name]).strip()
            if field_name == "bcr_patient_barcode":
                if normalize_barcode(left_value) != normalize_barcode(right_value):
                    grouping_matches_os = False
            elif field_name == "bcr_patient_uuid":
                if normalize_uuid(left_value) != normalize_uuid(right_value):
                    grouping_matches_os = False
            elif left_value != right_value:
                grouping_matches_os = False

        feature_row = confounder_lookup.get(str(grouping_row["feature_set_v1_row_index"]).strip())
        if feature_row is None:
            raise TreatmentGroupsDescriptiveV1Error(
                "Grouping row feature_set_v1_row_index not found in the validated baseline feature-set audit map: "
                f"{grouping_row['feature_set_v1_row_index']}"
            )
        analysis_row = baseline_analysis_lookup.get(str(grouping_row["baseline_analysis_v1_row_id"]).strip())
        if analysis_row is None:
            raise TreatmentGroupsDescriptiveV1Error(
                "Grouping row baseline_analysis_v1_row_id not found in baseline_analysis_v1.tsv: "
                f"{grouping_row['baseline_analysis_v1_row_id']}"
            )

        if str(feature_row["baseline_analysis_v1_row_id"]).strip() != str(
            grouping_row["baseline_analysis_v1_row_id"]
        ).strip():
            grouping_matches_feature_map = False
        if normalize_barcode(feature_row["bcr_patient_barcode"]) != barcode:
            grouping_matches_feature_map = False
        if normalize_uuid(feature_row["bcr_patient_uuid"]) != normalize_uuid(grouping_row["bcr_patient_uuid"]):
            grouping_matches_feature_map = False
        if str(feature_row["provisional_patient_row_id"]).strip() != str(
            grouping_row["provisional_patient_row_id"]
        ).strip():
            grouping_matches_feature_map = False
        for field_name in REQUIRED_CONFOUNDER_FIELDS:
            if str(feature_row.get(field_name, "")) != str(analysis_row.get(field_name, "")):
                grouping_matches_feature_map = False

        joined_rows.append(
            {
                **grouping_row,
                **{field_name: str(analysis_row.get(field_name, "")) for field_name in REQUIRED_CONFOUNDER_FIELDS},
            }
        )

    row_order_preserved = grouping_barcodes == os_barcodes
    all_os_patients_accounted_for = set(grouping_barcodes) == set(os_barcodes)
    no_duplicate_grouping_barcodes = len(grouping_barcodes) == len(set(grouping_barcodes))
    return (
        joined_rows,
        row_order_preserved,
        all_os_patients_accounted_for,
        no_duplicate_grouping_barcodes,
        grouping_matches_os and grouping_matches_feature_map,
    )


def build_group_metrics(joined_rows: list[dict[str, str]]) -> GroupMetrics:
    group_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    os_event_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    manual_review_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    radiation_overlap_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    regimen_context_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    timing_coverage_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    coverage_non_missing_counts: dict[tuple[str, str], int] = {}
    coverage_missing_like_counts: dict[tuple[str, str], int] = {}
    coverage_distinct_non_missing_counts: dict[tuple[str, str], int] = {}
    coverage_non_missing_fractions: dict[tuple[str, str], float] = {}

    rows_by_group = {group_name: [] for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    for row in joined_rows:
        group_name = str(row["treatment_group_v1"]).strip()
        if group_name not in group_counts:
            raise TreatmentGroupsDescriptiveV1Error(
                f"Unexpected treatment_group_v1 encountered in grouped descriptive review: {group_name}"
            )
        group_counts[group_name] += 1
        rows_by_group[group_name].append(row)
        if str(row["os_event"]).strip() == "1":
            os_event_counts[group_name] += 1
        if str(row["treatment_group_v1_requires_manual_review"]).strip() == "yes":
            manual_review_counts[group_name] += 1
        if str(row["has_any_radiation_row"]).strip() == "yes":
            radiation_overlap_counts[group_name] += 1
        if has_any_regimen_context(row):
            regimen_context_counts[group_name] += 1
        if str(row["has_any_treatment_timing"]).strip() == "yes":
            timing_coverage_counts[group_name] += 1

    present_groups = [group_name for group_name in APPROVED_TREATMENT_GROUP_ORDER if group_counts[group_name] > 0]
    large_groups = [
        group_name for group_name in APPROVED_TREATMENT_GROUP_ORDER if group_counts[group_name] >= MIN_GROUP_SIZE_FOR_FREEZE_REVIEW
    ]

    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        group_rows = rows_by_group[group_name]
        row_count = len(group_rows)
        for field_name in REQUIRED_CONFOUNDER_FIELDS:
            non_missing_values = [
                str(row.get(field_name, ""))
                for row in group_rows
                if not is_missing_like(row.get(field_name, ""))
            ]
            non_missing_count = len(non_missing_values)
            missing_like_count = row_count - non_missing_count
            distinct_non_missing_count = len(set(non_missing_values))
            coverage_non_missing_counts[(group_name, field_name)] = non_missing_count
            coverage_missing_like_counts[(group_name, field_name)] = missing_like_count
            coverage_distinct_non_missing_counts[(group_name, field_name)] = distinct_non_missing_count
            coverage_non_missing_fractions[(group_name, field_name)] = (
                non_missing_count / row_count if row_count > 0 else 0.0
            )

    groups_manual_review_within_tolerance = [
        group_name
        for group_name in present_groups
        if (
            manual_review_counts[group_name] / group_counts[group_name]
            if group_counts[group_name] > 0
            else 0.0
        )
        <= MANUAL_REVIEW_FRACTION_TOLERANCE
    ]
    groups_all_ready_fields_within_tolerance = [
        group_name
        for group_name in present_groups
        if all(
            coverage_non_missing_fractions[(group_name, field_name)] >= CONFOUNDER_NON_MISSING_TOLERANCE
            for field_name in COVERAGE_READY_FIELDS
        )
    ]

    return GroupMetrics(
        total_cohort_size=len(joined_rows),
        group_counts=group_counts,
        os_event_counts=os_event_counts,
        manual_review_counts=manual_review_counts,
        radiation_overlap_counts=radiation_overlap_counts,
        regimen_context_counts=regimen_context_counts,
        timing_coverage_counts=timing_coverage_counts,
        present_groups=present_groups,
        large_groups=large_groups,
        groups_manual_review_within_tolerance=groups_manual_review_within_tolerance,
        groups_all_ready_fields_within_tolerance=groups_all_ready_fields_within_tolerance,
        coverage_non_missing_counts=coverage_non_missing_counts,
        coverage_missing_like_counts=coverage_missing_like_counts,
        coverage_distinct_non_missing_counts=coverage_distinct_non_missing_counts,
        coverage_non_missing_fractions=coverage_non_missing_fractions,
    )


def build_group_summary_rows(
    run_id: str,
    metrics: GroupMetrics,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for group_order, group_name in enumerate(APPROVED_TREATMENT_GROUP_ORDER, start=1):
        patient_count = metrics.group_counts[group_name]
        manual_review_count = metrics.manual_review_counts[group_name]
        notes = ""
        if patient_count == 0:
            notes = "No patients assigned to this provisional treatment group in the current grouping run."
        elif manual_review_count == patient_count:
            notes = "Every patient in this provisional group still requires manual review."
        elif manual_review_count > 0:
            notes = "This provisional group mixes reviewed and manual-review-required patients."
        elif patient_count < MIN_GROUP_SIZE_FOR_FREEZE_REVIEW:
            notes = "Small provisional group; keep interpretation descriptive only."

        rows.append(
            {
                "treatment_groups_descriptive_v1_run_id": run_id,
                "group_order": str(group_order),
                "treatment_group_v1": group_name,
                "patient_count": str(patient_count),
                "patient_fraction": format_fraction(patient_count, metrics.total_cohort_size),
                "os_event_count": str(metrics.os_event_counts[group_name]),
                "os_event_fraction": format_fraction(metrics.os_event_counts[group_name], patient_count),
                "manual_review_count": str(manual_review_count),
                "manual_review_fraction": format_fraction(manual_review_count, patient_count),
                "radiation_overlap_count": str(metrics.radiation_overlap_counts[group_name]),
                "radiation_overlap_fraction": format_fraction(
                    metrics.radiation_overlap_counts[group_name],
                    patient_count,
                ),
                "regimen_context_count": str(metrics.regimen_context_counts[group_name]),
                "regimen_context_fraction": format_fraction(
                    metrics.regimen_context_counts[group_name],
                    patient_count,
                ),
                "timing_coverage_count": str(metrics.timing_coverage_counts[group_name]),
                "timing_coverage_fraction": format_fraction(
                    metrics.timing_coverage_counts[group_name],
                    patient_count,
                ),
                "notes": notes,
            }
        )
    return rows


def build_confounder_coverage_rows(
    run_id: str,
    metrics: GroupMetrics,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        row_count = metrics.group_counts[group_name]
        for field_name in REQUIRED_CONFOUNDER_FIELDS:
            notes = ""
            if row_count == 0:
                notes = "No patients assigned to this provisional group in the current run."
            elif (
                field_name in COVERAGE_READY_FIELDS
                and metrics.coverage_non_missing_fractions[(group_name, field_name)]
                < CONFOUNDER_NON_MISSING_TOLERANCE
            ):
                notes = (
                    "Non-missing coverage falls below the descriptive threshold "
                    f"({CONFOUNDER_NON_MISSING_TOLERANCE:.2f})."
                )
            rows.append(
                {
                    "treatment_groups_descriptive_v1_run_id": run_id,
                    "treatment_group_v1": group_name,
                    "field_name": field_name,
                    "row_count": str(row_count),
                    "non_missing_count": str(metrics.coverage_non_missing_counts[(group_name, field_name)]),
                    "missing_like_count": str(metrics.coverage_missing_like_counts[(group_name, field_name)]),
                    "missing_like_fraction": format_fraction(
                        metrics.coverage_missing_like_counts[(group_name, field_name)],
                        row_count,
                    ),
                    "distinct_non_missing_count": str(
                        metrics.coverage_distinct_non_missing_counts[(group_name, field_name)]
                    ),
                    "notes": notes,
                }
            )
    return rows


def build_value_composition_rows(
    run_id: str,
    joined_rows: list[dict[str, str]],
    metrics: GroupMetrics,
) -> list[dict[str, str]]:
    rows_by_group = {group_name: [] for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    for row in joined_rows:
        rows_by_group[row["treatment_group_v1"]].append(row)

    value_rows: list[dict[str, str]] = []
    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        group_rows = rows_by_group[group_name]
        patient_count = len(group_rows)
        if patient_count == 0:
            continue

        age_values = [
            age_value
            for age_value in (parse_float_or_none(row.get("age_at_diagnosis", "")) for row in group_rows)
            if age_value is not None
        ]
        value_rows.append(
            {
                "treatment_groups_descriptive_v1_run_id": run_id,
                "treatment_group_v1": group_name,
                "field_name": "age_at_diagnosis",
                "value_summary_type": "numeric_summary",
                "group_patient_count": str(patient_count),
                "non_missing_count": str(len(age_values)),
                "missing_like_count": str(patient_count - len(age_values)),
                "distinct_non_missing_count": str(len({value for value in age_values})),
                "value_rank": "",
                "value_label": "",
                "value_count": "",
                "value_fraction": "",
                "numeric_min": format_numeric(min(age_values) if age_values else None),
                "numeric_median": format_numeric(statistics.median(age_values) if age_values else None),
                "numeric_max": format_numeric(max(age_values) if age_values else None),
                "notes": "Age summary uses non-missing values only.",
            }
        )

        for field_name in CATEGORICAL_COMPOSITION_FIELDS:
            non_missing_values = [
                str(row.get(field_name, ""))
                for row in group_rows
                if not is_missing_like(row.get(field_name, ""))
            ]
            value_counts = Counter(non_missing_values)
            sorted_values = sorted(value_counts.items(), key=lambda item: (-item[1], item[0]))
            if not sorted_values:
                value_rows.append(
                    {
                        "treatment_groups_descriptive_v1_run_id": run_id,
                        "treatment_group_v1": group_name,
                        "field_name": field_name,
                        "value_summary_type": "categorical_no_non_missing_values",
                        "group_patient_count": str(patient_count),
                        "non_missing_count": "0",
                        "missing_like_count": str(patient_count),
                        "distinct_non_missing_count": "0",
                        "value_rank": "",
                        "value_label": "",
                        "value_count": "",
                        "value_fraction": "",
                        "numeric_min": "",
                        "numeric_median": "",
                        "numeric_max": "",
                        "notes": "All values are missing-like within this provisional group.",
                    }
                )
                continue

            distinct_count = len(sorted_values)
            truncated = distinct_count > MAX_TOP_VALUES
            for rank, (value_label, value_count) in enumerate(sorted_values[:MAX_TOP_VALUES], start=1):
                notes = (
                    f"Showing the top {MAX_TOP_VALUES} values out of {distinct_count} distinct non-missing values."
                    if truncated
                    else "Showing all distinct non-missing values for this field within the provisional group."
                )
                value_rows.append(
                    {
                        "treatment_groups_descriptive_v1_run_id": run_id,
                        "treatment_group_v1": group_name,
                        "field_name": field_name,
                        "value_summary_type": "categorical_top_value",
                        "group_patient_count": str(patient_count),
                        "non_missing_count": str(len(non_missing_values)),
                        "missing_like_count": str(patient_count - len(non_missing_values)),
                        "distinct_non_missing_count": str(distinct_count),
                        "value_rank": str(rank),
                        "value_label": value_label,
                        "value_count": str(value_count),
                        "value_fraction": format_fraction(value_count, patient_count),
                        "numeric_min": "",
                        "numeric_median": "",
                        "numeric_max": "",
                        "notes": notes,
                    }
                )
    return value_rows


def build_manual_review_summary_rows(
    run_id: str,
    metrics: GroupMetrics,
    grouping_conflict_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    conflict_counter: Counter[tuple[str, str]] = Counter()
    for conflict_row in grouping_conflict_rows:
        group_name = str(conflict_row["treatment_group_v1"]).strip()
        conflict_type = str(conflict_row["conflict_type"]).strip()
        conflict_counter[(group_name, conflict_type)] += 1

    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        manual_review_count = metrics.manual_review_counts[group_name]
        group_patient_count = metrics.group_counts[group_name]
        group_notes = ""
        if manual_review_count == 0:
            group_notes = "No manual-review-required patients remain in this provisional group."
        elif manual_review_count == group_patient_count:
            group_notes = "Manual-review burden spans the full provisional group."
        rows.append(
            {
                "treatment_groups_descriptive_v1_run_id": run_id,
                "treatment_group_v1": group_name,
                "summary_level": "group_total",
                "group_patient_count": str(group_patient_count),
                "manual_review_count": str(manual_review_count),
                "manual_review_fraction": format_fraction(manual_review_count, group_patient_count),
                "conflict_type": "",
                "conflict_type_count": "",
                "conflict_type_fraction_within_manual_review": "",
                "notes": group_notes,
            }
        )

        conflict_items = [
            (conflict_type, count)
            for (counter_group_name, conflict_type), count in conflict_counter.items()
            if counter_group_name == group_name
        ]
        for conflict_type, count in sorted(conflict_items, key=lambda item: (-item[1], item[0])):
            rows.append(
                {
                    "treatment_groups_descriptive_v1_run_id": run_id,
                    "treatment_group_v1": group_name,
                    "summary_level": "conflict_type",
                    "group_patient_count": str(group_patient_count),
                    "manual_review_count": str(manual_review_count),
                    "manual_review_fraction": format_fraction(manual_review_count, group_patient_count),
                    "conflict_type": conflict_type,
                    "conflict_type_count": str(count),
                    "conflict_type_fraction_within_manual_review": format_fraction(count, manual_review_count),
                    "notes": "Conflict subtype count carried forward from patient_treatment_grouping_v1_conflict_audit.tsv.",
                }
            )
    return rows


def determine_readiness_interpretation(
    workflow_inputs: WorkflowInputs,
    metrics: GroupMetrics,
) -> tuple[str, str]:
    if metrics.total_cohort_size <= 0 or not metrics.present_groups:
        return (
            READINESS_NOT_READY,
            "The grouped descriptive review could not establish a non-empty provisional treatment grouping cohort.",
        )

    all_large_groups_within_manual_review_tolerance = all(
        group_name in metrics.groups_manual_review_within_tolerance for group_name in metrics.large_groups
    )
    all_large_groups_within_coverage_tolerance = all(
        group_name in metrics.groups_all_ready_fields_within_tolerance for group_name in metrics.large_groups
    )
    no_manual_review_remaining = sum(metrics.manual_review_counts.values()) == 0
    no_drug_name_normalization_pending = not bool(
        workflow_inputs.grouping_run_log.get("rules", {}).get("no_drug_name_normalization", False)
    )

    if (
        metrics.large_groups
        and all_large_groups_within_manual_review_tolerance
        and all_large_groups_within_coverage_tolerance
        and no_manual_review_remaining
        and no_drug_name_normalization_pending
    ):
        return (
            READINESS_TREATMENT_ARM_FREEZE_REVIEW_NEXT,
            "All large provisional groups satisfy the descriptive manual-review and confounder-coverage checks, "
            "and the saved upstream blockers are no longer present.",
        )

    return (
        READINESS_DESCRIPTIVE_COMPLETED,
        "Grouped descriptive review is complete, but treatment-arm freeze review remains blocked by the current "
        "upstream state: drug names are still unnormalized and manual review remains concentrated in unresolved "
        "provisional groups.",
    )


def build_summary_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
    metrics: GroupMetrics,
) -> list[dict[str, str]]:
    summary_rows: list[dict[str, str]] = []

    def add(section: str, metric: str, value: str, notes: str) -> None:
        summary_rows.append(
            {
                "treatment_groups_descriptive_v1_run_id": run_id,
                "summary_section": section,
                "summary_metric": metric,
                "summary_value": value,
                "notes": notes,
            }
        )

    readiness_interpretation, readiness_notes = determine_readiness_interpretation(workflow_inputs, metrics)
    total_manual_review = sum(metrics.manual_review_counts.values())
    large_groups_meeting_both = [
        group_name
        for group_name in metrics.large_groups
        if group_name in metrics.groups_manual_review_within_tolerance
        and group_name in metrics.groups_all_ready_fields_within_tolerance
    ]

    add(
        "design",
        "review_scope",
        "grouped_descriptive_review_only",
        "This workflow remains descriptive only and does not freeze treatment arms or estimate effects.",
    )
    add(
        "design",
        "input_layer",
        "patient_treatment_grouping_v1_plus_os_endpoint_v1_plus_baseline_analysis_v1_validated_against_feature_set_lineage",
        "The review uses only saved downstream tables on disk and validates the feature-set/model-input lineage.",
    )
    add(
        "design",
        "no_drug_name_normalization_yet",
        "true"
        if workflow_inputs.grouping_run_log.get("rules", {}).get("no_drug_name_normalization", False)
        else "false",
        "Drug-name normalization remains out of scope for this descriptive review layer.",
    )
    add(
        "thresholds",
        "large_group_patient_count_threshold",
        str(MIN_GROUP_SIZE_FOR_FREEZE_REVIEW),
        "Large-group threshold used for descriptive readiness triage only.",
    )
    add(
        "thresholds",
        "manual_review_fraction_tolerance",
        format_float(MANUAL_REVIEW_FRACTION_TOLERANCE),
        "A provisional group is within tolerance when its manual-review fraction is less than or equal to this value.",
    )
    add(
        "thresholds",
        "adequate_non_missing_coverage_threshold",
        format_float(CONFOUNDER_NON_MISSING_TOLERANCE),
        "Adequate receptor/stage coverage uses non-missing fraction greater than or equal to this value.",
    )
    add(
        "counts",
        "total_cohort_size",
        str(metrics.total_cohort_size),
        "Total OS-linked patient rows reviewed in treatment_groups_descriptive_v1.",
    )
    add(
        "counts",
        "present_treatment_group_count",
        str(len(metrics.present_groups)),
        "Number of provisional treatment groups with at least one patient in the current grouping run.",
    )
    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        add(
            "counts",
            f"{group_name}_count",
            str(metrics.group_counts[group_name]),
            "Patient count carried forward from patient_treatment_grouping_v1.tsv.",
        )
    add(
        "counts",
        "patients_requiring_manual_review",
        str(total_manual_review),
        "Manual-review count carried forward from treatment_group_v1_requires_manual_review.",
    )
    add(
        "counts",
        "patients_requiring_manual_review_fraction",
        format_fraction(total_manual_review, metrics.total_cohort_size),
        "Manual-review fraction across the full OS-linked grouping cohort.",
    )
    add(
        "coverage",
        "groups_with_patient_count_ge_50",
        str(len(metrics.large_groups)),
        "Count of provisional groups with at least 50 patients.",
    )
    add(
        "coverage",
        "groups_with_manual_review_fraction_le_tolerance",
        str(len(metrics.groups_manual_review_within_tolerance)),
        "Count of non-empty provisional groups whose manual-review fraction is within the descriptive tolerance.",
    )
    for field_name in COVERAGE_READY_FIELDS:
        adequate_groups = [
            group_name
            for group_name in metrics.present_groups
            if metrics.coverage_non_missing_fractions[(group_name, field_name)]
            >= CONFOUNDER_NON_MISSING_TOLERANCE
        ]
        add(
            "coverage",
            f"groups_with_adequate_{field_name}_coverage",
            str(len(adequate_groups)),
            "Count of non-empty provisional groups meeting the descriptive non-missing coverage threshold.",
        )
    add(
        "coverage",
        "groups_with_adequate_all_required_coverage",
        str(len(metrics.groups_all_ready_fields_within_tolerance)),
        "Count of non-empty provisional groups meeting the descriptive threshold for ER, PR, HER2, and pathologic stage.",
    )
    add(
        "coverage",
        "large_groups_meeting_both_manual_review_and_coverage_tolerances",
        str(len(large_groups_meeting_both)),
        "Large provisional groups that meet both descriptive tolerances at the same time.",
    )
    add(
        "readiness",
        "descriptive_readiness_interpretation",
        readiness_interpretation,
        readiness_notes,
    )
    add(
        "readiness",
        "treatment_arm_freeze_review_blocked",
        "true" if readiness_interpretation != READINESS_TREATMENT_ARM_FREEZE_REVIEW_NEXT else "false",
        "Treatment-arm freeze remains blocked unless the descriptive evidence and upstream blockers are resolved.",
    )
    add(
        "readiness",
        "treatment_recommendation_modeling_status",
        TREATMENT_RECOMMENDATION_MODELING_STATUS,
        "Treatment recommendation modeling remains out of scope.",
    )
    add(
        "readiness",
        "largest_manual_review_group",
        max(metrics.group_counts, key=lambda name: metrics.manual_review_counts[name]),
        "Provisional group with the highest absolute manual-review count.",
    )
    return summary_rows


def validate_outputs(
    workflow_inputs: WorkflowInputs,
    joined_rows: list[dict[str, str]],
    group_summary_rows: list[dict[str, str]],
    confounder_coverage_rows: list[dict[str, str]],
    value_composition_rows: list[dict[str, str]],
    manual_review_summary_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    metrics: GroupMetrics,
    row_order_preserved: bool,
    all_os_patients_accounted_for: bool,
    no_duplicate_grouping_barcodes: bool,
    grouping_lineage_matches_inputs: bool,
) -> dict[str, Any]:
    group_summary_lookup = {row["treatment_group_v1"]: row for row in group_summary_rows}
    summary_metric_lookup = {row["summary_metric"]: row["summary_value"] for row in summary_rows}
    manual_review_group_rows = [
        row for row in manual_review_summary_rows if row["summary_level"] == "group_total"
    ]
    manual_review_group_lookup = {row["treatment_group_v1"]: row for row in manual_review_group_rows}

    group_counts_reconcile = all(
        group_name in group_summary_lookup
        and parse_int(group_summary_lookup[group_name]["patient_count"], "patient_count")
        == metrics.group_counts[group_name]
        and parse_int(group_summary_lookup[group_name]["os_event_count"], "os_event_count")
        == metrics.os_event_counts[group_name]
        and parse_int(group_summary_lookup[group_name]["manual_review_count"], "manual_review_count")
        == metrics.manual_review_counts[group_name]
        and parse_int(group_summary_lookup[group_name]["radiation_overlap_count"], "radiation_overlap_count")
        == metrics.radiation_overlap_counts[group_name]
        and parse_int(group_summary_lookup[group_name]["regimen_context_count"], "regimen_context_count")
        == metrics.regimen_context_counts[group_name]
        and parse_int(group_summary_lookup[group_name]["timing_coverage_count"], "timing_coverage_count")
        == metrics.timing_coverage_counts[group_name]
        for group_name in APPROVED_TREATMENT_GROUP_ORDER
    )

    confounder_coverage_row_count_expected = len(APPROVED_TREATMENT_GROUP_ORDER) * len(REQUIRED_CONFOUNDER_FIELDS)
    confounder_coverage_reconciles = (
        len(confounder_coverage_rows) == confounder_coverage_row_count_expected
        and all(
            parse_int(row["row_count"], "row_count") == metrics.group_counts[row["treatment_group_v1"]]
            and parse_int(row["non_missing_count"], "non_missing_count")
            == metrics.coverage_non_missing_counts[(row["treatment_group_v1"], row["field_name"])]
            and parse_int(row["missing_like_count"], "missing_like_count")
            == metrics.coverage_missing_like_counts[(row["treatment_group_v1"], row["field_name"])]
            and parse_int(row["distinct_non_missing_count"], "distinct_non_missing_count")
            == metrics.coverage_distinct_non_missing_counts[(row["treatment_group_v1"], row["field_name"])]
            for row in confounder_coverage_rows
        )
    )

    manual_review_group_totals_reconcile = len(manual_review_group_rows) == len(APPROVED_TREATMENT_GROUP_ORDER) and all(
        group_name in manual_review_group_lookup
        and parse_int(manual_review_group_lookup[group_name]["group_patient_count"], "group_patient_count")
        == metrics.group_counts[group_name]
        and parse_int(manual_review_group_lookup[group_name]["manual_review_count"], "manual_review_count")
        == metrics.manual_review_counts[group_name]
        for group_name in APPROVED_TREATMENT_GROUP_ORDER
    )

    conflict_subtype_rows = [
        row for row in manual_review_summary_rows if row["summary_level"] == "conflict_type"
    ]
    conflict_subtype_counts_sum_to_manual_review = True
    subtype_counter: Counter[str] = Counter()
    for row in conflict_subtype_rows:
        subtype_counter[row["treatment_group_v1"]] += parse_int(row["conflict_type_count"], "conflict_type_count")
    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        if subtype_counter.get(group_name, 0) != metrics.manual_review_counts[group_name]:
            if metrics.manual_review_counts[group_name] != 0:
                conflict_subtype_counts_sum_to_manual_review = False

    summary_counts_reconcile = (
        summary_metric_lookup.get("total_cohort_size") == str(metrics.total_cohort_size)
        and summary_metric_lookup.get("patients_requiring_manual_review")
        == str(sum(metrics.manual_review_counts.values()))
        and summary_metric_lookup.get("groups_with_patient_count_ge_50") == str(len(metrics.large_groups))
        and all(
            summary_metric_lookup.get(f"{group_name}_count") == str(metrics.group_counts[group_name])
            for group_name in APPROVED_TREATMENT_GROUP_ORDER
        )
    )

    present_group_set = set(metrics.present_groups)
    value_composition_nonzero_age_rows = {
        row["treatment_group_v1"]
        for row in value_composition_rows
        if row["field_name"] == "age_at_diagnosis" and row["value_summary_type"] == "numeric_summary"
    }
    value_composition_has_age_row_for_present_groups = value_composition_nonzero_age_rows == present_group_set

    readiness_allowed = summary_metric_lookup.get("descriptive_readiness_interpretation", "") in {
        READINESS_NOT_READY,
        READINESS_DESCRIPTIVE_COMPLETED,
        READINESS_TREATMENT_ARM_FREEZE_REVIEW_NEXT,
    }
    modeling_scope_stays_out = (
        summary_metric_lookup.get("treatment_recommendation_modeling_status", "")
        == TREATMENT_RECOMMENDATION_MODELING_STATUS
    )

    passed = all(
        [
            len(workflow_inputs.grouping_rows) > 0,
            len(workflow_inputs.os_endpoint_rows) > 0,
            len(workflow_inputs.baseline_feature_set_rows) > 0,
            len(workflow_inputs.baseline_feature_set_audit_map_rows) > 0,
            len(workflow_inputs.baseline_analysis_rows) > 0,
            len(joined_rows) == len(workflow_inputs.grouping_rows) == len(workflow_inputs.os_endpoint_rows),
            row_order_preserved,
            all_os_patients_accounted_for,
            no_duplicate_grouping_barcodes,
            grouping_lineage_matches_inputs,
            group_counts_reconcile,
            confounder_coverage_reconciles,
            manual_review_group_totals_reconcile,
            conflict_subtype_counts_sum_to_manual_review,
            summary_counts_reconcile,
            value_composition_has_age_row_for_present_groups,
            readiness_allowed,
            modeling_scope_stays_out,
        ]
    )

    return {
        "passed": passed,
        "required_upstream_pointers_found": True,
        "required_source_tables_found": True,
        "grouping_run_log_completed": True,
        "os_endpoint_run_log_completed": True,
        "baseline_model_input_run_log_completed": True,
        "baseline_feature_set_run_log_completed": True,
        "baseline_analysis_run_log_completed": True,
        "grouping_row_count_positive": len(workflow_inputs.grouping_rows) > 0,
        "os_endpoint_row_count_positive": len(workflow_inputs.os_endpoint_rows) > 0,
        "baseline_feature_set_row_count_positive": len(workflow_inputs.baseline_feature_set_rows) > 0,
        "baseline_analysis_row_count_positive": len(workflow_inputs.baseline_analysis_rows) > 0,
        "grouping_row_count_matches_os_cohort": len(joined_rows)
        == len(workflow_inputs.grouping_rows)
        == len(workflow_inputs.os_endpoint_rows),
        "row_order_preserved_from_os_endpoint": row_order_preserved,
        "all_os_cohort_patients_accounted_once_via_grouping_linkage": all_os_patients_accounted_for,
        "no_grouping_duplicate_barcodes": no_duplicate_grouping_barcodes,
        "grouping_lineage_matches_os_and_baseline_inputs": grouping_lineage_matches_inputs,
        "group_counts_reconcile_to_patient_treatment_grouping_table": group_counts_reconcile,
        "confounder_coverage_row_counts_reconcile_to_group_counts": confounder_coverage_reconciles,
        "manual_review_group_totals_reconcile_to_group_summary": manual_review_group_totals_reconcile,
        "manual_review_conflict_subtype_counts_sum_to_manual_review": conflict_subtype_counts_sum_to_manual_review,
        "summary_counts_reconcile_to_group_summary": summary_counts_reconcile,
        "value_composition_has_age_row_for_each_present_group": value_composition_has_age_row_for_present_groups,
        "summary_readiness_interpretation_allowed": readiness_allowed,
        "treatment_recommendation_modeling_remains_out_of_scope": modeling_scope_stays_out,
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
        "grouping_row_count": len(workflow_inputs.grouping_rows),
        "os_endpoint_row_count": len(workflow_inputs.os_endpoint_rows),
        "baseline_feature_set_row_count": len(workflow_inputs.baseline_feature_set_rows),
        "baseline_feature_set_audit_map_row_count": len(workflow_inputs.baseline_feature_set_audit_map_rows),
        "baseline_analysis_row_count": len(workflow_inputs.baseline_analysis_rows),
        "group_summary_row_count": len(group_summary_rows),
        "confounder_coverage_row_count": len(confounder_coverage_rows),
        "value_composition_row_count": len(value_composition_rows),
        "manual_review_summary_row_count": len(manual_review_summary_rows),
        "summary_row_count": len(summary_rows),
    }


def build_latest_pointer_payload(
    *,
    run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, str]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "treatment_groups_descriptive_v1_run_id": run_id,
        "patient_treatment_grouping_v1_run_id": str(
            workflow_inputs.grouping_latest_pointer["patient_treatment_grouping_v1_run_id"]
        ),
        "patient_treatment_profile_v1_run_id": str(
            workflow_inputs.grouping_latest_pointer["patient_treatment_profile_v1_run_id"]
        ),
        "treatment_os_overlap_v1_run_id": str(
            workflow_inputs.grouping_latest_pointer["treatment_os_overlap_v1_run_id"]
        ),
        "os_endpoint_v1_run_id": str(workflow_inputs.os_endpoint_latest_pointer["os_endpoint_v1_run_id"]),
        "baseline_model_input_v1_run_id": str(
            workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
        ),
        "baseline_feature_set_v1_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["baseline_feature_set_v1_run_id"]
        ),
        "baseline_analysis_v1_run_id": str(
            workflow_inputs.baseline_analysis_latest_pointer["baseline_analysis_v1_run_id"]
        ),
        "cohort_v1_build_id": str(workflow_inputs.grouping_latest_pointer["cohort_v1_build_id"]),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "treatment_groups_descriptive_v1_group_summary_tsv": repo_relative(
            output_paths["group_summary_tsv"],
            paths.repo_root,
        ),
        "treatment_groups_descriptive_v1_confounder_coverage_tsv": repo_relative(
            output_paths["confounder_coverage_tsv"],
            paths.repo_root,
        ),
        "treatment_groups_descriptive_v1_value_composition_tsv": repo_relative(
            output_paths["value_composition_tsv"],
            paths.repo_root,
        ),
        "treatment_groups_descriptive_v1_manual_review_summary_tsv": repo_relative(
            output_paths["manual_review_summary_tsv"],
            paths.repo_root,
        ),
        "treatment_groups_descriptive_v1_summary_tsv": repo_relative(
            output_paths["summary_tsv"],
            paths.repo_root,
        ),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "patient_treatment_grouping_v1_latest_json": repo_relative(paths.grouping_latest_pointer, paths.repo_root),
        "os_endpoint_v1_latest_json": repo_relative(paths.os_endpoint_latest_pointer, paths.repo_root),
        "baseline_model_input_v1_latest_json": repo_relative(
            paths.baseline_model_input_latest_pointer,
            paths.repo_root,
        ),
        "baseline_feature_set_v1_latest_json": repo_relative(
            paths.baseline_feature_set_latest_pointer,
            paths.repo_root,
        ),
        "baseline_analysis_v1_latest_json": repo_relative(
            paths.baseline_analysis_latest_pointer,
            paths.repo_root,
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    audit_run_dir = paths.audit_runs_root / run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)
        helper_module.create_run_directory(audit_run_dir)

        output_paths = {
            "audit_run_directory": audit_run_dir,
            "group_summary_tsv": audit_run_dir / "treatment_groups_descriptive_v1_group_summary.tsv",
            "confounder_coverage_tsv": audit_run_dir / "treatment_groups_descriptive_v1_confounder_coverage.tsv",
            "value_composition_tsv": audit_run_dir / "treatment_groups_descriptive_v1_value_composition.tsv",
            "manual_review_summary_tsv": audit_run_dir / "treatment_groups_descriptive_v1_manual_review_summary.tsv",
            "summary_tsv": audit_run_dir / "treatment_groups_descriptive_v1_summary.tsv",
            "run_log_json": run_log_path,
        }

        (
            joined_rows,
            row_order_preserved,
            all_os_patients_accounted_for,
            no_duplicate_grouping_barcodes,
            grouping_lineage_matches_inputs,
        ) = build_joined_rows(workflow_inputs)
        metrics = build_group_metrics(joined_rows)
        group_summary_rows = build_group_summary_rows(run_id, metrics)
        confounder_coverage_rows = build_confounder_coverage_rows(run_id, metrics)
        value_composition_rows = build_value_composition_rows(run_id, joined_rows, metrics)
        manual_review_summary_rows = build_manual_review_summary_rows(
            run_id,
            metrics,
            workflow_inputs.grouping_conflict_rows,
        )
        summary_rows = build_summary_rows(run_id, workflow_inputs, metrics)

        validation = validate_outputs(
            workflow_inputs=workflow_inputs,
            joined_rows=joined_rows,
            group_summary_rows=group_summary_rows,
            confounder_coverage_rows=confounder_coverage_rows,
            value_composition_rows=value_composition_rows,
            manual_review_summary_rows=manual_review_summary_rows,
            summary_rows=summary_rows,
            metrics=metrics,
            row_order_preserved=row_order_preserved,
            all_os_patients_accounted_for=all_os_patients_accounted_for,
            no_duplicate_grouping_barcodes=no_duplicate_grouping_barcodes,
            grouping_lineage_matches_inputs=grouping_lineage_matches_inputs,
        )
        if not validation["passed"]:
            raise TreatmentGroupsDescriptiveV1Error(
                "Output validation did not pass. Failed checks: "
                + str({key: value for key, value in validation.items() if value is False})
            )

        helper_module.write_dict_rows_tsv(
            output_paths["group_summary_tsv"],
            GROUP_SUMMARY_FIELDNAMES,
            group_summary_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_paths["confounder_coverage_tsv"],
            CONFOUNDER_COVERAGE_FIELDNAMES,
            confounder_coverage_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_paths["value_composition_tsv"],
            VALUE_COMPOSITION_FIELDNAMES,
            value_composition_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_paths["manual_review_summary_tsv"],
            MANUAL_REVIEW_SUMMARY_FIELDNAMES,
            manual_review_summary_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_paths["summary_tsv"],
            SUMMARY_FIELDNAMES,
            summary_rows,
        )

        latest_pointer_payload = build_latest_pointer_payload(
            run_id=run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )
        completed_at = utc_now()
        run_log_payload: dict[str, Any] = {
            "status": "completed",
            "treatment_groups_descriptive_v1_run_id": run_id,
            "patient_treatment_grouping_v1_run_id": str(
                workflow_inputs.grouping_latest_pointer["patient_treatment_grouping_v1_run_id"]
            ),
            "patient_treatment_profile_v1_run_id": str(
                workflow_inputs.grouping_latest_pointer["patient_treatment_profile_v1_run_id"]
            ),
            "treatment_os_overlap_v1_run_id": str(
                workflow_inputs.grouping_latest_pointer["treatment_os_overlap_v1_run_id"]
            ),
            "os_endpoint_v1_run_id": str(workflow_inputs.os_endpoint_latest_pointer["os_endpoint_v1_run_id"]),
            "baseline_model_input_v1_run_id": str(
                workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
            ),
            "baseline_feature_set_v1_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["baseline_feature_set_v1_run_id"]
            ),
            "baseline_analysis_v1_run_id": str(
                workflow_inputs.baseline_analysis_latest_pointer["baseline_analysis_v1_run_id"]
            ),
            "cohort_v1_build_id": str(workflow_inputs.grouping_latest_pointer["cohort_v1_build_id"]),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "patient_treatment_grouping_v1_latest_json": repo_relative(
                    paths.grouping_latest_pointer,
                    paths.repo_root,
                ),
                "os_endpoint_v1_latest_json": repo_relative(paths.os_endpoint_latest_pointer, paths.repo_root),
                "baseline_model_input_v1_latest_json": repo_relative(
                    paths.baseline_model_input_latest_pointer,
                    paths.repo_root,
                ),
                "baseline_feature_set_v1_latest_json": repo_relative(
                    paths.baseline_feature_set_latest_pointer,
                    paths.repo_root,
                ),
                "baseline_analysis_v1_latest_json": repo_relative(
                    paths.baseline_analysis_latest_pointer,
                    paths.repo_root,
                ),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "treatment_groups_descriptive_v1_group_summary_tsv": repo_relative(
                    output_paths["group_summary_tsv"],
                    paths.repo_root,
                ),
                "treatment_groups_descriptive_v1_confounder_coverage_tsv": repo_relative(
                    output_paths["confounder_coverage_tsv"],
                    paths.repo_root,
                ),
                "treatment_groups_descriptive_v1_value_composition_tsv": repo_relative(
                    output_paths["value_composition_tsv"],
                    paths.repo_root,
                ),
                "treatment_groups_descriptive_v1_manual_review_summary_tsv": repo_relative(
                    output_paths["manual_review_summary_tsv"],
                    paths.repo_root,
                ),
                "treatment_groups_descriptive_v1_summary_tsv": repo_relative(
                    output_paths["summary_tsv"],
                    paths.repo_root,
                ),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation,
            "rules": {
                "unit_of_analysis": "patient/case",
                "output_layer": "treatment_groups_descriptive_v1",
                "review_scope": "grouped_descriptive_review_only",
                "no_drug_name_normalization": True,
                "no_treatment_group_redefinition": True,
                "no_silent_exclusion_of_mixed_patients": True,
                "no_modeling": True,
                "no_treatment_recommendation_outputs": True,
                "no_statistical_testing_beyond_descriptive_counts_and_fractions": True,
                "manual_review_fraction_tolerance": format_float(MANUAL_REVIEW_FRACTION_TOLERANCE),
                "adequate_non_missing_coverage_threshold": format_float(CONFOUNDER_NON_MISSING_TOLERANCE),
                "large_group_patient_count_threshold": str(MIN_GROUP_SIZE_FOR_FREEZE_REVIEW),
            },
            "counts": {
                "total_cohort_size": metrics.total_cohort_size,
                "present_treatment_group_count": len(metrics.present_groups),
                "groups_with_patient_count_ge_50": len(metrics.large_groups),
                "groups_with_manual_review_fraction_le_tolerance": len(
                    metrics.groups_manual_review_within_tolerance
                ),
                "groups_with_adequate_all_required_coverage": len(
                    metrics.groups_all_ready_fields_within_tolerance
                ),
                "patients_requiring_manual_review": sum(metrics.manual_review_counts.values()),
                **{
                    f"group_{group_name}_count": metrics.group_counts[group_name]
                    for group_name in APPROVED_TREATMENT_GROUP_ORDER
                },
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "patient_treatment_grouping_v1_latest_pointer": workflow_inputs.grouping_latest_pointer,
                "os_endpoint_v1_latest_pointer": workflow_inputs.os_endpoint_latest_pointer,
                "baseline_model_input_v1_latest_pointer": workflow_inputs.baseline_model_input_latest_pointer,
                "baseline_feature_set_v1_latest_pointer": workflow_inputs.baseline_feature_set_latest_pointer,
                "baseline_analysis_v1_latest_pointer": workflow_inputs.baseline_analysis_latest_pointer,
            },
        }
        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        completed_at = utc_now()
        failure_payload = {
            "status": "failed",
            "treatment_groups_descriptive_v1_run_id": run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(completed_at),
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "outputs": {
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
            },
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    counts = run_log.get("counts", {})
    outputs = run_log.get("outputs", {})
    summary_lookup = {}
    summary_path = Path(run_log["repo_root"]) / outputs["treatment_groups_descriptive_v1_summary_tsv"]
    if summary_path.exists():
        with summary_path.open("r", encoding="utf-8") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            for raw_line in handle:
                values = raw_line.rstrip("\n").split("\t")
                row = dict(zip(header, values))
                summary_lookup[row["summary_metric"]] = row["summary_value"]
    print("TCGA-BRCA treatment groups descriptive review v1 completed")
    print(f"  Run ID                  : {run_log['treatment_groups_descriptive_v1_run_id']}")
    print(f"  Grouping run ID         : {run_log['patient_treatment_grouping_v1_run_id']}")
    print(f"  OS endpoint run ID      : {run_log['os_endpoint_v1_run_id']}")
    print(f"  Cohort size             : {counts.get('total_cohort_size', 0)}")
    print(f"  Manual review patients  : {counts.get('patients_requiring_manual_review', 0)}")
    print(f"  Large groups >= 50      : {counts.get('groups_with_patient_count_ge_50', 0)}")
    print(
        "  Readiness interpretation: "
        f"{summary_lookup.get('descriptive_readiness_interpretation', '[not available]')}"
    )
    print(f"  Audit run dir           : {outputs.get('audit_run_directory', '')}")
    print(f"  Latest pointer          : {outputs.get('latest_pointer_json', '')}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
