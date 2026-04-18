#!/usr/bin/env python
"""Build an auditable TCGA-BRCA provisional patient-level treatment grouping v1."""

from __future__ import annotations

import importlib.util
import json
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
READINESS_GROUPED_DESCRIPTIVE_REVIEW = "ready_for_grouped_descriptive_review"
TREATMENT_ARM_FREEZE_REVIEW_STATUS = "blocked"
TREATMENT_RECOMMENDATION_MODELING_STATUS = "out_of_scope"
GROUPING_FIELDNAMES = [
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
    "treatment_group_v1_rule",
    "treatment_group_v1_requires_manual_review",
    "has_any_drug_row",
    "has_any_radiation_row",
    "drug_row_count",
    "radiation_row_count",
    "drug_therapy_type_values_json",
    "drug_therapy_type_single_or_mixed",
    "dominant_therapy_type_if_any",
    "dominant_therapy_type_rule",
    "regimen_context_values_json",
    "has_any_treatment_timing",
    "treatment_profile_status",
    "treatment_profile_requires_manual_review",
    "treatment_profile_flags_json",
    "grouping_flags_json",
]
SPEC_FIELDNAMES = [
    "patient_treatment_grouping_v1_run_id",
    "field_name",
    "source_origin",
    "field_category",
    "decision_rule",
    "notes",
]
ARM_SUMMARY_FIELDNAMES = [
    "patient_treatment_grouping_v1_run_id",
    "group_order",
    "treatment_group_v1",
    "patient_count",
    "patient_fraction_of_cohort",
    "os_event_count",
    "os_event_fraction_within_group",
    "manual_review_count",
    "manual_review_fraction_within_group",
    "radiation_overlap_count",
    "radiation_overlap_fraction_within_group",
    "regimen_context_count",
    "regimen_context_fraction_within_group",
    "treatment_timing_count",
    "treatment_timing_fraction_within_group",
    "notes",
]
CONFLICT_AUDIT_FIELDNAMES = [
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
    "treatment_group_v1",
    "treatment_group_v1_rule",
    "conflict_type",
    "source_evidence_summary",
    "recommended_next_handling",
    "review_priority",
    "drug_therapy_type_single_or_mixed",
    "dominant_therapy_type_if_any",
    "dominant_therapy_type_rule",
    "regimen_context_values_json",
    "has_any_treatment_timing",
    "treatment_profile_status",
    "treatment_profile_flags_json",
    "grouping_flags_json",
    "treatment_group_v1_requires_manual_review",
]
SUMMARY_FIELDNAMES = [
    "patient_treatment_grouping_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]


class PatientTreatmentGroupingV1Error(RuntimeError):
    """Raised when the patient treatment grouping v1 workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the patient treatment grouping v1 workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    patient_treatment_profile_latest_pointer: Path
    os_endpoint_v1_latest_pointer: Path
    baseline_model_input_v1_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    patient_treatment_profile_latest_pointer: dict[str, Any]
    patient_treatment_profile_run_log: dict[str, Any]
    os_endpoint_v1_latest_pointer: dict[str, Any]
    os_endpoint_v1_run_log: dict[str, Any]
    baseline_model_input_v1_latest_pointer: dict[str, Any]
    baseline_model_input_v1_run_log: dict[str, Any]
    patient_treatment_profile_rows: list[dict[str, str]]
    os_endpoint_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class GroupingMetrics:
    """Scalar and grouped metrics derived from the treatment grouping table."""

    total: int
    group_counts: dict[str, int]
    group_event_counts: dict[str, int]
    group_manual_review_counts: dict[str, int]
    group_radiation_counts: dict[str, int]
    group_regimen_context_counts: dict[str, int]
    group_timing_counts: dict[str, int]
    patients_requiring_manual_review: int
    patients_with_any_radiation: int
    patients_with_any_regimen_context: int
    patients_with_any_treatment_timing: int
    patients_with_compound_raw_therapy_type_present: int


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise PatientTreatmentGroupingV1Error(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise PatientTreatmentGroupingV1Error(f"Unable to create an import spec for: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)  # type: ignore[union-attr]
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


def parse_json_list(raw_value: str) -> list[str]:
    if not str(raw_value).strip():
        return []
    parsed = json.loads(raw_value)
    if not isinstance(parsed, list):
        raise PatientTreatmentGroupingV1Error(f"Expected a JSON list, got: {raw_value}")
    return [str(value) for value in parsed]


def json_list_contains(row: dict[str, str], field_name: str, expected_value: str) -> bool:
    raw_value = str(row.get(field_name, "")).strip()
    if not raw_value:
        return False
    return expected_value in parse_json_list(raw_value)


def format_fraction(numerator: int, denominator: int) -> str:
    if denominator <= 0:
        return "0.0000"
    return f"{numerator / denominator:.4f}"


def parse_int_or_none(value: str) -> int | None:
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return int(stripped)
    except ValueError:
        return None


def normalize_barcode(value: str) -> str:
    return value.strip().upper()


def normalize_uuid(value: str) -> str:
    return value.strip().upper()


def require_columns(rows: list[dict[str, str]], required_columns: set[str], label: str) -> None:
    if not rows:
        raise PatientTreatmentGroupingV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise PatientTreatmentGroupingV1Error(f"{label} is missing required columns: {sorted(missing)}")


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
        raise PatientTreatmentGroupingV1Error(f"{label} does not report status == completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise PatientTreatmentGroupingV1Error(f"{label} does not report validation.passed == true.")
    return run_log


def build_unique_lookup_by_barcode(
    rows: list[dict[str, str]],
    *,
    barcode_field: str,
    label: str,
) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        barcode = normalize_barcode(str(row.get(barcode_field, "")))
        if not barcode:
            raise PatientTreatmentGroupingV1Error(f"Blank {barcode_field} encountered in {label}.")
        if barcode in lookup:
            raise PatientTreatmentGroupingV1Error(
                f"Duplicate {barcode_field} '{barcode}' encountered in {label}."
            )
        lookup[barcode] = row
    return lookup


def has_any_regimen_context(row: dict[str, str]) -> bool:
    return str(row["regimen_context_values_json"]).strip() not in {"", "[]"}


def determine_treatment_group(profile_row: dict[str, str]) -> tuple[str, str]:
    has_any_drug_row = str(profile_row["has_any_drug_row"]).strip()
    single_or_mixed = str(profile_row["drug_therapy_type_single_or_mixed"]).strip()
    dominant = str(profile_row["dominant_therapy_type_if_any"]).strip()

    if has_any_drug_row == "no":
        return ("no_drug_record", "has_any_drug_row_no")
    if single_or_mixed == "missing_type_only":
        return ("missing_type_only", "drug_therapy_type_single_or_mixed_missing_type_only")
    if single_or_mixed == "mixed_type":
        return ("mixed_multi_type", "drug_therapy_type_single_or_mixed_mixed_type")
    if single_or_mixed == "single_type":
        if dominant == "Chemotherapy":
            return ("single_chemotherapy", "single_type_exact_raw_dominant_Chemotherapy")
        if dominant == "Hormone Therapy":
            return ("single_hormone_therapy", "single_type_exact_raw_dominant_Hormone_Therapy")
        if dominant == "Targeted Molecular therapy":
            return (
                "single_targeted_therapy",
                "single_type_exact_raw_dominant_Targeted_Molecular_therapy",
            )
        if dominant == "Immunotherapy":
            return ("single_immunotherapy", "single_type_exact_raw_dominant_Immunotherapy")
        return ("single_ancillary_or_other", "single_type_other_exact_raw_dominant_type")

    raise PatientTreatmentGroupingV1Error(
        "Unsupported patient treatment profile combination for grouping: "
        f"barcode={profile_row['bcr_patient_barcode']}, "
        f"has_any_drug_row={has_any_drug_row}, "
        f"drug_therapy_type_single_or_mixed={single_or_mixed}, "
        f"dominant_therapy_type_if_any={dominant or '[blank]'}"
    )


def build_grouping_flags(profile_row: dict[str, str], treatment_group: str) -> list[str]:
    flags: list[str] = []

    if treatment_group == "no_drug_record":
        flags.append("no_drug_record")
    if profile_row["treatment_profile_requires_manual_review"] == "yes":
        flags.append("upstream_manual_review_required")
    if treatment_group == "mixed_multi_type":
        flags.append("mixed_multi_type")
    if treatment_group == "missing_type_only":
        flags.append("missing_type_only")
    if treatment_group == "single_ancillary_or_other":
        flags.append("single_ancillary_or_other_fallback")
    if json_list_contains(profile_row, "treatment_profile_flags_json", "compound_raw_therapy_type_present"):
        flags.append("compound_raw_therapy_type_present")
    if profile_row["has_any_radiation_row"] == "yes":
        flags.append("has_any_radiation_overlap")
    if has_any_regimen_context(profile_row):
        flags.append("has_any_regimen_context")
    if profile_row["has_any_treatment_timing"] == "yes":
        flags.append("has_any_treatment_timing")
    if treatment_group == "no_drug_record" and profile_row["has_any_radiation_row"] == "yes":
        flags.append("no_drug_record_with_radiation_overlap")
    if json_list_contains(profile_row, "treatment_profile_flags_json", "drug_timing_window_aggregated_inverted"):
        flags.append("drug_timing_window_aggregated_inverted")
    if treatment_group == "mixed_multi_type" and profile_row["dominant_therapy_type_rule"] == "blank_due_to_tie":
        flags.append("mixed_multi_type_tied_dominant")
    if treatment_group == "mixed_multi_type" and str(profile_row["dominant_therapy_type_if_any"]).strip():
        flags.append("mixed_multi_type_with_provisional_dominant")

    return ordered_unique(flags)


def build_grouping_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
) -> tuple[list[dict[str, str]], GroupingMetrics]:
    profile_by_barcode = build_unique_lookup_by_barcode(
        workflow_inputs.patient_treatment_profile_rows,
        barcode_field="bcr_patient_barcode",
        label="patient_treatment_profile_v1.tsv",
    )

    grouping_rows: list[dict[str, str]] = []
    for os_row in workflow_inputs.os_endpoint_rows:
        barcode = normalize_barcode(str(os_row["bcr_patient_barcode"]))
        os_uuid = normalize_uuid(str(os_row["bcr_patient_uuid"]))
        profile_row = profile_by_barcode.get(barcode)

        if profile_row is None:
            raise PatientTreatmentGroupingV1Error(
                f"Missing patient_treatment_profile_v1.tsv row for barcode {barcode}."
            )
        if normalize_uuid(str(profile_row["bcr_patient_uuid"])) != os_uuid:
            raise PatientTreatmentGroupingV1Error(
                f"UUID mismatch between os_endpoint_v1.tsv and patient_treatment_profile_v1.tsv for barcode {barcode}."
            )
        for field_name in ["provisional_patient_row_id", "baseline_analysis_v1_row_id", "feature_set_v1_row_index"]:
            if str(profile_row[field_name]).strip() != str(os_row[field_name]).strip():
                raise PatientTreatmentGroupingV1Error(
                    f"{field_name} mismatch between os_endpoint_v1.tsv and patient_treatment_profile_v1.tsv "
                    f"for barcode {barcode}."
                )

        treatment_group, treatment_group_rule = determine_treatment_group(profile_row)
        grouping_flags = build_grouping_flags(profile_row, treatment_group)

        grouping_rows.append(
            {
                "patient_treatment_grouping_v1_run_id": run_id,
                "patient_treatment_profile_v1_run_id": str(
                    workflow_inputs.patient_treatment_profile_latest_pointer["patient_treatment_profile_v1_run_id"]
                ),
                "treatment_os_overlap_v1_run_id": str(profile_row["treatment_os_overlap_v1_run_id"]),
                "os_endpoint_v1_run_id": str(profile_row["os_endpoint_v1_run_id"]),
                "baseline_model_input_v1_run_id": str(profile_row["baseline_model_input_v1_run_id"]),
                "cohort_v1_build_id": str(profile_row["cohort_v1_build_id"]),
                "bcr_patient_barcode": str(os_row["bcr_patient_barcode"]),
                "bcr_patient_uuid": str(os_row["bcr_patient_uuid"]),
                "provisional_patient_row_id": str(os_row["provisional_patient_row_id"]),
                "baseline_analysis_v1_row_id": str(os_row["baseline_analysis_v1_row_id"]),
                "feature_set_v1_row_index": str(os_row["feature_set_v1_row_index"]),
                "os_event": str(os_row["os_event"]),
                "os_time_days": str(os_row["os_time_days"]),
                "treatment_group_v1": treatment_group,
                "treatment_group_v1_rule": treatment_group_rule,
                "treatment_group_v1_requires_manual_review": str(
                    profile_row["treatment_profile_requires_manual_review"]
                ),
                "has_any_drug_row": str(profile_row["has_any_drug_row"]),
                "has_any_radiation_row": str(profile_row["has_any_radiation_row"]),
                "drug_row_count": str(profile_row["drug_row_count"]),
                "radiation_row_count": str(profile_row["radiation_row_count"]),
                "drug_therapy_type_values_json": str(profile_row["drug_therapy_type_values_json"]),
                "drug_therapy_type_single_or_mixed": str(profile_row["drug_therapy_type_single_or_mixed"]),
                "dominant_therapy_type_if_any": str(profile_row["dominant_therapy_type_if_any"]),
                "dominant_therapy_type_rule": str(profile_row["dominant_therapy_type_rule"]),
                "regimen_context_values_json": str(profile_row["regimen_context_values_json"]),
                "has_any_treatment_timing": str(profile_row["has_any_treatment_timing"]),
                "treatment_profile_status": str(profile_row["treatment_profile_status"]),
                "treatment_profile_requires_manual_review": str(
                    profile_row["treatment_profile_requires_manual_review"]
                ),
                "treatment_profile_flags_json": str(profile_row["treatment_profile_flags_json"]),
                "grouping_flags_json": json_list(grouping_flags),
            }
        )

    metrics = build_grouping_metrics(grouping_rows)
    return grouping_rows, metrics


def build_grouping_metrics(grouping_rows: list[dict[str, str]]) -> GroupingMetrics:
    group_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    group_event_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    group_manual_review_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    group_radiation_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    group_regimen_context_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}
    group_timing_counts = {group_name: 0 for group_name in APPROVED_TREATMENT_GROUP_ORDER}

    for row in grouping_rows:
        group_name = row["treatment_group_v1"]
        group_counts[group_name] += 1
        group_event_counts[group_name] += int(row["os_event"])
        if row["treatment_group_v1_requires_manual_review"] == "yes":
            group_manual_review_counts[group_name] += 1
        if row["has_any_radiation_row"] == "yes":
            group_radiation_counts[group_name] += 1
        if has_any_regimen_context(row):
            group_regimen_context_counts[group_name] += 1
        if row["has_any_treatment_timing"] == "yes":
            group_timing_counts[group_name] += 1

    return GroupingMetrics(
        total=len(grouping_rows),
        group_counts=group_counts,
        group_event_counts=group_event_counts,
        group_manual_review_counts=group_manual_review_counts,
        group_radiation_counts=group_radiation_counts,
        group_regimen_context_counts=group_regimen_context_counts,
        group_timing_counts=group_timing_counts,
        patients_requiring_manual_review=sum(
            1 for row in grouping_rows if row["treatment_group_v1_requires_manual_review"] == "yes"
        ),
        patients_with_any_radiation=sum(1 for row in grouping_rows if row["has_any_radiation_row"] == "yes"),
        patients_with_any_regimen_context=sum(1 for row in grouping_rows if has_any_regimen_context(row)),
        patients_with_any_treatment_timing=sum(
            1 for row in grouping_rows if row["has_any_treatment_timing"] == "yes"
        ),
        patients_with_compound_raw_therapy_type_present=sum(
            1
            for row in grouping_rows
            if json_list_contains(row, "treatment_profile_flags_json", "compound_raw_therapy_type_present")
        ),
    )


def determine_conflict_type(grouping_row: dict[str, str]) -> str:
    if grouping_row["treatment_group_v1"] == "missing_type_only":
        return "missing_type_only"
    if grouping_row["treatment_group_v1"] == "mixed_multi_type":
        if grouping_row["dominant_therapy_type_rule"] == "blank_due_to_tie":
            return "mixed_multi_type_tied_dominant"
        if str(grouping_row["dominant_therapy_type_if_any"]).strip():
            return "mixed_multi_type_with_provisional_dominant"
        return "mixed_multi_type_without_single_dominant"
    raise PatientTreatmentGroupingV1Error(
        "Conflict audit requested for a row outside the v1 manual-review scope: "
        f"{grouping_row['bcr_patient_barcode']}"
    )


def build_source_evidence_summary(grouping_row: dict[str, str]) -> str:
    parts = [
        f"treatment_group_v1={grouping_row['treatment_group_v1']}",
        f"drug_row_count={grouping_row['drug_row_count']}",
        f"radiation_row_count={grouping_row['radiation_row_count']}",
        f"drug_therapy_type_single_or_mixed={grouping_row['drug_therapy_type_single_or_mixed']}",
        f"drug_therapy_type_values_json={grouping_row['drug_therapy_type_values_json']}",
        f"dominant_therapy_type_if_any={grouping_row['dominant_therapy_type_if_any'] or '[blank]'}",
        f"regimen_context_values_json={grouping_row['regimen_context_values_json']}",
        f"treatment_profile_flags_json={grouping_row['treatment_profile_flags_json']}",
        f"grouping_flags_json={grouping_row['grouping_flags_json']}",
    ]
    return "; ".join(parts)


def recommended_handling_for_conflict(conflict_type: str) -> tuple[str, str]:
    if conflict_type == "missing_type_only":
        return (
            "Keep the patient in missing_type_only, inspect the source-backed therapy-type documentation, and do not infer untreated status or assign a single-type arm.",
            "high",
        )
    if conflict_type == "mixed_multi_type_tied_dominant":
        return (
            "Keep the patient in mixed_multi_type, review the upstream mixed therapy evidence, and do not collapse to a dominant single-type arm before arm-freeze review.",
            "high",
        )
    if conflict_type == "mixed_multi_type_with_provisional_dominant":
        return (
            "Keep the patient in mixed_multi_type for grouped descriptive review, retain the provisional dominant raw type only as visible evidence, and do not reclassify to a single-type arm.",
            "medium",
        )
    if conflict_type == "mixed_multi_type_without_single_dominant":
        return (
            "Keep the patient in mixed_multi_type, review the saved profile evidence, and defer any tighter arm assignment to later manual review.",
            "high",
        )
    raise PatientTreatmentGroupingV1Error(f"Unhandled conflict_type: {conflict_type}")


def build_conflict_audit_rows(grouping_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for grouping_row in grouping_rows:
        if grouping_row["treatment_group_v1_requires_manual_review"] != "yes":
            continue
        conflict_type = determine_conflict_type(grouping_row)
        recommended_next_handling, review_priority = recommended_handling_for_conflict(conflict_type)
        rows.append(
            {
                "patient_treatment_grouping_v1_run_id": grouping_row["patient_treatment_grouping_v1_run_id"],
                "patient_treatment_profile_v1_run_id": grouping_row["patient_treatment_profile_v1_run_id"],
                "treatment_os_overlap_v1_run_id": grouping_row["treatment_os_overlap_v1_run_id"],
                "os_endpoint_v1_run_id": grouping_row["os_endpoint_v1_run_id"],
                "baseline_model_input_v1_run_id": grouping_row["baseline_model_input_v1_run_id"],
                "cohort_v1_build_id": grouping_row["cohort_v1_build_id"],
                "bcr_patient_barcode": grouping_row["bcr_patient_barcode"],
                "bcr_patient_uuid": grouping_row["bcr_patient_uuid"],
                "provisional_patient_row_id": grouping_row["provisional_patient_row_id"],
                "baseline_analysis_v1_row_id": grouping_row["baseline_analysis_v1_row_id"],
                "feature_set_v1_row_index": grouping_row["feature_set_v1_row_index"],
                "treatment_group_v1": grouping_row["treatment_group_v1"],
                "treatment_group_v1_rule": grouping_row["treatment_group_v1_rule"],
                "conflict_type": conflict_type,
                "source_evidence_summary": build_source_evidence_summary(grouping_row),
                "recommended_next_handling": recommended_next_handling,
                "review_priority": review_priority,
                "drug_therapy_type_single_or_mixed": grouping_row["drug_therapy_type_single_or_mixed"],
                "dominant_therapy_type_if_any": grouping_row["dominant_therapy_type_if_any"],
                "dominant_therapy_type_rule": grouping_row["dominant_therapy_type_rule"],
                "regimen_context_values_json": grouping_row["regimen_context_values_json"],
                "has_any_treatment_timing": grouping_row["has_any_treatment_timing"],
                "treatment_profile_status": grouping_row["treatment_profile_status"],
                "treatment_profile_flags_json": grouping_row["treatment_profile_flags_json"],
                "grouping_flags_json": grouping_row["grouping_flags_json"],
                "treatment_group_v1_requires_manual_review": grouping_row[
                    "treatment_group_v1_requires_manual_review"
                ],
            }
        )
    return rows


def build_arm_summary_rows(
    run_id: str,
    grouping_rows: list[dict[str, str]],
    metrics: GroupingMetrics,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for group_order, group_name in enumerate(APPROVED_TREATMENT_GROUP_ORDER, start=1):
        patient_count = metrics.group_counts[group_name]
        rows.append(
            {
                "patient_treatment_grouping_v1_run_id": run_id,
                "group_order": str(group_order),
                "treatment_group_v1": group_name,
                "patient_count": str(patient_count),
                "patient_fraction_of_cohort": format_fraction(patient_count, metrics.total),
                "os_event_count": str(metrics.group_event_counts[group_name]),
                "os_event_fraction_within_group": format_fraction(
                    metrics.group_event_counts[group_name],
                    patient_count,
                ),
                "manual_review_count": str(metrics.group_manual_review_counts[group_name]),
                "manual_review_fraction_within_group": format_fraction(
                    metrics.group_manual_review_counts[group_name],
                    patient_count,
                ),
                "radiation_overlap_count": str(metrics.group_radiation_counts[group_name]),
                "radiation_overlap_fraction_within_group": format_fraction(
                    metrics.group_radiation_counts[group_name],
                    patient_count,
                ),
                "regimen_context_count": str(metrics.group_regimen_context_counts[group_name]),
                "regimen_context_fraction_within_group": format_fraction(
                    metrics.group_regimen_context_counts[group_name],
                    patient_count,
                ),
                "treatment_timing_count": str(metrics.group_timing_counts[group_name]),
                "treatment_timing_fraction_within_group": format_fraction(
                    metrics.group_timing_counts[group_name],
                    patient_count,
                ),
                "notes": (
                    "No patients mapped to this provisional treatment group in the reviewed run."
                    if patient_count == 0
                    else ""
                ),
            }
        )
    return rows


def build_spec_rows(run_id: str) -> list[dict[str, str]]:
    spec_meta: dict[str, dict[str, str]] = {
        "patient_treatment_grouping_v1_run_id": {
            "source_origin": "workflow_generated",
            "field_category": "provenance",
            "decision_rule": "UTC run identifier assigned once per grouping workflow execution.",
            "notes": "Written to every patient-level grouping row for audit linkage.",
        },
        "patient_treatment_profile_v1_run_id": {
            "source_origin": "patient_treatment_profile_v1_latest_pointer",
            "field_category": "provenance",
            "decision_rule": "Carried forward from the saved patient_treatment_profile_v1 pointer.",
            "notes": "Links the grouping layer directly to the upstream patient-level treatment profile.",
        },
        "treatment_os_overlap_v1_run_id": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "provenance",
            "decision_rule": "Carried forward unchanged from the saved patient treatment profile row.",
            "notes": "Preserves linkage to the earlier treatment-OS overlap audit.",
        },
        "os_endpoint_v1_run_id": {
            "source_origin": "patient_treatment_profile_v1.tsv|os_endpoint_v1_latest_pointer",
            "field_category": "provenance",
            "decision_rule": "Must reconcile across the saved patient treatment profile pointer and OS endpoint pointer.",
            "notes": "Identifies the saved OS cohort used as the unit-of-analysis backbone.",
        },
        "baseline_model_input_v1_run_id": {
            "source_origin": "patient_treatment_profile_v1.tsv|baseline_model_input_v1_latest_pointer",
            "field_category": "provenance",
            "decision_rule": "Carried forward from the upstream saved profile and validated against the latest baseline model-input pointer.",
            "notes": "Used for lineage confirmation only, not as a grouping input matrix.",
        },
        "cohort_v1_build_id": {
            "source_origin": "shared upstream pointers",
            "field_category": "provenance",
            "decision_rule": "Must reconcile across the saved profile, OS endpoint, and baseline model-input pointers.",
            "notes": "Represents the shared OS cohort build carried through grouping v1.",
        },
        "bcr_patient_barcode": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row and used as the patient join key.",
            "notes": "Exactly one output row is written per OS cohort patient barcode.",
        },
        "bcr_patient_uuid": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row after UUID reconciliation against the upstream treatment profile row.",
            "notes": "Preserves patient UUID-level audit linkage.",
        },
        "provisional_patient_row_id": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row in upstream row order.",
            "notes": "Preserves linkage to prior patient-level audit layers.",
        },
        "baseline_analysis_v1_row_id": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied unchanged from the saved OS endpoint row.",
            "notes": "Supports downstream joins back to prior baseline analysis outputs.",
        },
        "feature_set_v1_row_index": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied unchanged from the saved OS endpoint row.",
            "notes": "Supports linkage to earlier feature-set and model-input outputs.",
        },
        "os_event": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "audit_outcome",
            "decision_rule": "Copied from the saved OS endpoint row for audit convenience only.",
            "notes": "Included for descriptive audit joins and not to define an analysis-ready treatment dataset.",
        },
        "os_time_days": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "audit_outcome",
            "decision_rule": "Copied from the saved OS endpoint row for audit convenience only.",
            "notes": "Included for descriptive audit joins and not to define an analysis-ready treatment dataset.",
        },
        "treatment_group_v1": {
            "source_origin": "derived_from_patient_treatment_profile_v1",
            "field_category": "grouping",
            "decision_rule": "Allowed values are restricted exactly to the approved v1 provisional grouping categories.",
            "notes": "Mixed cases remain mixed and are not collapsed into dominant-type single arms.",
        },
        "treatment_group_v1_rule": {
            "source_origin": "workflow_generated",
            "field_category": "grouping",
            "decision_rule": "Stores the exact rule path used to assign treatment_group_v1 from the saved treatment profile fields.",
            "notes": "Makes the provisional grouping logic explicit for audit review.",
        },
        "treatment_group_v1_requires_manual_review": {
            "source_origin": "patient_treatment_profile_v1.tsv.treatment_profile_requires_manual_review",
            "field_category": "review",
            "decision_rule": "Mirrors the upstream narrow-core manual-review flag exactly without expanding review scope in grouping v1.",
            "notes": "Preserves the current 357-patient manual-review burden.",
        },
        "has_any_drug_row": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "coverage",
            "decision_rule": "Carried forward unchanged from the saved patient treatment profile row.",
            "notes": "no_drug_record remains explicit and is not interpreted as untreated.",
        },
        "has_any_radiation_row": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "coverage",
            "decision_rule": "Carried forward unchanged from the saved patient treatment profile row.",
            "notes": "Radiation overlap stays visible as audit context rather than a separate final arm.",
        },
        "drug_row_count": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "coverage",
            "decision_rule": "Carried forward unchanged from the saved patient treatment profile row.",
            "notes": "Counts underlying drug-row evidence without re-reading raw tables.",
        },
        "radiation_row_count": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "coverage",
            "decision_rule": "Carried forward unchanged from the saved patient treatment profile row.",
            "notes": "Counts underlying radiation-row evidence without re-reading raw tables.",
        },
        "drug_therapy_type_values_json": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "therapy_type",
            "decision_rule": "Carried forward as the saved raw therapy-type evidence array from the patient treatment profile.",
            "notes": "Raw therapy-type labels remain visible, unsplit, and unnormalized in grouping v1.",
        },
        "drug_therapy_type_single_or_mixed": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "therapy_type",
            "decision_rule": "Uses the upstream patient-level therapy-type structure label: no_drug_record, single_type, mixed_type, or missing_type_only.",
            "notes": "Grouping v1 depends on the saved profile rather than re-inferring row-by-row treatment structure.",
        },
        "dominant_therapy_type_if_any": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "therapy_type",
            "decision_rule": "Carries forward the exact raw dominant therapy type from the saved profile when present.",
            "notes": "Compound raw values such as Chemotherapy|Hormone Therapy stay unsplit and may fall into single_ancillary_or_other.",
        },
        "dominant_therapy_type_rule": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "therapy_type",
            "decision_rule": "Carries forward the saved upstream explanation of why dominant_therapy_type_if_any is populated or blank.",
            "notes": "Used to distinguish tied-dominant mixed cases in the conflict audit.",
        },
        "regimen_context_values_json": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "regimen_context",
            "decision_rule": "Carried forward as the saved distinct regimen-context values JSON array.",
            "notes": "Regimen context is visible as coverage evidence only and does not define final arms here.",
        },
        "has_any_treatment_timing": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "timing",
            "decision_rule": "Carried forward unchanged from the saved patient treatment profile row.",
            "notes": "Summarizes presence of any aggregated drug or radiation timing evidence.",
        },
        "treatment_profile_status": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "review",
            "decision_rule": "Carries forward the upstream patient-treatment-profile status without modification.",
            "notes": "Maintains the explicit upstream no-drug, single-type, mixed-type, and missing-type-only distinctions.",
        },
        "treatment_profile_requires_manual_review": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "review",
            "decision_rule": "Carries forward the saved upstream manual-review flag for direct comparison with treatment_group_v1_requires_manual_review.",
            "notes": "Used to validate that grouping v1 keeps the upstream review scope unchanged.",
        },
        "treatment_profile_flags_json": {
            "source_origin": "patient_treatment_profile_v1.tsv",
            "field_category": "review",
            "decision_rule": "Carries forward the ordered upstream structural treatment-profile flags as saved.",
            "notes": "Retains structural evidence such as compound raw therapy types and timing inversion flags.",
        },
        "grouping_flags_json": {
            "source_origin": "workflow_generated",
            "field_category": "review",
            "decision_rule": "Stores ordered grouping-relevant flags derived from the saved profile evidence and provisional grouping assignment.",
            "notes": "Includes manual-review carry-forward, mixed/missing structure, ancillary fallback, compound raw labels, radiation overlap, regimen context, timing, no-drug-with-radiation, and timing inversion when present.",
        },
    }

    missing_meta_fields = [field_name for field_name in GROUPING_FIELDNAMES if field_name not in spec_meta]
    if missing_meta_fields:
        raise PatientTreatmentGroupingV1Error(
            f"Spec metadata missing for grouping fields: {missing_meta_fields}"
        )

    return [
        {
            "patient_treatment_grouping_v1_run_id": run_id,
            "field_name": field_name,
            "source_origin": spec_meta[field_name]["source_origin"],
            "field_category": spec_meta[field_name]["field_category"],
            "decision_rule": spec_meta[field_name]["decision_rule"],
            "notes": spec_meta[field_name]["notes"],
        }
        for field_name in GROUPING_FIELDNAMES
    ]


def build_summary_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
    metrics: GroupingMetrics,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    def add(summary_section: str, summary_metric: str, summary_value: Any, notes: str = "") -> None:
        rows.append(
            {
                "patient_treatment_grouping_v1_run_id": run_id,
                "summary_section": summary_section,
                "summary_metric": summary_metric,
                "summary_value": str(summary_value),
                "notes": notes,
            }
        )

    add(
        "inputs",
        "patient_treatment_profile_v1_run_id",
        workflow_inputs.patient_treatment_profile_latest_pointer["patient_treatment_profile_v1_run_id"],
    )
    add(
        "inputs",
        "treatment_os_overlap_v1_run_id",
        workflow_inputs.patient_treatment_profile_latest_pointer["treatment_os_overlap_v1_run_id"],
    )
    add("inputs", "os_endpoint_v1_run_id", workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"])
    add(
        "inputs",
        "baseline_model_input_v1_run_id",
        workflow_inputs.baseline_model_input_v1_latest_pointer["baseline_model_input_v1_run_id"],
    )
    add("inputs", "cohort_v1_build_id", workflow_inputs.os_endpoint_v1_latest_pointer["cohort_v1_build_id"])

    add("coverage", "total_os_cohort_size", metrics.total)
    add("coverage", "patients_requiring_manual_review", metrics.patients_requiring_manual_review)
    add(
        "coverage",
        "patients_requiring_manual_review_fraction",
        format_fraction(metrics.patients_requiring_manual_review, metrics.total),
        "Manual-review scope mirrors patient_treatment_profile_v1 exactly in grouping v1.",
    )
    add("coverage", "patients_with_any_radiation", metrics.patients_with_any_radiation)
    add(
        "coverage",
        "patients_with_any_radiation_fraction",
        format_fraction(metrics.patients_with_any_radiation, metrics.total),
    )
    add("coverage", "patients_with_any_regimen_context", metrics.patients_with_any_regimen_context)
    add(
        "coverage",
        "patients_with_any_regimen_context_fraction",
        format_fraction(metrics.patients_with_any_regimen_context, metrics.total),
    )
    add("coverage", "patients_with_any_treatment_timing", metrics.patients_with_any_treatment_timing)
    add(
        "coverage",
        "patients_with_any_treatment_timing_fraction",
        format_fraction(metrics.patients_with_any_treatment_timing, metrics.total),
    )
    add(
        "coverage",
        "patients_with_compound_raw_therapy_type_present",
        metrics.patients_with_compound_raw_therapy_type_present,
        "Compound raw therapy-type labels remain preserved exactly and are not split in grouping v1.",
    )
    add(
        "coverage",
        "patients_with_compound_raw_therapy_type_present_fraction",
        format_fraction(metrics.patients_with_compound_raw_therapy_type_present, metrics.total),
    )

    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        add("treatment_group", f"{group_name}_count", metrics.group_counts[group_name])
        add(
            "treatment_group",
            f"{group_name}_fraction",
            format_fraction(metrics.group_counts[group_name], metrics.total),
        )

    add(
        "readiness",
        "provisional_readiness_interpretation",
        READINESS_GROUPED_DESCRIPTIVE_REVIEW,
        "The saved patient-level grouping table is cohort-complete and suitable for grouped descriptive review only.",
    )
    add(
        "readiness",
        "treatment_arm_freeze_review_status",
        TREATMENT_ARM_FREEZE_REVIEW_STATUS,
        "Final treatment-arm freeze remains blocked because drug names remain unnormalized and manual review remains unresolved.",
    )
    add(
        "readiness",
        "remaining_blocker",
        "drug_name_not_normalized_and_manual_review_remaining",
        "This workflow remains provisional grouping/prep only.",
    )
    add(
        "readiness",
        "treatment_recommendation_modeling_status",
        TREATMENT_RECOMMENDATION_MODELING_STATUS,
        "Treatment recommendation modeling remains out of scope for this phase.",
    )
    add(
        "readiness",
        "grouping_scope",
        "grouping_prep_only",
        "This workflow does not normalize drug names, freeze treatment arms, or perform modeling.",
    )

    return rows


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_config = (
        repo_root / "09-trials" / "01-tcga-only-source-audited" / "04-config" / "trial_config.yaml"
    )
    if not trial_config.exists():
        raise PatientTreatmentGroupingV1Error(f"Required trial config not found: {trial_config}")

    trial_config_data = helper_module.load_yaml(trial_config)
    processed_root = repo_root / str(trial_config_data.get("processed_data_root", "01-data/processed"))
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    treatment_prep_root = audit_root / "tcga-brca" / "treatment-prep"
    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        processed_runs_root=processed_root / "tcga-brca" / "treatment-prep" / "patient_treatment_grouping_v1_runs",
        audit_runs_root=treatment_prep_root / "patient_treatment_grouping_v1_runs",
        latest_pointer=treatment_prep_root / "tcga_brca_patient_treatment_grouping_v1_latest.json",
        patient_treatment_profile_latest_pointer=treatment_prep_root
        / "tcga_brca_patient_treatment_profile_v1_latest.json",
        os_endpoint_v1_latest_pointer=audit_root
        / "tcga-brca"
        / "endpoint-prep"
        / "tcga_brca_os_endpoint_v1_latest.json",
        baseline_model_input_v1_latest_pointer=audit_root
        / "tcga-brca"
        / "model-input"
        / "tcga_brca_baseline_model_input_v1_latest.json",
    )


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    for pointer_path, label in [
        (paths.patient_treatment_profile_latest_pointer, "patient treatment profile v1 latest pointer"),
        (paths.os_endpoint_v1_latest_pointer, "OS endpoint v1 latest pointer"),
        (paths.baseline_model_input_v1_latest_pointer, "baseline model-input v1 latest pointer"),
    ]:
        if not pointer_path.exists():
            raise PatientTreatmentGroupingV1Error(f"Required {label} not found: {pointer_path}")

    profile_pointer = helper_module.load_json(paths.patient_treatment_profile_latest_pointer)
    helper_module.require_keys(
        profile_pointer,
        {
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "patient_treatment_profile_v1_tsv",
            "run_log_json",
        },
        "patient treatment profile v1 latest pointer",
        paths.patient_treatment_profile_latest_pointer,
    )
    os_pointer = helper_module.load_json(paths.os_endpoint_v1_latest_pointer)
    helper_module.require_keys(
        os_pointer,
        {
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "os_endpoint_v1_tsv",
            "run_log_json",
        },
        "OS endpoint v1 latest pointer",
        paths.os_endpoint_v1_latest_pointer,
    )
    baseline_pointer = helper_module.load_json(paths.baseline_model_input_v1_latest_pointer)
    helper_module.require_keys(
        baseline_pointer,
        {
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "run_log_json",
        },
        "baseline model-input v1 latest pointer",
        paths.baseline_model_input_v1_latest_pointer,
    )

    profile_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(profile_pointer["run_log_json"]),
        label="patient treatment profile v1 run log",
        helper_module=helper_module,
    )
    os_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(os_pointer["run_log_json"]),
        label="OS endpoint v1 run log",
        helper_module=helper_module,
    )
    baseline_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(baseline_pointer["run_log_json"]),
        label="baseline model-input v1 run log",
        helper_module=helper_module,
    )

    if str(profile_pointer["os_endpoint_v1_run_id"]) != str(os_pointer["os_endpoint_v1_run_id"]):
        raise PatientTreatmentGroupingV1Error(
            "Mismatch between patient treatment profile and OS endpoint pointers for os_endpoint_v1_run_id."
        )
    if str(profile_pointer["baseline_model_input_v1_run_id"]) != str(
        os_pointer["baseline_model_input_v1_run_id"]
    ):
        raise PatientTreatmentGroupingV1Error(
            "Mismatch between patient treatment profile and OS endpoint pointers for baseline_model_input_v1_run_id."
        )
    if str(profile_pointer["baseline_model_input_v1_run_id"]) != str(
        baseline_pointer["baseline_model_input_v1_run_id"]
    ):
        raise PatientTreatmentGroupingV1Error(
            "Mismatch between patient treatment profile and baseline model-input pointers for baseline_model_input_v1_run_id."
        )
    if str(profile_pointer["cohort_v1_build_id"]) != str(os_pointer["cohort_v1_build_id"]):
        raise PatientTreatmentGroupingV1Error(
            "Mismatch between patient treatment profile and OS endpoint pointers for cohort_v1_build_id."
        )
    if str(profile_pointer["cohort_v1_build_id"]) != str(baseline_pointer["cohort_v1_build_id"]):
        raise PatientTreatmentGroupingV1Error(
            "Mismatch between patient treatment profile and baseline model-input pointers for cohort_v1_build_id."
        )

    profile_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(profile_pointer["patient_treatment_profile_v1_tsv"]),
        "patient_treatment_profile_v1.tsv",
    )
    os_endpoint_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(os_pointer["os_endpoint_v1_tsv"]),
        "os_endpoint_v1.tsv",
    )

    profile_rows = helper_module.read_tsv_dict_rows(profile_tsv)
    os_endpoint_rows = helper_module.read_tsv_dict_rows(os_endpoint_tsv)

    require_columns(
        profile_rows,
        {
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
            "has_any_drug_row",
            "has_any_radiation_row",
            "drug_row_count",
            "radiation_row_count",
            "drug_therapy_type_values_json",
            "drug_therapy_type_single_or_mixed",
            "dominant_therapy_type_if_any",
            "dominant_therapy_type_rule",
            "regimen_context_values_json",
            "has_any_treatment_timing",
            "treatment_profile_status",
            "treatment_profile_requires_manual_review",
            "treatment_profile_flags_json",
        },
        "patient_treatment_profile_v1.tsv",
    )
    require_columns(
        os_endpoint_rows,
        {
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

    if not profile_rows:
        raise PatientTreatmentGroupingV1Error("patient_treatment_profile_v1.tsv has no data rows.")
    if not os_endpoint_rows:
        raise PatientTreatmentGroupingV1Error("os_endpoint_v1.tsv has no data rows.")

    return WorkflowInputs(
        patient_treatment_profile_latest_pointer=profile_pointer,
        patient_treatment_profile_run_log=profile_run_log,
        os_endpoint_v1_latest_pointer=os_pointer,
        os_endpoint_v1_run_log=os_run_log,
        baseline_model_input_v1_latest_pointer=baseline_pointer,
        baseline_model_input_v1_run_log=baseline_run_log,
        patient_treatment_profile_rows=profile_rows,
        os_endpoint_rows=os_endpoint_rows,
        input_paths={
            "patient_treatment_profile_v1_tsv": profile_tsv,
            "os_endpoint_v1_tsv": os_endpoint_tsv,
        },
    )


def validate_outputs(
    workflow_inputs: WorkflowInputs,
    grouping_rows: list[dict[str, str]],
    spec_rows: list[dict[str, str]],
    arm_summary_rows: list[dict[str, str]],
    conflict_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    metrics: GroupingMetrics,
) -> dict[str, Any]:
    os_barcodes = [normalize_barcode(row["bcr_patient_barcode"]) for row in workflow_inputs.os_endpoint_rows]
    grouping_barcodes = [normalize_barcode(row["bcr_patient_barcode"]) for row in grouping_rows]
    profile_lookup = build_unique_lookup_by_barcode(
        workflow_inputs.patient_treatment_profile_rows,
        barcode_field="bcr_patient_barcode",
        label="patient_treatment_profile_v1.tsv",
    )

    row_order_preserved = grouping_barcodes == os_barcodes
    all_os_patients_in_grouping = set(grouping_barcodes) == set(os_barcodes)
    no_grouping_duplicate_barcodes = len(grouping_barcodes) == len(set(grouping_barcodes))
    group_categories_restricted = all(
        row["treatment_group_v1"] in APPROVED_TREATMENT_GROUP_ORDER for row in grouping_rows
    )
    manual_review_scope_matches_upstream = all(
        row["treatment_group_v1_requires_manual_review"] == row["treatment_profile_requires_manual_review"]
        for row in grouping_rows
    )

    profile_continuity_grouping_vs_profile = True
    for grouping_row in grouping_rows:
        profile_row = profile_lookup.get(normalize_barcode(grouping_row["bcr_patient_barcode"]))
        if profile_row is None:
            profile_continuity_grouping_vs_profile = False
            break
        for field_name in [
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "has_any_drug_row",
            "has_any_radiation_row",
            "drug_row_count",
            "radiation_row_count",
            "drug_therapy_type_values_json",
            "drug_therapy_type_single_or_mixed",
            "dominant_therapy_type_if_any",
            "dominant_therapy_type_rule",
            "regimen_context_values_json",
            "has_any_treatment_timing",
            "treatment_profile_status",
            "treatment_profile_requires_manual_review",
            "treatment_profile_flags_json",
        ]:
            if grouping_row[field_name] != profile_row[field_name]:
                profile_continuity_grouping_vs_profile = False
                break
        if grouping_row["treatment_group_v1_requires_manual_review"] != profile_row[
            "treatment_profile_requires_manual_review"
        ]:
            profile_continuity_grouping_vs_profile = False
        if not profile_continuity_grouping_vs_profile:
            break

    group_counts_from_rows = Counter(row["treatment_group_v1"] for row in grouping_rows)
    group_counts_reconcile = all(
        group_counts_from_rows.get(group_name, 0) == metrics.group_counts[group_name]
        for group_name in APPROVED_TREATMENT_GROUP_ORDER
    )
    patient_counts_reconcile_to_os_cohort = sum(metrics.group_counts.values()) == metrics.total == len(
        workflow_inputs.os_endpoint_rows
    )

    arm_summary_lookup = {row["treatment_group_v1"]: row for row in arm_summary_rows}
    arm_summary_counts_reconcile = len(arm_summary_rows) == len(APPROVED_TREATMENT_GROUP_ORDER) and all(
        group_name in arm_summary_lookup
        and parse_int_or_none(arm_summary_lookup[group_name]["patient_count"]) == metrics.group_counts[group_name]
        and parse_int_or_none(arm_summary_lookup[group_name]["os_event_count"])
        == metrics.group_event_counts[group_name]
        and parse_int_or_none(arm_summary_lookup[group_name]["manual_review_count"])
        == metrics.group_manual_review_counts[group_name]
        and parse_int_or_none(arm_summary_lookup[group_name]["radiation_overlap_count"])
        == metrics.group_radiation_counts[group_name]
        and parse_int_or_none(arm_summary_lookup[group_name]["regimen_context_count"])
        == metrics.group_regimen_context_counts[group_name]
        and parse_int_or_none(arm_summary_lookup[group_name]["treatment_timing_count"])
        == metrics.group_timing_counts[group_name]
        for group_name in APPROVED_TREATMENT_GROUP_ORDER
    )

    summary_metric_lookup = {row["summary_metric"]: row["summary_value"] for row in summary_rows}

    def summary_int(metric_name: str) -> int | None:
        return parse_int_or_none(summary_metric_lookup.get(metric_name, ""))

    summary_counts_reconcile = (
        summary_int("total_os_cohort_size") == metrics.total
        and summary_int("patients_requiring_manual_review") == metrics.patients_requiring_manual_review
        and summary_int("patients_with_any_radiation") == metrics.patients_with_any_radiation
        and summary_int("patients_with_any_regimen_context") == metrics.patients_with_any_regimen_context
        and summary_int("patients_with_any_treatment_timing") == metrics.patients_with_any_treatment_timing
        and summary_int("patients_with_compound_raw_therapy_type_present")
        == metrics.patients_with_compound_raw_therapy_type_present
        and all(
            summary_int(f"{group_name}_count") == metrics.group_counts[group_name]
            for group_name in APPROVED_TREATMENT_GROUP_ORDER
        )
    )

    readiness_allowed = (
        summary_metric_lookup.get("provisional_readiness_interpretation", "")
        == READINESS_GROUPED_DESCRIPTIVE_REVIEW
        and summary_metric_lookup.get("treatment_arm_freeze_review_status", "")
        == TREATMENT_ARM_FREEZE_REVIEW_STATUS
        and summary_metric_lookup.get("treatment_recommendation_modeling_status", "")
        == TREATMENT_RECOMMENDATION_MODELING_STATUS
    )

    spec_field_names_match = [row["field_name"] for row in spec_rows] == GROUPING_FIELDNAMES
    conflict_audit_row_count_matches_manual_review = len(conflict_rows) == metrics.patients_requiring_manual_review

    passed = all(
        [
            len(workflow_inputs.patient_treatment_profile_rows) > 0,
            len(workflow_inputs.os_endpoint_rows) > 0,
            len(grouping_rows) > 0,
            len(grouping_rows) == len(workflow_inputs.os_endpoint_rows),
            row_order_preserved,
            all_os_patients_in_grouping,
            no_grouping_duplicate_barcodes,
            group_categories_restricted,
            manual_review_scope_matches_upstream,
            profile_continuity_grouping_vs_profile,
            patient_counts_reconcile_to_os_cohort,
            group_counts_reconcile,
            arm_summary_counts_reconcile,
            summary_counts_reconcile,
            conflict_audit_row_count_matches_manual_review,
            len(spec_rows) == len(GROUPING_FIELDNAMES),
            spec_field_names_match,
            readiness_allowed,
        ]
    )

    return {
        "passed": passed,
        "required_upstream_pointers_found": True,
        "required_source_tables_found": True,
        "patient_treatment_profile_run_log_completed": True,
        "os_endpoint_run_log_completed": True,
        "baseline_model_input_run_log_completed": True,
        "patient_treatment_profile_row_count_positive": len(workflow_inputs.patient_treatment_profile_rows) > 0,
        "os_endpoint_row_count_positive": len(workflow_inputs.os_endpoint_rows) > 0,
        "output_rows_positive": len(grouping_rows) > 0,
        "grouping_row_count_matches_os_cohort": len(grouping_rows) == len(workflow_inputs.os_endpoint_rows),
        "row_order_preserved_from_os_endpoint": row_order_preserved,
        "all_os_patients_in_grouping": all_os_patients_in_grouping,
        "no_grouping_duplicate_barcodes": no_grouping_duplicate_barcodes,
        "group_categories_restricted_to_approved_set": group_categories_restricted,
        "manual_review_scope_matches_upstream_profile": manual_review_scope_matches_upstream,
        "profile_continuity_grouping_vs_profile": profile_continuity_grouping_vs_profile,
        "patient_counts_reconcile_to_os_cohort": patient_counts_reconcile_to_os_cohort,
        "group_counts_reconcile_to_grouping_table": group_counts_reconcile,
        "arm_summary_counts_reconcile_to_grouping_table": arm_summary_counts_reconcile,
        "summary_counts_reconcile_to_grouping_table": summary_counts_reconcile,
        "conflict_audit_row_count_matches_manual_review": conflict_audit_row_count_matches_manual_review,
        "spec_row_count_matches_grouping_columns": len(spec_rows) == len(GROUPING_FIELDNAMES),
        "spec_field_names_match_grouping_columns": spec_field_names_match,
        "summary_readiness_interpretation_allowed": readiness_allowed,
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
        "os_cohort_row_count": len(workflow_inputs.os_endpoint_rows),
        "patient_treatment_profile_row_count": len(workflow_inputs.patient_treatment_profile_rows),
        "grouping_row_count": len(grouping_rows),
        "manual_review_row_count": len(conflict_rows),
        "spec_row_count": len(spec_rows),
        "arm_summary_row_count": len(arm_summary_rows),
        "summary_row_count": len(summary_rows),
    }


def build_latest_pointer_payload(
    *,
    run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "patient_treatment_grouping_v1_run_id": run_id,
        "patient_treatment_profile_v1_run_id": str(
            workflow_inputs.patient_treatment_profile_latest_pointer["patient_treatment_profile_v1_run_id"]
        ),
        "treatment_os_overlap_v1_run_id": str(
            workflow_inputs.patient_treatment_profile_latest_pointer["treatment_os_overlap_v1_run_id"]
        ),
        "os_endpoint_v1_run_id": str(workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"]),
        "baseline_model_input_v1_run_id": str(
            workflow_inputs.baseline_model_input_v1_latest_pointer["baseline_model_input_v1_run_id"]
        ),
        "cohort_v1_build_id": str(workflow_inputs.os_endpoint_v1_latest_pointer["cohort_v1_build_id"]),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "patient_treatment_grouping_v1_tsv": repo_relative(
            output_paths["patient_treatment_grouping_v1_tsv"],
            paths.repo_root,
        ),
        "patient_treatment_grouping_v1_spec_tsv": repo_relative(
            output_paths["patient_treatment_grouping_v1_spec_tsv"],
            paths.repo_root,
        ),
        "patient_treatment_grouping_v1_arm_summary_tsv": repo_relative(
            output_paths["patient_treatment_grouping_v1_arm_summary_tsv"],
            paths.repo_root,
        ),
        "patient_treatment_grouping_v1_conflict_audit_tsv": repo_relative(
            output_paths["patient_treatment_grouping_v1_conflict_audit_tsv"],
            paths.repo_root,
        ),
        "patient_treatment_grouping_v1_summary_tsv": repo_relative(
            output_paths["patient_treatment_grouping_v1_summary_tsv"],
            paths.repo_root,
        ),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "patient_treatment_profile_v1_latest_json": repo_relative(
            paths.patient_treatment_profile_latest_pointer,
            paths.repo_root,
        ),
        "os_endpoint_v1_latest_json": repo_relative(paths.os_endpoint_v1_latest_pointer, paths.repo_root),
        "baseline_model_input_v1_latest_json": repo_relative(
            paths.baseline_model_input_v1_latest_pointer,
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
    processed_run_dir = paths.processed_runs_root / run_id
    audit_run_dir = paths.audit_runs_root / run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)

        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        grouping_path = processed_run_dir / "patient_treatment_grouping_v1.tsv"
        spec_path = audit_run_dir / "patient_treatment_grouping_v1_spec.tsv"
        arm_summary_path = audit_run_dir / "patient_treatment_grouping_v1_arm_summary.tsv"
        conflict_audit_path = audit_run_dir / "patient_treatment_grouping_v1_conflict_audit.tsv"
        summary_path = audit_run_dir / "patient_treatment_grouping_v1_summary.tsv"

        grouping_rows, metrics = build_grouping_rows(run_id, workflow_inputs)
        spec_rows = build_spec_rows(run_id)
        arm_summary_rows = build_arm_summary_rows(run_id, grouping_rows, metrics)
        conflict_rows = build_conflict_audit_rows(grouping_rows)
        summary_rows = build_summary_rows(run_id, workflow_inputs, metrics)

        validation = validate_outputs(
            workflow_inputs=workflow_inputs,
            grouping_rows=grouping_rows,
            spec_rows=spec_rows,
            arm_summary_rows=arm_summary_rows,
            conflict_rows=conflict_rows,
            summary_rows=summary_rows,
            metrics=metrics,
        )
        if not validation["passed"]:
            raise PatientTreatmentGroupingV1Error(
                "Output validation did not pass. Failed checks: "
                + str({key: value for key, value in validation.items() if value is False})
            )

        helper_module.write_dict_rows_tsv(grouping_path, GROUPING_FIELDNAMES, grouping_rows)
        helper_module.write_dict_rows_tsv(spec_path, SPEC_FIELDNAMES, spec_rows)
        helper_module.write_dict_rows_tsv(arm_summary_path, ARM_SUMMARY_FIELDNAMES, arm_summary_rows)
        helper_module.write_dict_rows_tsv(conflict_audit_path, CONFLICT_AUDIT_FIELDNAMES, conflict_rows)
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "patient_treatment_grouping_v1_tsv": grouping_path,
            "patient_treatment_grouping_v1_spec_tsv": spec_path,
            "patient_treatment_grouping_v1_arm_summary_tsv": arm_summary_path,
            "patient_treatment_grouping_v1_conflict_audit_tsv": conflict_audit_path,
            "patient_treatment_grouping_v1_summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }

        latest_pointer_payload = build_latest_pointer_payload(
            run_id=run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload: dict[str, Any] = {
            "status": "completed",
            "patient_treatment_grouping_v1_run_id": run_id,
            "patient_treatment_profile_v1_run_id": str(
                workflow_inputs.patient_treatment_profile_latest_pointer["patient_treatment_profile_v1_run_id"]
            ),
            "treatment_os_overlap_v1_run_id": str(
                workflow_inputs.patient_treatment_profile_latest_pointer["treatment_os_overlap_v1_run_id"]
            ),
            "os_endpoint_v1_run_id": str(workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"]),
            "baseline_model_input_v1_run_id": str(
                workflow_inputs.baseline_model_input_v1_latest_pointer["baseline_model_input_v1_run_id"]
            ),
            "cohort_v1_build_id": str(workflow_inputs.os_endpoint_v1_latest_pointer["cohort_v1_build_id"]),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "patient_treatment_profile_v1_latest_json": repo_relative(
                    paths.patient_treatment_profile_latest_pointer,
                    paths.repo_root,
                ),
                "os_endpoint_v1_latest_json": repo_relative(paths.os_endpoint_v1_latest_pointer, paths.repo_root),
                "baseline_model_input_v1_latest_json": repo_relative(
                    paths.baseline_model_input_v1_latest_pointer,
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
                "patient_treatment_grouping_v1_tsv": repo_relative(grouping_path, paths.repo_root),
                "patient_treatment_grouping_v1_spec_tsv": repo_relative(spec_path, paths.repo_root),
                "patient_treatment_grouping_v1_arm_summary_tsv": repo_relative(
                    arm_summary_path,
                    paths.repo_root,
                ),
                "patient_treatment_grouping_v1_conflict_audit_tsv": repo_relative(
                    conflict_audit_path,
                    paths.repo_root,
                ),
                "patient_treatment_grouping_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation,
            "rules": {
                "unit_of_analysis": "patient/case",
                "input_layer": "saved_patient_treatment_profile_v1_plus_saved_os_endpoint_v1",
                "output_layer": "patient_treatment_grouping_v1",
                "allowed_treatment_group_categories": APPROVED_TREATMENT_GROUP_ORDER,
                "manual_review_scope": "mirrored_from_patient_treatment_profile_v1_narrow_core",
                "mixed_group_handling": "retain_as_mixed_multi_type",
                "single_other_handling": "single_ancillary_or_other_exact_raw_dominant_type",
                "compound_raw_therapy_type_preserved": True,
                "baseline_model_input_used_for_lineage_confirmation_only": True,
                "no_drug_name_normalization": True,
                "no_treatment_arm_freeze": True,
                "no_modeling": True,
                "no_causal_analysis": True,
                "grouping_prep_only": True,
            },
            "counts": {
                "os_endpoint_row_count": len(workflow_inputs.os_endpoint_rows),
                "patient_treatment_profile_row_count": len(workflow_inputs.patient_treatment_profile_rows),
                "patient_treatment_grouping_row_count": len(grouping_rows),
                "patient_treatment_grouping_spec_row_count": len(spec_rows),
                "patient_treatment_grouping_arm_summary_row_count": len(arm_summary_rows),
                "patient_treatment_grouping_conflict_audit_row_count": len(conflict_rows),
                "patient_treatment_grouping_summary_row_count": len(summary_rows),
                "patients_requiring_manual_review": metrics.patients_requiring_manual_review,
                "patients_with_any_radiation": metrics.patients_with_any_radiation,
                "patients_with_any_regimen_context": metrics.patients_with_any_regimen_context,
                "patients_with_any_treatment_timing": metrics.patients_with_any_treatment_timing,
                "patients_with_compound_raw_therapy_type_present": (
                    metrics.patients_with_compound_raw_therapy_type_present
                ),
                **{f"group_{group_name}_count": metrics.group_counts[group_name] for group_name in APPROVED_TREATMENT_GROUP_ORDER},
                "provisional_readiness_interpretation": READINESS_GROUPED_DESCRIPTIVE_REVIEW,
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "patient_treatment_profile_v1_latest_pointer": workflow_inputs.patient_treatment_profile_latest_pointer,
                "os_endpoint_v1_latest_pointer": workflow_inputs.os_endpoint_v1_latest_pointer,
                "baseline_model_input_v1_latest_pointer": workflow_inputs.baseline_model_input_v1_latest_pointer,
            },
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        completed_at = utc_now()
        failure_payload = {
            "status": "failed",
            "patient_treatment_grouping_v1_run_id": run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(completed_at),
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
            },
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    counts = run_log.get("counts", {})
    print("TCGA-BRCA patient treatment grouping v1 completed")
    print(f"  Run ID                : {run_log['patient_treatment_grouping_v1_run_id']}")
    print(f"  Profile run ID        : {run_log['patient_treatment_profile_v1_run_id']}")
    print(f"  OS endpoint run ID    : {run_log['os_endpoint_v1_run_id']}")
    print(f"  Cohort size           : {counts.get('patient_treatment_grouping_row_count', 0)}")
    for group_name in APPROVED_TREATMENT_GROUP_ORDER:
        print(f"  {group_name:<21}: {counts.get(f'group_{group_name}_count', 0)}")
    print(f"  Manual review rows    : {counts.get('patients_requiring_manual_review', 0)}")
    print(f"  Processed run dir     : {run_log['outputs']['processed_run_directory']}")
    print(f"  Audit run dir         : {run_log['outputs']['audit_run_directory']}")
    print(f"  Latest pointer        : {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
