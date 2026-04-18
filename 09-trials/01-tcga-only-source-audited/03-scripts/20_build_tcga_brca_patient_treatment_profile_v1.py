#!/usr/bin/env python
"""Build an auditable patient-level TCGA-BRCA treatment profile v1 from saved treatment tables."""

from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

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
PROFILE_FIELDNAMES = [
    "patient_treatment_profile_v1_run_id",
    "treatment_os_overlap_v1_run_id",
    "os_endpoint_v1_run_id",
    "clinical_biotab_parse_run_id",
    "baseline_model_input_v1_run_id",
    "cohort_v1_build_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "provisional_patient_row_id",
    "baseline_analysis_v1_row_id",
    "feature_set_v1_row_index",
    "has_any_drug_row",
    "drug_row_count",
    "has_any_radiation_row",
    "radiation_row_count",
    "drug_therapy_type_values_json",
    "drug_therapy_type_distinct_count",
    "drug_therapy_type_single_or_mixed",
    "dominant_therapy_type_if_any",
    "dominant_therapy_type_rule",
    "therapy_type_counts_json",
    "has_any_regimen_context",
    "regimen_context_values_json",
    "regimen_context_distinct_count",
    "regimen_context_single_or_mixed",
    "drug_name_values_json",
    "distinct_drug_name_count",
    "earliest_drug_start_days",
    "latest_drug_end_days",
    "earliest_radiation_start_days",
    "latest_radiation_end_days",
    "has_any_treatment_timing",
    "treatment_profile_requires_manual_review",
    "treatment_profile_status",
    "treatment_profile_flags_json",
]
SPEC_FIELDNAMES = [
    "patient_treatment_profile_v1_run_id",
    "field_name",
    "source_origin",
    "field_category",
    "decision_rule",
    "notes",
]
SUMMARY_FIELDNAMES = [
    "patient_treatment_profile_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
CONFLICT_AUDIT_FIELDNAMES = [
    "patient_treatment_profile_v1_run_id",
    "treatment_os_overlap_v1_run_id",
    "os_endpoint_v1_run_id",
    "clinical_biotab_parse_run_id",
    "baseline_model_input_v1_run_id",
    "cohort_v1_build_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "provisional_patient_row_id",
    "baseline_analysis_v1_row_id",
    "feature_set_v1_row_index",
    "conflict_type",
    "source_evidence_summary",
    "recommended_next_handling",
    "review_priority",
    "drug_therapy_type_single_or_mixed",
    "dominant_therapy_type_if_any",
    "dominant_therapy_type_rule",
    "treatment_profile_status",
    "treatment_profile_flags_json",
    "treatment_profile_requires_manual_review",
]
READINESS_NOT_READY = "not_ready"
READINESS_COARSE_GROUPING = "ready_for_coarse_exploratory_grouping_next"
READINESS_ARM_FREEZE_REVIEW = "ready_for_treatment_arm_freeze_review"
STATUS_NO_DRUG_RECORD = "no_drug_record"
STATUS_SINGLE = "single_drug_therapy_type"
STATUS_MIXED = "mixed_drug_therapy_types"
STATUS_MISSING_TYPE_ONLY = "drug_rows_missing_therapy_type_only"
SINGLE_OR_MIXED_NO_DRUG = "no_drug_record"
SINGLE_OR_MIXED_SINGLE = "single_type"
SINGLE_OR_MIXED_MIXED = "mixed_type"
SINGLE_OR_MIXED_MISSING = "missing_type_only"


class PatientTreatmentProfileV1Error(RuntimeError):
    """Raised when the patient treatment profile v1 workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the patient treatment profile v1 workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    treatment_os_overlap_latest_pointer: Path
    os_endpoint_v1_latest_pointer: Path
    clinical_biotabs_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    treatment_os_overlap_latest_pointer: dict[str, Any]
    treatment_os_overlap_run_log: dict[str, Any]
    os_endpoint_v1_latest_pointer: dict[str, Any]
    os_endpoint_v1_run_log: dict[str, Any]
    clinical_biotabs_latest_pointer: dict[str, Any]
    clinical_biotabs_run_log: dict[str, Any]
    os_endpoint_rows: list[dict[str, str]]
    treatment_os_overlap_rows: list[dict[str, str]]
    clinical_drug_rows: list[dict[str, str]]
    clinical_radiation_rows: list[dict[str, str]]
    clinical_patient_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class ProfileMetrics:
    """Scalar metrics derived from the patient treatment profile table."""

    total: int
    patients_with_drug: int
    patients_without_drug: int
    patients_with_radiation: int
    patients_without_radiation: int
    patients_single_therapy_type: int
    patients_mixed_therapy_type: int
    patients_missing_therapy_type_only: int
    patients_with_dominant_therapy_type_identifiable: int
    patients_requiring_manual_review: int
    patients_with_any_regimen_context: int
    patients_with_any_treatment_timing: int
    patients_no_drug_or_radiation_rows: int
    patients_radiation_without_drug_row: int
    patients_compound_raw_therapy_type_present: int
    patients_drug_timing_window_aggregated_inverted: int


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise PatientTreatmentProfileV1Error(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise PatientTreatmentProfileV1Error(f"Unable to create an import spec for: {script_path}")
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


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def normalize_barcode(value: str) -> str:
    return value.strip().upper()


def normalize_uuid(value: str) -> str:
    return value.strip().upper()


def normalize_missing_like(value: str) -> str:
    return value.strip().lower()


def is_missing_like(value: str) -> bool:
    return normalize_missing_like(value) in MISSING_LIKE_TOKENS


def parse_int_or_none(value: str) -> int | None:
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return int(stripped)
    except ValueError:
        return None


def require_columns(rows: list[dict[str, str]], required_columns: set[str], label: str) -> None:
    if not rows:
        raise PatientTreatmentProfileV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise PatientTreatmentProfileV1Error(f"{label} is missing required columns: {sorted(missing)}")


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
        raise PatientTreatmentProfileV1Error(f"{label} is not completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise PatientTreatmentProfileV1Error(f"{label} does not report validation.passed == true.")
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
            raise PatientTreatmentProfileV1Error(f"{label} contains an empty barcode value.")
        if barcode in lookup:
            raise PatientTreatmentProfileV1Error(f"{label} contains a duplicate barcode: {barcode}")
        lookup[barcode] = row
    return lookup


def group_rows_by_barcode(
    rows: list[dict[str, str]],
    *,
    barcode_field: str,
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        barcode = normalize_barcode(str(row.get(barcode_field, "")))
        if barcode:
            grouped[barcode].append(row)
    return grouped


def collect_raw_values(rows: list[dict[str, str]], field_name: str) -> list[str]:
    return [str(row.get(field_name, "")) for row in rows]


def ordered_non_missing_distinct_values(raw_values: list[str]) -> list[str]:
    normalized_values = [value.strip() for value in raw_values if not is_missing_like(value.strip())]
    return ordered_unique(normalized_values)


def ordered_count_objects(raw_values: list[str]) -> list[dict[str, Any]]:
    counter: Counter[str] = Counter()
    ordered_values: list[str] = []
    for raw_value in raw_values:
        normalized_value = raw_value.strip()
        if is_missing_like(normalized_value):
            continue
        if normalized_value not in counter:
            ordered_values.append(normalized_value)
        counter[normalized_value] += 1
    return [{"value": value, "count": counter[value]} for value in ordered_values]


def collect_parseable_ints(rows: list[dict[str, str]], field_name: str) -> list[int]:
    values: list[int] = []
    for row in rows:
        raw_value = str(row.get(field_name, ""))
        if is_missing_like(raw_value):
            continue
        parsed = parse_int_or_none(raw_value)
        if parsed is not None:
            values.append(parsed)
    return values


def collect_regimen_context_values(
    drug_rows: list[dict[str, str]],
    radiation_rows: list[dict[str, str]],
) -> list[str]:
    raw_values: list[str] = []
    for row in drug_rows:
        raw_values.append(str(row.get("therapy_regimen", "")))
        raw_values.append(str(row.get("pharm_regimen", "")))
    for row in radiation_rows:
        raw_values.append(str(row.get("therapy_regimen", "")))
    return ordered_non_missing_distinct_values(raw_values)


def single_or_mixed_from_overlap(row: dict[str, str]) -> str:
    has_drug_row = str(row.get("has_drug_row", "")).strip().lower()
    overlap_value = str(row.get("drug_therapy_type_single_or_mixed", "")).strip().lower()
    if has_drug_row == "no":
        return SINGLE_OR_MIXED_NO_DRUG
    if overlap_value == "single":
        return SINGLE_OR_MIXED_SINGLE
    if overlap_value == "mixed":
        return SINGLE_OR_MIXED_MIXED
    if overlap_value == "missing":
        return SINGLE_OR_MIXED_MISSING
    raise PatientTreatmentProfileV1Error(
        "Unable to map overlap drug_therapy_type_single_or_mixed value "
        f"{row.get('drug_therapy_type_single_or_mixed', '')!r} for barcode {row.get('bcr_patient_barcode', '')!r}."
    )


def json_list_contains(row: dict[str, str], field_name: str, expected_value: str) -> bool:
    raw_value = str(row.get(field_name, "")).strip()
    if not raw_value:
        return False
    parsed = json.loads(raw_value)
    if not isinstance(parsed, list):
        raise PatientTreatmentProfileV1Error(f"Expected JSON list in {field_name}: {raw_value!r}")
    return expected_value in parsed


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_config = (
        repo_root / "09-trials" / "01-tcga-only-source-audited" / "04-config" / "trial_config.yaml"
    )
    if not trial_config.exists():
        raise PatientTreatmentProfileV1Error(f"Required trial config not found: {trial_config}")
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
        processed_runs_root=processed_root / "tcga-brca" / "treatment-prep" / "patient_treatment_profile_v1_runs",
        audit_runs_root=treatment_prep_root / "patient_treatment_profile_v1_runs",
        latest_pointer=treatment_prep_root / "tcga_brca_patient_treatment_profile_v1_latest.json",
        treatment_os_overlap_latest_pointer=treatment_prep_root / "tcga_brca_treatment_os_overlap_v1_latest.json",
        os_endpoint_v1_latest_pointer=audit_root / "tcga-brca" / "endpoint-prep" / "tcga_brca_os_endpoint_v1_latest.json",
        clinical_biotabs_latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_biotabs_latest.json",
    )


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    for pointer_path, label in [
        (paths.treatment_os_overlap_latest_pointer, "treatment-OS overlap v1 latest pointer"),
        (paths.os_endpoint_v1_latest_pointer, "OS endpoint v1 latest pointer"),
        (paths.clinical_biotabs_latest_pointer, "clinical biotabs latest pointer"),
    ]:
        if not pointer_path.exists():
            raise PatientTreatmentProfileV1Error(f"Required {label} not found: {pointer_path}")

    overlap_pointer = helper_module.load_json(paths.treatment_os_overlap_latest_pointer)
    helper_module.require_keys(
        overlap_pointer,
        {
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "clinical_biotab_parse_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "treatment_os_overlap_summary_tsv",
            "run_log_json",
            "os_endpoint_v1_latest_json",
            "clinical_biotabs_latest_json",
        },
        "treatment-OS overlap v1 latest pointer",
        paths.treatment_os_overlap_latest_pointer,
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
    biotabs_pointer = helper_module.load_json(paths.clinical_biotabs_latest_pointer)
    helper_module.require_keys(
        biotabs_pointer,
        {"parse_run_id", "source_run_id", "processed_run_directory", "run_log_json"},
        "clinical biotabs latest pointer",
        paths.clinical_biotabs_latest_pointer,
    )

    overlap_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(overlap_pointer["run_log_json"]),
        label="treatment-OS overlap v1 run log",
        helper_module=helper_module,
    )
    os_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(os_pointer["run_log_json"]),
        label="OS endpoint v1 run log",
        helper_module=helper_module,
    )
    biotabs_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(biotabs_pointer["run_log_json"]),
        label="clinical biotabs run log",
        helper_module=helper_module,
    )

    if str(overlap_pointer["os_endpoint_v1_run_id"]) != str(os_pointer["os_endpoint_v1_run_id"]):
        raise PatientTreatmentProfileV1Error(
            "Mismatch between treatment-OS overlap and OS endpoint pointers for os_endpoint_v1_run_id."
        )
    if str(overlap_pointer["clinical_biotab_parse_run_id"]) != str(biotabs_pointer["parse_run_id"]):
        raise PatientTreatmentProfileV1Error(
            "Mismatch between treatment-OS overlap and clinical biotab pointers for parse_run_id."
        )
    if str(overlap_pointer["baseline_model_input_v1_run_id"]) != str(os_pointer["baseline_model_input_v1_run_id"]):
        raise PatientTreatmentProfileV1Error(
            "Mismatch between treatment-OS overlap and OS endpoint pointers for baseline_model_input_v1_run_id."
        )
    if str(overlap_pointer["cohort_v1_build_id"]) != str(os_pointer["cohort_v1_build_id"]):
        raise PatientTreatmentProfileV1Error(
            "Mismatch between treatment-OS overlap and OS endpoint pointers for cohort_v1_build_id."
        )

    os_endpoint_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(os_pointer["os_endpoint_v1_tsv"]),
        "os_endpoint_v1.tsv",
    )
    overlap_summary_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(overlap_pointer["treatment_os_overlap_summary_tsv"]),
        "treatment_os_overlap_summary.tsv",
    )
    biotabs_processed_dir = Path(str(biotabs_pointer["processed_run_directory"]))
    clinical_drug_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (biotabs_processed_dir / "clinical_drug.tsv").as_posix(),
        "clinical_drug.tsv",
    )
    clinical_radiation_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (biotabs_processed_dir / "clinical_radiation.tsv").as_posix(),
        "clinical_radiation.tsv",
    )
    clinical_patient_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (biotabs_processed_dir / "clinical_patient.tsv").as_posix(),
        "clinical_patient.tsv",
    )

    os_endpoint_rows = helper_module.read_tsv_dict_rows(os_endpoint_tsv)
    overlap_rows = helper_module.read_tsv_dict_rows(overlap_summary_tsv)
    clinical_drug_rows = helper_module.read_tsv_dict_rows(clinical_drug_tsv)
    clinical_radiation_rows = helper_module.read_tsv_dict_rows(clinical_radiation_tsv)
    clinical_patient_rows = helper_module.read_tsv_dict_rows(clinical_patient_tsv)

    require_columns(
        os_endpoint_rows,
        {
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
        },
        "os_endpoint_v1.tsv",
    )
    require_columns(
        overlap_rows,
        {
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "has_drug_row",
            "drug_row_count",
            "has_radiation_row",
            "radiation_row_count",
            "drug_therapy_type_values_json",
            "drug_therapy_type_distinct_count",
            "drug_therapy_type_single_or_mixed",
            "dominant_therapy_type_if_any",
        },
        "treatment_os_overlap_summary.tsv",
    )
    require_columns(
        clinical_drug_rows,
        {
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "pharmaceutical_therapy_type",
            "pharmaceutical_therapy_drug_name",
            "pharmaceutical_tx_started_days_to",
            "pharmaceutical_tx_ended_days_to",
            "therapy_regimen",
            "pharm_regimen",
        },
        "clinical_drug.tsv",
    )
    require_columns(
        clinical_radiation_rows,
        {
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "radiation_therapy_started_days_to",
            "radiation_therapy_ended_days_to",
            "therapy_regimen",
        },
        "clinical_radiation.tsv",
    )
    require_columns(
        clinical_patient_rows,
        {"bcr_patient_barcode", "bcr_patient_uuid"},
        "clinical_patient.tsv",
    )

    if not os_endpoint_rows:
        raise PatientTreatmentProfileV1Error("os_endpoint_v1.tsv has no data rows.")
    if not overlap_rows:
        raise PatientTreatmentProfileV1Error("treatment_os_overlap_summary.tsv has no data rows.")
    if not clinical_drug_rows:
        raise PatientTreatmentProfileV1Error("clinical_drug.tsv has no data rows.")
    if not clinical_radiation_rows:
        raise PatientTreatmentProfileV1Error("clinical_radiation.tsv has no data rows.")
    if not clinical_patient_rows:
        raise PatientTreatmentProfileV1Error("clinical_patient.tsv has no data rows.")

    return WorkflowInputs(
        treatment_os_overlap_latest_pointer=overlap_pointer,
        treatment_os_overlap_run_log=overlap_run_log,
        os_endpoint_v1_latest_pointer=os_pointer,
        os_endpoint_v1_run_log=os_run_log,
        clinical_biotabs_latest_pointer=biotabs_pointer,
        clinical_biotabs_run_log=biotabs_run_log,
        os_endpoint_rows=os_endpoint_rows,
        treatment_os_overlap_rows=overlap_rows,
        clinical_drug_rows=clinical_drug_rows,
        clinical_radiation_rows=clinical_radiation_rows,
        clinical_patient_rows=clinical_patient_rows,
        input_paths={
            "os_endpoint_v1_tsv": os_endpoint_tsv,
            "treatment_os_overlap_summary_tsv": overlap_summary_tsv,
            "clinical_drug_tsv": clinical_drug_tsv,
            "clinical_radiation_tsv": clinical_radiation_tsv,
            "clinical_patient_tsv": clinical_patient_tsv,
        },
    )


def build_profile_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
) -> tuple[list[dict[str, str]], ProfileMetrics]:
    drug_by_barcode = group_rows_by_barcode(workflow_inputs.clinical_drug_rows, barcode_field="bcr_patient_barcode")
    radiation_by_barcode = group_rows_by_barcode(
        workflow_inputs.clinical_radiation_rows,
        barcode_field="bcr_patient_barcode",
    )
    clinical_patient_by_barcode = build_unique_lookup_by_barcode(
        workflow_inputs.clinical_patient_rows,
        barcode_field="bcr_patient_barcode",
        label="clinical_patient.tsv",
    )
    overlap_by_barcode = build_unique_lookup_by_barcode(
        workflow_inputs.treatment_os_overlap_rows,
        barcode_field="bcr_patient_barcode",
        label="treatment_os_overlap_summary.tsv",
    )

    profile_rows: list[dict[str, str]] = []
    for os_row in workflow_inputs.os_endpoint_rows:
        barcode = normalize_barcode(str(os_row["bcr_patient_barcode"]))
        os_uuid = normalize_uuid(str(os_row["bcr_patient_uuid"]))
        overlap_row = overlap_by_barcode.get(barcode)
        clinical_patient_row = clinical_patient_by_barcode.get(barcode)

        if overlap_row is None:
            raise PatientTreatmentProfileV1Error(
                f"Missing treatment_os_overlap_summary.tsv row for barcode {barcode}."
            )
        if clinical_patient_row is None:
            raise PatientTreatmentProfileV1Error(f"Missing clinical_patient.tsv row for barcode {barcode}.")
        if normalize_uuid(str(overlap_row["bcr_patient_uuid"])) != os_uuid:
            raise PatientTreatmentProfileV1Error(
                f"UUID mismatch between os_endpoint_v1.tsv and treatment_os_overlap_summary.tsv for barcode {barcode}."
            )
        if normalize_uuid(str(clinical_patient_row["bcr_patient_uuid"])) != os_uuid:
            raise PatientTreatmentProfileV1Error(
                f"UUID mismatch between os_endpoint_v1.tsv and clinical_patient.tsv for barcode {barcode}."
            )
        if str(overlap_row["provisional_patient_row_id"]).strip() != str(os_row["provisional_patient_row_id"]).strip():
            raise PatientTreatmentProfileV1Error(
                "provisional_patient_row_id mismatch between os_endpoint_v1.tsv and "
                f"treatment_os_overlap_summary.tsv for barcode {barcode}."
            )

        drug_rows = drug_by_barcode.get(barcode, [])
        radiation_rows = radiation_by_barcode.get(barcode, [])

        drug_therapy_raw_values = collect_raw_values(drug_rows, "pharmaceutical_therapy_type")
        non_missing_therapy_values = ordered_non_missing_distinct_values(drug_therapy_raw_values)
        therapy_type_counts = ordered_count_objects(drug_therapy_raw_values)
        drug_name_raw_values = collect_raw_values(drug_rows, "pharmaceutical_therapy_drug_name")
        distinct_drug_names = ordered_non_missing_distinct_values(drug_name_raw_values)

        if not drug_rows:
            single_or_mixed = SINGLE_OR_MIXED_NO_DRUG
            dominant_therapy_type = ""
            dominant_therapy_type_rule = "no_drug_record"
            treatment_profile_status = STATUS_NO_DRUG_RECORD
        elif not non_missing_therapy_values:
            single_or_mixed = SINGLE_OR_MIXED_MISSING
            dominant_therapy_type = ""
            dominant_therapy_type_rule = "blank_due_to_no_usable_type"
            treatment_profile_status = STATUS_MISSING_TYPE_ONLY
        elif len(non_missing_therapy_values) == 1:
            single_or_mixed = SINGLE_OR_MIXED_SINGLE
            dominant_therapy_type = non_missing_therapy_values[0]
            dominant_therapy_type_rule = "single_non_missing_type"
            treatment_profile_status = STATUS_SINGLE
        else:
            single_or_mixed = SINGLE_OR_MIXED_MIXED
            sorted_counts = sorted(
                (
                    {"value": str(item["value"]), "count": int(item["count"])}
                    for item in therapy_type_counts
                ),
                key=lambda item: (-item["count"], item["value"]),
            )
            if len(sorted_counts) == 1 or sorted_counts[0]["count"] > sorted_counts[1]["count"]:
                dominant_therapy_type = sorted_counts[0]["value"]
                dominant_therapy_type_rule = "most_frequent_non_missing_type"
            else:
                dominant_therapy_type = ""
                dominant_therapy_type_rule = "blank_due_to_tie"
            treatment_profile_status = STATUS_MIXED

        regimen_context_values = collect_regimen_context_values(drug_rows, radiation_rows)
        if not regimen_context_values:
            regimen_context_single_or_mixed = "no_regimen_context"
        elif len(regimen_context_values) == 1:
            regimen_context_single_or_mixed = "single_context"
        else:
            regimen_context_single_or_mixed = "mixed_context"

        drug_start_values = collect_parseable_ints(drug_rows, "pharmaceutical_tx_started_days_to")
        drug_end_values = collect_parseable_ints(drug_rows, "pharmaceutical_tx_ended_days_to")
        radiation_start_values = collect_parseable_ints(radiation_rows, "radiation_therapy_started_days_to")
        radiation_end_values = collect_parseable_ints(radiation_rows, "radiation_therapy_ended_days_to")

        earliest_drug_start = min(drug_start_values) if drug_start_values else None
        latest_drug_end = max(drug_end_values) if drug_end_values else None
        earliest_radiation_start = min(radiation_start_values) if radiation_start_values else None
        latest_radiation_end = max(radiation_end_values) if radiation_end_values else None

        flags: list[str] = []
        if radiation_rows:
            flags.append("has_any_radiation_row")
        if drug_rows and radiation_rows:
            flags.append("has_both_drug_and_radiation_rows")
        if not drug_rows and radiation_rows:
            flags.append("radiation_without_drug_row")
        if len(regimen_context_values) > 1:
            flags.append("multiple_regimen_context_values")
        if dominant_therapy_type_rule == "blank_due_to_tie":
            flags.append("dominant_therapy_type_tie")
        if any("|" in value for value in non_missing_therapy_values):
            flags.append("compound_raw_therapy_type_present")
        if earliest_drug_start is not None or latest_drug_end is not None:
            flags.append("drug_timing_window_present")
        if earliest_radiation_start is not None or latest_radiation_end is not None:
            flags.append("radiation_timing_window_present")
        if earliest_drug_start is not None and latest_drug_end is not None and earliest_drug_start > latest_drug_end:
            flags.append("drug_timing_window_aggregated_inverted")
        if not drug_rows and not radiation_rows:
            flags.append("no_drug_or_radiation_rows")
        if single_or_mixed == SINGLE_OR_MIXED_MIXED:
            flags.append("mixed_therapy_types")
        if single_or_mixed == SINGLE_OR_MIXED_MISSING:
            flags.append("drug_rows_missing_therapy_type_only")
        if single_or_mixed == SINGLE_OR_MIXED_MIXED and dominant_therapy_type:
            flags.append("provisional_dominant_therapy_type_identified")
        flags = ordered_unique(flags)

        requires_manual_review = single_or_mixed in {SINGLE_OR_MIXED_MIXED, SINGLE_OR_MIXED_MISSING}
        has_any_treatment_timing = any(
            value is not None
            for value in [
                earliest_drug_start,
                latest_drug_end,
                earliest_radiation_start,
                latest_radiation_end,
            ]
        )

        profile_rows.append(
            {
                "patient_treatment_profile_v1_run_id": run_id,
                "treatment_os_overlap_v1_run_id": str(
                    workflow_inputs.treatment_os_overlap_latest_pointer["treatment_os_overlap_v1_run_id"]
                ),
                "os_endpoint_v1_run_id": str(
                    workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"]
                ),
                "clinical_biotab_parse_run_id": str(
                    workflow_inputs.clinical_biotabs_latest_pointer["parse_run_id"]
                ),
                "baseline_model_input_v1_run_id": str(
                    workflow_inputs.treatment_os_overlap_latest_pointer["baseline_model_input_v1_run_id"]
                ),
                "cohort_v1_build_id": str(
                    workflow_inputs.treatment_os_overlap_latest_pointer["cohort_v1_build_id"]
                ),
                "bcr_patient_barcode": str(os_row["bcr_patient_barcode"]),
                "bcr_patient_uuid": str(os_row["bcr_patient_uuid"]),
                "provisional_patient_row_id": str(os_row["provisional_patient_row_id"]),
                "baseline_analysis_v1_row_id": str(os_row["baseline_analysis_v1_row_id"]),
                "feature_set_v1_row_index": str(os_row["feature_set_v1_row_index"]),
                "has_any_drug_row": yes_no(bool(drug_rows)),
                "drug_row_count": str(len(drug_rows)),
                "has_any_radiation_row": yes_no(bool(radiation_rows)),
                "radiation_row_count": str(len(radiation_rows)),
                "drug_therapy_type_values_json": json_list(drug_therapy_raw_values),
                "drug_therapy_type_distinct_count": str(len(non_missing_therapy_values)),
                "drug_therapy_type_single_or_mixed": single_or_mixed,
                "dominant_therapy_type_if_any": dominant_therapy_type,
                "dominant_therapy_type_rule": dominant_therapy_type_rule,
                "therapy_type_counts_json": json.dumps(therapy_type_counts, ensure_ascii=True),
                "has_any_regimen_context": yes_no(bool(regimen_context_values)),
                "regimen_context_values_json": json_list(regimen_context_values),
                "regimen_context_distinct_count": str(len(regimen_context_values)),
                "regimen_context_single_or_mixed": regimen_context_single_or_mixed,
                "drug_name_values_json": json_list(drug_name_raw_values),
                "distinct_drug_name_count": str(len(distinct_drug_names)),
                "earliest_drug_start_days": "" if earliest_drug_start is None else str(earliest_drug_start),
                "latest_drug_end_days": "" if latest_drug_end is None else str(latest_drug_end),
                "earliest_radiation_start_days": (
                    "" if earliest_radiation_start is None else str(earliest_radiation_start)
                ),
                "latest_radiation_end_days": "" if latest_radiation_end is None else str(latest_radiation_end),
                "has_any_treatment_timing": yes_no(has_any_treatment_timing),
                "treatment_profile_requires_manual_review": yes_no(requires_manual_review),
                "treatment_profile_status": treatment_profile_status,
                "treatment_profile_flags_json": json_list(flags),
            }
        )

    metrics = build_profile_metrics(profile_rows)
    return profile_rows, metrics


def build_profile_metrics(profile_rows: list[dict[str, str]]) -> ProfileMetrics:
    return ProfileMetrics(
        total=len(profile_rows),
        patients_with_drug=sum(1 for row in profile_rows if row["has_any_drug_row"] == "yes"),
        patients_without_drug=sum(1 for row in profile_rows if row["has_any_drug_row"] == "no"),
        patients_with_radiation=sum(1 for row in profile_rows if row["has_any_radiation_row"] == "yes"),
        patients_without_radiation=sum(1 for row in profile_rows if row["has_any_radiation_row"] == "no"),
        patients_single_therapy_type=sum(
            1 for row in profile_rows if row["drug_therapy_type_single_or_mixed"] == SINGLE_OR_MIXED_SINGLE
        ),
        patients_mixed_therapy_type=sum(
            1 for row in profile_rows if row["drug_therapy_type_single_or_mixed"] == SINGLE_OR_MIXED_MIXED
        ),
        patients_missing_therapy_type_only=sum(
            1 for row in profile_rows if row["drug_therapy_type_single_or_mixed"] == SINGLE_OR_MIXED_MISSING
        ),
        patients_with_dominant_therapy_type_identifiable=sum(
            1 for row in profile_rows if str(row["dominant_therapy_type_if_any"]).strip()
        ),
        patients_requiring_manual_review=sum(
            1 for row in profile_rows if row["treatment_profile_requires_manual_review"] == "yes"
        ),
        patients_with_any_regimen_context=sum(
            1 for row in profile_rows if row["has_any_regimen_context"] == "yes"
        ),
        patients_with_any_treatment_timing=sum(
            1 for row in profile_rows if row["has_any_treatment_timing"] == "yes"
        ),
        patients_no_drug_or_radiation_rows=sum(
            1 for row in profile_rows if json_list_contains(row, "treatment_profile_flags_json", "no_drug_or_radiation_rows")
        ),
        patients_radiation_without_drug_row=sum(
            1 for row in profile_rows if json_list_contains(row, "treatment_profile_flags_json", "radiation_without_drug_row")
        ),
        patients_compound_raw_therapy_type_present=sum(
            1
            for row in profile_rows
            if json_list_contains(row, "treatment_profile_flags_json", "compound_raw_therapy_type_present")
        ),
        patients_drug_timing_window_aggregated_inverted=sum(
            1
            for row in profile_rows
            if json_list_contains(row, "treatment_profile_flags_json", "drug_timing_window_aggregated_inverted")
        ),
    )


def determine_conflict_type(profile_row: dict[str, str]) -> str:
    if profile_row["drug_therapy_type_single_or_mixed"] == SINGLE_OR_MIXED_MISSING:
        return "drug_rows_missing_therapy_type_only"
    if profile_row["drug_therapy_type_single_or_mixed"] == SINGLE_OR_MIXED_MIXED:
        if profile_row["dominant_therapy_type_rule"] == "blank_due_to_tie":
            return "mixed_therapy_types_tied_dominant"
        return "mixed_therapy_types_with_provisional_dominant"
    raise PatientTreatmentProfileV1Error(
        "Conflict audit requested for a row that does not require manual review: "
        f"{profile_row['bcr_patient_barcode']}"
    )


def build_source_evidence_summary(profile_row: dict[str, str]) -> str:
    parts = [
        f"drug_row_count={profile_row['drug_row_count']}",
        f"radiation_row_count={profile_row['radiation_row_count']}",
        f"drug_therapy_type_single_or_mixed={profile_row['drug_therapy_type_single_or_mixed']}",
        f"therapy_type_counts_json={profile_row['therapy_type_counts_json']}",
        f"dominant_therapy_type_if_any={profile_row['dominant_therapy_type_if_any'] or '[blank]'}",
        f"regimen_context_values_json={profile_row['regimen_context_values_json']}",
        f"distinct_drug_name_count={profile_row['distinct_drug_name_count']}",
        f"treatment_profile_flags_json={profile_row['treatment_profile_flags_json']}",
    ]
    return "; ".join(parts)


def recommended_handling_for_conflict(conflict_type: str) -> tuple[str, str]:
    if conflict_type == "mixed_therapy_types_tied_dominant":
        return (
            "Keep patient in the mixed/manual-review bucket and inspect raw drug rows before any coarse grouping assignment.",
            "high",
        )
    if conflict_type == "mixed_therapy_types_with_provisional_dominant":
        return (
            "Retain the provisional dominant therapy type for coarse exploratory grouping only, but keep the mixed/manual-review flag and do not freeze a final arm.",
            "medium",
        )
    if conflict_type == "drug_rows_missing_therapy_type_only":
        return (
            "Inspect raw drug rows and source forms, keep exposure unresolved, and do not infer untreated status or assign a therapy type.",
            "high",
        )
    raise PatientTreatmentProfileV1Error(f"Unhandled conflict_type: {conflict_type}")


def build_conflict_audit_rows(
    profile_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for profile_row in profile_rows:
        if profile_row["treatment_profile_requires_manual_review"] != "yes":
            continue
        conflict_type = determine_conflict_type(profile_row)
        recommended_handling, review_priority = recommended_handling_for_conflict(conflict_type)
        rows.append(
            {
                "patient_treatment_profile_v1_run_id": profile_row["patient_treatment_profile_v1_run_id"],
                "treatment_os_overlap_v1_run_id": profile_row["treatment_os_overlap_v1_run_id"],
                "os_endpoint_v1_run_id": profile_row["os_endpoint_v1_run_id"],
                "clinical_biotab_parse_run_id": profile_row["clinical_biotab_parse_run_id"],
                "baseline_model_input_v1_run_id": profile_row["baseline_model_input_v1_run_id"],
                "cohort_v1_build_id": profile_row["cohort_v1_build_id"],
                "bcr_patient_barcode": profile_row["bcr_patient_barcode"],
                "bcr_patient_uuid": profile_row["bcr_patient_uuid"],
                "provisional_patient_row_id": profile_row["provisional_patient_row_id"],
                "baseline_analysis_v1_row_id": profile_row["baseline_analysis_v1_row_id"],
                "feature_set_v1_row_index": profile_row["feature_set_v1_row_index"],
                "conflict_type": conflict_type,
                "source_evidence_summary": build_source_evidence_summary(profile_row),
                "recommended_next_handling": recommended_handling,
                "review_priority": review_priority,
                "drug_therapy_type_single_or_mixed": profile_row["drug_therapy_type_single_or_mixed"],
                "dominant_therapy_type_if_any": profile_row["dominant_therapy_type_if_any"],
                "dominant_therapy_type_rule": profile_row["dominant_therapy_type_rule"],
                "treatment_profile_status": profile_row["treatment_profile_status"],
                "treatment_profile_flags_json": profile_row["treatment_profile_flags_json"],
                "treatment_profile_requires_manual_review": profile_row["treatment_profile_requires_manual_review"],
            }
        )
    return rows


def build_spec_rows(run_id: str) -> list[dict[str, str]]:
    spec_meta: dict[str, dict[str, str]] = {
        "patient_treatment_profile_v1_run_id": {
            "source_origin": "workflow_generated",
            "field_category": "provenance",
            "decision_rule": "UTC run identifier assigned once per workflow execution.",
            "notes": "Written to every output row for audit linkage.",
        },
        "treatment_os_overlap_v1_run_id": {
            "source_origin": "treatment_os_overlap_v1_latest_pointer",
            "field_category": "provenance",
            "decision_rule": "Carried forward from the saved treatment-OS overlap v1 pointer.",
            "notes": "Preserves direct linkage to the upstream overlap audit.",
        },
        "os_endpoint_v1_run_id": {
            "source_origin": "os_endpoint_v1_latest_pointer",
            "field_category": "provenance",
            "decision_rule": "Carried forward from the saved OS endpoint v1 pointer.",
            "notes": "Links the treatment profile back to the frozen OS cohort source.",
        },
        "clinical_biotab_parse_run_id": {
            "source_origin": "clinical_biotabs_latest_pointer",
            "field_category": "provenance",
            "decision_rule": "Carried forward from the saved clinical biotab parse pointer.",
            "notes": "Identifies the parsed treatment-table run used here.",
        },
        "baseline_model_input_v1_run_id": {
            "source_origin": "treatment_os_overlap_v1_latest_pointer",
            "field_category": "provenance",
            "decision_rule": "Carried forward from the upstream overlap pointer without reloading the model-input layer.",
            "notes": "Included for cross-run provenance only.",
        },
        "cohort_v1_build_id": {
            "source_origin": "shared upstream pointers",
            "field_category": "provenance",
            "decision_rule": "Must reconcile across the saved overlap and OS endpoint pointers.",
            "notes": "Represents the OS cohort build carried into this patient-level profile.",
        },
        "bcr_patient_barcode": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row and used as the patient join key.",
            "notes": "One output row per OS cohort barcode.",
        },
        "bcr_patient_uuid": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row after cross-checking against overlap and clinical_patient linkage.",
            "notes": "Preserves UUID-level audit linkage.",
        },
        "provisional_patient_row_id": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row in existing cohort order.",
            "notes": "Used to preserve prior row-level audit linkage.",
        },
        "baseline_analysis_v1_row_id": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row without modification.",
            "notes": "Supports linkage to prior baseline analysis outputs.",
        },
        "feature_set_v1_row_index": {
            "source_origin": "os_endpoint_v1.tsv",
            "field_category": "identifier",
            "decision_rule": "Copied from the saved OS endpoint row without modification.",
            "notes": "Supports linkage to prior feature-set and model-input outputs.",
        },
        "has_any_drug_row": {
            "source_origin": "clinical_drug.tsv",
            "field_category": "coverage",
            "decision_rule": "yes when at least one clinical_drug row matches the patient barcode; otherwise no.",
            "notes": "No-drug-record remains explicit and is not relabeled as untreated.",
        },
        "drug_row_count": {
            "source_origin": "clinical_drug.tsv",
            "field_category": "coverage",
            "decision_rule": "Count all clinical_drug rows matching the patient barcode.",
            "notes": "Raw one-to-many treatment rows remain untouched upstream.",
        },
        "has_any_radiation_row": {
            "source_origin": "clinical_radiation.tsv",
            "field_category": "coverage",
            "decision_rule": "yes when at least one clinical_radiation row matches the patient barcode; otherwise no.",
            "notes": "Radiation is summarized alongside drug exposure but not merged into a final arm.",
        },
        "radiation_row_count": {
            "source_origin": "clinical_radiation.tsv",
            "field_category": "coverage",
            "decision_rule": "Count all clinical_radiation rows matching the patient barcode.",
            "notes": "Preserves patient-level radiation complexity without collapsing it into a treatment arm.",
        },
        "drug_therapy_type_values_json": {
            "source_origin": "clinical_drug.tsv.pharmaceutical_therapy_type",
            "field_category": "therapy_type",
            "decision_rule": "Store raw row-order therapy-type values, including duplicates and missing-like source values, as a JSON array.",
            "notes": "Raw labels stay visible; no normalization or splitting is applied.",
        },
        "drug_therapy_type_distinct_count": {
            "source_origin": "clinical_drug.tsv.pharmaceutical_therapy_type",
            "field_category": "therapy_type",
            "decision_rule": "Count distinct non-missing therapy-type labels after value.strip().",
            "notes": "Missing-like values do not contribute to the distinct count.",
        },
        "drug_therapy_type_single_or_mixed": {
            "source_origin": "derived_from_clinical_drug.tsv.pharmaceutical_therapy_type",
            "field_category": "therapy_type",
            "decision_rule": "Allowed values: no_drug_record, single_type, mixed_type, missing_type_only.",
            "notes": "Mixed and missing-type-only cases stay explicit instead of being forced into a final arm.",
        },
        "dominant_therapy_type_if_any": {
            "source_origin": "derived_from_clinical_drug.tsv.pharmaceutical_therapy_type",
            "field_category": "therapy_type",
            "decision_rule": "Use the most frequent non-missing exact raw therapy type when uniquely identifiable; otherwise leave blank.",
            "notes": "This is a provisional summary field only and is not a frozen treatment arm.",
        },
        "dominant_therapy_type_rule": {
            "source_origin": "workflow_generated",
            "field_category": "therapy_type",
            "decision_rule": "Allowed values: no_drug_record, single_non_missing_type, most_frequent_non_missing_type, blank_due_to_tie, blank_due_to_no_usable_type.",
            "notes": "Documents why dominant_therapy_type_if_any is populated or blank.",
        },
        "therapy_type_counts_json": {
            "source_origin": "clinical_drug.tsv.pharmaceutical_therapy_type",
            "field_category": "therapy_type",
            "decision_rule": "Store ordered count objects {value, count} for non-missing stripped therapy-type labels.",
            "notes": "Order follows first appearance in the raw drug rows.",
        },
        "has_any_regimen_context": {
            "source_origin": "combined drug and radiation regimen fields",
            "field_category": "regimen_context",
            "decision_rule": "yes when at least one non-missing regimen-context value is present across drug therapy_regimen, drug pharm_regimen, or radiation therapy_regimen.",
            "notes": "Uses the combined-minimal regimen-context scope chosen for v1.",
        },
        "regimen_context_values_json": {
            "source_origin": "clinical_drug.tsv.therapy_regimen|pharm_regimen and clinical_radiation.tsv.therapy_regimen",
            "field_category": "regimen_context",
            "decision_rule": "Store distinct non-missing raw regimen-context values in first-seen order as a JSON array.",
            "notes": "No prioritization is applied beyond preserving first-seen order.",
        },
        "regimen_context_distinct_count": {
            "source_origin": "combined regimen context fields",
            "field_category": "regimen_context",
            "decision_rule": "Count distinct non-missing regimen-context values after value.strip().",
            "notes": "The currently empty drug pharm_regimen field remains included for forward compatibility.",
        },
        "regimen_context_single_or_mixed": {
            "source_origin": "derived_from_combined_regimen_context",
            "field_category": "regimen_context",
            "decision_rule": "Allowed values: no_regimen_context, single_context, mixed_context.",
            "notes": "Context is summarized but not used to freeze a final treatment arm.",
        },
        "drug_name_values_json": {
            "source_origin": "clinical_drug.tsv.pharmaceutical_therapy_drug_name",
            "field_category": "drug_name",
            "decision_rule": "Store raw row-order drug-name values, including duplicates and missing-like source values, as a JSON array.",
            "notes": "Drug names remain completely unnormalized in v1.",
        },
        "distinct_drug_name_count": {
            "source_origin": "clinical_drug.tsv.pharmaceutical_therapy_drug_name",
            "field_category": "drug_name",
            "decision_rule": "Count distinct non-missing exact raw drug-name labels after value.strip().",
            "notes": "Used only as a structural complexity summary, not as an arm-definition rule.",
        },
        "earliest_drug_start_days": {
            "source_origin": "clinical_drug.tsv.pharmaceutical_tx_started_days_to",
            "field_category": "timing",
            "decision_rule": "Minimum parseable integer among the patient's drug start-day values.",
            "notes": "Left blank if no parseable drug start-day value exists.",
        },
        "latest_drug_end_days": {
            "source_origin": "clinical_drug.tsv.pharmaceutical_tx_ended_days_to",
            "field_category": "timing",
            "decision_rule": "Maximum parseable integer among the patient's drug end-day values.",
            "notes": "Left blank if no parseable drug end-day value exists.",
        },
        "earliest_radiation_start_days": {
            "source_origin": "clinical_radiation.tsv.radiation_therapy_started_days_to",
            "field_category": "timing",
            "decision_rule": "Minimum parseable integer among the patient's radiation start-day values.",
            "notes": "Left blank if no parseable radiation start-day value exists.",
        },
        "latest_radiation_end_days": {
            "source_origin": "clinical_radiation.tsv.radiation_therapy_ended_days_to",
            "field_category": "timing",
            "decision_rule": "Maximum parseable integer among the patient's radiation end-day values.",
            "notes": "Left blank if no parseable radiation end-day value exists.",
        },
        "has_any_treatment_timing": {
            "source_origin": "derived_from_aggregated_timing_fields",
            "field_category": "timing",
            "decision_rule": "yes when any aggregated drug or radiation timing field is populated; otherwise no.",
            "notes": "Timing coverage is summarized structurally only.",
        },
        "treatment_profile_requires_manual_review": {
            "source_origin": "workflow_generated",
            "field_category": "review",
            "decision_rule": "yes for mixed_type and missing_type_only patients; otherwise no.",
            "notes": "This is the narrow-core manual-review scope selected for v1.",
        },
        "treatment_profile_status": {
            "source_origin": "workflow_generated",
            "field_category": "review",
            "decision_rule": "Allowed values: no_drug_record, single_drug_therapy_type, mixed_drug_therapy_types, drug_rows_missing_therapy_type_only.",
            "notes": "Separates no-drug, single-type, mixed-type, and unusable-type cases without freezing a final arm.",
        },
        "treatment_profile_flags_json": {
            "source_origin": "workflow_generated",
            "field_category": "review",
            "decision_rule": "Store ordered unique structural flags as a JSON array.",
            "notes": "Includes raw-complexity markers such as radiation_without_drug_row, dominant_therapy_type_tie, and drug_timing_window_aggregated_inverted.",
        },
    }

    missing_meta_fields = [field_name for field_name in PROFILE_FIELDNAMES if field_name not in spec_meta]
    if missing_meta_fields:
        raise PatientTreatmentProfileV1Error(
            f"Spec metadata missing for profile fields: {missing_meta_fields}"
        )

    return [
        {
            "patient_treatment_profile_v1_run_id": run_id,
            "field_name": field_name,
            "source_origin": spec_meta[field_name]["source_origin"],
            "field_category": spec_meta[field_name]["field_category"],
            "decision_rule": spec_meta[field_name]["decision_rule"],
            "notes": spec_meta[field_name]["notes"],
        }
        for field_name in PROFILE_FIELDNAMES
    ]


def determine_readiness(metrics: ProfileMetrics) -> tuple[str, str]:
    if metrics.total == 0 or metrics.patients_with_drug == 0:
        return (
            READINESS_NOT_READY,
            "No usable patient-level drug coverage was generated, so coarse exploratory grouping would still be blocked.",
        )
    return (
        READINESS_COARSE_GROUPING,
        "Patient-level treatment structure is cohort-complete, preserves no-drug and mixed cases explicitly, and remains compatible with later arm-freeze review.",
    )


def build_summary_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
    metrics: ProfileMetrics,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    def add(summary_section: str, summary_metric: str, summary_value: Any, notes: str = "") -> None:
        rows.append(
            {
                "patient_treatment_profile_v1_run_id": run_id,
                "summary_section": summary_section,
                "summary_metric": summary_metric,
                "summary_value": str(summary_value),
                "notes": notes,
            }
        )

    readiness_interpretation, readiness_notes = determine_readiness(metrics)

    add("inputs", "treatment_os_overlap_v1_run_id", workflow_inputs.treatment_os_overlap_latest_pointer["treatment_os_overlap_v1_run_id"])
    add("inputs", "os_endpoint_v1_run_id", workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"])
    add("inputs", "clinical_biotab_parse_run_id", workflow_inputs.clinical_biotabs_latest_pointer["parse_run_id"])
    add("inputs", "baseline_model_input_v1_run_id", workflow_inputs.treatment_os_overlap_latest_pointer["baseline_model_input_v1_run_id"])
    add("inputs", "cohort_v1_build_id", workflow_inputs.treatment_os_overlap_latest_pointer["cohort_v1_build_id"])

    add("coverage", "total_os_cohort_size", metrics.total)
    add("coverage", "patients_with_any_drug_row", metrics.patients_with_drug)
    add(
        "coverage",
        "patients_with_no_drug_row",
        metrics.patients_without_drug,
        "No-drug-record remains explicit and is not interpreted as untreated.",
    )
    add("coverage", "patients_with_any_radiation_row", metrics.patients_with_radiation)
    add("coverage", "patients_with_no_radiation_row", metrics.patients_without_radiation)
    add("coverage", "patients_with_no_drug_or_radiation_rows", metrics.patients_no_drug_or_radiation_rows)
    add("coverage", "patients_with_radiation_without_drug_row", metrics.patients_radiation_without_drug_row)

    add("therapy_type", "patients_with_single_therapy_type", metrics.patients_single_therapy_type)
    add("therapy_type", "patients_with_mixed_therapy_types", metrics.patients_mixed_therapy_type)
    add("therapy_type", "patients_with_missing_therapy_type_only", metrics.patients_missing_therapy_type_only)
    add(
        "therapy_type",
        "patients_with_dominant_therapy_type_identifiable",
        metrics.patients_with_dominant_therapy_type_identifiable,
        "Dominant therapy type remains provisional and is not a frozen final arm.",
    )
    add(
        "therapy_type",
        "patients_with_compound_raw_therapy_type_present",
        metrics.patients_compound_raw_therapy_type_present,
        "Compound raw labels are preserved exactly and not split in v1.",
    )

    add("regimen_context", "patients_with_any_regimen_context", metrics.patients_with_any_regimen_context)
    add("timing", "patients_with_any_treatment_timing", metrics.patients_with_any_treatment_timing)
    add(
        "timing",
        "patients_with_drug_timing_window_aggregated_inverted",
        metrics.patients_drug_timing_window_aggregated_inverted,
        "These are flagged structurally for review but are not part of the narrow-core manual-review trigger in v1.",
    )

    add("manual_review", "patients_requiring_manual_review", metrics.patients_requiring_manual_review)
    add(
        "manual_review",
        "manual_review_scope",
        "narrow_core",
        "Manual review includes mixed therapy types, dominant-type ties, and drug rows with unusable therapy type.",
    )

    add("readiness", "provisional_readiness_interpretation", readiness_interpretation, readiness_notes)
    add(
        "readiness",
        "treatment_arm_freeze_review_status",
        "blocked",
        "Drug names remain unnormalized and manual-review cases remain unresolved.",
    )
    add(
        "readiness",
        "remaining_blocker",
        "drug_name_not_normalized_and_manual_review_remaining",
        "This workflow is intentionally patient-level treatment profiling only.",
    )
    add(
        "readiness",
        "treatment_recommendation_modeling_status",
        "out_of_scope",
        "Treatment recommendation modeling remains out of scope for this phase.",
    )
    add(
        "readiness",
        "treatment_profile_scope",
        "construction_only",
        "This workflow does not freeze treatment arms, normalize drug names, or perform causal analysis.",
    )

    return rows


def validate_outputs(
    workflow_inputs: WorkflowInputs,
    profile_rows: list[dict[str, str]],
    spec_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    conflict_rows: list[dict[str, str]],
    metrics: ProfileMetrics,
) -> dict[str, Any]:
    os_barcodes = [normalize_barcode(row["bcr_patient_barcode"]) for row in workflow_inputs.os_endpoint_rows]
    overlap_lookup = build_unique_lookup_by_barcode(
        workflow_inputs.treatment_os_overlap_rows,
        barcode_field="bcr_patient_barcode",
        label="treatment_os_overlap_summary.tsv",
    )
    profile_barcodes = [normalize_barcode(row["bcr_patient_barcode"]) for row in profile_rows]
    clinical_patient_lookup = build_unique_lookup_by_barcode(
        workflow_inputs.clinical_patient_rows,
        barcode_field="bcr_patient_barcode",
        label="clinical_patient.tsv",
    )

    row_order_preserved = profile_barcodes == os_barcodes
    all_os_patients_in_profile = set(profile_barcodes) == set(os_barcodes)
    no_profile_duplicate_barcodes = len(profile_barcodes) == len(set(profile_barcodes))
    clinical_patient_linkage_complete = all(barcode in clinical_patient_lookup for barcode in os_barcodes)

    overlap_continuity_by_patient = True
    for profile_row in profile_rows:
        overlap_row = overlap_lookup.get(normalize_barcode(profile_row["bcr_patient_barcode"]))
        if overlap_row is None:
            overlap_continuity_by_patient = False
            break
        if profile_row["has_any_drug_row"] != overlap_row["has_drug_row"]:
            overlap_continuity_by_patient = False
            break
        if profile_row["drug_row_count"] != overlap_row["drug_row_count"]:
            overlap_continuity_by_patient = False
            break
        if profile_row["has_any_radiation_row"] != overlap_row["has_radiation_row"]:
            overlap_continuity_by_patient = False
            break
        if profile_row["radiation_row_count"] != overlap_row["radiation_row_count"]:
            overlap_continuity_by_patient = False
            break
        if profile_row["drug_therapy_type_distinct_count"] != overlap_row["drug_therapy_type_distinct_count"]:
            overlap_continuity_by_patient = False
            break
        if profile_row["drug_therapy_type_single_or_mixed"] != single_or_mixed_from_overlap(overlap_row):
            overlap_continuity_by_patient = False
            break
        if profile_row["dominant_therapy_type_if_any"] != overlap_row["dominant_therapy_type_if_any"]:
            overlap_continuity_by_patient = False
            break

    overlap_category_counts_profile = Counter(
        row["drug_therapy_type_single_or_mixed"] for row in profile_rows
    )
    overlap_category_counts_source = Counter(
        single_or_mixed_from_overlap(row) for row in workflow_inputs.treatment_os_overlap_rows
    )
    overlap_therapy_type_summary_reconciles = (
        overlap_category_counts_profile == overlap_category_counts_source
        and metrics.patients_with_dominant_therapy_type_identifiable
        == sum(1 for row in workflow_inputs.treatment_os_overlap_rows if str(row["dominant_therapy_type_if_any"]).strip())
    )

    summary_metric_lookup = {
        row["summary_metric"]: row["summary_value"]
        for row in summary_rows
    }

    def summary_int(metric_name: str) -> int | None:
        raw_value = summary_metric_lookup.get(metric_name, "")
        return parse_int_or_none(raw_value)

    summary_counts_reconcile = (
        summary_int("total_os_cohort_size") == metrics.total
        and summary_int("patients_with_any_drug_row") == metrics.patients_with_drug
        and summary_int("patients_with_no_drug_row") == metrics.patients_without_drug
        and summary_int("patients_with_any_radiation_row") == metrics.patients_with_radiation
        and summary_int("patients_with_no_radiation_row") == metrics.patients_without_radiation
        and summary_int("patients_with_single_therapy_type") == metrics.patients_single_therapy_type
        and summary_int("patients_with_mixed_therapy_types") == metrics.patients_mixed_therapy_type
        and summary_int("patients_with_missing_therapy_type_only") == metrics.patients_missing_therapy_type_only
        and summary_int("patients_with_dominant_therapy_type_identifiable")
        == metrics.patients_with_dominant_therapy_type_identifiable
        and summary_int("patients_requiring_manual_review") == metrics.patients_requiring_manual_review
        and summary_int("patients_with_any_regimen_context") == metrics.patients_with_any_regimen_context
        and summary_int("patients_with_any_treatment_timing") == metrics.patients_with_any_treatment_timing
    )

    readiness_value = summary_metric_lookup.get("provisional_readiness_interpretation", "")
    readiness_allowed = readiness_value in {
        READINESS_NOT_READY,
        READINESS_COARSE_GROUPING,
        READINESS_ARM_FREEZE_REVIEW,
    }

    patient_counts_reconcile = (
        metrics.patients_with_drug + metrics.patients_without_drug == metrics.total
        and metrics.patients_with_radiation + metrics.patients_without_radiation == metrics.total
        and metrics.patients_single_therapy_type
        + metrics.patients_mixed_therapy_type
        + metrics.patients_missing_therapy_type_only
        + metrics.patients_without_drug
        == metrics.total
    )

    spec_field_names_match = [row["field_name"] for row in spec_rows] == PROFILE_FIELDNAMES
    conflict_row_count_matches_manual_review = len(conflict_rows) == metrics.patients_requiring_manual_review

    passed = all(
        [
            len(workflow_inputs.os_endpoint_rows) > 0,
            len(workflow_inputs.treatment_os_overlap_rows) > 0,
            len(workflow_inputs.clinical_drug_rows) > 0,
            len(workflow_inputs.clinical_radiation_rows) > 0,
            len(profile_rows) > 0,
            len(profile_rows) == len(workflow_inputs.os_endpoint_rows),
            row_order_preserved,
            all_os_patients_in_profile,
            no_profile_duplicate_barcodes,
            clinical_patient_linkage_complete,
            patient_counts_reconcile,
            summary_counts_reconcile,
            conflict_row_count_matches_manual_review,
            len(spec_rows) == len(PROFILE_FIELDNAMES),
            spec_field_names_match,
            overlap_continuity_by_patient,
            overlap_therapy_type_summary_reconciles,
            readiness_allowed,
        ]
    )

    return {
        "passed": passed,
        "required_upstream_pointers_found": True,
        "required_source_tables_found": True,
        "treatment_os_overlap_run_log_completed": True,
        "os_endpoint_run_log_completed": True,
        "clinical_biotabs_run_log_completed": True,
        "os_cohort_row_count_positive": len(workflow_inputs.os_endpoint_rows) > 0,
        "overlap_summary_row_count_positive": len(workflow_inputs.treatment_os_overlap_rows) > 0,
        "clinical_drug_row_count_positive": len(workflow_inputs.clinical_drug_rows) > 0,
        "clinical_radiation_row_count_positive": len(workflow_inputs.clinical_radiation_rows) > 0,
        "clinical_patient_row_count_positive": len(workflow_inputs.clinical_patient_rows) > 0,
        "output_rows_positive": len(profile_rows) > 0,
        "profile_row_count_matches_os_cohort": len(profile_rows) == len(workflow_inputs.os_endpoint_rows),
        "row_order_preserved_from_os_endpoint": row_order_preserved,
        "all_os_patients_in_profile": all_os_patients_in_profile,
        "no_profile_duplicate_barcodes": no_profile_duplicate_barcodes,
        "clinical_patient_linkage_complete": clinical_patient_linkage_complete,
        "patient_counts_reconcile_to_os_cohort": patient_counts_reconcile,
        "summary_counts_reconcile_to_patient_table": summary_counts_reconcile,
        "conflict_audit_row_count_matches_manual_review": conflict_row_count_matches_manual_review,
        "spec_row_count_matches_profile_columns": len(spec_rows) == len(PROFILE_FIELDNAMES),
        "spec_field_names_match_profile_columns": spec_field_names_match,
        "overlap_row_count_matches_os_cohort": (
            len(workflow_inputs.treatment_os_overlap_rows) == len(workflow_inputs.os_endpoint_rows)
        ),
        "overlap_continuity_profile_vs_overlap": overlap_continuity_by_patient,
        "overlap_therapy_type_summary_reconciles": overlap_therapy_type_summary_reconciles,
        "summary_readiness_interpretation_allowed": readiness_allowed,
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
        "os_cohort_row_count": len(workflow_inputs.os_endpoint_rows),
        "profile_row_count": len(profile_rows),
        "manual_review_row_count": len(conflict_rows),
        "spec_row_count": len(spec_rows),
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
        "patient_treatment_profile_v1_run_id": run_id,
        "treatment_os_overlap_v1_run_id": str(
            workflow_inputs.treatment_os_overlap_latest_pointer["treatment_os_overlap_v1_run_id"]
        ),
        "os_endpoint_v1_run_id": str(
            workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"]
        ),
        "clinical_biotab_parse_run_id": str(
            workflow_inputs.clinical_biotabs_latest_pointer["parse_run_id"]
        ),
        "baseline_model_input_v1_run_id": str(
            workflow_inputs.treatment_os_overlap_latest_pointer["baseline_model_input_v1_run_id"]
        ),
        "cohort_v1_build_id": str(
            workflow_inputs.treatment_os_overlap_latest_pointer["cohort_v1_build_id"]
        ),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "patient_treatment_profile_v1_tsv": repo_relative(
            output_paths["patient_treatment_profile_v1_tsv"],
            paths.repo_root,
        ),
        "patient_treatment_profile_v1_spec_tsv": repo_relative(
            output_paths["patient_treatment_profile_v1_spec_tsv"],
            paths.repo_root,
        ),
        "patient_treatment_profile_v1_summary_tsv": repo_relative(
            output_paths["patient_treatment_profile_v1_summary_tsv"],
            paths.repo_root,
        ),
        "patient_treatment_profile_v1_conflict_audit_tsv": repo_relative(
            output_paths["patient_treatment_profile_v1_conflict_audit_tsv"],
            paths.repo_root,
        ),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "treatment_os_overlap_v1_latest_json": repo_relative(
            paths.treatment_os_overlap_latest_pointer,
            paths.repo_root,
        ),
        "os_endpoint_v1_latest_json": repo_relative(
            paths.os_endpoint_v1_latest_pointer,
            paths.repo_root,
        ),
        "clinical_biotabs_latest_json": repo_relative(
            paths.clinical_biotabs_latest_pointer,
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

        patient_treatment_profile_path = processed_run_dir / "patient_treatment_profile_v1.tsv"
        spec_path = audit_run_dir / "patient_treatment_profile_v1_spec.tsv"
        summary_path = audit_run_dir / "patient_treatment_profile_v1_summary.tsv"
        conflict_audit_path = audit_run_dir / "patient_treatment_profile_v1_conflict_audit.tsv"

        profile_rows, metrics = build_profile_rows(run_id, workflow_inputs)
        spec_rows = build_spec_rows(run_id)
        summary_rows = build_summary_rows(run_id, workflow_inputs, metrics)
        conflict_rows = build_conflict_audit_rows(profile_rows)

        validation = validate_outputs(
            workflow_inputs=workflow_inputs,
            profile_rows=profile_rows,
            spec_rows=spec_rows,
            summary_rows=summary_rows,
            conflict_rows=conflict_rows,
            metrics=metrics,
        )
        if not validation["passed"]:
            raise PatientTreatmentProfileV1Error(
                "Output validation did not pass. Failed checks: "
                + str({key: value for key, value in validation.items() if value is False})
            )

        helper_module.write_dict_rows_tsv(
            patient_treatment_profile_path,
            PROFILE_FIELDNAMES,
            profile_rows,
        )
        helper_module.write_dict_rows_tsv(spec_path, SPEC_FIELDNAMES, spec_rows)
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)
        helper_module.write_dict_rows_tsv(conflict_audit_path, CONFLICT_AUDIT_FIELDNAMES, conflict_rows)

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "patient_treatment_profile_v1_tsv": patient_treatment_profile_path,
            "patient_treatment_profile_v1_spec_tsv": spec_path,
            "patient_treatment_profile_v1_summary_tsv": summary_path,
            "patient_treatment_profile_v1_conflict_audit_tsv": conflict_audit_path,
            "run_log_json": run_log_path,
        }

        latest_pointer_payload = build_latest_pointer_payload(
            run_id=run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        readiness_value = next(
            row["summary_value"]
            for row in summary_rows
            if row["summary_metric"] == "provisional_readiness_interpretation"
        )

        completed_at = utc_now()
        run_log_payload: dict[str, Any] = {
            "status": "completed",
            "patient_treatment_profile_v1_run_id": run_id,
            "treatment_os_overlap_v1_run_id": str(
                workflow_inputs.treatment_os_overlap_latest_pointer["treatment_os_overlap_v1_run_id"]
            ),
            "os_endpoint_v1_run_id": str(
                workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"]
            ),
            "clinical_biotab_parse_run_id": str(
                workflow_inputs.clinical_biotabs_latest_pointer["parse_run_id"]
            ),
            "baseline_model_input_v1_run_id": str(
                workflow_inputs.treatment_os_overlap_latest_pointer["baseline_model_input_v1_run_id"]
            ),
            "cohort_v1_build_id": str(
                workflow_inputs.treatment_os_overlap_latest_pointer["cohort_v1_build_id"]
            ),
            "clinical_source_run_id": str(
                workflow_inputs.clinical_biotabs_latest_pointer["source_run_id"]
            ),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "treatment_os_overlap_v1_latest_json": repo_relative(
                    paths.treatment_os_overlap_latest_pointer,
                    paths.repo_root,
                ),
                "os_endpoint_v1_latest_json": repo_relative(
                    paths.os_endpoint_v1_latest_pointer,
                    paths.repo_root,
                ),
                "clinical_biotabs_latest_json": repo_relative(
                    paths.clinical_biotabs_latest_pointer,
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
                "patient_treatment_profile_v1_tsv": repo_relative(
                    patient_treatment_profile_path,
                    paths.repo_root,
                ),
                "patient_treatment_profile_v1_spec_tsv": repo_relative(spec_path, paths.repo_root),
                "patient_treatment_profile_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "patient_treatment_profile_v1_conflict_audit_tsv": repo_relative(
                    conflict_audit_path,
                    paths.repo_root,
                ),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation,
            "rules": {
                "unit_of_analysis": "patient/case",
                "input_layer": "saved_os_endpoint_and_treatment_overlap_plus_clinical_treatment_biotabs",
                "output_layer": "patient_treatment_profile_v1",
                "no_drug_name_normalization": True,
                "no_treatment_arm_freeze": True,
                "no_modeling": True,
                "no_causal_analysis": True,
                "treatment_profile_only": True,
                "regimen_context_scope": "combined_minimal",
                "manual_review_scope": "narrow_core",
                "raw_compound_therapy_type_preserved": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "missing_like_tokens": sorted(MISSING_LIKE_TOKENS),
            },
            "counts": {
                "os_endpoint_row_count": len(workflow_inputs.os_endpoint_rows),
                "treatment_os_overlap_row_count": len(workflow_inputs.treatment_os_overlap_rows),
                "clinical_drug_row_count": len(workflow_inputs.clinical_drug_rows),
                "clinical_radiation_row_count": len(workflow_inputs.clinical_radiation_rows),
                "clinical_patient_row_count": len(workflow_inputs.clinical_patient_rows),
                "patient_treatment_profile_row_count": len(profile_rows),
                "patient_treatment_profile_spec_row_count": len(spec_rows),
                "patient_treatment_profile_summary_row_count": len(summary_rows),
                "patient_treatment_profile_conflict_audit_row_count": len(conflict_rows),
                "patients_with_drug": metrics.patients_with_drug,
                "patients_without_drug": metrics.patients_without_drug,
                "patients_with_radiation": metrics.patients_with_radiation,
                "patients_without_radiation": metrics.patients_without_radiation,
                "patients_single_therapy_type": metrics.patients_single_therapy_type,
                "patients_mixed_therapy_type": metrics.patients_mixed_therapy_type,
                "patients_missing_therapy_type_only": metrics.patients_missing_therapy_type_only,
                "patients_with_dominant_therapy_type_identifiable": (
                    metrics.patients_with_dominant_therapy_type_identifiable
                ),
                "patients_requiring_manual_review": metrics.patients_requiring_manual_review,
                "patients_with_any_regimen_context": metrics.patients_with_any_regimen_context,
                "patients_with_any_treatment_timing": metrics.patients_with_any_treatment_timing,
                "patients_no_drug_or_radiation_rows": metrics.patients_no_drug_or_radiation_rows,
                "patients_compound_raw_therapy_type_present": metrics.patients_compound_raw_therapy_type_present,
                "patients_drug_timing_window_aggregated_inverted": (
                    metrics.patients_drug_timing_window_aggregated_inverted
                ),
                "provisional_readiness_interpretation": readiness_value,
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "treatment_os_overlap_v1_latest_pointer": workflow_inputs.treatment_os_overlap_latest_pointer,
                "os_endpoint_v1_latest_pointer": workflow_inputs.os_endpoint_v1_latest_pointer,
                "clinical_biotabs_latest_pointer": workflow_inputs.clinical_biotabs_latest_pointer,
            },
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "patient_treatment_profile_v1_run_id": run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_patient_treatment_profile_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    counts = run_log.get("counts", {})
    print("TCGA-BRCA patient treatment profile v1 workflow complete.")
    print(f"  Run ID                : {run_log['patient_treatment_profile_v1_run_id']}")
    print(f"  OS cohort size        : {counts.get('os_endpoint_row_count', '?')}")
    print(f"  Patients with drug    : {counts.get('patients_with_drug', '?')}")
    print(f"  Patients with rad     : {counts.get('patients_with_radiation', '?')}")
    print(f"  Single therapy type   : {counts.get('patients_single_therapy_type', '?')}")
    print(f"  Mixed therapy types   : {counts.get('patients_mixed_therapy_type', '?')}")
    print(f"  Manual review         : {counts.get('patients_requiring_manual_review', '?')}")
    print(
        "  Readiness             : "
        f"{counts.get('provisional_readiness_interpretation', '?')}"
    )
    print(f"  Processed run dir     : {run_log['outputs']['processed_run_directory']}")
    print(f"  Audit run dir         : {run_log['outputs']['audit_run_directory']}")
    print(f"  Latest pointer        : {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
