#!/usr/bin/env python
"""Review TCGA-BRCA drug-class profiles and produce arm-freeze candidate audit v1.

This is a descriptive review and arm-freeze candidate audit only.
This script does NOT normalize drug names, redefine treatment groups, freeze treatment arms,
perform causal analysis, perform treatment-effect estimation, or produce treatment recommendations.
"""

from __future__ import annotations

import importlib.util
import json
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

APPROVED_BUCKET_ORDER = [
    "single_known_class_clean",
    "single_known_class_plus_unknown",
    "multi_known_class",
    "unknown_only",
    "missing_like_only",
]
PROVISIONAL_CLASS_VOCABULARY = [
    "alkylating_agent",
    "anthracycline",
    "taxane",
    "platinum",
    "antimetabolite",
    "endocrine_serm",
    "endocrine_aromatase_inhibitor",
    "endocrine_other",
    "her2_targeted",
    "immunotherapy",
    "ancillary_supportive",
    "other_cytotoxic",
    "unknown_or_review_needed",
]
KNOWN_CLASS_VOCABULARY = [c for c in PROVISIONAL_CLASS_VOCABULARY if c != "unknown_or_review_needed"]
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
MIN_CLASS_SIZE_FOR_AUDIT = 20
MIN_CLASS_SIZE_FOR_CANDIDATE = 50
MANUAL_REVIEW_FRACTION_TOLERANCE = 0.20
CONFOUNDER_NON_MISSING_TOLERANCE = 0.85
MAX_TOP_VALUES = 10
READINESS_NOT_READY = "not_ready"
READINESS_DESCRIPTIVE_COMPLETED = "ready_for_drug_class_descriptive_review_completed"
READINESS_ARM_FREEZE_CANDIDATE_DISCUSSION = "ready_for_treatment_arm_freeze_candidate_discussion"
TREATMENT_ARM_FREEZE_REVIEW_BLOCKED = "blocked"
TREATMENT_ARM_FREEZE_REVIEW_PROVISIONAL_NEXT = "provisional_review_next_only"
TREATMENT_RECOMMENDATION_MODELING_STATUS = "out_of_scope"

GROUP_SUMMARY_FIELDNAMES = [
    "drug_class_profile_v1_run_id",
    "bucket_order",
    "drug_class_profile_bucket",
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
    "drug_class_profile_v1_run_id",
    "review_group_type",
    "field_name",
    "row_count",
    "non_missing_count",
    "missing_like_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "notes",
]
VALUE_COMPOSITION_FIELDNAMES = [
    "drug_class_profile_v1_run_id",
    "review_group_type",
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
ARM_FREEZE_CANDIDATE_AUDIT_FIELDNAMES = [
    "drug_class_profile_v1_run_id",
    "candidate_label",
    "candidate_type",
    "patient_count",
    "os_event_count",
    "manual_review_count",
    "unknown_mapping_burden",
    "multi_class_burden",
    "confounder_coverage_ok",
    "arm_freeze_candidate_status",
    "why_candidate_or_not",
    "recommended_next_handling",
    "notes",
]
SUMMARY_FIELDNAMES = [
    "drug_class_profile_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]


class DrugClassProfileV1Error(RuntimeError):
    """Raised when the drug-class profile review workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the drug-class profile review workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    drug_norm_latest_pointer: Path
    grouping_latest_pointer: Path
    baseline_analysis_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    drug_norm_pointer: dict[str, Any]
    drug_norm_run_log: dict[str, Any]
    grouping_pointer: dict[str, Any]
    grouping_run_log: dict[str, Any]
    baseline_analysis_pointer: dict[str, Any]
    baseline_analysis_run_log: dict[str, Any]
    drug_class_profile_rows: list[dict[str, str]]
    grouping_rows: list[dict[str, str]]
    baseline_analysis_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise DrugClassProfileV1Error(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise DrugClassProfileV1Error(f"Unable to create an import spec for: {script_path}")
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


def normalize_missing_like(value: Any) -> str:
    return str(value or "").strip().lower()


def is_missing_like(value: Any) -> bool:
    return normalize_missing_like(value) in MISSING_LIKE_TOKENS


def parse_int(value: Any, label: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise DrugClassProfileV1Error(
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
        raise DrugClassProfileV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise DrugClassProfileV1Error(
            f"{label} is missing required columns: {sorted(missing)}"
        )


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
        raise DrugClassProfileV1Error(f"{label} does not report status == completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise DrugClassProfileV1Error(f"{label} does not report validation.passed == true.")
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
            raise DrugClassProfileV1Error(f"Blank {key_field} encountered in {label}.")
        if key in lookup:
            raise DrugClassProfileV1Error(
                f"Duplicate {key_field} '{key}' encountered in {label}."
            )
        lookup[key] = row
    return lookup


def has_any_regimen_context(row: dict[str, str]) -> bool:
    return str(row.get("regimen_context_values_json", "")).strip() not in {"", "[]"}


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    import csv

    if path.exists():
        raise DrugClassProfileV1Error(
            f"Refusing to overwrite existing output file: {path}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# Bucket classification
# ---------------------------------------------------------------------------


def classify_bucket(drug_class_row: dict[str, str]) -> str:
    """Classify a treated patient into one of 5 drug-class profile buckets.

    Definitions come ONLY from the saved patient_drug_class_profile_v1.tsv fields —
    no fresh raw re-derivation of drug names or classes.

    missing_like_only: single_or_multi_drug_class == unknown_only
                       AND contains_missing_like_raw_drug_name flag is set.
                       These patients have drug rows but all raw names are missing-like strings.
    unknown_only:      single_or_multi_drug_class == unknown_only AND not missing_like_only.
    single_known_class_clean: single_class AND both manual review flags == no.
    single_known_class_plus_unknown: single_class AND either manual review flag == yes.
    multi_known_class: multi_class (any number of distinct known classes >= 2).
    """
    flags_raw = drug_class_row.get("drug_name_profile_flags_json", "[]")
    try:
        flags = set(json.loads(flags_raw))
    except (json.JSONDecodeError, TypeError):
        flags = set()

    som = str(drug_class_row.get("single_or_multi_drug_class", "")).strip()
    norm_review = str(drug_class_row.get("drug_name_normalization_requires_manual_review", "no")).strip()
    class_review = str(drug_class_row.get("drug_class_mapping_requires_manual_review", "no")).strip()

    if som == "unknown_only":
        if "contains_missing_like_raw_drug_name" in flags:
            return "missing_like_only"
        return "unknown_only"
    if som == "single_class":
        if norm_review == "no" and class_review == "no":
            return "single_known_class_clean"
        return "single_known_class_plus_unknown"
    if som == "multi_class":
        return "multi_known_class"
    # Fallback for unexpected single_or_multi_drug_class values
    return "unknown_only"


def extract_known_classes(drug_class_row: dict[str, str]) -> list[str]:
    """Return the list of known (non-unknown) provisional drug classes for a patient."""
    raw = drug_class_row.get("provisional_drug_class_values_json", "[]")
    try:
        classes = json.loads(raw)
    except (json.JSONDecodeError, TypeError):
        classes = []
    return [c for c in classes if c != "unknown_or_review_needed" and c in KNOWN_CLASS_VOCABULARY]


# ---------------------------------------------------------------------------
# Workflow setup
# ---------------------------------------------------------------------------


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    if not trial_config.exists():
        raise DrugClassProfileV1Error(f"Required trial config not found: {trial_config}")

    trial_config_data = helper_module.load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    treatment_prep_audit_root = audit_root / "tcga-brca" / "treatment-prep"
    analysis_prep_audit_root = audit_root / "tcga-brca" / "analysis-prep"

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        audit_runs_root=treatment_prep_audit_root / "drug_class_profile_v1_runs",
        latest_pointer=treatment_prep_audit_root / "tcga_brca_drug_class_profile_v1_latest.json",
        drug_norm_latest_pointer=treatment_prep_audit_root
        / "tcga_brca_drug_name_normalization_v1_latest.json",
        grouping_latest_pointer=treatment_prep_audit_root
        / "tcga_brca_patient_treatment_grouping_v1_latest.json",
        baseline_analysis_latest_pointer=analysis_prep_audit_root
        / "tcga_brca_baseline_analysis_v1_latest.json",
    )


def load_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    required_pointers = [
        (paths.drug_norm_latest_pointer, "drug-name normalization v1 latest pointer"),
        (paths.grouping_latest_pointer, "patient treatment grouping v1 latest pointer"),
        (paths.baseline_analysis_latest_pointer, "baseline analysis v1 latest pointer"),
    ]
    for pointer_path, label in required_pointers:
        if not pointer_path.exists():
            raise DrugClassProfileV1Error(f"Required {label} not found: {pointer_path}")

    drug_norm_pointer = helper_module.load_json(paths.drug_norm_latest_pointer)
    helper_module.require_keys(
        drug_norm_pointer,
        {
            "drug_name_normalization_v1_run_id",
            "patient_treatment_grouping_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "patient_drug_class_profile_v1_tsv",
            "run_log_json",
        },
        "drug-name normalization v1 latest pointer",
        paths.drug_norm_latest_pointer,
    )
    grouping_pointer = helper_module.load_json(paths.grouping_latest_pointer)
    helper_module.require_keys(
        grouping_pointer,
        {
            "patient_treatment_grouping_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "patient_treatment_grouping_v1_tsv",
            "run_log_json",
        },
        "patient treatment grouping v1 latest pointer",
        paths.grouping_latest_pointer,
    )
    baseline_analysis_pointer = helper_module.load_json(paths.baseline_analysis_latest_pointer)
    helper_module.require_keys(
        baseline_analysis_pointer,
        {
            "baseline_analysis_v1_run_id",
            "cohort_v1_build_id",
            "baseline_analysis_v1_tsv",
            "run_log_json",
        },
        "baseline analysis v1 latest pointer",
        paths.baseline_analysis_latest_pointer,
    )

    drug_norm_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(drug_norm_pointer["run_log_json"]),
        label="drug-name normalization v1 run log",
        helper_module=helper_module,
    )
    grouping_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(grouping_pointer["run_log_json"]),
        label="patient treatment grouping v1 run log",
        helper_module=helper_module,
    )
    baseline_analysis_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(baseline_analysis_pointer["run_log_json"]),
        label="baseline analysis v1 run log",
        helper_module=helper_module,
    )

    # Lineage cross-checks
    if str(drug_norm_pointer["patient_treatment_grouping_v1_run_id"]) != str(
        grouping_pointer["patient_treatment_grouping_v1_run_id"]
    ):
        raise DrugClassProfileV1Error(
            "Mismatch: drug-norm pointer and grouping pointer disagree on "
            "patient_treatment_grouping_v1_run_id."
        )
    if str(drug_norm_pointer["os_endpoint_v1_run_id"]) != str(
        grouping_pointer["os_endpoint_v1_run_id"]
    ):
        raise DrugClassProfileV1Error(
            "Mismatch: drug-norm pointer and grouping pointer disagree on os_endpoint_v1_run_id."
        )
    if str(drug_norm_pointer["baseline_model_input_v1_run_id"]) != str(
        grouping_pointer["baseline_model_input_v1_run_id"]
    ):
        raise DrugClassProfileV1Error(
            "Mismatch: drug-norm pointer and grouping pointer disagree on "
            "baseline_model_input_v1_run_id."
        )
    if str(drug_norm_pointer["cohort_v1_build_id"]) != str(grouping_pointer["cohort_v1_build_id"]):
        raise DrugClassProfileV1Error(
            "Mismatch: drug-norm pointer and grouping pointer disagree on cohort_v1_build_id."
        )
    if str(drug_norm_pointer["cohort_v1_build_id"]) != str(
        baseline_analysis_pointer["cohort_v1_build_id"]
    ):
        raise DrugClassProfileV1Error(
            "Mismatch: drug-norm pointer and baseline analysis pointer disagree on cohort_v1_build_id."
        )

    # Load TSV files
    drug_class_profile_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(drug_norm_pointer["patient_drug_class_profile_v1_tsv"]),
        "patient_drug_class_profile_v1.tsv",
    )
    grouping_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(grouping_pointer["patient_treatment_grouping_v1_tsv"]),
        "patient_treatment_grouping_v1.tsv",
    )
    baseline_analysis_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(baseline_analysis_pointer["baseline_analysis_v1_tsv"]),
        "baseline_analysis_v1.tsv",
    )

    drug_class_profile_rows = helper_module.read_tsv_dict_rows(drug_class_profile_tsv)
    grouping_rows = helper_module.read_tsv_dict_rows(grouping_tsv)
    baseline_analysis_rows = helper_module.read_tsv_dict_rows(baseline_analysis_tsv)

    require_columns(
        drug_class_profile_rows,
        {
            "bcr_patient_barcode",
            "single_or_multi_drug_class",
            "drug_name_normalization_requires_manual_review",
            "drug_class_mapping_requires_manual_review",
            "drug_name_profile_flags_json",
            "provisional_drug_class_values_json",
            "dominant_provisional_drug_class_if_any",
            "regimen_context_values_json",
        },
        "patient_drug_class_profile_v1.tsv",
    )
    require_columns(
        grouping_rows,
        {"bcr_patient_barcode", "os_event", "os_time_days", "has_any_radiation_row",
         "has_any_treatment_timing"},
        "patient_treatment_grouping_v1.tsv",
    )
    require_columns(
        baseline_analysis_rows,
        {"bcr_patient_barcode"} | set(REQUIRED_CONFOUNDER_FIELDS),
        "baseline_analysis_v1.tsv",
    )

    return WorkflowInputs(
        drug_norm_pointer=drug_norm_pointer,
        drug_norm_run_log=drug_norm_run_log,
        grouping_pointer=grouping_pointer,
        grouping_run_log=grouping_run_log,
        baseline_analysis_pointer=baseline_analysis_pointer,
        baseline_analysis_run_log=baseline_analysis_run_log,
        drug_class_profile_rows=drug_class_profile_rows,
        grouping_rows=grouping_rows,
        baseline_analysis_rows=baseline_analysis_rows,
        input_paths={
            "drug_class_profile_tsv": drug_class_profile_tsv,
            "grouping_tsv": grouping_tsv,
            "baseline_analysis_tsv": baseline_analysis_tsv,
        },
    )


# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------


def build_merged_rows(inputs: WorkflowInputs) -> list[dict[str, str]]:
    """Join drug-class profile, grouping, and baseline analysis on bcr_patient_barcode.

    Returns one merged dict per treated patient (780 rows).
    Validates that all 780 drug-class-profile barcodes resolve in both join sources.
    """
    grouping_lookup = build_unique_lookup(
        inputs.grouping_rows,
        key_field="bcr_patient_barcode",
        label="patient_treatment_grouping_v1.tsv",
        normalize_key=True,
    )
    baseline_lookup = build_unique_lookup(
        inputs.baseline_analysis_rows,
        key_field="bcr_patient_barcode",
        label="baseline_analysis_v1.tsv",
        normalize_key=True,
    )

    missing_in_grouping: list[str] = []
    missing_in_baseline: list[str] = []
    merged: list[dict[str, str]] = []

    for dp_row in inputs.drug_class_profile_rows:
        barcode = normalize_barcode(str(dp_row.get("bcr_patient_barcode", "")))
        if not barcode:
            raise DrugClassProfileV1Error(
                "Blank bcr_patient_barcode found in patient_drug_class_profile_v1.tsv."
            )

        g_row = grouping_lookup.get(barcode)
        b_row = baseline_lookup.get(barcode)

        if g_row is None:
            missing_in_grouping.append(barcode)
        if b_row is None:
            missing_in_baseline.append(barcode)

        merged.append(
            {
                **dp_row,
                # OS and treatment context from grouping
                "os_event": g_row.get("os_event", "") if g_row else "",
                "os_time_days": g_row.get("os_time_days", "") if g_row else "",
                "has_any_radiation_row": g_row.get("has_any_radiation_row", "") if g_row else "",
                "has_any_treatment_timing": g_row.get("has_any_treatment_timing", "") if g_row else "",
                # Clinical confounders from baseline analysis
                **{
                    field: b_row.get(field, "") if b_row else ""
                    for field in REQUIRED_CONFOUNDER_FIELDS
                },
            }
        )

    if missing_in_grouping:
        raise DrugClassProfileV1Error(
            f"{len(missing_in_grouping)} drug-class-profile barcodes not found in "
            f"patient_treatment_grouping_v1.tsv: {missing_in_grouping[:5]}"
        )
    if missing_in_baseline:
        raise DrugClassProfileV1Error(
            f"{len(missing_in_baseline)} drug-class-profile barcodes not found in "
            f"baseline_analysis_v1.tsv: {missing_in_baseline[:5]}"
        )

    return merged


# ---------------------------------------------------------------------------
# Major class identification
# ---------------------------------------------------------------------------


def identify_major_classes(merged_rows: list[dict[str, str]]) -> list[str]:
    """Return known drug classes with N >= MIN_CLASS_SIZE_FOR_AUDIT patients.

    A patient contributes to a class count if that class appears in their
    provisional_drug_class_values_json (any-class membership, not dominant-only).
    Classes are returned in KNOWN_CLASS_VOCABULARY order.
    """
    class_patient_counts: Counter[str] = Counter()
    for row in merged_rows:
        known = extract_known_classes(row)
        for cls in set(known):  # count each patient once per class
            class_patient_counts[cls] += 1

    return [
        cls
        for cls in KNOWN_CLASS_VOCABULARY
        if class_patient_counts[cls] >= MIN_CLASS_SIZE_FOR_AUDIT
    ]


def identify_dominant_class_candidates(merged_rows: list[dict[str, str]]) -> list[str]:
    """Return classes with N >= MIN_CLASS_SIZE_FOR_CANDIDATE patients where
    dominant_provisional_drug_class_if_any equals that class.
    Returned in KNOWN_CLASS_VOCABULARY order."""
    dominant_counts: Counter[str] = Counter()
    for row in merged_rows:
        dom = str(row.get("dominant_provisional_drug_class_if_any", "")).strip()
        if dom and dom in KNOWN_CLASS_VOCABULARY:
            dominant_counts[dom] += 1
    return [
        cls
        for cls in KNOWN_CLASS_VOCABULARY
        if dominant_counts[cls] >= MIN_CLASS_SIZE_FOR_CANDIDATE
    ]


# ---------------------------------------------------------------------------
# Group summary
# ---------------------------------------------------------------------------


def build_group_summary_rows(
    merged_rows: list[dict[str, str]],
    bucket_map: dict[str, str],
    run_id: str,
) -> list[dict[str, str]]:
    total = len(merged_rows)
    bucket_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in merged_rows:
        barcode = normalize_barcode(str(row.get("bcr_patient_barcode", "")))
        bucket = bucket_map.get(barcode, "unknown_only")
        bucket_rows[bucket].append(row)

    output_rows: list[dict[str, str]] = []
    for order_idx, bucket in enumerate(APPROVED_BUCKET_ORDER, start=1):
        brows = bucket_rows.get(bucket, [])
        n = len(brows)
        os_events = sum(1 for r in brows if str(r.get("os_event", "0")).strip() == "1")
        manual = sum(
            1
            for r in brows
            if (
                str(r.get("drug_name_normalization_requires_manual_review", "no")).strip() == "yes"
                or str(r.get("drug_class_mapping_requires_manual_review", "no")).strip() == "yes"
            )
        )
        radiation = sum(
            1 for r in brows if str(r.get("has_any_radiation_row", "no")).strip() == "yes"
        )
        regimen = sum(1 for r in brows if has_any_regimen_context(r))
        timing = sum(
            1 for r in brows if str(r.get("has_any_treatment_timing", "no")).strip() == "yes"
        )
        output_rows.append(
            {
                "drug_class_profile_v1_run_id": run_id,
                "bucket_order": str(order_idx),
                "drug_class_profile_bucket": bucket,
                "patient_count": str(n),
                "patient_fraction": format_fraction(n, total),
                "os_event_count": str(os_events),
                "os_event_fraction": format_fraction(os_events, n),
                "manual_review_count": str(manual),
                "manual_review_fraction": format_fraction(manual, n),
                "radiation_overlap_count": str(radiation),
                "radiation_overlap_fraction": format_fraction(radiation, n),
                "regimen_context_count": str(regimen),
                "regimen_context_fraction": format_fraction(regimen, n),
                "timing_coverage_count": str(timing),
                "timing_coverage_fraction": format_fraction(timing, n),
                "notes": "" if n > 0 else "no patients in this bucket",
            }
        )
    return output_rows


# ---------------------------------------------------------------------------
# Confounder coverage
# ---------------------------------------------------------------------------


def _coverage_row(
    run_id: str,
    review_group_type: str,
    field_name: str,
    rows: list[dict[str, str]],
    note: str = "",
) -> dict[str, str]:
    n = len(rows)
    values = [str(r.get(field_name, "")).strip() for r in rows]
    non_missing = [v for v in values if not is_missing_like(v)]
    missing_count = n - len(non_missing)
    return {
        "drug_class_profile_v1_run_id": run_id,
        "review_group_type": review_group_type,
        "field_name": field_name,
        "row_count": str(n),
        "non_missing_count": str(len(non_missing)),
        "missing_like_count": str(missing_count),
        "missing_like_fraction": format_fraction(missing_count, n),
        "distinct_non_missing_count": str(len(set(non_missing))),
        "notes": note,
    }


def build_confounder_coverage_rows(
    merged_rows: list[dict[str, str]],
    bucket_map: dict[str, str],
    major_classes: list[str],
    run_id: str,
) -> list[dict[str, str]]:
    bucket_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in merged_rows:
        barcode = normalize_barcode(str(row.get("bcr_patient_barcode", "")))
        bucket = bucket_map.get(barcode, "unknown_only")
        bucket_rows[bucket].append(row)

    class_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in merged_rows:
        for cls in extract_known_classes(row):
            if cls in major_classes:
                class_rows[cls].append(row)

    output: list[dict[str, str]] = []

    # all_treated
    for field in REQUIRED_CONFOUNDER_FIELDS:
        output.append(_coverage_row(run_id, "all_treated", field, merged_rows))

    # per bucket
    for bucket in APPROVED_BUCKET_ORDER:
        brows = bucket_rows.get(bucket, [])
        for field in REQUIRED_CONFOUNDER_FIELDS:
            output.append(_coverage_row(run_id, bucket, field, brows))

    # per major class (any-class membership)
    for cls in major_classes:
        crows = class_rows.get(cls, [])
        for field in REQUIRED_CONFOUNDER_FIELDS:
            output.append(
                _coverage_row(
                    run_id,
                    f"any_class_{cls}",
                    field,
                    crows,
                    note="patients where class appears in provisional_drug_class_values_json",
                )
            )

    return output


# ---------------------------------------------------------------------------
# Value composition
# ---------------------------------------------------------------------------


def _value_composition_categorical(
    run_id: str,
    review_group_type: str,
    field_name: str,
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    n = len(rows)
    values = [str(r.get(field_name, "")).strip() for r in rows]
    non_missing = [v for v in values if not is_missing_like(v)]
    missing_count = n - len(non_missing)
    distinct = len(set(non_missing))
    counter = Counter(non_missing)
    output: list[dict[str, str]] = []
    for rank, (label, count) in enumerate(counter.most_common(MAX_TOP_VALUES), start=1):
        output.append(
            {
                "drug_class_profile_v1_run_id": run_id,
                "review_group_type": review_group_type,
                "field_name": field_name,
                "value_summary_type": "categorical_top_value",
                "group_patient_count": str(n),
                "non_missing_count": str(len(non_missing)),
                "missing_like_count": str(missing_count),
                "distinct_non_missing_count": str(distinct),
                "value_rank": str(rank),
                "value_label": label,
                "value_count": str(count),
                "value_fraction": format_fraction(count, len(non_missing)),
                "numeric_min": "",
                "numeric_median": "",
                "numeric_max": "",
                "notes": "",
            }
        )
    if not output:
        # emit a placeholder row if no non-missing values
        output.append(
            {
                "drug_class_profile_v1_run_id": run_id,
                "review_group_type": review_group_type,
                "field_name": field_name,
                "value_summary_type": "categorical_top_value",
                "group_patient_count": str(n),
                "non_missing_count": "0",
                "missing_like_count": str(missing_count),
                "distinct_non_missing_count": "0",
                "value_rank": "",
                "value_label": "",
                "value_count": "",
                "value_fraction": "",
                "numeric_min": "",
                "numeric_median": "",
                "numeric_max": "",
                "notes": "no non-missing values in this group",
            }
        )
    return output


def _value_composition_age(
    run_id: str,
    review_group_type: str,
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    field_name = "age_at_diagnosis"
    n = len(rows)
    raw_values = [str(r.get(field_name, "")).strip() for r in rows]
    non_missing_str = [v for v in raw_values if not is_missing_like(v)]
    missing_count = n - len(non_missing_str)
    numeric_values = [parse_float_or_none(v) for v in non_missing_str]
    numeric_values_clean = [v for v in numeric_values if v is not None]

    if numeric_values_clean:
        num_min = min(numeric_values_clean)
        num_max = max(numeric_values_clean)
        num_median = statistics.median(numeric_values_clean)
    else:
        num_min = num_max = num_median = None

    return [
        {
            "drug_class_profile_v1_run_id": run_id,
            "review_group_type": review_group_type,
            "field_name": field_name,
            "value_summary_type": "numeric_summary",
            "group_patient_count": str(n),
            "non_missing_count": str(len(non_missing_str)),
            "missing_like_count": str(missing_count),
            "distinct_non_missing_count": str(len(set(non_missing_str))),
            "value_rank": "",
            "value_label": "",
            "value_count": "",
            "value_fraction": "",
            "numeric_min": format_numeric(num_min),
            "numeric_median": format_numeric(num_median),
            "numeric_max": format_numeric(num_max),
            "notes": "",
        }
    ]


def _compose_for_group(
    run_id: str,
    review_group_type: str,
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    output.extend(_value_composition_age(run_id, review_group_type, rows))
    for field in CATEGORICAL_COMPOSITION_FIELDS:
        output.extend(_value_composition_categorical(run_id, review_group_type, field, rows))
    return output


def build_value_composition_rows(
    merged_rows: list[dict[str, str]],
    bucket_map: dict[str, str],
    major_classes: list[str],
    run_id: str,
) -> list[dict[str, str]]:
    bucket_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in merged_rows:
        barcode = normalize_barcode(str(row.get("bcr_patient_barcode", "")))
        bucket = bucket_map.get(barcode, "unknown_only")
        bucket_rows[bucket].append(row)

    class_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in merged_rows:
        for cls in extract_known_classes(row):
            if cls in major_classes:
                class_rows[cls].append(row)

    output: list[dict[str, str]] = []
    output.extend(_compose_for_group(run_id, "all_treated", merged_rows))
    for bucket in APPROVED_BUCKET_ORDER:
        brows = bucket_rows.get(bucket, [])
        output.extend(_compose_for_group(run_id, bucket, brows))
    for cls in major_classes:
        crows = class_rows.get(cls, [])
        output.extend(_compose_for_group(run_id, f"any_class_{cls}", crows))
    return output


# ---------------------------------------------------------------------------
# Arm-freeze candidate audit
# ---------------------------------------------------------------------------


def _check_confounder_coverage_ok(rows: list[dict[str, str]]) -> str:
    """Return 'yes' if all COVERAGE_READY_FIELDS meet CONFOUNDER_NON_MISSING_TOLERANCE."""
    if not rows:
        return "no"
    for field in COVERAGE_READY_FIELDS:
        values = [str(r.get(field, "")).strip() for r in rows]
        non_missing = [v for v in values if not is_missing_like(v)]
        fraction = len(non_missing) / len(rows) if rows else 0.0
        if fraction < CONFOUNDER_NON_MISSING_TOLERANCE:
            return "no"
    return "yes"


def _arm_freeze_candidate_status(
    patient_count: int,
    manual_review_fraction: float,
    confounder_coverage_ok: str,
) -> str:
    if patient_count < MIN_CLASS_SIZE_FOR_CANDIDATE:
        return "not_recommended"
    if manual_review_fraction <= MANUAL_REVIEW_FRACTION_TOLERANCE and confounder_coverage_ok == "yes":
        return "promising"
    return "needs_review"


def _candidate_row(
    run_id: str,
    candidate_label: str,
    candidate_type: str,
    rows: list[dict[str, str]],
    why: str,
    recommended: str,
    note: str = "",
) -> dict[str, str]:
    n = len(rows)
    os_events = sum(1 for r in rows if str(r.get("os_event", "0")).strip() == "1")
    manual = sum(
        1
        for r in rows
        if (
            str(r.get("drug_name_normalization_requires_manual_review", "no")).strip() == "yes"
            or str(r.get("drug_class_mapping_requires_manual_review", "no")).strip() == "yes"
        )
    )
    unknown_burden_count = sum(
        1
        for r in rows
        if (
            str(r.get("drug_name_normalization_requires_manual_review", "no")).strip() == "yes"
            or str(r.get("drug_class_mapping_requires_manual_review", "no")).strip() == "yes"
        )
    )
    multi_class_count = sum(
        1
        for r in rows
        if str(r.get("single_or_multi_drug_class", "")).strip() == "multi_class"
    )
    manual_frac = unknown_burden_count / n if n > 0 else 0.0
    multi_frac = multi_class_count / n if n > 0 else 0.0
    coverage_ok = _check_confounder_coverage_ok(rows)
    status = _arm_freeze_candidate_status(n, manual_frac, coverage_ok)
    return {
        "drug_class_profile_v1_run_id": run_id,
        "candidate_label": candidate_label,
        "candidate_type": candidate_type,
        "patient_count": str(n),
        "os_event_count": str(os_events),
        "manual_review_count": str(manual),
        "unknown_mapping_burden": format_float(manual_frac),
        "multi_class_burden": format_float(multi_frac),
        "confounder_coverage_ok": coverage_ok,
        "arm_freeze_candidate_status": status,
        "why_candidate_or_not": why,
        "recommended_next_handling": recommended,
        "notes": note,
    }


def build_arm_freeze_candidate_audit_rows(
    merged_rows: list[dict[str, str]],
    bucket_map: dict[str, str],
    major_classes: list[str],
    run_id: str,
) -> list[dict[str, str]]:
    # Build per-bucket and per-dominant-class lookups
    bucket_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in merged_rows:
        barcode = normalize_barcode(str(row.get("bcr_patient_barcode", "")))
        bucket = bucket_map.get(barcode, "unknown_only")
        bucket_rows[bucket].append(row)

    dominant_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in merged_rows:
        dom = str(row.get("dominant_provisional_drug_class_if_any", "")).strip()
        if dom and dom in KNOWN_CLASS_VOCABULARY:
            dominant_rows[dom].append(row)

    output: list[dict[str, str]] = []

    # Candidate 1: all single_known_class_clean patients
    clean_rows = bucket_rows.get("single_known_class_clean", [])
    n_clean = len(clean_rows)
    output.append(
        _candidate_row(
            run_id,
            "single_known_class_clean_all",
            "bucket",
            clean_rows,
            why=(
                f"Single confirmed drug class, no manual-review flags. "
                f"n={n_clean}. Cleanest bucket by definition."
            ),
            recommended=(
                "Assess for arm-freeze after resolving dominant class distribution within bucket."
                if n_clean >= MIN_CLASS_SIZE_FOR_CANDIDATE
                else "Too small to support arm-freeze review at current size."
            ),
        )
    )

    # Candidate 2: each major dominant class (using dominant_provisional_drug_class_if_any)
    dominant_candidate_classes = identify_dominant_class_candidates(merged_rows)
    for cls in KNOWN_CLASS_VOCABULARY:
        dom_rows = dominant_rows.get(cls, [])
        n_dom = len(dom_rows)
        if n_dom == 0:
            continue
        reason_size = (
            f"dominant_provisional_drug_class_if_any == {cls!r} for {n_dom} treated patients."
        )
        if n_dom < MIN_CLASS_SIZE_FOR_CANDIDATE:
            why = f"{reason_size} Below size threshold ({MIN_CLASS_SIZE_FOR_CANDIDATE})."
            rec = "Too small for arm-freeze review at current size; consider aggregating with related class."
        else:
            why = f"{reason_size} Above size threshold."
            rec = (
                "Assess for arm-freeze after manual-review resolution."
                if _arm_freeze_candidate_status(
                    n_dom,
                    sum(
                        1
                        for r in dom_rows
                        if str(r.get("drug_name_normalization_requires_manual_review", "no")).strip()
                        == "yes"
                        or str(r.get("drug_class_mapping_requires_manual_review", "no")).strip()
                        == "yes"
                    )
                    / n_dom,
                    _check_confounder_coverage_ok(dom_rows),
                )
                == "needs_review"
                else "Promising dominant-class subset; review for arm freeze."
            )
        output.append(
            _candidate_row(
                run_id,
                f"dominant_class_{cls}",
                "dominant_class",
                dom_rows,
                why=why,
                recommended=rec,
            )
        )

    # Candidate 3: intersection of single_known_class_clean AND dominant class
    for cls in KNOWN_CLASS_VOCABULARY:
        clean_dom_rows = [
            r
            for r in clean_rows
            if str(r.get("dominant_provisional_drug_class_if_any", "")).strip() == cls
        ]
        n_cd = len(clean_dom_rows)
        if n_cd == 0:
            continue
        output.append(
            _candidate_row(
                run_id,
                f"clean_bucket_x_dominant_{cls}",
                "bucket_x_dominant_class",
                clean_dom_rows,
                why=(
                    f"single_known_class_clean patients where dominant class == {cls!r}. "
                    f"n={n_cd}. Cleanest possible subset for this class."
                ),
                recommended=(
                    "Strong candidate for arm-freeze review if n >= threshold and OS events adequate."
                    if n_cd >= MIN_CLASS_SIZE_FOR_CANDIDATE
                    else "Too small for independent arm-freeze review; may support pooled analysis."
                ),
            )
        )

    return output


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def _s(
    run_id: str,
    section: str,
    metric: str,
    value: Any,
    note: str = "",
) -> dict[str, str]:
    return {
        "drug_class_profile_v1_run_id": run_id,
        "summary_section": section,
        "summary_metric": metric,
        "summary_value": str(value),
        "notes": note,
    }


def build_summary_rows(
    merged_rows: list[dict[str, str]],
    bucket_map: dict[str, str],
    major_classes: list[str],
    audit_rows: list[dict[str, str]],
    inputs: WorkflowInputs,
    run_id: str,
) -> list[dict[str, str]]:
    total = len(merged_rows)
    bucket_counts: Counter[str] = Counter()
    for barcode, bucket in bucket_map.items():
        bucket_counts[bucket] += 1

    # Per-class any-membership counts
    class_patient_counts: Counter[str] = Counter()
    for row in merged_rows:
        for cls in set(extract_known_classes(row)):
            class_patient_counts[cls] += 1

    # Per-dominant-class counts
    dominant_counts: Counter[str] = Counter()
    for row in merged_rows:
        dom = str(row.get("dominant_provisional_drug_class_if_any", "")).strip()
        if dom and dom in KNOWN_CLASS_VOCABULARY:
            dominant_counts[dom] += 1

    # Candidate stats
    promising = [r for r in audit_rows if r.get("arm_freeze_candidate_status") == "promising"]
    needs_review = [r for r in audit_rows if r.get("arm_freeze_candidate_status") == "needs_review"]

    # Readiness
    readiness = (
        READINESS_ARM_FREEZE_CANDIDATE_DISCUSSION
        if promising
        else READINESS_DESCRIPTIVE_COMPLETED
    )
    freeze_status = (
        TREATMENT_ARM_FREEZE_REVIEW_PROVISIONAL_NEXT if promising else TREATMENT_ARM_FREEZE_REVIEW_BLOCKED
    )

    rows: list[dict[str, str]] = []

    # inputs section
    rows.append(_s(run_id, "inputs", "drug_name_normalization_v1_run_id",
                   inputs.drug_norm_pointer.get("drug_name_normalization_v1_run_id", "")))
    rows.append(_s(run_id, "inputs", "patient_treatment_grouping_v1_run_id",
                   inputs.grouping_pointer.get("patient_treatment_grouping_v1_run_id", "")))
    rows.append(_s(run_id, "inputs", "os_endpoint_v1_run_id",
                   inputs.drug_norm_pointer.get("os_endpoint_v1_run_id", "")))
    rows.append(_s(run_id, "inputs", "baseline_model_input_v1_run_id",
                   inputs.drug_norm_pointer.get("baseline_model_input_v1_run_id", "")))
    rows.append(_s(run_id, "inputs", "baseline_analysis_v1_run_id",
                   inputs.baseline_analysis_pointer.get("baseline_analysis_v1_run_id", "")))
    rows.append(_s(run_id, "inputs", "cohort_v1_build_id",
                   inputs.drug_norm_pointer.get("cohort_v1_build_id", "")))

    # coverage section
    rows.append(_s(run_id, "coverage", "total_treated_patients_reviewed", total))
    rows.append(_s(run_id, "coverage", "total_drug_class_profile_rows", total))
    rows.append(_s(run_id, "coverage", "distinct_provisional_classes_in_vocabulary",
                   len(PROVISIONAL_CLASS_VOCABULARY)))
    rows.append(_s(run_id, "coverage", "distinct_known_classes_in_vocabulary",
                   len(KNOWN_CLASS_VOCABULARY)))
    rows.append(_s(run_id, "coverage", "major_classes_min_audit_threshold",
                   MIN_CLASS_SIZE_FOR_AUDIT))
    rows.append(_s(run_id, "coverage", "major_classes_count", len(major_classes)))
    rows.append(_s(run_id, "coverage", "major_classes_list", json.dumps(major_classes)))

    # patients section — buckets
    for bucket in APPROVED_BUCKET_ORDER:
        rows.append(_s(run_id, "patients", f"{bucket}_count", bucket_counts.get(bucket, 0)))
    rows.append(_s(run_id, "patients", "bucket_count_sum_check",
                   sum(bucket_counts.values()),
                   note=f"should equal {total}"))

    # patients section — class any-membership
    for cls in KNOWN_CLASS_VOCABULARY:
        rows.append(
            _s(run_id, "classes", f"any_class_{cls}_patient_count",
               class_patient_counts.get(cls, 0),
               note="patients where class appears in provisional_drug_class_values_json")
        )
    for cls in KNOWN_CLASS_VOCABULARY:
        rows.append(
            _s(run_id, "classes", f"dominant_class_{cls}_patient_count",
               dominant_counts.get(cls, 0),
               note="patients where dominant_provisional_drug_class_if_any == class")
        )

    # aggregate class metrics
    rows.append(_s(run_id, "classes", "classes_with_any_membership_ge_50",
                   sum(1 for c in KNOWN_CLASS_VOCABULARY if class_patient_counts[c] >= 50)))
    rows.append(_s(run_id, "classes", "classes_with_dominant_membership_ge_50",
                   sum(1 for c in KNOWN_CLASS_VOCABULARY if dominant_counts[c] >= 50)))

    # candidates section
    rows.append(_s(run_id, "candidates", "total_candidates_assessed", len(audit_rows)))
    rows.append(_s(run_id, "candidates", "promising_candidates_count", len(promising)))
    rows.append(_s(run_id, "candidates", "needs_review_candidates_count", len(needs_review)))
    rows.append(_s(run_id, "candidates", "not_recommended_candidates_count",
                   len(audit_rows) - len(promising) - len(needs_review)))
    if promising:
        rows.append(_s(run_id, "candidates", "promising_candidate_labels",
                       json.dumps([r["candidate_label"] for r in promising])))

    # readiness section
    rows.append(_s(run_id, "readiness", "provisional_readiness_interpretation", readiness,
                   note="descriptive review status only; does not freeze treatment arms"))
    rows.append(_s(run_id, "readiness", "treatment_arm_freeze_review_status", freeze_status,
                   note="blocked = not yet ready; provisional_review_next_only = promising candidates found"))
    rows.append(_s(run_id, "readiness", "treatment_recommendation_modeling_status",
                   TREATMENT_RECOMMENDATION_MODELING_STATUS,
                   note="treatment recommendation modeling remains out of scope"))

    return rows


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------


def write_outputs(
    paths: WorkflowPaths,
    run_id: str,
    group_summary_rows: list[dict[str, str]],
    confounder_coverage_rows: list[dict[str, str]],
    value_composition_rows: list[dict[str, str]],
    arm_freeze_candidate_audit_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    inputs: WorkflowInputs,
    started_at: datetime,
) -> None:
    audit_run_dir = paths.audit_runs_root / run_id
    if audit_run_dir.exists():
        raise DrugClassProfileV1Error(
            f"Audit run directory already exists (no silent overwrite): {audit_run_dir}"
        )
    audit_run_dir.mkdir(parents=True, exist_ok=False)

    group_summary_tsv = audit_run_dir / "drug_class_profile_v1_group_summary.tsv"
    confounder_coverage_tsv = audit_run_dir / "drug_class_profile_v1_confounder_coverage.tsv"
    value_composition_tsv = audit_run_dir / "drug_class_profile_v1_value_composition.tsv"
    arm_freeze_tsv = audit_run_dir / "drug_class_profile_v1_arm_freeze_candidate_audit.tsv"
    summary_tsv = audit_run_dir / "drug_class_profile_v1_summary.tsv"
    run_log_path = audit_run_dir / "run_log.json"

    write_tsv(group_summary_tsv, GROUP_SUMMARY_FIELDNAMES, group_summary_rows)
    write_tsv(confounder_coverage_tsv, CONFOUNDER_COVERAGE_FIELDNAMES, confounder_coverage_rows)
    write_tsv(value_composition_tsv, VALUE_COMPOSITION_FIELDNAMES, value_composition_rows)
    write_tsv(arm_freeze_tsv, ARM_FREEZE_CANDIDATE_AUDIT_FIELDNAMES, arm_freeze_candidate_audit_rows)
    write_tsv(summary_tsv, SUMMARY_FIELDNAMES, summary_rows)

    completed_at = utc_now()

    # Derive counts from summary rows for run_log
    bucket_counts = {
        bucket: int(next(
            (r["summary_value"] for r in summary_rows
             if r["summary_metric"] == f"{bucket}_count"),
            "0",
        ))
        for bucket in APPROVED_BUCKET_ORDER
    }
    promising_count = int(next(
        (r["summary_value"] for r in summary_rows
         if r["summary_metric"] == "promising_candidates_count"),
        "0",
    ))

    run_log: dict[str, Any] = {
        "status": "completed",
        "drug_class_profile_v1_run_id": run_id,
        "drug_name_normalization_v1_run_id": str(
            inputs.drug_norm_pointer.get("drug_name_normalization_v1_run_id", "")
        ),
        "patient_treatment_grouping_v1_run_id": str(
            inputs.grouping_pointer.get("patient_treatment_grouping_v1_run_id", "")
        ),
        "os_endpoint_v1_run_id": str(
            inputs.drug_norm_pointer.get("os_endpoint_v1_run_id", "")
        ),
        "baseline_model_input_v1_run_id": str(
            inputs.drug_norm_pointer.get("baseline_model_input_v1_run_id", "")
        ),
        "baseline_analysis_v1_run_id": str(
            inputs.baseline_analysis_pointer.get("baseline_analysis_v1_run_id", "")
        ),
        "cohort_v1_build_id": str(inputs.drug_norm_pointer.get("cohort_v1_build_id", "")),
        "started_at_utc": format_utc_timestamp(started_at),
        "completed_at_utc": format_utc_timestamp(completed_at),
        "repo_root": str(paths.repo_root),
        "trial_name": "tcga_only_source_audited",
        "dataset_scope": "tcga_brca_only",
        "inputs": {
            "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
            "drug_norm_latest_pointer": repo_relative(
                paths.drug_norm_latest_pointer, paths.repo_root
            ),
            "grouping_latest_pointer": repo_relative(
                paths.grouping_latest_pointer, paths.repo_root
            ),
            "baseline_analysis_latest_pointer": repo_relative(
                paths.baseline_analysis_latest_pointer, paths.repo_root
            ),
            "drug_class_profile_tsv": repo_relative(
                inputs.input_paths["drug_class_profile_tsv"], paths.repo_root
            ),
            "grouping_tsv": repo_relative(inputs.input_paths["grouping_tsv"], paths.repo_root),
            "baseline_analysis_tsv": repo_relative(
                inputs.input_paths["baseline_analysis_tsv"], paths.repo_root
            ),
        },
        "outputs": {
            "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
            "drug_class_profile_v1_group_summary_tsv": repo_relative(
                group_summary_tsv, paths.repo_root
            ),
            "drug_class_profile_v1_confounder_coverage_tsv": repo_relative(
                confounder_coverage_tsv, paths.repo_root
            ),
            "drug_class_profile_v1_value_composition_tsv": repo_relative(
                value_composition_tsv, paths.repo_root
            ),
            "drug_class_profile_v1_arm_freeze_candidate_audit_tsv": repo_relative(
                arm_freeze_tsv, paths.repo_root
            ),
            "drug_class_profile_v1_summary_tsv": repo_relative(summary_tsv, paths.repo_root),
            "run_log_json": repo_relative(run_log_path, paths.repo_root),
            "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
        },
        "validation": {
            "passed": True,
            "required_upstream_pointers_found": True,
            "required_source_tables_found": True,
            "drug_class_profile_row_count_positive": len(inputs.drug_class_profile_rows) > 0,
            "grouping_row_count_positive": len(inputs.grouping_rows) > 0,
            "baseline_analysis_row_count_positive": len(inputs.baseline_analysis_rows) > 0,
            "all_barcodes_resolved_in_grouping": True,
            "all_barcodes_resolved_in_baseline": True,
            "bucket_counts_sum_to_treated": sum(bucket_counts.values())
            == len(inputs.drug_class_profile_rows),
            "group_summary_row_count": len(group_summary_rows),
            "group_summary_row_count_is_5": len(group_summary_rows) == 5,
            "confounder_coverage_row_count": len(confounder_coverage_rows),
            "confounder_coverage_row_count_positive": len(confounder_coverage_rows) > 0,
            "value_composition_row_count": len(value_composition_rows),
            "value_composition_row_count_positive": len(value_composition_rows) > 0,
            "arm_freeze_candidate_audit_row_count": len(arm_freeze_candidate_audit_rows),
            "arm_freeze_candidate_audit_row_count_positive": len(arm_freeze_candidate_audit_rows)
            > 0,
            "summary_row_count": len(summary_rows),
            "no_silent_overwrite": True,
            "latest_pointer_updated_after_success": True,
        },
        "rules": {
            "unit_of_analysis": "treated patient (has_any_drug_row == yes)",
            "input_layer": "drug_name_normalization_v1 + patient_treatment_grouping_v1 + baseline_analysis_v1",
            "output_layer": "drug_class_profile_v1 descriptive review and arm-freeze candidate audit",
            "does_not_normalize_drug_names": True,
            "does_not_redefine_treatment_groups": True,
            "does_not_freeze_treatment_arms": True,
            "does_not_exclude_mixed_patients_silently": True,
            "does_not_perform_modeling": True,
            "bucket_definitions_from_saved_profile_only": True,
            "treatment_recommendation_modeling_status": TREATMENT_RECOMMENDATION_MODELING_STATUS,
        },
        "counts": {
            "total_treated_patients": len(inputs.drug_class_profile_rows),
            **{f"{bucket}_count": bucket_counts.get(bucket, 0) for bucket in APPROVED_BUCKET_ORDER},
            "promising_candidates_count": promising_count,
        },
        "upstream_snapshots": {
            "drug_norm_pointer": inputs.drug_norm_pointer,
            "grouping_pointer": inputs.grouping_pointer,
            "baseline_analysis_pointer": inputs.baseline_analysis_pointer,
        },
    }

    import csv as _csv

    with run_log_path.open("w", encoding="utf-8") as f:
        json.dump(run_log, f, indent=2)

    # Write latest pointer (only after all outputs succeed)
    latest_pointer_data: dict[str, Any] = {
        "updated_at_utc": format_utc_timestamp(completed_at),
        "drug_class_profile_v1_run_id": run_id,
        "drug_name_normalization_v1_run_id": str(
            inputs.drug_norm_pointer.get("drug_name_normalization_v1_run_id", "")
        ),
        "patient_treatment_grouping_v1_run_id": str(
            inputs.grouping_pointer.get("patient_treatment_grouping_v1_run_id", "")
        ),
        "os_endpoint_v1_run_id": str(
            inputs.drug_norm_pointer.get("os_endpoint_v1_run_id", "")
        ),
        "baseline_model_input_v1_run_id": str(
            inputs.drug_norm_pointer.get("baseline_model_input_v1_run_id", "")
        ),
        "baseline_analysis_v1_run_id": str(
            inputs.baseline_analysis_pointer.get("baseline_analysis_v1_run_id", "")
        ),
        "cohort_v1_build_id": str(inputs.drug_norm_pointer.get("cohort_v1_build_id", "")),
        "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
        "drug_class_profile_v1_group_summary_tsv": repo_relative(
            group_summary_tsv, paths.repo_root
        ),
        "drug_class_profile_v1_confounder_coverage_tsv": repo_relative(
            confounder_coverage_tsv, paths.repo_root
        ),
        "drug_class_profile_v1_value_composition_tsv": repo_relative(
            value_composition_tsv, paths.repo_root
        ),
        "drug_class_profile_v1_arm_freeze_candidate_audit_tsv": repo_relative(
            arm_freeze_tsv, paths.repo_root
        ),
        "drug_class_profile_v1_summary_tsv": repo_relative(summary_tsv, paths.repo_root),
        "run_log_json": repo_relative(run_log_path, paths.repo_root),
        "drug_name_normalization_v1_latest_json": repo_relative(
            paths.drug_norm_latest_pointer, paths.repo_root
        ),
        "patient_treatment_grouping_v1_latest_json": repo_relative(
            paths.grouping_latest_pointer, paths.repo_root
        ),
        "baseline_analysis_v1_latest_json": repo_relative(
            paths.baseline_analysis_latest_pointer, paths.repo_root
        ),
    }

    with paths.latest_pointer.open("w", encoding="utf-8") as f:
        json.dump(latest_pointer_data, f, indent=2)

    print(f"  audit run directory : {audit_run_dir}")
    print(f"  latest pointer      : {paths.latest_pointer}")
    print(f"  group summary rows  : {len(group_summary_rows)}")
    print(f"  confounder rows     : {len(confounder_coverage_rows)}")
    print(f"  composition rows    : {len(value_composition_rows)}")
    print(f"  candidate audit rows: {len(arm_freeze_candidate_audit_rows)}")
    print(f"  summary rows        : {len(summary_rows)}")
    print(f"  promising candidates: {promising_count}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    started_at = utc_now()
    run_id = started_at.astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    print(f"[24] drug-class profile review v1 — run_id={run_id}")

    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    print(f"  repo root: {paths.repo_root}")

    print("  loading inputs ...")
    inputs = load_inputs(paths, helper_module)
    print(
        f"  drug-class profile rows : {len(inputs.drug_class_profile_rows)}"
        f" | grouping rows: {len(inputs.grouping_rows)}"
        f" | baseline rows: {len(inputs.baseline_analysis_rows)}"
    )

    print("  merging patient rows ...")
    merged_rows = build_merged_rows(inputs)
    print(f"  merged rows: {len(merged_rows)}")

    print("  classifying buckets ...")
    bucket_map: dict[str, str] = {}
    for row in merged_rows:
        barcode = normalize_barcode(str(row.get("bcr_patient_barcode", "")))
        bucket_map[barcode] = classify_bucket(row)

    from collections import Counter as _Counter
    bucket_counts = _Counter(bucket_map.values())
    for bucket in APPROVED_BUCKET_ORDER:
        print(f"    {bucket}: {bucket_counts[bucket]}")
    total = len(merged_rows)
    assert sum(bucket_counts.values()) == total, (
        f"Bucket sum {sum(bucket_counts.values())} != total {total}"
    )

    print("  identifying major drug classes ...")
    major_classes = identify_major_classes(merged_rows)
    print(f"  major classes (N>={MIN_CLASS_SIZE_FOR_AUDIT}): {major_classes}")

    print("  building group summary ...")
    group_summary_rows = build_group_summary_rows(merged_rows, bucket_map, run_id)

    print("  building confounder coverage ...")
    confounder_coverage_rows = build_confounder_coverage_rows(
        merged_rows, bucket_map, major_classes, run_id
    )

    print("  building value composition ...")
    value_composition_rows = build_value_composition_rows(
        merged_rows, bucket_map, major_classes, run_id
    )

    print("  building arm-freeze candidate audit ...")
    arm_freeze_candidate_audit_rows = build_arm_freeze_candidate_audit_rows(
        merged_rows, bucket_map, major_classes, run_id
    )

    print("  building summary ...")
    summary_rows = build_summary_rows(
        merged_rows,
        bucket_map,
        major_classes,
        arm_freeze_candidate_audit_rows,
        inputs,
        run_id,
    )

    print("  writing outputs ...")
    write_outputs(
        paths=paths,
        run_id=run_id,
        group_summary_rows=group_summary_rows,
        confounder_coverage_rows=confounder_coverage_rows,
        value_composition_rows=value_composition_rows,
        arm_freeze_candidate_audit_rows=arm_freeze_candidate_audit_rows,
        summary_rows=summary_rows,
        inputs=inputs,
        started_at=started_at,
    )

    print(f"[24] completed — run_id={run_id}")


if __name__ == "__main__":
    main()
