#!/usr/bin/env python
"""Profile TCGA-BRCA baseline-analysis-prep v1 fields for descriptive review."""

from __future__ import annotations

import importlib.util
import json
import math
import statistics
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

FIELD_SUMMARY_FIELDNAMES = [
    "field_name",
    "source_origin",
    "field_category",
    "value_container_type",
    "row_count",
    "non_missing_count",
    "missing_like_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "dominant_value",
    "dominant_value_fraction",
    "example_values_small_sample",
    "profiling_interpretation",
]
MISSINGNESS_RANKED_FIELDNAMES = [
    "missingness_rank",
    "field_name",
    "source_origin",
    "field_category",
    "value_container_type",
    "row_count",
    "non_missing_count",
    "missing_like_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "dominant_value",
    "dominant_value_fraction",
    "profiling_interpretation",
]
VALUE_SUMMARY_FIELDNAMES = [
    "field_name",
    "source_origin",
    "field_category",
    "value_summary_type",
    "non_missing_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "numeric_min",
    "numeric_median",
    "numeric_max",
    "top_values_json",
    "notes",
]
CANDIDATE_FIELDS_FIELDNAMES = [
    "field_name",
    "source_origin",
    "field_category",
    "value_container_type",
    "non_missing_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "dominant_value_fraction",
    "shortlist_bucket",
    "shortlist_rule",
    "shortlist_manual_review_priority",
    "candidate_bucket",
    "candidate_rule",
    "excluded_reason",
    "reason",
    "needs_manual_review",
    "notes_placeholder",
]
EXCLUDED_FIELDS_FIELDNAMES = [
    "field_name",
    "source_origin",
    "field_category",
    "value_container_type",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "dominant_value_fraction",
    "candidate_rule",
    "excluded_reason",
    "reason",
    "needs_manual_review",
    "notes_placeholder",
]
SUMMARY_FIELDNAMES = [
    "baseline_profile_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
NOTES_PLACEHOLDER = "[fill in during baseline profile v1 review]"
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
BALANCED_MAX_MISSING_FRACTION = 0.20
REVIEW_MAX_MISSING_FRACTION = 0.50
BALANCED_MAX_DOMINANT_FRACTION = 0.95
MIN_DISTINCT_NON_MISSING_COUNT = 2
MAX_EXAMPLE_VALUES = 5
MAX_TOP_VALUE_COUNT = 10
FIELD_CATEGORY_ORDER = [
    "audit_id_field",
    "baseline_clinical_field",
    "followup_join_evidence_field",
    "endpoint_candidate_field",
    "biospecimen_sample_anchor_evidence_field",
    "derived_prep_flag",
]
CANDIDATE_BUCKET_ORDER = [
    "candidate_for_baseline_modeling_prep",
    "candidate_but_review_needed",
    "exclude_for_now",
]


class BaselineProfileV1Error(RuntimeError):
    """Raised when the baseline profile workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the baseline profile workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    baseline_latest_pointer: Path
    minimal_cohort_latest_pointer: Path
    clinical_shortlist_latest_pointer: Path
    ambiguity_resolution_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    baseline_latest_pointer: dict[str, Any]
    baseline_run_log: dict[str, Any]
    baseline_rows: list[dict[str, str]]
    baseline_spec_rows: list[dict[str, str]]
    baseline_summary_rows: list[dict[str, str]]
    minimal_cohort_latest_pointer: dict[str, Any]
    minimal_cohort_run_log: dict[str, Any]
    minimal_cohort_summary_rows: list[dict[str, str]]
    clinical_shortlist_latest_pointer: dict[str, Any]
    clinical_shortlist_run_log: dict[str, Any]
    clinical_shortlist_rows: list[dict[str, str]]
    ambiguity_resolution_latest_pointer: dict[str, Any]
    ambiguity_resolution_run_log: dict[str, Any]
    ambiguity_resolution_summary_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise BaselineProfileV1Error(f"Required helper script not found: {script_path}")

    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise BaselineProfileV1Error(f"Unable to create an import spec for: {script_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def parse_int(value: Any, label: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise BaselineProfileV1Error(f"Expected integer-like value for {label}: {value!r}") from exc


def parse_float(value: Any, label: str) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError) as exc:
        raise BaselineProfileV1Error(f"Expected float-like value for {label}: {value!r}") from exc


def ordered_unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def build_summary_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        metric = str(row.get("summary_metric") or "")
        if metric in lookup:
            raise BaselineProfileV1Error(f"Duplicate summary_metric detected: {metric}")
        lookup[metric] = row
    return lookup


def json_list(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=True)


def normalize_value(value: Any) -> str:
    return str(value or "").strip().lower()


def is_scalar_missing_like(value: Any) -> bool:
    return normalize_value(value) in MISSING_LIKE_TOKENS


def is_json_array_cell_missing_like(raw_value: Any) -> bool:
    raw_text = str(raw_value or "")
    if is_scalar_missing_like(raw_text):
        return True

    try:
        parsed = json.loads(raw_text)
    except json.JSONDecodeError:
        return is_scalar_missing_like(raw_text)

    if isinstance(parsed, list):
        if not parsed:
            return True
        return all(is_scalar_missing_like(item) for item in parsed)

    return is_scalar_missing_like(parsed)


def cell_is_missing_like(field_name: str, raw_value: Any) -> bool:
    if field_name.endswith("_json"):
        return is_json_array_cell_missing_like(raw_value)
    return is_scalar_missing_like(raw_value)


def format_fraction(value: float) -> str:
    return f"{value:.6f}"


def format_number(value: float) -> str:
    if not math.isfinite(value):
        raise BaselineProfileV1Error(f"Cannot format non-finite numeric value: {value!r}")
    if float(value).is_integer():
        return str(int(value))
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    if text in {"", "-0"}:
        return "0"
    return text


def metric_triplet_text(
    *,
    missing_like_fraction: float,
    distinct_non_missing_count: int,
    dominant_value_fraction: float,
) -> str:
    return (
        f"missing_like_fraction={missing_like_fraction:.3f}, "
        f"distinct_non_missing_count={distinct_non_missing_count}, "
        f"dominant_value_fraction={dominant_value_fraction:.3f}"
    )


def candidate_bucket_sort_key(bucket: str) -> tuple[int, str]:
    if bucket in CANDIDATE_BUCKET_ORDER:
        return (CANDIDATE_BUCKET_ORDER.index(bucket), bucket)
    return (len(CANDIDATE_BUCKET_ORDER), bucket)


def field_category_sort_key(field_category: str) -> tuple[int, str]:
    if field_category in FIELD_CATEGORY_ORDER:
        return (FIELD_CATEGORY_ORDER.index(field_category), field_category)
    return (len(FIELD_CATEGORY_ORDER), field_category)


def project_rows(rows: list[dict[str, str]], fieldnames: list[str]) -> list[dict[str, str]]:
    return [{fieldname: str(row.get(fieldname, "")) for fieldname in fieldnames} for row in rows]


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    helper_paths = helper_module.build_workflow_paths()
    trial_config_data = helper_module.load_yaml(helper_paths.trial_config)

    processed_root = helper_paths.repo_root / str(
        trial_config_data.get("processed_data_root", "01-data/processed")
    )
    audit_root = helper_paths.repo_root / str(
        trial_config_data.get("audit_root", "01-data/audit")
    )
    results_root = helper_paths.repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    analysis_prep_root = audit_root / "tcga-brca" / "analysis-prep"

    return WorkflowPaths(
        repo_root=helper_paths.repo_root,
        trial_config=helper_paths.trial_config,
        results_root=results_root,
        processed_runs_root=processed_root / "tcga-brca" / "analysis-prep" / "baseline_profile_v1_runs",
        audit_runs_root=analysis_prep_root / "baseline_profile_v1_runs",
        latest_pointer=analysis_prep_root / "tcga_brca_baseline_profile_v1_latest.json",
        baseline_latest_pointer=analysis_prep_root / "tcga_brca_baseline_analysis_v1_latest.json",
        minimal_cohort_latest_pointer=(
            audit_root / "tcga-brca" / "cohort" / "tcga_brca_minimal_cohort_v1_latest.json"
        ),
        clinical_shortlist_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_shortlist_latest.json"
        ),
        ambiguity_resolution_latest_pointer=(
            audit_root
            / "tcga-brca"
            / "cohort"
            / "tcga_brca_blueprint_ambiguity_resolution_latest.json"
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
        raise BaselineProfileV1Error(f"Required {label} not found: {pointer_path}")
    pointer = helper_module.load_json(pointer_path)
    helper_module.require_keys(pointer, required_keys, label, pointer_path)
    return pointer


def require_completed_run_log(
    repo_root: Path,
    run_log_relative_path: str,
    *,
    label: str,
    helper_module: Any,
) -> tuple[Path, dict[str, Any]]:
    run_log_path = helper_module.resolve_existing_path(repo_root, run_log_relative_path, label)
    run_log = helper_module.load_json(run_log_path)
    if run_log.get("status") != "completed":
        raise BaselineProfileV1Error(f"{label} is not completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise BaselineProfileV1Error(f"{label} does not report validation.passed == true.")
    return run_log_path, run_log


def build_input_paths(
    paths: WorkflowPaths,
    *,
    baseline_latest_pointer: dict[str, Any],
    minimal_cohort_latest_pointer: dict[str, Any],
    clinical_shortlist_latest_pointer: dict[str, Any],
    ambiguity_resolution_latest_pointer: dict[str, Any],
    helper_module: Any,
) -> dict[str, Path]:
    return {
        "baseline_analysis_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_latest_pointer["baseline_analysis_v1_tsv"]),
            "baseline analysis v1 TSV",
        ),
        "baseline_analysis_v1_spec_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_latest_pointer["baseline_analysis_v1_spec_tsv"]),
            "baseline analysis v1 spec TSV",
        ),
        "baseline_analysis_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_latest_pointer["baseline_analysis_v1_summary_tsv"]),
            "baseline analysis v1 summary TSV",
        ),
        "minimal_cohort_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["minimal_cohort_v1_summary_tsv"]),
            "minimal cohort v1 summary TSV",
        ),
        "clinical_shortlist_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(clinical_shortlist_latest_pointer["clinical_shortlist_tsv"]),
            "clinical shortlist TSV",
        ),
        "clinical_shortlist_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(clinical_shortlist_latest_pointer["clinical_shortlist_summary_tsv"]),
            "clinical shortlist summary TSV",
        ),
        "ambiguity_resolution_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["ambiguity_resolution_summary_tsv"]),
            "ambiguity-resolution summary TSV",
        ),
    }


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    baseline_latest_pointer = require_pointer(
        paths.baseline_latest_pointer,
        label="Baseline-analysis-prep latest pointer",
        required_keys={
            "baseline_analysis_v1_run_id",
            "cohort_v1_build_id",
            "dry_run_build_id",
            "blueprint_run_id",
            "ambiguity_resolution_run_id",
            "shortlist_run_id",
            "core_audit_run_id",
            "clinical_parse_run_id",
            "endpoint_crosswalk_run_id",
            "biospecimen_crosswalk_run_id",
            "biospecimen_parse_run_id",
            "clinical_source_run_id",
            "biospecimen_source_run_id",
            "baseline_analysis_v1_tsv",
            "baseline_analysis_v1_spec_tsv",
            "baseline_analysis_v1_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )
    minimal_cohort_latest_pointer = require_pointer(
        paths.minimal_cohort_latest_pointer,
        label="Minimal cohort v1 latest pointer",
        required_keys={
            "cohort_v1_build_id",
            "dry_run_build_id",
            "blueprint_run_id",
            "ambiguity_resolution_run_id",
            "shortlist_run_id",
            "core_audit_run_id",
            "clinical_parse_run_id",
            "endpoint_crosswalk_run_id",
            "biospecimen_crosswalk_run_id",
            "biospecimen_parse_run_id",
            "clinical_source_run_id",
            "biospecimen_source_run_id",
            "minimal_cohort_v1_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )
    clinical_shortlist_latest_pointer = require_pointer(
        paths.clinical_shortlist_latest_pointer,
        label="Clinical shortlist latest pointer",
        required_keys={
            "shortlist_run_id",
            "core_audit_run_id",
            "parse_run_id",
            "source_run_id",
            "clinical_shortlist_tsv",
            "clinical_shortlist_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )
    ambiguity_resolution_latest_pointer = require_pointer(
        paths.ambiguity_resolution_latest_pointer,
        label="Ambiguity-resolution latest pointer",
        required_keys={
            "ambiguity_resolution_run_id",
            "blueprint_run_id",
            "shortlist_run_id",
            "core_audit_run_id",
            "clinical_parse_run_id",
            "endpoint_crosswalk_run_id",
            "biospecimen_crosswalk_run_id",
            "biospecimen_parse_run_id",
            "clinical_source_run_id",
            "biospecimen_source_run_id",
            "ambiguity_resolution_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )

    expected_pairings = {
        "cohort_v1_build_id": str(minimal_cohort_latest_pointer["cohort_v1_build_id"]),
        "dry_run_build_id": str(minimal_cohort_latest_pointer["dry_run_build_id"]),
        "blueprint_run_id": str(minimal_cohort_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(clinical_shortlist_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(clinical_shortlist_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(clinical_shortlist_latest_pointer["parse_run_id"]),
        "clinical_source_run_id": str(clinical_shortlist_latest_pointer["source_run_id"]),
        "endpoint_crosswalk_run_id": str(minimal_cohort_latest_pointer["endpoint_crosswalk_run_id"]),
        "biospecimen_crosswalk_run_id": str(
            minimal_cohort_latest_pointer["biospecimen_crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(minimal_cohort_latest_pointer["biospecimen_parse_run_id"]),
        "biospecimen_source_run_id": str(minimal_cohort_latest_pointer["biospecimen_source_run_id"]),
    }
    for key, expected_value in expected_pairings.items():
        observed_value = str(baseline_latest_pointer[key])
        if observed_value != expected_value:
            raise BaselineProfileV1Error(
                "Baseline-analysis-prep latest pointer is out of sync with current upstream latest pointers: "
                f"{key}={observed_value!r} vs expected {expected_value!r}."
            )

    if str(minimal_cohort_latest_pointer["shortlist_run_id"]) != str(
        clinical_shortlist_latest_pointer["shortlist_run_id"]
    ):
        raise BaselineProfileV1Error(
            "Minimal cohort v1 latest pointer is out of sync with the current clinical shortlist latest pointer."
        )
    if str(minimal_cohort_latest_pointer["ambiguity_resolution_run_id"]) != str(
        ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
    ):
        raise BaselineProfileV1Error(
            "Minimal cohort v1 latest pointer is out of sync with the current ambiguity-resolution latest pointer."
        )

    baseline_run_log_path, baseline_run_log = require_completed_run_log(
        paths.repo_root,
        str(baseline_latest_pointer["run_log_json"]),
        label="Baseline-analysis-prep run log",
        helper_module=helper_module,
    )
    minimal_cohort_run_log_path, minimal_cohort_run_log = require_completed_run_log(
        paths.repo_root,
        str(minimal_cohort_latest_pointer["run_log_json"]),
        label="Minimal cohort v1 run log",
        helper_module=helper_module,
    )
    clinical_shortlist_run_log_path, clinical_shortlist_run_log = require_completed_run_log(
        paths.repo_root,
        str(clinical_shortlist_latest_pointer["run_log_json"]),
        label="Clinical shortlist run log",
        helper_module=helper_module,
    )
    ambiguity_resolution_run_log_path, ambiguity_resolution_run_log = require_completed_run_log(
        paths.repo_root,
        str(ambiguity_resolution_latest_pointer["run_log_json"]),
        label="Ambiguity-resolution run log",
        helper_module=helper_module,
    )

    input_paths = build_input_paths(
        paths,
        baseline_latest_pointer=baseline_latest_pointer,
        minimal_cohort_latest_pointer=minimal_cohort_latest_pointer,
        clinical_shortlist_latest_pointer=clinical_shortlist_latest_pointer,
        ambiguity_resolution_latest_pointer=ambiguity_resolution_latest_pointer,
        helper_module=helper_module,
    )
    input_paths["baseline_analysis_v1_run_log_json"] = baseline_run_log_path
    input_paths["minimal_cohort_v1_run_log_json"] = minimal_cohort_run_log_path
    input_paths["clinical_shortlist_run_log_json"] = clinical_shortlist_run_log_path
    input_paths["ambiguity_resolution_run_log_json"] = ambiguity_resolution_run_log_path

    baseline_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_analysis_v1_tsv"])
    baseline_spec_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_analysis_v1_spec_tsv"])
    baseline_summary_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_analysis_v1_summary_tsv"])
    minimal_cohort_summary_rows = helper_module.read_tsv_dict_rows(input_paths["minimal_cohort_v1_summary_tsv"])
    clinical_shortlist_rows = helper_module.read_tsv_dict_rows(input_paths["clinical_shortlist_tsv"])
    ambiguity_resolution_summary_rows = helper_module.read_tsv_dict_rows(
        input_paths["ambiguity_resolution_summary_tsv"]
    )

    if not baseline_rows:
        raise BaselineProfileV1Error("The referenced baseline_analysis_v1.tsv contains no rows.")
    if not baseline_spec_rows:
        raise BaselineProfileV1Error("The referenced baseline_analysis_v1_spec.tsv contains no rows.")
    if not clinical_shortlist_rows:
        raise BaselineProfileV1Error("The referenced clinical_shortlist.tsv contains no rows.")

    baseline_summary_lookup = build_summary_lookup(baseline_summary_rows)
    if "final_row_count" not in baseline_summary_lookup:
        raise BaselineProfileV1Error(
            "Baseline analysis v1 summary TSV is missing summary_metric=final_row_count."
        )
    baseline_summary_row_count = parse_int(
        baseline_summary_lookup["final_row_count"]["summary_value"],
        "baseline analysis v1 final_row_count",
    )
    if baseline_summary_row_count != len(baseline_rows):
        raise BaselineProfileV1Error(
            "Baseline analysis v1 summary row count does not match baseline_analysis_v1.tsv row count: "
            f"{baseline_summary_row_count} vs {len(baseline_rows)}."
        )

    minimal_cohort_summary_lookup = build_summary_lookup(minimal_cohort_summary_rows)
    if "final_v1_row_count" not in minimal_cohort_summary_lookup:
        raise BaselineProfileV1Error(
            "Minimal cohort v1 summary TSV is missing summary_metric=final_v1_row_count."
        )
    minimal_cohort_summary_row_count = parse_int(
        minimal_cohort_summary_lookup["final_v1_row_count"]["summary_value"],
        "minimal cohort v1 final_v1_row_count",
    )
    if minimal_cohort_summary_row_count != len(baseline_rows):
        raise BaselineProfileV1Error(
            "Minimal cohort v1 summary row count does not match baseline_analysis_v1.tsv row count: "
            f"{minimal_cohort_summary_row_count} vs {len(baseline_rows)}."
        )

    retained_spec_field_names = ordered_unique(
        [
            str(row.get("field_name") or "")
            for row in baseline_spec_rows
            if str(row.get("retained_in_baseline_analysis_v1") or "") == "yes"
        ]
    )
    baseline_field_names = list(baseline_rows[0].keys())
    if len(retained_spec_field_names) != len(baseline_field_names) or set(retained_spec_field_names) != set(
        baseline_field_names
    ):
        raise BaselineProfileV1Error(
            "Retained baseline-analysis-prep spec rows do not reconcile to baseline_analysis_v1.tsv columns."
        )

    return WorkflowInputs(
        baseline_latest_pointer=baseline_latest_pointer,
        baseline_run_log=baseline_run_log,
        baseline_rows=baseline_rows,
        baseline_spec_rows=baseline_spec_rows,
        baseline_summary_rows=baseline_summary_rows,
        minimal_cohort_latest_pointer=minimal_cohort_latest_pointer,
        minimal_cohort_run_log=minimal_cohort_run_log,
        minimal_cohort_summary_rows=minimal_cohort_summary_rows,
        clinical_shortlist_latest_pointer=clinical_shortlist_latest_pointer,
        clinical_shortlist_run_log=clinical_shortlist_run_log,
        clinical_shortlist_rows=clinical_shortlist_rows,
        ambiguity_resolution_latest_pointer=ambiguity_resolution_latest_pointer,
        ambiguity_resolution_run_log=ambiguity_resolution_run_log,
        ambiguity_resolution_summary_rows=ambiguity_resolution_summary_rows,
        input_paths=input_paths,
    )


def build_retained_spec_lookup(
    baseline_spec_rows: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in baseline_spec_rows:
        if str(row.get("retained_in_baseline_analysis_v1") or "") != "yes":
            continue
        field_name = str(row.get("field_name") or "")
        if not field_name:
            raise BaselineProfileV1Error("Encountered an empty field_name in retained baseline spec rows.")
        if field_name in lookup:
            raise BaselineProfileV1Error(f"Duplicate retained baseline spec row detected: {field_name}")
        lookup[field_name] = row
    return lookup


def build_shortlist_lookup(
    shortlist_rows: list[dict[str, str]],
) -> dict[tuple[str, str], dict[str, str]]:
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for row in shortlist_rows:
        table_name = str(row.get("table_name") or "")
        field_name = str(row.get("field_name") or "")
        key = (table_name, field_name)
        if key in lookup:
            raise BaselineProfileV1Error(
                f"Duplicate clinical shortlist row detected for table_name={table_name!r}, field_name={field_name!r}."
            )
        lookup[key] = row
    return lookup


def value_container_type(field_name: str) -> str:
    return "json_array_cell" if field_name.endswith("_json") else "scalar"


def example_values_small_sample(non_missing_values: list[str]) -> str:
    return json_list(ordered_unique(non_missing_values)[:MAX_EXAMPLE_VALUES])


def dominant_value_stats(non_missing_values: list[str]) -> tuple[str, float]:
    if not non_missing_values:
        return ("", 0.0)

    counts: dict[str, int] = {}
    for value in non_missing_values:
        counts[value] = counts.get(value, 0) + 1

    dominant_value, dominant_count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
    return (dominant_value, dominant_count / len(non_missing_values))


def build_field_profiles(
    workflow_inputs: WorkflowInputs,
) -> list[dict[str, str]]:
    retained_spec_lookup = build_retained_spec_lookup(workflow_inputs.baseline_spec_rows)
    shortlist_lookup = build_shortlist_lookup(workflow_inputs.clinical_shortlist_rows)
    row_count = len(workflow_inputs.baseline_rows)
    field_profiles: list[dict[str, str]] = []

    for field_name in workflow_inputs.baseline_rows[0].keys():
        if field_name not in retained_spec_lookup:
            raise BaselineProfileV1Error(
                f"baseline_analysis_v1.tsv column is missing from the retained baseline spec lookup: {field_name}"
            )

        spec_row = retained_spec_lookup[field_name]
        source_origin = str(spec_row.get("source_origin") or "")
        field_category = str(spec_row.get("field_category") or "")
        raw_values = [str(row.get(field_name, "")) for row in workflow_inputs.baseline_rows]
        non_missing_values: list[str] = []
        missing_like_count = 0
        for raw_value in raw_values:
            if cell_is_missing_like(field_name, raw_value):
                missing_like_count += 1
                continue
            non_missing_values.append(raw_value)

        distinct_non_missing_count = len(set(non_missing_values))
        dominant_value, dominant_value_fraction = dominant_value_stats(non_missing_values)
        missing_like_fraction = (missing_like_count / row_count) if row_count else 0.0
        shortlist_row = shortlist_lookup.get((source_origin, field_name), {})

        profile_row = {
            "field_name": field_name,
            "source_origin": source_origin,
            "field_category": field_category,
            "value_container_type": value_container_type(field_name),
            "row_count": str(row_count),
            "non_missing_count": str(len(non_missing_values)),
            "missing_like_count": str(missing_like_count),
            "missing_like_fraction": format_fraction(missing_like_fraction),
            "distinct_non_missing_count": str(distinct_non_missing_count),
            "dominant_value": dominant_value,
            "dominant_value_fraction": format_fraction(dominant_value_fraction),
            "example_values_small_sample": example_values_small_sample(non_missing_values),
            "profiling_interpretation": "",
            "shortlist_bucket": str(shortlist_row.get("shortlist_bucket") or ""),
            "shortlist_rule": str(shortlist_row.get("shortlist_rule") or ""),
            "shortlist_manual_review_priority": str(shortlist_row.get("manual_review_priority") or ""),
            "ambiguity_carried_forward": str(spec_row.get("ambiguity_carried_forward") or ""),
        }
        field_profiles.append(profile_row)

    return field_profiles


def classify_candidate_row(profile_row: dict[str, str]) -> dict[str, str]:
    field_name = profile_row["field_name"]
    field_category = profile_row["field_category"]
    missing_like_fraction = parse_float(profile_row["missing_like_fraction"], f"{field_name} missing_like_fraction")
    distinct_non_missing_count = parse_int(
        profile_row["distinct_non_missing_count"], f"{field_name} distinct_non_missing_count"
    )
    dominant_value_fraction = parse_float(
        profile_row["dominant_value_fraction"], f"{field_name} dominant_value_fraction"
    )
    metrics_text = metric_triplet_text(
        missing_like_fraction=missing_like_fraction,
        distinct_non_missing_count=distinct_non_missing_count,
        dominant_value_fraction=dominant_value_fraction,
    )

    if field_category == "audit_id_field":
        excluded_reason = (
            "ambiguity_carried_forward" if field_name.endswith("_flags_json") else "identifier_only"
        )
        return {
            "candidate_bucket": "exclude_for_now",
            "candidate_rule": "audit_or_identifier_exclude",
            "excluded_reason": excluded_reason,
            "reason": (
                "Audit and identifier tracking fields are preserved for provenance and validation only; "
                f"they are not direct first-pass baseline modeling inputs ({metrics_text})."
            ),
            "needs_manual_review": "no",
        }

    if field_category == "derived_prep_flag":
        return {
            "candidate_bucket": "exclude_for_now",
            "candidate_rule": "derived_flag_exclude",
            "excluded_reason": "derived_flag_only",
            "reason": (
                "Derived prep flags are useful for audit coverage and review only; "
                f"they are not direct first-pass baseline modeling inputs ({metrics_text})."
            ),
            "needs_manual_review": "no",
        }

    if field_category == "followup_join_evidence_field":
        return {
            "candidate_bucket": "candidate_but_review_needed",
            "candidate_rule": "followup_join_evidence_review_only",
            "excluded_reason": "",
            "reason": (
                "Grouped follow-up join evidence is retained for audit traceability and linkage review only; "
                f"it should stay out of first-pass baseline modeling prep ({metrics_text})."
            ),
            "needs_manual_review": "yes",
        }

    if field_category == "endpoint_candidate_field":
        return {
            "candidate_bucket": "candidate_but_review_needed",
            "candidate_rule": "endpoint_candidate_not_frozen_review",
            "excluded_reason": "",
            "reason": (
                "Endpoint candidate fields remain provisional because the endpoint freeze is still blocked; "
                f"keep them out of first-pass baseline modeling prep ({metrics_text})."
            ),
            "needs_manual_review": "yes",
        }

    if field_category == "biospecimen_sample_anchor_evidence_field":
        return {
            "candidate_bucket": "candidate_but_review_needed",
            "candidate_rule": "biospecimen_evidence_container_review_only",
            "excluded_reason": "",
            "reason": (
                "Biospecimen sample-anchor columns remain evidence containers for later review rather than direct "
                f"first-pass baseline inputs ({metrics_text})."
            ),
            "needs_manual_review": "yes",
        }

    if field_category != "baseline_clinical_field":
        raise BaselineProfileV1Error(f"Unexpected field_category for candidate classification: {field_category}")

    if distinct_non_missing_count < MIN_DISTINCT_NON_MISSING_COUNT:
        return {
            "candidate_bucket": "exclude_for_now",
            "candidate_rule": "constant_or_low_information_exclude",
            "excluded_reason": "near_constant_low_information",
            "reason": (
                "The retained baseline clinical field does not show enough non-missing variation to support "
                f"first-pass baseline modeling prep ({metrics_text})."
            ),
            "needs_manual_review": "no",
        }

    if missing_like_fraction > REVIEW_MAX_MISSING_FRACTION:
        return {
            "candidate_bucket": "exclude_for_now",
            "candidate_rule": "high_missingness_exclude",
            "excluded_reason": "too_sparse",
            "reason": (
                "The retained baseline clinical field is too sparse for first-pass baseline modeling prep under the "
                f"balanced profile ({metrics_text})."
            ),
            "needs_manual_review": "no",
        }

    if (
        missing_like_fraction <= BALANCED_MAX_MISSING_FRACTION
        and dominant_value_fraction <= BALANCED_MAX_DOMINANT_FRACTION
    ):
        return {
            "candidate_bucket": "candidate_for_baseline_modeling_prep",
            "candidate_rule": "baseline_balanced_complete_variant",
            "excluded_reason": "",
            "reason": (
                "The retained baseline clinical field meets the balanced completeness and variation thresholds for "
                f"first-pass baseline modeling prep ({metrics_text})."
            ),
            "needs_manual_review": "no",
        }

    if dominant_value_fraction > BALANCED_MAX_DOMINANT_FRACTION:
        return {
            "candidate_bucket": "candidate_but_review_needed",
            "candidate_rule": "baseline_balanced_near_constant_review",
            "excluded_reason": "",
            "reason": (
                "The retained baseline clinical field still varies, but its dominant category is unusually strong and "
                f"should be reviewed before first-pass baseline modeling prep ({metrics_text})."
            ),
            "needs_manual_review": "yes",
        }

    return {
        "candidate_bucket": "candidate_but_review_needed",
        "candidate_rule": "baseline_balanced_moderate_missingness_review",
        "excluded_reason": "",
        "reason": (
            "The retained baseline clinical field still varies, but its completeness falls outside the direct-candidate "
            f"balanced threshold and should be reviewed before first-pass baseline modeling prep ({metrics_text})."
        ),
        "needs_manual_review": "yes",
    }


def profiling_interpretation_for_candidate(
    candidate_row: dict[str, str],
    field_category: str,
) -> str:
    candidate_bucket = candidate_row["candidate_bucket"]
    candidate_rule = candidate_row["candidate_rule"]

    if field_category == "audit_id_field":
        return "audit_or_identifier_tracking_field"
    if field_category == "derived_prep_flag":
        return "derived_review_flag"
    if field_category == "followup_join_evidence_field":
        return "followup_join_evidence_container"
    if field_category == "endpoint_candidate_field":
        return "provisional_endpoint_candidate"
    if field_category == "biospecimen_sample_anchor_evidence_field":
        return "biospecimen_sample_anchor_evidence_container"
    if candidate_bucket == "candidate_for_baseline_modeling_prep":
        return "baseline_clinical_balanced_candidate"
    if candidate_rule == "baseline_balanced_near_constant_review":
        return "baseline_clinical_review_near_constant"
    if candidate_rule == "baseline_balanced_moderate_missingness_review":
        return "baseline_clinical_review_moderate_missingness"
    if candidate_rule == "high_missingness_exclude":
        return "baseline_clinical_exclude_high_missingness"
    if candidate_rule == "constant_or_low_information_exclude":
        return "baseline_clinical_exclude_low_information"
    raise BaselineProfileV1Error(
        f"Unable to derive profiling_interpretation for field_category={field_category!r}, "
        f"candidate_rule={candidate_rule!r}."
    )


def build_candidate_rows(field_profiles: list[dict[str, str]]) -> list[dict[str, str]]:
    candidate_rows: list[dict[str, str]] = []
    for profile_row in field_profiles:
        classification = classify_candidate_row(profile_row)
        profile_row["profiling_interpretation"] = profiling_interpretation_for_candidate(
            classification,
            profile_row["field_category"],
        )
        candidate_rows.append(
            {
                "field_name": profile_row["field_name"],
                "source_origin": profile_row["source_origin"],
                "field_category": profile_row["field_category"],
                "value_container_type": profile_row["value_container_type"],
                "non_missing_count": profile_row["non_missing_count"],
                "missing_like_fraction": profile_row["missing_like_fraction"],
                "distinct_non_missing_count": profile_row["distinct_non_missing_count"],
                "dominant_value_fraction": profile_row["dominant_value_fraction"],
                "shortlist_bucket": profile_row["shortlist_bucket"],
                "shortlist_rule": profile_row["shortlist_rule"],
                "shortlist_manual_review_priority": profile_row["shortlist_manual_review_priority"],
                "candidate_bucket": classification["candidate_bucket"],
                "candidate_rule": classification["candidate_rule"],
                "excluded_reason": classification["excluded_reason"],
                "reason": classification["reason"],
                "needs_manual_review": classification["needs_manual_review"],
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )

    candidate_rows.sort(
        key=lambda row: (
            candidate_bucket_sort_key(str(row["candidate_bucket"])),
            field_category_sort_key(str(row["field_category"])),
            str(row["field_name"]),
        )
    )
    return candidate_rows


def build_missingness_ranked_rows(field_profiles: list[dict[str, str]]) -> list[dict[str, str]]:
    sorted_profiles = sorted(
        field_profiles,
        key=lambda row: (
            -parse_float(row["missing_like_fraction"], f"{row['field_name']} missing_like_fraction"),
            str(row["field_name"]),
        ),
    )
    missingness_rows: list[dict[str, str]] = []
    for rank, profile_row in enumerate(sorted_profiles, start=1):
        missingness_rows.append(
            {
                "missingness_rank": str(rank),
                "field_name": profile_row["field_name"],
                "source_origin": profile_row["source_origin"],
                "field_category": profile_row["field_category"],
                "value_container_type": profile_row["value_container_type"],
                "row_count": profile_row["row_count"],
                "non_missing_count": profile_row["non_missing_count"],
                "missing_like_count": profile_row["missing_like_count"],
                "missing_like_fraction": profile_row["missing_like_fraction"],
                "distinct_non_missing_count": profile_row["distinct_non_missing_count"],
                "dominant_value": profile_row["dominant_value"],
                "dominant_value_fraction": profile_row["dominant_value_fraction"],
                "profiling_interpretation": profile_row["profiling_interpretation"],
            }
        )
    return missingness_rows


def try_parse_numeric_values(non_missing_values: list[str]) -> list[float] | None:
    numeric_values: list[float] = []
    for value in non_missing_values:
        try:
            numeric_values.append(float(value))
        except ValueError:
            return None
    return numeric_values


def top_values_json(non_missing_values: list[str]) -> str:
    counts: dict[str, int] = {}
    for value in non_missing_values:
        counts[value] = counts.get(value, 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:MAX_TOP_VALUE_COUNT]
    payload = [
        {
            "value": value,
            "count": count,
            "fraction": format_fraction(count / len(non_missing_values)),
        }
        for value, count in ordered
    ]
    return json.dumps(payload, ensure_ascii=True)


def build_value_summary_rows(
    field_profiles: list[dict[str, str]],
    baseline_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    row_count = len(baseline_rows)

    for profile_row in field_profiles:
        if profile_row["field_category"] != "baseline_clinical_field":
            continue
        if profile_row["value_container_type"] != "scalar":
            continue

        field_name = profile_row["field_name"]
        non_missing_values = [
            str(row.get(field_name, ""))
            for row in baseline_rows
            if not cell_is_missing_like(field_name, row.get(field_name, ""))
        ]
        numeric_values = try_parse_numeric_values(non_missing_values)
        if numeric_values is not None and non_missing_values:
            numeric_min = min(numeric_values)
            numeric_median = statistics.median(numeric_values)
            numeric_max = max(numeric_values)
            rows.append(
                {
                    "field_name": field_name,
                    "source_origin": profile_row["source_origin"],
                    "field_category": profile_row["field_category"],
                    "value_summary_type": "numeric_like_scalar",
                    "non_missing_count": profile_row["non_missing_count"],
                    "missing_like_fraction": profile_row["missing_like_fraction"],
                    "distinct_non_missing_count": profile_row["distinct_non_missing_count"],
                    "numeric_min": format_number(float(numeric_min)),
                    "numeric_median": format_number(float(numeric_median)),
                    "numeric_max": format_number(float(numeric_max)),
                    "top_values_json": "",
                    "notes": "All non-missing retained scalar values parsed safely as numeric.",
                }
            )
            continue

        rows.append(
            {
                "field_name": field_name,
                "source_origin": profile_row["source_origin"],
                "field_category": profile_row["field_category"],
                "value_summary_type": "categorical_like_scalar",
                "non_missing_count": profile_row["non_missing_count"],
                "missing_like_fraction": profile_row["missing_like_fraction"],
                "distinct_non_missing_count": profile_row["distinct_non_missing_count"],
                "numeric_min": "",
                "numeric_median": "",
                "numeric_max": "",
                "top_values_json": top_values_json(non_missing_values) if row_count else "[]",
                "notes": "Retained scalar field summarized with top categorical values only.",
            }
        )

    rows.sort(key=lambda row: str(row["field_name"]))
    return rows


def build_excluded_rows(candidate_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    excluded_rows = [
        {
            "field_name": row["field_name"],
            "source_origin": row["source_origin"],
            "field_category": row["field_category"],
            "value_container_type": row["value_container_type"],
            "missing_like_fraction": row["missing_like_fraction"],
            "distinct_non_missing_count": row["distinct_non_missing_count"],
            "dominant_value_fraction": row["dominant_value_fraction"],
            "candidate_rule": row["candidate_rule"],
            "excluded_reason": row["excluded_reason"],
            "reason": row["reason"],
            "needs_manual_review": row["needs_manual_review"],
            "notes_placeholder": row["notes_placeholder"],
        }
        for row in candidate_rows
        if row["candidate_bucket"] == "exclude_for_now"
    ]
    excluded_rows.sort(key=lambda row: (str(row["excluded_reason"]), str(row["field_name"])))
    return excluded_rows


def count_yes(baseline_rows: list[dict[str, str]], field_name: str) -> int:
    return sum(str(row.get(field_name) or "") == "yes" for row in baseline_rows)


def build_summary_rows(
    baseline_profile_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
    field_profiles: list[dict[str, str]],
    candidate_rows: list[dict[str, str]],
    value_summary_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    category_counts: dict[str, int] = {}
    for row in field_profiles:
        field_category = str(row["field_category"])
        category_counts[field_category] = category_counts.get(field_category, 0) + 1

    candidate_count = sum(
        row["candidate_bucket"] == "candidate_for_baseline_modeling_prep" for row in candidate_rows
    )
    review_needed_count = sum(
        row["candidate_bucket"] == "candidate_but_review_needed" for row in candidate_rows
    )
    excluded_count = sum(row["candidate_bucket"] == "exclude_for_now" for row in candidate_rows)
    high_missingness_count = sum(
        parse_float(row["missing_like_fraction"], f"{row['field_name']} missing_like_fraction")
        > REVIEW_MAX_MISSING_FRACTION
        for row in field_profiles
    )
    near_constant_count = sum(
        parse_float(row["dominant_value_fraction"], f"{row['field_name']} dominant_value_fraction")
        > BALANCED_MAX_DOMINANT_FRACTION
        for row in field_profiles
    )
    rows_with_followup_evidence = count_yes(workflow_inputs.baseline_rows, "row_has_followup_match")
    rows_with_biospecimen_sample_evidence = count_yes(
        workflow_inputs.baseline_rows, "row_has_biospecimen_sample_match"
    )
    final_row_count = len(workflow_inputs.baseline_rows)

    return [
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "design",
            "summary_metric": "baseline_profile_v1_status",
            "summary_value": "descriptive_baseline_profile_v1",
            "notes": (
                "This workflow is a descriptive baseline profiling and candidate-field selection layer only."
            ),
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "design",
            "summary_metric": "profiling_source_table",
            "summary_value": "baseline_analysis_v1",
            "notes": "The profile starts only from the saved baseline_analysis_v1 TSV on disk.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "design",
            "summary_metric": "candidate_selection_profile",
            "summary_value": "balanced",
            "notes": (
                "Direct candidate threshold: missing<=0.20, dominant<=0.95, distinct>=2; "
                "review threshold extends to missing<=0.50 for varying baseline clinical fields."
            ),
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "design",
            "summary_metric": "proposed_unit_of_analysis",
            "summary_value": "patient/case",
            "notes": "The profile preserves one row per patient/case from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "design",
            "summary_metric": "endpoint_policy",
            "summary_value": "review_only_not_frozen",
            "notes": "Endpoint candidates remain review-only and are not collapsed into a final endpoint.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "design",
            "summary_metric": "treatment_policy",
            "summary_value": "excluded",
            "notes": "Treatment detail remains excluded from this profiling layer.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "baseline_analysis_v1_run_id",
            "summary_value": str(workflow_inputs.baseline_latest_pointer["baseline_analysis_v1_run_id"]),
            "notes": "Baseline-analysis-prep run profiled by this workflow.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "cohort_v1_build_id",
            "summary_value": str(workflow_inputs.baseline_latest_pointer["cohort_v1_build_id"]),
            "notes": "Source minimal cohort v1 build carried forward into baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "ambiguity_resolution_run_id",
            "summary_value": str(workflow_inputs.baseline_latest_pointer["ambiguity_resolution_run_id"]),
            "notes": "Current ambiguity-resolution run aligned to the profile inputs.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "shortlist_run_id",
            "summary_value": str(workflow_inputs.baseline_latest_pointer["shortlist_run_id"]),
            "notes": "Current clinical shortlist run aligned to the profile inputs.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "input_baseline_analysis_v1_row_count",
            "summary_value": str(final_row_count),
            "notes": "Row count read directly from baseline_analysis_v1.tsv.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "final_row_count",
            "summary_value": str(final_row_count),
            "notes": "Final row count represented in the baseline profile outputs.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_field_count_total",
            "summary_value": str(len(field_profiles)),
            "notes": "Total retained fields profiled from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_audit_id_field_count",
            "summary_value": str(category_counts.get("audit_id_field", 0)),
            "notes": "Retained audit/id fields carried from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_baseline_clinical_field_count",
            "summary_value": str(category_counts.get("baseline_clinical_field", 0)),
            "notes": "Retained baseline clinical fields carried from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_followup_join_evidence_field_count",
            "summary_value": str(category_counts.get("followup_join_evidence_field", 0)),
            "notes": "Retained follow-up join-evidence fields carried from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_endpoint_candidate_field_count",
            "summary_value": str(category_counts.get("endpoint_candidate_field", 0)),
            "notes": "Retained endpoint candidate fields carried side by side from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_biospecimen_sample_anchor_evidence_field_count",
            "summary_value": str(category_counts.get("biospecimen_sample_anchor_evidence_field", 0)),
            "notes": "Retained biospecimen sample-anchor evidence fields carried from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_derived_prep_flag_count",
            "summary_value": str(category_counts.get("derived_prep_flag", 0)),
            "notes": "Retained derived prep flags carried from baseline_analysis_v1.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "candidate_selection",
            "summary_metric": "candidate_field_count",
            "summary_value": str(candidate_count),
            "notes": "Count of retained fields suggested as direct first-pass baseline candidates.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "candidate_selection",
            "summary_metric": "review_needed_field_count",
            "summary_value": str(review_needed_count),
            "notes": "Count of retained fields kept for manual review before first-pass baseline prep.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "candidate_selection",
            "summary_metric": "excluded_for_now_field_count",
            "summary_value": str(excluded_count),
            "notes": "Count of retained fields excluded for now from first-pass baseline prep.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "candidate_selection",
            "summary_metric": "high_missingness_field_count",
            "summary_value": str(high_missingness_count),
            "notes": "Count of retained fields with missing_like_fraction > 0.50.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "candidate_selection",
            "summary_metric": "near_constant_field_count",
            "summary_value": str(near_constant_count),
            "notes": "Count of retained fields with dominant_value_fraction > 0.95.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "candidate_selection",
            "summary_metric": "value_summary_field_count",
            "summary_value": str(len(value_summary_rows)),
            "notes": "Count of scalar baseline clinical fields summarized in the value-summary output.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "coverage",
            "summary_metric": "followup_evidence_coverage_count",
            "summary_value": str(rows_with_followup_evidence),
            "notes": "Rows where row_has_followup_match == yes.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "coverage",
            "summary_metric": "biospecimen_sample_evidence_coverage_count",
            "summary_value": str(rows_with_biospecimen_sample_evidence),
            "notes": "Rows where row_has_biospecimen_sample_match == yes.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "readiness_for_descriptive_reporting",
            "summary_value": "ready",
            "notes": "The saved field summaries and missingness outputs are suitable for descriptive review.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "readiness_for_baseline_modeling_prep",
            "summary_value": "ready_for_first_pass_feature_selection_only",
            "notes": (
                "The saved outputs are suitable for first-pass baseline feature selection only while endpoint and "
                "treatment constraints remain explicit."
            ),
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "endpoint_freeze_status",
            "summary_value": "blocked",
            "notes": "Endpoint freeze remains blocked; endpoint candidates stay out of first-pass baseline prep.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "treatment_inclusion_status",
            "summary_value": "excluded",
            "notes": "Treatment detail remains excluded from this profiling workflow.",
        },
        {
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "overall_readiness_interpretation",
            "summary_value": "suitable_for_descriptive_reporting_and_first_pass_baseline_feature_selection_only",
            "notes": (
                "Use this layer for descriptive baseline reporting, missingness review, and candidate-field selection "
                "only; do not treat it as model training, final feature freeze, endpoint freeze, or treatment modeling."
            ),
        },
    ]


def build_latest_pointer_payload(
    *,
    baseline_profile_v1_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
        "baseline_analysis_v1_run_id": str(workflow_inputs.baseline_latest_pointer["baseline_analysis_v1_run_id"]),
        "cohort_v1_build_id": str(workflow_inputs.baseline_latest_pointer["cohort_v1_build_id"]),
        "dry_run_build_id": str(workflow_inputs.baseline_latest_pointer["dry_run_build_id"]),
        "blueprint_run_id": str(workflow_inputs.baseline_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            workflow_inputs.baseline_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(workflow_inputs.baseline_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(workflow_inputs.baseline_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(workflow_inputs.baseline_latest_pointer["clinical_parse_run_id"]),
        "endpoint_crosswalk_run_id": str(
            workflow_inputs.baseline_latest_pointer["endpoint_crosswalk_run_id"]
        ),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.baseline_latest_pointer["biospecimen_crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(workflow_inputs.baseline_latest_pointer["biospecimen_parse_run_id"]),
        "clinical_source_run_id": str(workflow_inputs.baseline_latest_pointer["clinical_source_run_id"]),
        "biospecimen_source_run_id": str(workflow_inputs.baseline_latest_pointer["biospecimen_source_run_id"]),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "baseline_profile_v1_field_summary_tsv": repo_relative(
            output_paths["field_summary_tsv"], paths.repo_root
        ),
        "baseline_profile_v1_missingness_ranked_tsv": repo_relative(
            output_paths["missingness_ranked_tsv"], paths.repo_root
        ),
        "baseline_profile_v1_value_summary_tsv": repo_relative(
            output_paths["value_summary_tsv"], paths.repo_root
        ),
        "baseline_profile_v1_candidate_fields_tsv": repo_relative(
            output_paths["candidate_fields_tsv"], paths.repo_root
        ),
        "baseline_profile_v1_excluded_fields_tsv": repo_relative(
            output_paths["excluded_fields_tsv"], paths.repo_root
        ),
        "baseline_profile_v1_summary_tsv": repo_relative(output_paths["summary_tsv"], paths.repo_root),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "baseline_analysis_v1_latest_json": repo_relative(paths.baseline_latest_pointer, paths.repo_root),
        "minimal_cohort_v1_latest_json": repo_relative(paths.minimal_cohort_latest_pointer, paths.repo_root),
        "blueprint_ambiguity_resolution_latest_json": repo_relative(
            paths.ambiguity_resolution_latest_pointer, paths.repo_root
        ),
        "clinical_shortlist_latest_json": repo_relative(
            paths.clinical_shortlist_latest_pointer, paths.repo_root
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    baseline_profile_v1_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    processed_run_dir = paths.processed_runs_root / baseline_profile_v1_run_id
    audit_run_dir = paths.audit_runs_root / baseline_profile_v1_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)

        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        field_profiles = build_field_profiles(workflow_inputs)
        candidate_rows = build_candidate_rows(field_profiles)
        missingness_rows = build_missingness_ranked_rows(field_profiles)
        value_summary_rows = build_value_summary_rows(field_profiles, workflow_inputs.baseline_rows)
        excluded_rows = build_excluded_rows(candidate_rows)
        summary_rows = build_summary_rows(
            baseline_profile_v1_run_id=baseline_profile_v1_run_id,
            workflow_inputs=workflow_inputs,
            field_profiles=field_profiles,
            candidate_rows=candidate_rows,
            value_summary_rows=value_summary_rows,
        )

        if not field_profiles:
            raise BaselineProfileV1Error("baseline_profile_v1_field_summary.tsv rows were not generated.")
        if not missingness_rows:
            raise BaselineProfileV1Error(
                "baseline_profile_v1_missingness_ranked.tsv rows were not generated."
            )
        if not value_summary_rows:
            raise BaselineProfileV1Error("baseline_profile_v1_value_summary.tsv rows were not generated.")
        if not candidate_rows:
            raise BaselineProfileV1Error("baseline_profile_v1_candidate_fields.tsv rows were not generated.")
        if not summary_rows:
            raise BaselineProfileV1Error("baseline_profile_v1_summary.tsv rows were not generated.")

        field_summary_path = processed_run_dir / "baseline_profile_v1_field_summary.tsv"
        missingness_ranked_path = processed_run_dir / "baseline_profile_v1_missingness_ranked.tsv"
        value_summary_path = processed_run_dir / "baseline_profile_v1_value_summary.tsv"
        candidate_fields_path = audit_run_dir / "baseline_profile_v1_candidate_fields.tsv"
        excluded_fields_path = audit_run_dir / "baseline_profile_v1_excluded_fields.tsv"
        summary_path = audit_run_dir / "baseline_profile_v1_summary.tsv"

        helper_module.write_dict_rows_tsv(
            field_summary_path, FIELD_SUMMARY_FIELDNAMES, project_rows(field_profiles, FIELD_SUMMARY_FIELDNAMES)
        )
        helper_module.write_dict_rows_tsv(
            missingness_ranked_path, MISSINGNESS_RANKED_FIELDNAMES, missingness_rows
        )
        helper_module.write_dict_rows_tsv(value_summary_path, VALUE_SUMMARY_FIELDNAMES, value_summary_rows)
        helper_module.write_dict_rows_tsv(candidate_fields_path, CANDIDATE_FIELDS_FIELDNAMES, candidate_rows)
        helper_module.write_dict_rows_tsv(excluded_fields_path, EXCLUDED_FIELDS_FIELDNAMES, excluded_rows)
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        summary_lookup = build_summary_lookup(summary_rows)
        field_summary_row_count_matches_retained_fields = len(field_profiles) == len(workflow_inputs.baseline_rows[0])
        candidate_bucket_counts_match_retained_fields = len(candidate_rows) == len(field_profiles)
        excluded_table_matches_exclude_bucket = len(excluded_rows) == sum(
            row["candidate_bucket"] == "exclude_for_now" for row in candidate_rows
        )
        output_rows_positive = len(field_profiles) > 0
        required_source_table_found = workflow_inputs.input_paths["baseline_analysis_v1_tsv"].exists()
        summary_contains_endpoint_block = (
            summary_lookup["endpoint_freeze_status"]["summary_value"] == "blocked"
        )
        summary_contains_treatment_excluded = (
            summary_lookup["treatment_inclusion_status"]["summary_value"] == "excluded"
        )
        baseline_row_count_matches_input = len(workflow_inputs.baseline_rows) == parse_int(
            summary_lookup["final_row_count"]["summary_value"], "baseline profile final_row_count"
        )

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "field_summary_tsv": field_summary_path,
            "missingness_ranked_tsv": missingness_ranked_path,
            "value_summary_tsv": value_summary_path,
            "candidate_fields_tsv": candidate_fields_path,
            "excluded_fields_tsv": excluded_fields_path,
            "summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            baseline_profile_v1_run_id=baseline_profile_v1_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "baseline_analysis_v1_run_id": str(
                workflow_inputs.baseline_latest_pointer["baseline_analysis_v1_run_id"]
            ),
            "cohort_v1_build_id": str(workflow_inputs.baseline_latest_pointer["cohort_v1_build_id"]),
            "dry_run_build_id": str(workflow_inputs.baseline_latest_pointer["dry_run_build_id"]),
            "blueprint_run_id": str(workflow_inputs.baseline_latest_pointer["blueprint_run_id"]),
            "ambiguity_resolution_run_id": str(
                workflow_inputs.baseline_latest_pointer["ambiguity_resolution_run_id"]
            ),
            "shortlist_run_id": str(workflow_inputs.baseline_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(workflow_inputs.baseline_latest_pointer["core_audit_run_id"]),
            "clinical_parse_run_id": str(workflow_inputs.baseline_latest_pointer["clinical_parse_run_id"]),
            "endpoint_crosswalk_run_id": str(
                workflow_inputs.baseline_latest_pointer["endpoint_crosswalk_run_id"]
            ),
            "biospecimen_crosswalk_run_id": str(
                workflow_inputs.baseline_latest_pointer["biospecimen_crosswalk_run_id"]
            ),
            "biospecimen_parse_run_id": str(workflow_inputs.baseline_latest_pointer["biospecimen_parse_run_id"]),
            "clinical_source_run_id": str(workflow_inputs.baseline_latest_pointer["clinical_source_run_id"]),
            "biospecimen_source_run_id": str(workflow_inputs.baseline_latest_pointer["biospecimen_source_run_id"]),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "baseline_analysis_v1_latest_json": repo_relative(
                    paths.baseline_latest_pointer, paths.repo_root
                ),
                "minimal_cohort_v1_latest_json": repo_relative(
                    paths.minimal_cohort_latest_pointer, paths.repo_root
                ),
                "clinical_shortlist_latest_json": repo_relative(
                    paths.clinical_shortlist_latest_pointer, paths.repo_root
                ),
                "blueprint_ambiguity_resolution_latest_json": repo_relative(
                    paths.ambiguity_resolution_latest_pointer, paths.repo_root
                ),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "baseline_profile_v1_field_summary_tsv": repo_relative(
                    field_summary_path, paths.repo_root
                ),
                "baseline_profile_v1_missingness_ranked_tsv": repo_relative(
                    missingness_ranked_path, paths.repo_root
                ),
                "baseline_profile_v1_value_summary_tsv": repo_relative(
                    value_summary_path, paths.repo_root
                ),
                "baseline_profile_v1_candidate_fields_tsv": repo_relative(
                    candidate_fields_path, paths.repo_root
                ),
                "baseline_profile_v1_excluded_fields_tsv": repo_relative(
                    excluded_fields_path, paths.repo_root
                ),
                "baseline_profile_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    output_rows_positive
                    and required_source_table_found
                    and field_summary_row_count_matches_retained_fields
                    and candidate_bucket_counts_match_retained_fields
                    and excluded_table_matches_exclude_bucket
                    and baseline_row_count_matches_input
                    and summary_contains_endpoint_block
                    and summary_contains_treatment_excluded
                ),
                "required_pointers_found": True,
                "baseline_analysis_v1_latest_pointer_found": True,
                "baseline_analysis_run_log_completed": True,
                "baseline_analysis_validation_passed": True,
                "minimal_cohort_v1_latest_pointer_found": True,
                "minimal_cohort_run_log_completed": True,
                "minimal_cohort_validation_passed": True,
                "clinical_shortlist_latest_pointer_found": True,
                "clinical_shortlist_run_log_completed": True,
                "clinical_shortlist_validation_passed": True,
                "ambiguity_resolution_latest_pointer_found": True,
                "ambiguity_resolution_run_log_completed": True,
                "ambiguity_resolution_validation_passed": True,
                "required_source_table_found": required_source_table_found,
                "output_rows_positive": output_rows_positive,
                "field_summary_row_count_matches_retained_fields": (
                    field_summary_row_count_matches_retained_fields
                ),
                "candidate_bucket_counts_match_retained_fields": (
                    candidate_bucket_counts_match_retained_fields
                ),
                "excluded_table_matches_exclude_bucket": excluded_table_matches_exclude_bucket,
                "baseline_row_count_matches_input": baseline_row_count_matches_input,
                "summary_contains_endpoint_block": summary_contains_endpoint_block,
                "summary_contains_treatment_excluded": summary_contains_treatment_excluded,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "unit_of_analysis": "patient/case",
                "input_layer": "baseline_analysis_v1_only",
                "output_layer": "baseline_profile_v1",
                "candidate_selection_profile": "balanced",
                "direct_candidate_missing_like_fraction_max": BALANCED_MAX_MISSING_FRACTION,
                "direct_candidate_dominant_value_fraction_max": BALANCED_MAX_DOMINANT_FRACTION,
                "direct_candidate_distinct_non_missing_count_min": MIN_DISTINCT_NON_MISSING_COUNT,
                "review_candidate_missing_like_fraction_max": REVIEW_MAX_MISSING_FRACTION,
                "no_missing_value_imputation": True,
                "no_endpoint_freeze": True,
                "no_treatment_reintegration": True,
                "no_modeling": True,
                "no_metabric": True,
                "no_raw_xml_or_ssf_parsing": True,
                "one_row_per_patient": True,
                "followup_join_evidence_separate_from_endpoint_candidates": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "missing_like_tokens_json": json_list(sorted(MISSING_LIKE_TOKENS)),
            },
            "counts": {
                "input_baseline_analysis_v1_row_count": len(workflow_inputs.baseline_rows),
                "field_summary_row_count": len(field_profiles),
                "missingness_ranked_row_count": len(missingness_rows),
                "value_summary_row_count": len(value_summary_rows),
                "candidate_fields_row_count": len(candidate_rows),
                "excluded_fields_row_count": len(excluded_rows),
                "summary_row_count": len(summary_rows),
                "candidate_for_baseline_modeling_prep_count": sum(
                    row["candidate_bucket"] == "candidate_for_baseline_modeling_prep"
                    for row in candidate_rows
                ),
                "candidate_but_review_needed_count": sum(
                    row["candidate_bucket"] == "candidate_but_review_needed"
                    for row in candidate_rows
                ),
                "exclude_for_now_count": sum(
                    row["candidate_bucket"] == "exclude_for_now" for row in candidate_rows
                ),
                "rows_with_followup_evidence": count_yes(
                    workflow_inputs.baseline_rows, "row_has_followup_match"
                ),
                "rows_with_biospecimen_sample_evidence": count_yes(
                    workflow_inputs.baseline_rows, "row_has_biospecimen_sample_match"
                ),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "baseline_analysis_v1_latest_pointer": workflow_inputs.baseline_latest_pointer,
                "baseline_analysis_v1_run_log": workflow_inputs.baseline_run_log,
                "minimal_cohort_v1_latest_pointer": workflow_inputs.minimal_cohort_latest_pointer,
                "minimal_cohort_v1_run_log": workflow_inputs.minimal_cohort_run_log,
                "clinical_shortlist_latest_pointer": workflow_inputs.clinical_shortlist_latest_pointer,
                "clinical_shortlist_run_log": workflow_inputs.clinical_shortlist_run_log,
                "ambiguity_resolution_latest_pointer": workflow_inputs.ambiguity_resolution_latest_pointer,
                "ambiguity_resolution_run_log": workflow_inputs.ambiguity_resolution_run_log,
            },
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "baseline_profile_v1_run_id": baseline_profile_v1_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_baseline_profile_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA baseline profile v1 workflow complete.")
    print(f"Baseline profile run ID: {run_log['baseline_profile_v1_run_id']}")
    print(f"Baseline analysis run ID: {run_log['baseline_analysis_v1_run_id']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"Field summary TSV: {run_log['outputs']['baseline_profile_v1_field_summary_tsv']}")
    print(f"Candidate fields TSV: {run_log['outputs']['baseline_profile_v1_candidate_fields_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['baseline_profile_v1_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
