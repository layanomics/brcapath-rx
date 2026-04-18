#!/usr/bin/env python
"""Build a reproducible TCGA-BRCA baseline feature-set v1 from saved baseline-profile outputs."""

from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SPEC_FIELDNAMES = [
    "field_name",
    "source_origin",
    "field_category",
    "candidate_bucket",
    "feature_set_v1_decision",
    "decision_rule",
    "reason",
    "notes_placeholder",
]
MISSINGNESS_FIELDNAMES = [
    "field_name",
    "row_count",
    "non_missing_count",
    "missing_like_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "notes",
]
AUDIT_MAP_FIELDNAMES = [
    "feature_set_v1_row_index",
    "baseline_analysis_v1_row_id",
    "provisional_patient_row_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
]
SUMMARY_FIELDNAMES = [
    "baseline_feature_set_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
SPEC_NOTES_PLACEHOLDER = "[fill in during baseline feature-set v1 review]"
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
REVIEW_FIELD_DECISIONS = {
    "ethnicity": "include_in_feature_set_v1",
    "gender": "include_in_feature_set_v1",
    "icd_10": "exclude_from_feature_set_v1",
    "icd_o_3_site": "exclude_from_feature_set_v1",
}
REVIEW_FIELD_ORDER = ["ethnicity", "gender", "icd_10", "icd_o_3_site"]
EXPECTED_RETAINED_BASELINE_CLINICAL_FIELD_COUNT = 27
EXPECTED_DIRECT_CANDIDATE_COUNT = 23
EXPECTED_REVIEW_NEEDED_COUNT = 18
EXPECTED_EXCLUDED_FOR_NOW_COUNT = 13
EXPECTED_TOTAL_RETAINED_FIELD_COUNT = 54
FINAL_FEATURE_COLUMNS = [
    "age_at_diagnosis",
    "ajcc_metastasis_pathologic_pm",
    "ajcc_nodes_pathologic_pn",
    "ajcc_pathologic_tumor_stage",
    "ajcc_staging_edition",
    "ajcc_tumor_pathologic_pt",
    "anatomic_neoplasm_subdivision",
    "axillary_staging_method",
    "birth_days_to",
    "er_status_by_ihc",
    "ethnicity",
    "gender",
    "her2_status_by_ihc",
    "histological_type",
    "history_other_malignancy",
    "icd_o_3_histology",
    "initial_pathologic_dx_year",
    "lymph_nodes_examined_count",
    "lymph_nodes_examined_he_count",
    "margin_status",
    "menopause_status",
    "method_initial_path_dx",
    "pr_status_by_ihc",
    "race",
    "surgical_procedure_first",
]


class BaselineFeatureSetV1Error(RuntimeError):
    """Raised when the baseline feature-set workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the baseline feature-set workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    baseline_latest_pointer: Path
    baseline_profile_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved baseline layers."""

    baseline_latest_pointer: dict[str, Any]
    baseline_run_log: dict[str, Any]
    baseline_rows: list[dict[str, str]]
    baseline_profile_latest_pointer: dict[str, Any]
    baseline_profile_run_log: dict[str, Any]
    candidate_rows: list[dict[str, str]]
    excluded_rows: list[dict[str, str]]
    profile_summary_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise BaselineFeatureSetV1Error(f"Required helper script not found: {script_path}")

    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise BaselineFeatureSetV1Error(f"Unable to create an import spec for: {script_path}")

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
        raise BaselineFeatureSetV1Error(f"Expected integer-like value for {label}: {value!r}") from exc


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
            raise BaselineFeatureSetV1Error(f"Duplicate summary_metric detected: {metric}")
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


def require_pointer(
    pointer_path: Path,
    *,
    label: str,
    required_keys: set[str],
    helper_module: Any,
) -> dict[str, Any]:
    if not pointer_path.exists():
        raise BaselineFeatureSetV1Error(f"Required {label} not found: {pointer_path}")
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
        raise BaselineFeatureSetV1Error(f"{label} is not completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise BaselineFeatureSetV1Error(f"{label} does not report validation.passed == true.")
    return run_log_path, run_log


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
        processed_runs_root=(
            processed_root / "tcga-brca" / "analysis-prep" / "baseline_feature_set_v1_runs"
        ),
        audit_runs_root=analysis_prep_root / "baseline_feature_set_v1_runs",
        latest_pointer=analysis_prep_root / "tcga_brca_baseline_feature_set_v1_latest.json",
        baseline_latest_pointer=analysis_prep_root / "tcga_brca_baseline_analysis_v1_latest.json",
        baseline_profile_latest_pointer=analysis_prep_root / "tcga_brca_baseline_profile_v1_latest.json",
    )


def build_input_paths(
    paths: WorkflowPaths,
    *,
    baseline_latest_pointer: dict[str, Any],
    baseline_profile_latest_pointer: dict[str, Any],
    helper_module: Any,
) -> dict[str, Path]:
    return {
        "baseline_analysis_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_latest_pointer["baseline_analysis_v1_tsv"]),
            "baseline analysis v1 TSV",
        ),
        "baseline_analysis_run_log_json": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_latest_pointer["run_log_json"]),
            "baseline analysis v1 run log",
        ),
        "baseline_profile_v1_candidate_fields_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_profile_latest_pointer["baseline_profile_v1_candidate_fields_tsv"]),
            "baseline profile v1 candidate fields TSV",
        ),
        "baseline_profile_v1_excluded_fields_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_profile_latest_pointer["baseline_profile_v1_excluded_fields_tsv"]),
            "baseline profile v1 excluded fields TSV",
        ),
        "baseline_profile_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_profile_latest_pointer["baseline_profile_v1_summary_tsv"]),
            "baseline profile v1 summary TSV",
        ),
        "baseline_profile_v1_run_log_json": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_profile_latest_pointer["run_log_json"]),
            "baseline profile v1 run log",
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
            "run_log_json",
        },
        helper_module=helper_module,
    )
    baseline_profile_latest_pointer = require_pointer(
        paths.baseline_profile_latest_pointer,
        label="Baseline-profile latest pointer",
        required_keys={
            "baseline_profile_v1_run_id",
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
            "baseline_profile_v1_candidate_fields_tsv",
            "baseline_profile_v1_excluded_fields_tsv",
            "baseline_profile_v1_summary_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )

    expected_pairings = {
        "baseline_analysis_v1_run_id": str(baseline_latest_pointer["baseline_analysis_v1_run_id"]),
        "cohort_v1_build_id": str(baseline_latest_pointer["cohort_v1_build_id"]),
        "dry_run_build_id": str(baseline_latest_pointer["dry_run_build_id"]),
        "blueprint_run_id": str(baseline_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(baseline_latest_pointer["ambiguity_resolution_run_id"]),
        "shortlist_run_id": str(baseline_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(baseline_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(baseline_latest_pointer["clinical_parse_run_id"]),
        "endpoint_crosswalk_run_id": str(baseline_latest_pointer["endpoint_crosswalk_run_id"]),
        "biospecimen_crosswalk_run_id": str(baseline_latest_pointer["biospecimen_crosswalk_run_id"]),
        "biospecimen_parse_run_id": str(baseline_latest_pointer["biospecimen_parse_run_id"]),
        "clinical_source_run_id": str(baseline_latest_pointer["clinical_source_run_id"]),
        "biospecimen_source_run_id": str(baseline_latest_pointer["biospecimen_source_run_id"]),
    }
    for key, expected_value in expected_pairings.items():
        observed_value = str(baseline_profile_latest_pointer[key])
        if observed_value != expected_value:
            raise BaselineFeatureSetV1Error(
                "Baseline-profile latest pointer is out of sync with the current baseline-analysis latest pointer: "
                f"{key}={observed_value!r} vs expected {expected_value!r}."
            )

    input_paths = build_input_paths(
        paths,
        baseline_latest_pointer=baseline_latest_pointer,
        baseline_profile_latest_pointer=baseline_profile_latest_pointer,
        helper_module=helper_module,
    )

    _, baseline_run_log = require_completed_run_log(
        paths.repo_root,
        str(baseline_latest_pointer["run_log_json"]),
        label="Baseline-analysis-prep run log",
        helper_module=helper_module,
    )
    _, baseline_profile_run_log = require_completed_run_log(
        paths.repo_root,
        str(baseline_profile_latest_pointer["run_log_json"]),
        label="Baseline-profile run log",
        helper_module=helper_module,
    )

    baseline_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_analysis_v1_tsv"])
    candidate_rows = helper_module.read_tsv_dict_rows(
        input_paths["baseline_profile_v1_candidate_fields_tsv"]
    )
    excluded_rows = helper_module.read_tsv_dict_rows(
        input_paths["baseline_profile_v1_excluded_fields_tsv"]
    )
    profile_summary_rows = helper_module.read_tsv_dict_rows(
        input_paths["baseline_profile_v1_summary_tsv"]
    )

    if not baseline_rows:
        raise BaselineFeatureSetV1Error("The referenced baseline_analysis_v1.tsv contains no rows.")
    if not candidate_rows:
        raise BaselineFeatureSetV1Error(
            "The referenced baseline_profile_v1_candidate_fields.tsv contains no rows."
        )
    if not profile_summary_rows:
        raise BaselineFeatureSetV1Error(
            "The referenced baseline_profile_v1_summary.tsv contains no rows."
        )

    return WorkflowInputs(
        baseline_latest_pointer=baseline_latest_pointer,
        baseline_run_log=baseline_run_log,
        baseline_rows=baseline_rows,
        baseline_profile_latest_pointer=baseline_profile_latest_pointer,
        baseline_profile_run_log=baseline_profile_run_log,
        candidate_rows=candidate_rows,
        excluded_rows=excluded_rows,
        profile_summary_rows=profile_summary_rows,
        input_paths=input_paths,
    )


def build_candidate_lookup(candidate_rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in candidate_rows:
        field_name = str(row.get("field_name") or "")
        if not field_name:
            raise BaselineFeatureSetV1Error("Encountered an empty field_name in baseline profile candidate rows.")
        if field_name in lookup:
            raise BaselineFeatureSetV1Error(f"Duplicate candidate field row detected: {field_name}")
        lookup[field_name] = row
    return lookup


def sort_baseline_rows(baseline_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    sorted_rows = sorted(
        baseline_rows,
        key=lambda row: parse_int(row["baseline_analysis_v1_row_id"], "baseline_analysis_v1_row_id"),
    )
    observed_ids = [
        parse_int(row["baseline_analysis_v1_row_id"], "baseline_analysis_v1_row_id")
        for row in sorted_rows
    ]
    if len(set(observed_ids)) != len(observed_ids):
        raise BaselineFeatureSetV1Error(
            "baseline_analysis_v1.tsv does not preserve one row per baseline_analysis_v1_row_id."
        )
    if observed_ids != list(range(1, len(sorted_rows) + 1)):
        raise BaselineFeatureSetV1Error(
            "baseline_analysis_v1.tsv does not preserve sequential baseline_analysis_v1_row_id values."
        )
    return sorted_rows


def validate_profile_state(workflow_inputs: WorkflowInputs) -> None:
    baseline_column_order = list(workflow_inputs.baseline_rows[0].keys())
    candidate_lookup = build_candidate_lookup(workflow_inputs.candidate_rows)

    if len(candidate_lookup) != EXPECTED_TOTAL_RETAINED_FIELD_COUNT:
        raise BaselineFeatureSetV1Error(
            "baseline_profile_v1_candidate_fields.tsv did not contain the expected retained field count: "
            f"{len(candidate_lookup)} vs expected {EXPECTED_TOTAL_RETAINED_FIELD_COUNT}."
        )
    if len(baseline_column_order) != EXPECTED_TOTAL_RETAINED_FIELD_COUNT:
        raise BaselineFeatureSetV1Error(
            "baseline_analysis_v1.tsv did not contain the expected retained column count: "
            f"{len(baseline_column_order)} vs expected {EXPECTED_TOTAL_RETAINED_FIELD_COUNT}."
        )
    if set(candidate_lookup) != set(baseline_column_order):
        raise BaselineFeatureSetV1Error(
            "baseline_profile_v1_candidate_fields.tsv did not reconcile exactly to baseline_analysis_v1.tsv columns."
        )

    baseline_clinical_rows = [
        row for row in workflow_inputs.candidate_rows if row["field_category"] == "baseline_clinical_field"
    ]
    direct_candidate_rows = [
        row
        for row in baseline_clinical_rows
        if row["candidate_bucket"] == "candidate_for_baseline_modeling_prep"
    ]
    review_needed_baseline_rows = [
        row
        for row in baseline_clinical_rows
        if row["candidate_bucket"] == "candidate_but_review_needed"
    ]
    baseline_excluded_rows = [
        row for row in baseline_clinical_rows if row["candidate_bucket"] == "exclude_for_now"
    ]

    if len(baseline_clinical_rows) != EXPECTED_RETAINED_BASELINE_CLINICAL_FIELD_COUNT:
        raise BaselineFeatureSetV1Error(
            "Unexpected retained baseline clinical field count in the saved profile outputs: "
            f"{len(baseline_clinical_rows)} vs expected {EXPECTED_RETAINED_BASELINE_CLINICAL_FIELD_COUNT}."
        )
    if len(direct_candidate_rows) != EXPECTED_DIRECT_CANDIDATE_COUNT:
        raise BaselineFeatureSetV1Error(
            "Unexpected direct baseline candidate count in the saved profile outputs: "
            f"{len(direct_candidate_rows)} vs expected {EXPECTED_DIRECT_CANDIDATE_COUNT}."
        )
    if len(review_needed_baseline_rows) != len(REVIEW_FIELD_DECISIONS):
        raise BaselineFeatureSetV1Error(
            "Unexpected review-needed baseline clinical field count in the saved profile outputs: "
            f"{len(review_needed_baseline_rows)} vs expected {len(REVIEW_FIELD_DECISIONS)}."
        )
    if baseline_excluded_rows:
        raise BaselineFeatureSetV1Error(
            "Encountered unexpected baseline clinical fields already marked exclude_for_now in the saved profile "
            f"outputs: {sorted(row['field_name'] for row in baseline_excluded_rows)}"
        )

    observed_review_fields = ordered_unique([row["field_name"] for row in review_needed_baseline_rows])
    if observed_review_fields != REVIEW_FIELD_ORDER:
        raise BaselineFeatureSetV1Error(
            "Unexpected review-needed baseline clinical fields in the saved profile outputs: "
            f"{observed_review_fields!r} vs expected {REVIEW_FIELD_ORDER!r}."
        )

    exclude_bucket_rows = [
        row for row in workflow_inputs.candidate_rows if row["candidate_bucket"] == "exclude_for_now"
    ]
    excluded_lookup = build_candidate_lookup(workflow_inputs.excluded_rows)
    if len(exclude_bucket_rows) != EXPECTED_EXCLUDED_FOR_NOW_COUNT:
        raise BaselineFeatureSetV1Error(
            "Unexpected exclude_for_now field count in the saved profile outputs: "
            f"{len(exclude_bucket_rows)} vs expected {EXPECTED_EXCLUDED_FOR_NOW_COUNT}."
        )
    if len(workflow_inputs.excluded_rows) != len(exclude_bucket_rows):
        raise BaselineFeatureSetV1Error(
            "baseline_profile_v1_excluded_fields.tsv row count did not match candidate_bucket=exclude_for_now rows."
        )
    if set(excluded_lookup) != {row["field_name"] for row in exclude_bucket_rows}:
        raise BaselineFeatureSetV1Error(
            "baseline_profile_v1_excluded_fields.tsv did not reconcile exactly to candidate_bucket=exclude_for_now "
            "fields."
        )

    summary_lookup = build_summary_lookup(workflow_inputs.profile_summary_rows)
    expected_summary_counts = {
        "final_row_count": len(workflow_inputs.baseline_rows),
        "retained_baseline_clinical_field_count": EXPECTED_RETAINED_BASELINE_CLINICAL_FIELD_COUNT,
        "candidate_field_count": EXPECTED_DIRECT_CANDIDATE_COUNT,
        "review_needed_field_count": EXPECTED_REVIEW_NEEDED_COUNT,
        "excluded_for_now_field_count": EXPECTED_EXCLUDED_FOR_NOW_COUNT,
    }
    for metric, expected_value in expected_summary_counts.items():
        if metric not in summary_lookup:
            raise BaselineFeatureSetV1Error(
                f"baseline_profile_v1_summary.tsv is missing summary_metric={metric}."
            )
        observed_value = parse_int(summary_lookup[metric]["summary_value"], metric)
        if observed_value != expected_value:
            raise BaselineFeatureSetV1Error(
                "baseline_profile_v1_summary.tsv did not match the expected saved profile state: "
                f"{metric}={observed_value} vs expected {expected_value}."
            )


def fixed_review_decision_reason(field_name: str, candidate_row: dict[str, str]) -> str:
    base_reason = str(candidate_row["reason"])
    if field_name == "ethnicity":
        return (
            "Fixed review-field decision: keep ethnicity in baseline feature-set v1 because the saved profile "
            "missingness remains acceptable for first-pass baseline prep even though the non-missing distribution is "
            f"near-constant. {base_reason}"
        )
    if field_name == "gender":
        return (
            "Fixed review-field decision: keep gender in baseline feature-set v1 because it is complete and remains an "
            f"explicitly retained baseline demographic field. {base_reason}"
        )
    if field_name == "icd_10":
        return (
            "Fixed review-field decision: keep icd_10 out of baseline feature-set v1 because the saved profile shows a "
            f"strongly dominant code distribution for this first-pass prep layer. {base_reason}"
        )
    if field_name == "icd_o_3_site":
        return (
            "Fixed review-field decision: keep icd_o_3_site out of baseline feature-set v1 because the saved profile "
            f"shows a strongly dominant site-code distribution for this first-pass prep layer. {base_reason}"
        )
    raise BaselineFeatureSetV1Error(f"Unexpected review field for fixed decision reasoning: {field_name}")


def build_feature_set_rows(
    baseline_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    sorted_rows = sort_baseline_rows(baseline_rows)
    feature_rows: list[dict[str, str]] = []
    audit_map_rows: list[dict[str, str]] = []

    for row_index, source_row in enumerate(sorted_rows, start=1):
        feature_rows.append(
            {field_name: str(source_row.get(field_name, "")) for field_name in FINAL_FEATURE_COLUMNS}
        )
        audit_map_rows.append(
            {
                "feature_set_v1_row_index": str(row_index),
                "baseline_analysis_v1_row_id": str(source_row["baseline_analysis_v1_row_id"]),
                "provisional_patient_row_id": str(source_row["provisional_patient_row_id"]),
                "bcr_patient_barcode": str(source_row["bcr_patient_barcode"]),
                "bcr_patient_uuid": str(source_row["bcr_patient_uuid"]),
            }
        )

    return feature_rows, audit_map_rows


def build_spec_rows(
    *,
    baseline_column_order: list[str],
    candidate_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    candidate_lookup = build_candidate_lookup(candidate_rows)
    spec_rows: list[dict[str, str]] = []

    for field_name in baseline_column_order:
        if field_name not in candidate_lookup:
            raise BaselineFeatureSetV1Error(
                f"Saved profile candidate rows are missing field_name={field_name!r}."
            )
        candidate_row = candidate_lookup[field_name]
        field_category = str(candidate_row["field_category"])
        candidate_bucket = str(candidate_row["candidate_bucket"])

        if field_name in REVIEW_FIELD_DECISIONS:
            feature_decision = REVIEW_FIELD_DECISIONS[field_name]
            decision_rule = "fixed_review_field_decision"
            reason = fixed_review_decision_reason(field_name, candidate_row)
        elif field_category == "baseline_clinical_field" and candidate_bucket == "candidate_for_baseline_modeling_prep":
            feature_decision = "include_in_feature_set_v1"
            decision_rule = "auto_include_candidate_bucket"
            reason = (
                "Auto-included because the saved baseline-profile output placed this field in "
                f"candidate_for_baseline_modeling_prep. {candidate_row['reason']}"
            )
        elif field_category == "baseline_clinical_field":
            raise BaselineFeatureSetV1Error(
                "Encountered an unexpected baseline clinical field outside the approved direct-candidate and fixed "
                f"review-decision paths: {field_name}"
            )
        else:
            feature_decision = "exclude_from_feature_set_v1"
            decision_rule = "excluded_non_feature_category"
            reason = str(candidate_row["reason"])

        spec_rows.append(
            {
                "field_name": field_name,
                "source_origin": str(candidate_row["source_origin"]),
                "field_category": field_category,
                "candidate_bucket": candidate_bucket,
                "feature_set_v1_decision": feature_decision,
                "decision_rule": decision_rule,
                "reason": reason,
                "notes_placeholder": SPEC_NOTES_PLACEHOLDER,
            }
        )

    if len(spec_rows) != EXPECTED_TOTAL_RETAINED_FIELD_COUNT:
        raise BaselineFeatureSetV1Error(
            "baseline_feature_set_v1_spec.tsv row count did not reconcile to the retained baseline profile fields."
        )

    return spec_rows


def build_missingness_rows(
    feature_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    row_count = len(feature_rows)
    missingness_rows: list[dict[str, str]] = []

    for field_name in FINAL_FEATURE_COLUMNS:
        non_missing_values: list[str] = []
        missing_like_count = 0
        for row in feature_rows:
            raw_value = str(row.get(field_name, ""))
            if cell_is_missing_like(field_name, raw_value):
                missing_like_count += 1
                continue
            non_missing_values.append(raw_value)

        missingness_rows.append(
            {
                "field_name": field_name,
                "row_count": str(row_count),
                "non_missing_count": str(len(non_missing_values)),
                "missing_like_count": str(missing_like_count),
                "missing_like_fraction": f"{(missing_like_count / row_count) if row_count else 0.0:.6f}",
                "distinct_non_missing_count": str(len(set(non_missing_values))),
                "notes": "scalar_feature",
            }
        )

    return missingness_rows


def build_review_decision_outcomes(
    candidate_rows: list[dict[str, str]],
    spec_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    candidate_lookup = build_candidate_lookup(candidate_rows)
    spec_lookup = {row["field_name"]: row for row in spec_rows}
    outcomes: list[dict[str, str]] = []

    for field_name in REVIEW_FIELD_ORDER:
        candidate_row = candidate_lookup[field_name]
        spec_row = spec_lookup[field_name]
        outcomes.append(
            {
                "field_name": field_name,
                "feature_set_v1_decision": str(spec_row["feature_set_v1_decision"]),
                "candidate_bucket": str(candidate_row["candidate_bucket"]),
                "missing_like_fraction": str(candidate_row["missing_like_fraction"]),
                "dominant_value_fraction": str(candidate_row["dominant_value_fraction"]),
                "distinct_non_missing_count": str(candidate_row["distinct_non_missing_count"]),
                "reason": str(spec_row["reason"]),
            }
        )

    return outcomes


def build_summary_rows(
    *,
    baseline_feature_set_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
    spec_rows: list[dict[str, str]],
    feature_rows: list[dict[str, str]],
    review_decision_outcomes: list[dict[str, str]],
) -> list[dict[str, str]]:
    included_rows = [row for row in spec_rows if row["feature_set_v1_decision"] == "include_in_feature_set_v1"]
    excluded_rows = [row for row in spec_rows if row["feature_set_v1_decision"] == "exclude_from_feature_set_v1"]
    auto_included_rows = [row for row in included_rows if row["decision_rule"] == "auto_include_candidate_bucket"]
    review_decided_included_rows = [
        row
        for row in included_rows
        if row["decision_rule"] == "fixed_review_field_decision"
    ]
    review_decided_excluded_rows = [
        row
        for row in excluded_rows
        if row["decision_rule"] == "fixed_review_field_decision"
    ]

    return [
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "design",
            "summary_metric": "baseline_feature_set_v1_status",
            "summary_value": "baseline_feature_set_v1_complete",
            "notes": "This workflow builds a first clean baseline feature-set table only; it does not perform modeling.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "design",
            "summary_metric": "feature_source_table",
            "summary_value": "baseline_analysis_v1",
            "notes": "The feature-set starts only from the saved baseline_analysis_v1 TSV on disk.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "design",
            "summary_metric": "feature_selection_scope",
            "summary_value": "baseline_clinical_fields_only",
            "notes": "Endpoint, follow-up evidence, biospecimen evidence, audit fields, and derived flags are kept out.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "design",
            "summary_metric": "unit_of_analysis",
            "summary_value": "patient/case",
            "notes": "The feature-set preserves one row per patient/case from baseline_analysis_v1.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "design",
            "summary_metric": "review_field_decision_policy",
            "summary_value": json.dumps(REVIEW_FIELD_DECISIONS, ensure_ascii=True, sort_keys=True),
            "notes": "Only the four approved review-needed baseline clinical fields receive fixed include/exclude decisions.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "design",
            "summary_metric": "endpoint_policy",
            "summary_value": "excluded_not_frozen",
            "notes": "Endpoint candidates remain blocked and are not part of baseline feature-set v1.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "design",
            "summary_metric": "treatment_policy",
            "summary_value": "excluded",
            "notes": "Treatment detail remains excluded from this feature-set workflow.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "baseline_analysis_v1_run_id",
            "summary_value": str(workflow_inputs.baseline_latest_pointer["baseline_analysis_v1_run_id"]),
            "notes": "Baseline-analysis-prep run used as the patient-level source table.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "baseline_profile_v1_run_id",
            "summary_value": str(workflow_inputs.baseline_profile_latest_pointer["baseline_profile_v1_run_id"]),
            "notes": "Baseline-profile run used as the saved field-selection evidence layer.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "cohort_v1_build_id",
            "summary_value": str(workflow_inputs.baseline_profile_latest_pointer["cohort_v1_build_id"]),
            "notes": "Source minimal cohort v1 build carried forward through baseline_analysis_v1.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "ambiguity_resolution_run_id",
            "summary_value": str(workflow_inputs.baseline_profile_latest_pointer["ambiguity_resolution_run_id"]),
            "notes": "Current ambiguity-resolution run aligned to the saved baseline inputs.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "shortlist_run_id",
            "summary_value": str(workflow_inputs.baseline_profile_latest_pointer["shortlist_run_id"]),
            "notes": "Current clinical shortlist run aligned to the saved baseline inputs.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "input_baseline_analysis_v1_row_count",
            "summary_value": str(len(workflow_inputs.baseline_rows)),
            "notes": "Row count read directly from baseline_analysis_v1.tsv.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "final_row_count",
            "summary_value": str(len(feature_rows)),
            "notes": "Final patient-level row count represented in baseline feature-set v1.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "included_feature_count",
            "summary_value": str(len(included_rows)),
            "notes": "Total number of retained baseline clinical fields included in feature-set v1.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "auto_included_feature_count",
            "summary_value": str(len(auto_included_rows)),
            "notes": "Included fields auto-selected from candidate_for_baseline_modeling_prep.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "review_decided_included_feature_count",
            "summary_value": str(len(review_decided_included_rows)),
            "notes": "Included fields admitted only through the fixed four-field review-decision block.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "review_decided_excluded_feature_count",
            "summary_value": str(len(review_decided_excluded_rows)),
            "notes": "Excluded baseline clinical fields resolved only through the fixed four-field review-decision block.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "included_feature_names_json",
            "summary_value": json_list([row["field_name"] for row in included_rows]),
            "notes": "Feature names included in baseline_feature_set_v1.tsv, in saved column order.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "fields_deliberately_kept_out_json",
            "summary_value": json_list([row["field_name"] for row in excluded_rows]),
            "notes": "Retained baseline-profile fields deliberately excluded from baseline_feature_set_v1.tsv.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "review_field_decision_outcomes_json",
            "summary_value": json.dumps(review_decision_outcomes, ensure_ascii=True),
            "notes": "Fixed include/exclude outcomes for ethnicity, gender, icd_10, and icd_o_3_site.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "readiness_for_first_pass_baseline_model_input",
            "summary_value": "ready",
            "notes": "The saved feature-set is ready as a first-pass baseline model-input table only.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "endpoint_freeze_status",
            "summary_value": "blocked",
            "notes": "Endpoint freeze remains blocked and endpoint candidate fields stay out of feature-set v1.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "treatment_inclusion_status",
            "summary_value": "excluded",
            "notes": "Treatment detail remains excluded from this feature-set workflow.",
        },
        {
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "overall_readiness_interpretation",
            "summary_value": "suitable_for_first_pass_baseline_model_input_only",
            "notes": (
                "Use this layer for first-pass baseline model input preparation, descriptive review, and missingness "
                "review only; do not treat it as endpoint freeze, treatment modeling, or model training."
            ),
        },
    ]


def build_latest_pointer_payload(
    *,
    baseline_feature_set_v1_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
        "baseline_profile_v1_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["baseline_profile_v1_run_id"]
        ),
        "baseline_analysis_v1_run_id": str(
            workflow_inputs.baseline_latest_pointer["baseline_analysis_v1_run_id"]
        ),
        "cohort_v1_build_id": str(workflow_inputs.baseline_profile_latest_pointer["cohort_v1_build_id"]),
        "dry_run_build_id": str(workflow_inputs.baseline_profile_latest_pointer["dry_run_build_id"]),
        "blueprint_run_id": str(workflow_inputs.baseline_profile_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(workflow_inputs.baseline_profile_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(workflow_inputs.baseline_profile_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["clinical_parse_run_id"]
        ),
        "endpoint_crosswalk_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["endpoint_crosswalk_run_id"]
        ),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["biospecimen_crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["biospecimen_parse_run_id"]
        ),
        "clinical_source_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["clinical_source_run_id"]
        ),
        "biospecimen_source_run_id": str(
            workflow_inputs.baseline_profile_latest_pointer["biospecimen_source_run_id"]
        ),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "baseline_feature_set_v1_tsv": repo_relative(
            output_paths["baseline_feature_set_v1_tsv"], paths.repo_root
        ),
        "baseline_feature_set_v1_spec_tsv": repo_relative(
            output_paths["baseline_feature_set_v1_spec_tsv"], paths.repo_root
        ),
        "baseline_feature_set_v1_missingness_tsv": repo_relative(
            output_paths["baseline_feature_set_v1_missingness_tsv"], paths.repo_root
        ),
        "baseline_feature_set_v1_summary_tsv": repo_relative(
            output_paths["baseline_feature_set_v1_summary_tsv"], paths.repo_root
        ),
        "baseline_feature_set_v1_audit_map_tsv": repo_relative(
            output_paths["baseline_feature_set_v1_audit_map_tsv"], paths.repo_root
        ),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "baseline_analysis_v1_latest_json": repo_relative(paths.baseline_latest_pointer, paths.repo_root),
        "baseline_profile_v1_latest_json": repo_relative(
            paths.baseline_profile_latest_pointer, paths.repo_root
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    baseline_feature_set_v1_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    processed_run_dir = paths.processed_runs_root / baseline_feature_set_v1_run_id
    audit_run_dir = paths.audit_runs_root / baseline_feature_set_v1_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)
        validate_profile_state(workflow_inputs)

        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        feature_set_tsv_path = processed_run_dir / "baseline_feature_set_v1.tsv"
        spec_path = audit_run_dir / "baseline_feature_set_v1_spec.tsv"
        missingness_path = audit_run_dir / "baseline_feature_set_v1_missingness.tsv"
        summary_path = audit_run_dir / "baseline_feature_set_v1_summary.tsv"
        audit_map_path = audit_run_dir / "baseline_feature_set_v1_audit_map.tsv"

        baseline_column_order = list(workflow_inputs.baseline_rows[0].keys())
        spec_rows = build_spec_rows(
            baseline_column_order=baseline_column_order,
            candidate_rows=workflow_inputs.candidate_rows,
        )
        feature_rows, audit_map_rows = build_feature_set_rows(workflow_inputs.baseline_rows)
        missingness_rows = build_missingness_rows(feature_rows)
        review_decision_outcomes = build_review_decision_outcomes(workflow_inputs.candidate_rows, spec_rows)
        summary_rows = build_summary_rows(
            baseline_feature_set_v1_run_id=baseline_feature_set_v1_run_id,
            workflow_inputs=workflow_inputs,
            spec_rows=spec_rows,
            feature_rows=feature_rows,
            review_decision_outcomes=review_decision_outcomes,
        )

        if not feature_rows:
            raise BaselineFeatureSetV1Error("baseline_feature_set_v1.tsv rows were not generated.")
        if not spec_rows:
            raise BaselineFeatureSetV1Error("baseline_feature_set_v1_spec.tsv rows were not generated.")
        if not missingness_rows:
            raise BaselineFeatureSetV1Error("baseline_feature_set_v1_missingness.tsv rows were not generated.")
        if not summary_rows:
            raise BaselineFeatureSetV1Error("baseline_feature_set_v1_summary.tsv rows were not generated.")
        if not audit_map_rows:
            raise BaselineFeatureSetV1Error("baseline_feature_set_v1_audit_map.tsv rows were not generated.")

        helper_module.write_dict_rows_tsv(feature_set_tsv_path, FINAL_FEATURE_COLUMNS, feature_rows)
        helper_module.write_dict_rows_tsv(spec_path, SPEC_FIELDNAMES, spec_rows)
        helper_module.write_dict_rows_tsv(missingness_path, MISSINGNESS_FIELDNAMES, missingness_rows)
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)
        helper_module.write_dict_rows_tsv(audit_map_path, AUDIT_MAP_FIELDNAMES, audit_map_rows)

        included_spec_rows = [row for row in spec_rows if row["feature_set_v1_decision"] == "include_in_feature_set_v1"]
        included_spec_field_names = [row["field_name"] for row in included_spec_rows]
        included_feature_columns_match_spec = included_spec_field_names == FINAL_FEATURE_COLUMNS
        feature_row_count_matches_baseline_analysis = len(feature_rows) == len(workflow_inputs.baseline_rows)
        audit_map_row_count_matches_feature_table = len(audit_map_rows) == len(feature_rows)
        feature_row_index_values = [
            parse_int(row["feature_set_v1_row_index"], "feature_set_v1_row_index") for row in audit_map_rows
        ]
        feature_row_index_sequential = feature_row_index_values == list(range(1, len(audit_map_rows) + 1))
        included_categories_baseline_only = all(
            row["field_category"] == "baseline_clinical_field" for row in included_spec_rows
        )
        required_source_tables_found = all(path.exists() for path in workflow_inputs.input_paths.values())
        output_rows_positive = len(feature_rows) > 0
        summary_lookup = build_summary_lookup(summary_rows)
        summary_endpoint_blocked = summary_lookup["endpoint_freeze_status"]["summary_value"] == "blocked"
        summary_treatment_excluded = (
            summary_lookup["treatment_inclusion_status"]["summary_value"] == "excluded"
        )
        auto_included_count = sum(row["decision_rule"] == "auto_include_candidate_bucket" for row in included_spec_rows)
        review_decided_included_count = sum(
            row["decision_rule"] == "fixed_review_field_decision" for row in included_spec_rows
        )
        review_decided_excluded_count = sum(
            row["decision_rule"] == "fixed_review_field_decision"
            and row["feature_set_v1_decision"] == "exclude_from_feature_set_v1"
            for row in spec_rows
        )
        final_feature_count_matches_expected = len(included_spec_rows) == len(FINAL_FEATURE_COLUMNS)
        auto_included_count_matches_expected = auto_included_count == EXPECTED_DIRECT_CANDIDATE_COUNT
        review_included_count_matches_expected = review_decided_included_count == 2
        review_excluded_count_matches_expected = review_decided_excluded_count == 2
        no_unexpected_review_fields = [row["field_name"] for row in review_decision_outcomes] == REVIEW_FIELD_ORDER

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "baseline_feature_set_v1_tsv": feature_set_tsv_path,
            "baseline_feature_set_v1_spec_tsv": spec_path,
            "baseline_feature_set_v1_missingness_tsv": missingness_path,
            "baseline_feature_set_v1_summary_tsv": summary_path,
            "baseline_feature_set_v1_audit_map_tsv": audit_map_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            baseline_feature_set_v1_run_id=baseline_feature_set_v1_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "baseline_profile_v1_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["baseline_profile_v1_run_id"]
            ),
            "baseline_analysis_v1_run_id": str(
                workflow_inputs.baseline_latest_pointer["baseline_analysis_v1_run_id"]
            ),
            "cohort_v1_build_id": str(workflow_inputs.baseline_profile_latest_pointer["cohort_v1_build_id"]),
            "dry_run_build_id": str(workflow_inputs.baseline_profile_latest_pointer["dry_run_build_id"]),
            "blueprint_run_id": str(workflow_inputs.baseline_profile_latest_pointer["blueprint_run_id"]),
            "ambiguity_resolution_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["ambiguity_resolution_run_id"]
            ),
            "shortlist_run_id": str(workflow_inputs.baseline_profile_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(workflow_inputs.baseline_profile_latest_pointer["core_audit_run_id"]),
            "clinical_parse_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["clinical_parse_run_id"]
            ),
            "endpoint_crosswalk_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["endpoint_crosswalk_run_id"]
            ),
            "biospecimen_crosswalk_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["biospecimen_crosswalk_run_id"]
            ),
            "biospecimen_parse_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["biospecimen_parse_run_id"]
            ),
            "clinical_source_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["clinical_source_run_id"]
            ),
            "biospecimen_source_run_id": str(
                workflow_inputs.baseline_profile_latest_pointer["biospecimen_source_run_id"]
            ),
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
                "baseline_profile_v1_latest_json": repo_relative(
                    paths.baseline_profile_latest_pointer, paths.repo_root
                ),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "baseline_feature_set_v1_tsv": repo_relative(feature_set_tsv_path, paths.repo_root),
                "baseline_feature_set_v1_spec_tsv": repo_relative(spec_path, paths.repo_root),
                "baseline_feature_set_v1_missingness_tsv": repo_relative(
                    missingness_path, paths.repo_root
                ),
                "baseline_feature_set_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "baseline_feature_set_v1_audit_map_tsv": repo_relative(audit_map_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    required_source_tables_found
                    and output_rows_positive
                    and feature_row_count_matches_baseline_analysis
                    and included_feature_columns_match_spec
                    and audit_map_row_count_matches_feature_table
                    and feature_row_index_sequential
                    and included_categories_baseline_only
                    and summary_endpoint_blocked
                    and summary_treatment_excluded
                    and final_feature_count_matches_expected
                    and auto_included_count_matches_expected
                    and review_included_count_matches_expected
                    and review_excluded_count_matches_expected
                    and no_unexpected_review_fields
                ),
                "required_pointers_found": True,
                "baseline_analysis_v1_latest_pointer_found": True,
                "baseline_profile_v1_latest_pointer_found": True,
                "baseline_analysis_run_log_completed": True,
                "baseline_analysis_validation_passed": True,
                "baseline_profile_run_log_completed": True,
                "baseline_profile_validation_passed": True,
                "required_source_tables_found": required_source_tables_found,
                "output_rows_positive": output_rows_positive,
                "feature_row_count_matches_baseline_analysis": feature_row_count_matches_baseline_analysis,
                "included_feature_columns_match_spec": included_feature_columns_match_spec,
                "audit_map_row_count_matches_feature_table": audit_map_row_count_matches_feature_table,
                "feature_row_index_sequential": feature_row_index_sequential,
                "included_categories_baseline_only": included_categories_baseline_only,
                "final_feature_count_matches_expected": final_feature_count_matches_expected,
                "auto_included_count_matches_expected": auto_included_count_matches_expected,
                "review_included_count_matches_expected": review_included_count_matches_expected,
                "review_excluded_count_matches_expected": review_excluded_count_matches_expected,
                "no_unexpected_review_fields": no_unexpected_review_fields,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "unit_of_analysis": "patient/case",
                "input_layer": "baseline_analysis_v1_and_baseline_profile_v1_only",
                "output_layer": "baseline_feature_set_v1",
                "fixed_review_field_decisions_json": json.dumps(
                    REVIEW_FIELD_DECISIONS, ensure_ascii=True, sort_keys=True
                ),
                "one_row_per_patient": True,
                "baseline_feature_columns_fixed_json": json_list(FINAL_FEATURE_COLUMNS),
                "no_missing_value_imputation": True,
                "no_endpoint_freeze": True,
                "no_treatment_reintegration": True,
                "no_modeling": True,
                "no_metabric": True,
                "no_raw_xml_or_ssf_parsing": True,
                "audit_map_sidecar_required": True,
                "feature_table_id_free": True,
                "followup_join_evidence_excluded": True,
                "biospecimen_sample_anchor_evidence_excluded": True,
                "audit_id_fields_excluded": True,
                "derived_prep_flags_excluded": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "missing_like_tokens_json": json_list(sorted(MISSING_LIKE_TOKENS)),
            },
            "counts": {
                "input_baseline_analysis_v1_row_count": len(workflow_inputs.baseline_rows),
                "final_row_count": len(feature_rows),
                "included_feature_count": len(included_spec_rows),
                "auto_included_feature_count": auto_included_count,
                "review_decided_included_feature_count": review_decided_included_count,
                "review_decided_excluded_feature_count": review_decided_excluded_count,
                "spec_row_count": len(spec_rows),
                "missingness_row_count": len(missingness_rows),
                "summary_row_count": len(summary_rows),
                "audit_map_row_count": len(audit_map_rows),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "baseline_analysis_v1_latest_pointer": workflow_inputs.baseline_latest_pointer,
                "baseline_analysis_v1_run_log": workflow_inputs.baseline_run_log,
                "baseline_profile_v1_latest_pointer": workflow_inputs.baseline_profile_latest_pointer,
                "baseline_profile_v1_run_log": workflow_inputs.baseline_profile_run_log,
            },
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "baseline_feature_set_v1_run_id": baseline_feature_set_v1_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_baseline_feature_set_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA baseline feature-set v1 workflow complete.")
    print(f"Baseline feature-set run ID: {run_log['baseline_feature_set_v1_run_id']}")
    print(f"Baseline profile run ID: {run_log['baseline_profile_v1_run_id']}")
    print(f"Baseline analysis run ID: {run_log['baseline_analysis_v1_run_id']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"Feature-set TSV: {run_log['outputs']['baseline_feature_set_v1_tsv']}")
    print(f"Spec TSV: {run_log['outputs']['baseline_feature_set_v1_spec_tsv']}")
    print(f"Missingness TSV: {run_log['outputs']['baseline_feature_set_v1_missingness_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['baseline_feature_set_v1_summary_tsv']}")
    print(f"Audit map TSV: {run_log['outputs']['baseline_feature_set_v1_audit_map_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
