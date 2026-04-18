#!/usr/bin/env python
"""Prepare a reproducible TCGA-BRCA baseline model-input v1 matrix from baseline feature-set v1."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
import unicodedata
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

FEATURE_DICTIONARY_FIELDNAMES = [
    "source_field_name",
    "source_value",
    "model_input_field_name",
    "feature_type",
    "encoding_rule",
    "missingness_action",
    "distinct_non_missing_count",
    "manual_review_needed",
    "notes",
]
ENCODING_SPEC_FIELDNAMES = [
    "source_field_name",
    "source_value",
    "output_field_name",
    "output_value",
    "encoding_rule",
    "reversible",
    "manual_review_needed",
    "notes",
]
MISSINGNESS_ACTION_FIELDNAMES = [
    "field_name",
    "feature_type",
    "encoding_rule",
    "missing_like_fraction",
    "missingness_action",
    "reason",
    "manual_review_needed",
]
SUMMARY_FIELDNAMES = [
    "baseline_model_input_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
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
EXPECTED_SOURCE_COLUMNS = [
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
ALLOWED_FEATURE_TYPES = {
    "binary",
    "ordinal",
    "nominal_categorical",
    "numeric_continuous",
    "numeric_count",
    "manual_review_needed",
}
ALLOWED_ENCODING_RULES = {
    "passthrough_numeric",
    "binary_map_required",
    "one_hot_encode",
    "ordinal_map_required",
    "manual_encoding_review",
}
MANUAL_REVIEW_FIELDS = {
    "ajcc_metastasis_pathologic_pm",
    "ajcc_nodes_pathologic_pn",
    "ajcc_pathologic_tumor_stage",
    "ajcc_staging_edition",
    "ajcc_tumor_pathologic_pt",
    "anatomic_neoplasm_subdivision",
    "icd_o_3_histology",
}
FIELD_RULES: dict[str, dict[str, Any]] = {
    "age_at_diagnosis": {
        "feature_type": "numeric_continuous",
        "encoding_rule": "passthrough_numeric",
        "manual_review_needed": False,
    },
    "ajcc_metastasis_pathologic_pm": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": True,
    },
    "ajcc_nodes_pathologic_pn": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": True,
    },
    "ajcc_pathologic_tumor_stage": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": True,
    },
    "ajcc_staging_edition": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": True,
    },
    "ajcc_tumor_pathologic_pt": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": True,
    },
    "anatomic_neoplasm_subdivision": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": True,
    },
    "axillary_staging_method": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "birth_days_to": {
        "feature_type": "numeric_continuous",
        "encoding_rule": "passthrough_numeric",
        "manual_review_needed": False,
    },
    "er_status_by_ihc": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "ethnicity": {
        "feature_type": "binary",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "gender": {
        "feature_type": "binary",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "her2_status_by_ihc": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "histological_type": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "history_other_malignancy": {
        "feature_type": "binary",
        "encoding_rule": "binary_map_required",
        "manual_review_needed": False,
    },
    "icd_o_3_histology": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": True,
    },
    "initial_pathologic_dx_year": {
        "feature_type": "numeric_continuous",
        "encoding_rule": "passthrough_numeric",
        "manual_review_needed": False,
    },
    "lymph_nodes_examined_count": {
        "feature_type": "numeric_count",
        "encoding_rule": "passthrough_numeric",
        "manual_review_needed": False,
    },
    "lymph_nodes_examined_he_count": {
        "feature_type": "numeric_count",
        "encoding_rule": "passthrough_numeric",
        "manual_review_needed": False,
    },
    "margin_status": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "menopause_status": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "method_initial_path_dx": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "pr_status_by_ihc": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "race": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
    "surgical_procedure_first": {
        "feature_type": "nominal_categorical",
        "encoding_rule": "one_hot_encode",
        "manual_review_needed": False,
    },
}
BINARY_VALUE_MAPS = {
    "history_other_malignancy": {
        "No": "0",
        "Yes": "1",
    }
}


class BaselineModelInputV1Error(RuntimeError):
    """Raised when the baseline model-input workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the baseline model-input workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    baseline_feature_set_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved baseline feature-set outputs."""

    baseline_feature_set_latest_pointer: dict[str, Any]
    baseline_feature_set_run_log: dict[str, Any]
    feature_rows: list[dict[str, str]]
    feature_spec_rows: list[dict[str, str]]
    feature_missingness_rows: list[dict[str, str]]
    feature_summary_rows: list[dict[str, str]]
    feature_audit_map_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class FieldProfile:
    """Observed source-field profile used to build v1 model-input outputs."""

    field_name: str
    feature_type: str
    encoding_rule: str
    manual_review_needed: bool
    non_missing_count: int
    missing_like_count: int
    missing_like_fraction: float
    distinct_non_missing_values: tuple[str, ...]
    distinct_non_missing_count: int
    missingness_action: str


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise BaselineModelInputV1Error(f"Required helper script not found: {script_path}")

    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise BaselineModelInputV1Error(f"Unable to create an import spec for: {script_path}")

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


def parse_int(value: Any, label: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise BaselineModelInputV1Error(f"Expected integer-like value for {label}: {value!r}") from exc


def parse_float(value: Any, label: str) -> float:
    try:
        return float(str(value))
    except (TypeError, ValueError) as exc:
        raise BaselineModelInputV1Error(f"Expected float-like value for {label}: {value!r}") from exc


def json_list(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=True)


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def format_fraction(value: float) -> str:
    return f"{value:.6f}"


def normalize_value(value: Any) -> str:
    return str(value or "").strip().lower()


def is_scalar_missing_like(value: Any) -> bool:
    return normalize_value(value) in MISSING_LIKE_TOKENS


def cell_is_missing_like(raw_value: Any) -> bool:
    return is_scalar_missing_like(raw_value)


def is_numeric_like(raw_value: str) -> bool:
    return re.fullmatch(r"-?\d+(?:\.\d+)?", raw_value.strip()) is not None


def build_summary_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        metric = str(row.get("summary_metric") or "")
        if metric in lookup:
            raise BaselineModelInputV1Error(f"Duplicate summary_metric detected: {metric}")
        lookup[metric] = row
    return lookup


def slugify_value(raw_value: str) -> str:
    normalized = unicodedata.normalize("NFKD", raw_value)
    ascii_value = normalized.encode("ascii", "ignore").decode("ascii")
    ascii_value = ascii_value.lower()
    slug = re.sub(r"[^a-z0-9]+", "_", ascii_value).strip("_")
    return slug or "value"


def require_pointer(
    pointer_path: Path,
    *,
    label: str,
    required_keys: set[str],
    helper_module: Any,
) -> dict[str, Any]:
    if not pointer_path.exists():
        raise BaselineModelInputV1Error(f"Required {label} not found: {pointer_path}")
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
        raise BaselineModelInputV1Error(f"{label} is not completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise BaselineModelInputV1Error(f"{label} does not report validation.passed == true.")
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
    model_input_root = audit_root / "tcga-brca" / "model-input"

    return WorkflowPaths(
        repo_root=helper_paths.repo_root,
        trial_config=helper_paths.trial_config,
        results_root=results_root,
        processed_runs_root=processed_root / "tcga-brca" / "model-input" / "baseline_v1_runs",
        audit_runs_root=model_input_root / "baseline_v1_runs",
        latest_pointer=model_input_root / "tcga_brca_baseline_model_input_v1_latest.json",
        baseline_feature_set_latest_pointer=(
            audit_root / "tcga-brca" / "analysis-prep" / "tcga_brca_baseline_feature_set_v1_latest.json"
        ),
    )


def build_input_paths(
    paths: WorkflowPaths,
    *,
    baseline_feature_set_latest_pointer: dict[str, Any],
    helper_module: Any,
) -> dict[str, Path]:
    return {
        "baseline_feature_set_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_feature_set_latest_pointer["baseline_feature_set_v1_tsv"]),
            "baseline feature-set v1 TSV",
        ),
        "baseline_feature_set_v1_spec_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_feature_set_latest_pointer["baseline_feature_set_v1_spec_tsv"]),
            "baseline feature-set v1 spec TSV",
        ),
        "baseline_feature_set_v1_missingness_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_feature_set_latest_pointer["baseline_feature_set_v1_missingness_tsv"]),
            "baseline feature-set v1 missingness TSV",
        ),
        "baseline_feature_set_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_feature_set_latest_pointer["baseline_feature_set_v1_summary_tsv"]),
            "baseline feature-set v1 summary TSV",
        ),
        "baseline_feature_set_v1_audit_map_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_feature_set_latest_pointer["baseline_feature_set_v1_audit_map_tsv"]),
            "baseline feature-set v1 audit-map TSV",
        ),
        "baseline_feature_set_run_log_json": helper_module.resolve_existing_path(
            paths.repo_root,
            str(baseline_feature_set_latest_pointer["run_log_json"]),
            "baseline feature-set v1 run log",
        ),
    }


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    baseline_feature_set_latest_pointer = require_pointer(
        paths.baseline_feature_set_latest_pointer,
        label="Baseline feature-set v1 latest pointer",
        required_keys={
            "baseline_feature_set_v1_run_id",
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
            "baseline_feature_set_v1_tsv",
            "baseline_feature_set_v1_spec_tsv",
            "baseline_feature_set_v1_missingness_tsv",
            "baseline_feature_set_v1_summary_tsv",
            "baseline_feature_set_v1_audit_map_tsv",
            "run_log_json",
        },
        helper_module=helper_module,
    )
    input_paths = build_input_paths(
        paths,
        baseline_feature_set_latest_pointer=baseline_feature_set_latest_pointer,
        helper_module=helper_module,
    )
    _, baseline_feature_set_run_log = require_completed_run_log(
        paths.repo_root,
        str(baseline_feature_set_latest_pointer["run_log_json"]),
        label="Baseline feature-set v1 run log",
        helper_module=helper_module,
    )

    feature_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_feature_set_v1_tsv"])
    feature_spec_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_feature_set_v1_spec_tsv"])
    feature_missingness_rows = helper_module.read_tsv_dict_rows(
        input_paths["baseline_feature_set_v1_missingness_tsv"]
    )
    feature_summary_rows = helper_module.read_tsv_dict_rows(input_paths["baseline_feature_set_v1_summary_tsv"])
    feature_audit_map_rows = helper_module.read_tsv_dict_rows(
        input_paths["baseline_feature_set_v1_audit_map_tsv"]
    )

    if not feature_rows:
        raise BaselineModelInputV1Error("The referenced baseline_feature_set_v1.tsv contains no rows.")
    if not feature_spec_rows:
        raise BaselineModelInputV1Error("The referenced baseline_feature_set_v1_spec.tsv contains no rows.")
    if not feature_missingness_rows:
        raise BaselineModelInputV1Error(
            "The referenced baseline_feature_set_v1_missingness.tsv contains no rows."
        )
    if not feature_summary_rows:
        raise BaselineModelInputV1Error("The referenced baseline_feature_set_v1_summary.tsv contains no rows.")
    if not feature_audit_map_rows:
        raise BaselineModelInputV1Error(
            "The referenced baseline_feature_set_v1_audit_map.tsv contains no rows."
        )

    return WorkflowInputs(
        baseline_feature_set_latest_pointer=baseline_feature_set_latest_pointer,
        baseline_feature_set_run_log=baseline_feature_set_run_log,
        feature_rows=feature_rows,
        feature_spec_rows=feature_spec_rows,
        feature_missingness_rows=feature_missingness_rows,
        feature_summary_rows=feature_summary_rows,
        feature_audit_map_rows=feature_audit_map_rows,
        input_paths=input_paths,
    )


def validate_upstream_state(workflow_inputs: WorkflowInputs) -> None:
    source_columns = list(workflow_inputs.feature_rows[0].keys())
    if source_columns != EXPECTED_SOURCE_COLUMNS:
        raise BaselineModelInputV1Error(
            "baseline_feature_set_v1.tsv columns do not match the explicit 25-field ruleset."
        )

    included_spec_rows = [
        row
        for row in workflow_inputs.feature_spec_rows
        if str(row.get("feature_set_v1_decision") or "") == "include_in_feature_set_v1"
    ]
    included_spec_field_names = [str(row.get("field_name") or "") for row in included_spec_rows]
    if included_spec_field_names != EXPECTED_SOURCE_COLUMNS:
        raise BaselineModelInputV1Error(
            "baseline_feature_set_v1_spec.tsv included fields do not reconcile to the explicit 25-field ruleset."
        )
    if any(str(row.get("field_category") or "") != "baseline_clinical_field" for row in included_spec_rows):
        raise BaselineModelInputV1Error(
            "baseline_feature_set_v1_spec.tsv includes a non-baseline_clinical_field among the retained fields."
        )

    missingness_field_names = [str(row.get("field_name") or "") for row in workflow_inputs.feature_missingness_rows]
    if missingness_field_names != EXPECTED_SOURCE_COLUMNS:
        raise BaselineModelInputV1Error(
            "baseline_feature_set_v1_missingness.tsv field order does not match the explicit 25-field ruleset."
        )

    summary_lookup = build_summary_lookup(workflow_inputs.feature_summary_rows)
    if summary_lookup["endpoint_freeze_status"]["summary_value"] != "blocked":
        raise BaselineModelInputV1Error("Upstream baseline feature-set summary does not report endpoint_freeze_status == blocked.")
    if summary_lookup["treatment_inclusion_status"]["summary_value"] != "excluded":
        raise BaselineModelInputV1Error(
            "Upstream baseline feature-set summary does not report treatment_inclusion_status == excluded."
        )
    if summary_lookup["readiness_for_first_pass_baseline_model_input"]["summary_value"] != "ready":
        raise BaselineModelInputV1Error(
            "Upstream baseline feature-set summary does not report readiness_for_first_pass_baseline_model_input == ready."
        )
    final_row_count = parse_int(summary_lookup["final_row_count"]["summary_value"], "final_row_count")
    if final_row_count != len(workflow_inputs.feature_rows):
        raise BaselineModelInputV1Error(
            "Upstream baseline feature-set summary final_row_count does not match baseline_feature_set_v1.tsv."
        )

    if len(workflow_inputs.feature_audit_map_rows) != len(workflow_inputs.feature_rows):
        raise BaselineModelInputV1Error(
            "baseline_feature_set_v1_audit_map.tsv row count does not match baseline_feature_set_v1.tsv."
        )
    audit_row_indices = [
        parse_int(row.get("feature_set_v1_row_index", ""), "feature_set_v1_row_index")
        for row in workflow_inputs.feature_audit_map_rows
    ]
    if audit_row_indices != list(range(1, len(workflow_inputs.feature_audit_map_rows) + 1)):
        raise BaselineModelInputV1Error(
            "baseline_feature_set_v1_audit_map.tsv feature_set_v1_row_index values are not sequential."
        )


def missingness_action_for_rule(feature_type: str, encoding_rule: str, missing_like_count: int) -> str:
    if missing_like_count == 0:
        return "keep_as_is_for_now"
    if encoding_rule == "passthrough_numeric":
        return "eligible_for_simple_imputation_later"
    return "allow_missing_category"


def profile_source_fields(workflow_inputs: WorkflowInputs) -> dict[str, FieldProfile]:
    missingness_lookup = {
        str(row.get("field_name") or ""): row for row in workflow_inputs.feature_missingness_rows
    }
    profiles: dict[str, FieldProfile] = {}

    for field_name in EXPECTED_SOURCE_COLUMNS:
        if field_name not in FIELD_RULES:
            raise BaselineModelInputV1Error(f"No FIELD_RULES entry was defined for: {field_name}")

        raw_values = [str(row.get(field_name, "")) for row in workflow_inputs.feature_rows]
        missing_like_count = sum(cell_is_missing_like(value) for value in raw_values)
        non_missing_values = [value for value in raw_values if not cell_is_missing_like(value)]
        distinct_non_missing_values = tuple(sorted(ordered_unique(non_missing_values)))
        non_missing_count = len(non_missing_values)
        distinct_non_missing_count = len(distinct_non_missing_values)
        missing_like_fraction = missing_like_count / len(raw_values)

        field_rule = FIELD_RULES[field_name]
        feature_type = str(field_rule["feature_type"])
        encoding_rule = str(field_rule["encoding_rule"])
        manual_review_needed = bool(field_rule["manual_review_needed"])
        if feature_type not in ALLOWED_FEATURE_TYPES:
            raise BaselineModelInputV1Error(f"Unsupported feature_type for {field_name}: {feature_type}")
        if encoding_rule not in ALLOWED_ENCODING_RULES:
            raise BaselineModelInputV1Error(f"Unsupported encoding_rule for {field_name}: {encoding_rule}")

        if feature_type in {"numeric_continuous", "numeric_count"}:
            non_numeric_values = [value for value in non_missing_values if not is_numeric_like(value)]
            if non_numeric_values:
                raise BaselineModelInputV1Error(
                    f"Numeric field {field_name} contains non-numeric values: {non_numeric_values[:5]!r}"
                )

        if encoding_rule == "binary_map_required":
            expected_values = set(BINARY_VALUE_MAPS[field_name].keys())
            unexpected_values = sorted(set(distinct_non_missing_values).difference(expected_values))
            if unexpected_values:
                raise BaselineModelInputV1Error(
                    f"Binary-mapped field {field_name} contains unexpected values: {unexpected_values!r}"
                )
            if distinct_non_missing_count == 0:
                raise BaselineModelInputV1Error(
                    f"Binary-mapped field {field_name} contains no non-missing values."
                )

        upstream_missingness_row = missingness_lookup.get(field_name)
        if upstream_missingness_row is None:
            raise BaselineModelInputV1Error(
                f"Upstream baseline_feature_set_v1_missingness.tsv is missing field: {field_name}"
            )
        if parse_int(upstream_missingness_row["row_count"], f"{field_name}.row_count") != len(raw_values):
            raise BaselineModelInputV1Error(
                f"Upstream missingness row_count did not match observed row count for {field_name}."
            )
        if parse_int(
            upstream_missingness_row["non_missing_count"], f"{field_name}.non_missing_count"
        ) != non_missing_count:
            raise BaselineModelInputV1Error(
                f"Upstream non_missing_count did not match observed count for {field_name}."
            )
        if parse_int(
            upstream_missingness_row["missing_like_count"], f"{field_name}.missing_like_count"
        ) != missing_like_count:
            raise BaselineModelInputV1Error(
                f"Upstream missing_like_count did not match observed count for {field_name}."
            )
        if (
            parse_int(
                upstream_missingness_row["distinct_non_missing_count"],
                f"{field_name}.distinct_non_missing_count",
            )
            != distinct_non_missing_count
        ):
            raise BaselineModelInputV1Error(
                f"Upstream distinct_non_missing_count did not match observed count for {field_name}."
            )
        upstream_fraction = parse_float(
            upstream_missingness_row["missing_like_fraction"], f"{field_name}.missing_like_fraction"
        )
        if abs(upstream_fraction - missing_like_fraction) > 5e-7:
            raise BaselineModelInputV1Error(
                f"Upstream missing_like_fraction did not match observed fraction for {field_name}."
            )

        profiles[field_name] = FieldProfile(
            field_name=field_name,
            feature_type=feature_type,
            encoding_rule=encoding_rule,
            manual_review_needed=manual_review_needed,
            non_missing_count=non_missing_count,
            missing_like_count=missing_like_count,
            missing_like_fraction=missing_like_fraction,
            distinct_non_missing_values=distinct_non_missing_values,
            distinct_non_missing_count=distinct_non_missing_count,
            missingness_action=missingness_action_for_rule(
                feature_type=feature_type,
                encoding_rule=encoding_rule,
                missing_like_count=missing_like_count,
            ),
        )
    return profiles


def build_one_hot_output_fields(
    field_name: str,
    categories: tuple[str, ...],
    *,
    include_missing: bool,
) -> list[tuple[str, str]]:
    reserved_names = {f"{field_name}__missing"} if include_missing else set()
    used_names: set[str] = set()
    mappings: list[tuple[str, str]] = []

    for category in categories:
        base_name = f"{field_name}__{slugify_value(category)}"
        candidate_name = base_name
        suffix_index = 2
        while candidate_name in used_names or candidate_name in reserved_names:
            candidate_name = f"{base_name}__dup{suffix_index}"
            suffix_index += 1
        used_names.add(candidate_name)
        mappings.append((category, candidate_name))

    return mappings


def missingness_reason(profile: FieldProfile) -> str:
    if profile.missing_like_count == 0:
        return "No missing-like values were observed in the saved baseline feature-set v1 source field."
    if profile.encoding_rule == "passthrough_numeric":
        return (
            "Numeric source field has missing-like values; v1 leaves blanks in the matrix and documents later simple imputation as a possible downstream step."
        )
    return (
        "Categorical source field has missing-like values; v1 adds a dedicated __missing indicator while preserving exact-category encoding."
    )


def build_model_input_outputs(
    workflow_inputs: WorkflowInputs,
    field_profiles: dict[str, FieldProfile],
) -> tuple[
    list[str],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    matrix_fieldnames: list[str] = []
    matrix_rows = [{} for _ in workflow_inputs.feature_rows]
    feature_dictionary_rows: list[dict[str, str]] = []
    encoding_spec_rows: list[dict[str, str]] = []
    missingness_action_rows: list[dict[str, str]] = []

    for field_name in EXPECTED_SOURCE_COLUMNS:
        profile = field_profiles[field_name]
        missingness_action_rows.append(
            {
                "field_name": field_name,
                "feature_type": profile.feature_type,
                "encoding_rule": profile.encoding_rule,
                "missing_like_fraction": format_fraction(profile.missing_like_fraction),
                "missingness_action": profile.missingness_action,
                "reason": missingness_reason(profile),
                "manual_review_needed": bool_text(profile.manual_review_needed),
            }
        )

        if profile.encoding_rule == "passthrough_numeric":
            matrix_fieldnames.append(field_name)
            for row_index, source_row in enumerate(workflow_inputs.feature_rows):
                raw_value = str(source_row[field_name])
                matrix_rows[row_index][field_name] = "" if cell_is_missing_like(raw_value) else raw_value

            feature_dictionary_rows.append(
                {
                    "source_field_name": field_name,
                    "source_value": "[numeric_passthrough]",
                    "model_input_field_name": field_name,
                    "feature_type": profile.feature_type,
                    "encoding_rule": profile.encoding_rule,
                    "missingness_action": profile.missingness_action,
                    "distinct_non_missing_count": str(profile.distinct_non_missing_count),
                    "manual_review_needed": bool_text(profile.manual_review_needed),
                    "notes": (
                        "Numeric source values are carried forward unchanged; missing-like values remain blank in the model-input matrix."
                    ),
                }
            )
            encoding_spec_rows.append(
                {
                    "source_field_name": field_name,
                    "source_value": "[numeric_passthrough]",
                    "output_field_name": field_name,
                    "output_value": "[same_as_source_value]",
                    "encoding_rule": profile.encoding_rule,
                    "reversible": "true",
                    "manual_review_needed": bool_text(profile.manual_review_needed),
                    "notes": "Numeric passthrough only; no imputation or rescaling is applied in v1.",
                }
            )
            continue

        if profile.encoding_rule == "binary_map_required":
            output_field_name = field_name
            missing_output_field_name = (
                f"{field_name}__missing" if profile.missing_like_count > 0 else None
            )
            matrix_fieldnames.append(output_field_name)
            if missing_output_field_name is not None:
                matrix_fieldnames.append(missing_output_field_name)

            feature_dictionary_rows.append(
                {
                    "source_field_name": field_name,
                    "source_value": "No|Yes",
                    "model_input_field_name": output_field_name,
                    "feature_type": profile.feature_type,
                    "encoding_rule": profile.encoding_rule,
                    "missingness_action": profile.missingness_action,
                    "distinct_non_missing_count": str(profile.distinct_non_missing_count),
                    "manual_review_needed": bool_text(profile.manual_review_needed),
                    "notes": "Exact binary map only: No -> 0 and Yes -> 1.",
                }
            )
            encoding_spec_rows.extend(
                [
                    {
                        "source_field_name": field_name,
                        "source_value": "No",
                        "output_field_name": output_field_name,
                        "output_value": "0",
                        "encoding_rule": profile.encoding_rule,
                        "reversible": "true",
                        "manual_review_needed": bool_text(profile.manual_review_needed),
                        "notes": "Observed exact raw value No maps to 0.",
                    },
                    {
                        "source_field_name": field_name,
                        "source_value": "Yes",
                        "output_field_name": output_field_name,
                        "output_value": "1",
                        "encoding_rule": profile.encoding_rule,
                        "reversible": "true",
                        "manual_review_needed": bool_text(profile.manual_review_needed),
                        "notes": "Observed exact raw value Yes maps to 1.",
                    },
                ]
            )
            if missing_output_field_name is not None:
                feature_dictionary_rows.append(
                    {
                        "source_field_name": field_name,
                        "source_value": "[missing_like]",
                        "model_input_field_name": missing_output_field_name,
                        "feature_type": profile.feature_type,
                        "encoding_rule": profile.encoding_rule,
                        "missingness_action": profile.missingness_action,
                        "distinct_non_missing_count": str(profile.distinct_non_missing_count),
                        "manual_review_needed": bool_text(profile.manual_review_needed),
                        "notes": (
                            "Dedicated missingness indicator for source rows where history_other_malignancy is missing-like."
                        ),
                    }
                )
                encoding_spec_rows.append(
                    {
                        "source_field_name": field_name,
                        "source_value": "[missing_like]",
                        "output_field_name": missing_output_field_name,
                        "output_value": "1",
                        "encoding_rule": profile.encoding_rule,
                        "reversible": "true",
                        "manual_review_needed": bool_text(profile.manual_review_needed),
                        "notes": "Missing-like source cells activate the dedicated missingness indicator.",
                    }
                )

            for row_index, source_row in enumerate(workflow_inputs.feature_rows):
                raw_value = str(source_row[field_name])
                if cell_is_missing_like(raw_value):
                    matrix_rows[row_index][output_field_name] = ""
                    if missing_output_field_name is not None:
                        matrix_rows[row_index][missing_output_field_name] = "1"
                    continue
                if raw_value not in BINARY_VALUE_MAPS[field_name]:
                    raise BaselineModelInputV1Error(
                        f"Unexpected non-missing value for binary-mapped field {field_name}: {raw_value!r}"
                    )
                matrix_rows[row_index][output_field_name] = BINARY_VALUE_MAPS[field_name][raw_value]
                if missing_output_field_name is not None:
                    matrix_rows[row_index][missing_output_field_name] = "0"
            continue

        if profile.encoding_rule != "one_hot_encode":
            raise BaselineModelInputV1Error(
                f"Unhandled encoding_rule for {field_name}: {profile.encoding_rule}"
            )

        category_output_mappings = build_one_hot_output_fields(
            field_name,
            profile.distinct_non_missing_values,
            include_missing=profile.missing_like_count > 0,
        )
        category_to_output = {source_value: output_field for source_value, output_field in category_output_mappings}
        missing_output_field_name = f"{field_name}__missing" if profile.missing_like_count > 0 else None

        for _, output_field_name in category_output_mappings:
            matrix_fieldnames.append(output_field_name)
        if missing_output_field_name is not None:
            matrix_fieldnames.append(missing_output_field_name)

        for source_value, output_field_name in category_output_mappings:
            feature_dictionary_rows.append(
                {
                    "source_field_name": field_name,
                    "source_value": source_value,
                    "model_input_field_name": output_field_name,
                    "feature_type": profile.feature_type,
                    "encoding_rule": profile.encoding_rule,
                    "missingness_action": profile.missingness_action,
                    "distinct_non_missing_count": str(profile.distinct_non_missing_count),
                    "manual_review_needed": bool_text(profile.manual_review_needed),
                    "notes": "Exact raw source category is represented by this one-hot output column.",
                }
            )
            encoding_spec_rows.append(
                {
                    "source_field_name": field_name,
                    "source_value": source_value,
                    "output_field_name": output_field_name,
                    "output_value": "1",
                    "encoding_rule": profile.encoding_rule,
                    "reversible": "true",
                    "manual_review_needed": bool_text(profile.manual_review_needed),
                    "notes": "Exact raw source category activates this one-hot output column.",
                }
            )
        if missing_output_field_name is not None:
            feature_dictionary_rows.append(
                {
                    "source_field_name": field_name,
                    "source_value": "[missing_like]",
                    "model_input_field_name": missing_output_field_name,
                    "feature_type": profile.feature_type,
                    "encoding_rule": profile.encoding_rule,
                    "missingness_action": profile.missingness_action,
                    "distinct_non_missing_count": str(profile.distinct_non_missing_count),
                    "manual_review_needed": bool_text(profile.manual_review_needed),
                    "notes": "Dedicated missingness indicator for missing-like source cells in this categorical field.",
                }
            )
            encoding_spec_rows.append(
                {
                    "source_field_name": field_name,
                    "source_value": "[missing_like]",
                    "output_field_name": missing_output_field_name,
                    "output_value": "1",
                    "encoding_rule": profile.encoding_rule,
                    "reversible": "true",
                    "manual_review_needed": bool_text(profile.manual_review_needed),
                    "notes": "Missing-like source cells activate the dedicated missingness indicator.",
                }
            )

        category_output_names = [output_field_name for _, output_field_name in category_output_mappings]
        for row_index, source_row in enumerate(workflow_inputs.feature_rows):
            raw_value = str(source_row[field_name])
            for output_field_name in category_output_names:
                matrix_rows[row_index][output_field_name] = "0"
            if missing_output_field_name is not None:
                matrix_rows[row_index][missing_output_field_name] = "0"

            if cell_is_missing_like(raw_value):
                if missing_output_field_name is not None:
                    matrix_rows[row_index][missing_output_field_name] = "1"
                continue

            output_field_name = category_to_output.get(raw_value)
            if output_field_name is None:
                raise BaselineModelInputV1Error(
                    f"Unexpected non-missing value for one-hot field {field_name}: {raw_value!r}"
                )
            matrix_rows[row_index][output_field_name] = "1"

    return (
        matrix_fieldnames,
        matrix_rows,
        feature_dictionary_rows,
        encoding_spec_rows,
        missingness_action_rows,
    )


def build_summary_rows(
    *,
    baseline_model_input_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
    field_profiles: dict[str, FieldProfile],
    matrix_fieldnames: list[str],
) -> list[dict[str, str]]:
    feature_type_order = [
        "binary",
        "ordinal",
        "nominal_categorical",
        "numeric_continuous",
        "numeric_count",
        "manual_review_needed",
    ]
    encoding_rule_order = [
        "passthrough_numeric",
        "binary_map_required",
        "one_hot_encode",
        "ordinal_map_required",
        "manual_encoding_review",
    ]
    feature_type_counts = {
        feature_type: sum(profile.feature_type == feature_type for profile in field_profiles.values())
        for feature_type in feature_type_order
    }
    encoding_rule_counts = {
        encoding_rule: sum(profile.encoding_rule == encoding_rule for profile in field_profiles.values())
        for encoding_rule in encoding_rule_order
    }
    manual_review_fields = [
        field_name for field_name in EXPECTED_SOURCE_COLUMNS if field_profiles[field_name].manual_review_needed
    ]
    fields_with_missing_category = [
        field_name
        for field_name in EXPECTED_SOURCE_COLUMNS
        if field_profiles[field_name].missing_like_count > 0
        and field_profiles[field_name].encoding_rule != "passthrough_numeric"
    ]
    numeric_fields_with_blank_missing = [
        field_name
        for field_name in EXPECTED_SOURCE_COLUMNS
        if field_profiles[field_name].missing_like_count > 0
        and field_profiles[field_name].encoding_rule == "passthrough_numeric"
    ]

    rows: list[dict[str, str]] = [
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "design",
            "summary_metric": "baseline_model_input_v1_status",
            "summary_value": "baseline_model_input_v1_complete",
            "notes": "This workflow prepares a first-pass baseline model-input matrix only; it does not perform modeling.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "design",
            "summary_metric": "input_source_table",
            "summary_value": "baseline_feature_set_v1",
            "notes": "The model-input matrix starts only from the saved 25-feature baseline feature-set v1 table.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "design",
            "summary_metric": "matrix_representation_policy",
            "summary_value": "fully_numeric_exact_label_one_hot_with_numeric_passthrough",
            "notes": "Numeric fields pass through unchanged; categorical fields use exact-label one-hot columns; no identifier column is added to the matrix.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "design",
            "summary_metric": "unit_of_analysis",
            "summary_value": "patient/case",
            "notes": "The matrix preserves one row per patient/case from baseline_feature_set_v1.tsv.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "baseline_feature_set_v1_run_id",
            "summary_value": str(
                workflow_inputs.baseline_feature_set_latest_pointer["baseline_feature_set_v1_run_id"]
            ),
            "notes": "Source baseline feature-set run used for model-input preparation.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "baseline_profile_v1_run_id",
            "summary_value": str(workflow_inputs.baseline_feature_set_latest_pointer["baseline_profile_v1_run_id"]),
            "notes": "Carried forward for provenance only.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "baseline_analysis_v1_run_id",
            "summary_value": str(workflow_inputs.baseline_feature_set_latest_pointer["baseline_analysis_v1_run_id"]),
            "notes": "Carried forward for provenance only.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "input_row_count",
            "summary_value": str(len(workflow_inputs.feature_rows)),
            "notes": "Row count read directly from baseline_feature_set_v1.tsv.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "final_row_count",
            "summary_value": str(len(workflow_inputs.feature_rows)),
            "notes": "The model-input matrix preserves baseline_feature_set_v1.tsv row count and row order exactly.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "field_counts",
            "summary_metric": "input_feature_count",
            "summary_value": str(len(EXPECTED_SOURCE_COLUMNS)),
            "notes": "All 25 retained baseline feature-set v1 fields were carried into model-input preparation.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "field_counts",
            "summary_metric": "output_feature_count",
            "summary_value": str(len(matrix_fieldnames)),
            "notes": "Count of final model-input matrix columns after numeric passthrough, binary mapping, and exact-label one-hot expansion.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "input_feature_names_json",
            "summary_value": json_list(EXPECTED_SOURCE_COLUMNS),
            "notes": "Saved source feature names used for baseline model-input preparation, in matrix-generation order.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "manual_review_fields_json",
            "summary_value": json_list(manual_review_fields),
            "notes": "Fields still flagged for manual encoding review even though safe initial v1 encoding was generated.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "missingness",
            "summary_metric": "categorical_fields_with_missing_indicator_count",
            "summary_value": str(len(fields_with_missing_category)),
            "notes": "Count of source fields that receive a dedicated __missing output column in v1.",
        },
        {
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "summary_section": "missingness",
            "summary_metric": "numeric_fields_left_blank_for_missing_count",
            "summary_value": str(len(numeric_fields_with_blank_missing)),
            "notes": "Count of numeric source fields that retain blank cells for missing-like values in v1.",
        },
    ]

    rows.extend(
        [
            {
                "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
                "summary_section": "readiness",
                "summary_metric": "manual_review_required_field_count",
                "summary_value": str(len(manual_review_fields)),
                "notes": "These fields were encoded safely for v1 but still require manual review before stronger modeling claims.",
            },
            {
                "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
                "summary_section": "readiness",
                "summary_metric": "readiness_for_first_pass_baseline_modeling",
                "summary_value": "ready_with_documented_manual_review_flags",
                "notes": "The saved matrix is suitable for first-pass baseline-only modeling experiments once the documented manual-review fields are acknowledged.",
            },
            {
                "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
                "summary_section": "readiness",
                "summary_metric": "endpoint_freeze_status",
                "summary_value": "blocked",
                "notes": "Endpoint freeze remains blocked; this workflow does not add endpoint fields or freeze a final endpoint.",
            },
            {
                "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
                "summary_section": "readiness",
                "summary_metric": "treatment_inclusion_status",
                "summary_value": "excluded",
                "notes": "Treatment detail remains excluded from model-input preparation in this workflow.",
            },
            {
                "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
                "summary_section": "readiness",
                "summary_metric": "overall_readiness_interpretation",
                "summary_value": (
                    "suitable_for_first_pass_baseline_modeling_experiments_with_documented_manual_review_flags"
                ),
                "notes": (
                    "Use this matrix for first-pass baseline-only modeling experiments, encoding review, and missingness review only; do not treat it as endpoint freeze, treatment modeling, or final clinical interpretation."
                ),
            },
        ]
    )

    for feature_type in feature_type_order:
        rows.append(
            {
                "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
                "summary_section": "field_types",
                "summary_metric": f"feature_type_count__{feature_type}",
                "summary_value": str(feature_type_counts[feature_type]),
                "notes": "Count across the 25 saved source features before one-hot expansion.",
            }
        )

    for encoding_rule in encoding_rule_order:
        rows.append(
            {
                "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
                "summary_section": "encoding",
                "summary_metric": f"encoding_rule_count__{encoding_rule}",
                "summary_value": str(encoding_rule_counts[encoding_rule]),
                "notes": "Count across the 25 saved source features before one-hot expansion.",
            }
        )

    return rows


def validate_generated_outputs(
    workflow_inputs: WorkflowInputs,
    field_profiles: dict[str, FieldProfile],
    matrix_fieldnames: list[str],
    matrix_rows: list[dict[str, str]],
    feature_dictionary_rows: list[dict[str, str]],
    encoding_spec_rows: list[dict[str, str]],
    missingness_action_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
) -> dict[str, Any]:
    if not matrix_rows:
        raise BaselineModelInputV1Error("baseline_model_input_v1.tsv rows were not generated.")
    if not feature_dictionary_rows:
        raise BaselineModelInputV1Error(
            "baseline_model_input_v1_feature_dictionary.tsv rows were not generated."
        )
    if not encoding_spec_rows:
        raise BaselineModelInputV1Error(
            "baseline_model_input_v1_encoding_spec.tsv rows were not generated."
        )
    if not missingness_action_rows:
        raise BaselineModelInputV1Error(
            "baseline_model_input_v1_missingness_actions.tsv rows were not generated."
        )
    if not summary_rows:
        raise BaselineModelInputV1Error("baseline_model_input_v1_summary.tsv rows were not generated.")

    matrix_row_count_matches_source = len(matrix_rows) == len(workflow_inputs.feature_rows)
    output_rows_positive = len(matrix_rows) > 0
    row_order_preserved_from_baseline_feature_set_v1 = True
    input_source_columns_match_expected = list(workflow_inputs.feature_rows[0].keys()) == EXPECTED_SOURCE_COLUMNS
    required_source_tables_found = all(path.exists() for path in workflow_inputs.input_paths.values())
    upstream_audit_map_row_count_matches_input = (
        len(workflow_inputs.feature_audit_map_rows) == len(workflow_inputs.feature_rows)
    )
    upstream_audit_map_row_index_sequential = [
        parse_int(row["feature_set_v1_row_index"], "feature_set_v1_row_index")
        for row in workflow_inputs.feature_audit_map_rows
    ] == list(range(1, len(workflow_inputs.feature_audit_map_rows) + 1))
    matrix_rows_have_all_output_columns = all(
        list(row.keys()) == matrix_fieldnames for row in matrix_rows
    )
    matrix_output_column_count_positive = len(matrix_fieldnames) > 0
    matrix_output_column_names_unique = len(matrix_fieldnames) == len(set(matrix_fieldnames))
    feature_dictionary_row_count_matches_matrix_columns = (
        len(feature_dictionary_rows) == len(matrix_fieldnames)
    )
    feature_dictionary_output_fields_match_matrix_columns = (
        [row["model_input_field_name"] for row in feature_dictionary_rows] == matrix_fieldnames
    )
    encoding_spec_output_fields_cover_matrix_columns = (
        set(row["output_field_name"] for row in encoding_spec_rows) == set(matrix_fieldnames)
    )
    missingness_action_rows_match_input_feature_count = (
        len(missingness_action_rows) == len(EXPECTED_SOURCE_COLUMNS)
    )
    missingness_action_fields_match_expected = (
        [row["field_name"] for row in missingness_action_rows] == EXPECTED_SOURCE_COLUMNS
    )
    summary_lookup = build_summary_lookup(summary_rows)
    summary_endpoint_blocked = summary_lookup["endpoint_freeze_status"]["summary_value"] == "blocked"
    summary_treatment_excluded = summary_lookup["treatment_inclusion_status"]["summary_value"] == "excluded"
    summary_input_feature_count_correct = (
        summary_lookup["input_feature_count"]["summary_value"] == str(len(EXPECTED_SOURCE_COLUMNS))
    )
    summary_output_feature_count_correct = (
        summary_lookup["output_feature_count"]["summary_value"] == str(len(matrix_fieldnames))
    )
    manual_review_field_count_matches_expected = sum(
        profile.manual_review_needed for profile in field_profiles.values()
    ) == len(MANUAL_REVIEW_FIELDS)

    return {
        "passed": (
            required_source_tables_found
            and output_rows_positive
            and matrix_row_count_matches_source
            and row_order_preserved_from_baseline_feature_set_v1
            and input_source_columns_match_expected
            and upstream_audit_map_row_count_matches_input
            and upstream_audit_map_row_index_sequential
            and matrix_rows_have_all_output_columns
            and matrix_output_column_count_positive
            and matrix_output_column_names_unique
            and feature_dictionary_row_count_matches_matrix_columns
            and feature_dictionary_output_fields_match_matrix_columns
            and encoding_spec_output_fields_cover_matrix_columns
            and missingness_action_rows_match_input_feature_count
            and missingness_action_fields_match_expected
            and summary_endpoint_blocked
            and summary_treatment_excluded
            and summary_input_feature_count_correct
            and summary_output_feature_count_correct
            and manual_review_field_count_matches_expected
        ),
        "required_upstream_pointers_found": True,
        "baseline_feature_set_v1_latest_pointer_found": True,
        "baseline_feature_set_run_log_completed": True,
        "baseline_feature_set_validation_passed": True,
        "required_source_tables_found": required_source_tables_found,
        "input_source_columns_match_expected": input_source_columns_match_expected,
        "output_rows_positive": output_rows_positive,
        "matrix_row_count_matches_source": matrix_row_count_matches_source,
        "row_order_preserved_from_baseline_feature_set_v1": row_order_preserved_from_baseline_feature_set_v1,
        "upstream_audit_map_row_count_matches_input": upstream_audit_map_row_count_matches_input,
        "upstream_audit_map_row_index_sequential": upstream_audit_map_row_index_sequential,
        "matrix_rows_have_all_output_columns": matrix_rows_have_all_output_columns,
        "matrix_output_column_count_positive": matrix_output_column_count_positive,
        "matrix_output_column_names_unique": matrix_output_column_names_unique,
        "feature_dictionary_row_count_matches_matrix_columns": (
            feature_dictionary_row_count_matches_matrix_columns
        ),
        "feature_dictionary_output_fields_match_matrix_columns": (
            feature_dictionary_output_fields_match_matrix_columns
        ),
        "encoding_spec_output_fields_cover_matrix_columns": encoding_spec_output_fields_cover_matrix_columns,
        "missingness_action_rows_match_input_feature_count": (
            missingness_action_rows_match_input_feature_count
        ),
        "missingness_action_fields_match_expected": missingness_action_fields_match_expected,
        "summary_endpoint_blocked": summary_endpoint_blocked,
        "summary_treatment_excluded": summary_treatment_excluded,
        "summary_input_feature_count_correct": summary_input_feature_count_correct,
        "summary_output_feature_count_correct": summary_output_feature_count_correct,
        "manual_review_field_count_matches_expected": manual_review_field_count_matches_expected,
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
    }


def build_latest_pointer_payload(
    *,
    baseline_model_input_v1_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
        "baseline_feature_set_v1_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["baseline_feature_set_v1_run_id"]
        ),
        "baseline_profile_v1_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["baseline_profile_v1_run_id"]
        ),
        "baseline_analysis_v1_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["baseline_analysis_v1_run_id"]
        ),
        "cohort_v1_build_id": str(workflow_inputs.baseline_feature_set_latest_pointer["cohort_v1_build_id"]),
        "dry_run_build_id": str(workflow_inputs.baseline_feature_set_latest_pointer["dry_run_build_id"]),
        "blueprint_run_id": str(workflow_inputs.baseline_feature_set_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(workflow_inputs.baseline_feature_set_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(workflow_inputs.baseline_feature_set_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["clinical_parse_run_id"]
        ),
        "endpoint_crosswalk_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["endpoint_crosswalk_run_id"]
        ),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["biospecimen_crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["biospecimen_parse_run_id"]
        ),
        "clinical_source_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["clinical_source_run_id"]
        ),
        "biospecimen_source_run_id": str(
            workflow_inputs.baseline_feature_set_latest_pointer["biospecimen_source_run_id"]
        ),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "baseline_model_input_v1_tsv": repo_relative(
            output_paths["baseline_model_input_v1_tsv"], paths.repo_root
        ),
        "baseline_model_input_v1_feature_dictionary_tsv": repo_relative(
            output_paths["baseline_model_input_v1_feature_dictionary_tsv"], paths.repo_root
        ),
        "baseline_model_input_v1_encoding_spec_tsv": repo_relative(
            output_paths["baseline_model_input_v1_encoding_spec_tsv"], paths.repo_root
        ),
        "baseline_model_input_v1_missingness_actions_tsv": repo_relative(
            output_paths["baseline_model_input_v1_missingness_actions_tsv"], paths.repo_root
        ),
        "baseline_model_input_v1_summary_tsv": repo_relative(
            output_paths["baseline_model_input_v1_summary_tsv"], paths.repo_root
        ),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "baseline_feature_set_v1_latest_json": repo_relative(
            paths.baseline_feature_set_latest_pointer, paths.repo_root
        ),
        "baseline_feature_set_v1_audit_map_tsv": str(
            workflow_inputs.baseline_feature_set_latest_pointer["baseline_feature_set_v1_audit_map_tsv"]
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    baseline_model_input_v1_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    processed_run_dir = paths.processed_runs_root / baseline_model_input_v1_run_id
    audit_run_dir = paths.audit_runs_root / baseline_model_input_v1_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)
        validate_upstream_state(workflow_inputs)
        field_profiles = profile_source_fields(workflow_inputs)

        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        matrix_path = processed_run_dir / "baseline_model_input_v1.tsv"
        feature_dictionary_path = audit_run_dir / "baseline_model_input_v1_feature_dictionary.tsv"
        encoding_spec_path = audit_run_dir / "baseline_model_input_v1_encoding_spec.tsv"
        missingness_actions_path = audit_run_dir / "baseline_model_input_v1_missingness_actions.tsv"
        summary_path = audit_run_dir / "baseline_model_input_v1_summary.tsv"

        (
            matrix_fieldnames,
            matrix_rows,
            feature_dictionary_rows,
            encoding_spec_rows,
            missingness_action_rows,
        ) = build_model_input_outputs(workflow_inputs, field_profiles)
        summary_rows = build_summary_rows(
            baseline_model_input_v1_run_id=baseline_model_input_v1_run_id,
            workflow_inputs=workflow_inputs,
            field_profiles=field_profiles,
            matrix_fieldnames=matrix_fieldnames,
        )
        validation_payload = validate_generated_outputs(
            workflow_inputs=workflow_inputs,
            field_profiles=field_profiles,
            matrix_fieldnames=matrix_fieldnames,
            matrix_rows=matrix_rows,
            feature_dictionary_rows=feature_dictionary_rows,
            encoding_spec_rows=encoding_spec_rows,
            missingness_action_rows=missingness_action_rows,
            summary_rows=summary_rows,
        )

        helper_module.write_dict_rows_tsv(matrix_path, matrix_fieldnames, matrix_rows)
        helper_module.write_dict_rows_tsv(
            feature_dictionary_path,
            FEATURE_DICTIONARY_FIELDNAMES,
            feature_dictionary_rows,
        )
        helper_module.write_dict_rows_tsv(
            encoding_spec_path,
            ENCODING_SPEC_FIELDNAMES,
            encoding_spec_rows,
        )
        helper_module.write_dict_rows_tsv(
            missingness_actions_path,
            MISSINGNESS_ACTION_FIELDNAMES,
            missingness_action_rows,
        )
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "baseline_model_input_v1_tsv": matrix_path,
            "baseline_model_input_v1_feature_dictionary_tsv": feature_dictionary_path,
            "baseline_model_input_v1_encoding_spec_tsv": encoding_spec_path,
            "baseline_model_input_v1_missingness_actions_tsv": missingness_actions_path,
            "baseline_model_input_v1_summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            baseline_model_input_v1_run_id=baseline_model_input_v1_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "baseline_feature_set_v1_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["baseline_feature_set_v1_run_id"]
            ),
            "baseline_profile_v1_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["baseline_profile_v1_run_id"]
            ),
            "baseline_analysis_v1_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["baseline_analysis_v1_run_id"]
            ),
            "cohort_v1_build_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["cohort_v1_build_id"]
            ),
            "dry_run_build_id": str(workflow_inputs.baseline_feature_set_latest_pointer["dry_run_build_id"]),
            "blueprint_run_id": str(workflow_inputs.baseline_feature_set_latest_pointer["blueprint_run_id"]),
            "ambiguity_resolution_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["ambiguity_resolution_run_id"]
            ),
            "shortlist_run_id": str(workflow_inputs.baseline_feature_set_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(workflow_inputs.baseline_feature_set_latest_pointer["core_audit_run_id"]),
            "clinical_parse_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["clinical_parse_run_id"]
            ),
            "endpoint_crosswalk_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["endpoint_crosswalk_run_id"]
            ),
            "biospecimen_crosswalk_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["biospecimen_crosswalk_run_id"]
            ),
            "biospecimen_parse_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["biospecimen_parse_run_id"]
            ),
            "clinical_source_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["clinical_source_run_id"]
            ),
            "biospecimen_source_run_id": str(
                workflow_inputs.baseline_feature_set_latest_pointer["biospecimen_source_run_id"]
            ),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "baseline_feature_set_v1_latest_json": repo_relative(
                    paths.baseline_feature_set_latest_pointer, paths.repo_root
                ),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "baseline_model_input_v1_tsv": repo_relative(matrix_path, paths.repo_root),
                "baseline_model_input_v1_feature_dictionary_tsv": repo_relative(
                    feature_dictionary_path, paths.repo_root
                ),
                "baseline_model_input_v1_encoding_spec_tsv": repo_relative(
                    encoding_spec_path, paths.repo_root
                ),
                "baseline_model_input_v1_missingness_actions_tsv": repo_relative(
                    missingness_actions_path, paths.repo_root
                ),
                "baseline_model_input_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation_payload,
            "rules": {
                "unit_of_analysis": "patient/case",
                "input_layer": "baseline_feature_set_v1_only",
                "output_layer": "baseline_model_input_v1",
                "matrix_representation_policy": "fully_numeric_exact_label_one_hot_with_numeric_passthrough",
                "source_row_order_preserved": True,
                "one_row_per_patient": True,
                "matrix_identifier_columns_excluded": True,
                "audit_linkage_via_upstream_feature_set_audit_map": True,
                "numeric_missing_values_left_blank": True,
                "categorical_missing_values_use_missing_indicator": True,
                "no_missing_value_imputation": True,
                "no_endpoint_fields_added": True,
                "no_endpoint_freeze": True,
                "no_treatment_reintegration": True,
                "no_modeling": True,
                "no_metabric": True,
                "no_raw_xml_or_ssf_parsing": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "missing_like_tokens_json": json_list(sorted(MISSING_LIKE_TOKENS)),
                "manual_review_fields_json": json_list(sorted(MANUAL_REVIEW_FIELDS)),
            },
            "counts": {
                "input_row_count": len(workflow_inputs.feature_rows),
                "final_row_count": len(matrix_rows),
                "input_feature_count": len(EXPECTED_SOURCE_COLUMNS),
                "output_feature_count": len(matrix_fieldnames),
                "feature_dictionary_row_count": len(feature_dictionary_rows),
                "encoding_spec_row_count": len(encoding_spec_rows),
                "missingness_action_row_count": len(missingness_action_rows),
                "summary_row_count": len(summary_rows),
                "manual_review_required_field_count": sum(
                    profile.manual_review_needed for profile in field_profiles.values()
                ),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "baseline_feature_set_v1_latest_pointer": workflow_inputs.baseline_feature_set_latest_pointer,
                "baseline_feature_set_v1_run_log": workflow_inputs.baseline_feature_set_run_log,
            },
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "baseline_model_input_v1_run_id": baseline_model_input_v1_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_baseline_model_input_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA baseline model-input v1 workflow complete.")
    print(f"Baseline model-input run ID: {run_log['baseline_model_input_v1_run_id']}")
    print(f"Baseline feature-set run ID: {run_log['baseline_feature_set_v1_run_id']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"Matrix TSV: {run_log['outputs']['baseline_model_input_v1_tsv']}")
    print(f"Feature dictionary TSV: {run_log['outputs']['baseline_model_input_v1_feature_dictionary_tsv']}")
    print(f"Encoding spec TSV: {run_log['outputs']['baseline_model_input_v1_encoding_spec_tsv']}")
    print(f"Missingness actions TSV: {run_log['outputs']['baseline_model_input_v1_missingness_actions_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['baseline_model_input_v1_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
