#!/usr/bin/env python
"""Prepare a provisional TCGA-BRCA baseline-analysis-prep v1 table from minimal cohort v1."""

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
    "retained_in_baseline_analysis_v1",
    "rationale",
    "ambiguity_carried_forward",
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
SUMMARY_FIELDNAMES = [
    "baseline_analysis_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
SPEC_NOTES_PLACEHOLDER = "[fill in during baseline analysis v1 review]"
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
FIELD_CATEGORY_MAP = {
    "included_audit_field": "audit_id_field",
    "included_baseline_field": "baseline_clinical_field",
    "included_followup_join_evidence_field": "followup_join_evidence_field",
    "included_endpoint_candidate_field": "endpoint_candidate_field",
    "included_biospecimen_sample_anchor_field": "biospecimen_sample_anchor_evidence_field",
    "excluded_treatment_field": "excluded_treatment_field",
    "excluded_sparse_timing_field": "excluded_sparse_timing_field",
    "excluded_biospecimen_sample_helper_field": "excluded_biospecimen_sample_helper_field",
    "excluded_child_or_side_biospecimen_field": "excluded_child_biospecimen_field",
}
DERIVED_PREP_FIELD_SPECS = [
    (
        "baseline_analysis_v1_row_id",
        "Assign a sequential row id after sorting by provisional_patient_row_id for review-only referencing.",
    ),
    (
        "row_has_followup_match",
        "Derive a reversible yes/no flag from followup_match_row_count > 0 without changing source values.",
    ),
    (
        "row_has_biospecimen_sample_match",
        "Derive a reversible yes/no flag from biospecimen_sample_match_row_count > 0 without changing source values.",
    ),
    (
        "row_has_both_patient_ids",
        "Derive a reversible yes/no flag when both patient identifiers are present in the retained row.",
    ),
    (
        "row_missing_one_patient_id",
        "Derive a reversible yes/no flag when exactly one of the retained patient identifiers is missing-like.",
    ),
]


class BaselineAnalysisV1Error(RuntimeError):
    """Raised when the baseline-analysis-prep workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the baseline-analysis-prep workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    minimal_cohort_latest_pointer: Path
    upstream_paths: Any


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved baseline-analysis-prep inputs loaded from saved audit layers."""

    upstream_inputs: Any
    minimal_cohort_latest_pointer: dict[str, Any]
    minimal_cohort_run_log: dict[str, Any]
    minimal_cohort_rows: list[dict[str, str]]
    minimal_cohort_spec_rows: list[dict[str, str]]
    minimal_cohort_summary_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise BaselineAnalysisV1Error(f"Required helper script not found: {script_path}")

    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise BaselineAnalysisV1Error(f"Unable to create an import spec for: {script_path}")

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
        raise BaselineAnalysisV1Error(f"Expected integer-like value for {label}: {value!r}") from exc


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
            raise BaselineAnalysisV1Error(f"Duplicate summary_metric detected: {metric}")
        lookup[metric] = row
    return lookup


def json_list(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=True)


def normalize_yes_no(value: bool) -> str:
    return "yes" if value else "no"


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


def missingness_note(field_name: str) -> str:
    if field_name.endswith("_json"):
        return "json_array_cell"
    if field_name in {
        "row_has_followup_match",
        "row_has_biospecimen_sample_match",
        "row_has_both_patient_ids",
        "row_missing_one_patient_id",
    }:
        return "derived_yes_no_flag"
    return "scalar"


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    upstream_paths = helper_module.build_workflow_paths()
    trial_config_data = helper_module.load_yaml(upstream_paths.trial_config)

    processed_root = upstream_paths.repo_root / str(
        trial_config_data.get("processed_data_root", "01-data/processed")
    )
    audit_root = upstream_paths.repo_root / str(
        trial_config_data.get("audit_root", "01-data/audit")
    )
    results_root = upstream_paths.repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    analysis_prep_root = audit_root / "tcga-brca" / "analysis-prep"

    return WorkflowPaths(
        repo_root=upstream_paths.repo_root,
        trial_config=upstream_paths.trial_config,
        results_root=results_root,
        processed_runs_root=processed_root / "tcga-brca" / "analysis-prep" / "baseline_v1_runs",
        audit_runs_root=analysis_prep_root / "baseline_v1_runs",
        latest_pointer=analysis_prep_root / "tcga_brca_baseline_analysis_v1_latest.json",
        minimal_cohort_latest_pointer=(
            audit_root / "tcga-brca" / "cohort" / "tcga_brca_minimal_cohort_v1_latest.json"
        ),
        upstream_paths=upstream_paths,
    )


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    upstream_inputs = helper_module.load_workflow_inputs(paths.upstream_paths)

    if not paths.minimal_cohort_latest_pointer.exists():
        raise BaselineAnalysisV1Error(
            f"Required minimal cohort v1 latest pointer not found: {paths.minimal_cohort_latest_pointer}"
        )

    minimal_cohort_latest_pointer = helper_module.load_json(paths.minimal_cohort_latest_pointer)
    helper_module.require_keys(
        minimal_cohort_latest_pointer,
        {
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
            "minimal_cohort_v1_tsv",
            "minimal_cohort_v1_spec_tsv",
            "minimal_cohort_v1_summary_tsv",
            "run_log_json",
        },
        "Minimal cohort v1 latest pointer",
        paths.minimal_cohort_latest_pointer,
    )

    expected_pairings = {
        "blueprint_run_id": str(upstream_inputs.blueprint_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            upstream_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(upstream_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(upstream_inputs.clinical_shortlist_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(upstream_inputs.clinical_biotab_latest_pointer["parse_run_id"]),
        "endpoint_crosswalk_run_id": str(
            upstream_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_crosswalk_run_id": str(
            upstream_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(upstream_inputs.biospecimen_biotab_latest_pointer["parse_run_id"]),
        "clinical_source_run_id": str(upstream_inputs.clinical_biotab_latest_pointer["source_run_id"]),
        "biospecimen_source_run_id": str(upstream_inputs.biospecimen_biotab_latest_pointer["source_run_id"]),
    }
    for key, expected_value in expected_pairings.items():
        observed_value = str(minimal_cohort_latest_pointer[key])
        if observed_value != expected_value:
            raise BaselineAnalysisV1Error(
                "Minimal cohort v1 latest pointer is out of sync with the current upstream latest pointers: "
                f"{key}={observed_value!r} vs expected {expected_value!r}."
            )

    input_paths = {
        "minimal_cohort_v1_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["minimal_cohort_v1_tsv"]),
            "minimal cohort v1 TSV",
        ),
        "minimal_cohort_v1_spec_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["minimal_cohort_v1_spec_tsv"]),
            "minimal cohort v1 spec TSV",
        ),
        "minimal_cohort_v1_summary_tsv": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["minimal_cohort_v1_summary_tsv"]),
            "minimal cohort v1 summary TSV",
        ),
        "minimal_cohort_v1_run_log_json": helper_module.resolve_existing_path(
            paths.repo_root,
            str(minimal_cohort_latest_pointer["run_log_json"]),
            "minimal cohort v1 run log",
        ),
    }

    minimal_cohort_run_log = helper_module.load_json(input_paths["minimal_cohort_v1_run_log_json"])
    if minimal_cohort_run_log.get("status") != "completed":
        raise BaselineAnalysisV1Error("The referenced minimal cohort v1 run log is not completed.")
    if not bool(minimal_cohort_run_log.get("validation", {}).get("passed", False)):
        raise BaselineAnalysisV1Error(
            "The referenced minimal cohort v1 run log does not report validation.passed == true."
        )

    minimal_cohort_rows = helper_module.read_tsv_dict_rows(input_paths["minimal_cohort_v1_tsv"])
    minimal_cohort_spec_rows = helper_module.read_tsv_dict_rows(input_paths["minimal_cohort_v1_spec_tsv"])
    minimal_cohort_summary_rows = helper_module.read_tsv_dict_rows(input_paths["minimal_cohort_v1_summary_tsv"])

    if not minimal_cohort_rows:
        raise BaselineAnalysisV1Error("The referenced minimal cohort v1 TSV contains no rows.")

    summary_lookup = build_summary_lookup(minimal_cohort_summary_rows)
    if "final_v1_row_count" not in summary_lookup:
        raise BaselineAnalysisV1Error(
            "Minimal cohort v1 summary TSV is missing summary_metric=final_v1_row_count."
        )
    summary_row_count = parse_int(
        summary_lookup["final_v1_row_count"]["summary_value"],
        "minimal cohort v1 final_v1_row_count",
    )
    if summary_row_count != len(minimal_cohort_rows):
        raise BaselineAnalysisV1Error(
            "Minimal cohort v1 summary row count does not match the cohort TSV row count: "
            f"{summary_row_count} vs {len(minimal_cohort_rows)}."
        )

    return WorkflowInputs(
        upstream_inputs=upstream_inputs,
        minimal_cohort_latest_pointer=minimal_cohort_latest_pointer,
        minimal_cohort_run_log=minimal_cohort_run_log,
        minimal_cohort_rows=minimal_cohort_rows,
        minimal_cohort_spec_rows=minimal_cohort_spec_rows,
        minimal_cohort_summary_rows=minimal_cohort_summary_rows,
        input_paths=input_paths,
    )


def require_mapped_field_category(source_category: str) -> str:
    mapped = FIELD_CATEGORY_MAP.get(source_category)
    if mapped is None:
        raise BaselineAnalysisV1Error(f"Unexpected minimal cohort v1 field_category: {source_category}")
    return mapped


def build_field_groups(
    minimal_cohort_spec_rows: list[dict[str, str]],
) -> dict[str, list[str]]:
    groups = {
        "audit_fields": [],
        "baseline_fields": [],
        "patient_endpoint_fields": [],
        "followup_join_evidence_fields": [],
        "followup_endpoint_fields": [],
        "biospecimen_anchor_fields": [],
    }

    for row in minimal_cohort_spec_rows:
        retained = str(row.get("include_in_v1") or "") == "yes"
        if not retained:
            continue

        source_table = str(row.get("source_table") or "")
        field_name = str(row.get("field_name") or "")
        field_category = require_mapped_field_category(str(row.get("field_category") or ""))

        if field_category == "audit_id_field":
            groups["audit_fields"].append(field_name)
        elif field_category == "baseline_clinical_field":
            groups["baseline_fields"].append(field_name)
        elif field_category == "followup_join_evidence_field":
            groups["followup_join_evidence_fields"].append(field_name)
        elif field_category == "endpoint_candidate_field" and source_table == "clinical_patient":
            groups["patient_endpoint_fields"].append(field_name)
        elif field_category == "endpoint_candidate_field" and source_table == "clinical_follow_up_v4_0":
            groups["followup_endpoint_fields"].append(field_name)
        elif field_category == "biospecimen_sample_anchor_evidence_field":
            groups["biospecimen_anchor_fields"].append(field_name)
        else:
            raise BaselineAnalysisV1Error(
                "Unable to place retained minimal cohort v1 field into a baseline-analysis-prep group: "
                f"field_name={field_name!r}, source_table={source_table!r}, field_category={field_category!r}."
            )

    for group_name, values in groups.items():
        groups[group_name] = ordered_unique(values)

    return groups


def build_output_column_order(field_groups: dict[str, list[str]]) -> list[str]:
    return [
        "cohort_v1_build_id",
        "baseline_analysis_v1_row_id",
        "provisional_patient_row_id",
        "bcr_patient_barcode",
        "bcr_patient_uuid",
        "ambiguity_note_flags_json",
        "join_status_flags_json",
        "followup_match_row_count",
        "biospecimen_sample_match_row_count",
        "row_has_followup_match",
        "row_has_biospecimen_sample_match",
        "row_has_both_patient_ids",
        "row_missing_one_patient_id",
        *field_groups["baseline_fields"],
        *field_groups["patient_endpoint_fields"],
        *field_groups["followup_join_evidence_fields"],
        *field_groups["followup_endpoint_fields"],
        *field_groups["biospecimen_anchor_fields"],
    ]


def sort_minimal_cohort_rows(
    minimal_cohort_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    sorted_rows = sorted(
        minimal_cohort_rows,
        key=lambda row: parse_int(row["provisional_patient_row_id"], "provisional_patient_row_id"),
    )

    observed_ids = [parse_int(row["provisional_patient_row_id"], "provisional_patient_row_id") for row in sorted_rows]
    if len(set(observed_ids)) != len(observed_ids):
        raise BaselineAnalysisV1Error("minimal_cohort_v1.tsv does not preserve one row per provisional_patient_row_id.")

    return sorted_rows


def build_prepared_rows(
    minimal_cohort_rows: list[dict[str, str]],
    output_column_order: list[str],
) -> list[dict[str, str]]:
    sorted_rows = sort_minimal_cohort_rows(minimal_cohort_rows)
    prepared_rows: list[dict[str, str]] = []

    for row_number, source_row in enumerate(sorted_rows, start=1):
        has_followup_match = parse_int(
            source_row["followup_match_row_count"], "followup_match_row_count"
        ) > 0
        has_biospecimen_sample_match = parse_int(
            source_row["biospecimen_sample_match_row_count"],
            "biospecimen_sample_match_row_count",
        ) > 0
        has_patient_barcode = not is_scalar_missing_like(source_row["bcr_patient_barcode"])
        has_patient_uuid = not is_scalar_missing_like(source_row["bcr_patient_uuid"])

        derived_values = {
            "baseline_analysis_v1_row_id": str(row_number),
            "row_has_followup_match": normalize_yes_no(has_followup_match),
            "row_has_biospecimen_sample_match": normalize_yes_no(has_biospecimen_sample_match),
            "row_has_both_patient_ids": normalize_yes_no(has_patient_barcode and has_patient_uuid),
            "row_missing_one_patient_id": normalize_yes_no(has_patient_barcode ^ has_patient_uuid),
        }

        prepared_row: dict[str, str] = {}
        for column_name in output_column_order:
            if column_name in derived_values:
                prepared_row[column_name] = derived_values[column_name]
            else:
                prepared_row[column_name] = str(source_row.get(column_name, ""))
        prepared_rows.append(prepared_row)

    return prepared_rows


def build_spec_rows(
    minimal_cohort_spec_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    spec_rows: list[dict[str, str]] = []
    for row in minimal_cohort_spec_rows:
        source_category = str(row.get("field_category") or "")
        mapped_category = require_mapped_field_category(source_category)
        spec_rows.append(
            {
                "field_name": str(row.get("field_name") or ""),
                "source_origin": str(row.get("source_table") or ""),
                "field_category": mapped_category,
                "retained_in_baseline_analysis_v1": (
                    "yes" if str(row.get("include_in_v1") or "") == "yes" else "no"
                ),
                "rationale": str(row.get("inclusion_rule") or ""),
                "ambiguity_carried_forward": str(row.get("ambiguity_carried_forward") or ""),
                "notes_placeholder": SPEC_NOTES_PLACEHOLDER,
            }
        )

    for field_name, rationale in DERIVED_PREP_FIELD_SPECS:
        spec_rows.append(
            {
                "field_name": field_name,
                "source_origin": "workflow_derived",
                "field_category": "derived_prep_flag",
                "retained_in_baseline_analysis_v1": "yes",
                "rationale": rationale,
                "ambiguity_carried_forward": "no",
                "notes_placeholder": SPEC_NOTES_PLACEHOLDER,
            }
        )

    if len(spec_rows) != len(minimal_cohort_spec_rows) + len(DERIVED_PREP_FIELD_SPECS):
        raise BaselineAnalysisV1Error("baseline_analysis_v1_spec.tsv row count did not reconcile.")

    return spec_rows


def build_missingness_rows(
    prepared_rows: list[dict[str, str]],
    output_column_order: list[str],
) -> list[dict[str, str]]:
    row_count = len(prepared_rows)
    missingness_rows: list[dict[str, str]] = []

    for field_name in output_column_order:
        non_missing_values: list[str] = []
        missing_like_count = 0
        for row in prepared_rows:
            raw_value = str(row.get(field_name, ""))
            if cell_is_missing_like(field_name, raw_value):
                missing_like_count += 1
                continue
            non_missing_values.append(raw_value)

        non_missing_count = len(non_missing_values)
        missingness_rows.append(
            {
                "field_name": field_name,
                "row_count": str(row_count),
                "non_missing_count": str(non_missing_count),
                "missing_like_count": str(missing_like_count),
                "missing_like_fraction": f"{(missing_like_count / row_count) if row_count else 0.0:.6f}",
                "distinct_non_missing_count": str(len(set(non_missing_values))),
                "notes": missingness_note(field_name),
            }
        )

    return missingness_rows


def distinct_json_flags(prepared_rows: list[dict[str, str]], field_name: str) -> list[str]:
    observed: list[str] = []
    seen: set[str] = set()
    for row in prepared_rows:
        raw_value = str(row.get(field_name, ""))
        if raw_value == "":
            continue
        parsed = json.loads(raw_value)
        if not isinstance(parsed, list):
            raise BaselineAnalysisV1Error(f"Expected {field_name} to contain JSON lists.")
        for item in parsed:
            item_str = str(item)
            if item_str in seen:
                continue
            seen.add(item_str)
            observed.append(item_str)
    return observed


def count_yes(prepared_rows: list[dict[str, str]], field_name: str) -> int:
    return sum(str(row.get(field_name) or "") == "yes" for row in prepared_rows)


def build_summary_rows(
    baseline_analysis_v1_run_id: str,
    workflow_inputs: WorkflowInputs,
    prepared_rows: list[dict[str, str]],
    spec_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    retained_rows = [row for row in spec_rows if row["retained_in_baseline_analysis_v1"] == "yes"]
    retained_category_counts: dict[str, int] = {}
    for row in retained_rows:
        retained_category_counts[row["field_category"]] = retained_category_counts.get(row["field_category"], 0) + 1

    carried_ambiguity_flags = distinct_json_flags(prepared_rows, "ambiguity_note_flags_json")
    final_row_count = len(prepared_rows)
    rows_with_followup_evidence = count_yes(prepared_rows, "row_has_followup_match")
    rows_with_biospecimen_sample_evidence = count_yes(
        prepared_rows, "row_has_biospecimen_sample_match"
    )
    rows_with_both_patient_ids = count_yes(prepared_rows, "row_has_both_patient_ids")
    rows_with_only_one_patient_id = count_yes(prepared_rows, "row_missing_one_patient_id")
    rows_with_no_patient_ids = sum(
        str(row["row_has_both_patient_ids"]) == "no" and str(row["row_missing_one_patient_id"]) == "no"
        for row in prepared_rows
    )

    return [
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "design",
            "summary_metric": "baseline_analysis_v1_status",
            "summary_value": "provisional_baseline_analysis_prep_v1",
            "notes": "This table is a baseline-analysis-prep layer only and is not a final modeling dataset.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "design",
            "summary_metric": "proposed_unit_of_analysis",
            "summary_value": "patient/case",
            "notes": "The output preserves one row per patient/case from minimal cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "design",
            "summary_metric": "baseline_analysis_source_table",
            "summary_value": "minimal_cohort_v1",
            "notes": "baseline_analysis_v1 starts only from the saved minimal cohort v1 output on disk.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "design",
            "summary_metric": "endpoint_policy",
            "summary_value": "candidate_columns_side_by_side_only",
            "notes": "Endpoint candidate columns remain separate, provenance-bearing, and provisional.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "design",
            "summary_metric": "treatment_policy",
            "summary_value": "excluded",
            "notes": "Treatment detail rows and treatment proxy fields remain excluded from baseline analysis prep v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "design",
            "summary_metric": "biospecimen_policy",
            "summary_value": "sample_anchor_only",
            "notes": "biospecimen evidence remains limited to grouped sample-anchor columns already present in cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "cohort_v1_build_id",
            "summary_value": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
            "notes": "Minimal cohort v1 build carried forward into this baseline-analysis-prep run.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "dry_run_build_id",
            "summary_value": str(workflow_inputs.minimal_cohort_latest_pointer["dry_run_build_id"]),
            "notes": "Original minimal dry-run build referenced by the current minimal cohort v1 pointer.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "blueprint_run_id",
            "summary_value": str(workflow_inputs.minimal_cohort_latest_pointer["blueprint_run_id"]),
            "notes": "Current cohort blueprint run aligned to the retained minimal cohort v1 pointer.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "ambiguity_resolution_run_id",
            "summary_value": str(workflow_inputs.minimal_cohort_latest_pointer["ambiguity_resolution_run_id"]),
            "notes": "Current ambiguity-resolution run aligned to the retained minimal cohort v1 pointer.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "shortlist_run_id",
            "summary_value": str(workflow_inputs.minimal_cohort_latest_pointer["shortlist_run_id"]),
            "notes": "Current clinical shortlist run aligned to the retained minimal cohort v1 pointer.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "inputs",
            "summary_metric": "endpoint_crosswalk_run_id",
            "summary_value": str(workflow_inputs.minimal_cohort_latest_pointer["endpoint_crosswalk_run_id"]),
            "notes": "Current endpoint crosswalk run aligned to the retained minimal cohort v1 pointer.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "input_v1_row_count",
            "summary_value": str(len(workflow_inputs.minimal_cohort_rows)),
            "notes": "Row count read directly from minimal_cohort_v1.tsv.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "row_counts",
            "summary_metric": "final_row_count",
            "summary_value": str(final_row_count),
            "notes": "Row count of baseline_analysis_v1.tsv after sorting and appending reversible prep columns.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_field_count_total",
            "summary_value": str(len(retained_rows)),
            "notes": "Total retained field count including carried cohort v1 fields plus derived prep flags.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_audit_id_field_count",
            "summary_value": str(retained_category_counts.get("audit_id_field", 0)),
            "notes": "Retained audit/id fields carried from minimal cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_baseline_clinical_field_count",
            "summary_value": str(retained_category_counts.get("baseline_clinical_field", 0)),
            "notes": "Retained baseline clinical fields carried from minimal cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_followup_join_evidence_field_count",
            "summary_value": str(retained_category_counts.get("followup_join_evidence_field", 0)),
            "notes": "Retained follow-up join-evidence fields carried from minimal cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_endpoint_candidate_field_count",
            "summary_value": str(retained_category_counts.get("endpoint_candidate_field", 0)),
            "notes": "Retained endpoint candidate fields carried side by side from minimal cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_biospecimen_sample_anchor_evidence_field_count",
            "summary_value": str(
                retained_category_counts.get("biospecimen_sample_anchor_evidence_field", 0)
            ),
            "notes": "Retained biospecimen sample-anchor evidence fields carried from minimal cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "field_selection",
            "summary_metric": "retained_derived_prep_flag_count",
            "summary_value": str(retained_category_counts.get("derived_prep_flag", 0)),
            "notes": "New reversible prep columns added by this workflow.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "coverage",
            "summary_metric": "rows_with_followup_evidence",
            "summary_value": str(rows_with_followup_evidence),
            "notes": "Rows where followup_match_row_count > 0.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "coverage",
            "summary_metric": "rows_with_biospecimen_sample_evidence",
            "summary_value": str(rows_with_biospecimen_sample_evidence),
            "notes": "Rows where biospecimen_sample_match_row_count > 0.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "coverage",
            "summary_metric": "rows_with_both_patient_ids",
            "summary_value": str(rows_with_both_patient_ids),
            "notes": "Rows where both retained patient identifiers are present.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "coverage",
            "summary_metric": "rows_with_only_one_patient_id",
            "summary_value": str(rows_with_only_one_patient_id),
            "notes": "Rows where exactly one retained patient identifier is missing-like.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "coverage",
            "summary_metric": "rows_with_no_patient_ids",
            "summary_value": str(rows_with_no_patient_ids),
            "notes": "Rows where both retained patient identifiers are missing-like.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "ambiguity",
            "summary_metric": "distinct_carried_ambiguity_flag_count",
            "summary_value": str(len(carried_ambiguity_flags)),
            "notes": "Distinct ambiguity flags still carried forward at the row level.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "ambiguity",
            "summary_metric": "carried_ambiguity_flags_json",
            "summary_value": json_list(carried_ambiguity_flags),
            "notes": "Distinct ambiguity flags retained from minimal cohort v1.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "readiness_for_descriptive_stats",
            "summary_value": "ready_provisional_baseline_analysis_v1",
            "notes": "This table is suitable for descriptive statistics review while remaining provisional.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "readiness_for_missingness_review",
            "summary_value": "ready",
            "notes": "This table is suitable for direct missingness review because values remain source-faithful.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "readiness_for_baseline_modeling_prep",
            "summary_value": "ready_provisional_only",
            "notes": "This table can support baseline modeling preparation only while endpoint and treatment remain unresolved.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "endpoint_freeze_status",
            "summary_value": "blocked",
            "notes": "Final endpoint freeze remains blocked and endpoint candidate columns remain side by side only.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "treatment_inclusion_status",
            "summary_value": "excluded",
            "notes": "Treatment detail remains excluded from this baseline-analysis-prep layer.",
        },
        {
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "summary_section": "readiness",
            "summary_metric": "overall_readiness_interpretation",
            "summary_value": (
                "suitable_for_descriptive_statistics_and_baseline_prep_only_not_for_endpoint_freeze_"
                "or_treatment_modeling"
            ),
            "notes": (
                "Use this table for descriptive statistics, missingness review, and baseline modeling preparation only; "
                "do not use it to claim a frozen endpoint, treatment modeling readiness, or stronger clinical-effect interpretation."
            ),
        },
    ]


def build_latest_pointer_payload(
    baseline_analysis_v1_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
        "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
        "dry_run_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["dry_run_build_id"]),
        "blueprint_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            workflow_inputs.minimal_cohort_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["clinical_parse_run_id"]),
        "endpoint_crosswalk_run_id": str(
            workflow_inputs.minimal_cohort_latest_pointer["endpoint_crosswalk_run_id"]
        ),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.minimal_cohort_latest_pointer["biospecimen_crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(
            workflow_inputs.minimal_cohort_latest_pointer["biospecimen_parse_run_id"]
        ),
        "clinical_source_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["clinical_source_run_id"]),
        "biospecimen_source_run_id": str(
            workflow_inputs.minimal_cohort_latest_pointer["biospecimen_source_run_id"]
        ),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "baseline_analysis_v1_tsv": repo_relative(output_paths["baseline_analysis_v1_tsv"], paths.repo_root),
        "baseline_analysis_v1_spec_tsv": repo_relative(
            output_paths["baseline_analysis_v1_spec_tsv"], paths.repo_root
        ),
        "baseline_analysis_v1_missingness_tsv": repo_relative(
            output_paths["baseline_analysis_v1_missingness_tsv"], paths.repo_root
        ),
        "baseline_analysis_v1_summary_tsv": repo_relative(
            output_paths["baseline_analysis_v1_summary_tsv"], paths.repo_root
        ),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "minimal_cohort_v1_latest_json": repo_relative(
            paths.minimal_cohort_latest_pointer, paths.repo_root
        ),
        "cohort_blueprint_latest_json": repo_relative(
            paths.upstream_paths.blueprint_latest_pointer,
            paths.repo_root,
        ),
        "blueprint_ambiguity_resolution_latest_json": repo_relative(
            paths.upstream_paths.ambiguity_resolution_latest_pointer,
            paths.repo_root,
        ),
        "clinical_shortlist_latest_json": repo_relative(
            paths.upstream_paths.clinical_shortlist_latest_pointer,
            paths.repo_root,
        ),
        "endpoint_crosswalk_latest_json": repo_relative(
            paths.upstream_paths.endpoint_crosswalk_latest_pointer,
            paths.repo_root,
        ),
        "biospecimen_crosswalk_latest_json": repo_relative(
            paths.upstream_paths.biospecimen_crosswalk_latest_pointer,
            paths.repo_root,
        ),
        "clinical_biotab_latest_json": repo_relative(
            paths.upstream_paths.clinical_biotab_latest_pointer,
            paths.repo_root,
        ),
        "biospecimen_biotab_latest_json": repo_relative(
            paths.upstream_paths.biospecimen_biotab_latest_pointer,
            paths.repo_root,
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    baseline_analysis_v1_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    processed_run_dir = paths.processed_runs_root / baseline_analysis_v1_run_id
    audit_run_dir = paths.audit_runs_root / baseline_analysis_v1_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)
        field_groups = build_field_groups(workflow_inputs.minimal_cohort_spec_rows)
        output_column_order = build_output_column_order(field_groups)

        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        baseline_analysis_tsv_path = processed_run_dir / "baseline_analysis_v1.tsv"
        spec_path = audit_run_dir / "baseline_analysis_v1_spec.tsv"
        missingness_path = audit_run_dir / "baseline_analysis_v1_missingness.tsv"
        summary_path = audit_run_dir / "baseline_analysis_v1_summary.tsv"

        prepared_rows = build_prepared_rows(workflow_inputs.minimal_cohort_rows, output_column_order)
        spec_rows = build_spec_rows(workflow_inputs.minimal_cohort_spec_rows)
        missingness_rows = build_missingness_rows(prepared_rows, output_column_order)
        summary_rows = build_summary_rows(
            baseline_analysis_v1_run_id=baseline_analysis_v1_run_id,
            workflow_inputs=workflow_inputs,
            prepared_rows=prepared_rows,
            spec_rows=spec_rows,
        )

        if not prepared_rows:
            raise BaselineAnalysisV1Error("baseline_analysis_v1.tsv rows were not generated.")
        if not spec_rows:
            raise BaselineAnalysisV1Error("baseline_analysis_v1_spec.tsv rows were not generated.")
        if not missingness_rows:
            raise BaselineAnalysisV1Error("baseline_analysis_v1_missingness.tsv rows were not generated.")
        if not summary_rows:
            raise BaselineAnalysisV1Error("baseline_analysis_v1_summary.tsv rows were not generated.")

        helper_module.write_dict_rows_tsv(baseline_analysis_tsv_path, output_column_order, prepared_rows)
        helper_module.write_dict_rows_tsv(spec_path, SPEC_FIELDNAMES, spec_rows)
        helper_module.write_dict_rows_tsv(missingness_path, MISSINGNESS_FIELDNAMES, missingness_rows)
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        retained_spec_field_names = ordered_unique(
            [
                row["field_name"]
                for row in spec_rows
                if row["retained_in_baseline_analysis_v1"] == "yes"
            ]
        )
        retained_spec_rows_match_output_columns = (
            len(retained_spec_field_names) == len(output_column_order)
            and set(retained_spec_field_names) == set(output_column_order)
        )
        output_row_count_matches_v1 = len(prepared_rows) == len(workflow_inputs.minimal_cohort_rows)
        provisional_patient_ids = [row["provisional_patient_row_id"] for row in prepared_rows]
        one_row_per_provisional_patient_row = len(set(provisional_patient_ids)) == len(provisional_patient_ids)
        baseline_row_ids = [
            parse_int(row["baseline_analysis_v1_row_id"], "baseline_analysis_v1_row_id")
            for row in prepared_rows
        ]
        sequential_row_ids = baseline_row_ids == list(range(1, len(prepared_rows) + 1))
        required_source_table_found = workflow_inputs.input_paths["minimal_cohort_v1_tsv"].exists()
        output_rows_positive = len(prepared_rows) > 0
        summary_lookup = build_summary_lookup(summary_rows)
        summary_endpoint_blocked = summary_lookup["endpoint_freeze_status"]["summary_value"] == "blocked"
        summary_treatment_excluded = (
            summary_lookup["treatment_inclusion_status"]["summary_value"] == "excluded"
        )

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "baseline_analysis_v1_tsv": baseline_analysis_tsv_path,
            "baseline_analysis_v1_spec_tsv": spec_path,
            "baseline_analysis_v1_missingness_tsv": missingness_path,
            "baseline_analysis_v1_summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            baseline_analysis_v1_run_id=baseline_analysis_v1_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "cohort_v1_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["cohort_v1_build_id"]),
            "dry_run_build_id": str(workflow_inputs.minimal_cohort_latest_pointer["dry_run_build_id"]),
            "blueprint_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["blueprint_run_id"]),
            "ambiguity_resolution_run_id": str(
                workflow_inputs.minimal_cohort_latest_pointer["ambiguity_resolution_run_id"]
            ),
            "shortlist_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["core_audit_run_id"]),
            "clinical_parse_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["clinical_parse_run_id"]),
            "endpoint_crosswalk_run_id": str(
                workflow_inputs.minimal_cohort_latest_pointer["endpoint_crosswalk_run_id"]
            ),
            "biospecimen_crosswalk_run_id": str(
                workflow_inputs.minimal_cohort_latest_pointer["biospecimen_crosswalk_run_id"]
            ),
            "biospecimen_parse_run_id": str(
                workflow_inputs.minimal_cohort_latest_pointer["biospecimen_parse_run_id"]
            ),
            "clinical_source_run_id": str(workflow_inputs.minimal_cohort_latest_pointer["clinical_source_run_id"]),
            "biospecimen_source_run_id": str(
                workflow_inputs.minimal_cohort_latest_pointer["biospecimen_source_run_id"]
            ),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "minimal_cohort_v1_latest_json": repo_relative(
                    paths.minimal_cohort_latest_pointer, paths.repo_root
                ),
                "cohort_blueprint_latest_json": repo_relative(
                    paths.upstream_paths.blueprint_latest_pointer, paths.repo_root
                ),
                "blueprint_ambiguity_resolution_latest_json": repo_relative(
                    paths.upstream_paths.ambiguity_resolution_latest_pointer,
                    paths.repo_root,
                ),
                "clinical_shortlist_latest_json": repo_relative(
                    paths.upstream_paths.clinical_shortlist_latest_pointer,
                    paths.repo_root,
                ),
                "endpoint_crosswalk_latest_json": repo_relative(
                    paths.upstream_paths.endpoint_crosswalk_latest_pointer,
                    paths.repo_root,
                ),
                "biospecimen_crosswalk_latest_json": repo_relative(
                    paths.upstream_paths.biospecimen_crosswalk_latest_pointer,
                    paths.repo_root,
                ),
                "clinical_biotab_latest_json": repo_relative(
                    paths.upstream_paths.clinical_biotab_latest_pointer,
                    paths.repo_root,
                ),
                "biospecimen_biotab_latest_json": repo_relative(
                    paths.upstream_paths.biospecimen_biotab_latest_pointer,
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
                "baseline_analysis_v1_tsv": repo_relative(
                    baseline_analysis_tsv_path, paths.repo_root
                ),
                "baseline_analysis_v1_spec_tsv": repo_relative(spec_path, paths.repo_root),
                "baseline_analysis_v1_missingness_tsv": repo_relative(
                    missingness_path, paths.repo_root
                ),
                "baseline_analysis_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    output_rows_positive
                    and output_row_count_matches_v1
                    and retained_spec_rows_match_output_columns
                    and sequential_row_ids
                    and one_row_per_provisional_patient_row
                    and required_source_table_found
                    and summary_endpoint_blocked
                    and summary_treatment_excluded
                ),
                "required_pointers_found": True,
                "minimal_cohort_v1_latest_pointer_found": True,
                "cohort_blueprint_latest_pointer_found": True,
                "ambiguity_resolution_latest_pointer_found": True,
                "clinical_shortlist_latest_pointer_found": True,
                "endpoint_crosswalk_latest_pointer_found": True,
                "biospecimen_crosswalk_latest_pointer_found": True,
                "minimal_cohort_run_log_completed": True,
                "minimal_cohort_validation_passed": True,
                "required_source_table_found": required_source_table_found,
                "output_rows_positive": output_rows_positive,
                "output_row_count_matches_v1": output_row_count_matches_v1,
                "retained_spec_rows_match_output_columns": retained_spec_rows_match_output_columns,
                "sequential_row_ids": sequential_row_ids,
                "one_row_per_provisional_patient_row": one_row_per_provisional_patient_row,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "unit_of_analysis": "patient/case",
                "input_layer": "minimal_cohort_v1_only",
                "output_layer": "baseline_analysis_prep_v1",
                "endpoint_policy": "candidate_columns_side_by_side_only",
                "treatment_policy": "excluded",
                "biospecimen_policy": "sample_anchor_only",
                "no_missing_value_imputation": True,
                "no_endpoint_freeze": True,
                "no_treatment_reintegration": True,
                "no_child_biospecimen_expansion": True,
                "no_modeling": True,
                "no_metabric": True,
                "no_raw_xml_or_ssf_parsing": True,
                "one_row_per_patient": True,
                "derived_columns_reversible_only": True,
                "existing_values_unchanged_for_carried_fields": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "missing_like_tokens_json": json_list(sorted(MISSING_LIKE_TOKENS)),
            },
            "counts": {
                "input_v1_row_count": len(workflow_inputs.minimal_cohort_rows),
                "final_row_count": len(prepared_rows),
                "retained_field_count_total": len(retained_spec_field_names),
                "retained_spec_row_count": len(retained_spec_field_names),
                "spec_row_count": len(spec_rows),
                "missingness_row_count": len(missingness_rows),
                "summary_row_count": len(summary_rows),
                "rows_with_followup_evidence": count_yes(prepared_rows, "row_has_followup_match"),
                "rows_with_biospecimen_sample_evidence": count_yes(
                    prepared_rows, "row_has_biospecimen_sample_match"
                ),
                "rows_with_both_patient_ids": count_yes(prepared_rows, "row_has_both_patient_ids"),
                "rows_with_only_one_patient_id": count_yes(
                    prepared_rows, "row_missing_one_patient_id"
                ),
                "rows_with_no_patient_ids": sum(
                    str(row["row_has_both_patient_ids"]) == "no"
                    and str(row["row_missing_one_patient_id"]) == "no"
                    for row in prepared_rows
                ),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "minimal_cohort_v1_latest_pointer": workflow_inputs.minimal_cohort_latest_pointer,
                "minimal_cohort_v1_run_log": workflow_inputs.minimal_cohort_run_log,
                "cohort_blueprint_latest_pointer": workflow_inputs.upstream_inputs.blueprint_latest_pointer,
                "cohort_blueprint_run_log": workflow_inputs.upstream_inputs.blueprint_run_log,
                "ambiguity_resolution_latest_pointer": workflow_inputs.upstream_inputs.ambiguity_resolution_latest_pointer,
                "ambiguity_resolution_run_log": workflow_inputs.upstream_inputs.ambiguity_resolution_run_log,
                "clinical_shortlist_latest_pointer": workflow_inputs.upstream_inputs.clinical_shortlist_latest_pointer,
                "clinical_shortlist_run_log": workflow_inputs.upstream_inputs.clinical_shortlist_run_log,
                "endpoint_crosswalk_latest_pointer": workflow_inputs.upstream_inputs.endpoint_crosswalk_latest_pointer,
                "endpoint_crosswalk_run_log": workflow_inputs.upstream_inputs.endpoint_crosswalk_run_log,
                "biospecimen_crosswalk_latest_pointer": workflow_inputs.upstream_inputs.biospecimen_crosswalk_latest_pointer,
                "biospecimen_crosswalk_run_log": workflow_inputs.upstream_inputs.biospecimen_crosswalk_run_log,
                "clinical_biotab_latest_pointer": workflow_inputs.upstream_inputs.clinical_biotab_latest_pointer,
                "clinical_biotab_run_log": workflow_inputs.upstream_inputs.clinical_biotab_run_log,
                "biospecimen_biotab_latest_pointer": workflow_inputs.upstream_inputs.biospecimen_biotab_latest_pointer,
                "biospecimen_biotab_run_log": workflow_inputs.upstream_inputs.biospecimen_biotab_run_log,
            },
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "baseline_analysis_v1_run_id": baseline_analysis_v1_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_baseline_analysis_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA baseline analysis prep v1 workflow complete.")
    print(f"Baseline analysis prep run ID: {run_log['baseline_analysis_v1_run_id']}")
    print(f"Cohort v1 build ID: {run_log['cohort_v1_build_id']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"Baseline analysis TSV: {run_log['outputs']['baseline_analysis_v1_tsv']}")
    print(f"Spec TSV: {run_log['outputs']['baseline_analysis_v1_spec_tsv']}")
    print(f"Missingness TSV: {run_log['outputs']['baseline_analysis_v1_missingness_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['baseline_analysis_v1_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
