#!/usr/bin/env python
"""Build a provisional minimal TCGA-BRCA cohort v1 from saved audit layers."""

from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

JOIN_AUDIT_FIELDNAMES = [
    "cohort_v1_build_id",
    "join_step_order",
    "join_step_name",
    "left_table_name",
    "right_table_name",
    "starting_left_row_count",
    "matched_left_row_count",
    "unmatched_left_row_count",
    "matched_right_row_count",
    "right_only_row_count",
    "left_rows_with_multiple_matches",
    "max_right_matches_per_left_row",
    "join_rule_used",
    "included_in_v1",
    "carried_forward_as_provisional",
    "excluded_from_v1",
    "ambiguity_carried_forward",
    "notes",
]
EXCLUSIONS_FIELDNAMES = [
    "cohort_v1_build_id",
    "exclusion_category",
    "scope_level",
    "source_table",
    "entity_id",
    "entity_value",
    "reason",
    "count_value",
    "details",
]
SUMMARY_FIELDNAMES = [
    "cohort_v1_build_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
SPEC_FIELDNAMES = [
    "field_name",
    "source_table",
    "field_category",
    "include_in_v1",
    "inclusion_rule",
    "ambiguity_carried_forward",
    "notes_placeholder",
]
NOTES_PLACEHOLDER = "[fill in during minimal cohort v1 review]"
EXPECTED_DRY_RUN_CONTRACT = {
    "final_row_count": 1097,
    "followup_matched_patient_count": 619,
    "followup_matched_row_count": 716,
    "biospecimen_matched_patient_count": 1097,
    "biospecimen_matched_row_count": 2293,
    "mismatch_signal": "1097 / 1098 / 1101",
}


class MinimalCohortV1Error(RuntimeError):
    """Raised when the minimal cohort v1 workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the minimal cohort v1 workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    dry_run_latest_pointer: Path
    upstream_paths: Any


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved minimal cohort v1 inputs loaded from saved audit layers."""

    upstream_inputs: Any
    dry_run_latest_pointer: dict[str, Any]
    dry_run_run_log: dict[str, Any]
    dry_run_join_audit_rows: list[dict[str, str]]
    dry_run_row_counts_rows: list[dict[str, str]]
    dry_run_exclusions_rows: list[dict[str, str]]
    dry_run_summary_rows: list[dict[str, str]]
    dry_run_contract: dict[str, Any]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_dry_run_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise MinimalCohortV1Error(f"Required dry-run script not found: {script_path}")

    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise MinimalCohortV1Error(f"Unable to create an import spec for: {script_path}")

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
        raise MinimalCohortV1Error(f"Expected integer-like value for {label}: {value!r}") from exc


def build_summary_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        metric = str(row.get("summary_metric") or "")
        if metric in lookup:
            raise MinimalCohortV1Error(f"Duplicate summary_metric detected in summary TSV: {metric}")
        lookup[metric] = row
    return lookup


def build_row_count_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        metric = str(row.get("count_metric") or "")
        if metric in lookup:
            raise MinimalCohortV1Error(f"Duplicate count_metric detected in row-count TSV: {metric}")
        lookup[metric] = row
    return lookup


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


def normalize_yes_no(value: bool) -> str:
    return "yes" if value else "no"


def required_field_names(
    rows: list[dict[str, str]],
    *,
    table_name: str | None = None,
    source_layer: str | None = None,
    proposed_role: str | None = None,
) -> list[str]:
    selected: list[str] = []
    for row in rows:
        if table_name is not None and str(row.get("table_name") or "") != table_name:
            continue
        if source_layer is not None and str(row.get("source_layer") or "") != source_layer:
            continue
        if proposed_role is not None and str(row.get("proposed_role") or "") != proposed_role:
            continue
        field_name = str(row.get("field_name") or "")
        if field_name:
            selected.append(field_name)
    return ordered_unique(selected)


def filter_rows(
    rows: list[dict[str, str]],
    *,
    table_name: str | None = None,
    source_layer: str | None = None,
    proposed_role: str | None = None,
) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    for row in rows:
        if table_name is not None and str(row.get("table_name") or "") != table_name:
            continue
        if source_layer is not None and str(row.get("source_layer") or "") != source_layer:
            continue
        if proposed_role is not None and str(row.get("proposed_role") or "") != proposed_role:
            continue
        selected.append(row)
    return selected


def adapt_dry_run_text(value: str) -> str:
    return (
        value.replace("minimal dry-run build", "minimal cohort v1 build")
        .replace("minimal dry-run cohort row build", "minimal cohort v1 row build")
        .replace("minimal dry-run cohort", "minimal cohort v1")
        .replace("dry-run join", "cohort v1 join")
        .replace("dry run", "cohort v1")
        .replace("dry-run", "cohort v1")
    )


def build_workflow_paths(dryrun_module: Any) -> WorkflowPaths:
    upstream_paths = dryrun_module.build_workflow_paths()
    trial_config_data = dryrun_module.load_yaml(upstream_paths.trial_config)

    processed_root = upstream_paths.repo_root / str(
        trial_config_data.get("processed_data_root", "01-data/processed")
    )
    results_root = upstream_paths.repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    audit_root = upstream_paths.repo_root / str(
        trial_config_data.get("audit_root", "01-data/audit")
    )

    return WorkflowPaths(
        repo_root=upstream_paths.repo_root,
        trial_config=upstream_paths.trial_config,
        results_root=results_root,
        processed_runs_root=processed_root / "tcga-brca" / "cohort" / "minimal_cohort_v1_runs",
        audit_runs_root=audit_root / "tcga-brca" / "cohort" / "minimal_cohort_v1_runs",
        latest_pointer=audit_root / "tcga-brca" / "cohort" / "tcga_brca_minimal_cohort_v1_latest.json",
        dry_run_latest_pointer=upstream_paths.latest_pointer,
        upstream_paths=upstream_paths,
    )


def extract_dry_run_contract(
    dry_run_row_counts_rows: list[dict[str, str]],
    dry_run_summary_rows: list[dict[str, str]],
) -> dict[str, Any]:
    row_count_lookup = build_row_count_lookup(dry_run_row_counts_rows)
    summary_lookup = build_summary_lookup(dry_run_summary_rows)

    return {
        "final_row_count": parse_int(
            row_count_lookup["final_cohort_row_count"]["count_value"],
            "dry-run final_cohort_row_count",
        ),
        "followup_matched_patient_count": parse_int(
            row_count_lookup["clinical_patient_rows_with_followup_match"]["count_value"],
            "dry-run clinical_patient_rows_with_followup_match",
        ),
        "followup_matched_row_count": parse_int(
            row_count_lookup["matched_followup_row_count"]["count_value"],
            "dry-run matched_followup_row_count",
        ),
        "biospecimen_matched_patient_count": parse_int(
            row_count_lookup["clinical_patient_rows_with_biospecimen_sample_match"]["count_value"],
            "dry-run clinical_patient_rows_with_biospecimen_sample_match",
        ),
        "biospecimen_matched_row_count": parse_int(
            row_count_lookup["matched_biospecimen_sample_row_count"]["count_value"],
            "dry-run matched_biospecimen_sample_row_count",
        ),
        "mismatch_signal": str(
            summary_lookup["clinical_patient_vs_source_vs_biospecimen_patient_counts"]["summary_value"]
        ),
    }


def validate_dry_run_contract(dry_run_contract: dict[str, Any]) -> None:
    for key, expected_value in EXPECTED_DRY_RUN_CONTRACT.items():
        observed_value = dry_run_contract[key]
        if observed_value != expected_value:
            raise MinimalCohortV1Error(
                "Current dry-run contract no longer matches the expected audited baseline: "
                f"{key}={observed_value!r} vs expected {expected_value!r}."
            )


def load_workflow_inputs(paths: WorkflowPaths, dryrun_module: Any) -> WorkflowInputs:
    upstream_inputs = dryrun_module.load_workflow_inputs(paths.upstream_paths)

    if not paths.dry_run_latest_pointer.exists():
        raise MinimalCohortV1Error(
            f"Required dry-run latest pointer not found: {paths.dry_run_latest_pointer}"
        )

    dry_run_latest_pointer = dryrun_module.load_json(paths.dry_run_latest_pointer)
    dryrun_module.require_keys(
        dry_run_latest_pointer,
        {
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
            "minimal_dry_run_cohort_tsv",
            "minimal_dry_run_join_audit_tsv",
            "minimal_dry_run_row_counts_tsv",
            "minimal_dry_run_exclusions_tsv",
            "minimal_dry_run_summary_tsv",
            "run_log_json",
        },
        "Minimal dry-run latest pointer",
        paths.dry_run_latest_pointer,
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
        observed_value = str(dry_run_latest_pointer[key])
        if observed_value != expected_value:
            raise MinimalCohortV1Error(
                "Dry-run latest pointer is out of sync with the current upstream latest pointers: "
                f"{key}={observed_value!r} vs expected {expected_value!r}."
            )

    input_paths = {
        **upstream_inputs.input_paths,
        "dry_run_cohort_tsv": dryrun_module.resolve_existing_path(
            paths.repo_root,
            str(dry_run_latest_pointer["minimal_dry_run_cohort_tsv"]),
            "minimal dry-run cohort TSV",
        ),
        "dry_run_join_audit_tsv": dryrun_module.resolve_existing_path(
            paths.repo_root,
            str(dry_run_latest_pointer["minimal_dry_run_join_audit_tsv"]),
            "minimal dry-run join-audit TSV",
        ),
        "dry_run_row_counts_tsv": dryrun_module.resolve_existing_path(
            paths.repo_root,
            str(dry_run_latest_pointer["minimal_dry_run_row_counts_tsv"]),
            "minimal dry-run row-count TSV",
        ),
        "dry_run_exclusions_tsv": dryrun_module.resolve_existing_path(
            paths.repo_root,
            str(dry_run_latest_pointer["minimal_dry_run_exclusions_tsv"]),
            "minimal dry-run exclusions TSV",
        ),
        "dry_run_summary_tsv": dryrun_module.resolve_existing_path(
            paths.repo_root,
            str(dry_run_latest_pointer["minimal_dry_run_summary_tsv"]),
            "minimal dry-run summary TSV",
        ),
        "dry_run_run_log_json": dryrun_module.resolve_existing_path(
            paths.repo_root,
            str(dry_run_latest_pointer["run_log_json"]),
            "minimal dry-run run log",
        ),
    }

    dry_run_run_log = dryrun_module.load_json(input_paths["dry_run_run_log_json"])
    if dry_run_run_log.get("status") != "completed":
        raise MinimalCohortV1Error("The referenced minimal dry-run run log is not completed.")
    if not bool(dry_run_run_log.get("validation", {}).get("passed", False)):
        raise MinimalCohortV1Error(
            "The referenced minimal dry-run run log does not report validation.passed == true."
        )

    dry_run_join_audit_rows = dryrun_module.read_tsv_dict_rows(input_paths["dry_run_join_audit_tsv"])
    dry_run_row_counts_rows = dryrun_module.read_tsv_dict_rows(input_paths["dry_run_row_counts_tsv"])
    dry_run_exclusions_rows = dryrun_module.read_tsv_dict_rows(input_paths["dry_run_exclusions_tsv"])
    dry_run_summary_rows = dryrun_module.read_tsv_dict_rows(input_paths["dry_run_summary_tsv"])

    dry_run_contract = extract_dry_run_contract(dry_run_row_counts_rows, dry_run_summary_rows)
    validate_dry_run_contract(dry_run_contract)

    return WorkflowInputs(
        upstream_inputs=upstream_inputs,
        dry_run_latest_pointer=dry_run_latest_pointer,
        dry_run_run_log=dry_run_run_log,
        dry_run_join_audit_rows=dry_run_join_audit_rows,
        dry_run_row_counts_rows=dry_run_row_counts_rows,
        dry_run_exclusions_rows=dry_run_exclusions_rows,
        dry_run_summary_rows=dry_run_summary_rows,
        dry_run_contract=dry_run_contract,
        input_paths=input_paths,
    )


def transform_cohort_rows(
    cohort_v1_build_id: str,
    dry_run_cohort_rows: list[dict[str, Any]],
    baseline_field_names: list[str],
    patient_endpoint_field_names: list[str],
    followup_candidate_field_names: list[str],
) -> list[dict[str, Any]]:
    transformed_rows: list[dict[str, Any]] = []
    for dry_run_row in dry_run_cohort_rows:
        cohort_row: dict[str, Any] = {
            "cohort_v1_build_id": cohort_v1_build_id,
            "provisional_patient_row_id": dry_run_row["provisional_patient_row_id"],
            "bcr_patient_barcode": dry_run_row["bcr_patient_barcode"],
            "bcr_patient_uuid": dry_run_row["bcr_patient_uuid"],
            "ambiguity_note_flags_json": dry_run_row["ambiguity_note_flags_json"],
            "join_status_flags_json": dry_run_row["join_status_flags_json"],
            "followup_match_row_count": dry_run_row["clinical_follow_up_v4_0__match_row_count"],
            "biospecimen_sample_match_row_count": dry_run_row["biospecimen_sample__match_row_count"],
        }

        for field_name in baseline_field_names:
            cohort_row[field_name] = dry_run_row[field_name]
        for field_name in patient_endpoint_field_names:
            cohort_row[f"clinical_patient__{field_name}"] = dry_run_row[
                f"clinical_patient__{field_name}"
            ]

        cohort_row["clinical_follow_up_v4_0__bcr_followup_barcode_json"] = dry_run_row[
            "clinical_follow_up_v4_0__bcr_followup_barcode_json"
        ]
        cohort_row["clinical_follow_up_v4_0__bcr_followup_uuid_json"] = dry_run_row[
            "clinical_follow_up_v4_0__bcr_followup_uuid_json"
        ]
        for field_name in followup_candidate_field_names:
            cohort_row[f"clinical_follow_up_v4_0__{field_name}_json"] = dry_run_row[
                f"clinical_follow_up_v4_0__{field_name}_json"
            ]

        cohort_row["biospecimen_sample__bcr_patient_uuid_json"] = dry_run_row[
            "biospecimen_sample__bcr_patient_uuid_json"
        ]
        cohort_row["biospecimen_sample__bcr_sample_barcode_json"] = dry_run_row[
            "biospecimen_sample__bcr_sample_barcode_json"
        ]
        cohort_row["biospecimen_sample__bcr_sample_uuid_json"] = dry_run_row[
            "biospecimen_sample__bcr_sample_uuid_json"
        ]
        transformed_rows.append(cohort_row)

    return transformed_rows


def build_spec_rows(
    workflow_inputs: WorkflowInputs,
    baseline_field_names: list[str],
    patient_endpoint_field_names: list[str],
    followup_candidate_field_names: list[str],
    derived: dict[str, Any],
) -> list[dict[str, Any]]:
    extras = derived["extras"]
    mismatch_flag = next(
        flag
        for flag in extras["ambiguity_note_flags"]
        if str(flag).startswith("case_count_mismatch_carried_forward_")
    )
    biospecimen_anchor_output_fields = [
        "biospecimen_sample__bcr_patient_uuid_json",
        "biospecimen_sample__bcr_sample_barcode_json",
        "biospecimen_sample__bcr_sample_uuid_json",
    ]
    sample_helper_rows = filter_rows(
        workflow_inputs.upstream_inputs.blueprint_optional_rows,
        table_name="biospecimen_sample",
        source_layer="biospecimen",
        proposed_role="biospecimen_sample_helper",
    )
    treatment_proxy_rows = [
        *filter_rows(
            workflow_inputs.upstream_inputs.blueprint_optional_rows,
            table_name="clinical_patient",
            source_layer="clinical",
            proposed_role="treatment_proxy_field",
        ),
        *filter_rows(
            workflow_inputs.upstream_inputs.blueprint_optional_rows,
            table_name="clinical_follow_up_v4_0",
            source_layer="clinical",
            proposed_role="treatment_proxy_field",
        ),
        *filter_rows(
            workflow_inputs.upstream_inputs.blueprint_optional_rows,
            table_name="clinical_drug",
            source_layer="clinical",
            proposed_role="treatment_proxy_field",
        ),
        *filter_rows(
            workflow_inputs.upstream_inputs.blueprint_optional_rows,
            table_name="clinical_radiation",
            source_layer="clinical",
            proposed_role="treatment_proxy_field",
        ),
    ]
    sparse_timing_rows = filter_rows(
        workflow_inputs.upstream_inputs.blueprint_deferred_rows,
        source_layer="endpoint",
        proposed_role="deferred_endpoint_candidate",
    )
    child_biospecimen_rows = filter_rows(
        workflow_inputs.upstream_inputs.blueprint_deferred_rows,
        source_layer="biospecimen",
    )

    rows: list[dict[str, Any]] = []

    def add_spec_row(
        field_name: str,
        source_table: str,
        field_category: str,
        include_in_v1: bool,
        inclusion_rule: str,
        ambiguity_flags: list[str],
    ) -> None:
        rows.append(
            {
                "field_name": field_name,
                "source_table": source_table,
                "field_category": field_category,
                "include_in_v1": normalize_yes_no(include_in_v1),
                "inclusion_rule": inclusion_rule,
                "ambiguity_carried_forward": json_list(ambiguity_flags),
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )

    derived_rows = [
        (
            "cohort_v1_build_id",
            "workflow_derived",
            "included_audit_field",
            True,
            "Populate every patient row with the UTC build id for this minimal cohort v1 run.",
            [],
        ),
        (
            "provisional_patient_row_id",
            "workflow_derived",
            "included_audit_field",
            True,
            "Assign 1-based row ids in clinical_patient source-row order and keep them provisional.",
            [],
        ),
        (
            "bcr_patient_barcode",
            "clinical_patient",
            "included_audit_field",
            True,
            "Preserve the clinical_patient barcode as a parallel patient identifier; do not collapse it into UUID.",
            ["parallel_patient_identifiers_retained"],
        ),
        (
            "bcr_patient_uuid",
            "clinical_patient",
            "included_audit_field",
            True,
            "Preserve the clinical_patient UUID as a parallel patient identifier; do not collapse it into barcode.",
            ["parallel_patient_identifiers_retained"],
        ),
        (
            "ambiguity_note_flags_json",
            "workflow_derived",
            "included_audit_field",
            True,
            "Carry global ambiguity flags into every patient row as a JSON array of provisional notes.",
            extras["ambiguity_note_flags"],
        ),
        (
            "join_status_flags_json",
            "workflow_derived",
            "included_audit_field",
            True,
            "Record join-state signals for follow-up and biospecimen sample attachment as a JSON array per row.",
            [],
        ),
        (
            "followup_match_row_count",
            "clinical_follow_up_v4_0",
            "included_audit_field",
            True,
            "Store the grouped count of exact (bcr_patient_uuid, bcr_patient_barcode) follow-up matches per patient row.",
            ["endpoint_overlap_fields_kept_side_by_side"],
        ),
        (
            "biospecimen_sample_match_row_count",
            "biospecimen_sample",
            "included_audit_field",
            True,
            "Store the grouped count of UUID-matched biospecimen_sample rows per patient row.",
            [
                "parallel_patient_identifiers_retained",
                mismatch_flag,
                "biospecimen_child_layer_expansion_deferred",
            ],
        ),
    ]
    for item in derived_rows:
        add_spec_row(*item)

    for field_name in baseline_field_names:
        add_spec_row(
            field_name,
            "clinical_patient",
            "included_baseline_field",
            True,
            "Include audited usable baseline covariates from clinical_patient exactly as accepted in the current blueprint required set.",
            [],
        )
    for field_name in patient_endpoint_field_names:
        add_spec_row(
            f"clinical_patient__{field_name}",
            "clinical_patient",
            "included_endpoint_candidate_field",
            True,
            "Carry patient-layer endpoint candidate values side by side with source-specific names and no endpoint collapse.",
            ["endpoint_overlap_fields_kept_side_by_side"],
        )

    add_spec_row(
        "clinical_follow_up_v4_0__bcr_followup_barcode_json",
        "clinical_follow_up_v4_0",
        "included_followup_join_evidence_field",
        True,
        "Retain grouped follow-up barcodes as JSON-array join evidence for exact-key matched patient rows.",
        [],
    )
    add_spec_row(
        "clinical_follow_up_v4_0__bcr_followup_uuid_json",
        "clinical_follow_up_v4_0",
        "included_followup_join_evidence_field",
        True,
        "Retain grouped follow-up UUIDs as JSON-array join evidence for exact-key matched patient rows.",
        [],
    )
    for field_name in followup_candidate_field_names:
        add_spec_row(
            f"clinical_follow_up_v4_0__{field_name}_json",
            "clinical_follow_up_v4_0",
            "included_endpoint_candidate_field",
            True,
            "Carry grouped follow-up endpoint candidate values as JSON arrays and keep them source-specific.",
            ["endpoint_overlap_fields_kept_side_by_side"],
        )

    for field_name in biospecimen_anchor_output_fields:
        add_spec_row(
            field_name,
            "biospecimen_sample",
            "included_biospecimen_sample_anchor_field",
            True,
            "Retain grouped biospecimen_sample anchor evidence as JSON arrays without expanding beyond the sample layer.",
            [
                "parallel_patient_identifiers_retained",
                mismatch_flag,
                "biospecimen_child_layer_expansion_deferred",
            ],
        )

    for row in sample_helper_rows:
        add_spec_row(
            str(row["field_name"]),
            "biospecimen_sample",
            "excluded_biospecimen_sample_helper_field",
            False,
            "Exclude optional biospecimen sample-helper fields because minimal cohort v1 keeps only the sample anchor evidence layer.",
            [],
        )
    for row in treatment_proxy_rows:
        add_spec_row(
            str(row["field_name"]),
            str(row["table_name"]),
            "excluded_treatment_field",
            False,
            "Exclude treatment proxy/detail fields from minimal cohort v1 to avoid unsupported patient-level treatment harmonization or aggregation.",
            ["treatment_detail_rows_excluded_from_minimal_build"],
        )
    for row in sparse_timing_rows:
        add_spec_row(
            str(row["field_name"]),
            str(row["table_name"]),
            "excluded_sparse_timing_field",
            False,
            "Exclude sparse death and progression timing fields from minimal cohort v1 and carry them forward only as unresolved endpoint audit evidence.",
            ["sparse_timing_fields_excluded_from_minimal_build"],
        )
    for row in child_biospecimen_rows:
        add_spec_row(
            str(row["field_name"]),
            str(row["table_name"]),
            "excluded_child_or_side_biospecimen_field",
            False,
            "Exclude deferred non-sample biospecimen fields from minimal cohort v1 because child and side-layer expansion remains out of scope.",
            ["biospecimen_child_layer_expansion_deferred"],
        )

    return rows


def build_join_audit_rows(
    cohort_v1_build_id: str,
    workflow_inputs: WorkflowInputs,
    derived: dict[str, Any],
    dryrun_module: Any,
) -> list[dict[str, Any]]:
    base_rows = dryrun_module.build_join_audit_rows(
        cohort_v1_build_id,
        workflow_inputs.upstream_inputs,
        derived,
    )
    transformed_rows: list[dict[str, Any]] = []
    for row in base_rows:
        transformed_rows.append(
            {
                "cohort_v1_build_id": row["dry_run_build_id"],
                "join_step_order": row["join_step_order"],
                "join_step_name": row["join_step_name"],
                "left_table_name": row["left_table_name"],
                "right_table_name": row["right_table_name"],
                "starting_left_row_count": row["starting_left_row_count"],
                "matched_left_row_count": row["matched_left_row_count"],
                "unmatched_left_row_count": row["unmatched_left_row_count"],
                "matched_right_row_count": row["matched_right_row_count"],
                "right_only_row_count": row["right_only_row_count"],
                "left_rows_with_multiple_matches": row["left_rows_with_multiple_matches"],
                "max_right_matches_per_left_row": row["max_right_matches_per_left_row"],
                "join_rule_used": adapt_dry_run_text(str(row["join_rule_used"])),
                "included_in_v1": "yes",
                "carried_forward_as_provisional": "yes",
                "excluded_from_v1": "no",
                "ambiguity_carried_forward": row["ambiguity_carried_forward"],
                "notes": adapt_dry_run_text(str(row["notes"])),
            }
        )
    return transformed_rows


def build_exclusions_rows(
    cohort_v1_build_id: str,
    workflow_inputs: WorkflowInputs,
    derived: dict[str, Any],
    dryrun_module: Any,
) -> list[dict[str, Any]]:
    base_rows = dryrun_module.build_exclusions_rows(
        cohort_v1_build_id,
        workflow_inputs.upstream_inputs,
        derived,
    )
    transformed_rows: list[dict[str, Any]] = []
    for row in base_rows:
        transformed_rows.append(
            {
                "cohort_v1_build_id": row["dry_run_build_id"],
                "exclusion_category": row["exclusion_category"],
                "scope_level": row["scope_level"],
                "source_table": row["source_table"],
                "entity_id": row["entity_id"],
                "entity_value": row["entity_value"],
                "reason": adapt_dry_run_text(str(row["reason"])),
                "count_value": row["count_value"],
                "details": adapt_dry_run_text(str(row["details"])),
            }
        )

    sample_helper_rows = filter_rows(
        workflow_inputs.upstream_inputs.blueprint_optional_rows,
        table_name="biospecimen_sample",
        source_layer="biospecimen",
        proposed_role="biospecimen_sample_helper",
    )
    for row in sample_helper_rows:
        transformed_rows.append(
            {
                "cohort_v1_build_id": cohort_v1_build_id,
                "exclusion_category": "excluded_biospecimen_sample_helper_field",
                "scope_level": "field",
                "source_table": "biospecimen_sample",
                "entity_id": "field_name",
                "entity_value": row["field_name"],
                "reason": (
                    "biospecimen_sample helper fields are excluded from minimal cohort v1 "
                    "because this build keeps only patient-linked sample anchor evidence."
                ),
                "count_value": 1,
                "details": "Excluded by v1 execution spec; no sample-helper columns are emitted into minimal_cohort_v1.tsv.",
            }
        )

    treatment_proxy_rows = [
        *filter_rows(
            workflow_inputs.upstream_inputs.blueprint_optional_rows,
            table_name="clinical_patient",
            source_layer="clinical",
            proposed_role="treatment_proxy_field",
        ),
        *filter_rows(
            workflow_inputs.upstream_inputs.blueprint_optional_rows,
            table_name="clinical_follow_up_v4_0",
            source_layer="clinical",
            proposed_role="treatment_proxy_field",
        ),
    ]
    for row in treatment_proxy_rows:
        transformed_rows.append(
            {
                "cohort_v1_build_id": cohort_v1_build_id,
                "exclusion_category": "excluded_treatment_proxy_field",
                "scope_level": "field",
                "source_table": row["table_name"],
                "entity_id": "field_name",
                "entity_value": row["field_name"],
                "reason": (
                    "Treatment proxy fields were left out of minimal cohort v1 so the build stays limited "
                    "to baseline fields, endpoint candidates, and biospecimen sample anchor evidence."
                ),
                "count_value": 1,
                "details": "Excluded by v1 execution spec; treatment harmonization remains out of scope.",
            }
        )

    for patient_uuid in derived["extras"]["extra_biospecimen_uuid_order"]:
        for row in derived["extras"]["extra_biospecimen_uuid_details"][patient_uuid]:
            transformed_rows.append(
                {
                    "cohort_v1_build_id": cohort_v1_build_id,
                    "exclusion_category": "unmatched_biospecimen_sample_row",
                    "scope_level": "sample_row",
                    "source_table": "biospecimen_sample",
                    "entity_id": "bcr_sample_barcode",
                    "entity_value": str(row.get("bcr_sample_barcode") or ""),
                    "reason": (
                        "This biospecimen_sample row was not attached because its patient UUID is absent "
                        "from the clinical_patient provisional patient universe."
                    ),
                    "count_value": 1,
                    "details": (
                        "bcr_patient_uuid="
                        + str(row.get("bcr_patient_uuid") or "")
                        + "; bcr_sample_uuid="
                        + str(row.get("bcr_sample_uuid") or "")
                    ),
                }
            )

    return transformed_rows


def build_summary_rows(
    cohort_v1_build_id: str,
    workflow_inputs: WorkflowInputs,
    baseline_field_names: list[str],
    patient_endpoint_field_names: list[str],
    followup_candidate_field_names: list[str],
    derived: dict[str, Any],
) -> list[dict[str, Any]]:
    metrics = derived["metrics"]
    extras = derived["extras"]
    ambiguity_summary_lookup = build_summary_lookup(
        workflow_inputs.upstream_inputs.ambiguity_summary_rows
    )
    sample_helper_fields = required_field_names(
        workflow_inputs.upstream_inputs.blueprint_optional_rows,
        table_name="biospecimen_sample",
        source_layer="biospecimen",
        proposed_role="biospecimen_sample_helper",
    )
    child_biospecimen_tables = ordered_unique(
        [
            str(row["table_name"])
            for row in filter_rows(
                workflow_inputs.upstream_inputs.blueprint_deferred_rows,
                source_layer="biospecimen",
            )
        ]
    )
    sparse_timing_fields = ordered_unique(
        [
            str(row["field_name"])
            for row in filter_rows(
                workflow_inputs.upstream_inputs.blueprint_deferred_rows,
                source_layer="endpoint",
                proposed_role="deferred_endpoint_candidate",
            )
        ]
    )
    mismatch_note_value = (
        f"{metrics['clinical_patient_count']} / "
        f"{metrics['clinical_source_case_coverage']} / "
        f"{metrics['biospecimen_patient_uuid_distinct_count']}"
    )
    biospecimen_anchor_fields = [
        "biospecimen_sample__bcr_patient_uuid_json",
        "biospecimen_sample__bcr_sample_barcode_json",
        "biospecimen_sample__bcr_sample_uuid_json",
    ]

    summary_rows = [
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "design",
            "summary_metric": "cohort_v1_status",
            "summary_value": "provisional_minimal_cohort_v1",
            "notes": "This build is the first real minimal cohort v1 and remains explicitly provisional.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "design",
            "summary_metric": "proposed_unit_of_analysis",
            "summary_value": "patient/case",
            "notes": "Minimal cohort v1 preserves patient/case as the row unit.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "design",
            "summary_metric": "provisional_patient_universe",
            "summary_value": "clinical_patient",
            "notes": "clinical_patient remains the provisional patient universe for cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "design",
            "summary_metric": "biospecimen_policy",
            "summary_value": "sample_anchor_only",
            "notes": "biospecimen is represented only through grouped biospecimen_sample anchor evidence.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "design",
            "summary_metric": "endpoint_policy",
            "summary_value": "candidate_columns_side_by_side_only",
            "notes": "Patient and follow-up endpoint-like fields remain side-by-side provenance-bearing candidates only.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "design",
            "summary_metric": "treatment_policy",
            "summary_value": "excluded_from_row_build",
            "notes": "Treatment proxy/detail fields remain outside the patient-level v1 row build.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "inputs",
            "summary_metric": "dry_run_build_id",
            "summary_value": workflow_inputs.dry_run_latest_pointer["dry_run_build_id"],
            "notes": "Current minimal dry-run build used as the execution contract reference for v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "inputs",
            "summary_metric": "blueprint_run_id",
            "summary_value": workflow_inputs.upstream_inputs.blueprint_latest_pointer["blueprint_run_id"],
            "notes": "Current cohort blueprint run used as the v1 field and join source.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "inputs",
            "summary_metric": "ambiguity_resolution_run_id",
            "summary_value": workflow_inputs.upstream_inputs.ambiguity_resolution_latest_pointer[
                "ambiguity_resolution_run_id"
            ],
            "notes": "Current ambiguity-resolution run used as the v1 containment source.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_baseline_field_count",
            "summary_value": len(baseline_field_names),
            "notes": "Baseline fields selected from cohort_blueprint_required_fields where table_name=clinical_patient, excluding the duplicated patient identifiers.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_patient_endpoint_candidate_field_count",
            "summary_value": len(patient_endpoint_field_names),
            "notes": "Patient-table endpoint candidate fields carried as explicit prefixed columns.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_followup_endpoint_candidate_field_count",
            "summary_value": len(followup_candidate_field_names),
            "notes": "Follow-up endpoint candidate fields carried as grouped JSON arrays.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_biospecimen_sample_anchor_field_count",
            "summary_value": len(biospecimen_anchor_fields),
            "notes": "Grouped biospecimen sample-anchor evidence fields carried into the v1 cohort table.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_baseline_fields_json",
            "summary_value": json_list(baseline_field_names),
            "notes": "Baseline clinical_patient fields included in minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_patient_endpoint_candidate_fields_json",
            "summary_value": json_list(
                [f"clinical_patient__{field_name}" for field_name in patient_endpoint_field_names]
            ),
            "notes": "Patient-layer endpoint candidate columns included in minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_followup_endpoint_candidate_fields_json",
            "summary_value": json_list(
                [
                    f"clinical_follow_up_v4_0__{field_name}_json"
                    for field_name in followup_candidate_field_names
                ]
            ),
            "notes": "Follow-up-layer endpoint candidate columns included in minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "field_selection",
            "summary_metric": "included_biospecimen_sample_anchor_fields_json",
            "summary_value": json_list(biospecimen_anchor_fields),
            "notes": "Biospecimen sample-anchor evidence columns included in minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "exclusions",
            "summary_metric": "excluded_sparse_timing_fields_json",
            "summary_value": json_list(sparse_timing_fields),
            "notes": "Sparse death and progression timing fields excluded from minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "exclusions",
            "summary_metric": "excluded_biospecimen_sample_helper_fields_json",
            "summary_value": json_list(sample_helper_fields),
            "notes": "Optional biospecimen sample-helper fields excluded from minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "exclusions",
            "summary_metric": "excluded_child_or_side_biospecimen_tables_json",
            "summary_value": json_list(child_biospecimen_tables),
            "notes": "Deferred non-sample biospecimen tables excluded from minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "row_counts",
            "summary_metric": "proposed_patient_universe_row_count",
            "summary_value": metrics["clinical_patient_count"],
            "notes": "clinical_patient source-row count used as the provisional patient universe.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "row_counts",
            "summary_metric": "final_v1_row_count",
            "summary_value": metrics["final_cohort_row_count"],
            "notes": "Final minimal cohort v1 row count.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_followup_candidate_count",
            "summary_value": metrics["followup_matched_patient_count"],
            "notes": "Clinical patient rows with at least one grouped follow-up evidence row.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_followup_row_count",
            "summary_value": metrics["followup_matched_right_row_count"],
            "notes": "Follow-up source rows grouped into matched patient-level evidence.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_biospecimen_sample_anchor_count",
            "summary_value": metrics["biospecimen_matched_patient_count"],
            "notes": "Clinical patient rows with at least one grouped biospecimen_sample evidence row.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_biospecimen_sample_row_count",
            "summary_value": metrics["biospecimen_matched_right_row_count"],
            "notes": "biospecimen_sample source rows grouped into matched patient-level evidence.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "identifier_signals",
            "summary_metric": "rows_with_both_patient_barcode_and_uuid",
            "summary_value": metrics["rows_with_both_ids"],
            "notes": "Clinical patient rows where both patient identifier forms are present.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "identifier_signals",
            "summary_metric": "rows_with_only_one_patient_identifier",
            "summary_value": metrics["rows_with_only_uuid"] + metrics["rows_with_only_barcode"],
            "notes": "Clinical patient rows with only one patient identifier form present.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "ambiguity",
            "summary_metric": "carried_global_ambiguity_note_count",
            "summary_value": metrics["carried_ambiguity_note_flag_count"],
            "notes": "Global ambiguity-note flags written into every minimal cohort v1 row.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "ambiguity",
            "summary_metric": "upstream_blueprint_ambiguity_count",
            "summary_value": len(workflow_inputs.upstream_inputs.blueprint_ambiguity_rows),
            "notes": "Current blueprint ambiguity count retained as context for minimal cohort v1.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "ambiguity",
            "summary_metric": "upstream_ambiguity_readiness_signal",
            "summary_value": ambiguity_summary_lookup["minimal_build_readiness_signal"]["summary_value"],
            "notes": "The upstream ambiguity-resolution layer still reports manual-review dependence for a final freeze.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "ambiguity",
            "summary_metric": "carried_ambiguity_note_flags_json",
            "summary_value": json_list(extras["ambiguity_note_flags"]),
            "notes": "Global ambiguity-note flags carried into every minimal cohort v1 row.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "clinical_patient_vs_source_vs_biospecimen_patient_counts",
            "summary_value": mismatch_note_value,
            "notes": "Minimal cohort v1 preserves the current count mismatch as an explicit audit note rather than resolving it.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "extra_source_case_id_count",
            "summary_value": metrics["extra_source_case_id_count"],
            "notes": "Case submitter IDs present in source metadata but not in clinical_patient.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "extra_biospecimen_patient_uuid_count",
            "summary_value": metrics["biospecimen_right_only_patient_uuid_count"],
            "notes": "biospecimen_sample patient UUIDs absent from clinical_patient.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "extra_biospecimen_sample_row_count",
            "summary_value": metrics["biospecimen_right_only_row_count"],
            "notes": "biospecimen_sample rows absent from the clinical_patient provisional patient universe.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "readiness",
            "summary_metric": "downstream_baseline_analysis_preparation_readiness",
            "summary_value": "ready_provisional_v1_only",
            "notes": "This build is suitable for downstream baseline-analysis preparation only while remaining explicitly provisional.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "readiness",
            "summary_metric": "endpoint_freeze_status",
            "summary_value": "blocked",
            "notes": "Endpoint freeze remains blocked by unresolved overlap and sparse timing ambiguities.",
        },
        {
            "cohort_v1_build_id": cohort_v1_build_id,
            "summary_section": "readiness",
            "summary_metric": "cohort_v1_readiness_interpretation",
            "summary_value": "sufficient_for_provisional_baseline_preparation_not_final_freeze",
            "notes": "The saved evidence is sufficient for a real provisional minimal cohort v1, but not for final cohort or endpoint freeze.",
        },
    ]

    if extras["extra_source_case_ids"]:
        summary_rows.append(
            {
                "cohort_v1_build_id": cohort_v1_build_id,
                "summary_section": "mismatch_signals",
                "summary_metric": "extra_source_case_ids_json",
                "summary_value": json_list(extras["extra_source_case_ids"]),
                "notes": "Source-only case submitter IDs carried forward as explicit audit notes.",
            }
        )
    if extras["extra_biospecimen_uuid_order"]:
        summary_rows.append(
            {
                "cohort_v1_build_id": cohort_v1_build_id,
                "summary_section": "mismatch_signals",
                "summary_metric": "extra_biospecimen_patient_uuids_json",
                "summary_value": json_list(extras["extra_biospecimen_uuid_order"]),
                "notes": "biospecimen-only patient UUIDs carried forward as explicit audit notes.",
            }
        )

    return summary_rows


def build_latest_pointer_payload(
    cohort_v1_build_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "cohort_v1_build_id": cohort_v1_build_id,
        "dry_run_build_id": str(workflow_inputs.dry_run_latest_pointer["dry_run_build_id"]),
        "blueprint_run_id": str(workflow_inputs.upstream_inputs.blueprint_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            workflow_inputs.upstream_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(workflow_inputs.upstream_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(
            workflow_inputs.upstream_inputs.clinical_shortlist_latest_pointer["core_audit_run_id"]
        ),
        "clinical_parse_run_id": str(workflow_inputs.upstream_inputs.clinical_biotab_latest_pointer["parse_run_id"]),
        "endpoint_crosswalk_run_id": str(
            workflow_inputs.upstream_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.upstream_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(
            workflow_inputs.upstream_inputs.biospecimen_biotab_latest_pointer["parse_run_id"]
        ),
        "clinical_source_run_id": str(workflow_inputs.upstream_inputs.clinical_biotab_latest_pointer["source_run_id"]),
        "biospecimen_source_run_id": str(
            workflow_inputs.upstream_inputs.biospecimen_biotab_latest_pointer["source_run_id"]
        ),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "minimal_cohort_v1_tsv": repo_relative(output_paths["cohort_tsv"], paths.repo_root),
        "minimal_cohort_v1_spec_tsv": repo_relative(output_paths["spec_tsv"], paths.repo_root),
        "minimal_cohort_v1_join_audit_tsv": repo_relative(output_paths["join_audit_tsv"], paths.repo_root),
        "minimal_cohort_v1_exclusions_tsv": repo_relative(output_paths["exclusions_tsv"], paths.repo_root),
        "minimal_cohort_v1_summary_tsv": repo_relative(output_paths["summary_tsv"], paths.repo_root),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "minimal_dry_run_latest_json": repo_relative(paths.dry_run_latest_pointer, paths.repo_root),
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


def write_failure_log(path: Path, payload: dict[str, Any], dryrun_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dryrun_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    cohort_v1_build_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    dryrun_module = load_dry_run_module()
    paths = build_workflow_paths(dryrun_module)
    processed_run_dir = paths.processed_runs_root / cohort_v1_build_id
    audit_run_dir = paths.audit_runs_root / cohort_v1_build_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = dryrun_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, dryrun_module)

        dryrun_module.create_run_directory(processed_run_dir)
        dryrun_module.create_run_directory(audit_run_dir)

        cohort_path = processed_run_dir / "minimal_cohort_v1.tsv"
        spec_path = audit_run_dir / "minimal_cohort_v1_spec.tsv"
        join_audit_path = audit_run_dir / "minimal_cohort_v1_join_audit.tsv"
        exclusions_path = audit_run_dir / "minimal_cohort_v1_exclusions.tsv"
        summary_path = audit_run_dir / "minimal_cohort_v1_summary.tsv"

        (
            dry_run_cohort_rows,
            baseline_field_names,
            patient_endpoint_field_names,
            followup_candidate_field_names,
            derived,
        ) = dryrun_module.build_cohort_rows(
            dry_run_build_id=cohort_v1_build_id,
            workflow_inputs=workflow_inputs.upstream_inputs,
        )
        if not dry_run_cohort_rows:
            raise MinimalCohortV1Error("Minimal cohort v1 rows were not generated.")

        cohort_rows = transform_cohort_rows(
            cohort_v1_build_id=cohort_v1_build_id,
            dry_run_cohort_rows=dry_run_cohort_rows,
            baseline_field_names=baseline_field_names,
            patient_endpoint_field_names=patient_endpoint_field_names,
            followup_candidate_field_names=followup_candidate_field_names,
        )
        spec_rows = build_spec_rows(
            workflow_inputs=workflow_inputs,
            baseline_field_names=baseline_field_names,
            patient_endpoint_field_names=patient_endpoint_field_names,
            followup_candidate_field_names=followup_candidate_field_names,
            derived=derived,
        )
        join_audit_rows = build_join_audit_rows(
            cohort_v1_build_id=cohort_v1_build_id,
            workflow_inputs=workflow_inputs,
            derived=derived,
            dryrun_module=dryrun_module,
        )
        exclusions_rows = build_exclusions_rows(
            cohort_v1_build_id=cohort_v1_build_id,
            workflow_inputs=workflow_inputs,
            derived=derived,
            dryrun_module=dryrun_module,
        )
        summary_rows = build_summary_rows(
            cohort_v1_build_id=cohort_v1_build_id,
            workflow_inputs=workflow_inputs,
            baseline_field_names=baseline_field_names,
            patient_endpoint_field_names=patient_endpoint_field_names,
            followup_candidate_field_names=followup_candidate_field_names,
            derived=derived,
        )

        if not spec_rows:
            raise MinimalCohortV1Error("Minimal cohort v1 spec rows were not generated.")
        if not join_audit_rows:
            raise MinimalCohortV1Error("Minimal cohort v1 join-audit rows were not generated.")
        if not exclusions_rows:
            raise MinimalCohortV1Error("Minimal cohort v1 exclusions rows were not generated.")
        if not summary_rows:
            raise MinimalCohortV1Error("Minimal cohort v1 summary rows were not generated.")

        cohort_fieldnames = list(cohort_rows[0].keys())
        dryrun_module.write_dict_rows_tsv(cohort_path, cohort_fieldnames, cohort_rows)
        dryrun_module.write_dict_rows_tsv(spec_path, SPEC_FIELDNAMES, spec_rows)
        dryrun_module.write_dict_rows_tsv(join_audit_path, JOIN_AUDIT_FIELDNAMES, join_audit_rows)
        dryrun_module.write_dict_rows_tsv(exclusions_path, EXCLUSIONS_FIELDNAMES, exclusions_rows)
        dryrun_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        summary_lookup = build_summary_lookup(summary_rows)
        computed_contract = {
            "final_row_count": parse_int(
                summary_lookup["final_v1_row_count"]["summary_value"],
                "v1 final_v1_row_count",
            ),
            "followup_matched_patient_count": parse_int(
                summary_lookup["matched_followup_candidate_count"]["summary_value"],
                "v1 matched_followup_candidate_count",
            ),
            "followup_matched_row_count": parse_int(
                summary_lookup["matched_followup_row_count"]["summary_value"],
                "v1 matched_followup_row_count",
            ),
            "biospecimen_matched_patient_count": parse_int(
                summary_lookup["matched_biospecimen_sample_anchor_count"]["summary_value"],
                "v1 matched_biospecimen_sample_anchor_count",
            ),
            "biospecimen_matched_row_count": parse_int(
                summary_lookup["matched_biospecimen_sample_row_count"]["summary_value"],
                "v1 matched_biospecimen_sample_row_count",
            ),
            "mismatch_signal": str(
                summary_lookup["clinical_patient_vs_source_vs_biospecimen_patient_counts"][
                    "summary_value"
                ]
            ),
        }
        if computed_contract != workflow_inputs.dry_run_contract:
            raise MinimalCohortV1Error(
                "Minimal cohort v1 contract does not match the current dry-run contract: "
                f"{computed_contract!r} vs {workflow_inputs.dry_run_contract!r}."
            )
        if computed_contract != EXPECTED_DRY_RUN_CONTRACT:
            raise MinimalCohortV1Error(
                "Minimal cohort v1 contract does not match the expected audited acceptance counts: "
                f"{computed_contract!r} vs {EXPECTED_DRY_RUN_CONTRACT!r}."
            )

        required_output_columns = [
            "cohort_v1_build_id",
            "provisional_patient_row_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "ambiguity_note_flags_json",
            "join_status_flags_json",
            "followup_match_row_count",
            "biospecimen_sample_match_row_count",
        ]
        required_output_columns.extend(baseline_field_names)
        required_output_columns.extend(
            [f"clinical_patient__{field_name}" for field_name in patient_endpoint_field_names]
        )
        required_output_columns.extend(
            [
                "clinical_follow_up_v4_0__bcr_followup_barcode_json",
                "clinical_follow_up_v4_0__bcr_followup_uuid_json",
            ]
        )
        required_output_columns.extend(
            [f"clinical_follow_up_v4_0__{field_name}_json" for field_name in followup_candidate_field_names]
        )
        required_output_columns.extend(
            [
                "biospecimen_sample__bcr_patient_uuid_json",
                "biospecimen_sample__bcr_sample_barcode_json",
                "biospecimen_sample__bcr_sample_uuid_json",
            ]
        )

        cohort_has_required_columns = all(
            column in cohort_fieldnames for column in required_output_columns
        )
        cohort_row_ids_are_sequential = [
            parse_int(row["provisional_patient_row_id"], "provisional_patient_row_id")
            for row in cohort_rows
        ] == list(range(1, len(cohort_rows) + 1))
        final_row_count_matches_patient_universe = (
            derived["metrics"]["final_cohort_row_count"] == derived["metrics"]["clinical_patient_count"]
        )
        join_audit_step_count_valid = len(join_audit_rows) == 4
        summary_contains_readiness = (
            summary_lookup["downstream_baseline_analysis_preparation_readiness"]["summary_value"]
            == "ready_provisional_v1_only"
        )
        summary_contains_endpoint_block = (
            summary_lookup["endpoint_freeze_status"]["summary_value"] == "blocked"
        )
        summary_contains_expected_mismatch = (
            summary_lookup["clinical_patient_vs_source_vs_biospecimen_patient_counts"]["summary_value"]
            == EXPECTED_DRY_RUN_CONTRACT["mismatch_signal"]
        )
        distinct_unmatched_sample_row_count = sum(
            row["exclusion_category"] == "unmatched_biospecimen_sample_row"
            for row in exclusions_rows
        )

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "cohort_tsv": cohort_path,
            "spec_tsv": spec_path,
            "join_audit_tsv": join_audit_path,
            "exclusions_tsv": exclusions_path,
            "summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            cohort_v1_build_id=cohort_v1_build_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "cohort_v1_build_id": cohort_v1_build_id,
            "dry_run_build_id": str(workflow_inputs.dry_run_latest_pointer["dry_run_build_id"]),
            "blueprint_run_id": str(workflow_inputs.upstream_inputs.blueprint_latest_pointer["blueprint_run_id"]),
            "ambiguity_resolution_run_id": str(
                workflow_inputs.upstream_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
            ),
            "shortlist_run_id": str(workflow_inputs.upstream_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(
                workflow_inputs.upstream_inputs.clinical_shortlist_latest_pointer["core_audit_run_id"]
            ),
            "clinical_parse_run_id": str(workflow_inputs.upstream_inputs.clinical_biotab_latest_pointer["parse_run_id"]),
            "endpoint_crosswalk_run_id": str(
                workflow_inputs.upstream_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
            ),
            "biospecimen_crosswalk_run_id": str(
                workflow_inputs.upstream_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
            ),
            "biospecimen_parse_run_id": str(
                workflow_inputs.upstream_inputs.biospecimen_biotab_latest_pointer["parse_run_id"]
            ),
            "clinical_source_run_id": str(workflow_inputs.upstream_inputs.clinical_biotab_latest_pointer["source_run_id"]),
            "biospecimen_source_run_id": str(
                workflow_inputs.upstream_inputs.biospecimen_biotab_latest_pointer["source_run_id"]
            ),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "minimal_dry_run_latest_json": repo_relative(paths.dry_run_latest_pointer, paths.repo_root),
                "cohort_blueprint_latest_json": repo_relative(paths.upstream_paths.blueprint_latest_pointer, paths.repo_root),
                "blueprint_ambiguity_resolution_latest_json": repo_relative(paths.upstream_paths.ambiguity_resolution_latest_pointer, paths.repo_root),
                "clinical_shortlist_latest_json": repo_relative(paths.upstream_paths.clinical_shortlist_latest_pointer, paths.repo_root),
                "endpoint_crosswalk_latest_json": repo_relative(paths.upstream_paths.endpoint_crosswalk_latest_pointer, paths.repo_root),
                "biospecimen_crosswalk_latest_json": repo_relative(paths.upstream_paths.biospecimen_crosswalk_latest_pointer, paths.repo_root),
                "clinical_biotab_latest_json": repo_relative(paths.upstream_paths.clinical_biotab_latest_pointer, paths.repo_root),
                "biospecimen_biotab_latest_json": repo_relative(paths.upstream_paths.biospecimen_biotab_latest_pointer, paths.repo_root),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "minimal_cohort_v1_tsv": repo_relative(cohort_path, paths.repo_root),
                "minimal_cohort_v1_spec_tsv": repo_relative(spec_path, paths.repo_root),
                "minimal_cohort_v1_join_audit_tsv": repo_relative(join_audit_path, paths.repo_root),
                "minimal_cohort_v1_exclusions_tsv": repo_relative(exclusions_path, paths.repo_root),
                "minimal_cohort_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    derived["metrics"]["clinical_patient_count"] > 0
                    and derived["metrics"]["final_cohort_row_count"] > 0
                    and final_row_count_matches_patient_universe
                    and cohort_has_required_columns
                    and cohort_row_ids_are_sequential
                    and join_audit_step_count_valid
                    and summary_contains_readiness
                    and summary_contains_endpoint_block
                    and summary_contains_expected_mismatch
                    and distinct_unmatched_sample_row_count == 9
                    and computed_contract == workflow_inputs.dry_run_contract
                ),
                "dry_run_latest_pointer_found": True,
                "dry_run_run_log_completed": True,
                "dry_run_validation_passed": True,
                "blueprint_latest_pointer_found": True,
                "ambiguity_resolution_latest_pointer_found": True,
                "clinical_shortlist_latest_pointer_found": True,
                "endpoint_crosswalk_latest_pointer_found": True,
                "biospecimen_crosswalk_latest_pointer_found": True,
                "clinical_biotab_latest_pointer_found": True,
                "biospecimen_biotab_latest_pointer_found": True,
                "required_source_tables_found": True,
                "processed_output_row_count_positive": derived["metrics"]["final_cohort_row_count"] > 0,
                "spec_row_count_positive": len(spec_rows) > 0,
                "join_audit_row_count_positive": len(join_audit_rows) > 0,
                "exclusions_row_count_positive": len(exclusions_rows) > 0,
                "summary_row_count_positive": len(summary_rows) > 0,
                "cohort_has_required_columns": cohort_has_required_columns,
                "cohort_row_ids_are_sequential": cohort_row_ids_are_sequential,
                "final_row_count_matches_clinical_patient": final_row_count_matches_patient_universe,
                "join_audit_step_count_valid": join_audit_step_count_valid,
                "summary_contains_readiness": summary_contains_readiness,
                "summary_contains_endpoint_block": summary_contains_endpoint_block,
                "summary_contains_expected_mismatch": summary_contains_expected_mismatch,
                "distinct_unmatched_biospecimen_sample_row_count_is_9": distinct_unmatched_sample_row_count == 9,
                "contract_matches_current_dry_run": computed_contract == workflow_inputs.dry_run_contract,
                "contract_matches_expected_acceptance_counts": computed_contract == EXPECTED_DRY_RUN_CONTRACT,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "unit_of_analysis": "patient/case",
                "provisional_patient_universe": "clinical_patient",
                "endpoint_policy": "candidate_columns_side_by_side_only",
                "treatment_policy": "exclude_one_to_many_treatment_detail",
                "biospecimen_policy": "sample_anchor_only",
                "no_canonical_patient_identifier": True,
                "no_endpoint_freeze": True,
                "no_modeling": True,
                "no_metabric": True,
                "no_raw_xml_or_ssf_parsing": True,
                "baseline_fields_selected_from_blueprint_required_fields": True,
                "biospecimen_sample_helper_fields_excluded": True,
            },
            "counts": {
                **derived["metrics"],
                "baseline_field_count": len(baseline_field_names),
                "patient_endpoint_candidate_field_count": len(patient_endpoint_field_names),
                "followup_endpoint_candidate_field_count": len(followup_candidate_field_names),
                "biospecimen_sample_anchor_field_count": 3,
                "spec_row_count": len(spec_rows),
                "join_audit_row_count": len(join_audit_rows),
                "exclusions_row_count": len(exclusions_rows),
                "summary_row_count": len(summary_rows),
            },
            "dry_run_contract": workflow_inputs.dry_run_contract,
            "computed_v1_contract": computed_contract,
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "minimal_dry_run_latest_pointer": workflow_inputs.dry_run_latest_pointer,
                "minimal_dry_run_run_log": workflow_inputs.dry_run_run_log,
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
                "source_run_log": workflow_inputs.upstream_inputs.source_run_log,
            },
        }

        dryrun_module.write_json(run_log_path, run_log_payload)
        dryrun_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "cohort_v1_build_id": cohort_v1_build_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_minimal_cohort_v1",
        }
        write_failure_log(run_log_path, failure_payload, dryrun_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA minimal cohort v1 workflow complete.")
    print(f"Cohort v1 build ID: {run_log['cohort_v1_build_id']}")
    print(f"Dry-run build ID: {run_log['dry_run_build_id']}")
    print(f"Blueprint run ID: {run_log['blueprint_run_id']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"Cohort TSV: {run_log['outputs']['minimal_cohort_v1_tsv']}")
    print(f"Spec TSV: {run_log['outputs']['minimal_cohort_v1_spec_tsv']}")
    print(f"Join audit TSV: {run_log['outputs']['minimal_cohort_v1_join_audit_tsv']}")
    print(f"Exclusions TSV: {run_log['outputs']['minimal_cohort_v1_exclusions_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['minimal_cohort_v1_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
