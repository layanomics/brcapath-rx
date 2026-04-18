#!/usr/bin/env python
"""Build an auditable minimal TCGA-BRCA dry-run cohort from saved audit layers."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover
    yaml = None

JOIN_AUDIT_FIELDNAMES = [
    "dry_run_build_id",
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
    "ambiguity_carried_forward",
    "notes",
]
ROW_COUNT_FIELDNAMES = [
    "dry_run_build_id",
    "count_stage",
    "count_metric",
    "count_value",
    "notes",
]
EXCLUSIONS_FIELDNAMES = [
    "dry_run_build_id",
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
    "dry_run_build_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
STATIC_AMBIGUITY_NOTE_FLAGS = [
    "parallel_patient_identifiers_retained",
    "endpoint_overlap_fields_kept_side_by_side",
    "sparse_timing_fields_excluded_from_minimal_build",
    "treatment_detail_rows_excluded_from_minimal_build",
    "biospecimen_child_layer_expansion_deferred",
]


class MinimalDryRunCohortError(RuntimeError):
    """Raised when the minimal dry-run cohort workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the minimal dry-run cohort workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    blueprint_latest_pointer: Path
    ambiguity_resolution_latest_pointer: Path
    clinical_shortlist_latest_pointer: Path
    endpoint_crosswalk_latest_pointer: Path
    biospecimen_crosswalk_latest_pointer: Path
    clinical_biotab_latest_pointer: Path
    biospecimen_biotab_latest_pointer: Path
    dry_run_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved minimal dry-run inputs loaded from saved audit layers."""

    blueprint_latest_pointer: dict[str, Any]
    blueprint_run_log: dict[str, Any]
    ambiguity_resolution_latest_pointer: dict[str, Any]
    ambiguity_resolution_run_log: dict[str, Any]
    clinical_shortlist_latest_pointer: dict[str, Any]
    clinical_shortlist_run_log: dict[str, Any]
    endpoint_crosswalk_latest_pointer: dict[str, Any]
    endpoint_crosswalk_run_log: dict[str, Any]
    biospecimen_crosswalk_latest_pointer: dict[str, Any]
    biospecimen_crosswalk_run_log: dict[str, Any]
    clinical_biotab_latest_pointer: dict[str, Any]
    clinical_biotab_run_log: dict[str, Any]
    biospecimen_biotab_latest_pointer: dict[str, Any]
    biospecimen_biotab_run_log: dict[str, Any]
    source_run_log: dict[str, Any]
    blueprint_required_rows: list[dict[str, str]]
    blueprint_optional_rows: list[dict[str, str]]
    blueprint_deferred_rows: list[dict[str, str]]
    blueprint_join_path_rows: list[dict[str, str]]
    blueprint_ambiguity_rows: list[dict[str, str]]
    blueprint_summary_rows: list[dict[str, str]]
    ambiguity_inventory_rows: list[dict[str, str]]
    ambiguity_priority_rows: list[dict[str, str]]
    ambiguity_actions_rows: list[dict[str, str]]
    ambiguity_summary_rows: list[dict[str, str]]
    clinical_shortlist_rows: list[dict[str, str]]
    clinical_shortlist_summary_rows: list[dict[str, str]]
    endpoint_crosswalk_rows: list[dict[str, str]]
    endpoint_summary_rows: list[dict[str, str]]
    biospecimen_crosswalk_rows: list[dict[str, str]]
    biospecimen_summary_rows: list[dict[str, str]]
    source_metadata_rows: list[dict[str, str]]
    clinical_table_manifest_rows: list[dict[str, str]]
    biospecimen_table_manifest_rows: list[dict[str, str]]
    clinical_patient_rows: list[dict[str, str]]
    clinical_follow_up_rows: list[dict[str, str]]
    biospecimen_sample_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise MinimalDryRunCohortError("Unable to locate the repository root from the script path.")


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def parse_yaml_scalar(value: str) -> str:
    value = value.strip()
    if value.startswith(("'", '"')) and value.endswith(("'", '"')) and len(value) >= 2:
        return value[1:-1]
    return value


def parse_simple_yaml_mapping(raw_text: str, source_path: Path) -> dict[str, Any]:
    parsed: dict[str, Any] = {}
    current_key: str | None = None

    for line_number, raw_line in enumerate(raw_text.splitlines(), start=1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        if raw_line.startswith((" ", "\t")):
            if current_key is None:
                raise MinimalDryRunCohortError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise MinimalDryRunCohortError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise MinimalDryRunCohortError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise MinimalDryRunCohortError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise MinimalDryRunCohortError(
                f"Expected a top-level key/value pair in YAML: {source_path}:{line_number}"
            )
        key, value = raw_line.split(":", 1)
        key = key.strip()
        value = value.strip()
        current_key = key
        parsed[key] = {} if value == "" else parse_yaml_scalar(value)

    return parsed


def load_yaml(path: Path) -> dict[str, Any]:
    raw_text = path.read_text(encoding="utf-8")
    if yaml is not None:
        data = yaml.safe_load(raw_text)
    else:
        data = parse_simple_yaml_mapping(raw_text, path)
    if not isinstance(data, dict):
        raise MinimalDryRunCohortError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise MinimalDryRunCohortError(f"Expected a JSON object in file: {path}")
    return data


def read_tsv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise MinimalDryRunCohortError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise MinimalDryRunCohortError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise MinimalDryRunCohortError(f"Run directory already exists: {path}")
    path.mkdir(parents=False, exist_ok=False)
    return path


def parse_int(value: Any, label: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise MinimalDryRunCohortError(f"Expected integer-like value for {label}: {value!r}") from exc


def require_keys(payload: dict[str, Any], required_keys: set[str], label: str, source_path: Path) -> None:
    missing_keys = required_keys.difference(payload.keys())
    if missing_keys:
        raise MinimalDryRunCohortError(
            f"{label} is missing required keys {sorted(missing_keys)}: {source_path}"
        )


def resolve_existing_path(repo_root: Path, relative_path: str, label: str) -> Path:
    path = repo_root / relative_path
    if not path.exists():
        raise MinimalDryRunCohortError(f"Required {label} not found: {path}")
    return path


def find_single_row(rows: list[dict[str, str]], label: str, **conditions: str) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if all(str(row.get(key, "")) == str(expected) for key, expected in conditions.items())
    ]
    if len(matches) != 1:
        raise MinimalDryRunCohortError(
            f"Expected exactly one {label} row matching {conditions}, found {len(matches)}."
        )
    return matches[0]


def build_summary_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        metric = str(row.get("summary_metric") or "")
        if metric in lookup:
            raise MinimalDryRunCohortError(f"Duplicate summary_metric detected in summary TSV: {metric}")
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


def truthy_string(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def build_workflow_paths() -> WorkflowPaths:
    repo_root = detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    trial_config_data = load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    cohort_root = audit_root / "tcga-brca" / "cohort"
    variable_root = audit_root / "tcga-brca" / "variables"

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        blueprint_latest_pointer=cohort_root / "tcga_brca_cohort_blueprint_latest.json",
        ambiguity_resolution_latest_pointer=(
            cohort_root / "tcga_brca_blueprint_ambiguity_resolution_latest.json"
        ),
        clinical_shortlist_latest_pointer=variable_root / "tcga_brca_clinical_shortlist_latest.json",
        endpoint_crosswalk_latest_pointer=variable_root / "tcga_brca_endpoint_crosswalk_latest.json",
        biospecimen_crosswalk_latest_pointer=(
            variable_root / "tcga_brca_biospecimen_identifier_crosswalk_latest.json"
        ),
        clinical_biotab_latest_pointer=variable_root / "tcga_brca_clinical_biotabs_latest.json",
        biospecimen_biotab_latest_pointer=variable_root / "tcga_brca_biospecimen_biotabs_latest.json",
        dry_run_runs_root=cohort_root / "minimal_dry_run_runs",
        latest_pointer=cohort_root / "tcga_brca_minimal_dry_run_cohort_latest.json",
    )


def unique_case_submitter_ids(
    source_metadata_rows: list[dict[str, str]],
    source_class: str,
) -> list[str]:
    observed: list[str] = []
    seen: set[str] = set()
    for row in source_metadata_rows:
        if str(row.get("source_class") or "") != source_class:
            continue
        raw_case_ids = str(row.get("case_submitter_ids") or "[]")
        case_ids = json.loads(raw_case_ids)
        if not isinstance(case_ids, list):
            raise MinimalDryRunCohortError(
                f"Expected case_submitter_ids JSON list in source metadata row: {raw_case_ids!r}"
            )
        for case_id in case_ids:
            case_id_str = str(case_id or "")
            if not case_id_str or case_id_str in seen:
                continue
            seen.add(case_id_str)
            observed.append(case_id_str)
    if not observed:
        raise MinimalDryRunCohortError(
            f"No case_submitter_ids were recovered for source_class={source_class!r} from source metadata."
        )
    return observed


def table_row_count_from_manifest(rows: list[dict[str, str]], table_name: str) -> int:
    table_row = find_single_row(rows, "table manifest", table_name=table_name)
    return parse_int(table_row["row_count"], f"{table_name} manifest row_count")


def load_workflow_inputs(paths: WorkflowPaths) -> WorkflowInputs:
    for pointer_path, label in [
        (paths.blueprint_latest_pointer, "cohort blueprint latest pointer"),
        (paths.ambiguity_resolution_latest_pointer, "ambiguity-resolution latest pointer"),
        (paths.clinical_shortlist_latest_pointer, "clinical shortlist latest pointer"),
        (paths.endpoint_crosswalk_latest_pointer, "endpoint crosswalk latest pointer"),
        (paths.biospecimen_crosswalk_latest_pointer, "biospecimen crosswalk latest pointer"),
        (paths.clinical_biotab_latest_pointer, "clinical biotab latest pointer"),
        (paths.biospecimen_biotab_latest_pointer, "biospecimen biotab latest pointer"),
    ]:
        if not pointer_path.exists():
            raise MinimalDryRunCohortError(f"Required {label} not found: {pointer_path}")

    blueprint_latest_pointer = load_json(paths.blueprint_latest_pointer)
    ambiguity_resolution_latest_pointer = load_json(paths.ambiguity_resolution_latest_pointer)
    clinical_shortlist_latest_pointer = load_json(paths.clinical_shortlist_latest_pointer)
    endpoint_crosswalk_latest_pointer = load_json(paths.endpoint_crosswalk_latest_pointer)
    biospecimen_crosswalk_latest_pointer = load_json(paths.biospecimen_crosswalk_latest_pointer)
    clinical_biotab_latest_pointer = load_json(paths.clinical_biotab_latest_pointer)
    biospecimen_biotab_latest_pointer = load_json(paths.biospecimen_biotab_latest_pointer)

    require_keys(
        blueprint_latest_pointer,
        {
            "blueprint_run_id",
            "shortlist_run_id",
            "core_audit_run_id",
            "clinical_parse_run_id",
            "endpoint_crosswalk_run_id",
            "biospecimen_crosswalk_run_id",
            "biospecimen_parse_run_id",
            "clinical_source_run_id",
            "biospecimen_source_run_id",
            "cohort_blueprint_required_fields_tsv",
            "cohort_blueprint_optional_fields_tsv",
            "cohort_blueprint_deferred_fields_tsv",
            "cohort_blueprint_join_path_tsv",
            "cohort_blueprint_ambiguities_tsv",
            "cohort_blueprint_summary_tsv",
            "run_log_json",
        },
        "Cohort blueprint latest pointer",
        paths.blueprint_latest_pointer,
    )
    require_keys(
        ambiguity_resolution_latest_pointer,
        {
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
            "ambiguity_resolution_inventory_tsv",
            "ambiguity_resolution_priority_tsv",
            "ambiguity_resolution_actions_tsv",
            "ambiguity_resolution_summary_tsv",
            "run_log_json",
        },
        "Ambiguity-resolution latest pointer",
        paths.ambiguity_resolution_latest_pointer,
    )
    require_keys(
        clinical_shortlist_latest_pointer,
        {
            "shortlist_run_id",
            "core_audit_run_id",
            "parse_run_id",
            "source_run_id",
            "clinical_shortlist_tsv",
            "clinical_shortlist_summary_tsv",
            "clinical_shortlist_by_table_tsv",
            "run_log_json",
        },
        "Clinical shortlist latest pointer",
        paths.clinical_shortlist_latest_pointer,
    )
    require_keys(
        endpoint_crosswalk_latest_pointer,
        {
            "crosswalk_run_id",
            "shortlist_run_id",
            "core_audit_run_id",
            "parse_run_id",
            "source_run_id",
            "endpoint_candidate_inventory_tsv",
            "endpoint_crosswalk_tsv",
            "endpoint_crosswalk_summary_tsv",
            "run_log_json",
        },
        "Endpoint crosswalk latest pointer",
        paths.endpoint_crosswalk_latest_pointer,
    )
    require_keys(
        biospecimen_crosswalk_latest_pointer,
        {
            "crosswalk_run_id",
            "parse_run_id",
            "source_run_id",
            "biospecimen_identifier_inventory_tsv",
            "biospecimen_identifier_crosswalk_tsv",
            "biospecimen_identifier_crosswalk_summary_tsv",
            "run_log_json",
        },
        "Biospecimen crosswalk latest pointer",
        paths.biospecimen_crosswalk_latest_pointer,
    )
    require_keys(
        clinical_biotab_latest_pointer,
        {
            "parse_run_id",
            "source_run_id",
            "processed_run_directory",
            "table_manifest_tsv",
            "run_log_json",
            "source_metadata_tsv",
        },
        "Clinical biotab latest pointer",
        paths.clinical_biotab_latest_pointer,
    )
    require_keys(
        biospecimen_biotab_latest_pointer,
        {
            "parse_run_id",
            "source_run_id",
            "processed_run_directory",
            "table_manifest_tsv",
            "run_log_json",
            "source_metadata_tsv",
        },
        "Biospecimen biotab latest pointer",
        paths.biospecimen_biotab_latest_pointer,
    )

    if str(blueprint_latest_pointer["shortlist_run_id"]) != str(
        clinical_shortlist_latest_pointer["shortlist_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Cohort blueprint latest pointer does not reference the current clinical shortlist run."
        )
    if str(blueprint_latest_pointer["endpoint_crosswalk_run_id"]) != str(
        endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Cohort blueprint latest pointer does not reference the current endpoint crosswalk run."
        )
    if str(blueprint_latest_pointer["biospecimen_crosswalk_run_id"]) != str(
        biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Cohort blueprint latest pointer does not reference the current biospecimen crosswalk run."
        )
    if str(ambiguity_resolution_latest_pointer["blueprint_run_id"]) != str(
        blueprint_latest_pointer["blueprint_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Ambiguity-resolution latest pointer does not reference the current cohort blueprint run."
        )
    if str(ambiguity_resolution_latest_pointer["shortlist_run_id"]) != str(
        clinical_shortlist_latest_pointer["shortlist_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Ambiguity-resolution latest pointer does not reference the current clinical shortlist run."
        )
    if str(ambiguity_resolution_latest_pointer["endpoint_crosswalk_run_id"]) != str(
        endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Ambiguity-resolution latest pointer does not reference the current endpoint crosswalk run."
        )
    if str(ambiguity_resolution_latest_pointer["biospecimen_crosswalk_run_id"]) != str(
        biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Ambiguity-resolution latest pointer does not reference the current biospecimen crosswalk run."
        )

    for label, left_value, right_value in [
        (
            "clinical shortlist parse_run_id vs clinical biotab parse_run_id",
            clinical_shortlist_latest_pointer["parse_run_id"],
            clinical_biotab_latest_pointer["parse_run_id"],
        ),
        (
            "endpoint crosswalk parse_run_id vs clinical biotab parse_run_id",
            endpoint_crosswalk_latest_pointer["parse_run_id"],
            clinical_biotab_latest_pointer["parse_run_id"],
        ),
        (
            "cohort blueprint clinical_parse_run_id vs clinical biotab parse_run_id",
            blueprint_latest_pointer["clinical_parse_run_id"],
            clinical_biotab_latest_pointer["parse_run_id"],
        ),
        (
            "ambiguity-resolution clinical_parse_run_id vs clinical biotab parse_run_id",
            ambiguity_resolution_latest_pointer["clinical_parse_run_id"],
            clinical_biotab_latest_pointer["parse_run_id"],
        ),
        (
            "biospecimen crosswalk parse_run_id vs biospecimen biotab parse_run_id",
            biospecimen_crosswalk_latest_pointer["parse_run_id"],
            biospecimen_biotab_latest_pointer["parse_run_id"],
        ),
        (
            "cohort blueprint biospecimen_parse_run_id vs biospecimen biotab parse_run_id",
            blueprint_latest_pointer["biospecimen_parse_run_id"],
            biospecimen_biotab_latest_pointer["parse_run_id"],
        ),
        (
            "ambiguity-resolution biospecimen_parse_run_id vs biospecimen biotab parse_run_id",
            ambiguity_resolution_latest_pointer["biospecimen_parse_run_id"],
            biospecimen_biotab_latest_pointer["parse_run_id"],
        ),
        (
            "clinical shortlist source_run_id vs clinical biotab source_run_id",
            clinical_shortlist_latest_pointer["source_run_id"],
            clinical_biotab_latest_pointer["source_run_id"],
        ),
        (
            "endpoint crosswalk source_run_id vs clinical biotab source_run_id",
            endpoint_crosswalk_latest_pointer["source_run_id"],
            clinical_biotab_latest_pointer["source_run_id"],
        ),
        (
            "cohort blueprint clinical_source_run_id vs clinical biotab source_run_id",
            blueprint_latest_pointer["clinical_source_run_id"],
            clinical_biotab_latest_pointer["source_run_id"],
        ),
        (
            "ambiguity-resolution clinical_source_run_id vs clinical biotab source_run_id",
            ambiguity_resolution_latest_pointer["clinical_source_run_id"],
            clinical_biotab_latest_pointer["source_run_id"],
        ),
        (
            "biospecimen crosswalk source_run_id vs biospecimen biotab source_run_id",
            biospecimen_crosswalk_latest_pointer["source_run_id"],
            biospecimen_biotab_latest_pointer["source_run_id"],
        ),
        (
            "cohort blueprint biospecimen_source_run_id vs biospecimen biotab source_run_id",
            blueprint_latest_pointer["biospecimen_source_run_id"],
            biospecimen_biotab_latest_pointer["source_run_id"],
        ),
        (
            "ambiguity-resolution biospecimen_source_run_id vs biospecimen biotab source_run_id",
            ambiguity_resolution_latest_pointer["biospecimen_source_run_id"],
            biospecimen_biotab_latest_pointer["source_run_id"],
        ),
    ]:
        if str(left_value) != str(right_value):
            raise MinimalDryRunCohortError(
                f"Upstream run-id mismatch for {label}: {left_value} vs {right_value}"
            )

    clinical_processed_run_directory = resolve_existing_path(
        paths.repo_root,
        str(clinical_biotab_latest_pointer["processed_run_directory"]),
        "clinical processed run directory",
    )
    biospecimen_processed_run_directory = resolve_existing_path(
        paths.repo_root,
        str(biospecimen_biotab_latest_pointer["processed_run_directory"]),
        "biospecimen processed run directory",
    )

    input_paths = {
        "blueprint_required_tsv": resolve_existing_path(
            paths.repo_root,
            str(blueprint_latest_pointer["cohort_blueprint_required_fields_tsv"]),
            "cohort blueprint required fields TSV",
        ),
        "blueprint_optional_tsv": resolve_existing_path(
            paths.repo_root,
            str(blueprint_latest_pointer["cohort_blueprint_optional_fields_tsv"]),
            "cohort blueprint optional fields TSV",
        ),
        "blueprint_deferred_tsv": resolve_existing_path(
            paths.repo_root,
            str(blueprint_latest_pointer["cohort_blueprint_deferred_fields_tsv"]),
            "cohort blueprint deferred fields TSV",
        ),
        "blueprint_join_path_tsv": resolve_existing_path(
            paths.repo_root,
            str(blueprint_latest_pointer["cohort_blueprint_join_path_tsv"]),
            "cohort blueprint join-path TSV",
        ),
        "blueprint_ambiguities_tsv": resolve_existing_path(
            paths.repo_root,
            str(blueprint_latest_pointer["cohort_blueprint_ambiguities_tsv"]),
            "cohort blueprint ambiguities TSV",
        ),
        "blueprint_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(blueprint_latest_pointer["cohort_blueprint_summary_tsv"]),
            "cohort blueprint summary TSV",
        ),
        "blueprint_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(blueprint_latest_pointer["run_log_json"]),
            "cohort blueprint run log",
        ),
        "ambiguity_inventory_tsv": resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["ambiguity_resolution_inventory_tsv"]),
            "ambiguity-resolution inventory TSV",
        ),
        "ambiguity_priority_tsv": resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["ambiguity_resolution_priority_tsv"]),
            "ambiguity-resolution priority TSV",
        ),
        "ambiguity_actions_tsv": resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["ambiguity_resolution_actions_tsv"]),
            "ambiguity-resolution actions TSV",
        ),
        "ambiguity_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["ambiguity_resolution_summary_tsv"]),
            "ambiguity-resolution summary TSV",
        ),
        "ambiguity_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(ambiguity_resolution_latest_pointer["run_log_json"]),
            "ambiguity-resolution run log",
        ),
        "clinical_shortlist_tsv": resolve_existing_path(
            paths.repo_root,
            str(clinical_shortlist_latest_pointer["clinical_shortlist_tsv"]),
            "clinical shortlist TSV",
        ),
        "clinical_shortlist_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(clinical_shortlist_latest_pointer["clinical_shortlist_summary_tsv"]),
            "clinical shortlist summary TSV",
        ),
        "clinical_shortlist_by_table_tsv": resolve_existing_path(
            paths.repo_root,
            str(clinical_shortlist_latest_pointer["clinical_shortlist_by_table_tsv"]),
            "clinical shortlist by-table TSV",
        ),
        "clinical_shortlist_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(clinical_shortlist_latest_pointer["run_log_json"]),
            "clinical shortlist run log",
        ),
        "endpoint_inventory_tsv": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_candidate_inventory_tsv"]),
            "endpoint candidate inventory TSV",
        ),
        "endpoint_crosswalk_tsv": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_crosswalk_tsv"]),
            "endpoint crosswalk TSV",
        ),
        "endpoint_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_crosswalk_summary_tsv"]),
            "endpoint crosswalk summary TSV",
        ),
        "endpoint_crosswalk_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["run_log_json"]),
            "endpoint crosswalk run log",
        ),
        "biospecimen_inventory_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_inventory_tsv"]),
            "biospecimen inventory TSV",
        ),
        "biospecimen_crosswalk_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_crosswalk_tsv"]),
            "biospecimen crosswalk TSV",
        ),
        "biospecimen_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_crosswalk_summary_tsv"]),
            "biospecimen crosswalk summary TSV",
        ),
        "biospecimen_crosswalk_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["run_log_json"]),
            "biospecimen crosswalk run log",
        ),
        "clinical_biotab_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(clinical_biotab_latest_pointer["run_log_json"]),
            "clinical biotab run log",
        ),
        "clinical_table_manifest_tsv": resolve_existing_path(
            paths.repo_root,
            str(clinical_biotab_latest_pointer["table_manifest_tsv"]),
            "clinical biotab table manifest TSV",
        ),
        "source_metadata_tsv": resolve_existing_path(
            paths.repo_root,
            str(clinical_biotab_latest_pointer["source_metadata_tsv"]),
            "source metadata TSV",
        ),
        "biospecimen_biotab_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_biotab_latest_pointer["run_log_json"]),
            "biospecimen biotab run log",
        ),
        "biospecimen_table_manifest_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_biotab_latest_pointer["table_manifest_tsv"]),
            "biospecimen biotab table manifest TSV",
        ),
        "clinical_patient_tsv": resolve_existing_path(
            clinical_processed_run_directory,
            "clinical_patient.tsv",
            "clinical_patient TSV",
        ),
        "clinical_follow_up_tsv": resolve_existing_path(
            clinical_processed_run_directory,
            "clinical_follow_up_v4_0.tsv",
            "clinical_follow_up_v4_0 TSV",
        ),
        "biospecimen_sample_tsv": resolve_existing_path(
            biospecimen_processed_run_directory,
            "biospecimen_sample.tsv",
            "biospecimen_sample TSV",
        ),
    }

    blueprint_run_log = load_json(input_paths["blueprint_run_log_json"])
    ambiguity_resolution_run_log = load_json(input_paths["ambiguity_run_log_json"])
    clinical_shortlist_run_log = load_json(input_paths["clinical_shortlist_run_log_json"])
    endpoint_crosswalk_run_log = load_json(input_paths["endpoint_crosswalk_run_log_json"])
    biospecimen_crosswalk_run_log = load_json(input_paths["biospecimen_crosswalk_run_log_json"])
    clinical_biotab_run_log = load_json(input_paths["clinical_biotab_run_log_json"])
    biospecimen_biotab_run_log = load_json(input_paths["biospecimen_biotab_run_log_json"])
    source_run_log_path_raw = str(
        clinical_biotab_run_log.get("inputs", {}).get("source_run_log_json") or ""
    )
    if not source_run_log_path_raw:
        raise MinimalDryRunCohortError(
            "Clinical biotab run log does not reference a source supplement run log."
        )
    input_paths["source_run_log_json"] = resolve_existing_path(
        paths.repo_root,
        source_run_log_path_raw,
        "source supplement run log",
    )
    source_run_log = load_json(input_paths["source_run_log_json"])

    for run_log, label in [
        (blueprint_run_log, "cohort blueprint"),
        (ambiguity_resolution_run_log, "ambiguity-resolution"),
        (clinical_shortlist_run_log, "clinical shortlist"),
        (endpoint_crosswalk_run_log, "endpoint crosswalk"),
        (biospecimen_crosswalk_run_log, "biospecimen crosswalk"),
        (clinical_biotab_run_log, "clinical biotab"),
        (biospecimen_biotab_run_log, "biospecimen biotab"),
        (source_run_log, "source supplement"),
    ]:
        if run_log.get("status") != "completed":
            raise MinimalDryRunCohortError(f"Upstream {label} run log is not completed.")

    if str(blueprint_run_log.get("blueprint_run_id") or "") != str(
        blueprint_latest_pointer["blueprint_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Cohort blueprint run log ID does not match the latest blueprint pointer."
        )
    if str(ambiguity_resolution_run_log.get("ambiguity_resolution_run_id") or "") != str(
        ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Ambiguity-resolution run log ID does not match the latest ambiguity-resolution pointer."
        )
    if str(clinical_shortlist_run_log.get("shortlist_run_id") or "") != str(
        clinical_shortlist_latest_pointer["shortlist_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Clinical shortlist run log ID does not match the latest clinical shortlist pointer."
        )
    if str(endpoint_crosswalk_run_log.get("crosswalk_run_id") or "") != str(
        endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Endpoint crosswalk run log ID does not match the latest endpoint crosswalk pointer."
        )
    if str(biospecimen_crosswalk_run_log.get("crosswalk_run_id") or "") != str(
        biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Biospecimen crosswalk run log ID does not match the latest biospecimen crosswalk pointer."
        )
    if str(clinical_biotab_run_log.get("parse_run_id") or "") != str(
        clinical_biotab_latest_pointer["parse_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Clinical biotab run log parse_run_id does not match the latest clinical biotab pointer."
        )
    if str(biospecimen_biotab_run_log.get("parse_run_id") or "") != str(
        biospecimen_biotab_latest_pointer["parse_run_id"]
    ):
        raise MinimalDryRunCohortError(
            "Biospecimen biotab run log parse_run_id does not match the latest biospecimen biotab pointer."
        )

    blueprint_required_rows = read_tsv_dict_rows(input_paths["blueprint_required_tsv"])
    blueprint_optional_rows = read_tsv_dict_rows(input_paths["blueprint_optional_tsv"])
    blueprint_deferred_rows = read_tsv_dict_rows(input_paths["blueprint_deferred_tsv"])
    blueprint_join_path_rows = read_tsv_dict_rows(input_paths["blueprint_join_path_tsv"])
    blueprint_ambiguity_rows = read_tsv_dict_rows(input_paths["blueprint_ambiguities_tsv"])
    blueprint_summary_rows = read_tsv_dict_rows(input_paths["blueprint_summary_tsv"])
    ambiguity_inventory_rows = read_tsv_dict_rows(input_paths["ambiguity_inventory_tsv"])
    ambiguity_priority_rows = read_tsv_dict_rows(input_paths["ambiguity_priority_tsv"])
    ambiguity_actions_rows = read_tsv_dict_rows(input_paths["ambiguity_actions_tsv"])
    ambiguity_summary_rows = read_tsv_dict_rows(input_paths["ambiguity_summary_tsv"])
    clinical_shortlist_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_tsv"])
    clinical_shortlist_summary_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_summary_tsv"])
    endpoint_crosswalk_rows = read_tsv_dict_rows(input_paths["endpoint_crosswalk_tsv"])
    endpoint_summary_rows = read_tsv_dict_rows(input_paths["endpoint_summary_tsv"])
    biospecimen_crosswalk_rows = read_tsv_dict_rows(input_paths["biospecimen_crosswalk_tsv"])
    biospecimen_summary_rows = read_tsv_dict_rows(input_paths["biospecimen_summary_tsv"])
    source_metadata_rows = read_tsv_dict_rows(input_paths["source_metadata_tsv"])
    clinical_table_manifest_rows = read_tsv_dict_rows(input_paths["clinical_table_manifest_tsv"])
    biospecimen_table_manifest_rows = read_tsv_dict_rows(input_paths["biospecimen_table_manifest_tsv"])
    clinical_patient_rows = read_tsv_dict_rows(input_paths["clinical_patient_tsv"])
    clinical_follow_up_rows = read_tsv_dict_rows(input_paths["clinical_follow_up_tsv"])
    biospecimen_sample_rows = read_tsv_dict_rows(input_paths["biospecimen_sample_tsv"])

    for rows, label in [
        (blueprint_required_rows, "cohort blueprint required fields"),
        (blueprint_optional_rows, "cohort blueprint optional fields"),
        (blueprint_deferred_rows, "cohort blueprint deferred fields"),
        (blueprint_join_path_rows, "cohort blueprint join path"),
        (blueprint_ambiguity_rows, "cohort blueprint ambiguities"),
        (blueprint_summary_rows, "cohort blueprint summary"),
        (ambiguity_inventory_rows, "ambiguity-resolution inventory"),
        (ambiguity_priority_rows, "ambiguity-resolution priority"),
        (ambiguity_actions_rows, "ambiguity-resolution actions"),
        (ambiguity_summary_rows, "ambiguity-resolution summary"),
        (clinical_shortlist_rows, "clinical shortlist"),
        (clinical_shortlist_summary_rows, "clinical shortlist summary"),
        (endpoint_crosswalk_rows, "endpoint crosswalk"),
        (endpoint_summary_rows, "endpoint crosswalk summary"),
        (biospecimen_crosswalk_rows, "biospecimen crosswalk"),
        (biospecimen_summary_rows, "biospecimen crosswalk summary"),
        (source_metadata_rows, "source metadata"),
        (clinical_table_manifest_rows, "clinical table manifest"),
        (biospecimen_table_manifest_rows, "biospecimen table manifest"),
        (clinical_patient_rows, "clinical_patient table"),
        (clinical_follow_up_rows, "clinical_follow_up_v4_0 table"),
        (biospecimen_sample_rows, "biospecimen_sample table"),
    ]:
        if not rows:
            raise MinimalDryRunCohortError(f"Required input table has no rows: {label}")

    if table_row_count_from_manifest(clinical_table_manifest_rows, "clinical_patient") != len(
        clinical_patient_rows
    ):
        raise MinimalDryRunCohortError(
            "clinical_patient row count does not match the clinical table manifest."
        )
    if table_row_count_from_manifest(clinical_table_manifest_rows, "clinical_follow_up_v4_0") != len(
        clinical_follow_up_rows
    ):
        raise MinimalDryRunCohortError(
            "clinical_follow_up_v4_0 row count does not match the clinical table manifest."
        )
    if table_row_count_from_manifest(biospecimen_table_manifest_rows, "biospecimen_sample") != len(
        biospecimen_sample_rows
    ):
        raise MinimalDryRunCohortError(
            "biospecimen_sample row count does not match the biospecimen table manifest."
        )

    return WorkflowInputs(
        blueprint_latest_pointer=blueprint_latest_pointer,
        blueprint_run_log=blueprint_run_log,
        ambiguity_resolution_latest_pointer=ambiguity_resolution_latest_pointer,
        ambiguity_resolution_run_log=ambiguity_resolution_run_log,
        clinical_shortlist_latest_pointer=clinical_shortlist_latest_pointer,
        clinical_shortlist_run_log=clinical_shortlist_run_log,
        endpoint_crosswalk_latest_pointer=endpoint_crosswalk_latest_pointer,
        endpoint_crosswalk_run_log=endpoint_crosswalk_run_log,
        biospecimen_crosswalk_latest_pointer=biospecimen_crosswalk_latest_pointer,
        biospecimen_crosswalk_run_log=biospecimen_crosswalk_run_log,
        clinical_biotab_latest_pointer=clinical_biotab_latest_pointer,
        clinical_biotab_run_log=clinical_biotab_run_log,
        biospecimen_biotab_latest_pointer=biospecimen_biotab_latest_pointer,
        biospecimen_biotab_run_log=biospecimen_biotab_run_log,
        source_run_log=source_run_log,
        blueprint_required_rows=blueprint_required_rows,
        blueprint_optional_rows=blueprint_optional_rows,
        blueprint_deferred_rows=blueprint_deferred_rows,
        blueprint_join_path_rows=blueprint_join_path_rows,
        blueprint_ambiguity_rows=blueprint_ambiguity_rows,
        blueprint_summary_rows=blueprint_summary_rows,
        ambiguity_inventory_rows=ambiguity_inventory_rows,
        ambiguity_priority_rows=ambiguity_priority_rows,
        ambiguity_actions_rows=ambiguity_actions_rows,
        ambiguity_summary_rows=ambiguity_summary_rows,
        clinical_shortlist_rows=clinical_shortlist_rows,
        clinical_shortlist_summary_rows=clinical_shortlist_summary_rows,
        endpoint_crosswalk_rows=endpoint_crosswalk_rows,
        endpoint_summary_rows=endpoint_summary_rows,
        biospecimen_crosswalk_rows=biospecimen_crosswalk_rows,
        biospecimen_summary_rows=biospecimen_summary_rows,
        source_metadata_rows=source_metadata_rows,
        clinical_table_manifest_rows=clinical_table_manifest_rows,
        biospecimen_table_manifest_rows=biospecimen_table_manifest_rows,
        clinical_patient_rows=clinical_patient_rows,
        clinical_follow_up_rows=clinical_follow_up_rows,
        biospecimen_sample_rows=biospecimen_sample_rows,
        input_paths=input_paths,
    )


def required_field_names(
    rows: list[dict[str, str]],
    *,
    table_name: str,
    proposed_role: str | None = None,
) -> list[str]:
    field_names: list[str] = []
    seen: set[str] = set()
    for row in rows:
        if str(row.get("table_name") or "") != table_name:
            continue
        if proposed_role is not None and str(row.get("proposed_role") or "") != proposed_role:
            continue
        field_name = str(row.get("field_name") or "")
        if not field_name or field_name in seen:
            continue
        seen.add(field_name)
        field_names.append(field_name)
    return field_names


def assert_required_columns(
    observed_columns: set[str],
    required_columns: list[str],
    label: str,
) -> None:
    missing_columns = [column for column in required_columns if column not in observed_columns]
    if missing_columns:
        raise MinimalDryRunCohortError(
            f"{label} is missing required columns: {missing_columns}"
        )


def assert_unique_identifier(rows: list[dict[str, str]], field_name: str, label: str) -> None:
    seen: set[str] = set()
    for row in rows:
        value = str(row.get(field_name) or "")
        if not value:
            raise MinimalDryRunCohortError(f"{label} contains an empty {field_name} value.")
        if value in seen:
            raise MinimalDryRunCohortError(f"{label} contains duplicate {field_name} value: {value}")
        seen.add(value)


def group_records(
    rows: list[dict[str, str]],
    key_fields: list[str],
    list_fields: list[str],
) -> dict[tuple[str, ...], dict[str, Any]]:
    grouped: dict[tuple[str, ...], dict[str, Any]] = {}
    for row in rows:
        key = tuple(str(row.get(field) or "") for field in key_fields)
        if key not in grouped:
            grouped[key] = {
                "row_count": 0,
                "lists": {field: [] for field in list_fields},
            }
        grouped[key]["row_count"] += 1
        for field in list_fields:
            grouped[key]["lists"][field].append(str(row.get(field) or ""))
    return grouped


def build_source_case_mismatch_flag(
    clinical_patient_count: int,
    clinical_source_case_coverage: int,
    biospecimen_patient_uuid_distinct_count: int,
) -> str:
    return (
        "case_count_mismatch_carried_forward_"
        f"{clinical_patient_count}_{clinical_source_case_coverage}_{biospecimen_patient_uuid_distinct_count}"
    )


def build_join_status_flags(
    followup_match_row_count: int,
    biospecimen_match_row_count: int,
) -> list[str]:
    flags: list[str] = []
    if followup_match_row_count > 0:
        flags.append("followup_matched")
    else:
        flags.append("followup_unmatched")
    if followup_match_row_count > 1:
        flags.append("followup_multirow")
    if biospecimen_match_row_count > 0:
        flags.append("biospecimen_sample_matched")
        flags.append("biospecimen_uuid_only_cross_layer_join")
    else:
        flags.append("biospecimen_sample_unmatched")
    if biospecimen_match_row_count > 1:
        flags.append("biospecimen_sample_multirow")
    return flags


def build_cohort_rows(
    dry_run_build_id: str,
    workflow_inputs: WorkflowInputs,
) -> tuple[
    list[dict[str, Any]],
    list[str],
    list[str],
    list[str],
    dict[str, Any],
]:
    clinical_patient_columns = set(workflow_inputs.clinical_patient_rows[0].keys())
    clinical_follow_up_columns = set(workflow_inputs.clinical_follow_up_rows[0].keys())
    biospecimen_sample_columns = set(workflow_inputs.biospecimen_sample_rows[0].keys())

    required_clinical_field_names = required_field_names(
        workflow_inputs.blueprint_required_rows,
        table_name="clinical_patient",
    )
    patient_endpoint_field_names = required_field_names(
        workflow_inputs.blueprint_optional_rows,
        table_name="clinical_patient",
        proposed_role="endpoint_candidate_join_preparation",
    )
    followup_candidate_field_names = required_field_names(
        workflow_inputs.blueprint_optional_rows,
        table_name="clinical_follow_up_v4_0",
        proposed_role="endpoint_candidate_join_preparation",
    )

    assert_required_columns(
        clinical_patient_columns,
        required_clinical_field_names,
        "clinical_patient",
    )
    assert_required_columns(
        clinical_patient_columns,
        patient_endpoint_field_names,
        "clinical_patient patient-layer endpoint candidates",
    )
    assert_required_columns(
        clinical_follow_up_columns,
        ["bcr_patient_uuid", "bcr_patient_barcode", "bcr_followup_barcode", "bcr_followup_uuid"],
        "clinical_follow_up_v4_0 follow-up identifiers",
    )
    assert_required_columns(
        clinical_follow_up_columns,
        followup_candidate_field_names,
        "clinical_follow_up_v4_0 follow-up candidates",
    )
    assert_required_columns(
        biospecimen_sample_columns,
        ["bcr_patient_uuid", "bcr_sample_barcode", "bcr_sample_uuid"],
        "biospecimen_sample sample-anchor identifiers",
    )

    assert_unique_identifier(workflow_inputs.clinical_patient_rows, "bcr_patient_uuid", "clinical_patient")
    assert_unique_identifier(workflow_inputs.clinical_patient_rows, "bcr_patient_barcode", "clinical_patient")

    baseline_field_names = [
        field_name
        for field_name in required_clinical_field_names
        if field_name not in {"bcr_patient_barcode", "bcr_patient_uuid"}
    ]

    clinical_patient_count = len(workflow_inputs.clinical_patient_rows)
    clinical_barcodes = [str(row["bcr_patient_barcode"]) for row in workflow_inputs.clinical_patient_rows]
    clinical_uuids = [str(row["bcr_patient_uuid"]) for row in workflow_inputs.clinical_patient_rows]
    clinical_barcode_set = set(clinical_barcodes)
    clinical_uuid_set = set(clinical_uuids)
    clinical_pair_set = {
        (str(row["bcr_patient_uuid"]), str(row["bcr_patient_barcode"]))
        for row in workflow_inputs.clinical_patient_rows
    }

    followup_grouped = group_records(
        workflow_inputs.clinical_follow_up_rows,
        key_fields=["bcr_patient_uuid", "bcr_patient_barcode"],
        list_fields=["bcr_followup_barcode", "bcr_followup_uuid", *followup_candidate_field_names],
    )
    biospecimen_grouped = group_records(
        workflow_inputs.biospecimen_sample_rows,
        key_fields=["bcr_patient_uuid"],
        list_fields=["bcr_patient_uuid", "bcr_sample_barcode", "bcr_sample_uuid"],
    )

    followup_matched_right_row_count = 0
    followup_right_only_rows: list[dict[str, str]] = []
    for row in workflow_inputs.clinical_follow_up_rows:
        key = (str(row["bcr_patient_uuid"]), str(row["bcr_patient_barcode"]))
        if key in clinical_pair_set:
            followup_matched_right_row_count += 1
        else:
            followup_right_only_rows.append(row)

    biospecimen_matched_right_row_count = 0
    biospecimen_right_only_rows: list[dict[str, str]] = []
    for row in workflow_inputs.biospecimen_sample_rows:
        if str(row["bcr_patient_uuid"]) in clinical_uuid_set:
            biospecimen_matched_right_row_count += 1
        else:
            biospecimen_right_only_rows.append(row)

    source_case_submitter_ids = unique_case_submitter_ids(
        workflow_inputs.source_metadata_rows,
        "clinical",
    )
    clinical_source_case_coverage = len(source_case_submitter_ids)
    extra_source_case_ids = [
        case_id for case_id in source_case_submitter_ids if case_id not in clinical_barcode_set
    ]
    clinical_cases_missing_from_source = [
        case_id for case_id in clinical_barcodes if case_id not in set(source_case_submitter_ids)
    ]

    extra_biospecimen_uuid_order = ordered_unique(
        [str(row["bcr_patient_uuid"]) for row in biospecimen_right_only_rows]
    )
    extra_biospecimen_uuid_details: dict[str, list[dict[str, str]]] = {
        patient_uuid: [] for patient_uuid in extra_biospecimen_uuid_order
    }
    for row in biospecimen_right_only_rows:
        extra_biospecimen_uuid_details[str(row["bcr_patient_uuid"])].append(row)

    biospecimen_patient_uuid_distinct_count = parse_int(
        find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_sample",
            field_name="bcr_patient_uuid",
        )["distinct_non_missing_count"],
        "biospecimen_sample bcr_patient_uuid distinct_non_missing_count",
    )
    mismatch_flag = build_source_case_mismatch_flag(
        clinical_patient_count=clinical_patient_count,
        clinical_source_case_coverage=clinical_source_case_coverage,
        biospecimen_patient_uuid_distinct_count=biospecimen_patient_uuid_distinct_count,
    )
    ambiguity_note_flags = [*STATIC_AMBIGUITY_NOTE_FLAGS[:1], mismatch_flag, *STATIC_AMBIGUITY_NOTE_FLAGS[1:]]

    cohort_rows: list[dict[str, Any]] = []
    followup_matched_patient_count = 0
    followup_patients_with_multiple_rows = 0
    biospecimen_matched_patient_count = 0
    biospecimen_patients_with_multiple_rows = 0
    both_id_count = 0
    uuid_only_count = 0
    barcode_only_count = 0
    neither_id_count = 0
    max_followup_rows_per_patient = 0
    max_biospecimen_rows_per_patient = 0

    for row_index, patient_row in enumerate(workflow_inputs.clinical_patient_rows, start=1):
        patient_uuid = str(patient_row["bcr_patient_uuid"])
        patient_barcode = str(patient_row["bcr_patient_barcode"])
        followup_key = (patient_uuid, patient_barcode)
        followup_payload = followup_grouped.get(followup_key)
        followup_match_row_count = int(followup_payload["row_count"]) if followup_payload else 0
        biospecimen_payload = biospecimen_grouped.get((patient_uuid,))
        biospecimen_match_row_count = int(biospecimen_payload["row_count"]) if biospecimen_payload else 0

        if patient_uuid and patient_barcode:
            both_id_count += 1
        elif patient_uuid:
            uuid_only_count += 1
        elif patient_barcode:
            barcode_only_count += 1
        else:
            neither_id_count += 1

        if followup_match_row_count > 0:
            followup_matched_patient_count += 1
        if followup_match_row_count > 1:
            followup_patients_with_multiple_rows += 1
        if biospecimen_match_row_count > 0:
            biospecimen_matched_patient_count += 1
        if biospecimen_match_row_count > 1:
            biospecimen_patients_with_multiple_rows += 1
        if followup_match_row_count > max_followup_rows_per_patient:
            max_followup_rows_per_patient = followup_match_row_count
        if biospecimen_match_row_count > max_biospecimen_rows_per_patient:
            max_biospecimen_rows_per_patient = biospecimen_match_row_count

        cohort_row: dict[str, Any] = {
            "dry_run_build_id": dry_run_build_id,
            "provisional_patient_row_id": row_index,
            "bcr_patient_barcode": patient_barcode,
            "bcr_patient_uuid": patient_uuid,
        }
        for field_name in baseline_field_names:
            cohort_row[field_name] = str(patient_row.get(field_name) or "")
        for field_name in patient_endpoint_field_names:
            cohort_row[f"clinical_patient__{field_name}"] = str(patient_row.get(field_name) or "")

        cohort_row["clinical_follow_up_v4_0__match_row_count"] = followup_match_row_count
        cohort_row["clinical_follow_up_v4_0__bcr_followup_barcode_json"] = (
            json_list(followup_payload["lists"]["bcr_followup_barcode"]) if followup_payload else "[]"
        )
        cohort_row["clinical_follow_up_v4_0__bcr_followup_uuid_json"] = (
            json_list(followup_payload["lists"]["bcr_followup_uuid"]) if followup_payload else "[]"
        )
        for field_name in followup_candidate_field_names:
            cohort_row[f"clinical_follow_up_v4_0__{field_name}_json"] = (
                json_list(followup_payload["lists"][field_name]) if followup_payload else "[]"
            )

        cohort_row["biospecimen_sample__match_row_count"] = biospecimen_match_row_count
        cohort_row["biospecimen_sample__bcr_patient_uuid_json"] = (
            json_list(biospecimen_payload["lists"]["bcr_patient_uuid"]) if biospecimen_payload else "[]"
        )
        cohort_row["biospecimen_sample__bcr_sample_barcode_json"] = (
            json_list(biospecimen_payload["lists"]["bcr_sample_barcode"]) if biospecimen_payload else "[]"
        )
        cohort_row["biospecimen_sample__bcr_sample_uuid_json"] = (
            json_list(biospecimen_payload["lists"]["bcr_sample_uuid"]) if biospecimen_payload else "[]"
        )

        cohort_row["ambiguity_note_flags_json"] = json_list(ambiguity_note_flags)
        cohort_row["join_status_flags_json"] = json_list(
            build_join_status_flags(
                followup_match_row_count=followup_match_row_count,
                biospecimen_match_row_count=biospecimen_match_row_count,
            )
        )
        cohort_rows.append(cohort_row)

    metrics = {
        "clinical_patient_count": clinical_patient_count,
        "clinical_source_case_coverage": clinical_source_case_coverage,
        "biospecimen_patient_uuid_distinct_count": biospecimen_patient_uuid_distinct_count,
        "followup_total_row_count": len(workflow_inputs.clinical_follow_up_rows),
        "followup_matched_patient_count": followup_matched_patient_count,
        "followup_unmatched_patient_count": clinical_patient_count - followup_matched_patient_count,
        "followup_matched_right_row_count": followup_matched_right_row_count,
        "followup_right_only_row_count": len(followup_right_only_rows),
        "followup_patients_with_multiple_rows": followup_patients_with_multiple_rows,
        "followup_max_rows_per_patient": max_followup_rows_per_patient,
        "biospecimen_total_row_count": len(workflow_inputs.biospecimen_sample_rows),
        "biospecimen_matched_patient_count": biospecimen_matched_patient_count,
        "biospecimen_unmatched_patient_count": clinical_patient_count - biospecimen_matched_patient_count,
        "biospecimen_matched_right_row_count": biospecimen_matched_right_row_count,
        "biospecimen_right_only_row_count": len(biospecimen_right_only_rows),
        "biospecimen_right_only_patient_uuid_count": len(extra_biospecimen_uuid_order),
        "biospecimen_patients_with_multiple_rows": biospecimen_patients_with_multiple_rows,
        "biospecimen_max_rows_per_patient": max_biospecimen_rows_per_patient,
        "rows_with_both_ids": both_id_count,
        "rows_with_only_uuid": uuid_only_count,
        "rows_with_only_barcode": barcode_only_count,
        "rows_with_neither_id": neither_id_count,
        "extra_source_case_id_count": len(extra_source_case_ids),
        "clinical_cases_missing_from_source_count": len(clinical_cases_missing_from_source),
        "carried_ambiguity_note_flag_count": len(ambiguity_note_flags),
        "final_cohort_row_count": len(cohort_rows),
    }
    extras = {
        "baseline_field_names": baseline_field_names,
        "patient_endpoint_field_names": patient_endpoint_field_names,
        "followup_candidate_field_names": followup_candidate_field_names,
        "ambiguity_note_flags": ambiguity_note_flags,
        "extra_source_case_ids": extra_source_case_ids,
        "clinical_cases_missing_from_source": clinical_cases_missing_from_source,
        "extra_biospecimen_uuid_order": extra_biospecimen_uuid_order,
        "extra_biospecimen_uuid_details": extra_biospecimen_uuid_details,
        "followup_right_only_rows": followup_right_only_rows,
    }
    return (
        cohort_rows,
        baseline_field_names,
        patient_endpoint_field_names,
        followup_candidate_field_names,
        {"metrics": metrics, "extras": extras},
    )


def build_join_audit_rows(
    dry_run_build_id: str,
    workflow_inputs: WorkflowInputs,
    derived: dict[str, Any],
) -> list[dict[str, Any]]:
    metrics = derived["metrics"]
    extras = derived["extras"]
    mismatch_flag = next(
        flag
        for flag in extras["ambiguity_note_flags"]
        if flag.startswith("case_count_mismatch_carried_forward_")
    )

    return [
        {
            "dry_run_build_id": dry_run_build_id,
            "join_step_order": 1,
            "join_step_name": "clinical_patient baseline build",
            "left_table_name": "patient/case",
            "right_table_name": "clinical_patient",
            "starting_left_row_count": metrics["clinical_patient_count"],
            "matched_left_row_count": metrics["clinical_patient_count"],
            "unmatched_left_row_count": 0,
            "matched_right_row_count": metrics["clinical_patient_count"],
            "right_only_row_count": 0,
            "left_rows_with_multiple_matches": 0,
            "max_right_matches_per_left_row": 1,
            "join_rule_used": (
                "Use clinical_patient as the provisional patient universe in source-row order; "
                "preserve bcr_patient_barcode and bcr_patient_uuid side by side."
            ),
            "ambiguity_carried_forward": "[]",
            "notes": (
                "Baseline fields come only from cohort_blueprint_required_fields rows with table_name=clinical_patient."
            ),
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "join_step_order": 2,
            "join_step_name": "attach patient-layer endpoint candidates",
            "left_table_name": "clinical_patient",
            "right_table_name": "clinical_patient",
            "starting_left_row_count": metrics["clinical_patient_count"],
            "matched_left_row_count": metrics["clinical_patient_count"],
            "unmatched_left_row_count": 0,
            "matched_right_row_count": metrics["clinical_patient_count"],
            "right_only_row_count": 0,
            "left_rows_with_multiple_matches": 0,
            "max_right_matches_per_left_row": 1,
            "join_rule_used": (
                "Carry patient-table endpoint candidates in place with source-specific prefixes and no endpoint collapse."
            ),
            "ambiguity_carried_forward": json_list(
                ["endpoint_overlap_fields_kept_side_by_side"]
            ),
            "notes": (
                "Attached patient-layer candidate fields: "
                + ", ".join(extras["patient_endpoint_field_names"])
                + "."
            ),
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "join_step_order": 3,
            "join_step_name": "attach follow-up-layer endpoint candidates",
            "left_table_name": "clinical_patient",
            "right_table_name": "clinical_follow_up_v4_0",
            "starting_left_row_count": metrics["clinical_patient_count"],
            "matched_left_row_count": metrics["followup_matched_patient_count"],
            "unmatched_left_row_count": metrics["followup_unmatched_patient_count"],
            "matched_right_row_count": metrics["followup_matched_right_row_count"],
            "right_only_row_count": metrics["followup_right_only_row_count"],
            "left_rows_with_multiple_matches": metrics["followup_patients_with_multiple_rows"],
            "max_right_matches_per_left_row": metrics["followup_max_rows_per_patient"],
            "join_rule_used": (
                "Group clinical_follow_up_v4_0 by exact (bcr_patient_uuid, bcr_patient_barcode) and retain JSON arrays per patient row."
            ),
            "ambiguity_carried_forward": json_list(
                ["endpoint_overlap_fields_kept_side_by_side"]
            ),
            "notes": (
                "Attached follow-up candidate fields as provenance-bearing JSON arrays: "
                + ", ".join(extras["followup_candidate_field_names"])
                + "."
            ),
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "join_step_order": 4,
            "join_step_name": "attach biospecimen sample anchor evidence",
            "left_table_name": "clinical_patient",
            "right_table_name": "biospecimen_sample",
            "starting_left_row_count": metrics["clinical_patient_count"],
            "matched_left_row_count": metrics["biospecimen_matched_patient_count"],
            "unmatched_left_row_count": metrics["biospecimen_unmatched_patient_count"],
            "matched_right_row_count": metrics["biospecimen_matched_right_row_count"],
            "right_only_row_count": metrics["biospecimen_right_only_row_count"],
            "left_rows_with_multiple_matches": metrics["biospecimen_patients_with_multiple_rows"],
            "max_right_matches_per_left_row": metrics["biospecimen_max_rows_per_patient"],
            "join_rule_used": (
                "Group biospecimen_sample by bcr_patient_uuid only and retain sample-anchor evidence as JSON arrays."
            ),
            "ambiguity_carried_forward": json_list(
                [
                    "parallel_patient_identifiers_retained",
                    mismatch_flag,
                    "biospecimen_child_layer_expansion_deferred",
                ]
            ),
            "notes": (
                "Biospecimen sample evidence remains UUID-anchored only; right-only patient UUIDs not attached to clinical_patient="
                + json_list(extras["extra_biospecimen_uuid_order"])
                + "."
            ),
        },
    ]


def build_row_counts_rows(
    dry_run_build_id: str,
    workflow_inputs: WorkflowInputs,
    derived: dict[str, Any],
) -> list[dict[str, Any]]:
    metrics = derived["metrics"]
    extras = derived["extras"]
    rows = [
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "patient_universe",
            "count_metric": "clinical_patient_row_count",
            "count_value": metrics["clinical_patient_count"],
            "notes": "Provisional patient universe from clinical_patient source rows.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "patient_universe",
            "count_metric": "final_cohort_row_count",
            "count_value": metrics["final_cohort_row_count"],
            "notes": "Final dry-run cohort row count; must equal clinical_patient row count.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "patient_universe",
            "count_metric": "rows_with_both_ids",
            "count_value": metrics["rows_with_both_ids"],
            "notes": "Rows with both bcr_patient_barcode and bcr_patient_uuid present.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "patient_universe",
            "count_metric": "rows_with_only_uuid",
            "count_value": metrics["rows_with_only_uuid"],
            "notes": "Rows with only bcr_patient_uuid present.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "patient_universe",
            "count_metric": "rows_with_only_barcode",
            "count_value": metrics["rows_with_only_barcode"],
            "notes": "Rows with only bcr_patient_barcode present.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "patient_universe",
            "count_metric": "rows_with_neither_id",
            "count_value": metrics["rows_with_neither_id"],
            "notes": "Rows with neither patient identifier present.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "followup_join",
            "count_metric": "clinical_follow_up_total_row_count",
            "count_value": metrics["followup_total_row_count"],
            "notes": "Raw clinical_follow_up_v4_0 source-row count.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "followup_join",
            "count_metric": "clinical_patient_rows_with_followup_match",
            "count_value": metrics["followup_matched_patient_count"],
            "notes": "Clinical patient rows with at least one exact (uuid, barcode) follow-up match.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "followup_join",
            "count_metric": "clinical_patient_rows_without_followup_match",
            "count_value": metrics["followup_unmatched_patient_count"],
            "notes": "Clinical patient rows without any exact (uuid, barcode) follow-up match.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "followup_join",
            "count_metric": "matched_followup_row_count",
            "count_value": metrics["followup_matched_right_row_count"],
            "notes": "Follow-up source rows matched into grouped patient-level evidence.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "followup_join",
            "count_metric": "right_only_followup_row_count",
            "count_value": metrics["followup_right_only_row_count"],
            "notes": "Follow-up source rows with no clinical_patient exact-key match.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "followup_join",
            "count_metric": "clinical_patient_rows_with_multiple_followup_rows",
            "count_value": metrics["followup_patients_with_multiple_rows"],
            "notes": "Clinical patient rows with more than one grouped follow-up row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "followup_join",
            "count_metric": "max_followup_rows_per_patient",
            "count_value": metrics["followup_max_rows_per_patient"],
            "notes": "Maximum grouped follow-up rows attached to a single clinical patient row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "biospecimen_sample_total_row_count",
            "count_value": metrics["biospecimen_total_row_count"],
            "notes": "Raw biospecimen_sample source-row count.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "clinical_patient_rows_with_biospecimen_sample_match",
            "count_value": metrics["biospecimen_matched_patient_count"],
            "notes": "Clinical patient rows with at least one biospecimen_sample UUID match.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "clinical_patient_rows_without_biospecimen_sample_match",
            "count_value": metrics["biospecimen_unmatched_patient_count"],
            "notes": "Clinical patient rows without any biospecimen_sample UUID match.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "matched_biospecimen_sample_row_count",
            "count_value": metrics["biospecimen_matched_right_row_count"],
            "notes": "biospecimen_sample rows matched into grouped patient-level evidence.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "right_only_biospecimen_sample_row_count",
            "count_value": metrics["biospecimen_right_only_row_count"],
            "notes": "biospecimen_sample rows with UUIDs absent from clinical_patient.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "right_only_biospecimen_patient_uuid_count",
            "count_value": metrics["biospecimen_right_only_patient_uuid_count"],
            "notes": "Distinct biospecimen patient UUIDs absent from clinical_patient.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "clinical_patient_rows_with_multiple_biospecimen_samples",
            "count_value": metrics["biospecimen_patients_with_multiple_rows"],
            "notes": "Clinical patient rows with more than one grouped biospecimen_sample row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "biospecimen_sample_join",
            "count_metric": "max_biospecimen_sample_rows_per_patient",
            "count_value": metrics["biospecimen_max_rows_per_patient"],
            "notes": "Maximum grouped biospecimen_sample rows attached to a single clinical patient row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "mismatch_signals",
            "count_metric": "clinical_source_case_submitter_coverage",
            "count_value": metrics["clinical_source_case_coverage"],
            "notes": "Unique clinical source case coverage from saved source metadata.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "mismatch_signals",
            "count_metric": "biospecimen_patient_uuid_distinct_count",
            "count_value": metrics["biospecimen_patient_uuid_distinct_count"],
            "notes": "Distinct patient UUID count from biospecimen_sample evidence.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "mismatch_signals",
            "count_metric": "extra_source_case_id_count",
            "count_value": metrics["extra_source_case_id_count"],
            "notes": "Case submitter IDs present in source metadata but absent from clinical_patient.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "mismatch_signals",
            "count_metric": "clinical_cases_missing_from_source_count",
            "count_value": metrics["clinical_cases_missing_from_source_count"],
            "notes": "Clinical patient barcodes absent from source metadata coverage.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "ambiguity",
            "count_metric": "carried_ambiguity_note_flag_count",
            "count_value": metrics["carried_ambiguity_note_flag_count"],
            "notes": "Static ambiguity-note flags carried into every dry-run cohort row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "ambiguity",
            "count_metric": "upstream_blueprint_ambiguity_count",
            "count_value": len(workflow_inputs.blueprint_ambiguity_rows),
            "notes": "Current cohort-blueprint ambiguity count carried forward as context.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "count_stage": "ambiguity",
            "count_metric": "upstream_actions_needed_before_minimal_build_count",
            "count_value": sum(
                truthy_string(row.get("needed_before_minimal_build"))
                for row in workflow_inputs.ambiguity_actions_rows
            ),
            "notes": "Count of grouped ambiguity-resolution actions still flagged before a minimal build.",
        },
    ]

    if extras["extra_source_case_ids"]:
        rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "count_stage": "mismatch_signals",
                "count_metric": "extra_source_case_ids_json",
                "count_value": json_list(extras["extra_source_case_ids"]),
                "notes": "Source-only case submitter IDs retained as an audit note.",
            }
        )

    return rows


def build_exclusions_rows(
    dry_run_build_id: str,
    workflow_inputs: WorkflowInputs,
    derived: dict[str, Any],
) -> list[dict[str, Any]]:
    extras = derived["extras"]
    metrics = derived["metrics"]
    exclusions_rows: list[dict[str, Any]] = []

    endpoint_crosswalk_lookup = {
        (str(row["table_name"]), str(row["field_name"])): row
        for row in workflow_inputs.endpoint_crosswalk_rows
    }
    for deferred_row in workflow_inputs.blueprint_deferred_rows:
        if str(deferred_row.get("source_layer") or "") != "endpoint":
            continue
        table_name = str(deferred_row["table_name"])
        field_name = str(deferred_row["field_name"])
        crosswalk_row = endpoint_crosswalk_lookup.get((table_name, field_name))
        details = (
            f"proposed_role={deferred_row['proposed_role']}; "
            f"rationale={deferred_row['rationale']}"
        )
        if crosswalk_row is not None:
            details = (
                details
                + "; crosswalk_role="
                + str(crosswalk_row.get("crosswalk_role") or "")
                + "; non_missing_count="
                + str(crosswalk_row.get("non_missing_count") or "")
            )
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "deferred_endpoint_field",
                "scope_level": "field",
                "source_table": table_name,
                "entity_id": "field_name",
                "entity_value": field_name,
                "reason": "Deferred from the minimal dry-run build and kept out of endpoint freeze logic.",
                "count_value": 1,
                "details": details,
            }
        )

    for action_row in workflow_inputs.ambiguity_actions_rows:
        action_id = str(action_row.get("action_id") or "")
        if action_id not in {"AR04", "AR05", "AR06"}:
            continue
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "ambiguity_resolution_action_applied",
                "scope_level": "workflow_rule",
                "source_table": str(action_row.get("action_scope") or ""),
                "entity_id": "action_id",
                "entity_value": action_id,
                "reason": "Documented dry-run containment action applied.",
                "count_value": 1,
                "details": str(action_row.get("provisional_rule_if_any") or ""),
            }
        )

    ambiguity_row = find_single_row(
        workflow_inputs.blueprint_ambiguity_rows,
        "cohort blueprint ambiguities",
        ambiguity_category="treatment_table_multiplicity",
    )
    for table_name in ["clinical_drug", "clinical_radiation"]:
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "excluded_treatment_detail_table",
                "scope_level": "table",
                "source_table": table_name,
                "entity_id": "table_name",
                "entity_value": table_name,
                "reason": "One-to-many treatment-detail rows are excluded from the minimal patient-level row build.",
                "count_value": table_row_count_from_manifest(workflow_inputs.clinical_table_manifest_rows, table_name),
                "details": str(ambiguity_row.get("current_saved_evidence") or ""),
            }
        )

    for table_name in [
        "biospecimen_portion",
        "biospecimen_analyte",
        "biospecimen_slide",
        "biospecimen_aliquot",
        "biospecimen_protocol",
        "biospecimen_shipment_portion",
        "biospecimen_diagnostic_slides",
    ]:
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "excluded_child_biospecimen_table",
                "scope_level": "table",
                "source_table": table_name,
                "entity_id": "table_name",
                "entity_value": table_name,
                "reason": "The minimal dry-run cohort stays at patient/case level and keeps biospecimen only through the sample anchor.",
                "count_value": table_row_count_from_manifest(workflow_inputs.biospecimen_table_manifest_rows, table_name),
                "details": "Child-layer or side-layer biospecimen table excluded from the row set.",
            }
        )

    for case_id in extras["extra_source_case_ids"]:
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "unmatched_source_case_submitter_id",
                "scope_level": "case",
                "source_table": "source_metadata",
                "entity_id": "case_submitter_id",
                "entity_value": case_id,
                "reason": "Source metadata includes this case submitter ID but clinical_patient remains the provisional patient universe.",
                "count_value": 1,
                "details": "Retained as an explicit audit mismatch note rather than appended to the dry-run cohort.",
            }
        )

    if metrics["followup_unmatched_patient_count"] > 0:
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "left_unmatched_join_group",
                "scope_level": "patient_rows",
                "source_table": "clinical_follow_up_v4_0",
                "entity_id": "clinical_patient_without_followup_match",
                "entity_value": "exact_uuid_barcode_join",
                "reason": "These clinical_patient rows had no matching grouped follow-up evidence.",
                "count_value": metrics["followup_unmatched_patient_count"],
                "details": "Count of clinical_patient rows with zero clinical_follow_up_v4_0 rows under the exact (bcr_patient_uuid, bcr_patient_barcode) rule.",
            }
        )

    if metrics["followup_right_only_row_count"] > 0:
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "right_only_followup_rows",
                "scope_level": "source_rows",
                "source_table": "clinical_follow_up_v4_0",
                "entity_id": "unmatched_followup_rows",
                "entity_value": "exact_uuid_barcode_join",
                "reason": "Follow-up rows existed but could not be matched to clinical_patient under the exact join rule.",
                "count_value": metrics["followup_right_only_row_count"],
                "details": "Right-only follow-up rows are not appended to the dry-run cohort.",
            }
        )

    for patient_uuid in extras["extra_biospecimen_uuid_order"]:
        rows_for_uuid = extras["extra_biospecimen_uuid_details"][patient_uuid]
        sample_barcodes = [str(row["bcr_sample_barcode"]) for row in rows_for_uuid]
        exclusions_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "exclusion_category": "unmatched_biospecimen_patient_uuid",
                "scope_level": "patient_uuid",
                "source_table": "biospecimen_sample",
                "entity_id": "bcr_patient_uuid",
                "entity_value": patient_uuid,
                "reason": "biospecimen_sample contains this patient UUID but it is absent from clinical_patient.",
                "count_value": len(rows_for_uuid),
                "details": "Unmatched sample barcodes=" + json_list(sample_barcodes),
            }
        )

    return exclusions_rows


def build_summary_rows(
    dry_run_build_id: str,
    workflow_inputs: WorkflowInputs,
    derived: dict[str, Any],
) -> list[dict[str, Any]]:
    metrics = derived["metrics"]
    extras = derived["extras"]
    ambiguity_summary_lookup = build_summary_lookup(workflow_inputs.ambiguity_summary_rows)
    readiness_signal = str(
        ambiguity_summary_lookup["minimal_build_readiness_signal"]["summary_value"]
    )
    mismatch_note_value = (
        f"{metrics['clinical_patient_count']} / "
        f"{metrics['clinical_source_case_coverage']} / "
        f"{metrics['biospecimen_patient_uuid_distinct_count']}"
    )
    summary_rows = [
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "design",
            "summary_metric": "proposed_unit_of_analysis",
            "summary_value": "patient/case",
            "notes": "The minimal dry-run cohort preserves patient/case as the row unit.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "design",
            "summary_metric": "provisional_patient_universe",
            "summary_value": "clinical_patient",
            "notes": "clinical_patient remains the provisional patient universe for this dry-run build.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "design",
            "summary_metric": "biospecimen_policy",
            "summary_value": "sample_anchor_only",
            "notes": "biospecimen is represented only through grouped biospecimen_sample anchor evidence.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "design",
            "summary_metric": "endpoint_policy",
            "summary_value": "candidate_columns_side_by_side_only",
            "notes": "Patient and follow-up endpoint-like fields remain side-by-side provenance-bearing candidates only.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "design",
            "summary_metric": "treatment_policy",
            "summary_value": "excluded_from_row_build",
            "notes": "One-to-many treatment-detail rows remain excluded from the patient-level dry-run row build.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "inputs",
            "summary_metric": "blueprint_run_id",
            "summary_value": workflow_inputs.blueprint_latest_pointer["blueprint_run_id"],
            "notes": "Current cohort blueprint run used as the dry-run field and join source.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "inputs",
            "summary_metric": "ambiguity_resolution_run_id",
            "summary_value": workflow_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"],
            "notes": "Current ambiguity-resolution run used as the dry-run containment source.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "field_selection",
            "summary_metric": "baseline_field_count",
            "summary_value": len(extras["baseline_field_names"]),
            "notes": "Baseline fields selected from cohort_blueprint_required_fields where table_name=clinical_patient, excluding duplicated patient IDs.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "field_selection",
            "summary_metric": "patient_endpoint_candidate_field_count",
            "summary_value": len(extras["patient_endpoint_field_names"]),
            "notes": "Patient-table endpoint candidate fields carried as explicit prefixed columns.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "field_selection",
            "summary_metric": "followup_endpoint_candidate_field_count",
            "summary_value": len(extras["followup_candidate_field_names"]),
            "notes": "Follow-up endpoint candidate fields carried as grouped JSON arrays.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "row_counts",
            "summary_metric": "proposed_patient_universe_row_count",
            "summary_value": metrics["clinical_patient_count"],
            "notes": "clinical_patient source-row count used as the provisional patient universe.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "row_counts",
            "summary_metric": "final_dry_run_row_count",
            "summary_value": metrics["final_cohort_row_count"],
            "notes": "Final dry-run cohort row count.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_followup_candidate_count",
            "summary_value": metrics["followup_matched_patient_count"],
            "notes": "Clinical patient rows with at least one grouped follow-up evidence row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_followup_row_count",
            "summary_value": metrics["followup_matched_right_row_count"],
            "notes": "Follow-up source rows grouped into matched patient-level evidence.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_biospecimen_sample_anchor_count",
            "summary_value": metrics["biospecimen_matched_patient_count"],
            "notes": "Clinical patient rows with at least one grouped biospecimen_sample evidence row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "join_results",
            "summary_metric": "matched_biospecimen_sample_row_count",
            "summary_value": metrics["biospecimen_matched_right_row_count"],
            "notes": "biospecimen_sample source rows grouped into matched patient-level evidence.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "identifier_signals",
            "summary_metric": "rows_with_both_patient_barcode_and_uuid",
            "summary_value": metrics["rows_with_both_ids"],
            "notes": "Clinical patient rows where both patient identifier forms are present.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "identifier_signals",
            "summary_metric": "rows_with_only_one_patient_identifier",
            "summary_value": metrics["rows_with_only_uuid"] + metrics["rows_with_only_barcode"],
            "notes": "Clinical patient rows with only one patient identifier form present.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "ambiguity",
            "summary_metric": "carried_global_ambiguity_note_count",
            "summary_value": metrics["carried_ambiguity_note_flag_count"],
            "notes": "Global ambiguity-note flags written into every dry-run cohort row.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "ambiguity",
            "summary_metric": "upstream_blueprint_ambiguity_count",
            "summary_value": len(workflow_inputs.blueprint_ambiguity_rows),
            "notes": "Current blueprint ambiguity count retained as context for this dry run.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "ambiguity",
            "summary_metric": "upstream_ambiguity_readiness_signal",
            "summary_value": readiness_signal,
            "notes": "The upstream ambiguity-resolution layer remains more conservative than this contained dry-run execution.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "clinical_patient_vs_source_vs_biospecimen_patient_counts",
            "summary_value": mismatch_note_value,
            "notes": "This dry run preserves the current count mismatch as an explicit audit note rather than resolving it.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "extra_source_case_id_count",
            "summary_value": metrics["extra_source_case_id_count"],
            "notes": "Case submitter IDs present in source metadata but not in clinical_patient.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "extra_biospecimen_patient_uuid_count",
            "summary_value": metrics["biospecimen_right_only_patient_uuid_count"],
            "notes": "biospecimen_sample patient UUIDs absent from clinical_patient.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "mismatch_signals",
            "summary_metric": "extra_biospecimen_sample_row_count",
            "summary_value": metrics["biospecimen_right_only_row_count"],
            "notes": "biospecimen_sample rows absent from the clinical_patient provisional patient universe.",
        },
        {
            "dry_run_build_id": dry_run_build_id,
            "summary_section": "readiness",
            "summary_metric": "dry_run_readiness_interpretation",
            "summary_value": "sufficient_for_join_test_not_final_freeze",
            "notes": (
                "The saved evidence is sufficient to execute and audit a minimal dry-run join test, "
                "but unresolved mismatch and endpoint ambiguities still block any final cohort freeze."
            ),
        },
    ]

    if extras["extra_source_case_ids"]:
        summary_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "summary_section": "mismatch_signals",
                "summary_metric": "extra_source_case_ids_json",
                "summary_value": json_list(extras["extra_source_case_ids"]),
                "notes": "Source-only case submitter IDs carried forward as explicit audit notes.",
            }
        )

    if extras["extra_biospecimen_uuid_order"]:
        summary_rows.append(
            {
                "dry_run_build_id": dry_run_build_id,
                "summary_section": "mismatch_signals",
                "summary_metric": "extra_biospecimen_patient_uuids_json",
                "summary_value": json_list(extras["extra_biospecimen_uuid_order"]),
                "notes": "biospecimen-only patient UUIDs carried forward as explicit audit notes.",
            }
        )

    return summary_rows


def build_latest_pointer_payload(
    dry_run_build_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "dry_run_build_id": dry_run_build_id,
        "blueprint_run_id": str(workflow_inputs.blueprint_latest_pointer["blueprint_run_id"]),
        "ambiguity_resolution_run_id": str(
            workflow_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
        ),
        "shortlist_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(workflow_inputs.clinical_biotab_latest_pointer["parse_run_id"]),
        "endpoint_crosswalk_run_id": str(
            workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(
            workflow_inputs.biospecimen_biotab_latest_pointer["parse_run_id"]
        ),
        "clinical_source_run_id": str(workflow_inputs.clinical_biotab_latest_pointer["source_run_id"]),
        "biospecimen_source_run_id": str(
            workflow_inputs.biospecimen_biotab_latest_pointer["source_run_id"]
        ),
        "dry_run_run_directory": repo_relative(output_paths["dry_run_run_directory"], paths.repo_root),
        "minimal_dry_run_cohort_tsv": repo_relative(output_paths["cohort_tsv"], paths.repo_root),
        "minimal_dry_run_join_audit_tsv": repo_relative(output_paths["join_audit_tsv"], paths.repo_root),
        "minimal_dry_run_row_counts_tsv": repo_relative(output_paths["row_counts_tsv"], paths.repo_root),
        "minimal_dry_run_exclusions_tsv": repo_relative(output_paths["exclusions_tsv"], paths.repo_root),
        "minimal_dry_run_summary_tsv": repo_relative(output_paths["summary_tsv"], paths.repo_root),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "cohort_blueprint_latest_json": repo_relative(paths.blueprint_latest_pointer, paths.repo_root),
        "blueprint_ambiguity_resolution_latest_json": repo_relative(
            paths.ambiguity_resolution_latest_pointer,
            paths.repo_root,
        ),
        "clinical_shortlist_latest_json": repo_relative(
            paths.clinical_shortlist_latest_pointer,
            paths.repo_root,
        ),
        "endpoint_crosswalk_latest_json": repo_relative(
            paths.endpoint_crosswalk_latest_pointer,
            paths.repo_root,
        ),
        "biospecimen_crosswalk_latest_json": repo_relative(
            paths.biospecimen_crosswalk_latest_pointer,
            paths.repo_root,
        ),
        "clinical_biotab_latest_json": repo_relative(
            paths.clinical_biotab_latest_pointer,
            paths.repo_root,
        ),
        "biospecimen_biotab_latest_json": repo_relative(
            paths.biospecimen_biotab_latest_pointer,
            paths.repo_root,
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    dry_run_build_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    dry_run_dir = paths.dry_run_runs_root / dry_run_build_id
    run_log_path = dry_run_dir / "run_log.json"

    try:
        trial_config = load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths)
        create_run_directory(dry_run_dir)

        cohort_path = dry_run_dir / "minimal_dry_run_cohort.tsv"
        join_audit_path = dry_run_dir / "minimal_dry_run_join_audit.tsv"
        row_counts_path = dry_run_dir / "minimal_dry_run_row_counts.tsv"
        exclusions_path = dry_run_dir / "minimal_dry_run_exclusions.tsv"
        summary_path = dry_run_dir / "minimal_dry_run_summary.tsv"

        (
            cohort_rows,
            baseline_field_names,
            patient_endpoint_field_names,
            followup_candidate_field_names,
            derived,
        ) = build_cohort_rows(dry_run_build_id=dry_run_build_id, workflow_inputs=workflow_inputs)
        if not cohort_rows:
            raise MinimalDryRunCohortError("Dry-run cohort rows were not generated.")

        join_audit_rows = build_join_audit_rows(
            dry_run_build_id=dry_run_build_id,
            workflow_inputs=workflow_inputs,
            derived=derived,
        )
        row_counts_rows = build_row_counts_rows(
            dry_run_build_id=dry_run_build_id,
            workflow_inputs=workflow_inputs,
            derived=derived,
        )
        exclusions_rows = build_exclusions_rows(
            dry_run_build_id=dry_run_build_id,
            workflow_inputs=workflow_inputs,
            derived=derived,
        )
        summary_rows = build_summary_rows(
            dry_run_build_id=dry_run_build_id,
            workflow_inputs=workflow_inputs,
            derived=derived,
        )

        if not join_audit_rows:
            raise MinimalDryRunCohortError("Join-audit rows were not generated.")
        if not row_counts_rows:
            raise MinimalDryRunCohortError("Row-count rows were not generated.")
        if not exclusions_rows:
            raise MinimalDryRunCohortError("Exclusions rows were not generated.")
        if not summary_rows:
            raise MinimalDryRunCohortError("Summary rows were not generated.")

        cohort_fieldnames = list(cohort_rows[0].keys())
        write_dict_rows_tsv(cohort_path, cohort_fieldnames, cohort_rows)
        write_dict_rows_tsv(join_audit_path, JOIN_AUDIT_FIELDNAMES, join_audit_rows)
        write_dict_rows_tsv(row_counts_path, ROW_COUNT_FIELDNAMES, row_counts_rows)
        write_dict_rows_tsv(exclusions_path, EXCLUSIONS_FIELDNAMES, exclusions_rows)
        write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        summary_lookup = build_summary_lookup(summary_rows)
        metrics = derived["metrics"]
        required_output_columns = [
            "dry_run_build_id",
            "provisional_patient_row_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "clinical_follow_up_v4_0__match_row_count",
            "biospecimen_sample__match_row_count",
            "ambiguity_note_flags_json",
            "join_status_flags_json",
        ]
        required_output_columns.extend(baseline_field_names)
        required_output_columns.extend(
            [f"clinical_patient__{field_name}" for field_name in patient_endpoint_field_names]
        )
        required_output_columns.extend(
            [f"clinical_follow_up_v4_0__{field_name}_json" for field_name in followup_candidate_field_names]
        )

        cohort_has_required_columns = all(
            column in cohort_fieldnames for column in required_output_columns
        )
        final_row_count_matches_patient_universe = (
            metrics["final_cohort_row_count"] == metrics["clinical_patient_count"]
        )
        summary_contains_readiness_interpretation = (
            summary_lookup["dry_run_readiness_interpretation"]["summary_value"]
            == "sufficient_for_join_test_not_final_freeze"
        )
        summary_contains_mismatch_note = (
            "clinical_patient_vs_source_vs_biospecimen_patient_counts" in summary_lookup
        )
        cohort_row_ids_are_sequential = [
            parse_int(row["provisional_patient_row_id"], "provisional_patient_row_id")
            for row in cohort_rows
        ] == list(range(1, len(cohort_rows) + 1))
        join_audit_step_count_valid = len(join_audit_rows) == 4

        output_paths = {
            "dry_run_run_directory": dry_run_dir,
            "cohort_tsv": cohort_path,
            "join_audit_tsv": join_audit_path,
            "row_counts_tsv": row_counts_path,
            "exclusions_tsv": exclusions_path,
            "summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            dry_run_build_id=dry_run_build_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "dry_run_build_id": dry_run_build_id,
            "blueprint_run_id": str(workflow_inputs.blueprint_latest_pointer["blueprint_run_id"]),
            "ambiguity_resolution_run_id": str(
                workflow_inputs.ambiguity_resolution_latest_pointer["ambiguity_resolution_run_id"]
            ),
            "shortlist_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["core_audit_run_id"]),
            "clinical_parse_run_id": str(workflow_inputs.clinical_biotab_latest_pointer["parse_run_id"]),
            "endpoint_crosswalk_run_id": str(
                workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
            ),
            "biospecimen_crosswalk_run_id": str(
                workflow_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
            ),
            "biospecimen_parse_run_id": str(
                workflow_inputs.biospecimen_biotab_latest_pointer["parse_run_id"]
            ),
            "clinical_source_run_id": str(workflow_inputs.clinical_biotab_latest_pointer["source_run_id"]),
            "biospecimen_source_run_id": str(
                workflow_inputs.biospecimen_biotab_latest_pointer["source_run_id"]
            ),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "cohort_blueprint_latest_json": repo_relative(paths.blueprint_latest_pointer, paths.repo_root),
                "blueprint_ambiguity_resolution_latest_json": repo_relative(
                    paths.ambiguity_resolution_latest_pointer,
                    paths.repo_root,
                ),
                "clinical_shortlist_latest_json": repo_relative(
                    paths.clinical_shortlist_latest_pointer,
                    paths.repo_root,
                ),
                "endpoint_crosswalk_latest_json": repo_relative(
                    paths.endpoint_crosswalk_latest_pointer,
                    paths.repo_root,
                ),
                "biospecimen_crosswalk_latest_json": repo_relative(
                    paths.biospecimen_crosswalk_latest_pointer,
                    paths.repo_root,
                ),
                "clinical_biotab_latest_json": repo_relative(
                    paths.clinical_biotab_latest_pointer,
                    paths.repo_root,
                ),
                "biospecimen_biotab_latest_json": repo_relative(
                    paths.biospecimen_biotab_latest_pointer,
                    paths.repo_root,
                ),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "dry_run_run_directory": repo_relative(dry_run_dir, paths.repo_root),
                "minimal_dry_run_cohort_tsv": repo_relative(cohort_path, paths.repo_root),
                "minimal_dry_run_join_audit_tsv": repo_relative(join_audit_path, paths.repo_root),
                "minimal_dry_run_row_counts_tsv": repo_relative(row_counts_path, paths.repo_root),
                "minimal_dry_run_exclusions_tsv": repo_relative(exclusions_path, paths.repo_root),
                "minimal_dry_run_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    metrics["clinical_patient_count"] > 0
                    and metrics["final_cohort_row_count"] > 0
                    and final_row_count_matches_patient_universe
                    and cohort_has_required_columns
                    and cohort_row_ids_are_sequential
                    and join_audit_step_count_valid
                    and summary_contains_readiness_interpretation
                    and summary_contains_mismatch_note
                ),
                "blueprint_latest_pointer_found": True,
                "ambiguity_resolution_latest_pointer_found": True,
                "clinical_shortlist_latest_pointer_found": True,
                "endpoint_crosswalk_latest_pointer_found": True,
                "biospecimen_crosswalk_latest_pointer_found": True,
                "clinical_biotab_latest_pointer_found": True,
                "biospecimen_biotab_latest_pointer_found": True,
                "blueprint_run_log_completed": True,
                "ambiguity_resolution_run_log_completed": True,
                "clinical_shortlist_run_log_completed": True,
                "endpoint_crosswalk_run_log_completed": True,
                "biospecimen_crosswalk_run_log_completed": True,
                "clinical_biotab_run_log_completed": True,
                "biospecimen_biotab_run_log_completed": True,
                "source_run_log_completed": True,
                "required_source_tables_found": True,
                "clinical_patient_row_count_positive": metrics["clinical_patient_count"] > 0,
                "final_cohort_row_count_positive": metrics["final_cohort_row_count"] > 0,
                "clinical_patient_identifiers_unique": True,
                "cohort_has_required_columns": cohort_has_required_columns,
                "cohort_row_ids_are_sequential": cohort_row_ids_are_sequential,
                "join_audit_row_count_positive": len(join_audit_rows) > 0,
                "row_counts_row_count_positive": len(row_counts_rows) > 0,
                "exclusions_row_count_positive": len(exclusions_rows) > 0,
                "summary_row_count_positive": len(summary_rows) > 0,
                "final_cohort_row_count_matches_clinical_patient": final_row_count_matches_patient_universe,
                "summary_contains_readiness_interpretation": summary_contains_readiness_interpretation,
                "summary_contains_mismatch_note": summary_contains_mismatch_note,
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
            },
            "counts": {
                **metrics,
                "baseline_field_count": len(baseline_field_names),
                "patient_endpoint_candidate_field_count": len(patient_endpoint_field_names),
                "followup_endpoint_candidate_field_count": len(followup_candidate_field_names),
                "join_audit_row_count": len(join_audit_rows),
                "row_counts_row_count": len(row_counts_rows),
                "exclusions_row_count": len(exclusions_rows),
                "summary_row_count": len(summary_rows),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "cohort_blueprint_latest_pointer": workflow_inputs.blueprint_latest_pointer,
                "cohort_blueprint_run_log": workflow_inputs.blueprint_run_log,
                "ambiguity_resolution_latest_pointer": workflow_inputs.ambiguity_resolution_latest_pointer,
                "ambiguity_resolution_run_log": workflow_inputs.ambiguity_resolution_run_log,
                "clinical_shortlist_latest_pointer": workflow_inputs.clinical_shortlist_latest_pointer,
                "clinical_shortlist_run_log": workflow_inputs.clinical_shortlist_run_log,
                "endpoint_crosswalk_latest_pointer": workflow_inputs.endpoint_crosswalk_latest_pointer,
                "endpoint_crosswalk_run_log": workflow_inputs.endpoint_crosswalk_run_log,
                "biospecimen_crosswalk_latest_pointer": workflow_inputs.biospecimen_crosswalk_latest_pointer,
                "biospecimen_crosswalk_run_log": workflow_inputs.biospecimen_crosswalk_run_log,
                "clinical_biotab_latest_pointer": workflow_inputs.clinical_biotab_latest_pointer,
                "clinical_biotab_run_log": workflow_inputs.clinical_biotab_run_log,
                "biospecimen_biotab_latest_pointer": workflow_inputs.biospecimen_biotab_latest_pointer,
                "biospecimen_biotab_run_log": workflow_inputs.biospecimen_biotab_run_log,
                "source_run_log": workflow_inputs.source_run_log,
            },
        }

        write_json(run_log_path, run_log_payload)
        write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "dry_run_build_id": dry_run_build_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_minimal_dry_run_cohort",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA minimal dry-run cohort workflow complete.")
    print(f"Dry-run build ID: {run_log['dry_run_build_id']}")
    print(f"Blueprint run ID: {run_log['blueprint_run_id']}")
    print(f"Ambiguity-resolution run ID: {run_log['ambiguity_resolution_run_id']}")
    print(f"Dry-run output directory: {run_log['outputs']['dry_run_run_directory']}")
    print(f"Cohort TSV: {run_log['outputs']['minimal_dry_run_cohort_tsv']}")
    print(f"Join audit TSV: {run_log['outputs']['minimal_dry_run_join_audit_tsv']}")
    print(f"Row counts TSV: {run_log['outputs']['minimal_dry_run_row_counts_tsv']}")
    print(f"Exclusions TSV: {run_log['outputs']['minimal_dry_run_exclusions_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['minimal_dry_run_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
