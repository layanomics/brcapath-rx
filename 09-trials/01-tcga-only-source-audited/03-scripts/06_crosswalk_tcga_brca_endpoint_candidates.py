#!/usr/bin/env python
"""Build an auditable endpoint crosswalk from the latest TCGA-BRCA clinical shortlist."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover
    yaml = None

TARGET_TABLES = (
    "clinical_patient",
    "clinical_follow_up_v4_0",
)
TARGET_TABLE_ORDER = {table_name: index for index, table_name in enumerate(TARGET_TABLES)}
PRIMARY_TABLE_RANK = {
    "clinical_follow_up_v4_0": 0,
    "clinical_patient": 1,
}
FAMILY_ORDER = (
    "survival_status_like",
    "last_contact_like",
    "death_time_like",
    "progression_or_tumor_status_like",
    "new_tumor_event_like",
    "followup_loss_like",
    "other_endpoint_like",
)
ROLE_ORDER = (
    "primary_candidate",
    "secondary_candidate",
    "overlapping_candidate",
    "ambiguous_candidate",
)
CANDIDATE_SOURCES = (
    "usable_endpoint_candidate",
    "escalated_manual_review",
)
SELECTION_RULES = (
    "select_usable_endpoint_candidate",
    "escalate_unclear_endpoint_like",
    "escalate_weak_endpoint_like",
)
CROSSWALK_RULES = (
    "role_primary_best_usable_completeness",
    "role_overlapping_same_canonical_signal",
    "role_secondary_distinct_usable_signal",
    "role_ambiguous_no_usable_primary_or_sparse_signal",
)
NOTES_PLACEHOLDER = "[fill in during endpoint crosswalk review]"
INVENTORY_FIELDNAMES = [
    "crosswalk_run_id",
    "shortlist_run_id",
    "core_audit_run_id",
    "parse_run_id",
    "source_run_id",
    "table_name",
    "field_name",
    "source_position",
    "alternate_column_name",
    "probable_field_group",
    "missing_like_fraction",
    "non_missing_count",
    "distinct_non_missing_count",
    "example_values_small_sample",
    "candidate_source",
    "original_shortlist_bucket",
    "original_shortlist_rule",
    "selection_rule",
    "endpoint_signal_family",
    "endpoint_signal_type",
    "endpoint_signal_rule",
    "canonical_signal_key",
    "manual_review_priority",
    "notes_placeholder",
]
CROSSWALK_FIELDNAMES = [
    "crosswalk_run_id",
    "shortlist_run_id",
    "core_audit_run_id",
    "parse_run_id",
    "source_run_id",
    "endpoint_signal_family",
    "endpoint_signal_type",
    "endpoint_signal_rule",
    "canonical_signal_key",
    "table_name",
    "field_name",
    "source_position",
    "alternate_column_name",
    "candidate_source",
    "original_shortlist_bucket",
    "original_shortlist_rule",
    "selection_rule",
    "missing_like_fraction",
    "non_missing_count",
    "distinct_non_missing_count",
    "example_values_small_sample",
    "crosswalk_role",
    "crosswalk_rule",
    "manual_review_priority",
    "notes_placeholder",
]
SUMMARY_FIELDNAMES = [
    "crosswalk_run_id",
    "shortlist_run_id",
    "core_audit_run_id",
    "parse_run_id",
    "source_run_id",
    "endpoint_signal_family",
    "field_count",
    "tables_present_json",
    "primary_candidates_json",
    "secondary_candidates_json",
    "overlapping_candidates_json",
    "ambiguous_candidates_json",
    "notes_placeholder",
]


class EndpointCrosswalkError(RuntimeError):
    """Raised when the endpoint crosswalk workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the endpoint crosswalk workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    shortlist_latest_pointer: Path
    core_audit_latest_pointer: Path
    crosswalk_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class ParsedTableContext:
    """One parsed clinical table plus schema/header context used for field validation."""

    table_name: str
    parsed_path: Path
    schema_path: Path
    row_count: int
    column_count: int
    header: list[str]
    schema_rows_by_position: dict[int, dict[str, Any]]


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from the latest upstream audit layers."""

    shortlist_latest_pointer: dict[str, Any]
    shortlist_run_log: dict[str, Any]
    core_audit_latest_pointer: dict[str, Any]
    core_audit_run_log: dict[str, Any]
    parse_latest_pointer: dict[str, Any]
    parse_run_log: dict[str, Any]
    shortlist_rows: list[dict[str, Any]]
    core_audit_rows: list[dict[str, Any]]
    core_audit_by_key: dict[tuple[str, str], dict[str, Any]]
    parsed_tables: dict[str, ParsedTableContext]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise EndpointCrosswalkError("Unable to locate the repository root from the script path.")


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
                raise EndpointCrosswalkError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise EndpointCrosswalkError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise EndpointCrosswalkError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise EndpointCrosswalkError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise EndpointCrosswalkError(
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
        raise EndpointCrosswalkError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise EndpointCrosswalkError(f"Expected a JSON object in file: {path}")
    return data


def read_tsv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def read_tsv_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as exc:
            raise EndpointCrosswalkError(f"Parsed TSV is empty: {path}") from exc


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise EndpointCrosswalkError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise EndpointCrosswalkError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise EndpointCrosswalkError(f"Run directory already exists: {path}")
    path.mkdir(parents=False, exist_ok=False)
    return path


def build_workflow_paths() -> WorkflowPaths:
    repo_root = detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    trial_config_data = load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        shortlist_latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_shortlist_latest.json",
        core_audit_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_core_field_audit_latest.json"
        ),
        crosswalk_runs_root=audit_root / "tcga-brca" / "variables" / "endpoint_crosswalk_runs",
        latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_endpoint_crosswalk_latest.json",
    )


def normalize_key(row: dict[str, Any]) -> tuple[str, str]:
    return str(row["table_name"]), str(row["field_name"])


def normalize_field_name(field_name: str) -> str:
    return field_name.strip().lower()


def parse_int(value: Any, label: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise EndpointCrosswalkError(f"Expected integer value for {label}: {value!r}") from exc


def parse_float(value: Any, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise EndpointCrosswalkError(f"Expected float value for {label}: {value!r}") from exc


def read_required_tsv(path: Path, label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise EndpointCrosswalkError(f"Required {label} does not exist: {path}")
    rows = read_tsv_dict_rows(path)
    if not rows:
        raise EndpointCrosswalkError(f"Required {label} has no rows: {path}")
    return rows


def validate_unique_field_keys(rows: list[dict[str, Any]], label: str) -> None:
    seen: set[tuple[str, str]] = set()
    duplicates: set[tuple[str, str]] = set()
    for row in rows:
        key = normalize_key(row)
        if key in seen:
            duplicates.add(key)
        seen.add(key)
    if duplicates:
        sample = sorted(f"{table_name}:{field_name}" for table_name, field_name in duplicates)[:10]
        raise EndpointCrosswalkError(f"Duplicate field keys found in {label}: {sample}")


def compare_shortlist_and_core_rows(
    shortlist_rows: list[dict[str, Any]],
    core_audit_by_key: dict[tuple[str, str], dict[str, Any]],
) -> None:
    fields_to_compare = (
        "source_position",
        "alternate_column_name",
        "probable_field_group",
        "missing_like_fraction",
        "non_missing_count",
        "distinct_non_missing_count",
        "example_values_small_sample",
    )

    for shortlist_row in shortlist_rows:
        key = normalize_key(shortlist_row)
        if key not in core_audit_by_key:
            raise EndpointCrosswalkError(f"Shortlist field is missing from core audit TSV: {key}")
        core_row = core_audit_by_key[key]
        for field_name in fields_to_compare:
            shortlist_value = shortlist_row[field_name]
            core_value = core_row[field_name]
            if shortlist_value != core_value:
                raise EndpointCrosswalkError(
                    "Shortlist and core audit disagree for "
                    f"{key} field '{field_name}': {shortlist_value!r} vs {core_value!r}"
                )


def build_core_audit_rows(raw_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    typed_rows: list[dict[str, Any]] = []
    for row in raw_rows:
        typed_rows.append(
            {
                **row,
                "source_position": parse_int(row["source_position"], "source_position"),
                "missing_like_fraction": round(
                    parse_float(row["missing_like_fraction"], "missing_like_fraction"),
                    6,
                ),
                "non_missing_count": parse_int(row["non_missing_count"], "non_missing_count"),
                "distinct_non_missing_count": parse_int(
                    row["distinct_non_missing_count"],
                    "distinct_non_missing_count",
                ),
            }
        )
    return typed_rows


def build_shortlist_rows(raw_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    typed_rows: list[dict[str, Any]] = []
    for row in raw_rows:
        typed_rows.append(
            {
                **row,
                "source_position": parse_int(row["source_position"], "source_position"),
                "missing_like_fraction": round(
                    parse_float(row["missing_like_fraction"], "missing_like_fraction"),
                    6,
                ),
                "non_missing_count": parse_int(row["non_missing_count"], "non_missing_count"),
                "distinct_non_missing_count": parse_int(
                    row["distinct_non_missing_count"],
                    "distinct_non_missing_count",
                ),
            }
        )
    return typed_rows


def build_parsed_table_context(
    *,
    repo_root: Path,
    manifest_row: dict[str, str],
    expected_parse_run_id: str,
    expected_source_run_id: str,
) -> ParsedTableContext:
    parse_run_id = str(manifest_row.get("parse_run_id") or "").strip()
    source_run_id = str(manifest_row.get("source_run_id") or "").strip()
    table_name = str(manifest_row.get("table_name") or "").strip()
    if parse_run_id != expected_parse_run_id:
        raise EndpointCrosswalkError(
            f"Parsed table manifest parse_run_id mismatch for {table_name}: {parse_run_id} vs {expected_parse_run_id}"
        )
    if source_run_id != expected_source_run_id:
        raise EndpointCrosswalkError(
            f"Parsed table manifest source_run_id mismatch for {table_name}: {source_run_id} vs {expected_source_run_id}"
        )

    parsed_path = repo_root / str(manifest_row.get("output_path") or "")
    schema_path = repo_root / str(manifest_row.get("schema_path") or "")
    row_count = parse_int(manifest_row.get("row_count"), f"{table_name} row_count")
    column_count = parse_int(manifest_row.get("column_count"), f"{table_name} column_count")

    for required_path, label in [
        (parsed_path, f"{table_name} parsed TSV"),
        (schema_path, f"{table_name} schema TSV"),
    ]:
        if not required_path.exists():
            raise EndpointCrosswalkError(f"Expected {label} does not exist: {required_path}")

    header = read_tsv_header(parsed_path)
    if len(header) != column_count:
        raise EndpointCrosswalkError(
            f"Parsed TSV header length does not match manifest column_count for {table_name}: "
            f"{len(header)} vs {column_count}"
        )

    schema_rows = read_required_tsv(schema_path, f"{table_name} schema TSV")
    if len(schema_rows) != column_count:
        raise EndpointCrosswalkError(
            f"Schema TSV row count does not match manifest column_count for {table_name}: "
            f"{len(schema_rows)} vs {column_count}"
        )

    schema_rows_by_position: dict[int, dict[str, Any]] = {}
    for schema_row in schema_rows:
        row_table_name = str(schema_row.get("table_name") or "").strip()
        if row_table_name != table_name:
            raise EndpointCrosswalkError(
                f"Schema TSV table_name mismatch for {table_name}: found {row_table_name}"
            )
        position = parse_int(schema_row["column_position"], f"{table_name} schema column_position")
        if position in schema_rows_by_position:
            raise EndpointCrosswalkError(f"Duplicate schema column_position {position} for {table_name}")
        schema_rows_by_position[position] = {
            **schema_row,
            "column_position": position,
        }

    expected_positions = list(range(1, column_count + 1))
    if sorted(schema_rows_by_position) != expected_positions:
        raise EndpointCrosswalkError(f"Schema TSV positions are incomplete for {table_name}")

    return ParsedTableContext(
        table_name=table_name,
        parsed_path=parsed_path,
        schema_path=schema_path,
        row_count=row_count,
        column_count=column_count,
        header=header,
        schema_rows_by_position=schema_rows_by_position,
    )


def load_workflow_inputs(paths: WorkflowPaths) -> WorkflowInputs:
    shortlist_latest_pointer = load_json(paths.shortlist_latest_pointer)
    core_audit_latest_pointer = load_json(paths.core_audit_latest_pointer)

    shortlist_required_keys = {
        "shortlist_run_id",
        "core_audit_run_id",
        "parse_run_id",
        "source_run_id",
        "clinical_shortlist_tsv",
        "run_log_json",
    }
    core_required_keys = {
        "audit_run_id",
        "parse_run_id",
        "source_run_id",
        "clinical_core_field_audit_tsv",
        "run_log_json",
        "clinical_biotab_latest_json",
    }
    missing_shortlist_keys = shortlist_required_keys.difference(shortlist_latest_pointer.keys())
    missing_core_keys = core_required_keys.difference(core_audit_latest_pointer.keys())
    if missing_shortlist_keys:
        raise EndpointCrosswalkError(
            f"Clinical shortlist latest pointer is missing required keys {sorted(missing_shortlist_keys)}: "
            f"{paths.shortlist_latest_pointer}"
        )
    if missing_core_keys:
        raise EndpointCrosswalkError(
            f"Clinical core audit latest pointer is missing required keys {sorted(missing_core_keys)}: "
            f"{paths.core_audit_latest_pointer}"
        )

    shortlist_run_id = str(shortlist_latest_pointer["shortlist_run_id"])
    core_audit_run_id = str(core_audit_latest_pointer["audit_run_id"])
    if str(shortlist_latest_pointer["core_audit_run_id"]) != core_audit_run_id:
        raise EndpointCrosswalkError(
            "Shortlist latest pointer references a different core audit run than the current core audit latest pointer: "
            f"{shortlist_latest_pointer['core_audit_run_id']} vs {core_audit_run_id}"
        )
    if str(shortlist_latest_pointer["parse_run_id"]) != str(core_audit_latest_pointer["parse_run_id"]):
        raise EndpointCrosswalkError(
            "Shortlist and core audit latest pointers disagree on parse_run_id: "
            f"{shortlist_latest_pointer['parse_run_id']} vs {core_audit_latest_pointer['parse_run_id']}"
        )
    if str(shortlist_latest_pointer["source_run_id"]) != str(core_audit_latest_pointer["source_run_id"]):
        raise EndpointCrosswalkError(
            "Shortlist and core audit latest pointers disagree on source_run_id: "
            f"{shortlist_latest_pointer['source_run_id']} vs {core_audit_latest_pointer['source_run_id']}"
        )

    input_paths = {
        "shortlist_latest_json": paths.shortlist_latest_pointer,
        "clinical_shortlist_tsv": paths.repo_root / shortlist_latest_pointer["clinical_shortlist_tsv"],
        "shortlist_run_log_json": paths.repo_root / shortlist_latest_pointer["run_log_json"],
        "core_audit_latest_json": paths.core_audit_latest_pointer,
        "clinical_core_field_audit_tsv": paths.repo_root / core_audit_latest_pointer["clinical_core_field_audit_tsv"],
        "core_audit_run_log_json": paths.repo_root / core_audit_latest_pointer["run_log_json"],
        "clinical_biotab_latest_json": paths.repo_root / core_audit_latest_pointer["clinical_biotab_latest_json"],
    }

    parse_latest_pointer = load_json(input_paths["clinical_biotab_latest_json"])
    parse_required_keys = {
        "parse_run_id",
        "source_run_id",
        "table_manifest_tsv",
        "run_log_json",
    }
    missing_parse_keys = parse_required_keys.difference(parse_latest_pointer.keys())
    if missing_parse_keys:
        raise EndpointCrosswalkError(
            f"Clinical biotab latest pointer is missing required keys {sorted(missing_parse_keys)}: "
            f"{input_paths['clinical_biotab_latest_json']}"
        )
    if str(parse_latest_pointer["parse_run_id"]) != str(core_audit_latest_pointer["parse_run_id"]):
        raise EndpointCrosswalkError(
            "Core audit latest pointer and clinical biotab latest pointer disagree on parse_run_id: "
            f"{core_audit_latest_pointer['parse_run_id']} vs {parse_latest_pointer['parse_run_id']}"
        )
    if str(parse_latest_pointer["source_run_id"]) != str(core_audit_latest_pointer["source_run_id"]):
        raise EndpointCrosswalkError(
            "Core audit latest pointer and clinical biotab latest pointer disagree on source_run_id: "
            f"{core_audit_latest_pointer['source_run_id']} vs {parse_latest_pointer['source_run_id']}"
        )

    input_paths["clinical_biotab_manifest_tsv"] = paths.repo_root / parse_latest_pointer["table_manifest_tsv"]
    input_paths["clinical_biotab_run_log_json"] = paths.repo_root / parse_latest_pointer["run_log_json"]

    shortlist_run_log = load_json(input_paths["shortlist_run_log_json"])
    core_audit_run_log = load_json(input_paths["core_audit_run_log_json"])
    parse_run_log = load_json(input_paths["clinical_biotab_run_log_json"])
    for run_log, label in [
        (shortlist_run_log, "shortlist"),
        (core_audit_run_log, "core audit"),
        (parse_run_log, "clinical biotab parse"),
    ]:
        if run_log.get("status") != "completed":
            raise EndpointCrosswalkError(f"{label.title()} run log is not marked completed.")

    if str(shortlist_run_log.get("shortlist_run_id") or "") != shortlist_run_id:
        raise EndpointCrosswalkError(
            f"Shortlist run log ID does not match shortlist latest pointer: {shortlist_run_log.get('shortlist_run_id')}"
        )
    if str(core_audit_run_log.get("audit_run_id") or "") != core_audit_run_id:
        raise EndpointCrosswalkError(
            f"Core audit run log ID does not match core audit latest pointer: {core_audit_run_log.get('audit_run_id')}"
        )
    if str(parse_run_log.get("parse_run_id") or "") != str(parse_latest_pointer["parse_run_id"]):
        raise EndpointCrosswalkError(
            f"Clinical biotab run log parse_run_id does not match latest pointer: {parse_run_log.get('parse_run_id')}"
        )

    shortlist_rows = build_shortlist_rows(
        read_required_tsv(input_paths["clinical_shortlist_tsv"], "clinical shortlist TSV")
    )
    core_audit_rows = build_core_audit_rows(
        read_required_tsv(input_paths["clinical_core_field_audit_tsv"], "clinical core field audit TSV")
    )
    validate_unique_field_keys(shortlist_rows, "clinical shortlist TSV")
    validate_unique_field_keys(core_audit_rows, "clinical core field audit TSV")

    shortlist_field_count = parse_int(shortlist_latest_pointer.get("field_count"), "shortlist pointer field_count")
    core_field_count = parse_int(core_audit_latest_pointer.get("field_count"), "core audit pointer field_count")
    if shortlist_field_count != len(shortlist_rows):
        raise EndpointCrosswalkError(
            f"Shortlist TSV row count {len(shortlist_rows)} does not match pointer field_count {shortlist_field_count}."
        )
    if core_field_count != len(core_audit_rows):
        raise EndpointCrosswalkError(
            f"Core audit TSV row count {len(core_audit_rows)} does not match pointer field_count {core_field_count}."
        )
    if len(shortlist_rows) != len(core_audit_rows):
        raise EndpointCrosswalkError(
            f"Shortlist TSV row count {len(shortlist_rows)} does not match core audit TSV row count {len(core_audit_rows)}."
        )

    core_audit_by_key = {normalize_key(row): row for row in core_audit_rows}
    compare_shortlist_and_core_rows(shortlist_rows, core_audit_by_key)

    manifest_rows = read_required_tsv(input_paths["clinical_biotab_manifest_tsv"], "clinical biotab manifest TSV")
    required_manifest_columns = {
        "parse_run_id",
        "source_run_id",
        "table_name",
        "output_path",
        "schema_path",
        "row_count",
        "column_count",
    }
    if not required_manifest_columns.issubset(manifest_rows[0].keys()):
        raise EndpointCrosswalkError(
            f"Clinical biotab manifest TSV is missing required columns {sorted(required_manifest_columns)}."
        )

    manifest_rows_by_table = {str(row["table_name"]): row for row in manifest_rows}
    parsed_tables: dict[str, ParsedTableContext] = {}
    for table_name in TARGET_TABLES:
        if table_name not in manifest_rows_by_table:
            raise EndpointCrosswalkError(f"Clinical biotab manifest does not include required table: {table_name}")
        parsed_tables[table_name] = build_parsed_table_context(
            repo_root=paths.repo_root,
            manifest_row=manifest_rows_by_table[table_name],
            expected_parse_run_id=str(parse_latest_pointer["parse_run_id"]),
            expected_source_run_id=str(parse_latest_pointer["source_run_id"]),
        )

    return WorkflowInputs(
        shortlist_latest_pointer=shortlist_latest_pointer,
        shortlist_run_log=shortlist_run_log,
        core_audit_latest_pointer=core_audit_latest_pointer,
        core_audit_run_log=core_audit_run_log,
        parse_latest_pointer=parse_latest_pointer,
        parse_run_log=parse_run_log,
        shortlist_rows=shortlist_rows,
        core_audit_rows=core_audit_rows,
        core_audit_by_key=core_audit_by_key,
        parsed_tables=parsed_tables,
        input_paths=input_paths,
    )


def classify_endpoint_signal(table_name: str, field_name: str) -> dict[str, str] | None:
    normalized_name = normalize_field_name(field_name)

    if "vital_status" in normalized_name:
        return {
            "endpoint_signal_family": "survival_status_like",
            "endpoint_signal_type": "survival_status_signal",
            "endpoint_signal_rule": "signal_vital_status_name",
            "canonical_signal_key": normalized_name,
        }
    if "last_contact" in normalized_name:
        return {
            "endpoint_signal_family": "last_contact_like",
            "endpoint_signal_type": "last_contact_time_signal",
            "endpoint_signal_rule": "signal_last_contact_name",
            "canonical_signal_key": normalized_name,
        }
    if normalized_name == "death_days_to" or ("death" in normalized_name and "days_to" in normalized_name):
        return {
            "endpoint_signal_family": "death_time_like",
            "endpoint_signal_type": "death_time_signal",
            "endpoint_signal_rule": "signal_death_days_to_name",
            "canonical_signal_key": normalized_name,
        }
    if "tumor_status" in normalized_name:
        return {
            "endpoint_signal_family": "progression_or_tumor_status_like",
            "endpoint_signal_type": "tumor_status_signal",
            "endpoint_signal_rule": "signal_tumor_status_name",
            "canonical_signal_key": normalized_name,
        }
    if "progression" in normalized_name:
        return {
            "endpoint_signal_family": "progression_or_tumor_status_like",
            "endpoint_signal_type": "progression_time_signal",
            "endpoint_signal_rule": "signal_progression_name",
            "canonical_signal_key": normalized_name,
        }
    if "new_tumor_event" in normalized_name:
        return {
            "endpoint_signal_family": "new_tumor_event_like",
            "endpoint_signal_type": "new_tumor_event_signal",
            "endpoint_signal_rule": "signal_new_tumor_event_name",
            "canonical_signal_key": normalized_name,
        }
    if table_name == "clinical_follow_up_v4_0" and "lost_to" in normalized_name:
        return {
            "endpoint_signal_family": "followup_loss_like",
            "endpoint_signal_type": "followup_loss_signal",
            "endpoint_signal_rule": "signal_followup_lost_to_followup_table",
            "canonical_signal_key": normalized_name,
        }
    if any(
        token in normalized_name
        for token in (
            "progression_free",
            "disease_free",
            "recurrence",
            "relapse",
            "event_free",
            "event_indicator",
        )
    ):
        return {
            "endpoint_signal_family": "other_endpoint_like",
            "endpoint_signal_type": "other_endpoint_signal",
            "endpoint_signal_rule": "signal_other_endpoint_like_fallback",
            "canonical_signal_key": normalized_name,
        }
    return None


def validate_selected_field_against_parsed_context(
    *,
    row: dict[str, Any],
    parsed_table: ParsedTableContext,
) -> None:
    source_position = int(row["source_position"])
    if source_position < 1 or source_position > len(parsed_table.header):
        raise EndpointCrosswalkError(
            f"source_position {source_position} is out of bounds for {parsed_table.table_name}"
        )

    header_field_name = parsed_table.header[source_position - 1]
    if header_field_name != str(row["field_name"]):
        raise EndpointCrosswalkError(
            f"Parsed TSV header mismatch for {parsed_table.table_name} position {source_position}: "
            f"{header_field_name!r} vs {row['field_name']!r}"
        )

    schema_row = parsed_table.schema_rows_by_position.get(source_position)
    if schema_row is None:
        raise EndpointCrosswalkError(
            f"Schema TSV is missing position {source_position} for {parsed_table.table_name}"
        )

    if str(schema_row["raw_column_name"]) != str(row["field_name"]):
        raise EndpointCrosswalkError(
            f"Schema raw_column_name mismatch for {parsed_table.table_name} position {source_position}: "
            f"{schema_row['raw_column_name']!r} vs {row['field_name']!r}"
        )

    if str(schema_row["alternate_column_name"]) != str(row["alternate_column_name"]):
        raise EndpointCrosswalkError(
            f"Schema alternate_column_name mismatch for {parsed_table.table_name} position {source_position}: "
            f"{schema_row['alternate_column_name']!r} vs {row['alternate_column_name']!r}"
        )


def selection_rule_for_bucket(shortlist_bucket: str) -> tuple[str, str] | None:
    if shortlist_bucket == "usable_endpoint_candidate":
        return "usable_endpoint_candidate", "select_usable_endpoint_candidate"
    if shortlist_bucket == "unclear_manual_review":
        return "escalated_manual_review", "escalate_unclear_endpoint_like"
    if shortlist_bucket == "weak_or_unusable":
        return "escalated_manual_review", "escalate_weak_endpoint_like"
    return None


def build_inventory_rows(
    *,
    crosswalk_run_id: str,
    workflow_inputs: WorkflowInputs,
) -> list[dict[str, Any]]:
    shortlist_run_id = str(workflow_inputs.shortlist_latest_pointer["shortlist_run_id"])
    core_audit_run_id = str(workflow_inputs.core_audit_latest_pointer["audit_run_id"])
    parse_run_id = str(workflow_inputs.core_audit_latest_pointer["parse_run_id"])
    source_run_id = str(workflow_inputs.core_audit_latest_pointer["source_run_id"])

    inventory_rows: list[dict[str, Any]] = []
    for shortlist_row in sorted(
        workflow_inputs.shortlist_rows,
        key=lambda item: (TARGET_TABLE_ORDER.get(str(item["table_name"]), 99), int(item["source_position"])),
    ):
        table_name = str(shortlist_row["table_name"])
        if table_name not in TARGET_TABLES:
            continue

        signal_metadata = classify_endpoint_signal(table_name, str(shortlist_row["field_name"]))
        if signal_metadata is None:
            continue

        selection = selection_rule_for_bucket(str(shortlist_row["shortlist_bucket"]))
        if selection is None:
            continue
        candidate_source, selection_rule = selection

        core_row = workflow_inputs.core_audit_by_key.get(normalize_key(shortlist_row))
        if core_row is None:
            raise EndpointCrosswalkError(
                f"Selected endpoint-like field is missing from the core audit TSV: {normalize_key(shortlist_row)}"
            )

        validate_selected_field_against_parsed_context(
            row=shortlist_row,
            parsed_table=workflow_inputs.parsed_tables[table_name],
        )

        inventory_rows.append(
            {
                "crosswalk_run_id": crosswalk_run_id,
                "shortlist_run_id": shortlist_run_id,
                "core_audit_run_id": core_audit_run_id,
                "parse_run_id": parse_run_id,
                "source_run_id": source_run_id,
                "table_name": table_name,
                "field_name": str(shortlist_row["field_name"]),
                "source_position": int(shortlist_row["source_position"]),
                "alternate_column_name": str(shortlist_row["alternate_column_name"]),
                "probable_field_group": str(core_row["probable_field_group"]),
                "missing_like_fraction": round(float(core_row["missing_like_fraction"]), 6),
                "non_missing_count": int(core_row["non_missing_count"]),
                "distinct_non_missing_count": int(core_row["distinct_non_missing_count"]),
                "example_values_small_sample": str(core_row["example_values_small_sample"]),
                "candidate_source": candidate_source,
                "original_shortlist_bucket": str(shortlist_row["shortlist_bucket"]),
                "original_shortlist_rule": str(shortlist_row["shortlist_rule"]),
                "selection_rule": selection_rule,
                "endpoint_signal_family": signal_metadata["endpoint_signal_family"],
                "endpoint_signal_type": signal_metadata["endpoint_signal_type"],
                "endpoint_signal_rule": signal_metadata["endpoint_signal_rule"],
                "canonical_signal_key": signal_metadata["canonical_signal_key"],
                "manual_review_priority": str(shortlist_row["manual_review_priority"]),
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )

    validate_unique_field_keys(inventory_rows, "endpoint candidate inventory")
    return inventory_rows


def primary_rank_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        float(row["missing_like_fraction"]),
        -int(row["non_missing_count"]),
        PRIMARY_TABLE_RANK[str(row["table_name"])],
        int(row["source_position"]),
    )


def crosswalk_sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        FAMILY_ORDER.index(str(row["endpoint_signal_family"])),
        ROLE_ORDER.index(str(row["crosswalk_role"])),
        float(row["missing_like_fraction"]),
        -int(row["non_missing_count"]),
        TARGET_TABLE_ORDER.get(str(row["table_name"]), 99),
        int(row["source_position"]),
    )


def assign_crosswalk_rows(
    *,
    inventory_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows_by_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in inventory_rows:
        rows_by_family[str(row["endpoint_signal_family"])].append(row)

    crosswalk_rows: list[dict[str, Any]] = []
    for endpoint_signal_family in sorted(rows_by_family, key=FAMILY_ORDER.index):
        family_rows = rows_by_family[endpoint_signal_family]
        usable_rows = [
            row for row in family_rows if str(row["original_shortlist_bucket"]) == "usable_endpoint_candidate"
        ]
        primary_key: tuple[str, str] | None = None
        primary_canonical_signal_key: str | None = None
        if usable_rows:
            primary_row = sorted(usable_rows, key=primary_rank_key)[0]
            primary_key = normalize_key(primary_row)
            primary_canonical_signal_key = str(primary_row["canonical_signal_key"])

        family_crosswalk_rows: list[dict[str, Any]] = []
        for row in family_rows:
            row_key = normalize_key(row)
            original_shortlist_bucket = str(row["original_shortlist_bucket"])

            if primary_key is not None and row_key == primary_key:
                crosswalk_role = "primary_candidate"
                crosswalk_rule = "role_primary_best_usable_completeness"
            elif primary_key is not None and str(row["canonical_signal_key"]) == primary_canonical_signal_key:
                crosswalk_role = "overlapping_candidate"
                crosswalk_rule = "role_overlapping_same_canonical_signal"
            elif primary_key is not None and original_shortlist_bucket == "usable_endpoint_candidate":
                crosswalk_role = "secondary_candidate"
                crosswalk_rule = "role_secondary_distinct_usable_signal"
            else:
                crosswalk_role = "ambiguous_candidate"
                crosswalk_rule = "role_ambiguous_no_usable_primary_or_sparse_signal"

            family_crosswalk_rows.append(
                {
                    "crosswalk_run_id": row["crosswalk_run_id"],
                    "shortlist_run_id": row["shortlist_run_id"],
                    "core_audit_run_id": row["core_audit_run_id"],
                    "parse_run_id": row["parse_run_id"],
                    "source_run_id": row["source_run_id"],
                    "endpoint_signal_family": row["endpoint_signal_family"],
                    "endpoint_signal_type": row["endpoint_signal_type"],
                    "endpoint_signal_rule": row["endpoint_signal_rule"],
                    "canonical_signal_key": row["canonical_signal_key"],
                    "table_name": row["table_name"],
                    "field_name": row["field_name"],
                    "source_position": row["source_position"],
                    "alternate_column_name": row["alternate_column_name"],
                    "candidate_source": row["candidate_source"],
                    "original_shortlist_bucket": row["original_shortlist_bucket"],
                    "original_shortlist_rule": row["original_shortlist_rule"],
                    "selection_rule": row["selection_rule"],
                    "missing_like_fraction": row["missing_like_fraction"],
                    "non_missing_count": row["non_missing_count"],
                    "distinct_non_missing_count": row["distinct_non_missing_count"],
                    "example_values_small_sample": row["example_values_small_sample"],
                    "crosswalk_role": crosswalk_role,
                    "crosswalk_rule": crosswalk_rule,
                    "manual_review_priority": "pending",
                    "notes_placeholder": NOTES_PLACEHOLDER,
                }
            )

        family_row_count = len(family_crosswalk_rows)
        for row in family_crosswalk_rows:
            crosswalk_role = str(row["crosswalk_role"])
            if crosswalk_role in {"overlapping_candidate", "ambiguous_candidate"}:
                row["manual_review_priority"] = "high"
            elif crosswalk_role == "secondary_candidate":
                row["manual_review_priority"] = "medium"
            elif family_row_count > 1:
                row["manual_review_priority"] = "medium"
            else:
                row["manual_review_priority"] = "low"

        crosswalk_rows.extend(sorted(family_crosswalk_rows, key=crosswalk_sort_key))

    return crosswalk_rows


def build_summary_rows(crosswalk_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows_by_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in crosswalk_rows:
        rows_by_family[str(row["endpoint_signal_family"])].append(row)

    summary_rows: list[dict[str, Any]] = []
    for endpoint_signal_family in sorted(rows_by_family, key=FAMILY_ORDER.index):
        family_rows = sorted(rows_by_family[endpoint_signal_family], key=crosswalk_sort_key)

        def rows_to_json_strings(role_name: str) -> str:
            values = sorted(
                f"{row['table_name']}.{row['field_name']}"
                for row in family_rows
                if str(row["crosswalk_role"]) == role_name
            )
            return json.dumps(values, ensure_ascii=True)

        summary_rows.append(
            {
                "crosswalk_run_id": str(family_rows[0]["crosswalk_run_id"]),
                "shortlist_run_id": str(family_rows[0]["shortlist_run_id"]),
                "core_audit_run_id": str(family_rows[0]["core_audit_run_id"]),
                "parse_run_id": str(family_rows[0]["parse_run_id"]),
                "source_run_id": str(family_rows[0]["source_run_id"]),
                "endpoint_signal_family": endpoint_signal_family,
                "field_count": len(family_rows),
                "tables_present_json": json.dumps(
                    sorted({str(row["table_name"]) for row in family_rows}),
                    ensure_ascii=True,
                ),
                "primary_candidates_json": rows_to_json_strings("primary_candidate"),
                "secondary_candidates_json": rows_to_json_strings("secondary_candidate"),
                "overlapping_candidates_json": rows_to_json_strings("overlapping_candidate"),
                "ambiguous_candidates_json": rows_to_json_strings("ambiguous_candidate"),
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )

    return summary_rows


def write_failure_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    crosswalk_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    crosswalk_run_dir = paths.crosswalk_runs_root / crosswalk_run_id
    run_log_path = crosswalk_run_dir / "run_log.json"

    try:
        trial_config = load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths)
        create_run_directory(crosswalk_run_dir)

        inventory_rows = build_inventory_rows(
            crosswalk_run_id=crosswalk_run_id,
            workflow_inputs=workflow_inputs,
        )
        if not inventory_rows:
            raise EndpointCrosswalkError("No endpoint-candidate fields were selected for the crosswalk inventory.")

        usable_candidate_count = sum(
            str(row["candidate_source"]) == "usable_endpoint_candidate" for row in inventory_rows
        )
        if usable_candidate_count <= 0:
            raise EndpointCrosswalkError(
                "At least one usable_endpoint_candidate field is required for the endpoint crosswalk workflow."
            )

        inventory_path = crosswalk_run_dir / "endpoint_candidate_inventory.tsv"
        write_dict_rows_tsv(inventory_path, INVENTORY_FIELDNAMES, inventory_rows)

        crosswalk_rows = assign_crosswalk_rows(inventory_rows=inventory_rows)
        if not crosswalk_rows:
            raise EndpointCrosswalkError("Endpoint crosswalk rows were not generated.")
        crosswalk_path = crosswalk_run_dir / "endpoint_crosswalk.tsv"
        write_dict_rows_tsv(crosswalk_path, CROSSWALK_FIELDNAMES, crosswalk_rows)

        summary_rows = build_summary_rows(crosswalk_rows)
        if not summary_rows:
            raise EndpointCrosswalkError("Endpoint crosswalk summary rows were not generated.")
        summary_path = crosswalk_run_dir / "endpoint_crosswalk_summary.tsv"
        write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        family_values = {str(row["endpoint_signal_family"]) for row in crosswalk_rows}
        summary_family_values = {str(row["endpoint_signal_family"]) for row in summary_rows}
        summary_total = sum(parse_int(row["field_count"], "summary field_count") for row in summary_rows)
        inventory_rules_valid = all(
            str(row["selection_rule"]) in SELECTION_RULES and str(row["candidate_source"]) in CANDIDATE_SOURCES
            for row in inventory_rows
        )
        crosswalk_rules_valid = all(
            str(row["crosswalk_rule"]) in CROSSWALK_RULES and str(row["crosswalk_role"]) in ROLE_ORDER
            for row in crosswalk_rows
        )
        inventory_field_headers_validated = True
        output_rows_positive = all(len(rows) > 0 for rows in (inventory_rows, crosswalk_rows, summary_rows))
        summary_matches_crosswalk = family_values == summary_family_values and summary_total == len(crosswalk_rows)

        shortlist_run_id = str(workflow_inputs.shortlist_latest_pointer["shortlist_run_id"])
        core_audit_run_id = str(workflow_inputs.core_audit_latest_pointer["audit_run_id"])
        parse_run_id = str(workflow_inputs.core_audit_latest_pointer["parse_run_id"])
        source_run_id = str(workflow_inputs.core_audit_latest_pointer["source_run_id"])

        latest_pointer_payload = {
            "updated_at_utc": format_utc_timestamp(utc_now()),
            "crosswalk_run_id": crosswalk_run_id,
            "shortlist_run_id": shortlist_run_id,
            "core_audit_run_id": core_audit_run_id,
            "parse_run_id": parse_run_id,
            "source_run_id": source_run_id,
            "crosswalk_run_directory": repo_relative(crosswalk_run_dir, paths.repo_root),
            "endpoint_candidate_inventory_tsv": repo_relative(inventory_path, paths.repo_root),
            "endpoint_crosswalk_tsv": repo_relative(crosswalk_path, paths.repo_root),
            "endpoint_crosswalk_summary_tsv": repo_relative(summary_path, paths.repo_root),
            "run_log_json": repo_relative(run_log_path, paths.repo_root),
            "shortlist_latest_json": repo_relative(paths.shortlist_latest_pointer, paths.repo_root),
            "core_audit_latest_json": repo_relative(paths.core_audit_latest_pointer, paths.repo_root),
            "field_count": len(inventory_rows),
        }

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "crosswalk_run_id": crosswalk_run_id,
            "shortlist_run_id": shortlist_run_id,
            "core_audit_run_id": core_audit_run_id,
            "parse_run_id": parse_run_id,
            "source_run_id": source_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "shortlist_latest_json": repo_relative(paths.shortlist_latest_pointer, paths.repo_root),
                "clinical_shortlist_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_tsv"],
                    paths.repo_root,
                ),
                "shortlist_run_log_json": repo_relative(
                    workflow_inputs.input_paths["shortlist_run_log_json"],
                    paths.repo_root,
                ),
                "core_audit_latest_json": repo_relative(paths.core_audit_latest_pointer, paths.repo_root),
                "clinical_core_field_audit_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_core_field_audit_tsv"],
                    paths.repo_root,
                ),
                "core_audit_run_log_json": repo_relative(
                    workflow_inputs.input_paths["core_audit_run_log_json"],
                    paths.repo_root,
                ),
                "clinical_biotab_latest_json": repo_relative(
                    workflow_inputs.input_paths["clinical_biotab_latest_json"],
                    paths.repo_root,
                ),
                "clinical_biotab_manifest_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_biotab_manifest_tsv"],
                    paths.repo_root,
                ),
                "clinical_biotab_run_log_json": repo_relative(
                    workflow_inputs.input_paths["clinical_biotab_run_log_json"],
                    paths.repo_root,
                ),
                "parsed_tables": {
                    table_name: {
                        "parsed_tsv": repo_relative(table_context.parsed_path, paths.repo_root),
                        "schema_tsv": repo_relative(table_context.schema_path, paths.repo_root),
                    }
                    for table_name, table_context in workflow_inputs.parsed_tables.items()
                },
            },
            "outputs": {
                "crosswalk_run_directory": repo_relative(crosswalk_run_dir, paths.repo_root),
                "endpoint_candidate_inventory_tsv": repo_relative(inventory_path, paths.repo_root),
                "endpoint_crosswalk_tsv": repo_relative(crosswalk_path, paths.repo_root),
                "endpoint_crosswalk_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    usable_candidate_count > 0
                    and output_rows_positive
                    and inventory_rules_valid
                    and crosswalk_rules_valid
                    and summary_matches_crosswalk
                    and inventory_field_headers_validated
                ),
                "shortlist_latest_pointer_found": True,
                "core_audit_latest_pointer_found": True,
                "shortlist_run_log_completed": True,
                "core_audit_run_log_completed": True,
                "clinical_biotab_run_log_completed": True,
                "shortlist_and_core_audit_run_ids_match": (
                    workflow_inputs.shortlist_latest_pointer["core_audit_run_id"]
                    == workflow_inputs.core_audit_latest_pointer["audit_run_id"]
                ),
                "patient_and_followup_parsed_tables_found": all(
                    table_name in workflow_inputs.parsed_tables for table_name in TARGET_TABLES
                ),
                "selected_fields_validated_against_headers_and_schema": inventory_field_headers_validated,
                "at_least_one_usable_endpoint_candidate_found": usable_candidate_count > 0,
                "inventory_row_count_positive": len(inventory_rows) > 0,
                "crosswalk_row_count_positive": len(crosswalk_rows) > 0,
                "summary_row_count_positive": len(summary_rows) > 0,
                "inventory_rules_valid": inventory_rules_valid,
                "crosswalk_rules_valid": crosswalk_rules_valid,
                "summary_family_counts_match_crosswalk": summary_matches_crosswalk,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "target_tables": list(TARGET_TABLES),
                "family_order": list(FAMILY_ORDER),
                "role_order": list(ROLE_ORDER),
                "candidate_sources": list(CANDIDATE_SOURCES),
                "selection_rules": list(SELECTION_RULES),
                "crosswalk_rules": list(CROSSWALK_RULES),
                "endpoint_signal_rule_order": [
                    {
                        "rule": "signal_vital_status_name",
                        "match": "field_name contains 'vital_status'",
                        "endpoint_signal_family": "survival_status_like",
                        "endpoint_signal_type": "survival_status_signal",
                    },
                    {
                        "rule": "signal_last_contact_name",
                        "match": "field_name contains 'last_contact'",
                        "endpoint_signal_family": "last_contact_like",
                        "endpoint_signal_type": "last_contact_time_signal",
                    },
                    {
                        "rule": "signal_death_days_to_name",
                        "match": "field_name == 'death_days_to' or contains both 'death' and 'days_to'",
                        "endpoint_signal_family": "death_time_like",
                        "endpoint_signal_type": "death_time_signal",
                    },
                    {
                        "rule": "signal_tumor_status_name",
                        "match": "field_name contains 'tumor_status'",
                        "endpoint_signal_family": "progression_or_tumor_status_like",
                        "endpoint_signal_type": "tumor_status_signal",
                    },
                    {
                        "rule": "signal_progression_name",
                        "match": "field_name contains 'progression'",
                        "endpoint_signal_family": "progression_or_tumor_status_like",
                        "endpoint_signal_type": "progression_time_signal",
                    },
                    {
                        "rule": "signal_new_tumor_event_name",
                        "match": "field_name contains 'new_tumor_event'",
                        "endpoint_signal_family": "new_tumor_event_like",
                        "endpoint_signal_type": "new_tumor_event_signal",
                    },
                    {
                        "rule": "signal_followup_lost_to_followup_table",
                        "match": "table_name == 'clinical_follow_up_v4_0' and field_name contains 'lost_to'",
                        "endpoint_signal_family": "followup_loss_like",
                        "endpoint_signal_type": "followup_loss_signal",
                    },
                    {
                        "rule": "signal_other_endpoint_like_fallback",
                        "match": "field_name contains disease-free, progression-free, recurrence, relapse, or event-free tokens",
                        "endpoint_signal_family": "other_endpoint_like",
                        "endpoint_signal_type": "other_endpoint_signal",
                    },
                ],
                "primary_candidate_ranking": [
                    "lowest_missing_like_fraction",
                    "highest_non_missing_count",
                    "clinical_follow_up_v4_0_before_clinical_patient",
                    "lowest_source_position",
                ],
                "canonical_signal_key_strategy": "normalized exact field_name",
                "manual_review_priority_rules": {
                    "high": ["overlapping_candidate", "ambiguous_candidate"],
                    "medium": ["secondary_candidate", "primary_candidate when family_row_count > 1"],
                    "low": ["primary_candidate when family_row_count == 1"],
                },
                "notes_placeholder": NOTES_PLACEHOLDER,
            },
            "inventory_counts_by_candidate_source": {
                source_name: sum(str(row["candidate_source"]) == source_name for row in inventory_rows)
                for source_name in CANDIDATE_SOURCES
            },
            "crosswalk_counts_by_role": {
                role_name: sum(str(row["crosswalk_role"]) == role_name for row in crosswalk_rows)
                for role_name in ROLE_ORDER
            },
            "crosswalk_counts_by_family": {
                family_name: sum(str(row["endpoint_signal_family"]) == family_name for row in crosswalk_rows)
                for family_name in family_values
            },
            "crosswalk_counts_by_family_and_role": {
                family_name: {
                    role_name: sum(
                        str(row["endpoint_signal_family"]) == family_name and str(row["crosswalk_role"]) == role_name
                        for row in crosswalk_rows
                    )
                    for role_name in ROLE_ORDER
                }
                for family_name in sorted(family_values, key=FAMILY_ORDER.index)
            },
            "tables": {
                table_name: {
                    "parsed_tsv": repo_relative(table_context.parsed_path, paths.repo_root),
                    "schema_tsv": repo_relative(table_context.schema_path, paths.repo_root),
                    "row_count": table_context.row_count,
                    "column_count": table_context.column_count,
                }
                for table_name, table_context in workflow_inputs.parsed_tables.items()
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "shortlist_latest_pointer": workflow_inputs.shortlist_latest_pointer,
                "shortlist_run_log": workflow_inputs.shortlist_run_log,
                "core_audit_latest_pointer": workflow_inputs.core_audit_latest_pointer,
                "core_audit_run_log": workflow_inputs.core_audit_run_log,
                "clinical_biotab_latest_pointer": workflow_inputs.parse_latest_pointer,
            },
        }

        write_json(run_log_path, run_log_payload)
        write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "crosswalk_run_id": crosswalk_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_endpoint_crosswalk",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA endpoint crosswalk audit complete.")
    print(f"Crosswalk run ID: {run_log['crosswalk_run_id']}")
    print(f"Shortlist run ID: {run_log['shortlist_run_id']}")
    print(f"Core audit run ID: {run_log['core_audit_run_id']}")
    print(f"Parse run ID: {run_log['parse_run_id']}")
    print(f"Source run ID: {run_log['source_run_id']}")
    print(f"Crosswalk output directory: {run_log['outputs']['crosswalk_run_directory']}")
    print(f"Inventory TSV: {run_log['outputs']['endpoint_candidate_inventory_tsv']}")
    print(f"Crosswalk TSV: {run_log['outputs']['endpoint_crosswalk_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['endpoint_crosswalk_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
