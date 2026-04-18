#!/usr/bin/env python
"""Build an auditable biospecimen identifier crosswalk from parsed TCGA-BRCA biospecimen biotabs."""

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

IDENTIFIER_TOKENS = (
    "patient",
    "sample",
    "portion",
    "analyte",
    "slide",
    "aliquot",
    "barcode",
    "uuid",
)
LEVEL_TOKENS = (
    "patient",
    "sample",
    "portion",
    "analyte",
    "slide",
    "aliquot",
)
GENERIC_IDENTIFIER_TOKENS = ("barcode", "uuid")
CENTRAL_TABLES = (
    "biospecimen_sample",
    "biospecimen_portion",
    "biospecimen_analyte",
    "biospecimen_slide",
    "biospecimen_aliquot",
)
OPTIONAL_SIDE_TABLES = (
    "biospecimen_protocol",
    "biospecimen_shipment_portion",
    "biospecimen_diagnostic_slides",
)
ALL_INCLUDED_TABLES = CENTRAL_TABLES + OPTIONAL_SIDE_TABLES
TABLE_ORDER = {table_name: index for index, table_name in enumerate(ALL_INCLUDED_TABLES)}
PATIENT_PRIMARY_TABLE_ORDER = {
    "biospecimen_sample": 0,
    "biospecimen_portion": 1,
    "biospecimen_analyte": 2,
    "biospecimen_slide": 3,
    "biospecimen_aliquot": 4,
    "biospecimen_diagnostic_slides": 5,
    "biospecimen_protocol": 6,
    "biospecimen_shipment_portion": 7,
}
CANONICAL_TABLE_BY_LEVEL = {
    "sample": "biospecimen_sample",
    "portion": "biospecimen_portion",
    "analyte": "biospecimen_analyte",
    "slide": "biospecimen_slide",
    "aliquot": "biospecimen_aliquot",
}
FAMILY_ORDER = (
    "patient_identifier_like",
    "sample_identifier_like",
    "portion_identifier_like",
    "analyte_identifier_like",
    "slide_identifier_like",
    "aliquot_identifier_like",
    "multi_level_or_unclear_identifier_like",
)
ROLE_ORDER = (
    "primary_link_candidate",
    "overlapping_link_candidate",
    "secondary_link_candidate",
    "side_table_candidate",
    "ambiguous_candidate",
)
CANDIDATE_SOURCES = (
    "field_name_level_and_identifier_token_match",
    "field_name_level_token_match_only",
    "field_name_identifier_token_match_only",
    "field_name_multi_level_token_match",
)
CROSSWALK_RULES = (
    "role_primary_patient_best_uuid_candidate",
    "role_primary_patient_best_barcode_candidate",
    "role_primary_canonical_level_identifier",
    "role_overlapping_repeated_strong_identifier",
    "role_secondary_level_helper_field",
    "role_side_table_strong_identifier",
    "role_ambiguous_generic_or_conflicted_identifier",
)
PATTERN_RULES = (
    "exact_pair_one_to_one_same_table",
    "barcode_prefix_same_row",
    "same_row_copresence_only",
    "not_applicable_no_pair_rule",
)
MISSING_LIKE_VALUES = {
    "",
    "[Not Available]",
    "[Not Applicable]",
    "[Unknown]",
    "[Not Evaluated]",
    "[Discrepancy]",
}
FIELD_SUBTYPE_STRONG = {"level_barcode", "level_uuid"}
FIELD_SUBTYPE_HELPER = {"level_helper", "level_context"}
NOTES_PLACEHOLDER = "[fill in during biospecimen identifier crosswalk review]"
INVENTORY_FIELDNAMES = [
    "crosswalk_run_id",
    "parse_run_id",
    "source_run_id",
    "table_name",
    "field_name",
    "source_position",
    "probable_identifier_level",
    "identifier_token_family",
    "missing_like_fraction",
    "non_missing_count",
    "distinct_non_missing_count",
    "example_values_small_sample",
    "candidate_source",
    "manual_review_priority",
    "notes_placeholder",
    "matched_identifier_tokens_json",
    "alternate_column_name",
    "cde_id_raw",
]
CROSSWALK_FIELDNAMES = [
    "crosswalk_run_id",
    "parse_run_id",
    "source_run_id",
    "identifier_linkage_family",
    "table_name",
    "field_name",
    "source_position",
    "probable_identifier_level",
    "identifier_token_family",
    "candidate_source",
    "missing_like_fraction",
    "non_missing_count",
    "distinct_non_missing_count",
    "example_values_small_sample",
    "crosswalk_role",
    "crosswalk_rule",
    "manual_review_priority",
    "notes_placeholder",
    "matched_identifier_tokens_json",
    "alternate_column_name",
    "cde_id_raw",
    "paired_field_name",
    "pattern_rule",
    "pattern_comparable_row_count",
    "pattern_pass_fraction",
    "pattern_nonmatch_examples_json",
]
SUMMARY_FIELDNAMES = [
    "crosswalk_run_id",
    "parse_run_id",
    "source_run_id",
    "identifier_linkage_family",
    "field_count",
    "tables_present_json",
    "primary_candidates_json",
    "secondary_candidates_json",
    "overlapping_candidates_json",
    "ambiguous_candidates_json",
    "side_table_candidates_json",
    "notes_placeholder",
]


class BiospecimenIdentifierCrosswalkError(RuntimeError):
    """Raised when the biospecimen identifier crosswalk workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the biospecimen identifier crosswalk workflow."""

    repo_root: Path
    trial_config: Path
    parse_latest_pointer: Path
    crosswalk_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class ParsedTableContext:
    """One parsed biospecimen table plus its schema and rows."""

    table_name: str
    source_filename: str
    parsed_path: Path
    schema_path: Path
    row_count: int
    column_count: int
    header: list[str]
    schema_rows_by_position: dict[int, dict[str, Any]]
    data_rows: list[dict[str, str]]


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from the latest biospecimen parse layer."""

    parse_latest_pointer: dict[str, Any]
    parse_run_log: dict[str, Any]
    source_inventory_rows: list[dict[str, str]]
    manifest_rows: list[dict[str, str]]
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
    raise BiospecimenIdentifierCrosswalkError("Unable to locate the repository root from the script path.")


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
                raise BiospecimenIdentifierCrosswalkError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise BiospecimenIdentifierCrosswalkError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise BiospecimenIdentifierCrosswalkError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise BiospecimenIdentifierCrosswalkError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise BiospecimenIdentifierCrosswalkError(
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
        raise BiospecimenIdentifierCrosswalkError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise BiospecimenIdentifierCrosswalkError(f"Expected a JSON object in file: {path}")
    return data


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise BiospecimenIdentifierCrosswalkError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise BiospecimenIdentifierCrosswalkError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise BiospecimenIdentifierCrosswalkError(f"TSV file is empty: {path}") from exc
    if not header:
        raise BiospecimenIdentifierCrosswalkError(f"TSV header is empty: {path}")
    return header


def read_tsv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def read_required_tsv(path: Path, label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise BiospecimenIdentifierCrosswalkError(f"Required {label} does not exist: {path}")
    rows = read_tsv_dict_rows(path)
    if not rows:
        raise BiospecimenIdentifierCrosswalkError(f"Required {label} has no rows: {path}")
    return rows


def parse_int(value: Any, label: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise BiospecimenIdentifierCrosswalkError(f"Expected integer value for {label}: {value!r}") from exc


def parse_float(value: Any, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise BiospecimenIdentifierCrosswalkError(f"Expected float value for {label}: {value!r}") from exc


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise BiospecimenIdentifierCrosswalkError(f"Run directory already exists: {path}")
    path.mkdir(parents=False, exist_ok=False)
    return path


def build_workflow_paths() -> WorkflowPaths:
    repo_root = detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    trial_config_data = load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        parse_latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_biospecimen_biotabs_latest.json",
        crosswalk_runs_root=audit_root / "tcga-brca" / "variables" / "biospecimen_identifier_crosswalk_runs",
        latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_biospecimen_identifier_crosswalk_latest.json"
        ),
    )


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
        raise BiospecimenIdentifierCrosswalkError(
            f"Manifest parse_run_id mismatch for {table_name}: {parse_run_id} vs {expected_parse_run_id}"
        )
    if source_run_id != expected_source_run_id:
        raise BiospecimenIdentifierCrosswalkError(
            f"Manifest source_run_id mismatch for {table_name}: {source_run_id} vs {expected_source_run_id}"
        )

    parsed_path = repo_root / str(manifest_row.get("output_path") or "")
    schema_path = repo_root / str(manifest_row.get("schema_path") or "")
    source_filename = str(manifest_row.get("source_filename") or "")
    row_count = parse_int(manifest_row.get("row_count"), f"{table_name} row_count")
    column_count = parse_int(manifest_row.get("column_count"), f"{table_name} column_count")

    for required_path, label in [
        (parsed_path, f"{table_name} parsed TSV"),
        (schema_path, f"{table_name} schema TSV"),
    ]:
        if not required_path.exists():
            raise BiospecimenIdentifierCrosswalkError(f"Expected {label} does not exist: {required_path}")

    header = read_tsv_header(parsed_path)
    if len(header) != column_count:
        raise BiospecimenIdentifierCrosswalkError(
            f"Parsed TSV header length does not match manifest column_count for {table_name}: "
            f"{len(header)} vs {column_count}"
        )

    data_rows = read_tsv_dict_rows(parsed_path)
    if len(data_rows) != row_count:
        raise BiospecimenIdentifierCrosswalkError(
            f"Parsed TSV row count does not match manifest row_count for {table_name}: {len(data_rows)} vs {row_count}"
        )

    schema_rows = read_required_tsv(schema_path, f"{table_name} schema TSV")
    if len(schema_rows) != column_count:
        raise BiospecimenIdentifierCrosswalkError(
            f"Schema TSV row count does not match manifest column_count for {table_name}: "
            f"{len(schema_rows)} vs {column_count}"
        )

    schema_rows_by_position: dict[int, dict[str, Any]] = {}
    for schema_row in schema_rows:
        row_table_name = str(schema_row.get("table_name") or "").strip()
        if row_table_name != table_name:
            raise BiospecimenIdentifierCrosswalkError(
                f"Schema TSV table_name mismatch for {table_name}: found {row_table_name}"
            )
        position = parse_int(schema_row.get("column_position"), f"{table_name} schema column_position")
        if position in schema_rows_by_position:
            raise BiospecimenIdentifierCrosswalkError(
                f"Duplicate schema column_position {position} for {table_name}"
            )
        schema_rows_by_position[position] = {
            **schema_row,
            "column_position": position,
        }

    expected_positions = list(range(1, column_count + 1))
    if sorted(schema_rows_by_position) != expected_positions:
        raise BiospecimenIdentifierCrosswalkError(f"Schema TSV positions are incomplete for {table_name}")

    return ParsedTableContext(
        table_name=table_name,
        source_filename=source_filename,
        parsed_path=parsed_path,
        schema_path=schema_path,
        row_count=row_count,
        column_count=column_count,
        header=header,
        schema_rows_by_position=schema_rows_by_position,
        data_rows=data_rows,
    )


def load_workflow_inputs(paths: WorkflowPaths) -> WorkflowInputs:
    parse_latest_pointer = load_json(paths.parse_latest_pointer)
    required_keys = {
        "parse_run_id",
        "source_run_id",
        "source_inventory_tsv",
        "table_manifest_tsv",
        "run_log_json",
    }
    missing_keys = required_keys.difference(parse_latest_pointer.keys())
    if missing_keys:
        raise BiospecimenIdentifierCrosswalkError(
            f"Biospecimen parse latest pointer is missing required keys {sorted(missing_keys)}: "
            f"{paths.parse_latest_pointer}"
        )

    input_paths = {
        "parse_latest_json": paths.parse_latest_pointer,
        "source_inventory_tsv": paths.repo_root / parse_latest_pointer["source_inventory_tsv"],
        "manifest_tsv": paths.repo_root / parse_latest_pointer["table_manifest_tsv"],
        "parse_run_log_json": paths.repo_root / parse_latest_pointer["run_log_json"],
    }

    parse_run_log = load_json(input_paths["parse_run_log_json"])
    if parse_run_log.get("status") != "completed":
        raise BiospecimenIdentifierCrosswalkError("Biospecimen parse run log is not marked completed.")
    if str(parse_run_log.get("parse_run_id") or "") != str(parse_latest_pointer["parse_run_id"]):
        raise BiospecimenIdentifierCrosswalkError(
            "Biospecimen parse run log parse_run_id does not match latest pointer: "
            f"{parse_run_log.get('parse_run_id')} vs {parse_latest_pointer['parse_run_id']}"
        )
    if str(parse_run_log.get("source_run_id") or "") != str(parse_latest_pointer["source_run_id"]):
        raise BiospecimenIdentifierCrosswalkError(
            "Biospecimen parse run log source_run_id does not match latest pointer: "
            f"{parse_run_log.get('source_run_id')} vs {parse_latest_pointer['source_run_id']}"
        )

    source_inventory_rows = read_required_tsv(input_paths["source_inventory_tsv"], "biospecimen source inventory TSV")
    manifest_rows = read_required_tsv(input_paths["manifest_tsv"], "biospecimen table manifest TSV")

    manifest_required_columns = {
        "parse_run_id",
        "source_run_id",
        "table_name",
        "source_filename",
        "output_path",
        "schema_path",
        "row_count",
        "column_count",
    }
    if not manifest_required_columns.issubset(manifest_rows[0].keys()):
        raise BiospecimenIdentifierCrosswalkError(
            f"Biospecimen manifest TSV is missing required columns {sorted(manifest_required_columns)}."
        )

    manifest_rows_by_table = {str(row["table_name"]): row for row in manifest_rows}
    missing_required_tables = [table_name for table_name in CENTRAL_TABLES if table_name not in manifest_rows_by_table]
    if missing_required_tables:
        raise BiospecimenIdentifierCrosswalkError(
            f"Biospecimen manifest TSV does not include required central tables: {missing_required_tables}"
        )

    parsed_tables: dict[str, ParsedTableContext] = {}
    for table_name in ALL_INCLUDED_TABLES:
        manifest_row = manifest_rows_by_table.get(table_name)
        if manifest_row is None:
            continue
        parsed_tables[table_name] = build_parsed_table_context(
            repo_root=paths.repo_root,
            manifest_row=manifest_row,
            expected_parse_run_id=str(parse_latest_pointer["parse_run_id"]),
            expected_source_run_id=str(parse_latest_pointer["source_run_id"]),
        )

    if not all(table_name in parsed_tables for table_name in CENTRAL_TABLES):
        raise BiospecimenIdentifierCrosswalkError(
            "Failed to load all required central biospecimen parsed tables after manifest validation."
        )

    return WorkflowInputs(
        parse_latest_pointer=parse_latest_pointer,
        parse_run_log=parse_run_log,
        source_inventory_rows=source_inventory_rows,
        manifest_rows=manifest_rows,
        parsed_tables=parsed_tables,
        input_paths=input_paths,
    )


def normalize_field_name(field_name: str) -> str:
    return field_name.strip().lower()


def find_matched_identifier_tokens(field_name: str) -> list[str]:
    normalized_name = normalize_field_name(field_name)
    return [token for token in IDENTIFIER_TOKENS if token in normalized_name]


def classify_probable_identifier_level(matched_tokens: list[str]) -> str:
    matched_level_tokens = [token for token in matched_tokens if token in LEVEL_TOKENS]
    if len(matched_level_tokens) == 1:
        return matched_level_tokens[0]
    return "multi_level_or_unclear"


def classify_identifier_token_family(matched_tokens: list[str]) -> str:
    level_tokens = [token for token in matched_tokens if token in LEVEL_TOKENS]
    generic_tokens = [token for token in matched_tokens if token in GENERIC_IDENTIFIER_TOKENS]
    if "barcode" in generic_tokens and "uuid" in generic_tokens:
        return "barcode_and_uuid"
    if "barcode" in generic_tokens:
        return "barcode"
    if "uuid" in generic_tokens:
        return "uuid"
    if len(level_tokens) > 1:
        return "multi_level"
    if level_tokens:
        return "level_only"
    return "identifier_like_other"


def classify_candidate_source(matched_tokens: list[str]) -> str:
    level_tokens = [token for token in matched_tokens if token in LEVEL_TOKENS]
    generic_tokens = [token for token in matched_tokens if token in GENERIC_IDENTIFIER_TOKENS]
    if len(level_tokens) > 1:
        return "field_name_multi_level_token_match"
    if len(level_tokens) == 1 and generic_tokens:
        return "field_name_level_and_identifier_token_match"
    if len(level_tokens) == 1:
        return "field_name_level_token_match_only"
    if generic_tokens:
        return "field_name_identifier_token_match_only"
    raise BiospecimenIdentifierCrosswalkError(
        f"Cannot classify candidate_source without identifier-like tokens: {matched_tokens}"
    )


def classify_field_subtype(field_name: str, matched_tokens: list[str], probable_level: str) -> str:
    normalized_name = normalize_field_name(field_name)
    generic_tokens = [token for token in matched_tokens if token in GENERIC_IDENTIFIER_TOKENS]
    matched_level_tokens = [token for token in matched_tokens if token in LEVEL_TOKENS]

    if len(matched_level_tokens) > 1:
        return "multi_level_identifier"
    if probable_level == "multi_level_or_unclear" and generic_tokens:
        return "generic_identifier"
    if probable_level == "multi_level_or_unclear":
        return "unclear_identifier_like"
    if generic_tokens == ["barcode"]:
        return "level_barcode"
    if generic_tokens == ["uuid"]:
        return "level_uuid"
    if any(token in normalized_name for token in ("number", "sequence", "type", "type_id")):
        return "level_helper"
    return "level_context"


def classify_identifier_linkage_family(probable_level: str) -> str:
    if probable_level in LEVEL_TOKENS:
        return f"{probable_level}_identifier_like"
    return "multi_level_or_unclear_identifier_like"


def table_category(table_name: str) -> str:
    return "central" if table_name in CENTRAL_TABLES else "side"


def first_non_missing_examples(values: list[str], limit: int = 5) -> str:
    examples: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value in MISSING_LIKE_VALUES or value in seen:
            continue
        seen.add(value)
        examples.append(value)
        if len(examples) >= limit:
            break
    return json.dumps(examples, ensure_ascii=True)


def preliminary_manual_review_priority(
    *,
    probable_level: str,
    candidate_source: str,
    field_subtype: str,
    table_name: str,
) -> str:
    if probable_level == "multi_level_or_unclear" or candidate_source in {
        "field_name_identifier_token_match_only",
        "field_name_multi_level_token_match",
    }:
        return "high"
    if field_subtype in FIELD_SUBTYPE_HELPER or table_category(table_name) == "side":
        return "medium"
    return "low"


def build_inventory_rows(
    *,
    crosswalk_run_id: str,
    workflow_inputs: WorkflowInputs,
) -> list[dict[str, Any]]:
    parse_run_id = str(workflow_inputs.parse_latest_pointer["parse_run_id"])
    source_run_id = str(workflow_inputs.parse_latest_pointer["source_run_id"])

    inventory_rows: list[dict[str, Any]] = []
    for table_name in sorted(workflow_inputs.parsed_tables, key=lambda name: TABLE_ORDER.get(name, 99)):
        table_context = workflow_inputs.parsed_tables[table_name]
        for source_position, field_name in enumerate(table_context.header, start=1):
            matched_tokens = find_matched_identifier_tokens(field_name)
            if not matched_tokens:
                continue

            probable_level = classify_probable_identifier_level(matched_tokens)
            candidate_source = classify_candidate_source(matched_tokens)
            identifier_token_family = classify_identifier_token_family(matched_tokens)
            field_subtype = classify_field_subtype(field_name, matched_tokens, probable_level)

            values = [str(row.get(field_name) or "").strip() for row in table_context.data_rows]
            non_missing_values = [value for value in values if value not in MISSING_LIKE_VALUES]
            non_missing_count = len(non_missing_values)
            distinct_non_missing_count = len(set(non_missing_values))
            row_count = len(values)
            missing_like_fraction = round(
                float((row_count - non_missing_count) / row_count) if row_count else 0.0,
                6,
            )
            schema_row = table_context.schema_rows_by_position[source_position]

            inventory_rows.append(
                {
                    "crosswalk_run_id": crosswalk_run_id,
                    "parse_run_id": parse_run_id,
                    "source_run_id": source_run_id,
                    "table_name": table_name,
                    "field_name": field_name,
                    "source_position": source_position,
                    "probable_identifier_level": probable_level,
                    "identifier_token_family": identifier_token_family,
                    "missing_like_fraction": missing_like_fraction,
                    "non_missing_count": non_missing_count,
                    "distinct_non_missing_count": distinct_non_missing_count,
                    "example_values_small_sample": first_non_missing_examples(non_missing_values),
                    "candidate_source": candidate_source,
                    "manual_review_priority": preliminary_manual_review_priority(
                        probable_level=probable_level,
                        candidate_source=candidate_source,
                        field_subtype=field_subtype,
                        table_name=table_name,
                    ),
                    "notes_placeholder": NOTES_PLACEHOLDER,
                    "matched_identifier_tokens_json": json.dumps(matched_tokens, ensure_ascii=True),
                    "alternate_column_name": str(schema_row.get("alternate_column_name") or ""),
                    "cde_id_raw": str(schema_row.get("cde_id_raw") or ""),
                    "_field_subtype": field_subtype,
                    "_table_category": table_category(table_name),
                    "_matched_tokens": matched_tokens,
                    "_identifier_linkage_family": classify_identifier_linkage_family(probable_level),
                }
            )

    return inventory_rows


def build_rows_by_table(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["table_name"])].append(row)
    return grouped


def pattern_strength_rank(pattern_rule: str) -> int:
    order = {
        "exact_pair_one_to_one_same_table": 0,
        "barcode_prefix_same_row": 1,
        "same_row_copresence_only": 2,
        "not_applicable_no_pair_rule": 99,
    }
    return order.get(pattern_rule, 99)


def compute_exact_pair_evidence(
    *,
    table_context: ParsedTableContext,
    field_name: str,
    paired_field_name: str,
) -> dict[str, Any]:
    comparable_pairs: list[tuple[str, str]] = []
    left_to_right: dict[str, set[str]] = defaultdict(set)
    right_to_left: dict[str, set[str]] = defaultdict(set)

    for row in table_context.data_rows:
        left_value = str(row.get(field_name) or "").strip()
        right_value = str(row.get(paired_field_name) or "").strip()
        if left_value in MISSING_LIKE_VALUES or right_value in MISSING_LIKE_VALUES:
            continue
        comparable_pairs.append((left_value, right_value))
        left_to_right[left_value].add(right_value)
        right_to_left[right_value].add(left_value)

    pass_count = 0
    nonmatch_examples: list[dict[str, str]] = []
    for left_value, right_value in comparable_pairs:
        passes = len(left_to_right[left_value]) == 1 and len(right_to_left[right_value]) == 1
        if passes:
            pass_count += 1
        elif len(nonmatch_examples) < 5:
            nonmatch_examples.append(
                {
                    "field_value": left_value,
                    "paired_field_value": right_value,
                }
            )

    comparable_row_count = len(comparable_pairs)
    return {
        "paired_field_name": paired_field_name,
        "pattern_rule": "exact_pair_one_to_one_same_table",
        "pattern_comparable_row_count": comparable_row_count,
        "pattern_pass_fraction": round(float(pass_count / comparable_row_count), 6) if comparable_row_count else 0.0,
        "pattern_nonmatch_examples_json": json.dumps(nonmatch_examples, ensure_ascii=True),
    }


def compute_prefix_evidence(
    *,
    table_context: ParsedTableContext,
    prefix_field_name: str,
    child_field_name: str,
    paired_field_name: str,
) -> dict[str, Any]:
    comparable_pairs: list[tuple[str, str]] = []
    nonmatch_examples: list[dict[str, str]] = []
    pass_count = 0

    for row in table_context.data_rows:
        prefix_value = str(row.get(prefix_field_name) or "").strip()
        child_value = str(row.get(child_field_name) or "").strip()
        if prefix_value in MISSING_LIKE_VALUES or child_value in MISSING_LIKE_VALUES:
            continue
        comparable_pairs.append((prefix_value, child_value))
        if child_value.startswith(prefix_value):
            pass_count += 1
        elif len(nonmatch_examples) < 5:
            nonmatch_examples.append(
                {
                    "prefix_value": prefix_value,
                    "child_value": child_value,
                }
            )

    comparable_row_count = len(comparable_pairs)
    return {
        "paired_field_name": paired_field_name,
        "pattern_rule": "barcode_prefix_same_row",
        "pattern_comparable_row_count": comparable_row_count,
        "pattern_pass_fraction": round(float(pass_count / comparable_row_count), 6) if comparable_row_count else 0.0,
        "pattern_nonmatch_examples_json": json.dumps(nonmatch_examples, ensure_ascii=True),
    }


def compute_copresence_evidence(
    *,
    table_context: ParsedTableContext,
    field_name: str,
    paired_field_name: str,
) -> dict[str, Any]:
    field_non_missing_count = 0
    paired_non_missing_count = 0
    nonmatch_examples: list[dict[str, str]] = []

    for row in table_context.data_rows:
        field_value = str(row.get(field_name) or "").strip()
        if field_value in MISSING_LIKE_VALUES:
            continue
        field_non_missing_count += 1
        paired_value = str(row.get(paired_field_name) or "").strip()
        if paired_value in MISSING_LIKE_VALUES:
            if len(nonmatch_examples) < 5:
                nonmatch_examples.append({"field_value": field_value})
            continue
        paired_non_missing_count += 1

    return {
        "paired_field_name": paired_field_name,
        "pattern_rule": "same_row_copresence_only",
        "pattern_comparable_row_count": field_non_missing_count,
        "pattern_pass_fraction": round(float(paired_non_missing_count / field_non_missing_count), 6)
        if field_non_missing_count
        else 0.0,
        "pattern_nonmatch_examples_json": json.dumps(nonmatch_examples, ensure_ascii=True),
    }


def find_exact_pair_field(table_rows: list[dict[str, Any]], current_row: dict[str, Any]) -> str | None:
    probable_level = str(current_row["probable_identifier_level"])
    field_subtype = str(current_row["_field_subtype"])
    if field_subtype not in FIELD_SUBTYPE_STRONG or probable_level == "multi_level_or_unclear":
        return None

    target_subtype = "level_uuid" if field_subtype == "level_barcode" else "level_barcode"
    for peer_row in table_rows:
        if peer_row is current_row:
            continue
        if str(peer_row["probable_identifier_level"]) != probable_level:
            continue
        if str(peer_row["_field_subtype"]) != target_subtype:
            continue
        return str(peer_row["field_name"])
    return None


def find_prefix_partner(table_rows: list[dict[str, Any]], current_row: dict[str, Any]) -> tuple[str, bool] | None:
    probable_level = str(current_row["probable_identifier_level"])
    field_subtype = str(current_row["_field_subtype"])
    if field_subtype != "level_barcode":
        return None

    strong_barcode_fields = [
        row
        for row in table_rows
        if str(row["_field_subtype"]) == "level_barcode" and str(row["field_name"]) != str(current_row["field_name"])
    ]
    strong_barcode_by_level = {
        str(row["probable_identifier_level"]): str(row["field_name"]) for row in strong_barcode_fields
    }

    if probable_level == "patient" and "sample" in strong_barcode_by_level:
        return strong_barcode_by_level["sample"], True
    if probable_level == "sample":
        if "patient" in strong_barcode_by_level:
            return strong_barcode_by_level["patient"], False
        for child_level in ("portion", "analyte", "slide", "aliquot"):
            if child_level in strong_barcode_by_level:
                return strong_barcode_by_level[child_level], True
        return None
    if probable_level in {"portion", "analyte", "slide", "aliquot"} and "sample" in strong_barcode_by_level:
        return strong_barcode_by_level["sample"], False
    return None


def find_copresence_partner(table_rows: list[dict[str, Any]], current_row: dict[str, Any]) -> str | None:
    probable_level = str(current_row["probable_identifier_level"])
    current_field_name = str(current_row["field_name"])

    same_level_strong_rows = [
        row
        for row in table_rows
        if str(row["probable_identifier_level"]) == probable_level
        and str(row["_field_subtype"]) in FIELD_SUBTYPE_STRONG
        and str(row["field_name"]) != current_field_name
    ]
    if same_level_strong_rows:
        return str(sorted(same_level_strong_rows, key=lambda row: int(row["source_position"]))[0]["field_name"])

    preferred_levels = []
    if probable_level == "patient":
        preferred_levels = ["sample"]
    elif probable_level == "sample":
        preferred_levels = ["patient", "portion", "analyte", "slide", "aliquot"]
    elif probable_level in {"portion", "analyte", "slide", "aliquot"}:
        preferred_levels = ["sample", "patient"]
    elif probable_level == "multi_level_or_unclear":
        preferred_levels = ["sample", "portion", "analyte", "slide", "aliquot", "patient"]

    for preferred_level in preferred_levels:
        preferred_rows = [
            row
            for row in table_rows
            if str(row["probable_identifier_level"]) == preferred_level and str(row["_field_subtype"]) in FIELD_SUBTYPE_STRONG
        ]
        if preferred_rows:
            return str(sorted(preferred_rows, key=lambda row: int(row["source_position"]))[0]["field_name"])
    return None


def add_pattern_evidence(
    inventory_rows: list[dict[str, Any]],
    parsed_tables: dict[str, ParsedTableContext],
) -> None:
    rows_by_table = build_rows_by_table(inventory_rows)
    for row in inventory_rows:
        table_name = str(row["table_name"])
        table_rows = rows_by_table[table_name]
        table_context = parsed_tables[table_name]
        field_name = str(row["field_name"])

        exact_pair_field = find_exact_pair_field(table_rows, row)
        if exact_pair_field is not None:
            evidence = compute_exact_pair_evidence(
                table_context=table_context,
                field_name=field_name,
                paired_field_name=exact_pair_field,
            )
        else:
            prefix_partner = find_prefix_partner(table_rows, row)
            if prefix_partner is not None:
                paired_field_name, current_is_prefix = prefix_partner
                if current_is_prefix:
                    evidence = compute_prefix_evidence(
                        table_context=table_context,
                        prefix_field_name=field_name,
                        child_field_name=paired_field_name,
                        paired_field_name=paired_field_name,
                    )
                else:
                    evidence = compute_prefix_evidence(
                        table_context=table_context,
                        prefix_field_name=paired_field_name,
                        child_field_name=field_name,
                        paired_field_name=paired_field_name,
                    )
            else:
                copresence_partner = find_copresence_partner(table_rows, row)
                if copresence_partner is not None:
                    evidence = compute_copresence_evidence(
                        table_context=table_context,
                        field_name=field_name,
                        paired_field_name=copresence_partner,
                    )
                else:
                    evidence = {
                        "paired_field_name": "",
                        "pattern_rule": "not_applicable_no_pair_rule",
                        "pattern_comparable_row_count": 0,
                        "pattern_pass_fraction": 0.0,
                        "pattern_nonmatch_examples_json": json.dumps([], ensure_ascii=True),
                    }

        row.update(evidence)


def has_strong_identifier_evidence(row: dict[str, Any]) -> bool:
    if str(row["_field_subtype"]) not in FIELD_SUBTYPE_STRONG:
        return False
    comparable_row_count = parse_int(row["pattern_comparable_row_count"], "pattern_comparable_row_count")
    if comparable_row_count <= 0:
        return False
    return parse_float(row["pattern_pass_fraction"], "pattern_pass_fraction") >= 0.95


def patient_primary_rank_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        -int(row["distinct_non_missing_count"]),
        float(row["missing_like_fraction"]),
        0 if str(row["_table_category"]) == "central" else 1,
        PATIENT_PRIMARY_TABLE_ORDER.get(str(row["table_name"]), 99),
        pattern_strength_rank(str(row["pattern_rule"])),
        -int(row["non_missing_count"]),
        int(row["source_position"]),
    )


def crosswalk_sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        FAMILY_ORDER.index(str(row["identifier_linkage_family"])),
        ROLE_ORDER.index(str(row["crosswalk_role"])),
        TABLE_ORDER.get(str(row["table_name"]), 99),
        float(row["missing_like_fraction"]),
        -int(row["non_missing_count"]),
        int(row["source_position"]),
    )


def crosswalk_manual_review_priority(row: dict[str, Any]) -> str:
    crosswalk_role = str(row["crosswalk_role"])
    if crosswalk_role == "ambiguous_candidate":
        return "high"
    if crosswalk_role in {"overlapping_link_candidate", "secondary_link_candidate", "side_table_candidate"}:
        return "medium"
    if str(row["_table_category"]) == "side":
        return "medium"
    return "low"


def assign_crosswalk_rows(inventory_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows_by_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in inventory_rows:
        rows_by_family[str(row["_identifier_linkage_family"])].append(row)

    crosswalk_rows: list[dict[str, Any]] = []
    for family_name in FAMILY_ORDER:
        family_rows = rows_by_family.get(family_name, [])
        if not family_rows:
            continue

        patient_primary_keys: set[tuple[str, str]] = set()
        if family_name == "patient_identifier_like":
            barcode_candidates = [
                row
                for row in family_rows
                if str(row["_field_subtype"]) == "level_barcode" and has_strong_identifier_evidence(row)
            ]
            uuid_candidates = [
                row
                for row in family_rows
                if str(row["_field_subtype"]) == "level_uuid" and has_strong_identifier_evidence(row)
            ]
            if barcode_candidates:
                top_barcode = sorted(barcode_candidates, key=patient_primary_rank_key)[0]
                patient_primary_keys.add((str(top_barcode["table_name"]), str(top_barcode["field_name"])))
            if uuid_candidates:
                top_uuid = sorted(uuid_candidates, key=patient_primary_rank_key)[0]
                patient_primary_keys.add((str(top_uuid["table_name"]), str(top_uuid["field_name"])))

        for row in family_rows:
            probable_level = str(row["probable_identifier_level"])
            field_subtype = str(row["_field_subtype"])
            table_name = str(row["table_name"])
            table_kind = str(row["_table_category"])
            row_key = (table_name, str(row["field_name"]))

            if family_name == "multi_level_or_unclear_identifier_like":
                crosswalk_role = "ambiguous_candidate"
                crosswalk_rule = "role_ambiguous_generic_or_conflicted_identifier"
            elif family_name == "patient_identifier_like" and row_key in patient_primary_keys:
                if field_subtype == "level_barcode":
                    crosswalk_role = "primary_link_candidate"
                    crosswalk_rule = "role_primary_patient_best_barcode_candidate"
                else:
                    crosswalk_role = "primary_link_candidate"
                    crosswalk_rule = "role_primary_patient_best_uuid_candidate"
            elif field_subtype in FIELD_SUBTYPE_STRONG and not has_strong_identifier_evidence(row):
                crosswalk_role = "ambiguous_candidate"
                crosswalk_rule = "role_ambiguous_generic_or_conflicted_identifier"
            elif probable_level in CANONICAL_TABLE_BY_LEVEL and table_name == CANONICAL_TABLE_BY_LEVEL[probable_level]:
                if field_subtype in FIELD_SUBTYPE_STRONG:
                    crosswalk_role = "primary_link_candidate"
                    crosswalk_rule = "role_primary_canonical_level_identifier"
                elif field_subtype in FIELD_SUBTYPE_HELPER:
                    crosswalk_role = "secondary_link_candidate"
                    crosswalk_rule = "role_secondary_level_helper_field"
                else:
                    crosswalk_role = "ambiguous_candidate"
                    crosswalk_rule = "role_ambiguous_generic_or_conflicted_identifier"
            elif field_subtype in FIELD_SUBTYPE_STRONG and table_kind == "side":
                crosswalk_role = "side_table_candidate"
                crosswalk_rule = "role_side_table_strong_identifier"
            elif field_subtype in FIELD_SUBTYPE_STRONG:
                crosswalk_role = "overlapping_link_candidate"
                crosswalk_rule = "role_overlapping_repeated_strong_identifier"
            elif field_subtype in FIELD_SUBTYPE_HELPER:
                crosswalk_role = "secondary_link_candidate"
                crosswalk_rule = "role_secondary_level_helper_field"
            else:
                crosswalk_role = "ambiguous_candidate"
                crosswalk_rule = "role_ambiguous_generic_or_conflicted_identifier"

            crosswalk_row = {
                "crosswalk_run_id": row["crosswalk_run_id"],
                "parse_run_id": row["parse_run_id"],
                "source_run_id": row["source_run_id"],
                "identifier_linkage_family": row["_identifier_linkage_family"],
                "table_name": row["table_name"],
                "field_name": row["field_name"],
                "source_position": row["source_position"],
                "probable_identifier_level": row["probable_identifier_level"],
                "identifier_token_family": row["identifier_token_family"],
                "candidate_source": row["candidate_source"],
                "missing_like_fraction": row["missing_like_fraction"],
                "non_missing_count": row["non_missing_count"],
                "distinct_non_missing_count": row["distinct_non_missing_count"],
                "example_values_small_sample": row["example_values_small_sample"],
                "crosswalk_role": crosswalk_role,
                "crosswalk_rule": crosswalk_rule,
                "manual_review_priority": "pending",
                "notes_placeholder": row["notes_placeholder"],
                "matched_identifier_tokens_json": row["matched_identifier_tokens_json"],
                "alternate_column_name": row["alternate_column_name"],
                "cde_id_raw": row["cde_id_raw"],
                "paired_field_name": row["paired_field_name"],
                "pattern_rule": row["pattern_rule"],
                "pattern_comparable_row_count": row["pattern_comparable_row_count"],
                "pattern_pass_fraction": row["pattern_pass_fraction"],
                "pattern_nonmatch_examples_json": row["pattern_nonmatch_examples_json"],
                "_field_subtype": row["_field_subtype"],
                "_table_category": row["_table_category"],
                "_matched_tokens": row["_matched_tokens"],
            }
            crosswalk_row["manual_review_priority"] = crosswalk_manual_review_priority(crosswalk_row)
            crosswalk_rows.append(crosswalk_row)

    return sorted(crosswalk_rows, key=crosswalk_sort_key)


def build_summary_rows(crosswalk_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows_by_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in crosswalk_rows:
        rows_by_family[str(row["identifier_linkage_family"])].append(row)

    summary_rows: list[dict[str, Any]] = []
    for family_name in FAMILY_ORDER:
        family_rows = rows_by_family.get(family_name, [])
        if not family_rows:
            continue

        def rows_to_json(role_name: str) -> str:
            values = sorted(
                f"{row['table_name']}.{row['field_name']}"
                for row in family_rows
                if str(row["crosswalk_role"]) == role_name
            )
            return json.dumps(values, ensure_ascii=True)

        summary_rows.append(
            {
                "crosswalk_run_id": str(family_rows[0]["crosswalk_run_id"]),
                "parse_run_id": str(family_rows[0]["parse_run_id"]),
                "source_run_id": str(family_rows[0]["source_run_id"]),
                "identifier_linkage_family": family_name,
                "field_count": len(family_rows),
                "tables_present_json": json.dumps(
                    sorted({str(row["table_name"]) for row in family_rows}),
                    ensure_ascii=True,
                ),
                "primary_candidates_json": rows_to_json("primary_link_candidate"),
                "secondary_candidates_json": rows_to_json("secondary_link_candidate"),
                "overlapping_candidates_json": rows_to_json("overlapping_link_candidate"),
                "ambiguous_candidates_json": rows_to_json("ambiguous_candidate"),
                "side_table_candidates_json": rows_to_json("side_table_candidate"),
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )

    return summary_rows


def sanitize_rows_for_tsv(rows: list[dict[str, Any]], fieldnames: list[str]) -> list[dict[str, Any]]:
    return [{field_name: row[field_name] for field_name in fieldnames} for row in rows]


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
            raise BiospecimenIdentifierCrosswalkError(
                "No identifier-like fields were found in the parsed biospecimen tables."
            )

        add_pattern_evidence(inventory_rows, workflow_inputs.parsed_tables)
        inventory_path = crosswalk_run_dir / "biospecimen_identifier_inventory.tsv"
        write_dict_rows_tsv(
            inventory_path,
            INVENTORY_FIELDNAMES,
            sanitize_rows_for_tsv(inventory_rows, INVENTORY_FIELDNAMES),
        )

        crosswalk_rows = assign_crosswalk_rows(inventory_rows)
        if not crosswalk_rows:
            raise BiospecimenIdentifierCrosswalkError("Biospecimen identifier crosswalk rows were not generated.")
        crosswalk_path = crosswalk_run_dir / "biospecimen_identifier_crosswalk.tsv"
        write_dict_rows_tsv(
            crosswalk_path,
            CROSSWALK_FIELDNAMES,
            sanitize_rows_for_tsv(crosswalk_rows, CROSSWALK_FIELDNAMES),
        )

        summary_rows = build_summary_rows(crosswalk_rows)
        if not summary_rows:
            raise BiospecimenIdentifierCrosswalkError(
                "Biospecimen identifier crosswalk summary rows were not generated."
            )
        summary_path = crosswalk_run_dir / "biospecimen_identifier_crosswalk_summary.tsv"
        write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        family_values = {str(row["identifier_linkage_family"]) for row in crosswalk_rows}
        summary_family_values = {str(row["identifier_linkage_family"]) for row in summary_rows}
        summary_total = sum(parse_int(row["field_count"], "summary field_count") for row in summary_rows)
        inventory_candidate_sources_valid = all(
            str(row["candidate_source"]) in CANDIDATE_SOURCES for row in inventory_rows
        )
        crosswalk_roles_valid = all(
            str(row["crosswalk_role"]) in ROLE_ORDER
            and str(row["crosswalk_rule"]) in CROSSWALK_RULES
            and str(row["pattern_rule"]) in PATTERN_RULES
            for row in crosswalk_rows
        )
        output_rows_positive = all(len(rows) > 0 for rows in (inventory_rows, crosswalk_rows, summary_rows))
        summary_matches_crosswalk = family_values == summary_family_values and summary_total == len(crosswalk_rows)

        parse_run_id = str(workflow_inputs.parse_latest_pointer["parse_run_id"])
        source_run_id = str(workflow_inputs.parse_latest_pointer["source_run_id"])
        latest_pointer_payload = {
            "updated_at_utc": format_utc_timestamp(utc_now()),
            "crosswalk_run_id": crosswalk_run_id,
            "parse_run_id": parse_run_id,
            "source_run_id": source_run_id,
            "crosswalk_run_directory": repo_relative(crosswalk_run_dir, paths.repo_root),
            "biospecimen_identifier_inventory_tsv": repo_relative(inventory_path, paths.repo_root),
            "biospecimen_identifier_crosswalk_tsv": repo_relative(crosswalk_path, paths.repo_root),
            "biospecimen_identifier_crosswalk_summary_tsv": repo_relative(summary_path, paths.repo_root),
            "run_log_json": repo_relative(run_log_path, paths.repo_root),
            "biospecimen_biotab_latest_json": repo_relative(paths.parse_latest_pointer, paths.repo_root),
            "field_count": len(inventory_rows),
        }

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "crosswalk_run_id": crosswalk_run_id,
            "parse_run_id": parse_run_id,
            "source_run_id": source_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "biospecimen_biotab_latest_json": repo_relative(paths.parse_latest_pointer, paths.repo_root),
                "biospecimen_source_inventory_tsv": repo_relative(
                    workflow_inputs.input_paths["source_inventory_tsv"],
                    paths.repo_root,
                ),
                "biospecimen_manifest_tsv": repo_relative(
                    workflow_inputs.input_paths["manifest_tsv"],
                    paths.repo_root,
                ),
                "biospecimen_parse_run_log_json": repo_relative(
                    workflow_inputs.input_paths["parse_run_log_json"],
                    paths.repo_root,
                ),
                "parsed_tables": {
                    table_name: {
                        "source_filename": table_context.source_filename,
                        "parsed_tsv": repo_relative(table_context.parsed_path, paths.repo_root),
                        "schema_tsv": repo_relative(table_context.schema_path, paths.repo_root),
                        "row_count": table_context.row_count,
                        "column_count": table_context.column_count,
                    }
                    for table_name, table_context in workflow_inputs.parsed_tables.items()
                },
            },
            "outputs": {
                "crosswalk_run_directory": repo_relative(crosswalk_run_dir, paths.repo_root),
                "biospecimen_identifier_inventory_tsv": repo_relative(inventory_path, paths.repo_root),
                "biospecimen_identifier_crosswalk_tsv": repo_relative(crosswalk_path, paths.repo_root),
                "biospecimen_identifier_crosswalk_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    output_rows_positive
                    and inventory_candidate_sources_valid
                    and crosswalk_roles_valid
                    and summary_matches_crosswalk
                ),
                "biospecimen_parse_latest_pointer_found": True,
                "biospecimen_parse_run_log_completed": True,
                "required_central_tables_found": all(
                    table_name in workflow_inputs.parsed_tables for table_name in CENTRAL_TABLES
                ),
                "optional_side_tables_loaded": sorted(
                    table_name
                    for table_name in workflow_inputs.parsed_tables
                    if table_name in OPTIONAL_SIDE_TABLES
                ),
                "identifier_like_field_count_positive": len(inventory_rows) > 0,
                "inventory_row_count_positive": len(inventory_rows) > 0,
                "crosswalk_row_count_positive": len(crosswalk_rows) > 0,
                "summary_row_count_positive": len(summary_rows) > 0,
                "inventory_candidate_sources_valid": inventory_candidate_sources_valid,
                "crosswalk_roles_and_rules_valid": crosswalk_roles_valid,
                "summary_family_counts_match_crosswalk": summary_matches_crosswalk,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "identifier_tokens": list(IDENTIFIER_TOKENS),
                "level_tokens": list(LEVEL_TOKENS),
                "central_tables": list(CENTRAL_TABLES),
                "optional_side_tables": list(OPTIONAL_SIDE_TABLES),
                "canonical_table_by_level": dict(CANONICAL_TABLE_BY_LEVEL),
                "candidate_sources": list(CANDIDATE_SOURCES),
                "family_order": list(FAMILY_ORDER),
                "role_order": list(ROLE_ORDER),
                "crosswalk_rules": list(CROSSWALK_RULES),
                "pattern_rules": [
                    {
                        "rule": "exact_pair_one_to_one_same_table",
                        "match": "same-level barcode/uuid fields in the same table map one-to-one",
                    },
                    {
                        "rule": "barcode_prefix_same_row",
                        "match": "child barcode starts with parent barcode in the same row",
                    },
                    {
                        "rule": "same_row_copresence_only",
                        "match": "field and paired field are both present in the same row without stronger normalization",
                    },
                    {
                        "rule": "not_applicable_no_pair_rule",
                        "match": "no reversible pairing rule was available from the saved parsed table only",
                    },
                ],
                "patient_primary_strategy": {
                    "barcode": "best strong patient barcode candidate by distinct coverage, completeness, table priority, and evidence rule",
                    "uuid": "best strong patient UUID candidate by distinct coverage, completeness, table priority, and evidence rule",
                },
                "notes_placeholder": NOTES_PLACEHOLDER,
            },
            "inventory_counts_by_candidate_source": {
                candidate_source: sum(str(row["candidate_source"]) == candidate_source for row in inventory_rows)
                for candidate_source in CANDIDATE_SOURCES
            },
            "crosswalk_counts_by_role": {
                role_name: sum(str(row["crosswalk_role"]) == role_name for row in crosswalk_rows)
                for role_name in ROLE_ORDER
            },
            "crosswalk_counts_by_family": {
                family_name: sum(str(row["identifier_linkage_family"]) == family_name for row in crosswalk_rows)
                for family_name in sorted(family_values, key=FAMILY_ORDER.index)
            },
            "crosswalk_counts_by_family_and_role": {
                family_name: {
                    role_name: sum(
                        str(row["identifier_linkage_family"]) == family_name and str(row["crosswalk_role"]) == role_name
                        for row in crosswalk_rows
                    )
                    for role_name in ROLE_ORDER
                }
                for family_name in sorted(family_values, key=FAMILY_ORDER.index)
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "biospecimen_biotab_latest_pointer": workflow_inputs.parse_latest_pointer,
                "biospecimen_biotab_run_log": workflow_inputs.parse_run_log,
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
            "workflow": "tcga_brca_biospecimen_identifier_crosswalk",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA biospecimen identifier crosswalk audit complete.")
    print(f"Crosswalk run ID: {run_log['crosswalk_run_id']}")
    print(f"Parse run ID: {run_log['parse_run_id']}")
    print(f"Source run ID: {run_log['source_run_id']}")
    print(f"Crosswalk output directory: {run_log['outputs']['crosswalk_run_directory']}")
    print(f"Inventory TSV: {run_log['outputs']['biospecimen_identifier_inventory_tsv']}")
    print(f"Crosswalk TSV: {run_log['outputs']['biospecimen_identifier_crosswalk_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['biospecimen_identifier_crosswalk_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
