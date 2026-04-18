#!/usr/bin/env python
"""Audit the 4 core TCGA-BRCA clinical biotab tables at the field level."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover
    yaml = None

CORE_TABLE_NAMES = (
    "clinical_patient",
    "clinical_drug",
    "clinical_radiation",
    "clinical_follow_up_v4_0",
)
FIELD_GROUP_LABELS = (
    "demographics",
    "diagnosis / pathology",
    "stage",
    "receptor / biomarker",
    "treatment / drug",
    "radiation",
    "follow-up / outcome-like",
    "identifier / admin",
    "other / unclear",
)
FIELD_AUDIT_FIELDNAMES = [
    "audit_run_id",
    "parse_run_id",
    "source_run_id",
    "table_name",
    "field_name",
    "source_position",
    "alternate_column_name",
    "cde_id_raw",
    "row_count",
    "non_missing_count",
    "missing_like_count",
    "missing_like_fraction",
    "distinct_non_missing_count",
    "example_values_small_sample",
    "probable_field_group",
    "group_assignment_rule",
    "blank_count",
    "not_available_count",
    "not_applicable_count",
    "unknown_count",
    "not_evaluated_count",
    "discrepancy_count",
    "na_count",
    "n_a_count",
    "null_count",
    "none_count",
    "nan_count",
]
GROUP_SUMMARY_FIELDNAMES = [
    "audit_run_id",
    "parse_run_id",
    "source_run_id",
    "table_name",
    "probable_field_group",
    "field_count",
    "field_fraction_of_table",
    "non_missing_any_field_count",
    "high_missingness_field_count",
    "field_names_json",
]
MISSINGNESS_SUMMARY_FIELDNAMES = [
    "audit_run_id",
    "parse_run_id",
    "source_run_id",
    "table_name",
    "field_name",
    "source_position",
    "probable_field_group",
    "row_count",
    "non_missing_count",
    "missing_like_count",
    "missing_like_fraction",
    "blank_count",
    "not_available_count",
    "not_applicable_count",
    "unknown_count",
    "not_evaluated_count",
    "discrepancy_count",
    "na_count",
    "n_a_count",
    "null_count",
    "none_count",
    "nan_count",
    "distinct_non_missing_count",
    "example_values_small_sample",
]
MISSING_TOKEN_COLUMNS = {
    "blank": "blank_count",
    "[not available]": "not_available_count",
    "[not applicable]": "not_applicable_count",
    "[unknown]": "unknown_count",
    "[not evaluated]": "not_evaluated_count",
    "[discrepancy]": "discrepancy_count",
    "na": "na_count",
    "n/a": "n_a_count",
    "null": "null_count",
    "none": "none_count",
    "nan": "nan_count",
}
IDENTIFIER_ADMIN_EXACT_NAMES = {
    "disease_code",
    "form_completion_date",
    "informed_consent_verified",
    "patient_id",
    "project_code",
    "prospective_collection",
    "retrospective_collection",
    "tissue_source_site",
}
IDENTIFIER_ADMIN_PREFIXES = ("bcr_",)
GROUP_KEYWORD_RULES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("demographics", ("gender", "race", "ethnicity", "menopause", "birth", "age_at_")),
    ("stage", ("ajcc", "clinical_t", "clinical_n", "clinical_m", "clinical_stage", "pathologic_", "stage_")),
    ("receptor / biomarker", ("er_", "pr_", "her2", "receptor", "ihc", "fish", "cent17", "cent_17")),
    ("radiation", ("radiation", "irradiat")),
    (
        "follow-up / outcome-like",
        (
            "followup",
            "follow_up",
            "vital_status",
            "death",
            "last_contact",
            "tumor_status",
            "progression",
            "new_tumor_event",
            "response",
            "lost_to",
        ),
    ),
    (
        "treatment / drug",
        (
            "drug",
            "pharm",
            "therapy",
            "treatment",
            "regimen",
            "dose",
            "route_of_administration",
            "stem_cell",
            "tx_",
            "tx_on_",
            "cycles",
        ),
    ),
    (
        "diagnosis / pathology",
        (
            "diagnosis",
            "histolog",
            "margin",
            "lymph",
            "metastasis_site",
            "axillary",
            "surgical",
            "anatomic",
            "neoplasm",
            "tumor_tissue_site",
            "site_of_primary_tumor",
            "icd_",
            "method_initial_path_dx",
            "metastatic_tumor",
            "subdivision",
        ),
    ),
)
DEFAULT_GROUP_BY_TABLE = {
    "clinical_drug": "treatment / drug",
    "clinical_radiation": "radiation",
    "clinical_follow_up_v4_0": "follow-up / outcome-like",
    "clinical_patient": "other / unclear",
}


class ClinicalCoreFieldAuditError(RuntimeError):
    """Raised when the clinical core field audit cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the clinical core field audit workflow."""

    repo_root: Path
    trial_config: Path
    clinical_biotab_latest_pointer: Path
    audit_runs_root: Path
    latest_pointer: Path
    results_root: Path


@dataclass(frozen=True)
class CoreTableInput:
    """One core parsed clinical table and its upstream context."""

    parse_run_id: str
    source_run_id: str
    table_name: str
    source_filename: str
    parsed_path: Path
    schema_path: Path
    row_count: int
    column_count: int


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise ClinicalCoreFieldAuditError("Unable to locate the repository root from the script path.")


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
                raise ClinicalCoreFieldAuditError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise ClinicalCoreFieldAuditError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise ClinicalCoreFieldAuditError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise ClinicalCoreFieldAuditError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise ClinicalCoreFieldAuditError(
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
        raise ClinicalCoreFieldAuditError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ClinicalCoreFieldAuditError(f"Expected a JSON object in file: {path}")
    return data


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise ClinicalCoreFieldAuditError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise ClinicalCoreFieldAuditError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ClinicalCoreFieldAuditError(f"Run directory already exists: {path}")
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
        clinical_biotab_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_biotabs_latest.json"
        ),
        audit_runs_root=audit_root / "tcga-brca" / "variables" / "clinical_core_field_audit_runs",
        latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_core_field_audit_latest.json",
        results_root=results_root,
    )


def read_tsv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def read_table_rows(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)

    if not rows:
        raise ClinicalCoreFieldAuditError(f"Parsed table is empty: {path}")

    header = rows[0]
    data_rows = rows[1:]
    if not header:
        raise ClinicalCoreFieldAuditError(f"Parsed table has an empty header row: {path}")
    if not data_rows:
        raise ClinicalCoreFieldAuditError(f"Parsed table has zero data rows: {path}")

    column_count = len(header)
    for line_number, row in enumerate(data_rows, start=2):
        if len(row) != column_count:
            raise ClinicalCoreFieldAuditError(
                f"Parsed table row {line_number} has {len(row)} columns, expected {column_count}: {path}"
            )

    return header, data_rows


def load_manifest_rows(manifest_path: Path) -> list[dict[str, str]]:
    manifest_rows = read_tsv_dict_rows(manifest_path)
    if not manifest_rows:
        raise ClinicalCoreFieldAuditError(f"Clinical biotab manifest has no rows: {manifest_path}")
    required_columns = {
        "parse_run_id",
        "source_run_id",
        "table_name",
        "source_filename",
        "output_path",
        "schema_path",
        "row_count",
        "column_count",
    }
    missing_columns = required_columns.difference(manifest_rows[0].keys())
    if missing_columns:
        raise ClinicalCoreFieldAuditError(
            f"Clinical biotab manifest is missing required columns {sorted(missing_columns)}: {manifest_path}"
        )
    return manifest_rows


def load_core_table_inputs(
    *,
    paths: WorkflowPaths,
    parse_pointer: dict[str, Any],
    parse_run_log: dict[str, Any],
) -> tuple[Path, list[CoreTableInput]]:
    manifest_path = paths.repo_root / str(parse_pointer.get("table_manifest_tsv") or "")
    if not manifest_path.exists():
        raise ClinicalCoreFieldAuditError(f"Clinical biotab manifest does not exist: {manifest_path}")

    manifest_rows = load_manifest_rows(manifest_path)
    rows_by_table: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in manifest_rows:
        rows_by_table[row["table_name"]].append(row)

    core_tables: list[CoreTableInput] = []
    parse_run_id = str(parse_pointer.get("parse_run_id") or "").strip()
    source_run_id = str(parse_pointer.get("source_run_id") or "").strip()
    if not parse_run_id or not source_run_id:
        raise ClinicalCoreFieldAuditError(
            f"Clinical biotab latest pointer is missing parse_run_id or source_run_id: {paths.clinical_biotab_latest_pointer}"
        )

    run_log_parse_run_id = str(parse_run_log.get("parse_run_id") or "").strip()
    if run_log_parse_run_id and run_log_parse_run_id != parse_run_id:
        raise ClinicalCoreFieldAuditError(
            f"Clinical biotab parse run mismatch between latest pointer and run log: {parse_run_id} vs {run_log_parse_run_id}"
        )

    for table_name in CORE_TABLE_NAMES:
        matching_rows = rows_by_table.get(table_name, [])
        if not matching_rows:
            raise ClinicalCoreFieldAuditError(f"Core clinical table is missing from the manifest: {table_name}")
        if len(matching_rows) != 1:
            raise ClinicalCoreFieldAuditError(f"Expected exactly one manifest row for {table_name}, found {len(matching_rows)}.")

        row = matching_rows[0]
        parsed_path = paths.repo_root / row["output_path"]
        schema_path = paths.repo_root / row["schema_path"]
        if not parsed_path.exists():
            raise ClinicalCoreFieldAuditError(f"Parsed clinical table does not exist on disk: {parsed_path}")
        if not schema_path.exists():
            raise ClinicalCoreFieldAuditError(f"Schema TSV does not exist on disk: {schema_path}")

        try:
            row_count = int(row["row_count"])
            column_count = int(row["column_count"])
        except ValueError as exc:
            raise ClinicalCoreFieldAuditError(
                f"Manifest row_count or column_count is not an integer for {table_name}: {row}"
            ) from exc

        if row_count <= 0:
            raise ClinicalCoreFieldAuditError(f"Core clinical table has zero rows in the manifest: {table_name}")
        if column_count <= 0:
            raise ClinicalCoreFieldAuditError(f"Core clinical table has zero columns in the manifest: {table_name}")

        core_tables.append(
            CoreTableInput(
                parse_run_id=parse_run_id,
                source_run_id=source_run_id,
                table_name=table_name,
                source_filename=row["source_filename"],
                parsed_path=parsed_path,
                schema_path=schema_path,
                row_count=row_count,
                column_count=column_count,
            )
        )

    return manifest_path, core_tables


def normalize_value(value: str) -> str:
    return value.strip().lower()


def classify_missing_like_token(value: str) -> str | None:
    normalized = normalize_value(value)
    if normalized == "":
        return "blank"
    if normalized in MISSING_TOKEN_COLUMNS:
        return normalized
    return None


def select_example_values(counter: Counter[str], limit: int = 5) -> str:
    ordered_values = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    example_values = [value for value, _count in ordered_values[:limit]]
    return json.dumps(example_values, ensure_ascii=True)


def assign_probable_field_group(table_name: str, field_name: str) -> tuple[str, str]:
    normalized = field_name.lower()

    if normalized in IDENTIFIER_ADMIN_EXACT_NAMES:
        return "identifier / admin", f"identifier_exact_name:{normalized}"

    for prefix in IDENTIFIER_ADMIN_PREFIXES:
        if normalized.startswith(prefix):
            return "identifier / admin", f"identifier_prefix:{prefix}"

    for group_name, patterns in GROUP_KEYWORD_RULES:
        for pattern in patterns:
            if pattern in normalized:
                return group_name, f"keyword:{pattern}"

    default_group = DEFAULT_GROUP_BY_TABLE.get(table_name)
    if default_group is None:
        raise ClinicalCoreFieldAuditError(f"No default field-group rule is defined for table: {table_name}")
    return default_group, f"default_table_group:{table_name}"


def load_schema_rows(schema_path: Path, expected_column_count: int) -> list[dict[str, str]]:
    schema_rows = read_tsv_dict_rows(schema_path)
    if len(schema_rows) != expected_column_count:
        raise ClinicalCoreFieldAuditError(
            f"Schema TSV row count {len(schema_rows)} does not match expected column count {expected_column_count}: {schema_path}"
        )

    required_columns = {"column_position", "raw_column_name", "alternate_column_name", "cde_id_raw"}
    if schema_rows:
        missing_columns = required_columns.difference(schema_rows[0].keys())
        if missing_columns:
            raise ClinicalCoreFieldAuditError(
                f"Schema TSV is missing required columns {sorted(missing_columns)}: {schema_path}"
            )

    for expected_position, row in enumerate(schema_rows, start=1):
        try:
            actual_position = int(row["column_position"])
        except ValueError as exc:
            raise ClinicalCoreFieldAuditError(
                f"Schema TSV column_position is not an integer in {schema_path}: {row}"
            ) from exc
        if actual_position != expected_position:
            raise ClinicalCoreFieldAuditError(
                f"Schema TSV column_position mismatch in {schema_path}: expected {expected_position}, found {actual_position}"
            )

    return schema_rows


def build_field_inventory_rows(
    *,
    audit_run_id: str,
    table_input: CoreTableInput,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    header, data_rows = read_table_rows(table_input.parsed_path)
    if len(header) != table_input.column_count:
        raise ClinicalCoreFieldAuditError(
            f"Parsed table header column count {len(header)} does not match manifest column count {table_input.column_count}: "
            f"{table_input.parsed_path}"
        )
    if len(data_rows) != table_input.row_count:
        raise ClinicalCoreFieldAuditError(
            f"Parsed table row count {len(data_rows)} does not match manifest row count {table_input.row_count}: "
            f"{table_input.parsed_path}"
        )

    schema_rows = load_schema_rows(table_input.schema_path, expected_column_count=table_input.column_count)

    field_inventory_rows: list[dict[str, Any]] = []
    for source_position, field_name in enumerate(header, start=1):
        schema_row = schema_rows[source_position - 1]
        schema_field_name = schema_row["raw_column_name"]
        if schema_field_name != field_name:
            raise ClinicalCoreFieldAuditError(
                f"Schema raw_column_name does not match parsed header for {table_input.table_name} position {source_position}: "
                f"{schema_field_name} vs {field_name}"
            )

        token_counts = {column_name: 0 for column_name in MISSING_TOKEN_COLUMNS.values()}
        non_missing_counter: Counter[str] = Counter()

        for row in data_rows:
            raw_value = row[source_position - 1]
            token_key = classify_missing_like_token(raw_value)
            if token_key is None:
                non_missing_counter[raw_value] += 1
            else:
                token_counts[MISSING_TOKEN_COLUMNS[token_key]] += 1

        missing_like_count = sum(token_counts.values())
        non_missing_count = len(data_rows) - missing_like_count
        missing_like_fraction = (missing_like_count / len(data_rows)) if data_rows else 0.0
        probable_field_group, group_assignment_rule = assign_probable_field_group(
            table_input.table_name,
            field_name,
        )

        field_inventory_rows.append(
            {
                "audit_run_id": audit_run_id,
                "parse_run_id": table_input.parse_run_id,
                "source_run_id": table_input.source_run_id,
                "table_name": table_input.table_name,
                "field_name": field_name,
                "source_position": source_position,
                "alternate_column_name": schema_row["alternate_column_name"],
                "cde_id_raw": schema_row["cde_id_raw"],
                "row_count": len(data_rows),
                "non_missing_count": non_missing_count,
                "missing_like_count": missing_like_count,
                "missing_like_fraction": round(missing_like_fraction, 6),
                "distinct_non_missing_count": len(non_missing_counter),
                "example_values_small_sample": select_example_values(non_missing_counter),
                "probable_field_group": probable_field_group,
                "group_assignment_rule": group_assignment_rule,
                **token_counts,
            }
        )

    if len(field_inventory_rows) != table_input.column_count:
        raise ClinicalCoreFieldAuditError(
            f"Generated field inventory row count {len(field_inventory_rows)} does not match column count "
            f"{table_input.column_count} for {table_input.table_name}"
        )

    table_summary = {
        "table_name": table_input.table_name,
        "source_filename": table_input.source_filename,
        "parsed_path": table_input.parsed_path,
        "schema_path": table_input.schema_path,
        "row_count": len(data_rows),
        "column_count": len(header),
        "schema_row_count": len(schema_rows),
        "field_inventory_row_count": len(field_inventory_rows),
        "non_missing_any_field_count": sum(1 for row in field_inventory_rows if int(row["non_missing_count"]) > 0),
        "high_missingness_field_count": sum(
            1 for row in field_inventory_rows if float(row["missing_like_fraction"]) >= 0.75
        ),
    }
    return field_inventory_rows, table_summary


def build_group_summary_rows(
    *,
    audit_run_id: str,
    parse_run_id: str,
    source_run_id: str,
    field_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    grouped_rows: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    table_field_counts: dict[str, int] = defaultdict(int)

    for row in field_rows:
        grouped_rows[(row["table_name"], row["probable_field_group"])].append(row)
        table_field_counts[row["table_name"]] += 1

    summary_rows: list[dict[str, Any]] = []
    for table_name in CORE_TABLE_NAMES:
        total_fields = table_field_counts[table_name]
        if total_fields <= 0:
            raise ClinicalCoreFieldAuditError(f"No field rows were collected for table: {table_name}")

        for group_name in FIELD_GROUP_LABELS:
            rows = sorted(
                grouped_rows.get((table_name, group_name), []),
                key=lambda item: int(item["source_position"]),
            )
            field_count = len(rows)
            summary_rows.append(
                {
                    "audit_run_id": audit_run_id,
                    "parse_run_id": parse_run_id,
                    "source_run_id": source_run_id,
                    "table_name": table_name,
                    "probable_field_group": group_name,
                    "field_count": field_count,
                    "field_fraction_of_table": round(field_count / total_fields, 6),
                    "non_missing_any_field_count": sum(1 for row in rows if int(row["non_missing_count"]) > 0),
                    "high_missingness_field_count": sum(
                        1 for row in rows if float(row["missing_like_fraction"]) >= 0.75
                    ),
                    "field_names_json": json.dumps(
                        [row["field_name"] for row in rows],
                        ensure_ascii=True,
                    ),
                }
            )

    return summary_rows


def build_missingness_summary_rows(field_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    projected_rows = [{field: row[field] for field in MISSINGNESS_SUMMARY_FIELDNAMES} for row in field_rows]
    return sorted(
        projected_rows,
        key=lambda row: (row["table_name"], -float(row["missing_like_fraction"]), int(row["source_position"])),
    )


def write_failure_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    audit_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    audit_run_dir = paths.audit_runs_root / audit_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = load_yaml(paths.trial_config)
        parse_pointer = load_json(paths.clinical_biotab_latest_pointer)

        parse_run_log_path = paths.repo_root / str(parse_pointer.get("run_log_json") or "")
        if not parse_run_log_path.exists():
            raise ClinicalCoreFieldAuditError(f"Clinical biotab parse run log does not exist: {parse_run_log_path}")
        parse_run_log = load_json(parse_run_log_path)

        manifest_path, core_tables = load_core_table_inputs(
            paths=paths,
            parse_pointer=parse_pointer,
            parse_run_log=parse_run_log,
        )

        create_run_directory(audit_run_dir)

        all_field_rows: list[dict[str, Any]] = []
        table_summaries: list[dict[str, Any]] = []
        field_inventory_paths: dict[str, str] = {}

        for table_input in core_tables:
            field_inventory_rows, table_summary = build_field_inventory_rows(
                audit_run_id=audit_run_id,
                table_input=table_input,
            )
            field_inventory_path = audit_run_dir / f"{table_input.table_name}_field_inventory.tsv"
            write_dict_rows_tsv(field_inventory_path, FIELD_AUDIT_FIELDNAMES, field_inventory_rows)
            field_inventory_paths[table_input.table_name] = repo_relative(field_inventory_path, paths.repo_root)
            all_field_rows.extend(field_inventory_rows)
            table_summaries.append(
                {
                    "table_name": table_input.table_name,
                    "source_filename": table_input.source_filename,
                    "parsed_path": repo_relative(table_summary["parsed_path"], paths.repo_root),
                    "schema_path": repo_relative(table_summary["schema_path"], paths.repo_root),
                    "row_count": table_summary["row_count"],
                    "column_count": table_summary["column_count"],
                    "schema_row_count": table_summary["schema_row_count"],
                    "field_inventory_path": repo_relative(field_inventory_path, paths.repo_root),
                    "field_inventory_row_count": table_summary["field_inventory_row_count"],
                    "non_missing_any_field_count": table_summary["non_missing_any_field_count"],
                    "high_missingness_field_count": table_summary["high_missingness_field_count"],
                }
            )

        combined_field_audit_path = audit_run_dir / "clinical_core_field_audit.tsv"
        write_dict_rows_tsv(combined_field_audit_path, FIELD_AUDIT_FIELDNAMES, all_field_rows)

        group_summary_rows = build_group_summary_rows(
            audit_run_id=audit_run_id,
            parse_run_id=core_tables[0].parse_run_id,
            source_run_id=core_tables[0].source_run_id,
            field_rows=all_field_rows,
        )
        group_summary_path = audit_run_dir / "clinical_core_field_group_summary.tsv"
        write_dict_rows_tsv(group_summary_path, GROUP_SUMMARY_FIELDNAMES, group_summary_rows)

        missingness_summary_rows = build_missingness_summary_rows(all_field_rows)
        missingness_summary_path = audit_run_dir / "clinical_core_missingness_summary.tsv"
        write_dict_rows_tsv(missingness_summary_path, MISSINGNESS_SUMMARY_FIELDNAMES, missingness_summary_rows)

        expected_field_count = sum(table.column_count for table in core_tables)
        combined_field_count = len(all_field_rows)
        combined_field_count_matches_expected = combined_field_count == expected_field_count
        field_inventory_counts_match = all(
            summary["field_inventory_row_count"] == summary["column_count"] for summary in table_summaries
        )
        row_count_positive = all(summary["row_count"] > 0 for summary in table_summaries)
        schema_row_counts_match = all(
            summary["schema_row_count"] == summary["column_count"] for summary in table_summaries
        )

        latest_pointer_payload = {
            "updated_at_utc": format_utc_timestamp(utc_now()),
            "audit_run_id": audit_run_id,
            "parse_run_id": core_tables[0].parse_run_id,
            "source_run_id": core_tables[0].source_run_id,
            "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
            "clinical_core_field_audit_tsv": repo_relative(combined_field_audit_path, paths.repo_root),
            "clinical_core_field_group_summary_tsv": repo_relative(group_summary_path, paths.repo_root),
            "clinical_core_missingness_summary_tsv": repo_relative(missingness_summary_path, paths.repo_root),
            "run_log_json": repo_relative(run_log_path, paths.repo_root),
            "clinical_biotab_latest_json": repo_relative(paths.clinical_biotab_latest_pointer, paths.repo_root),
            "field_inventory_tsvs": field_inventory_paths,
            "core_table_count": len(core_tables),
            "field_count": combined_field_count,
        }

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "audit_run_id": audit_run_id,
            "parse_run_id": core_tables[0].parse_run_id,
            "source_run_id": core_tables[0].source_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "clinical_biotab_latest_json": repo_relative(paths.clinical_biotab_latest_pointer, paths.repo_root),
                "clinical_biotab_manifest_tsv": repo_relative(manifest_path, paths.repo_root),
                "clinical_biotab_run_log_json": repo_relative(parse_run_log_path, paths.repo_root),
                "core_table_names": list(CORE_TABLE_NAMES),
            },
            "outputs": {
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "clinical_core_field_audit_tsv": repo_relative(combined_field_audit_path, paths.repo_root),
                "clinical_core_field_group_summary_tsv": repo_relative(group_summary_path, paths.repo_root),
                "clinical_core_missingness_summary_tsv": repo_relative(missingness_summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
                "field_inventory_tsvs": field_inventory_paths,
            },
            "validation": {
                "passed": (
                    row_count_positive
                    and schema_row_counts_match
                    and field_inventory_counts_match
                    and combined_field_count_matches_expected
                ),
                "required_core_table_count": len(CORE_TABLE_NAMES),
                "resolved_core_table_count": len(core_tables),
                "row_count_positive_for_all_tables": row_count_positive,
                "schema_row_counts_match_table_column_counts": schema_row_counts_match,
                "field_inventory_row_counts_match_table_column_counts": field_inventory_counts_match,
                "combined_field_count": combined_field_count,
                "expected_combined_field_count": expected_field_count,
                "combined_field_count_matches_expected": combined_field_count_matches_expected,
                "no_prior_run_overwrite": True,
            },
            "rules": {
                "missing_like_normalization": "value.strip().lower()",
                "missing_like_token_columns": MISSING_TOKEN_COLUMNS,
                "field_group_labels": list(FIELD_GROUP_LABELS),
                "identifier_admin_exact_names": sorted(IDENTIFIER_ADMIN_EXACT_NAMES),
                "identifier_admin_prefixes": list(IDENTIFIER_ADMIN_PREFIXES),
                "group_keyword_rules": [
                    {"probable_field_group": group_name, "patterns": list(patterns)}
                    for group_name, patterns in GROUP_KEYWORD_RULES
                ],
                "default_group_by_table": DEFAULT_GROUP_BY_TABLE,
                "high_missingness_threshold": 0.75,
                "example_value_sample_limit": 5,
            },
            "upstream_snapshots": {
                "clinical_biotab_latest_pointer": parse_pointer,
                "clinical_biotab_run_log": parse_run_log,
            },
            "tables": table_summaries,
            "latest_pointer": latest_pointer_payload,
        }

        write_json(run_log_path, run_log_payload)
        write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "audit_run_id": audit_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_clinical_core_field_audit",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA clinical core field audit complete.")
    print(f"Audit run ID: {run_log['audit_run_id']}")
    print(f"Parse run ID: {run_log['parse_run_id']}")
    print(f"Source run ID: {run_log['source_run_id']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"Combined field audit TSV: {run_log['outputs']['clinical_core_field_audit_tsv']}")
    print(f"Field group summary TSV: {run_log['outputs']['clinical_core_field_group_summary_tsv']}")
    print(f"Missingness summary TSV: {run_log['outputs']['clinical_core_missingness_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
