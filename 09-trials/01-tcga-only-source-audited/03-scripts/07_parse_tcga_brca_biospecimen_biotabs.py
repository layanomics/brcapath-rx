#!/usr/bin/env python
"""Parse TCGA-BRCA core Biospecimen Supplement biotab tables into auditable TSV outputs."""

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

SOURCE_CLASS = "biospecimen"
BIOTAB_FORMAT = "BCR Biotab"
HEADER_PREAMBLE_ROWS = 2
DATA_START_LINE = 3
EXPECTED_DOWNLOAD_SUBDIR = Path("01-data/raw/tcga-brca/gdc/downloads/biospecimen")
CORE_FILENAME_PREFIX = "nationwidechildrens.org_biospecimen_"
SSF_FILENAME_PREFIX = "nationwidechildrens.org_ssf_"
SOURCE_FILENAME_SUFFIX = "_brca.txt"
SOURCE_INVENTORY_FIELDNAMES = [
    "parse_run_id",
    "source_run_id",
    "source_class",
    "source_file_id",
    "source_filename",
    "source_path",
    "data_format",
    "case_count",
    "scope_status",
    "scope_reason",
    "table_name",
    "output_filename",
]
MANIFEST_FIELDNAMES = [
    "parse_run_id",
    "source_run_id",
    "source_class",
    "table_name",
    "source_file_id",
    "source_filename",
    "source_path",
    "output_filename",
    "output_path",
    "schema_path",
    "row_count",
    "column_count",
    "header_preamble_rows",
    "data_start_line",
    "parse_notes",
]
SCHEMA_FIELDNAMES = [
    "parse_run_id",
    "table_name",
    "column_position",
    "raw_column_name",
    "alternate_column_name",
    "cde_id_raw",
]


class BiospecimenBiotabParseError(RuntimeError):
    """Raised when the biospecimen biotab parsing workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the biospecimen biotab parsing workflow."""

    repo_root: Path
    trial_config: Path
    supplement_latest_pointer: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class DiscoveredBiotabSource:
    """One discovered biospecimen biotab source file and its scope status."""

    source_run_id: str
    source_file_id: str
    source_filename: str
    source_path: Path
    data_format: str
    case_count: str
    scope_status: str
    scope_reason: str
    output_filename: str

    @property
    def table_name(self) -> str:
        return Path(self.output_filename).stem if self.output_filename else ""


@dataclass(frozen=True)
class ParsedBiotabTable:
    """Parsed biospecimen table content and its two-row source preamble."""

    raw_header: list[str]
    alternate_header: list[str]
    cde_row: list[str]
    data_rows: list[list[str]]

    @property
    def row_count(self) -> int:
        return len(self.data_rows)

    @property
    def column_count(self) -> int:
        return len(self.raw_header)


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise BiospecimenBiotabParseError("Unable to locate the repository root from the script path.")


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
                raise BiospecimenBiotabParseError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise BiospecimenBiotabParseError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise BiospecimenBiotabParseError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise BiospecimenBiotabParseError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise BiospecimenBiotabParseError(
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
        raise BiospecimenBiotabParseError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise BiospecimenBiotabParseError(f"Expected a JSON object in file: {path}")
    return data


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise BiospecimenBiotabParseError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    if path.exists():
        raise BiospecimenBiotabParseError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise BiospecimenBiotabParseError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def build_workflow_paths() -> WorkflowPaths:
    repo_root = detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    supplement_latest_pointer = (
        repo_root / "01-data" / "audit" / "tcga-brca" / "source" / "tcga_brca_source_supplements_latest.json"
    )

    trial_config_data = load_yaml(trial_config)
    processed_root = repo_root / str(trial_config_data.get("processed_data_root", "01-data/processed"))
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        supplement_latest_pointer=supplement_latest_pointer,
        processed_runs_root=processed_root / "tcga-brca" / "biospecimen" / "biotab_runs",
        audit_runs_root=audit_root / "tcga-brca" / "variables" / "biospecimen_biotab_runs",
        latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_biospecimen_biotabs_latest.json",
    )


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise BiospecimenBiotabParseError(f"Run directory already exists: {path}")
    path.mkdir(parents=False, exist_ok=False)
    return path


def validate_supplement_pointer(
    pointer: dict[str, Any],
    pointer_path: Path,
    repo_root: Path,
) -> tuple[str, Path, Path, Path]:
    source_run_id = str(pointer.get("run_id") or "").strip()
    if not source_run_id:
        raise BiospecimenBiotabParseError(f"Supplement pointer is missing run_id: {pointer_path}")

    source_classes = pointer.get("source_classes")
    if not isinstance(source_classes, dict) or SOURCE_CLASS not in source_classes:
        raise BiospecimenBiotabParseError(
            f"Supplement pointer is missing source class '{SOURCE_CLASS}': {pointer_path}"
        )

    biospecimen_pointer = source_classes[SOURCE_CLASS]
    if not isinstance(biospecimen_pointer, dict):
        raise BiospecimenBiotabParseError(f"Biospecimen source-class payload is invalid in {pointer_path}")
    if not bool(biospecimen_pointer.get("download_completed")):
        raise BiospecimenBiotabParseError(
            "Biospecimen supplement download is not marked complete in the latest supplement pointer."
        )

    download_dir = repo_root / str(biospecimen_pointer.get("download_dir") or "")
    metadata_tsv = repo_root / str(pointer.get("metadata_tsv") or "")
    source_run_log = repo_root / str(pointer.get("run_log_json") or "")

    for required_path, label in [
        (download_dir, "biospecimen download directory"),
        (metadata_tsv, "supplement metadata TSV"),
        (source_run_log, "supplement run log"),
    ]:
        if not required_path.exists():
            raise BiospecimenBiotabParseError(f"Expected {label} does not exist: {required_path}")

    expected_download_root = repo_root / EXPECTED_DOWNLOAD_SUBDIR
    if expected_download_root not in [download_dir, *download_dir.parents]:
        raise BiospecimenBiotabParseError(
            f"Biospecimen download directory is outside the expected biospecimen raw root: {download_dir}"
        )

    return source_run_id, download_dir, metadata_tsv, source_run_log


def read_metadata_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required_columns = {"source_class", "file_id", "file_name", "data_format", "case_count"}
        if not required_columns.issubset(fieldnames):
            raise BiospecimenBiotabParseError(
                f"Supplement metadata TSV is missing required columns {sorted(required_columns)}: {path}"
            )
        return list(reader)


def build_output_filename(source_filename: str) -> str:
    prefix = "nationwidechildrens.org_"
    if not source_filename.startswith(prefix) or not source_filename.endswith(SOURCE_FILENAME_SUFFIX):
        raise BiospecimenBiotabParseError(
            f"Cannot derive an output filename from unexpected biospecimen source name: {source_filename}"
        )
    slug = source_filename[len(prefix) : -len(SOURCE_FILENAME_SUFFIX)].replace(".", "_")
    if not slug:
        raise BiospecimenBiotabParseError(f"Derived an empty table slug from source name: {source_filename}")
    return f"{slug}.tsv"


def classify_source_scope(source_filename: str) -> tuple[str, str]:
    if source_filename.startswith(CORE_FILENAME_PREFIX) and source_filename.endswith(SOURCE_FILENAME_SUFFIX):
        return (
            "parsed",
            "Core biospecimen hierarchy biotab table included in the first identifier-linkage parsing workflow.",
        )
    if source_filename.startswith(SSF_FILENAME_PREFIX) and source_filename.endswith(SOURCE_FILENAME_SUFFIX):
        return (
            "excluded",
            "Deferred in v1: SSF text biotab tables are outside the core biospecimen hierarchy scope.",
        )
    return ("excluded", "Excluded: biospecimen biotab filename family is outside the approved core v1 scope.")


def discover_biospecimen_biotab_sources(
    *,
    metadata_rows: list[dict[str, str]],
    source_run_id: str,
    download_dir: Path,
) -> tuple[list[DiscoveredBiotabSource], list[DiscoveredBiotabSource]]:
    discovered: list[DiscoveredBiotabSource] = []
    parsed: list[DiscoveredBiotabSource] = []
    output_filenames: set[str] = set()

    for row in metadata_rows:
        if row["source_class"] != SOURCE_CLASS:
            continue
        if row["data_format"] != BIOTAB_FORMAT:
            continue

        source_file_id = str(row.get("file_id") or "").strip()
        source_filename = str(row.get("file_name") or "").strip()
        if not source_file_id or not source_filename:
            raise BiospecimenBiotabParseError(
                "Encountered a biospecimen biotab metadata row missing file_id or file_name."
            )

        source_path = download_dir / source_file_id / source_filename
        if not source_path.exists():
            raise BiospecimenBiotabParseError(f"Biospecimen biotab source file does not exist on disk: {source_path}")

        scope_status, scope_reason = classify_source_scope(source_filename)
        output_filename = ""
        if scope_status == "parsed":
            output_filename = build_output_filename(source_filename)
            if output_filename in output_filenames:
                raise BiospecimenBiotabParseError(
                    f"Duplicate output filename derived from biospecimen biotab sources: {output_filename}"
                )
            output_filenames.add(output_filename)

        discovered_source = DiscoveredBiotabSource(
            source_run_id=source_run_id,
            source_file_id=source_file_id,
            source_filename=source_filename,
            source_path=source_path,
            data_format=str(row.get("data_format") or "").strip(),
            case_count=str(row.get("case_count") or "").strip(),
            scope_status=scope_status,
            scope_reason=scope_reason,
            output_filename=output_filename,
        )
        discovered.append(discovered_source)
        if scope_status == "parsed":
            parsed.append(discovered_source)

    if not discovered:
        raise BiospecimenBiotabParseError(
            "No biospecimen BCR Biotab source files were discovered in the supplement metadata TSV."
        )
    if not parsed:
        raise BiospecimenBiotabParseError(
            "No in-scope core biospecimen biotab source files were discovered in the supplement metadata TSV."
        )

    return (
        sorted(discovered, key=lambda item: item.source_filename),
        sorted(parsed, key=lambda item: item.source_filename),
    )


def parse_biotab_file(source: DiscoveredBiotabSource) -> ParsedBiotabTable:
    with source.source_path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    if len(rows) <= HEADER_PREAMBLE_ROWS:
        raise BiospecimenBiotabParseError(
            f"Biospecimen biotab source file does not contain data rows after the {HEADER_PREAMBLE_ROWS}-line preamble: "
            f"{source.source_path}"
        )

    raw_header = rows[0]
    cde_row = rows[1]
    data_rows = rows[HEADER_PREAMBLE_ROWS:]

    column_count = len(raw_header)
    if column_count <= 0:
        raise BiospecimenBiotabParseError(f"Biospecimen biotab source file has an empty header row: {source.source_path}")
    if len(cde_row) != column_count:
        raise BiospecimenBiotabParseError(
            f"CDE row has {len(cde_row)} columns, expected {column_count}: {source.source_path}"
        )
    if not data_rows:
        raise BiospecimenBiotabParseError(f"Biospecimen biotab source file has zero parsed data rows: {source.source_path}")

    for line_number, row in enumerate(data_rows, start=DATA_START_LINE):
        if len(row) != column_count:
            raise BiospecimenBiotabParseError(
                f"Data row {line_number} has {len(row)} columns, expected {column_count}: {source.source_path}"
            )

    # Core biospecimen biotabs expose exact source headers on line 1 and sparse CDE metadata on line 2.
    return ParsedBiotabTable(
        raw_header=raw_header,
        alternate_header=[""] * column_count,
        cde_row=cde_row,
        data_rows=data_rows,
    )


def build_source_inventory_rows(
    parse_run_id: str,
    repo_root: Path,
    discovered_sources: list[DiscoveredBiotabSource],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source in discovered_sources:
        rows.append(
            {
                "parse_run_id": parse_run_id,
                "source_run_id": source.source_run_id,
                "source_class": SOURCE_CLASS,
                "source_file_id": source.source_file_id,
                "source_filename": source.source_filename,
                "source_path": repo_relative(source.source_path, repo_root),
                "data_format": source.data_format,
                "case_count": source.case_count,
                "scope_status": source.scope_status,
                "scope_reason": source.scope_reason,
                "table_name": source.table_name,
                "output_filename": source.output_filename,
            }
        )
    return rows


def build_schema_rows(parse_run_id: str, table_name: str, parsed_table: ParsedBiotabTable) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for position, (raw_name, alternate_name, cde_raw) in enumerate(
        zip(parsed_table.raw_header, parsed_table.alternate_header, parsed_table.cde_row),
        start=1,
    ):
        rows.append(
            {
                "parse_run_id": parse_run_id,
                "table_name": table_name,
                "column_position": position,
                "raw_column_name": raw_name,
                "alternate_column_name": alternate_name,
                "cde_id_raw": cde_raw,
            }
        )
    return rows


def write_failure_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    parse_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    processed_run_dir = paths.processed_runs_root / parse_run_id
    audit_run_dir = paths.audit_runs_root / parse_run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        source_pointer = load_json(paths.supplement_latest_pointer)
        source_run_id, download_dir, metadata_tsv, source_run_log = validate_supplement_pointer(
            pointer=source_pointer,
            pointer_path=paths.supplement_latest_pointer,
            repo_root=paths.repo_root,
        )
        source_run_log_payload = load_json(source_run_log)
        trial_config = load_yaml(paths.trial_config)

        create_run_directory(processed_run_dir)
        create_run_directory(audit_run_dir)

        metadata_rows = read_metadata_rows(metadata_tsv)
        discovered_sources, parsed_sources = discover_biospecimen_biotab_sources(
            metadata_rows=metadata_rows,
            source_run_id=source_run_id,
            download_dir=download_dir,
        )

        source_inventory_rows = build_source_inventory_rows(
            parse_run_id=parse_run_id,
            repo_root=paths.repo_root,
            discovered_sources=discovered_sources,
        )
        source_inventory_path = audit_run_dir / "biospecimen_biotab_source_inventory.tsv"
        write_dict_rows_tsv(source_inventory_path, SOURCE_INVENTORY_FIELDNAMES, source_inventory_rows)

        manifest_rows: list[dict[str, Any]] = []
        parsed_table_summaries: list[dict[str, Any]] = []
        schema_files: list[str] = []
        parse_notes = (
            "Primary TSV preserves source line-1 headers exactly; source line-2 sparse CDE metadata is stored in the "
            "schema TSV. No alternate header row exists for this core biospecimen v1 workflow."
        )

        for source in parsed_sources:
            parsed = parse_biotab_file(source)

            output_path = processed_run_dir / source.output_filename
            write_tsv(output_path, parsed.raw_header, parsed.data_rows)

            schema_filename = f"schema_{source.table_name}.tsv"
            schema_path = audit_run_dir / schema_filename
            schema_rows = build_schema_rows(parse_run_id=parse_run_id, table_name=source.table_name, parsed_table=parsed)
            write_dict_rows_tsv(schema_path, SCHEMA_FIELDNAMES, schema_rows)
            schema_files.append(schema_filename)

            manifest_rows.append(
                {
                    "parse_run_id": parse_run_id,
                    "source_run_id": source.source_run_id,
                    "source_class": SOURCE_CLASS,
                    "table_name": source.table_name,
                    "source_file_id": source.source_file_id,
                    "source_filename": source.source_filename,
                    "source_path": repo_relative(source.source_path, paths.repo_root),
                    "output_filename": source.output_filename,
                    "output_path": repo_relative(output_path, paths.repo_root),
                    "schema_path": repo_relative(schema_path, paths.repo_root),
                    "row_count": parsed.row_count,
                    "column_count": parsed.column_count,
                    "header_preamble_rows": HEADER_PREAMBLE_ROWS,
                    "data_start_line": DATA_START_LINE,
                    "parse_notes": parse_notes,
                }
            )

            parsed_table_summaries.append(
                {
                    "table_name": source.table_name,
                    "source_file_id": source.source_file_id,
                    "source_filename": source.source_filename,
                    "source_path": repo_relative(source.source_path, paths.repo_root),
                    "output_filename": source.output_filename,
                    "output_path": repo_relative(output_path, paths.repo_root),
                    "schema_path": repo_relative(schema_path, paths.repo_root),
                    "row_count": parsed.row_count,
                    "column_count": parsed.column_count,
                }
            )

        manifest_path = audit_run_dir / "biospecimen_biotab_table_manifest.tsv"
        write_dict_rows_tsv(manifest_path, MANIFEST_FIELDNAMES, manifest_rows)

        latest_pointer_payload = {
            "updated_at_utc": format_utc_timestamp(utc_now()),
            "parse_run_id": parse_run_id,
            "source_run_id": source_run_id,
            "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
            "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
            "source_inventory_tsv": repo_relative(source_inventory_path, paths.repo_root),
            "table_manifest_tsv": repo_relative(manifest_path, paths.repo_root),
            "run_log_json": repo_relative(run_log_path, paths.repo_root),
            "source_supplement_latest_json": repo_relative(paths.supplement_latest_pointer, paths.repo_root),
            "source_metadata_tsv": repo_relative(metadata_tsv, paths.repo_root),
            "parsed_table_count": len(parsed_table_summaries),
            "discovered_biotab_file_count": len(discovered_sources),
            "excluded_biotab_file_count": len(
                [item for item in discovered_sources if item.scope_status != "parsed"]
            ),
        }

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "parse_run_id": parse_run_id,
            "source_run_id": source_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "supplement_latest_json": repo_relative(paths.supplement_latest_pointer, paths.repo_root),
                "source_metadata_tsv": repo_relative(metadata_tsv, paths.repo_root),
                "source_run_log_json": repo_relative(source_run_log, paths.repo_root),
                "biospecimen_download_dir": repo_relative(download_dir, paths.repo_root),
            },
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "source_inventory_tsv": repo_relative(source_inventory_path, paths.repo_root),
                "table_manifest_tsv": repo_relative(manifest_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": True,
                "biospecimen_download_completed": True,
                "discovered_biotab_file_count": len(discovered_sources),
                "parsed_table_count": len(parsed_table_summaries),
                "excluded_biotab_file_count": len(
                    [item for item in discovered_sources if item.scope_status != "parsed"]
                ),
                "header_preamble_rows": HEADER_PREAMBLE_ROWS,
                "data_start_line": DATA_START_LINE,
                "row_count_positive_for_all_tables": all(item["row_count"] > 0 for item in parsed_table_summaries),
                "unique_output_filenames": len({item["output_filename"] for item in parsed_table_summaries})
                == len(parsed_table_summaries),
                "schema_file_count": len(schema_files),
                "source_inventory_row_count": len(source_inventory_rows),
                "core_scope_only": all(
                    item["source_filename"].startswith(CORE_FILENAME_PREFIX)
                    for item in parsed_table_summaries
                ),
            },
            "source_pointer_snapshot": source_pointer,
            "source_run_log_snapshot": source_run_log_payload,
            "source_inventory": source_inventory_rows,
            "tables": parsed_table_summaries,
            "latest_pointer": latest_pointer_payload,
        }

        write_json(run_log_path, run_log_payload)
        write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "parse_run_id": parse_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_biospecimen_biotab_parsing",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA core biospecimen biotab parsing complete.")
    print(f"Parse run ID: {run_log['parse_run_id']}")
    print(f"Source run ID: {run_log['source_run_id']}")
    print(f"Discovered biospecimen biotab text files: {run_log['validation']['discovered_biotab_file_count']}")
    print(f"Parsed core biospecimen tables: {run_log['validation']['parsed_table_count']}")
    print(f"Excluded biospecimen biotab text files: {run_log['validation']['excluded_biotab_file_count']}")
    print(f"Processed output directory: {run_log['outputs']['processed_run_directory']}")
    print(f"Audit output directory: {run_log['outputs']['audit_run_directory']}")
    print(f"Source inventory: {run_log['outputs']['source_inventory_tsv']}")
    print(f"Table manifest: {run_log['outputs']['table_manifest_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
