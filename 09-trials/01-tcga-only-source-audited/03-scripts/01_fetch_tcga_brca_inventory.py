#!/usr/bin/env python
"""Fetch source-audited TCGA-BRCA inventory snapshots from the GDC API.

This script is intentionally limited to the inventory step of the
TCGA-only, source-audited BRCAPath-Rx rebuild. It:

1. Reads the repo's existing config and payload files.
2. Queries the GDC API for TCGA-BRCA case and file inventories.
3. Saves immutable raw page responses and effective payloads for each run.
4. Writes reviewer-friendly normalized TSV inventories.
5. Updates a stable "latest run" manifest for notebook review.

It does not download raw data files, parse XML, define cohorts, or choose
clinical endpoints.
"""

from __future__ import annotations

import csv
import json
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - exercised in environments without PyYAML
    yaml = None

GDC_API_BASE_URL = "https://api.gdc.cancer.gov"
PROJECT_ID = "TCGA-BRCA"
PAGE_SIZE = 2000
REQUEST_TIMEOUT_SECONDS = 120

CASE_TSV_COLUMNS = [
    "case_id",
    "submitter_id",
    "project_id",
    "primary_site",
    "disease_type",
]

FILE_TSV_COLUMNS = [
    "file_id",
    "file_name",
    "access",
    "state",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "platform",
    "file_size",
    "analysis_workflow_type",
    "case_ids",
    "case_submitter_ids",
    "project_ids",
    "sample_ids",
    "sample_submitter_ids",
    "sample_types",
    "tissue_types",
    "tumor_descriptors",
    "case_count",
    "sample_count",
]

LATEST_MANIFEST_NAME = "tcga_brca_inventory_latest.json"


class InventoryFetchError(RuntimeError):
    """Raised when the inventory workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the workflow."""

    repo_root: Path
    gdc_client_config: Path
    project_status_config: Path
    trial_config: Path
    cases_payload_config: Path
    files_payload_config: Path
    source_root: Path
    runs_root: Path
    latest_manifest: Path


def utc_now() -> datetime:
    """Return the current UTC time with timezone information."""

    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    """Return an ISO-like UTC timestamp without fractional seconds."""

    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    """Walk upward until the git root is found."""

    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise InventoryFetchError("Unable to locate the repository root from the script path.")


def build_workflow_paths() -> WorkflowPaths:
    """Build the repo-relative path map used by the inventory workflow."""

    repo_root = detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    source_root = repo_root / "01-data" / "audit" / "tcga-brca" / "source"
    return WorkflowPaths(
        repo_root=repo_root,
        gdc_client_config=repo_root / "07-config" / "gdc_client.yaml",
        project_status_config=repo_root / "07-config" / "project_status.yaml",
        trial_config=trial_root / "04-config" / "trial_config.yaml",
        cases_payload_config=trial_root / "04-config" / "gdc_cases_inventory_payload.json",
        files_payload_config=trial_root / "04-config" / "gdc_files_inventory_payload.json",
        source_root=source_root,
        runs_root=source_root / "runs",
        latest_manifest=source_root / LATEST_MANIFEST_NAME,
    )


def repo_relative(path: Path, repo_root: Path) -> str:
    """Convert a path to a forward-slash repo-relative string."""

    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def load_yaml(path: Path) -> dict[str, Any]:
    """Read a YAML configuration file."""

    with path.open("r", encoding="utf-8") as handle:
        raw_text = handle.read()
    if yaml is not None:
        data = yaml.safe_load(raw_text)
    else:
        data = parse_simple_yaml_mapping(raw_text, source_path=path)
    if not isinstance(data, dict):
        raise InventoryFetchError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    """Read a JSON file and confirm it contains an object."""

    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise InventoryFetchError(f"Expected a JSON object in file: {path}")
    return data


def parse_simple_yaml_mapping(raw_text: str, source_path: Path) -> dict[str, str]:
    """Parse the repo's simple top-level key/value YAML files without PyYAML."""

    parsed: dict[str, str] = {}
    for line_number, raw_line in enumerate(raw_text.splitlines(), start=1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if raw_line.startswith((" ", "\t")):
            raise InventoryFetchError(
                f"Nested YAML is not supported by the built-in fallback parser: {source_path}:{line_number}"
            )
        if ":" not in raw_line:
            raise InventoryFetchError(
                f"Expected a top-level key/value pair in YAML: {source_path}:{line_number}"
            )
        key, value = raw_line.split(":", 1)
        key = key.strip()
        value = value.strip()
        if value.startswith(("'", '"')) and value.endswith(("'", '"')) and len(value) >= 2:
            value = value[1:-1]
        parsed[key] = value
    return parsed


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    """Write structured JSON with explicit overwrite control."""

    if path.exists() and not overwrite:
        raise InventoryFetchError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_text(path: Path, text: str, overwrite: bool = False) -> None:
    """Write text with explicit overwrite control."""

    if path.exists() and not overwrite:
        raise InventoryFetchError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def create_run_directory(runs_root: Path, run_id: str) -> Path:
    """Create a new immutable run directory."""

    runs_root.mkdir(parents=True, exist_ok=True)
    run_dir = runs_root / run_id
    if run_dir.exists():
        raise InventoryFetchError(f"Run directory already exists: {run_dir}")
    run_dir.mkdir(parents=False, exist_ok=False)
    return run_dir


def build_effective_payload(base_payload: dict[str, Any], sort_value: str, offset: int) -> dict[str, Any]:
    """Apply runtime inventory overrides to the checked-in payload template."""

    payload = deepcopy(base_payload)
    payload["format"] = "JSON"
    payload["size"] = str(PAGE_SIZE)
    payload["from"] = offset
    payload["sort"] = sort_value
    return payload


def get_gdc_status(session: requests.Session) -> tuple[dict[str, Any], str]:
    """Fetch the GDC API status payload and raw text."""

    response = session.get(f"{GDC_API_BASE_URL}/status", timeout=REQUEST_TIMEOUT_SECONDS)
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        raise InventoryFetchError(f"GDC status request failed with HTTP {response.status_code}.") from exc
    try:
        return response.json(), response.text
    except json.JSONDecodeError as exc:
        raise InventoryFetchError("Unable to decode the GDC status response as JSON.") from exc


def post_gdc_payload(
    session: requests.Session, endpoint: str, payload: dict[str, Any]
) -> tuple[dict[str, Any], str]:
    """POST an inventory query to a GDC API endpoint."""

    response = session.post(
        f"{GDC_API_BASE_URL}/{endpoint}",
        json=payload,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        body_preview = response.text[:500]
        raise InventoryFetchError(
            f"GDC {endpoint} request failed with HTTP {response.status_code}. Body preview: {body_preview}"
        ) from exc
    try:
        return response.json(), response.text
    except json.JSONDecodeError as exc:
        raise InventoryFetchError(f"Unable to decode the GDC {endpoint} response as JSON.") from exc


def unique_preserving_order(values: list[str]) -> list[str]:
    """Drop blanks and duplicates while preserving the first-seen order."""

    seen: set[str] = set()
    result: list[str] = []
    for value in values:
        normalized = str(value).strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        result.append(normalized)
    return result


def to_json_array_cell(values: list[str]) -> str:
    """Serialize list-like TSV columns as compact JSON arrays."""

    return json.dumps(values, ensure_ascii=True, separators=(",", ":"))


def normalize_case_hit(hit: dict[str, Any]) -> dict[str, str]:
    """Convert one raw case hit into the normalized case inventory row."""

    project_id = (hit.get("project") or {}).get("project_id")
    return {
        "case_id": str(hit.get("case_id") or hit.get("id") or "").strip(),
        "submitter_id": str(hit.get("submitter_id") or "").strip(),
        "project_id": str(project_id or "").strip(),
        "primary_site": str(hit.get("primary_site") or "").strip(),
        "disease_type": str(hit.get("disease_type") or "").strip(),
    }


def normalize_file_hit(hit: dict[str, Any]) -> dict[str, str]:
    """Convert one raw file hit into the normalized file inventory row."""

    case_ids: list[str] = []
    case_submitter_ids: list[str] = []
    project_ids: list[str] = []
    sample_records: list[dict[str, str]] = []
    sample_record_keys: set[tuple[str, str, str, str, str]] = set()

    for case in hit.get("cases") or []:
        case_ids.extend(unique_preserving_order([str(case.get("case_id") or "").strip()]))
        case_submitter_ids.extend(unique_preserving_order([str(case.get("submitter_id") or "").strip()]))
        project_id = str((case.get("project") or {}).get("project_id") or "").strip()
        project_ids.extend(unique_preserving_order([project_id]))

        for sample in case.get("samples") or []:
            sample_record = {
                "sample_id": str(sample.get("sample_id") or "").strip(),
                "sample_submitter_id": str(sample.get("submitter_id") or "").strip(),
                "sample_type": str(sample.get("sample_type") or "").strip(),
                "tissue_type": str(sample.get("tissue_type") or "").strip(),
                "tumor_descriptor": str(sample.get("tumor_descriptor") or "").strip(),
            }
            sample_record_key = (
                sample_record["sample_id"],
                sample_record["sample_submitter_id"],
                sample_record["sample_type"],
                sample_record["tissue_type"],
                sample_record["tumor_descriptor"],
            )
            if sample_record_key in sample_record_keys:
                continue
            sample_record_keys.add(sample_record_key)
            sample_records.append(sample_record)

    unique_case_ids = unique_preserving_order(case_ids)
    unique_case_submitter_ids = unique_preserving_order(case_submitter_ids)
    unique_project_ids = unique_preserving_order(project_ids)
    sample_ids = unique_preserving_order([record["sample_id"] for record in sample_records])
    sample_submitter_ids = unique_preserving_order([record["sample_submitter_id"] for record in sample_records])
    sample_types = unique_preserving_order([record["sample_type"] for record in sample_records])
    tissue_types = unique_preserving_order([record["tissue_type"] for record in sample_records])
    tumor_descriptors = unique_preserving_order([record["tumor_descriptor"] for record in sample_records])

    return {
        "file_id": str(hit.get("file_id") or hit.get("id") or "").strip(),
        "file_name": str(hit.get("file_name") or "").strip(),
        "access": str(hit.get("access") or "").strip(),
        "state": str(hit.get("state") or "").strip(),
        "data_category": str(hit.get("data_category") or "").strip(),
        "data_type": str(hit.get("data_type") or "").strip(),
        "data_format": str(hit.get("data_format") or "").strip(),
        "experimental_strategy": str(hit.get("experimental_strategy") or "").strip(),
        "platform": str(hit.get("platform") or "").strip(),
        "file_size": str(hit.get("file_size") or "").strip(),
        "analysis_workflow_type": str((hit.get("analysis") or {}).get("workflow_type") or "").strip(),
        "case_ids": to_json_array_cell(unique_case_ids),
        "case_submitter_ids": to_json_array_cell(unique_case_submitter_ids),
        "project_ids": to_json_array_cell(unique_project_ids),
        "sample_ids": to_json_array_cell(sample_ids),
        "sample_submitter_ids": to_json_array_cell(sample_submitter_ids),
        "sample_types": to_json_array_cell(sample_types),
        "tissue_types": to_json_array_cell(tissue_types),
        "tumor_descriptors": to_json_array_cell(tumor_descriptors),
        "case_count": str(len(unique_case_ids)),
        "sample_count": str(len(sample_records)),
    }


def validate_response_page(
    endpoint: str,
    payload_offset: int,
    payload_size: int,
    response_json: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Validate a raw API page before normalizing it."""

    data = response_json.get("data")
    if not isinstance(data, dict):
        raise InventoryFetchError(f"GDC {endpoint} response is missing the top-level data object.")

    hits = data.get("hits")
    pagination = data.get("pagination")
    if not isinstance(hits, list) or not isinstance(pagination, dict):
        raise InventoryFetchError(f"GDC {endpoint} response is missing hits or pagination.")

    try:
        total = int(pagination["total"])
        count = int(pagination["count"])
        current_from = int(pagination["from"])
        page = int(pagination["page"])
        pages = int(pagination["pages"])
    except (KeyError, TypeError, ValueError) as exc:
        raise InventoryFetchError(f"GDC {endpoint} pagination metadata is incomplete: {pagination}") from exc

    if current_from != payload_offset:
        raise InventoryFetchError(
            f"GDC {endpoint} returned pagination.from={current_from} for requested offset {payload_offset}."
        )
    if count != len(hits):
        raise InventoryFetchError(
            f"GDC {endpoint} pagination.count={count} does not match the number of hits ({len(hits)})."
        )
    if count > payload_size:
        raise InventoryFetchError(
            f"GDC {endpoint} returned {count} rows for page size {payload_size}, which is unexpected."
        )
    return hits, {
        "total": total,
        "count": count,
        "from": current_from,
        "page": page,
        "pages": pages,
    }


def write_inventory_tsv(
    session: requests.Session,
    endpoint: str,
    base_payload: dict[str, Any],
    sort_value: str,
    payload_output_path: Path,
    response_prefix: str,
    tsv_output_path: Path,
    fieldnames: list[str],
    normalizer: Any,
) -> dict[str, Any]:
    """Fetch one endpoint with pagination and write a normalized TSV."""

    effective_payload = build_effective_payload(base_payload=base_payload, sort_value=sort_value, offset=0)
    write_json(payload_output_path, effective_payload)

    response_paths: list[str] = []
    rows_written = 0
    total_expected: int | None = None
    page_index = 1
    offset = 0
    pages_expected: int | None = None

    with tsv_output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()

        while True:
            page_payload = build_effective_payload(base_payload=base_payload, sort_value=sort_value, offset=offset)
            response_json, response_text = post_gdc_payload(session=session, endpoint=endpoint, payload=page_payload)

            page_path = payload_output_path.parent / f"{response_prefix}_page_{page_index:04d}.json"
            write_text(page_path, response_text)
            response_paths.append(page_path.name)

            hits, pagination = validate_response_page(
                endpoint=endpoint,
                payload_offset=offset,
                payload_size=PAGE_SIZE,
                response_json=response_json,
            )
            if total_expected is None:
                total_expected = pagination["total"]
                pages_expected = pagination["pages"]
            elif total_expected != pagination["total"]:
                raise InventoryFetchError(
                    f"GDC {endpoint} changed pagination.total during paging ({total_expected} -> {pagination['total']})."
                )

            for hit in hits:
                row = normalizer(hit)
                writer.writerow(row)
                rows_written += 1

            if rows_written == total_expected:
                break

            if pagination["count"] == 0:
                raise InventoryFetchError(
                    f"GDC {endpoint} returned an empty page before the expected total ({rows_written}/{total_expected})."
                )
            offset += pagination["count"]
            page_index += 1

    if total_expected is None:
        raise InventoryFetchError(f"GDC {endpoint} did not return any pagination metadata.")
    if rows_written != total_expected:
        raise InventoryFetchError(
            f"Normalized {endpoint} row count {rows_written} does not match expected total {total_expected}."
        )

    return {
        "endpoint": endpoint,
        "rows_written": rows_written,
        "total_expected": total_expected,
        "pages_written": page_index,
        "pages_expected": pages_expected,
        "response_files": response_paths,
        "sort": sort_value,
        "page_size": PAGE_SIZE,
        "tsv_output_path": str(tsv_output_path),
    }


def validate_case_tsv(path: Path) -> dict[str, Any]:
    """Run post-write validation on the normalized case inventory TSV."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing_columns = [column for column in CASE_TSV_COLUMNS if column not in fieldnames]
        if missing_columns:
            raise InventoryFetchError(f"Case inventory is missing required columns: {missing_columns}")

        row_count = 0
        unique_submitters: set[str] = set()
        for row in reader:
            row_count += 1
            if row["project_id"] != PROJECT_ID:
                raise InventoryFetchError(
                    f"Case inventory contains a non-{PROJECT_ID} row for case {row['case_id']}: {row['project_id']}"
                )
            if not row["case_id"] or not row["submitter_id"]:
                raise InventoryFetchError("Case inventory contains a row without case_id or submitter_id.")
            unique_submitters.add(row["submitter_id"])

    if row_count <= 0:
        raise InventoryFetchError("Case inventory validation failed because the row count is zero.")

    return {
        "row_count": row_count,
        "unique_case_submitter_count": len(unique_submitters),
        "required_columns": CASE_TSV_COLUMNS,
    }


def parse_json_array_cell(value: str, column_name: str) -> list[str]:
    """Parse a TSV cell that stores a JSON array."""

    if value is None or str(value).strip() == "":
        return []
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise InventoryFetchError(f"Unable to parse JSON array in column {column_name}: {value}") from exc
    if not isinstance(parsed, list):
        raise InventoryFetchError(f"Expected a JSON array in column {column_name}: {value}")
    return [str(item).strip() for item in parsed if str(item).strip()]


def validate_file_tsv(path: Path, expected_total: int) -> dict[str, Any]:
    """Run post-write validation on the normalized file inventory TSV."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing_columns = [column for column in FILE_TSV_COLUMNS if column not in fieldnames]
        if missing_columns:
            raise InventoryFetchError(f"File inventory is missing required columns: {missing_columns}")

        row_count = 0
        unique_file_ids: set[str] = set()
        for row in reader:
            row_count += 1
            file_id = row["file_id"]
            if not file_id:
                raise InventoryFetchError("File inventory contains a row without file_id.")
            if file_id in unique_file_ids:
                raise InventoryFetchError(f"File inventory contains a duplicate file_id: {file_id}")
            unique_file_ids.add(file_id)

            project_ids = parse_json_array_cell(row["project_ids"], column_name="project_ids")
            if not project_ids:
                raise InventoryFetchError(f"File inventory row {file_id} has no associated project_ids.")
            if any(project_id != PROJECT_ID for project_id in project_ids):
                raise InventoryFetchError(
                    f"File inventory row {file_id} contains non-{PROJECT_ID} project_ids: {project_ids}"
                )

            case_ids = parse_json_array_cell(row["case_ids"], column_name="case_ids")
            if int(row["case_count"]) != len(case_ids):
                raise InventoryFetchError(
                    f"File inventory row {file_id} has case_count={row['case_count']} but {len(case_ids)} case_ids."
                )

    if row_count <= 0:
        raise InventoryFetchError("File inventory validation failed because the row count is zero.")
    if row_count != expected_total:
        raise InventoryFetchError(
            f"File inventory row count {row_count} does not match expected API total {expected_total}."
        )

    return {
        "row_count": row_count,
        "unique_file_id_count": len(unique_file_ids),
        "required_columns": FILE_TSV_COLUMNS,
        "expected_total": expected_total,
    }


def write_latest_manifest(
    latest_manifest_path: Path,
    repo_root: Path,
    run_id: str,
    run_dir: Path,
    case_tsv_path: Path,
    file_tsv_path: Path,
    run_log_path: Path,
    gdc_status: dict[str, Any],
) -> dict[str, Any]:
    """Update the stable pointer file used by the review notebook."""

    latest_manifest = {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "run_id": run_id,
        "gdc_data_release": gdc_status.get("data_release"),
        "gdc_tag": gdc_status.get("tag"),
        "run_directory": repo_relative(run_dir, repo_root),
        "case_inventory_tsv": repo_relative(case_tsv_path, repo_root),
        "file_inventory_tsv": repo_relative(file_tsv_path, repo_root),
        "run_log_json": repo_relative(run_log_path, repo_root),
    }
    write_json(latest_manifest_path, latest_manifest, overwrite=True)
    return latest_manifest


def build_run_log(
    *,
    repo_root: Path,
    run_id: str,
    started_at: datetime,
    completed_at: datetime,
    paths: WorkflowPaths,
    run_dir: Path,
    gdc_status: dict[str, Any],
    gdc_client_config: dict[str, Any],
    project_status_config: dict[str, Any],
    trial_config: dict[str, Any],
    case_fetch_metadata: dict[str, Any],
    file_fetch_metadata: dict[str, Any],
    case_validation: dict[str, Any],
    file_validation: dict[str, Any],
    latest_manifest: dict[str, Any],
) -> dict[str, Any]:
    """Assemble the structured run log written at the end of the workflow."""

    return {
        "status": "completed",
        "run_id": run_id,
        "started_at_utc": format_utc_timestamp(started_at),
        "completed_at_utc": format_utc_timestamp(completed_at),
        "repo_root": str(repo_root.resolve()),
        "trial_name": trial_config.get("trial_name"),
        "project_status": project_status_config.get("status"),
        "gdc_api": {
            "base_url": GDC_API_BASE_URL,
            "status": gdc_status.get("status"),
            "tag": gdc_status.get("tag"),
            "version": gdc_status.get("version"),
            "data_release": gdc_status.get("data_release"),
            "data_release_version": gdc_status.get("data_release_version"),
        },
        "gdc_client_provenance": {
            "configured_path": gdc_client_config.get("gdc_client_path"),
            "configured_version": gdc_client_config.get("gdc_client_version"),
            "configured_project": gdc_client_config.get("gdc_project"),
        },
        "configs_used": {
            "gdc_client_yaml": repo_relative(paths.gdc_client_config, repo_root),
            "project_status_yaml": repo_relative(paths.project_status_config, repo_root),
            "trial_config_yaml": repo_relative(paths.trial_config, repo_root),
            "cases_payload_json": repo_relative(paths.cases_payload_config, repo_root),
            "files_payload_json": repo_relative(paths.files_payload_config, repo_root),
        },
        "outputs": {
            "run_directory": repo_relative(run_dir, repo_root),
            "case_inventory_tsv": repo_relative(run_dir / "tcga_brca_case_inventory.tsv", repo_root),
            "file_inventory_tsv": repo_relative(run_dir / "tcga_brca_file_inventory.tsv", repo_root),
            "latest_manifest_json": repo_relative(paths.latest_manifest, repo_root),
            "gdc_status_json": repo_relative(run_dir / "gdc_status.json", repo_root),
            "cases_payload_effective_json": repo_relative(run_dir / "cases_payload_effective.json", repo_root),
            "files_payload_effective_json": repo_relative(run_dir / "files_payload_effective.json", repo_root),
            "case_response_files": case_fetch_metadata["response_files"],
            "file_response_files": file_fetch_metadata["response_files"],
        },
        "fetch": {
            "cases": case_fetch_metadata,
            "files": file_fetch_metadata,
        },
        "validation": {
            "passed": True,
            "cases": case_validation,
            "files": file_validation,
        },
        "latest_manifest": latest_manifest,
    }


def write_failure_log(
    *,
    path: Path,
    run_id: str,
    started_at: datetime,
    error_message: str,
    repo_root: Path,
) -> None:
    """Persist a failure log when the workflow aborts after creating a run directory."""

    payload = {
        "status": "failed",
        "run_id": run_id,
        "started_at_utc": format_utc_timestamp(started_at),
        "failed_at_utc": format_utc_timestamp(utc_now()),
        "repo_root": str(repo_root.resolve()),
        "error": error_message,
    }
    write_json(path, payload, overwrite=True)


def run_inventory_workflow() -> dict[str, Any]:
    """Execute the full inventory workflow and return the final run summary."""

    started_at = utc_now()
    run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    run_dir = create_run_directory(paths.runs_root, run_id)
    run_log_path = run_dir / "run_log.json"

    try:
        gdc_client_config = load_yaml(paths.gdc_client_config)
        project_status_config = load_yaml(paths.project_status_config)
        trial_config = load_yaml(paths.trial_config)
        cases_payload_config = load_json(paths.cases_payload_config)
        files_payload_config = load_json(paths.files_payload_config)

        session = requests.Session()
        session.headers.update({"Content-Type": "application/json", "Accept": "application/json"})

        gdc_status, gdc_status_text = get_gdc_status(session)
        write_text(run_dir / "gdc_status.json", gdc_status_text)

        case_tsv_path = run_dir / "tcga_brca_case_inventory.tsv"
        file_tsv_path = run_dir / "tcga_brca_file_inventory.tsv"

        case_fetch_metadata = write_inventory_tsv(
            session=session,
            endpoint="cases",
            base_payload=cases_payload_config,
            sort_value="submitter_id:asc",
            payload_output_path=run_dir / "cases_payload_effective.json",
            response_prefix="cases_response",
            tsv_output_path=case_tsv_path,
            fieldnames=CASE_TSV_COLUMNS,
            normalizer=normalize_case_hit,
        )
        file_fetch_metadata = write_inventory_tsv(
            session=session,
            endpoint="files",
            base_payload=files_payload_config,
            sort_value="file_id:asc",
            payload_output_path=run_dir / "files_payload_effective.json",
            response_prefix="files_response",
            tsv_output_path=file_tsv_path,
            fieldnames=FILE_TSV_COLUMNS,
            normalizer=normalize_file_hit,
        )

        case_validation = validate_case_tsv(case_tsv_path)
        file_validation = validate_file_tsv(file_tsv_path, expected_total=file_fetch_metadata["total_expected"])

        latest_manifest = write_latest_manifest(
            latest_manifest_path=paths.latest_manifest,
            repo_root=paths.repo_root,
            run_id=run_id,
            run_dir=run_dir,
            case_tsv_path=case_tsv_path,
            file_tsv_path=file_tsv_path,
            run_log_path=run_log_path,
            gdc_status=gdc_status,
        )

        completed_at = utc_now()
        run_log = build_run_log(
            repo_root=paths.repo_root,
            run_id=run_id,
            started_at=started_at,
            completed_at=completed_at,
            paths=paths,
            run_dir=run_dir,
            gdc_status=gdc_status,
            gdc_client_config=gdc_client_config,
            project_status_config=project_status_config,
            trial_config=trial_config,
            case_fetch_metadata=case_fetch_metadata,
            file_fetch_metadata=file_fetch_metadata,
            case_validation=case_validation,
            file_validation=file_validation,
            latest_manifest=latest_manifest,
        )
        write_json(run_log_path, run_log)
        return run_log

    except Exception as exc:
        write_failure_log(
            path=run_log_path,
            run_id=run_id,
            started_at=started_at,
            error_message=str(exc),
            repo_root=paths.repo_root,
        )
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    """Print a concise human-readable summary for researchers."""

    print("TCGA-BRCA inventory fetch complete.")
    print(f"Run ID: {run_log['run_id']}")
    print(f"GDC data release: {run_log['gdc_api'].get('data_release')}")
    print(f"Case count: {run_log['validation']['cases']['row_count']}")
    print(f"Unique case submitters: {run_log['validation']['cases']['unique_case_submitter_count']}")
    print(f"File count: {run_log['validation']['files']['row_count']}")
    print(f"Validation passed: {run_log['validation']['passed']}")
    print(f"Latest manifest: {run_log['outputs']['latest_manifest_json']}")


def main() -> int:
    """Run the inventory workflow from the command line."""

    run_log = run_inventory_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
