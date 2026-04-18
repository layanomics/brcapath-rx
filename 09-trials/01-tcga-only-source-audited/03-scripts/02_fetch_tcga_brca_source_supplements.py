#!/usr/bin/env python
"""Fetch TCGA-BRCA Clinical and Biospecimen Supplement source manifests."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover
    yaml = None

GDC_API_BASE_URL = "https://api.gdc.cancer.gov"
PROJECT_ID = "TCGA-BRCA"
PAGE_SIZE = 2000
REQUEST_TIMEOUT_SECONDS = 120
MANIFEST_COLUMNS = ["id", "filename", "md5", "size", "state"]
METADATA_TSV_COLUMNS = [
    "source_class",
    "file_id",
    "file_name",
    "md5sum",
    "file_size",
    "access",
    "state",
    "data_category",
    "data_type",
    "data_format",
    "case_ids",
    "case_submitter_ids",
    "project_ids",
    "case_count",
]


class SupplementFetchError(RuntimeError):
    """Raised when the supplement acquisition workflow cannot complete safely."""


@dataclass(frozen=True)
class SourceClassConfig:
    """Configuration for one supplement source class."""

    source_class: str
    data_type: str
    data_category: str
    payload_config: Path
    manifest_dir: Path
    download_dir: Path
    log_dir: Path
    manifest_filename: str
    download_log_filename: str = "gdc_client_download.log"


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the supplement workflow."""

    repo_root: Path
    gdc_client_config: Path
    project_status_config: Path
    trial_config: Path
    download_scope_config: Path
    source_root: Path
    supplement_audit_root: Path
    supplement_runs_root: Path
    latest_manifest: Path
    source_classes: dict[str, SourceClassConfig]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise SupplementFetchError("Unable to locate the repository root from the script path.")


def build_workflow_paths() -> WorkflowPaths:
    repo_root = detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    raw_root = repo_root / "01-data" / "raw" / "tcga-brca" / "gdc"
    source_root = repo_root / "01-data" / "audit" / "tcga-brca" / "source"

    clinical = SourceClassConfig(
        source_class="clinical",
        data_type="Clinical Supplement",
        data_category="Clinical",
        payload_config=trial_root / "04-config" / "gdc_clinical_supplement_payload.json",
        manifest_dir=raw_root / "manifests" / "clinical",
        download_dir=raw_root / "downloads" / "clinical",
        log_dir=raw_root / "logs" / "clinical",
        manifest_filename="clinical_supplement_manifest.tsv",
    )
    biospecimen = SourceClassConfig(
        source_class="biospecimen",
        data_type="Biospecimen Supplement",
        data_category="Biospecimen",
        payload_config=trial_root / "04-config" / "gdc_biospecimen_supplement_payload.json",
        manifest_dir=raw_root / "manifests" / "biospecimen",
        download_dir=raw_root / "downloads" / "biospecimen",
        log_dir=raw_root / "logs" / "biospecimen",
        manifest_filename="biospecimen_supplement_manifest.tsv",
    )

    return WorkflowPaths(
        repo_root=repo_root,
        gdc_client_config=repo_root / "07-config" / "gdc_client.yaml",
        project_status_config=repo_root / "07-config" / "project_status.yaml",
        trial_config=trial_root / "04-config" / "trial_config.yaml",
        download_scope_config=trial_root / "04-config" / "gdc_download_scope.yaml",
        source_root=source_root,
        supplement_audit_root=source_root / "supplements",
        supplement_runs_root=source_root / "supplements" / "runs",
        latest_manifest=source_root / "tcga_brca_source_supplements_latest.json",
        source_classes={"clinical": clinical, "biospecimen": biospecimen},
    )


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
                raise SupplementFetchError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise SupplementFetchError(
                        f"Cannot mix list and scalar/mapping values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise SupplementFetchError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise SupplementFetchError(
                    f"Cannot mix mapping and scalar/list values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise SupplementFetchError(f"Expected a top-level key/value pair in YAML: {source_path}:{line_number}")
        key, value = raw_line.split(":", 1)
        key = key.strip()
        value = value.strip()
        current_key = key
        if value == "":
            parsed[key] = {}
        else:
            parsed[key] = parse_yaml_scalar(value)

    return parsed


def load_yaml(path: Path) -> dict[str, Any]:
    raw_text = path.read_text(encoding="utf-8")
    if yaml is not None:
        data = yaml.safe_load(raw_text)
    else:
        data = parse_simple_yaml_mapping(raw_text, path)
    if not isinstance(data, dict):
        raise SupplementFetchError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SupplementFetchError(f"Expected a JSON object in file: {path}")
    return data


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise SupplementFetchError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_text(path: Path, text: str, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise SupplementFetchError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def build_metadata_payload(base_payload: dict[str, Any], offset: int) -> dict[str, Any]:
    payload = deepcopy(base_payload)
    payload["format"] = "JSON"
    payload["size"] = str(PAGE_SIZE)
    payload["from"] = offset
    payload["sort"] = "file_id:asc"
    return payload


def build_manifest_payload(base_payload: dict[str, Any]) -> dict[str, Any]:
    payload = deepcopy(base_payload)
    payload["sort"] = "file_id:asc"
    return payload


def create_run_directory(runs_root: Path, run_id: str) -> Path:
    runs_root.mkdir(parents=True, exist_ok=True)
    run_dir = runs_root / run_id
    if run_dir.exists():
        raise SupplementFetchError(f"Run directory already exists: {run_dir}")
    run_dir.mkdir(parents=False, exist_ok=False)
    return run_dir


def get_gdc_status(session: requests.Session) -> tuple[dict[str, Any], str]:
    response = session.get(f"{GDC_API_BASE_URL}/status", timeout=REQUEST_TIMEOUT_SECONDS)
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        raise SupplementFetchError(f"GDC status request failed with HTTP {response.status_code}.") from exc
    try:
        return response.json(), response.text
    except json.JSONDecodeError as exc:
        raise SupplementFetchError("Unable to decode the GDC status response as JSON.") from exc


def post_json(session: requests.Session, endpoint: str, payload: dict[str, Any]) -> tuple[dict[str, Any], str]:
    response = session.post(f"{GDC_API_BASE_URL}/{endpoint}", json=payload, timeout=REQUEST_TIMEOUT_SECONDS)
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        body_preview = response.text[:500]
        raise SupplementFetchError(
            f"GDC {endpoint} request failed with HTTP {response.status_code}. Body preview: {body_preview}"
        ) from exc
    try:
        return response.json(), response.text
    except json.JSONDecodeError as exc:
        raise SupplementFetchError(f"Unable to decode the GDC {endpoint} response as JSON.") from exc


def post_manifest(session: requests.Session, payload: dict[str, Any]) -> str:
    response = session.post(
        f"{GDC_API_BASE_URL}/files?return_type=manifest",
        json=payload,
        timeout=REQUEST_TIMEOUT_SECONDS,
    )
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        body_preview = response.text[:500]
        raise SupplementFetchError(
            f"GDC manifest request failed with HTTP {response.status_code}. Body preview: {body_preview}"
        ) from exc
    return response.text


def validate_response_page(
    *, source_class: str, offset: int, response_json: dict[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    data = response_json.get("data")
    if not isinstance(data, dict):
        raise SupplementFetchError(f"{source_class} response is missing the top-level data object.")

    hits = data.get("hits")
    pagination = data.get("pagination")
    if not isinstance(hits, list) or not isinstance(pagination, dict):
        raise SupplementFetchError(f"{source_class} response is missing hits or pagination.")

    try:
        total = int(pagination["total"])
        count = int(pagination["count"])
        current_from = int(pagination["from"])
        pages = int(pagination["pages"])
        page = int(pagination["page"])
    except (KeyError, TypeError, ValueError) as exc:
        raise SupplementFetchError(f"{source_class} pagination metadata is incomplete: {pagination}") from exc

    if current_from != offset:
        raise SupplementFetchError(
            f"{source_class} response returned pagination.from={current_from} for requested offset {offset}."
        )
    if count != len(hits):
        raise SupplementFetchError(
            f"{source_class} pagination.count={count} does not match returned hits {len(hits)}."
        )
    if count > PAGE_SIZE:
        raise SupplementFetchError(f"{source_class} response exceeded page size {PAGE_SIZE}.")

    return hits, {"total": total, "count": count, "from": current_from, "page": page, "pages": pages}


def unique_preserving_order(values: list[str]) -> list[str]:
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
    return json.dumps(values, ensure_ascii=True, separators=(",", ":"))


def normalize_hit(hit: dict[str, Any], source_class_config: SourceClassConfig) -> dict[str, str]:
    case_ids: list[str] = []
    case_submitter_ids: list[str] = []
    project_ids: list[str] = []
    for case in hit.get("cases") or []:
        case_ids.extend(unique_preserving_order([str(case.get("case_id") or "").strip()]))
        case_submitter_ids.extend(unique_preserving_order([str(case.get("submitter_id") or "").strip()]))
        project_id = str((case.get("project") or {}).get("project_id") or "").strip()
        project_ids.extend(unique_preserving_order([project_id]))

    unique_case_ids = unique_preserving_order(case_ids)
    unique_case_submitters = unique_preserving_order(case_submitter_ids)
    unique_project_ids = unique_preserving_order(project_ids)

    return {
        "source_class": source_class_config.source_class,
        "file_id": str(hit.get("file_id") or hit.get("id") or "").strip(),
        "file_name": str(hit.get("file_name") or "").strip(),
        "md5sum": str(hit.get("md5sum") or "").strip(),
        "file_size": str(hit.get("file_size") or "").strip(),
        "access": str(hit.get("access") or "").strip(),
        "state": str(hit.get("state") or "").strip(),
        "data_category": str(hit.get("data_category") or "").strip(),
        "data_type": str(hit.get("data_type") or "").strip(),
        "data_format": str(hit.get("data_format") or "").strip(),
        "case_ids": to_json_array_cell(unique_case_ids),
        "case_submitter_ids": to_json_array_cell(unique_case_submitters),
        "project_ids": to_json_array_cell(unique_project_ids),
        "case_count": str(len(unique_case_ids)),
    }


def fetch_metadata_for_source_class(
    *,
    session: requests.Session,
    source_class_config: SourceClassConfig,
    base_payload: dict[str, Any],
    audit_run_dir: Path,
) -> dict[str, Any]:
    effective_payload_path = audit_run_dir / f"{source_class_config.source_class}_metadata_payload_effective.json"
    response_prefix = f"{source_class_config.source_class}_response"
    first_payload = build_metadata_payload(base_payload, offset=0)
    write_json(effective_payload_path, first_payload)

    rows: list[dict[str, str]] = []
    response_files: list[str] = []
    offset = 0
    page_index = 1
    total_expected: int | None = None
    pages_expected: int | None = None

    while True:
        payload = build_metadata_payload(base_payload, offset=offset)
        response_json, response_text = post_json(session, "files", payload)
        page_path = audit_run_dir / f"{response_prefix}_page_{page_index:04d}.json"
        write_text(page_path, response_text)
        response_files.append(page_path.name)

        hits, pagination = validate_response_page(
            source_class=source_class_config.source_class,
            offset=offset,
            response_json=response_json,
        )
        if total_expected is None:
            total_expected = pagination["total"]
            pages_expected = pagination["pages"]
        elif total_expected != pagination["total"]:
            raise SupplementFetchError(
                f"{source_class_config.source_class} pagination.total changed during paging "
                f"({total_expected} -> {pagination['total']})."
            )

        rows.extend(normalize_hit(hit, source_class_config) for hit in hits)
        if len(rows) == total_expected:
            break
        if pagination["count"] == 0:
            raise SupplementFetchError(
                f"{source_class_config.source_class} returned an empty page before expected total {total_expected}."
            )
        offset += pagination["count"]
        page_index += 1

    if total_expected is None:
        raise SupplementFetchError(f"{source_class_config.source_class} did not return pagination metadata.")
    if len(rows) != total_expected:
        raise SupplementFetchError(
            f"{source_class_config.source_class} normalized row count {len(rows)} "
            f"does not match expected total {total_expected}."
        )

    return {
        "rows": rows,
        "rows_written": len(rows),
        "total_expected": total_expected,
        "pages_written": page_index,
        "pages_expected": pages_expected,
        "response_files": response_files,
        "payload_path": effective_payload_path,
    }


def write_manifest_for_source_class(
    *,
    session: requests.Session,
    source_class_config: SourceClassConfig,
    base_payload: dict[str, Any],
    audit_run_dir: Path,
    run_id: str,
) -> dict[str, Any]:
    manifest_payload = build_manifest_payload(base_payload)
    manifest_payload_path = audit_run_dir / f"{source_class_config.source_class}_manifest_payload_effective.json"
    write_json(manifest_payload_path, manifest_payload)

    manifest_text = post_manifest(session, manifest_payload)
    manifest_run_dir = source_class_config.manifest_dir / run_id
    if manifest_run_dir.exists():
        raise SupplementFetchError(f"Manifest run directory already exists: {manifest_run_dir}")
    manifest_run_dir.mkdir(parents=True, exist_ok=False)

    manifest_path = manifest_run_dir / source_class_config.manifest_filename
    write_text(manifest_path, manifest_text)

    return {
        "manifest_payload_path": manifest_payload_path,
        "manifest_path": manifest_path,
    }


def read_manifest_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if fieldnames != MANIFEST_COLUMNS:
            raise SupplementFetchError(f"Manifest columns do not match expected schema in {path}: {fieldnames}")
        return list(reader)


def write_combined_metadata_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if path.exists():
        raise SupplementFetchError(f"Refusing to overwrite existing metadata TSV: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_TSV_COLUMNS, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def parse_json_array_cell(value: str, column_name: str) -> list[str]:
    if value is None or str(value).strip() == "":
        return []
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise SupplementFetchError(f"Unable to parse JSON array in column {column_name}: {value}") from exc
    if not isinstance(parsed, list):
        raise SupplementFetchError(f"Expected a JSON array in column {column_name}: {value}")
    return [str(item).strip() for item in parsed if str(item).strip()]


def validate_metadata_rows(rows: list[dict[str, str]], source_classes: dict[str, SourceClassConfig]) -> dict[str, Any]:
    if len(rows) <= 0:
        raise SupplementFetchError("Supplement metadata validation failed because the row count is zero.")

    observed_source_classes = {row["source_class"] for row in rows}
    expected_source_classes = set(source_classes.keys())
    if observed_source_classes != expected_source_classes:
        raise SupplementFetchError(
            f"Observed source classes {sorted(observed_source_classes)} do not match expected {sorted(expected_source_classes)}."
        )

    counts_by_source_class: dict[str, int] = {key: 0 for key in source_classes}
    for row in rows:
        source_class = row["source_class"]
        config = source_classes[source_class]
        if row["data_type"] != config.data_type:
            raise SupplementFetchError(
                f"{source_class} row {row['file_id']} has data_type={row['data_type']} but expected {config.data_type}."
            )
        if row["data_category"] != config.data_category:
            raise SupplementFetchError(
                f"{source_class} row {row['file_id']} has data_category={row['data_category']} "
                f"but expected {config.data_category}."
            )
        project_ids = parse_json_array_cell(row["project_ids"], "project_ids")
        if not project_ids:
            raise SupplementFetchError(f"{source_class} row {row['file_id']} has no project_ids.")
        if any(project_id != PROJECT_ID for project_id in project_ids):
            raise SupplementFetchError(
                f"{source_class} row {row['file_id']} contains non-{PROJECT_ID} project_ids: {project_ids}"
            )
        case_ids = parse_json_array_cell(row["case_ids"], "case_ids")
        if int(row["case_count"]) != len(case_ids):
            raise SupplementFetchError(
                f"{source_class} row {row['file_id']} has case_count={row['case_count']} but {len(case_ids)} case_ids."
            )
        counts_by_source_class[source_class] += 1

    return {
        "row_count": len(rows),
        "counts_by_source_class": counts_by_source_class,
        "required_columns": METADATA_TSV_COLUMNS,
    }


def validate_manifest_against_rows(
    *,
    source_class_config: SourceClassConfig,
    manifest_rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
) -> dict[str, Any]:
    if len(manifest_rows) <= 0:
        raise SupplementFetchError(f"{source_class_config.source_class} manifest row count is zero.")
    if len(manifest_rows) != len(metadata_rows):
        raise SupplementFetchError(
            f"{source_class_config.source_class} manifest row count {len(manifest_rows)} "
            f"does not match metadata row count {len(metadata_rows)}."
        )

    metadata_by_file_id = {row["file_id"]: row for row in metadata_rows}
    for manifest_row in manifest_rows:
        file_id = manifest_row["id"]
        if file_id not in metadata_by_file_id:
            raise SupplementFetchError(
                f"{source_class_config.source_class} manifest contains file_id not found in metadata TSV: {file_id}"
            )
        metadata_row = metadata_by_file_id[file_id]
        if manifest_row["filename"] != metadata_row["file_name"]:
            raise SupplementFetchError(
                f"{source_class_config.source_class} filename mismatch for {file_id}: "
                f"{manifest_row['filename']} vs {metadata_row['file_name']}"
            )
        if manifest_row["md5"] != metadata_row["md5sum"]:
            raise SupplementFetchError(
                f"{source_class_config.source_class} md5 mismatch for {file_id}: "
                f"{manifest_row['md5']} vs {metadata_row['md5sum']}"
            )
        if manifest_row["size"] != metadata_row["file_size"]:
            raise SupplementFetchError(
                f"{source_class_config.source_class} size mismatch for {file_id}: "
                f"{manifest_row['size']} vs {metadata_row['file_size']}"
            )
        if manifest_row["state"] != metadata_row["state"]:
            raise SupplementFetchError(
                f"{source_class_config.source_class} state mismatch for {file_id}: "
                f"{manifest_row['state']} vs {metadata_row['state']}"
            )

    return {"row_count": len(manifest_rows), "columns": MANIFEST_COLUMNS}


def run_gdc_client_download(
    *,
    source_class_config: SourceClassConfig,
    gdc_client_path: str,
    manifest_path: Path,
    download_dir: Path,
    log_path: Path,
    token_file: str | None,
) -> dict[str, Any]:
    if download_dir.exists():
        raise SupplementFetchError(f"Download directory already exists: {download_dir}")
    download_dir.mkdir(parents=True, exist_ok=False)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    command = [
        gdc_client_path,
        "download",
        "--manifest",
        str(manifest_path),
        "--dir",
        str(download_dir),
        "--log-file",
        str(log_path),
    ]
    if token_file:
        command.extend(["--token-file", token_file])

    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise SupplementFetchError(
            f"gdc-client download failed for {source_class_config.source_class} with exit code {result.returncode}. "
            f"stderr: {result.stderr[:500]}"
        )

    return {
        "command": command,
        "return_code": result.returncode,
        "stdout_preview": result.stdout[:500],
        "stderr_preview": result.stderr[:500],
    }


def validate_download_outputs(manifest_rows: list[dict[str, str]], download_dir: Path, source_class: str) -> dict[str, Any]:
    missing_files: list[str] = []
    for row in manifest_rows:
        expected_path = download_dir / row["id"] / row["filename"]
        if not expected_path.exists():
            missing_files.append(str(expected_path))
    if missing_files:
        preview = missing_files[:5]
        raise SupplementFetchError(
            f"{source_class} download validation failed; missing {len(missing_files)} files. Examples: {preview}"
        )
    return {"validated_downloaded_file_count": len(manifest_rows)}


def write_failure_log(path: Path, run_id: str, started_at: datetime, repo_root: Path, error_message: str) -> None:
    payload = {
        "status": "failed",
        "run_id": run_id,
        "started_at_utc": format_utc_timestamp(started_at),
        "failed_at_utc": format_utc_timestamp(utc_now()),
        "repo_root": str(repo_root.resolve()),
        "error": error_message,
    }
    write_json(path, payload, overwrite=True)


def write_latest_manifest(
    *,
    latest_manifest_path: Path,
    repo_root: Path,
    run_id: str,
    run_dir: Path,
    metadata_tsv_path: Path,
    source_class_results: dict[str, dict[str, Any]],
    run_log_path: Path,
    gdc_status: dict[str, Any],
) -> dict[str, Any]:
    latest_manifest = {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "run_id": run_id,
        "gdc_data_release": gdc_status.get("data_release"),
        "gdc_tag": gdc_status.get("tag"),
        "audit_run_directory": repo_relative(run_dir, repo_root),
        "metadata_tsv": repo_relative(metadata_tsv_path, repo_root),
        "run_log_json": repo_relative(run_log_path, repo_root),
        "source_classes": {},
    }
    for source_class, result in source_class_results.items():
        latest_manifest["source_classes"][source_class] = {
            "manifest_path": repo_relative(result["manifest_path"], repo_root),
            "download_dir": repo_relative(result["download_dir"], repo_root) if result["download_dir"] else None,
            "download_log_path": repo_relative(result["download_log_path"], repo_root)
            if result["download_log_path"]
            else None,
            "download_requested": result["download_requested"],
            "download_completed": result["download_completed"],
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
    gdc_status: dict[str, Any],
    gdc_client_config: dict[str, Any],
    project_status_config: dict[str, Any],
    trial_config: dict[str, Any],
    download_scope_config: dict[str, Any],
    run_dir: Path,
    metadata_tsv_path: Path,
    validation: dict[str, Any],
    source_class_results: dict[str, dict[str, Any]],
    latest_manifest: dict[str, Any],
) -> dict[str, Any]:
    outputs = {
        "audit_run_directory": repo_relative(run_dir, repo_root),
        "metadata_tsv": repo_relative(metadata_tsv_path, repo_root),
        "latest_manifest_json": repo_relative(paths.latest_manifest, repo_root),
    }
    for source_class, result in source_class_results.items():
        outputs[f"{source_class}_manifest_path"] = repo_relative(result["manifest_path"], repo_root)
        outputs[f"{source_class}_download_dir"] = (
            repo_relative(result["download_dir"], repo_root) if result["download_dir"] else None
        )
        outputs[f"{source_class}_download_log_path"] = (
            repo_relative(result["download_log_path"], repo_root) if result["download_log_path"] else None
        )

    serialized_source_class_results: dict[str, dict[str, Any]] = {}
    for source_class, result in source_class_results.items():
        serialized_source_class_results[source_class] = {
            **result,
            "manifest_path": repo_relative(result["manifest_path"], repo_root),
            "download_dir": repo_relative(result["download_dir"], repo_root) if result["download_dir"] else None,
            "download_log_path": (
                repo_relative(result["download_log_path"], repo_root) if result["download_log_path"] else None
            ),
        }

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
        "download_scope_config": download_scope_config,
        "configs_used": {
            "gdc_client_yaml": repo_relative(paths.gdc_client_config, repo_root),
            "project_status_yaml": repo_relative(paths.project_status_config, repo_root),
            "trial_config_yaml": repo_relative(paths.trial_config, repo_root),
            "gdc_download_scope_yaml": repo_relative(paths.download_scope_config, repo_root),
            "clinical_payload_json": repo_relative(paths.source_classes["clinical"].payload_config, repo_root),
            "biospecimen_payload_json": repo_relative(paths.source_classes["biospecimen"].payload_config, repo_root),
        },
        "outputs": outputs,
        "source_classes": serialized_source_class_results,
        "validation": validation,
        "latest_manifest": latest_manifest,
    }


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA supplement acquisition complete.")
    print(f"Run ID: {run_log['run_id']}")
    print(f"GDC data release: {run_log['gdc_api'].get('data_release')}")
    print(
        "Clinical files: "
        f"{run_log['validation']['metadata']['counts_by_source_class']['clinical']}"
    )
    print(
        "Biospecimen files: "
        f"{run_log['validation']['metadata']['counts_by_source_class']['biospecimen']}"
    )
    print(f"Metadata TSV: {run_log['outputs']['metadata_tsv']}")
    print(f"Clinical manifest: {run_log['outputs']['clinical_manifest_path']}")
    print(f"Biospecimen manifest: {run_log['outputs']['biospecimen_manifest_path']}")
    download_completed = all(
        result["download_completed"] for result in run_log["source_classes"].values()
    )
    print(f"Download completed: {download_completed}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch TCGA-BRCA Clinical and Biospecimen Supplement manifests and metadata."
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="After manifest generation, download the raw supplement files with gdc-client.",
    )
    parser.add_argument(
        "--token-file",
        help="Optional GDC token file path passed through to gdc-client when --download is used.",
    )
    return parser.parse_args()


def run_workflow(download_requested: bool, token_file: str | None) -> dict[str, Any]:
    started_at = utc_now()
    run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    run_dir = create_run_directory(paths.supplement_runs_root, run_id)
    run_log_path = run_dir / "run_log.json"

    try:
        gdc_client_config = load_yaml(paths.gdc_client_config)
        project_status_config = load_yaml(paths.project_status_config)
        trial_config = load_yaml(paths.trial_config)
        download_scope_config = load_yaml(paths.download_scope_config)

        session = requests.Session()
        session.headers.update({"Content-Type": "application/json", "Accept": "application/json"})

        gdc_status, gdc_status_text = get_gdc_status(session)
        write_text(run_dir / "gdc_status.json", gdc_status_text)

        all_rows: list[dict[str, str]] = []
        source_class_results: dict[str, dict[str, Any]] = {}

        for source_class, source_class_config in paths.source_classes.items():
            base_payload = load_json(source_class_config.payload_config)

            metadata_result = fetch_metadata_for_source_class(
                session=session,
                source_class_config=source_class_config,
                base_payload=base_payload,
                audit_run_dir=run_dir,
            )
            manifest_result = write_manifest_for_source_class(
                session=session,
                source_class_config=source_class_config,
                base_payload=base_payload,
                audit_run_dir=run_dir,
                run_id=run_id,
            )
            manifest_rows = read_manifest_rows(manifest_result["manifest_path"])
            metadata_rows = metadata_result["rows"]
            manifest_validation = validate_manifest_against_rows(
                source_class_config=source_class_config,
                manifest_rows=manifest_rows,
                metadata_rows=metadata_rows,
            )

            download_dir = None
            download_log_path = None
            download_validation = None
            gdc_client_result = None
            download_completed = False
            if download_requested:
                download_dir = source_class_config.download_dir / run_id
                download_log_path = source_class_config.log_dir / run_id / source_class_config.download_log_filename
                gdc_client_result = run_gdc_client_download(
                    source_class_config=source_class_config,
                    gdc_client_path=str(gdc_client_config.get("gdc_client_path")),
                    manifest_path=manifest_result["manifest_path"],
                    download_dir=download_dir,
                    log_path=download_log_path,
                    token_file=token_file,
                )
                download_validation = validate_download_outputs(
                    manifest_rows=manifest_rows,
                    download_dir=download_dir,
                    source_class=source_class,
                )
                download_completed = True

            all_rows.extend(metadata_rows)
            source_class_results[source_class] = {
                "source_class": source_class,
                "data_type": source_class_config.data_type,
                "data_category": source_class_config.data_category,
                "count": len(metadata_rows),
                "metadata_payload_path": repo_relative(metadata_result["payload_path"], paths.repo_root),
                "response_files": metadata_result["response_files"],
                "manifest_payload_path": repo_relative(manifest_result["manifest_payload_path"], paths.repo_root),
                "manifest_path": manifest_result["manifest_path"],
                "manifest_validation": manifest_validation,
                "download_requested": download_requested,
                "download_completed": download_completed,
                "download_dir": download_dir,
                "download_log_path": download_log_path,
                "download_validation": download_validation,
                "gdc_client_result": gdc_client_result,
            }

        metadata_tsv_path = run_dir / "tcga_brca_source_supplements_metadata.tsv"
        write_combined_metadata_tsv(metadata_tsv_path, all_rows)
        metadata_validation = validate_metadata_rows(all_rows, paths.source_classes)

        validation = {
            "passed": True,
            "metadata": metadata_validation,
            "manifests": {
                source_class: result["manifest_validation"] for source_class, result in source_class_results.items()
            },
            "downloads": {
                source_class: result["download_validation"]
                for source_class, result in source_class_results.items()
                if result["download_requested"]
            },
        }

        latest_manifest = write_latest_manifest(
            latest_manifest_path=paths.latest_manifest,
            repo_root=paths.repo_root,
            run_id=run_id,
            run_dir=run_dir,
            metadata_tsv_path=metadata_tsv_path,
            source_class_results=source_class_results,
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
            gdc_status=gdc_status,
            gdc_client_config=gdc_client_config,
            project_status_config=project_status_config,
            trial_config=trial_config,
            download_scope_config=download_scope_config,
            run_dir=run_dir,
            metadata_tsv_path=metadata_tsv_path,
            validation=validation,
            source_class_results=source_class_results,
            latest_manifest=latest_manifest,
        )
        write_json(run_log_path, run_log)
        return run_log

    except Exception as exc:
        write_failure_log(
            path=run_log_path,
            run_id=run_id,
            started_at=started_at,
            repo_root=paths.repo_root,
            error_message=str(exc),
        )
        raise


def main() -> int:
    args = parse_args()
    run_log = run_workflow(download_requested=args.download, token_file=args.token_file)
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
