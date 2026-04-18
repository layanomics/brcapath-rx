#!/usr/bin/env python
"""Build an auditable TCGA-BRCA cohort-construction blueprint from saved audit layers."""

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

FIELD_BUCKET_ORDER = (
    "required_for_first_baseline",
    "optional_for_first_baseline",
    "deferred_for_later_expansion",
)
SOURCE_LAYER_ORDER = ("clinical", "endpoint", "biospecimen")
JOIN_STRENGTHS = (
    "strong_current_evidence",
    "workable_but_needs_review",
    "deferred",
)
LINK_TIMINGS = ("link_now", "link_later", "prepare_only")
SEVERITY_ORDER = ("high", "medium", "low")
NOTES_PLACEHOLDER = "[fill in during cohort blueprint review]"
FIELD_TABLE_FIELDNAMES = [
    "blueprint_run_id",
    "field_bucket",
    "source_layer",
    "table_name",
    "field_name",
    "proposed_role",
    "rationale",
    "dependency_for_join",
    "ambiguity_flag",
    "notes_placeholder",
    "upstream_source_type",
    "upstream_bucket_or_role",
    "upstream_rule",
    "manual_review_priority",
    "missing_like_fraction",
    "non_missing_count",
    "distinct_non_missing_count",
    "paired_field_name",
    "evidence_rule",
    "evidence_rate",
]
JOIN_PATH_FIELDNAMES = [
    "blueprint_run_id",
    "step_order",
    "left_node",
    "right_node",
    "left_table_name",
    "right_table_name",
    "left_field_name",
    "right_field_name",
    "join_keys_or_evidence",
    "join_strength",
    "link_timing",
    "rationale",
    "ambiguity_flag",
    "notes_placeholder",
]
AMBIGUITY_FIELDNAMES = [
    "blueprint_run_id",
    "ambiguity_id",
    "ambiguity_category",
    "source_layer",
    "table_name",
    "field_name_or_step",
    "severity",
    "why_unresolved",
    "current_saved_evidence",
    "proposed_blueprint_handling",
    "blocks_final_cohort_construction",
    "notes_placeholder",
]
SUMMARY_FIELDNAMES = [
    "blueprint_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]


class CohortBlueprintError(RuntimeError):
    """Raised when the cohort blueprint workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the cohort blueprint workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    clinical_shortlist_latest_pointer: Path
    endpoint_crosswalk_latest_pointer: Path
    biospecimen_crosswalk_latest_pointer: Path
    blueprint_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved blueprint inputs loaded from saved audit layers."""

    clinical_shortlist_latest_pointer: dict[str, Any]
    clinical_shortlist_run_log: dict[str, Any]
    endpoint_crosswalk_latest_pointer: dict[str, Any]
    endpoint_crosswalk_run_log: dict[str, Any]
    biospecimen_crosswalk_latest_pointer: dict[str, Any]
    biospecimen_crosswalk_run_log: dict[str, Any]
    core_audit_latest_pointer: dict[str, Any]
    core_audit_run_log: dict[str, Any]
    clinical_biotab_run_log: dict[str, Any]
    source_run_log: dict[str, Any]
    clinical_shortlist_rows: list[dict[str, str]]
    clinical_shortlist_summary_rows: list[dict[str, str]]
    clinical_shortlist_by_table_rows: list[dict[str, str]]
    clinical_core_field_audit_rows: list[dict[str, str]]
    endpoint_inventory_rows: list[dict[str, str]]
    endpoint_crosswalk_rows: list[dict[str, str]]
    endpoint_summary_rows: list[dict[str, str]]
    biospecimen_inventory_rows: list[dict[str, str]]
    biospecimen_crosswalk_rows: list[dict[str, str]]
    biospecimen_summary_rows: list[dict[str, str]]
    source_metadata_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise CohortBlueprintError("Unable to locate the repository root from the script path.")


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
                raise CohortBlueprintError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise CohortBlueprintError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise CohortBlueprintError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise CohortBlueprintError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise CohortBlueprintError(
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
        raise CohortBlueprintError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise CohortBlueprintError(f"Expected a JSON object in file: {path}")
    return data


def read_tsv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise CohortBlueprintError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise CohortBlueprintError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise CohortBlueprintError(f"Run directory already exists: {path}")
    path.mkdir(parents=False, exist_ok=False)
    return path


def parse_int(value: Any, label: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise CohortBlueprintError(f"Expected integer-like value for {label}: {value!r}") from exc


def require_keys(payload: dict[str, Any], required_keys: set[str], label: str, source_path: Path) -> None:
    missing_keys = required_keys.difference(payload.keys())
    if missing_keys:
        raise CohortBlueprintError(
            f"{label} is missing required keys {sorted(missing_keys)}: {source_path}"
        )


def resolve_existing_path(repo_root: Path, relative_path: str, label: str) -> Path:
    path = repo_root / relative_path
    if not path.exists():
        raise CohortBlueprintError(f"Required {label} not found: {path}")
    return path


def unique_case_submitter_count(source_metadata_rows: list[dict[str, str]], source_class: str) -> int:
    observed: set[str] = set()
    for row in source_metadata_rows:
        if str(row.get("source_class") or "") != source_class:
            continue
        raw_case_ids = str(row.get("case_submitter_ids") or "[]")
        case_ids = json.loads(raw_case_ids)
        if not isinstance(case_ids, list):
            raise CohortBlueprintError(
                f"Expected case_submitter_ids JSON list in source metadata row: {raw_case_ids!r}"
            )
        for case_id in case_ids:
            if case_id:
                observed.add(str(case_id))
    if not observed:
        raise CohortBlueprintError(
            f"No case_submitter_ids were recovered for source_class={source_class!r} from source metadata."
        )
    return len(observed)


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

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        clinical_shortlist_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_shortlist_latest.json"
        ),
        endpoint_crosswalk_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_endpoint_crosswalk_latest.json"
        ),
        biospecimen_crosswalk_latest_pointer=(
            audit_root
            / "tcga-brca"
            / "variables"
            / "tcga_brca_biospecimen_identifier_crosswalk_latest.json"
        ),
        blueprint_runs_root=cohort_root / "cohort_blueprint_runs",
        latest_pointer=cohort_root / "tcga_brca_cohort_blueprint_latest.json",
    )


def load_workflow_inputs(paths: WorkflowPaths) -> WorkflowInputs:
    for pointer_path, label in [
        (paths.clinical_shortlist_latest_pointer, "clinical shortlist latest pointer"),
        (paths.endpoint_crosswalk_latest_pointer, "endpoint crosswalk latest pointer"),
        (paths.biospecimen_crosswalk_latest_pointer, "biospecimen crosswalk latest pointer"),
    ]:
        if not pointer_path.exists():
            raise CohortBlueprintError(f"Required {label} not found: {pointer_path}")

    clinical_shortlist_latest_pointer = load_json(paths.clinical_shortlist_latest_pointer)
    endpoint_crosswalk_latest_pointer = load_json(paths.endpoint_crosswalk_latest_pointer)
    biospecimen_crosswalk_latest_pointer = load_json(paths.biospecimen_crosswalk_latest_pointer)

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
            "core_audit_latest_json",
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

    if str(endpoint_crosswalk_latest_pointer["shortlist_run_id"]) != str(
        clinical_shortlist_latest_pointer["shortlist_run_id"]
    ):
        raise CohortBlueprintError(
            "Endpoint crosswalk latest pointer does not reference the current clinical shortlist run: "
            f"{endpoint_crosswalk_latest_pointer['shortlist_run_id']} vs "
            f"{clinical_shortlist_latest_pointer['shortlist_run_id']}"
        )
    if str(endpoint_crosswalk_latest_pointer["core_audit_run_id"]) != str(
        clinical_shortlist_latest_pointer["core_audit_run_id"]
    ):
        raise CohortBlueprintError(
            "Endpoint crosswalk latest pointer does not reference the current core audit run: "
            f"{endpoint_crosswalk_latest_pointer['core_audit_run_id']} vs "
            f"{clinical_shortlist_latest_pointer['core_audit_run_id']}"
        )
    if str(endpoint_crosswalk_latest_pointer["parse_run_id"]) != str(
        clinical_shortlist_latest_pointer["parse_run_id"]
    ):
        raise CohortBlueprintError(
            "Endpoint crosswalk latest pointer parse_run_id does not match the clinical shortlist parse_run_id: "
            f"{endpoint_crosswalk_latest_pointer['parse_run_id']} vs "
            f"{clinical_shortlist_latest_pointer['parse_run_id']}"
        )
    if str(endpoint_crosswalk_latest_pointer["source_run_id"]) != str(
        clinical_shortlist_latest_pointer["source_run_id"]
    ):
        raise CohortBlueprintError(
            "Endpoint crosswalk latest pointer source_run_id does not match the clinical shortlist source_run_id: "
            f"{endpoint_crosswalk_latest_pointer['source_run_id']} vs "
            f"{clinical_shortlist_latest_pointer['source_run_id']}"
        )

    input_paths = {
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
        "core_audit_latest_json": resolve_existing_path(
            paths.repo_root,
            str(clinical_shortlist_latest_pointer["core_audit_latest_json"]),
            "clinical core audit latest pointer",
        ),
        "endpoint_candidate_inventory_tsv": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_candidate_inventory_tsv"]),
            "endpoint candidate inventory TSV",
        ),
        "endpoint_crosswalk_tsv": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_crosswalk_tsv"]),
            "endpoint crosswalk TSV",
        ),
        "endpoint_crosswalk_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["endpoint_crosswalk_summary_tsv"]),
            "endpoint crosswalk summary TSV",
        ),
        "endpoint_crosswalk_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(endpoint_crosswalk_latest_pointer["run_log_json"]),
            "endpoint crosswalk run log",
        ),
        "biospecimen_identifier_inventory_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_inventory_tsv"]),
            "biospecimen identifier inventory TSV",
        ),
        "biospecimen_identifier_crosswalk_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_crosswalk_tsv"]),
            "biospecimen identifier crosswalk TSV",
        ),
        "biospecimen_identifier_crosswalk_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_crosswalk_summary_tsv"]),
            "biospecimen identifier crosswalk summary TSV",
        ),
        "biospecimen_identifier_crosswalk_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["run_log_json"]),
            "biospecimen identifier crosswalk run log",
        ),
    }

    clinical_shortlist_run_log = load_json(input_paths["clinical_shortlist_run_log_json"])
    endpoint_crosswalk_run_log = load_json(input_paths["endpoint_crosswalk_run_log_json"])
    biospecimen_crosswalk_run_log = load_json(input_paths["biospecimen_identifier_crosswalk_run_log_json"])
    for run_log, label in [
        (clinical_shortlist_run_log, "clinical shortlist"),
        (endpoint_crosswalk_run_log, "endpoint crosswalk"),
        (biospecimen_crosswalk_run_log, "biospecimen crosswalk"),
    ]:
        if run_log.get("status") != "completed":
            raise CohortBlueprintError(f"Upstream {label} run log is not completed.")

    core_audit_latest_pointer = load_json(input_paths["core_audit_latest_json"])
    require_keys(
        core_audit_latest_pointer,
        {
            "audit_run_id",
            "parse_run_id",
            "source_run_id",
            "clinical_core_field_audit_tsv",
            "run_log_json",
        },
        "Clinical core audit latest pointer",
        input_paths["core_audit_latest_json"],
    )
    if str(core_audit_latest_pointer["audit_run_id"]) != str(
        clinical_shortlist_latest_pointer["core_audit_run_id"]
    ):
        raise CohortBlueprintError(
            "Clinical shortlist latest pointer and core audit latest pointer disagree on audit_run_id: "
            f"{clinical_shortlist_latest_pointer['core_audit_run_id']} vs "
            f"{core_audit_latest_pointer['audit_run_id']}"
        )
    input_paths["clinical_core_field_audit_tsv"] = resolve_existing_path(
        paths.repo_root,
        str(core_audit_latest_pointer["clinical_core_field_audit_tsv"]),
        "clinical core field audit TSV",
    )
    input_paths["core_audit_run_log_json"] = resolve_existing_path(
        paths.repo_root,
        str(core_audit_latest_pointer["run_log_json"]),
        "clinical core audit run log",
    )
    core_audit_run_log = load_json(input_paths["core_audit_run_log_json"])
    if core_audit_run_log.get("status") != "completed":
        raise CohortBlueprintError("Upstream clinical core audit run log is not completed.")

    clinical_biotab_run_log_path_raw = str(
        core_audit_run_log.get("inputs", {}).get("clinical_biotab_run_log_json") or ""
    )
    if not clinical_biotab_run_log_path_raw:
        raise CohortBlueprintError("Clinical core audit run log does not reference a clinical biotab run log.")
    input_paths["clinical_biotab_run_log_json"] = resolve_existing_path(
        paths.repo_root,
        clinical_biotab_run_log_path_raw,
        "clinical biotab run log",
    )
    clinical_biotab_run_log = load_json(input_paths["clinical_biotab_run_log_json"])
    if clinical_biotab_run_log.get("status") != "completed":
        raise CohortBlueprintError("Upstream clinical biotab run log is not completed.")

    source_metadata_tsv_path_raw = str(
        clinical_biotab_run_log.get("inputs", {}).get("source_metadata_tsv") or ""
    )
    source_run_log_path_raw = str(clinical_biotab_run_log.get("inputs", {}).get("source_run_log_json") or "")
    if not source_metadata_tsv_path_raw or not source_run_log_path_raw:
        raise CohortBlueprintError(
            "Clinical biotab run log does not reference source metadata TSV and source run log."
        )
    input_paths["source_metadata_tsv"] = resolve_existing_path(
        paths.repo_root,
        source_metadata_tsv_path_raw,
        "source metadata TSV",
    )
    input_paths["source_run_log_json"] = resolve_existing_path(
        paths.repo_root,
        source_run_log_path_raw,
        "source supplement run log",
    )
    source_run_log = load_json(input_paths["source_run_log_json"])
    if source_run_log.get("status") != "completed":
        raise CohortBlueprintError("Upstream source supplement run log is not completed.")

    clinical_shortlist_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_tsv"])
    clinical_shortlist_summary_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_summary_tsv"])
    clinical_shortlist_by_table_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_by_table_tsv"])
    clinical_core_field_audit_rows = read_tsv_dict_rows(input_paths["clinical_core_field_audit_tsv"])
    endpoint_inventory_rows = read_tsv_dict_rows(input_paths["endpoint_candidate_inventory_tsv"])
    endpoint_crosswalk_rows = read_tsv_dict_rows(input_paths["endpoint_crosswalk_tsv"])
    endpoint_summary_rows = read_tsv_dict_rows(input_paths["endpoint_crosswalk_summary_tsv"])
    biospecimen_inventory_rows = read_tsv_dict_rows(input_paths["biospecimen_identifier_inventory_tsv"])
    biospecimen_crosswalk_rows = read_tsv_dict_rows(input_paths["biospecimen_identifier_crosswalk_tsv"])
    biospecimen_summary_rows = read_tsv_dict_rows(input_paths["biospecimen_identifier_crosswalk_summary_tsv"])
    source_metadata_rows = read_tsv_dict_rows(input_paths["source_metadata_tsv"])

    for label, rows in [
        ("clinical shortlist", clinical_shortlist_rows),
        ("clinical shortlist summary", clinical_shortlist_summary_rows),
        ("clinical shortlist by-table", clinical_shortlist_by_table_rows),
        ("clinical core audit", clinical_core_field_audit_rows),
        ("endpoint inventory", endpoint_inventory_rows),
        ("endpoint crosswalk", endpoint_crosswalk_rows),
        ("endpoint summary", endpoint_summary_rows),
        ("biospecimen inventory", biospecimen_inventory_rows),
        ("biospecimen crosswalk", biospecimen_crosswalk_rows),
        ("biospecimen summary", biospecimen_summary_rows),
        ("source metadata", source_metadata_rows),
    ]:
        if not rows:
            raise CohortBlueprintError(f"Required input table has no rows: {label}")

    return WorkflowInputs(
        clinical_shortlist_latest_pointer=clinical_shortlist_latest_pointer,
        clinical_shortlist_run_log=clinical_shortlist_run_log,
        endpoint_crosswalk_latest_pointer=endpoint_crosswalk_latest_pointer,
        endpoint_crosswalk_run_log=endpoint_crosswalk_run_log,
        biospecimen_crosswalk_latest_pointer=biospecimen_crosswalk_latest_pointer,
        biospecimen_crosswalk_run_log=biospecimen_crosswalk_run_log,
        core_audit_latest_pointer=core_audit_latest_pointer,
        core_audit_run_log=core_audit_run_log,
        clinical_biotab_run_log=clinical_biotab_run_log,
        source_run_log=source_run_log,
        clinical_shortlist_rows=clinical_shortlist_rows,
        clinical_shortlist_summary_rows=clinical_shortlist_summary_rows,
        clinical_shortlist_by_table_rows=clinical_shortlist_by_table_rows,
        clinical_core_field_audit_rows=clinical_core_field_audit_rows,
        endpoint_inventory_rows=endpoint_inventory_rows,
        endpoint_crosswalk_rows=endpoint_crosswalk_rows,
        endpoint_summary_rows=endpoint_summary_rows,
        biospecimen_inventory_rows=biospecimen_inventory_rows,
        biospecimen_crosswalk_rows=biospecimen_crosswalk_rows,
        biospecimen_summary_rows=biospecimen_summary_rows,
        source_metadata_rows=source_metadata_rows,
        input_paths=input_paths,
    )


def build_row_index(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = (str(row.get("table_name") or ""), str(row.get("field_name") or ""))
        if key in index:
            raise CohortBlueprintError(f"Duplicate row key found in indexed input table: {key}")
        index[key] = row
    return index


def clinical_field_row(
    blueprint_run_id: str,
    bucket: str,
    source_row: dict[str, str],
    upstream_source_type: str,
    upstream_bucket_or_role: str,
    upstream_rule: str,
    proposed_role: str,
    rationale: str,
    dependency_for_join: str,
    ambiguity_flag: bool,
) -> dict[str, Any]:
    return {
        "blueprint_run_id": blueprint_run_id,
        "field_bucket": bucket,
        "source_layer": "clinical",
        "table_name": str(source_row["table_name"]),
        "field_name": str(source_row["field_name"]),
        "proposed_role": proposed_role,
        "rationale": rationale,
        "dependency_for_join": dependency_for_join,
        "ambiguity_flag": ambiguity_flag,
        "notes_placeholder": NOTES_PLACEHOLDER,
        "upstream_source_type": upstream_source_type,
        "upstream_bucket_or_role": upstream_bucket_or_role,
        "upstream_rule": upstream_rule,
        "manual_review_priority": str(source_row.get("manual_review_priority") or ""),
        "missing_like_fraction": str(source_row.get("missing_like_fraction") or ""),
        "non_missing_count": str(source_row.get("non_missing_count") or ""),
        "distinct_non_missing_count": str(source_row.get("distinct_non_missing_count") or ""),
        "paired_field_name": "",
        "evidence_rule": "",
        "evidence_rate": "",
    }


def endpoint_field_row(
    blueprint_run_id: str,
    bucket: str,
    source_row: dict[str, str],
    proposed_role: str,
    rationale: str,
    dependency_for_join: str,
    ambiguity_flag: bool,
) -> dict[str, Any]:
    return {
        "blueprint_run_id": blueprint_run_id,
        "field_bucket": bucket,
        "source_layer": "endpoint",
        "table_name": str(source_row["table_name"]),
        "field_name": str(source_row["field_name"]),
        "proposed_role": proposed_role,
        "rationale": rationale,
        "dependency_for_join": dependency_for_join,
        "ambiguity_flag": ambiguity_flag,
        "notes_placeholder": NOTES_PLACEHOLDER,
        "upstream_source_type": "endpoint_crosswalk",
        "upstream_bucket_or_role": str(source_row.get("crosswalk_role") or ""),
        "upstream_rule": str(source_row.get("crosswalk_rule") or ""),
        "manual_review_priority": str(source_row.get("manual_review_priority") or ""),
        "missing_like_fraction": str(source_row.get("missing_like_fraction") or ""),
        "non_missing_count": str(source_row.get("non_missing_count") or ""),
        "distinct_non_missing_count": str(source_row.get("distinct_non_missing_count") or ""),
        "paired_field_name": "",
        "evidence_rule": str(source_row.get("crosswalk_rule") or ""),
        "evidence_rate": "",
    }


def biospecimen_field_row(
    blueprint_run_id: str,
    bucket: str,
    source_row: dict[str, str],
    proposed_role: str,
    rationale: str,
    dependency_for_join: str,
    ambiguity_flag: bool,
) -> dict[str, Any]:
    return {
        "blueprint_run_id": blueprint_run_id,
        "field_bucket": bucket,
        "source_layer": "biospecimen",
        "table_name": str(source_row["table_name"]),
        "field_name": str(source_row["field_name"]),
        "proposed_role": proposed_role,
        "rationale": rationale,
        "dependency_for_join": dependency_for_join,
        "ambiguity_flag": ambiguity_flag,
        "notes_placeholder": NOTES_PLACEHOLDER,
        "upstream_source_type": "biospecimen_identifier_crosswalk",
        "upstream_bucket_or_role": str(source_row.get("crosswalk_role") or ""),
        "upstream_rule": str(source_row.get("crosswalk_rule") or ""),
        "manual_review_priority": str(source_row.get("manual_review_priority") or ""),
        "missing_like_fraction": str(source_row.get("missing_like_fraction") or ""),
        "non_missing_count": str(source_row.get("non_missing_count") or ""),
        "distinct_non_missing_count": str(source_row.get("distinct_non_missing_count") or ""),
        "paired_field_name": str(source_row.get("paired_field_name") or ""),
        "evidence_rule": str(source_row.get("pattern_rule") or ""),
        "evidence_rate": str(source_row.get("pattern_pass_fraction") or ""),
    }


def add_field_row(
    row: dict[str, Any],
    target_rows: list[dict[str, Any]],
    seen_keys: set[tuple[str, str, str]],
) -> None:
    key = (str(row["source_layer"]), str(row["table_name"]), str(row["field_name"]))
    if key in seen_keys:
        return
    seen_keys.add(key)
    target_rows.append(row)


def build_field_tables(
    blueprint_run_id: str,
    workflow_inputs: WorkflowInputs,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    shortlist_by_key = build_row_index(workflow_inputs.clinical_shortlist_rows)
    core_audit_by_key = build_row_index(workflow_inputs.clinical_core_field_audit_rows)

    required_rows: list[dict[str, Any]] = []
    optional_rows: list[dict[str, Any]] = []
    deferred_rows: list[dict[str, Any]] = []
    seen_keys: set[tuple[str, str, str]] = set()

    for table_name, field_name in [
        ("clinical_patient", "bcr_patient_uuid"),
        ("clinical_patient", "bcr_patient_barcode"),
    ]:
        source_row = core_audit_by_key.get((table_name, field_name))
        if source_row is None:
            raise CohortBlueprintError(f"Required core audit field not found: {table_name}.{field_name}")
        add_field_row(
            clinical_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="required_for_first_baseline",
                source_row=source_row,
                upstream_source_type="clinical_core_field_audit",
                upstream_bucket_or_role="manual_blueprint_backbone",
                upstream_rule="cohort_blueprint_required_case_identifier",
                proposed_role="case_join_identifier",
                rationale="Retain both audited clinical patient identifiers as the case backbone for the first baseline blueprint.",
                dependency_for_join="case_to_clinical_and_case_to_followup_preparation",
                ambiguity_flag=False,
            ),
            required_rows,
            seen_keys,
        )

    for source_row in workflow_inputs.clinical_shortlist_rows:
        if str(source_row.get("shortlist_bucket") or "") != "usable_baseline":
            continue
        add_field_row(
            clinical_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="required_for_first_baseline",
                source_row=source_row,
                upstream_source_type="clinical_shortlist",
                upstream_bucket_or_role=str(source_row.get("shortlist_bucket") or ""),
                upstream_rule=str(source_row.get("shortlist_rule") or ""),
                proposed_role="baseline_covariate",
                rationale="Carry forward audited usable baseline fields as required case-level inputs for the first TCGA-only baseline study blueprint.",
                dependency_for_join="case_level_baseline_core",
                ambiguity_flag=False,
            ),
            required_rows,
            seen_keys,
        )

    for table_name, field_name in [
        ("biospecimen_sample", "bcr_patient_uuid"),
        ("biospecimen_sample", "bcr_sample_barcode"),
        ("biospecimen_sample", "bcr_sample_uuid"),
    ]:
        matching_rows = [
            row
            for row in workflow_inputs.biospecimen_crosswalk_rows
            if str(row.get("table_name") or "") == table_name
            and str(row.get("field_name") or "") == field_name
        ]
        if not matching_rows:
            raise CohortBlueprintError(f"Required biospecimen blueprint field not found: {table_name}.{field_name}")
        add_field_row(
            biospecimen_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="required_for_first_baseline",
                source_row=matching_rows[0],
                proposed_role="biospecimen_sample_anchor",
                rationale="Use the audited biospecimen sample identifiers as the required specimen anchor beneath the patient/case unit.",
                dependency_for_join="case_to_sample_anchor",
                ambiguity_flag=False,
            ),
            required_rows,
            seen_keys,
        )

    for source_row in workflow_inputs.clinical_shortlist_rows:
        if str(source_row.get("shortlist_bucket") or "") != "usable_treatment_proxy":
            continue
        add_field_row(
            clinical_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="optional_for_first_baseline",
                source_row=source_row,
                upstream_source_type="clinical_shortlist",
                upstream_bucket_or_role=str(source_row.get("shortlist_bucket") or ""),
                upstream_rule=str(source_row.get("shortlist_rule") or ""),
                proposed_role="treatment_proxy_field",
                rationale="Preserve audited treatment-proxy fields as optional first-study context without turning them into harmonized exposure variables.",
                dependency_for_join="optional_case_level_treatment_proxy",
                ambiguity_flag=True,
            ),
            optional_rows,
            seen_keys,
        )

    for table_name, field_name in [
        ("clinical_drug", "bcr_patient_uuid"),
        ("clinical_drug", "bcr_patient_barcode"),
        ("clinical_radiation", "bcr_patient_uuid"),
        ("clinical_radiation", "bcr_patient_barcode"),
        ("clinical_follow_up_v4_0", "bcr_patient_uuid"),
        ("clinical_follow_up_v4_0", "bcr_patient_barcode"),
        ("clinical_follow_up_v4_0", "bcr_followup_barcode"),
        ("clinical_follow_up_v4_0", "bcr_followup_uuid"),
    ]:
        source_row = core_audit_by_key.get((table_name, field_name))
        if source_row is None:
            raise CohortBlueprintError(f"Expected clinical core audit field not found: {table_name}.{field_name}")
        add_field_row(
            clinical_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="optional_for_first_baseline",
                source_row=source_row,
                upstream_source_type="clinical_core_field_audit",
                upstream_bucket_or_role="manual_blueprint_join_support",
                upstream_rule="cohort_blueprint_optional_join_identifier",
                proposed_role=(
                    "followup_join_identifier" if table_name == "clinical_follow_up_v4_0" else "case_join_identifier"
                ),
                rationale=(
                    "Keep the audited patient and follow-up identifiers needed to attach optional treatment and endpoint-preparation layers later."
                ),
                dependency_for_join=(
                    "case_to_followup_or_case_to_treatment_optional_link"
                    if table_name == "clinical_follow_up_v4_0"
                    else "case_to_treatment_optional_link"
                ),
                ambiguity_flag=(table_name == "clinical_follow_up_v4_0"),
            ),
            optional_rows,
            seen_keys,
        )

    for source_row in workflow_inputs.endpoint_crosswalk_rows:
        role = str(source_row.get("crosswalk_role") or "")
        if role not in {"primary_candidate", "overlapping_candidate"}:
            continue
        add_field_row(
            endpoint_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="optional_for_first_baseline",
                source_row=source_row,
                proposed_role="endpoint_candidate_join_preparation",
                rationale="Carry audited endpoint candidate fields only as join-preparation inputs; the blueprint does not freeze a final endpoint definition.",
                dependency_for_join="case_to_endpoint_candidate_preparation",
                ambiguity_flag=(role == "overlapping_candidate"),
            ),
            optional_rows,
            seen_keys,
        )

    for source_row in workflow_inputs.biospecimen_crosswalk_rows:
        if (
            str(source_row.get("identifier_linkage_family") or "") == "sample_identifier_like"
            and str(source_row.get("crosswalk_role") or "") == "secondary_link_candidate"
            and str(source_row.get("table_name") or "") == "biospecimen_sample"
        ):
            add_field_row(
                biospecimen_field_row(
                    blueprint_run_id=blueprint_run_id,
                    bucket="optional_for_first_baseline",
                    source_row=source_row,
                    proposed_role="biospecimen_sample_helper",
                    rationale="Retain audited sample helper fields as optional context around the required sample anchor without making them baseline requirements.",
                    dependency_for_join="sample_anchor_context_only",
                    ambiguity_flag=True,
                ),
                optional_rows,
                seen_keys,
            )

    for source_row in workflow_inputs.endpoint_crosswalk_rows:
        if str(source_row.get("crosswalk_role") or "") != "ambiguous_candidate":
            continue
        add_field_row(
            endpoint_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="deferred_for_later_expansion",
                source_row=source_row,
                proposed_role="deferred_endpoint_candidate",
                rationale="Defer endpoint candidates that remain sparse or unresolved until later endpoint review rather than using them automatically in first-baseline construction.",
                dependency_for_join="deferred_endpoint_resolution",
                ambiguity_flag=True,
            ),
            deferred_rows,
            seen_keys,
        )

    biospecimen_family_deferred = {
        "portion_identifier_like",
        "analyte_identifier_like",
        "slide_identifier_like",
        "aliquot_identifier_like",
        "multi_level_or_unclear_identifier_like",
    }
    for source_row in workflow_inputs.biospecimen_crosswalk_rows:
        family = str(source_row.get("identifier_linkage_family") or "")
        role = str(source_row.get("crosswalk_role") or "")
        if family in biospecimen_family_deferred:
            add_field_row(
                biospecimen_field_row(
                    blueprint_run_id=blueprint_run_id,
                    bucket="deferred_for_later_expansion",
                    source_row=source_row,
                    proposed_role=(
                        "ambiguous_biospecimen_identifier"
                        if family == "multi_level_or_unclear_identifier_like"
                        else "child_biospecimen_identifier"
                    ),
                    rationale=(
                        "Keep deeper biospecimen identifiers for later expansion because the first baseline blueprint requires only the sample anchor, not mandatory child-layer joins."
                    ),
                    dependency_for_join=(
                        "later_child_layer_expansion" if family != "multi_level_or_unclear_identifier_like" else "manual_identifier_review"
                    ),
                    ambiguity_flag=True,
                ),
                deferred_rows,
                seen_keys,
            )
            continue

        if family in {"patient_identifier_like", "sample_identifier_like"} and role in {
            "overlapping_link_candidate",
            "side_table_candidate",
            "ambiguous_candidate",
        }:
            add_field_row(
                biospecimen_field_row(
                    blueprint_run_id=blueprint_run_id,
                    bucket="deferred_for_later_expansion",
                    source_row=source_row,
                    proposed_role="later_biospecimen_side_or_overlap_identifier",
                    rationale="Keep overlapping and side-table biospecimen identifiers for later expansion rather than making them automatic first-baseline requirements.",
                    dependency_for_join="later_side_link_review",
                    ambiguity_flag=True,
                ),
                deferred_rows,
                seen_keys,
            )

    for source_row in workflow_inputs.biospecimen_crosswalk_rows:
        add_field_row(
            biospecimen_field_row(
                blueprint_run_id=blueprint_run_id,
                bucket="deferred_for_later_expansion",
                source_row=source_row,
                proposed_role="later_biospecimen_side_or_overlap_identifier",
                rationale="Carry forward any remaining audited biospecimen identifier evidence as later-expansion material instead of dropping it from the blueprint.",
                dependency_for_join="later_side_link_review",
                ambiguity_flag=True,
            ),
            deferred_rows,
            seen_keys,
        )

    required_rows.sort(key=lambda row: (SOURCE_LAYER_ORDER.index(str(row["source_layer"])), str(row["table_name"]), str(row["field_name"])))
    optional_rows.sort(key=lambda row: (SOURCE_LAYER_ORDER.index(str(row["source_layer"])), str(row["table_name"]), str(row["field_name"])))
    deferred_rows.sort(key=lambda row: (SOURCE_LAYER_ORDER.index(str(row["source_layer"])), str(row["table_name"]), str(row["field_name"])))
    return required_rows, optional_rows, deferred_rows


def build_join_path_rows(blueprint_run_id: str, workflow_inputs: WorkflowInputs) -> list[dict[str, Any]]:
    endpoint_by_key = build_row_index(workflow_inputs.endpoint_crosswalk_rows)
    biospecimen_by_key = build_row_index(workflow_inputs.biospecimen_crosswalk_rows)

    patient_vital = endpoint_by_key[("clinical_patient", "vital_status")]
    followup_vital = endpoint_by_key[("clinical_follow_up_v4_0", "vital_status")]
    followup_last_contact = endpoint_by_key[("clinical_follow_up_v4_0", "last_contact_days_to")]
    portion_sample = biospecimen_by_key[("biospecimen_portion", "bcr_sample_barcode")]
    protocol_analyte = biospecimen_by_key[("biospecimen_protocol", "bcr_analyte_barcode")]
    shipment_portion = biospecimen_by_key[("biospecimen_shipment_portion", "bcr_shipment_portion_uuid")]
    diagnostic_slide = biospecimen_by_key[("biospecimen_diagnostic_slides", "ffpe_slide_uuid")]

    return [
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 1,
            "left_node": "patient/case",
            "right_node": "clinical_patient",
            "left_table_name": "",
            "right_table_name": "clinical_patient",
            "left_field_name": "case_identifier",
            "right_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "join_keys_or_evidence": "clinical_patient carries both audited case identifiers with 1097 non-missing values each in the core audit.",
            "join_strength": "strong_current_evidence",
            "link_timing": "link_now",
            "rationale": "Use clinical_patient as the central case-level baseline table for the first blueprint.",
            "ambiguity_flag": False,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 2,
            "left_node": "patient/case",
            "right_node": "clinical_drug",
            "left_table_name": "clinical_patient",
            "right_table_name": "clinical_drug",
            "left_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "right_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "join_keys_or_evidence": "clinical_drug contains repeated patient identifiers but remains one-to-many with 2406 rows over 780 distinct patients.",
            "join_strength": "workable_but_needs_review",
            "link_timing": "link_later",
            "rationale": "Treatment rows can be attached later as optional treatment-proxy evidence, not as required baseline rows.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 3,
            "left_node": "patient/case",
            "right_node": "clinical_radiation",
            "left_table_name": "clinical_patient",
            "right_table_name": "clinical_radiation",
            "left_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "right_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "join_keys_or_evidence": "clinical_radiation contains repeated patient identifiers but remains one-to-many with 618 rows over 528 distinct patients.",
            "join_strength": "workable_but_needs_review",
            "link_timing": "link_later",
            "rationale": "Radiation rows should remain optional treatment-proxy attachments instead of required first-baseline rows.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 4,
            "left_node": "patient/case",
            "right_node": "clinical_patient endpoint-like fields",
            "left_table_name": "clinical_patient",
            "right_table_name": "clinical_patient",
            "left_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "right_field_name": "vital_status / last_contact_days_to / tumor_status / new_tumor_event_dx_indicator",
            "join_keys_or_evidence": (
                f"clinical_patient endpoint-like fields already sit on the case table; vital_status has "
                f"{patient_vital['non_missing_count']} non-missing rows."
            ),
            "join_strength": "strong_current_evidence",
            "link_timing": "prepare_only",
            "rationale": "Patient-table endpoint candidates can be carried as case-level preparation fields without creating a final endpoint.",
            "ambiguity_flag": False,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 5,
            "left_node": "patient/case",
            "right_node": "clinical_follow_up_v4_0 endpoint-like fields",
            "left_table_name": "clinical_patient",
            "right_table_name": "clinical_follow_up_v4_0",
            "left_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "right_field_name": "bcr_patient_barcode / bcr_patient_uuid",
            "join_keys_or_evidence": (
                f"follow-up endpoint-like fields reuse patient identifiers, but overlap with patient-table signals; "
                f"follow-up vital_status has {followup_vital['non_missing_count']} non-missing rows and "
                f"last_contact_days_to has {followup_last_contact['non_missing_count']}."
            ),
            "join_strength": "workable_but_needs_review",
            "link_timing": "prepare_only",
            "rationale": "Follow-up endpoint candidates are usable for later reconciliation but should not be treated as automatically equivalent to patient-table fields.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 6,
            "left_node": "patient/case",
            "right_node": "biospecimen_sample",
            "left_table_name": "clinical_patient",
            "right_table_name": "biospecimen_sample",
            "left_field_name": "bcr_patient_uuid / bcr_patient_barcode",
            "right_field_name": "bcr_patient_uuid / bcr_sample_barcode / bcr_sample_uuid",
            "join_keys_or_evidence": "biospecimen_sample provides the main specimen anchor, but the saved evidence is stronger inside the biospecimen chain than across clinical-to-biospecimen reconciliation.",
            "join_strength": "workable_but_needs_review",
            "link_timing": "link_now",
            "rationale": "The blueprint can anchor on biospecimen_sample now while preserving cross-layer patient identifier review for later cohort construction.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 7,
            "left_node": "sample",
            "right_node": "portion",
            "left_table_name": "biospecimen_sample",
            "right_table_name": "biospecimen_portion",
            "left_field_name": "bcr_sample_barcode / bcr_sample_uuid",
            "right_field_name": "bcr_sample_barcode / bcr_portion_barcode / bcr_portion_uuid",
            "join_keys_or_evidence": (
                f"biospecimen_portion repeats bcr_sample_barcode with barcode-prefix evidence to portion barcode at "
                f"{portion_sample['pattern_pass_fraction']}."
            ),
            "join_strength": "strong_current_evidence",
            "link_timing": "link_later",
            "rationale": "Sample-to-portion is the strongest saved child-layer handoff and can be retained for later expansion.",
            "ambiguity_flag": False,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 8,
            "left_node": "portion",
            "right_node": "analyte",
            "left_table_name": "biospecimen_portion",
            "right_table_name": "biospecimen_analyte",
            "left_field_name": "bcr_portion_barcode / bcr_portion_uuid",
            "right_field_name": "bcr_analyte_barcode / bcr_analyte_uuid / subportion_sequence",
            "join_keys_or_evidence": "analyte rows retain sample barcode plus subportion_sequence, but no direct saved portion barcode appears in biospecimen_analyte.",
            "join_strength": "workable_but_needs_review",
            "link_timing": "link_later",
            "rationale": "Portion-to-analyte looks workable but still needs careful review because the saved bridge is indirect.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 9,
            "left_node": "analyte",
            "right_node": "slide",
            "left_table_name": "biospecimen_analyte",
            "right_table_name": "biospecimen_slide",
            "left_field_name": "bcr_analyte_barcode / bcr_analyte_uuid",
            "right_field_name": "bcr_slide_barcode / bcr_slide_uuid / bcr_sample_barcode",
            "join_keys_or_evidence": "No direct saved analyte identifier appears in biospecimen_slide; slide rows carry sample barcode and slide ids only.",
            "join_strength": "deferred",
            "link_timing": "link_later",
            "rationale": "An analyte-to-slide step should remain deferred until stronger saved linkage evidence is reviewed.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 10,
            "left_node": "analyte",
            "right_node": "aliquot",
            "left_table_name": "biospecimen_analyte",
            "right_table_name": "biospecimen_aliquot",
            "left_field_name": "bcr_analyte_barcode / bcr_analyte_uuid",
            "right_field_name": "bcr_aliquot_barcode / bcr_aliquot_uuid / bcr_sample_barcode",
            "join_keys_or_evidence": "No direct saved analyte identifier appears in biospecimen_aliquot; aliquot rows carry sample barcode and aliquot ids only.",
            "join_strength": "deferred",
            "link_timing": "link_later",
            "rationale": "An analyte-to-aliquot step should remain deferred until a stronger saved bridge is available.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 11,
            "left_node": "sample/analyte",
            "right_node": "protocol",
            "left_table_name": "biospecimen_sample / biospecimen_analyte",
            "right_table_name": "biospecimen_protocol",
            "left_field_name": "bcr_sample_barcode / bcr_analyte_barcode",
            "right_field_name": "bcr_sample_barcode / bcr_analyte_barcode",
            "join_keys_or_evidence": (
                f"biospecimen_protocol carries both sample and analyte side-table candidates; protocol analyte barcode "
                f"shows barcode-prefix evidence at {protocol_analyte['pattern_pass_fraction']}."
            ),
            "join_strength": "workable_but_needs_review",
            "link_timing": "link_later",
            "rationale": "Protocol remains a side-link layer that looks usable later but should not be mandatory for the first baseline cohort.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 12,
            "left_node": "sample/portion",
            "right_node": "shipment_portion",
            "left_table_name": "biospecimen_sample / biospecimen_portion",
            "right_table_name": "biospecimen_shipment_portion",
            "left_field_name": "bcr_sample_barcode / portion_number / portion_sequence",
            "right_field_name": "bcr_shipment_portion_uuid / bcr_sample_barcode / portion_number / portion_sequence",
            "join_keys_or_evidence": (
                f"shipment_portion retains sample barcode plus portion helper fields and a side-table portion UUID with "
                f"{shipment_portion['pattern_pass_fraction']} same-row evidence."
            ),
            "join_strength": "workable_but_needs_review",
            "link_timing": "link_later",
            "rationale": "Shipment data can support later specimen-tracking review but is not required for the first baseline blueprint.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "step_order": 13,
            "left_node": "sample/slide",
            "right_node": "diagnostic_slides",
            "left_table_name": "biospecimen_sample / biospecimen_slide",
            "right_table_name": "biospecimen_diagnostic_slides",
            "left_field_name": "bcr_sample_barcode / bcr_slide_barcode",
            "right_field_name": "bcr_sample_barcode / ffpe_slide_uuid / bcr_patient_barcode",
            "join_keys_or_evidence": (
                f"diagnostic_slides provides sample and patient barcode evidence plus ffpe_slide_uuid as a side-table slide key with "
                f"{diagnostic_slide['pattern_pass_fraction']} saved evidence."
            ),
            "join_strength": "deferred",
            "link_timing": "link_later",
            "rationale": "Diagnostic slide linkage should remain deferred because its saved linkage is side-table evidence rather than part of the required first-baseline backbone.",
            "ambiguity_flag": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
    ]


def build_ambiguity_rows(
    blueprint_run_id: str,
    workflow_inputs: WorkflowInputs,
) -> list[dict[str, Any]]:
    endpoint_by_key = build_row_index(workflow_inputs.endpoint_crosswalk_rows)
    core_audit_by_key = build_row_index(workflow_inputs.clinical_core_field_audit_rows)
    biospecimen_by_key = build_row_index(workflow_inputs.biospecimen_crosswalk_rows)

    clinical_patient_row_count = parse_int(
        core_audit_by_key[("clinical_patient", "bcr_patient_barcode")]["non_missing_count"],
        "clinical_patient bcr_patient_barcode non_missing_count",
    )
    clinical_source_case_coverage = unique_case_submitter_count(
        workflow_inputs.source_metadata_rows,
        "clinical",
    )
    biospecimen_patient_uuid_distinct = parse_int(
        biospecimen_by_key[("biospecimen_sample", "bcr_patient_uuid")]["distinct_non_missing_count"],
        "biospecimen_sample bcr_patient_uuid distinct_non_missing_count",
    )

    return [
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A01",
            "ambiguity_category": "endpoint_overlap_patient_vs_followup",
            "source_layer": "endpoint",
            "table_name": "clinical_patient | clinical_follow_up_v4_0",
            "field_name_or_step": "vital_status",
            "severity": "high",
            "why_unresolved": "The same endpoint-like field exists in both patient and follow-up tables and should not be treated as automatically equivalent.",
            "current_saved_evidence": (
                f"clinical_patient.vital_status has {endpoint_by_key[('clinical_patient', 'vital_status')]['non_missing_count']} non-missing rows; "
                f"clinical_follow_up_v4_0.vital_status has {endpoint_by_key[('clinical_follow_up_v4_0', 'vital_status')]['non_missing_count']}."
            ),
            "proposed_blueprint_handling": "Carry both rows as endpoint join-preparation fields and require later reconciliation before endpoint freeze.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A02",
            "ambiguity_category": "endpoint_overlap_patient_vs_followup",
            "source_layer": "endpoint",
            "table_name": "clinical_patient | clinical_follow_up_v4_0",
            "field_name_or_step": "last_contact_days_to",
            "severity": "high",
            "why_unresolved": "Patient and follow-up tables both carry last-contact timing, but their later use still needs explicit reconciliation.",
            "current_saved_evidence": (
                f"clinical_patient.last_contact_days_to has {endpoint_by_key[('clinical_patient', 'last_contact_days_to')]['non_missing_count']} non-missing rows; "
                f"clinical_follow_up_v4_0.last_contact_days_to has {endpoint_by_key[('clinical_follow_up_v4_0', 'last_contact_days_to')]['non_missing_count']}."
            ),
            "proposed_blueprint_handling": "Keep both rows optional and endpoint-preparatory only; do not choose one as the final censoring source here.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A03",
            "ambiguity_category": "endpoint_overlap_patient_vs_followup",
            "source_layer": "endpoint",
            "table_name": "clinical_patient | clinical_follow_up_v4_0",
            "field_name_or_step": "tumor_status",
            "severity": "high",
            "why_unresolved": "Tumor-status signals overlap across patient and follow-up tables and may reflect related but not guaranteed identical source semantics.",
            "current_saved_evidence": (
                f"clinical_patient.tumor_status has {endpoint_by_key[('clinical_patient', 'tumor_status')]['non_missing_count']} non-missing rows; "
                f"clinical_follow_up_v4_0.tumor_status has {endpoint_by_key[('clinical_follow_up_v4_0', 'tumor_status')]['non_missing_count']}."
            ),
            "proposed_blueprint_handling": "Preserve both rows for later endpoint review and keep them out of any frozen progression definition.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A04",
            "ambiguity_category": "endpoint_overlap_patient_vs_followup",
            "source_layer": "endpoint",
            "table_name": "clinical_patient | clinical_follow_up_v4_0",
            "field_name_or_step": "new_tumor_event_dx_indicator",
            "severity": "high",
            "why_unresolved": "New-tumor-event indicators overlap across patient and follow-up, but their relative timing and preferred later use remain unsettled.",
            "current_saved_evidence": (
                f"clinical_patient.new_tumor_event_dx_indicator has {endpoint_by_key[('clinical_patient', 'new_tumor_event_dx_indicator')]['non_missing_count']} non-missing rows; "
                f"clinical_follow_up_v4_0.new_tumor_event_dx_indicator has {endpoint_by_key[('clinical_follow_up_v4_0', 'new_tumor_event_dx_indicator')]['non_missing_count']}."
            ),
            "proposed_blueprint_handling": "Keep both rows as later endpoint-reconciliation targets rather than automatic event indicators.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A05",
            "ambiguity_category": "sparse_endpoint_timing",
            "source_layer": "endpoint",
            "table_name": "clinical_patient | clinical_follow_up_v4_0",
            "field_name_or_step": "death_days_to",
            "severity": "high",
            "why_unresolved": "Death timing exists in both tables but remains too sparse to promote directly into a finalized endpoint workflow here.",
            "current_saved_evidence": (
                f"clinical_patient.death_days_to has {endpoint_by_key[('clinical_patient', 'death_days_to')]['non_missing_count']} non-missing rows; "
                f"clinical_follow_up_v4_0.death_days_to has {endpoint_by_key[('clinical_follow_up_v4_0', 'death_days_to')]['non_missing_count']}."
            ),
            "proposed_blueprint_handling": "Record death timing as deferred endpoint evidence only and revisit it during endpoint freeze.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A06",
            "ambiguity_category": "all_missing_progression_timing",
            "source_layer": "endpoint",
            "table_name": "clinical_patient",
            "field_name_or_step": "days_to_patient_progression_free | days_to_tumor_progression",
            "severity": "high",
            "why_unresolved": "Progression timing fields were surfaced only as ambiguous candidates and remain unusable from the saved audit outputs.",
            "current_saved_evidence": "Both progression timing fields have 0 non-missing rows in the saved endpoint crosswalk.",
            "proposed_blueprint_handling": "Keep both fields deferred and exclude them from any automatic first-baseline endpoint preparation.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A07",
            "ambiguity_category": "clinical_vs_biospecimen_patient_identifiers",
            "source_layer": "clinical | biospecimen",
            "table_name": "clinical_patient | biospecimen_sample | biospecimen_diagnostic_slides",
            "field_name_or_step": "bcr_patient_barcode | bcr_patient_uuid",
            "severity": "high",
            "why_unresolved": "Clinical and biospecimen layers expose both patient barcodes and UUIDs, but the blueprint should not assume a single harmonized patient identifier without later reconciliation.",
            "current_saved_evidence": (
                f"clinical_patient has {clinical_patient_row_count} non-missing values for both patient identifiers; "
                f"biospecimen_sample.bcr_patient_uuid shows {biospecimen_patient_uuid_distinct} distinct patients; "
                f"biospecimen_diagnostic_slides.bcr_patient_barcode is present in 1142 rows."
            ),
            "proposed_blueprint_handling": "Use both identifier forms as evidence-bearing join keys now and leave final harmonization rules for later cohort construction.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A08",
            "ambiguity_category": "case_count_mismatch",
            "source_layer": "clinical | source | biospecimen",
            "table_name": "clinical_patient | source_metadata | biospecimen_sample",
            "field_name_or_step": "case coverage counts",
            "severity": "high",
            "why_unresolved": "The saved clinical, source-acquisition, and biospecimen layers do not expose the exact same patient-count signal, so final cohort construction still needs explicit reconciliation logic.",
            "current_saved_evidence": (
                f"clinical_patient row count is {clinical_patient_row_count}; "
                f"clinical source unique case coverage is {clinical_source_case_coverage}; "
                f"biospecimen_sample patient UUID distinct count is {biospecimen_patient_uuid_distinct}."
            ),
            "proposed_blueprint_handling": "Flag cross-layer case counts as a blocking reconciliation task before any final cohort table is frozen.",
            "blocks_final_cohort_construction": True,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A09",
            "ambiguity_category": "treatment_table_multiplicity",
            "source_layer": "clinical",
            "table_name": "clinical_drug | clinical_radiation",
            "field_name_or_step": "patient-linked treatment rows",
            "severity": "medium",
            "why_unresolved": "Treatment tables are one-to-many and include optional proxy signals, so later aggregation rules still need explicit design choices.",
            "current_saved_evidence": "clinical_drug has 2406 rows over 780 distinct patients; clinical_radiation has 618 rows over 528 distinct patients.",
            "proposed_blueprint_handling": "Keep treatment fields optional and proxy-only until later cohort construction defines aggregation or prioritization rules.",
            "blocks_final_cohort_construction": False,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A10",
            "ambiguity_category": "missing_direct_child_link",
            "source_layer": "biospecimen",
            "table_name": "biospecimen_analyte | biospecimen_slide",
            "field_name_or_step": "analyte -> slide",
            "severity": "medium",
            "why_unresolved": "The saved biospecimen crosswalk does not expose a direct analyte identifier in the slide table.",
            "current_saved_evidence": "biospecimen_slide retains sample barcode and slide ids, while analyte ids appear only in biospecimen_analyte and protocol side links.",
            "proposed_blueprint_handling": "Keep analyte-to-slide deferred and use sample-level or side-table evidence only during later expansion review.",
            "blocks_final_cohort_construction": False,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "ambiguity_id": "A11",
            "ambiguity_category": "missing_direct_child_link",
            "source_layer": "biospecimen",
            "table_name": "biospecimen_analyte | biospecimen_aliquot",
            "field_name_or_step": "analyte -> aliquot",
            "severity": "medium",
            "why_unresolved": "The saved biospecimen crosswalk does not expose a direct analyte identifier in the aliquot table.",
            "current_saved_evidence": "biospecimen_aliquot carries sample barcode and aliquot identifiers only; analyte identifiers remain upstream.",
            "proposed_blueprint_handling": "Keep analyte-to-aliquot deferred and revisit it only during later child-layer expansion.",
            "blocks_final_cohort_construction": False,
            "notes_placeholder": NOTES_PLACEHOLDER,
        },
    ]


def biospecimen_by_key_for_summary(
    workflow_inputs: WorkflowInputs,
    key: tuple[str, str],
    field_name: str,
) -> str:
    row = build_row_index(workflow_inputs.biospecimen_crosswalk_rows)[key]
    return str(row[field_name])


def build_summary_rows(
    blueprint_run_id: str,
    workflow_inputs: WorkflowInputs,
    required_rows: list[dict[str, Any]],
    optional_rows: list[dict[str, Any]],
    deferred_rows: list[dict[str, Any]],
    join_path_rows: list[dict[str, Any]],
    ambiguity_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    source_run_id_set = {
        str(workflow_inputs.clinical_shortlist_latest_pointer["source_run_id"]),
        str(workflow_inputs.biospecimen_crosswalk_latest_pointer["source_run_id"]),
    }
    clinical_source_case_coverage = unique_case_submitter_count(workflow_inputs.source_metadata_rows, "clinical")
    source_summary_rows = [
        row for row in workflow_inputs.source_metadata_rows if str(row.get("source_class") or "") == "clinical"
    ]
    source_metadata_file_count = len(source_summary_rows)
    summary_rows = [
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "design",
            "summary_metric": "proposed_unit_of_analysis",
            "summary_value": "patient/case",
            "notes": "The first baseline blueprint keeps the analysis unit at patient/case level.",
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "design",
            "summary_metric": "proposed_biospecimen_anchor",
            "summary_value": "sample",
            "notes": "biospecimen_sample is the required specimen anchor beneath the case unit.",
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "design",
            "summary_metric": "endpoint_policy",
            "summary_value": "candidate_join_preparation_only",
            "notes": "Endpoint rows remain preparatory only; the final endpoint is not frozen here.",
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "design",
            "summary_metric": "treatment_policy",
            "summary_value": "proxy_only_optional",
            "notes": "Treatment rows remain optional proxy evidence only.",
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "field_counts",
            "summary_metric": "required_field_count",
            "summary_value": len(required_rows),
            "notes": "Required first-baseline fields across clinical and biospecimen sample anchor layers.",
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "field_counts",
            "summary_metric": "optional_field_count",
            "summary_value": len(optional_rows),
            "notes": "Optional attach-later first-study fields.",
        },
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "field_counts",
            "summary_metric": "deferred_field_count",
            "summary_value": len(deferred_rows),
            "notes": "Deferred later-expansion or unresolved fields.",
        },
    ]

    for bucket_label, rows in [
        ("required_for_first_baseline", required_rows),
        ("optional_for_first_baseline", optional_rows),
        ("deferred_for_later_expansion", deferred_rows),
    ]:
        for source_layer in SOURCE_LAYER_ORDER:
            layer_count = sum(str(row["source_layer"]) == source_layer for row in rows)
            summary_rows.append(
                {
                    "blueprint_run_id": blueprint_run_id,
                    "summary_section": "field_counts",
                    "summary_metric": f"{bucket_label}_{source_layer}_count",
                    "summary_value": layer_count,
                    "notes": f"{bucket_label} rows carried from the {source_layer} layer.",
                }
            )

    summary_rows.append(
        {
            "blueprint_run_id": blueprint_run_id,
            "summary_section": "join_path",
            "summary_metric": "join_step_count",
            "summary_value": len(join_path_rows),
            "notes": "Total fixed join/path blueprint steps.",
        }
    )
    for join_strength in JOIN_STRENGTHS:
        summary_rows.append(
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "join_path",
                "summary_metric": f"{join_strength}_count",
                "summary_value": sum(str(row["join_strength"]) == join_strength for row in join_path_rows),
                "notes": f"Join steps marked as {join_strength}.",
            }
        )
    for link_timing in LINK_TIMINGS:
        summary_rows.append(
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "join_path",
                "summary_metric": f"{link_timing}_count",
                "summary_value": sum(str(row["link_timing"]) == link_timing for row in join_path_rows),
                "notes": f"Join steps marked as {link_timing}.",
            }
        )

    summary_rows.extend(
        [
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "ambiguities",
                "summary_metric": "ambiguity_count",
                "summary_value": len(ambiguity_rows),
                "notes": "Total unresolved issues recorded by the blueprint.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "ambiguities",
                "summary_metric": "blocking_ambiguity_count",
                "summary_value": sum(bool(row["blocks_final_cohort_construction"]) for row in ambiguity_rows),
                "notes": "Ambiguities currently marked as blocking final cohort construction.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "upstream_runs",
                "summary_metric": "clinical_shortlist_run_id",
                "summary_value": str(workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
                "notes": "Current clinical shortlist run used by the blueprint.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "upstream_runs",
                "summary_metric": "clinical_core_audit_run_id",
                "summary_value": str(workflow_inputs.core_audit_latest_pointer["audit_run_id"]),
                "notes": "Clinical core audit run referenced by the current shortlist.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "upstream_runs",
                "summary_metric": "clinical_parse_run_id",
                "summary_value": str(workflow_inputs.core_audit_latest_pointer["parse_run_id"]),
                "notes": "Clinical biotab parse run behind the current shortlist and endpoint layers.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "upstream_runs",
                "summary_metric": "endpoint_crosswalk_run_id",
                "summary_value": str(workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]),
                "notes": "Endpoint crosswalk run used by the blueprint.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "upstream_runs",
                "summary_metric": "biospecimen_crosswalk_run_id",
                "summary_value": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]),
                "notes": "Biospecimen identifier crosswalk run used by the blueprint.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "upstream_runs",
                "summary_metric": "biospecimen_parse_run_id",
                "summary_value": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["parse_run_id"]),
                "notes": "Biospecimen biotab parse run behind the biospecimen crosswalk.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "upstream_runs",
                "summary_metric": "shared_source_run_id_count",
                "summary_value": len(source_run_id_set),
                "notes": "Count of unique source_run_id values across the clinical and biospecimen upstream layers.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "count_signals",
                "summary_metric": "clinical_patient_case_count",
                "summary_value": next(
                    table["row_count"]
                    for table in workflow_inputs.clinical_biotab_run_log["tables"]
                    if str(table.get("table_name") or "") == "clinical_patient"
                ),
                "notes": "clinical_patient row count from the saved clinical biotab parse run log.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "count_signals",
                "summary_metric": "clinical_source_case_submitter_coverage",
                "summary_value": clinical_source_case_coverage,
                "notes": f"Unique clinical source case coverage computed from {source_metadata_file_count} source-metadata rows.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "count_signals",
                "summary_metric": "biospecimen_patient_uuid_distinct_count",
                "summary_value": biospecimen_by_key_for_summary(workflow_inputs, ("biospecimen_sample", "bcr_patient_uuid"), "distinct_non_missing_count"),
                "notes": "Distinct patient UUID count from the biospecimen sample anchor row.",
            },
            {
                "blueprint_run_id": blueprint_run_id,
                "summary_section": "count_signals",
                "summary_metric": "biospecimen_sample_anchor_count",
                "summary_value": biospecimen_by_key_for_summary(workflow_inputs, ("biospecimen_sample", "bcr_sample_barcode"), "non_missing_count"),
                "notes": "Sample-anchor row count from the biospecimen crosswalk.",
            },
        ]
    )
    return summary_rows


def validate_field_assignments(
    required_rows: list[dict[str, Any]],
    optional_rows: list[dict[str, Any]],
    deferred_rows: list[dict[str, Any]],
) -> bool:
    seen: set[tuple[str, str, str]] = set()
    for row_set in (required_rows, optional_rows, deferred_rows):
        for row in row_set:
            key = (str(row["source_layer"]), str(row["table_name"]), str(row["field_name"]))
            if key in seen:
                return False
            seen.add(key)
    return True


def build_latest_pointer_payload(
    blueprint_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "blueprint_run_id": blueprint_run_id,
        "shortlist_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(workflow_inputs.core_audit_latest_pointer["audit_run_id"]),
        "clinical_parse_run_id": str(workflow_inputs.core_audit_latest_pointer["parse_run_id"]),
        "endpoint_crosswalk_run_id": str(workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["parse_run_id"]),
        "clinical_source_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["source_run_id"]),
        "biospecimen_source_run_id": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["source_run_id"]),
        "blueprint_run_directory": repo_relative(output_paths["blueprint_run_directory"], paths.repo_root),
        "cohort_blueprint_required_fields_tsv": repo_relative(
            output_paths["required_fields_tsv"], paths.repo_root
        ),
        "cohort_blueprint_optional_fields_tsv": repo_relative(
            output_paths["optional_fields_tsv"], paths.repo_root
        ),
        "cohort_blueprint_deferred_fields_tsv": repo_relative(
            output_paths["deferred_fields_tsv"], paths.repo_root
        ),
        "cohort_blueprint_join_path_tsv": repo_relative(output_paths["join_path_tsv"], paths.repo_root),
        "cohort_blueprint_ambiguities_tsv": repo_relative(
            output_paths["ambiguities_tsv"], paths.repo_root
        ),
        "cohort_blueprint_summary_tsv": repo_relative(output_paths["summary_tsv"], paths.repo_root),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "clinical_shortlist_latest_json": repo_relative(
            paths.clinical_shortlist_latest_pointer, paths.repo_root
        ),
        "endpoint_crosswalk_latest_json": repo_relative(
            paths.endpoint_crosswalk_latest_pointer, paths.repo_root
        ),
        "biospecimen_crosswalk_latest_json": repo_relative(
            paths.biospecimen_crosswalk_latest_pointer, paths.repo_root
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    blueprint_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    blueprint_run_dir = paths.blueprint_runs_root / blueprint_run_id
    run_log_path = blueprint_run_dir / "run_log.json"

    try:
        trial_config = load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths)
        create_run_directory(blueprint_run_dir)

        required_fields_path = blueprint_run_dir / "cohort_blueprint_required_fields.tsv"
        optional_fields_path = blueprint_run_dir / "cohort_blueprint_optional_fields.tsv"
        deferred_fields_path = blueprint_run_dir / "cohort_blueprint_deferred_fields.tsv"
        join_path_path = blueprint_run_dir / "cohort_blueprint_join_path.tsv"
        ambiguities_path = blueprint_run_dir / "cohort_blueprint_ambiguities.tsv"
        summary_path = blueprint_run_dir / "cohort_blueprint_summary.tsv"

        required_rows, optional_rows, deferred_rows = build_field_tables(
            blueprint_run_id=blueprint_run_id,
            workflow_inputs=workflow_inputs,
        )
        if not required_rows:
            raise CohortBlueprintError("Required cohort blueprint field rows were not generated.")
        if not optional_rows:
            raise CohortBlueprintError("Optional cohort blueprint field rows were not generated.")
        if not deferred_rows:
            raise CohortBlueprintError("Deferred cohort blueprint field rows were not generated.")

        join_path_rows = build_join_path_rows(blueprint_run_id, workflow_inputs)
        ambiguity_rows = build_ambiguity_rows(blueprint_run_id, workflow_inputs)
        summary_rows = build_summary_rows(
            blueprint_run_id=blueprint_run_id,
            workflow_inputs=workflow_inputs,
            required_rows=required_rows,
            optional_rows=optional_rows,
            deferred_rows=deferred_rows,
            join_path_rows=join_path_rows,
            ambiguity_rows=ambiguity_rows,
        )

        write_dict_rows_tsv(required_fields_path, FIELD_TABLE_FIELDNAMES, required_rows)
        write_dict_rows_tsv(optional_fields_path, FIELD_TABLE_FIELDNAMES, optional_rows)
        write_dict_rows_tsv(deferred_fields_path, FIELD_TABLE_FIELDNAMES, deferred_rows)
        write_dict_rows_tsv(join_path_path, JOIN_PATH_FIELDNAMES, join_path_rows)
        write_dict_rows_tsv(ambiguities_path, AMBIGUITY_FIELDNAMES, ambiguity_rows)
        write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        output_row_count_positive = all(
            len(rows) > 0
            for rows in (
                required_rows,
                optional_rows,
                deferred_rows,
                join_path_rows,
                ambiguity_rows,
                summary_rows,
            )
        )
        field_assignment_disjoint = validate_field_assignments(required_rows, optional_rows, deferred_rows)
        join_strength_values_valid = all(
            str(row["join_strength"]) in JOIN_STRENGTHS and str(row["link_timing"]) in LINK_TIMINGS
            for row in join_path_rows
        )
        ambiguity_severity_values_valid = all(
            str(row["severity"]) in SEVERITY_ORDER for row in ambiguity_rows
        )
        summary_rows_positive = len(summary_rows) > 0

        output_paths = {
            "blueprint_run_directory": blueprint_run_dir,
            "required_fields_tsv": required_fields_path,
            "optional_fields_tsv": optional_fields_path,
            "deferred_fields_tsv": deferred_fields_path,
            "join_path_tsv": join_path_path,
            "ambiguities_tsv": ambiguities_path,
            "summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            blueprint_run_id=blueprint_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "blueprint_run_id": blueprint_run_id,
            "shortlist_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(workflow_inputs.core_audit_latest_pointer["audit_run_id"]),
            "clinical_parse_run_id": str(workflow_inputs.core_audit_latest_pointer["parse_run_id"]),
            "endpoint_crosswalk_run_id": str(workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]),
            "biospecimen_crosswalk_run_id": str(
                workflow_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
            ),
            "biospecimen_parse_run_id": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["parse_run_id"]),
            "clinical_source_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["source_run_id"]),
            "biospecimen_source_run_id": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["source_run_id"]),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "results_root": repo_relative(paths.results_root, paths.repo_root),
                "clinical_shortlist_latest_json": repo_relative(
                    paths.clinical_shortlist_latest_pointer, paths.repo_root
                ),
                "clinical_shortlist_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_tsv"], paths.repo_root
                ),
                "clinical_shortlist_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_summary_tsv"], paths.repo_root
                ),
                "clinical_shortlist_by_table_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_by_table_tsv"], paths.repo_root
                ),
                "clinical_shortlist_run_log_json": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_run_log_json"], paths.repo_root
                ),
                "clinical_core_audit_latest_json": repo_relative(
                    workflow_inputs.input_paths["core_audit_latest_json"], paths.repo_root
                ),
                "clinical_core_field_audit_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_core_field_audit_tsv"], paths.repo_root
                ),
                "clinical_core_audit_run_log_json": repo_relative(
                    workflow_inputs.input_paths["core_audit_run_log_json"], paths.repo_root
                ),
                "clinical_biotab_run_log_json": repo_relative(
                    workflow_inputs.input_paths["clinical_biotab_run_log_json"], paths.repo_root
                ),
                "source_metadata_tsv": repo_relative(
                    workflow_inputs.input_paths["source_metadata_tsv"], paths.repo_root
                ),
                "source_run_log_json": repo_relative(
                    workflow_inputs.input_paths["source_run_log_json"], paths.repo_root
                ),
                "endpoint_crosswalk_latest_json": repo_relative(
                    paths.endpoint_crosswalk_latest_pointer, paths.repo_root
                ),
                "endpoint_candidate_inventory_tsv": repo_relative(
                    workflow_inputs.input_paths["endpoint_candidate_inventory_tsv"], paths.repo_root
                ),
                "endpoint_crosswalk_tsv": repo_relative(
                    workflow_inputs.input_paths["endpoint_crosswalk_tsv"], paths.repo_root
                ),
                "endpoint_crosswalk_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["endpoint_crosswalk_summary_tsv"], paths.repo_root
                ),
                "endpoint_crosswalk_run_log_json": repo_relative(
                    workflow_inputs.input_paths["endpoint_crosswalk_run_log_json"], paths.repo_root
                ),
                "biospecimen_crosswalk_latest_json": repo_relative(
                    paths.biospecimen_crosswalk_latest_pointer, paths.repo_root
                ),
                "biospecimen_identifier_inventory_tsv": repo_relative(
                    workflow_inputs.input_paths["biospecimen_identifier_inventory_tsv"], paths.repo_root
                ),
                "biospecimen_identifier_crosswalk_tsv": repo_relative(
                    workflow_inputs.input_paths["biospecimen_identifier_crosswalk_tsv"], paths.repo_root
                ),
                "biospecimen_identifier_crosswalk_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["biospecimen_identifier_crosswalk_summary_tsv"], paths.repo_root
                ),
                "biospecimen_identifier_crosswalk_run_log_json": repo_relative(
                    workflow_inputs.input_paths["biospecimen_identifier_crosswalk_run_log_json"],
                    paths.repo_root,
                ),
            },
            "outputs": {
                "blueprint_run_directory": repo_relative(blueprint_run_dir, paths.repo_root),
                "cohort_blueprint_required_fields_tsv": repo_relative(required_fields_path, paths.repo_root),
                "cohort_blueprint_optional_fields_tsv": repo_relative(optional_fields_path, paths.repo_root),
                "cohort_blueprint_deferred_fields_tsv": repo_relative(deferred_fields_path, paths.repo_root),
                "cohort_blueprint_join_path_tsv": repo_relative(join_path_path, paths.repo_root),
                "cohort_blueprint_ambiguities_tsv": repo_relative(ambiguities_path, paths.repo_root),
                "cohort_blueprint_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    output_row_count_positive
                    and field_assignment_disjoint
                    and join_strength_values_valid
                    and ambiguity_severity_values_valid
                    and summary_rows_positive
                ),
                "clinical_shortlist_latest_pointer_found": True,
                "endpoint_crosswalk_latest_pointer_found": True,
                "biospecimen_crosswalk_latest_pointer_found": True,
                "clinical_shortlist_run_log_completed": True,
                "clinical_core_audit_run_log_completed": True,
                "clinical_biotab_run_log_completed": True,
                "source_run_log_completed": True,
                "endpoint_crosswalk_run_log_completed": True,
                "biospecimen_crosswalk_run_log_completed": True,
                "required_fields_row_count_positive": len(required_rows) > 0,
                "optional_fields_row_count_positive": len(optional_rows) > 0,
                "deferred_fields_row_count_positive": len(deferred_rows) > 0,
                "join_path_row_count_positive": len(join_path_rows) > 0,
                "ambiguities_row_count_positive": len(ambiguity_rows) > 0,
                "summary_row_count_positive": len(summary_rows) > 0,
                "field_assignment_disjoint": field_assignment_disjoint,
                "join_strength_values_valid": join_strength_values_valid,
                "ambiguity_severity_values_valid": ambiguity_severity_values_valid,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "field_bucket_order": list(FIELD_BUCKET_ORDER),
                "field_bucket_precedence": list(FIELD_BUCKET_ORDER),
                "source_layer_order": list(SOURCE_LAYER_ORDER),
                "join_strengths": list(JOIN_STRENGTHS),
                "link_timings": list(LINK_TIMINGS),
                "severity_order": list(SEVERITY_ORDER),
                "proposed_unit_of_analysis": "patient/case",
                "proposed_biospecimen_anchor": "sample",
                "treatment_policy": "treatment_proxy_only_optional",
                "endpoint_policy": "candidate_join_preparation_only",
            },
            "counts": {
                "required_field_count": len(required_rows),
                "optional_field_count": len(optional_rows),
                "deferred_field_count": len(deferred_rows),
                "join_step_count": len(join_path_rows),
                "ambiguity_count": len(ambiguity_rows),
                "summary_count": len(summary_rows),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "clinical_shortlist_latest_pointer": workflow_inputs.clinical_shortlist_latest_pointer,
                "clinical_shortlist_run_log": workflow_inputs.clinical_shortlist_run_log,
                "core_audit_latest_pointer": workflow_inputs.core_audit_latest_pointer,
                "core_audit_run_log": workflow_inputs.core_audit_run_log,
                "clinical_biotab_run_log": workflow_inputs.clinical_biotab_run_log,
                "source_run_log": workflow_inputs.source_run_log,
                "endpoint_crosswalk_latest_pointer": workflow_inputs.endpoint_crosswalk_latest_pointer,
                "endpoint_crosswalk_run_log": workflow_inputs.endpoint_crosswalk_run_log,
                "biospecimen_crosswalk_latest_pointer": workflow_inputs.biospecimen_crosswalk_latest_pointer,
                "biospecimen_crosswalk_run_log": workflow_inputs.biospecimen_crosswalk_run_log,
            },
        }

        write_json(run_log_path, run_log_payload)
        write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "blueprint_run_id": blueprint_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_cohort_blueprint",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA cohort blueprint workflow complete.")
    print(f"Blueprint run ID: {run_log['blueprint_run_id']}")
    print(f"Clinical shortlist run ID: {run_log['shortlist_run_id']}")
    print(f"Endpoint crosswalk run ID: {run_log['endpoint_crosswalk_run_id']}")
    print(f"Biospecimen crosswalk run ID: {run_log['biospecimen_crosswalk_run_id']}")
    print(f"Blueprint output directory: {run_log['outputs']['blueprint_run_directory']}")
    print(f"Required fields TSV: {run_log['outputs']['cohort_blueprint_required_fields_tsv']}")
    print(f"Optional fields TSV: {run_log['outputs']['cohort_blueprint_optional_fields_tsv']}")
    print(f"Deferred fields TSV: {run_log['outputs']['cohort_blueprint_deferred_fields_tsv']}")
    print(f"Join path TSV: {run_log['outputs']['cohort_blueprint_join_path_tsv']}")
    print(f"Ambiguities TSV: {run_log['outputs']['cohort_blueprint_ambiguities_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['cohort_blueprint_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
