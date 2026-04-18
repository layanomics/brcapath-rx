#!/usr/bin/env python
"""Build a focused shortlist audit from the latest TCGA-BRCA clinical core field audit."""

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

SHORTLIST_BUCKETS = (
    "usable_baseline",
    "usable_treatment_proxy",
    "usable_endpoint_candidate",
    "weak_or_unusable",
    "unclear_manual_review",
)
MANUAL_REVIEW_PRIORITIES = ("high", "medium", "low")
SHORTLIST_PROFILE = "balanced"
NOTES_PLACEHOLDER = "[fill in during shortlist review]"
BASELINE_GROUPS = {
    "demographics",
    "diagnosis / pathology",
    "stage",
    "receptor / biomarker",
}
DIRECT_TREATMENT_TABLES = {"clinical_drug", "clinical_radiation"}
ENDPOINT_SIGNAL_PATTERNS = (
    "vital_status",
    "last_contact",
    "death",
    "tumor_status",
    "followup_lost_to",
    "days_to_patient_progression_free",
    "days_to_tumor_progression",
    "new_tumor_event",
)
TREATMENT_SIGNAL_PATTERNS = (
    "drug",
    "pharm",
    "therapy",
    "radiation",
    "regimen",
    "dose",
    "adjuvant",
    "started",
    "ended",
    "ongoing",
)
RESPONSE_SIGNAL_PATTERNS = ("response",)
SHORTLIST_FIELDNAMES = [
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
    "support_high_completeness",
    "support_high_missingness",
    "support_likely_high_value",
    "support_likely_weak",
    "shortlist_bucket",
    "shortlist_rule",
    "manual_review_priority",
    "notes_placeholder",
]
SHORTLIST_SUMMARY_FIELDNAMES = [
    "shortlist_run_id",
    "core_audit_run_id",
    "parse_run_id",
    "source_run_id",
    "shortlist_bucket",
    "field_count",
    "field_fraction",
    "high_priority_count",
    "medium_priority_count",
    "low_priority_count",
    "field_names_json",
]
SHORTLIST_BY_TABLE_FIELDNAMES = [
    "shortlist_run_id",
    "core_audit_run_id",
    "parse_run_id",
    "source_run_id",
    "table_name",
    "shortlist_bucket",
    "field_count",
    "field_fraction_of_table",
    "high_priority_count",
    "medium_priority_count",
    "low_priority_count",
    "field_names_json",
]
REVIEW_SUPPORT_TABLES = {
    "support_high_completeness": "21_clinical_core_high_completeness_fields.tsv",
    "support_high_missingness": "22_clinical_core_high_missingness_fields.tsv",
    "support_likely_high_value": "23_clinical_core_likely_high_value_fields.tsv",
    "support_likely_weak": "24_clinical_core_likely_weak_fields.tsv",
}
SHORTLIST_RULES = {
    "weak_all_missing",
    "weak_extreme_missingness",
    "weak_near_constant",
    "manual_identifier_admin",
    "manual_other_unclear_group",
    "manual_nte_prefixed",
    "manual_response_like_ambiguous",
    "manual_overlapping_candidate_signals",
    "endpoint_candidate_time_or_status_signal",
    "treatment_proxy_direct_treatment_table",
    "treatment_proxy_treatment_signal",
    "baseline_candidate_structured_core_field",
    "manual_default_fallback",
}


class ClinicalShortlistError(RuntimeError):
    """Raised when the clinical shortlist workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the clinical shortlist workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    core_audit_latest_pointer: Path
    shortlist_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved shortlist inputs from the latest core audit and saved review tables."""

    core_audit_latest_pointer: dict[str, Any]
    core_audit_run_log: dict[str, Any]
    shortlist_rows: list[dict[str, Any]]
    group_summary_rows: list[dict[str, Any]]
    missingness_rows: list[dict[str, Any]]
    support_key_sets: dict[str, set[tuple[str, str]]]
    input_paths: dict[str, Path]


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise ClinicalShortlistError("Unable to locate the repository root from the script path.")


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
                raise ClinicalShortlistError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise ClinicalShortlistError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise ClinicalShortlistError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise ClinicalShortlistError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise ClinicalShortlistError(
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
        raise ClinicalShortlistError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ClinicalShortlistError(f"Expected a JSON object in file: {path}")
    return data


def read_tsv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise ClinicalShortlistError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise ClinicalShortlistError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise ClinicalShortlistError(f"Run directory already exists: {path}")
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
        core_audit_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_core_field_audit_latest.json"
        ),
        shortlist_runs_root=audit_root / "tcga-brca" / "variables" / "clinical_shortlist_runs",
        latest_pointer=audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_shortlist_latest.json",
    )


def normalize_key(row: dict[str, Any]) -> tuple[str, str]:
    return str(row["table_name"]), str(row["field_name"])


def read_required_tsv(path: Path, label: str) -> list[dict[str, str]]:
    if not path.exists():
        raise ClinicalShortlistError(f"Required {label} does not exist: {path}")
    rows = read_tsv_dict_rows(path)
    if not rows:
        raise ClinicalShortlistError(f"Required {label} has no rows: {path}")
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
        raise ClinicalShortlistError(f"Duplicate field keys found in {label}: {sample}")


def parse_int(value: Any, label: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ClinicalShortlistError(f"Expected integer value for {label}: {value!r}") from exc


def parse_float(value: Any, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ClinicalShortlistError(f"Expected float value for {label}: {value!r}") from exc


def build_support_key_set(
    *,
    path: Path,
    label: str,
    audit_keys: set[tuple[str, str]],
) -> set[tuple[str, str]]:
    if not path.exists():
        raise ClinicalShortlistError(f"Required {label} does not exist: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if not {"table_name", "field_name"}.issubset(fieldnames):
            raise ClinicalShortlistError(f"Required {label} is missing table_name or field_name columns: {path}")
        rows = list(reader)

    support_keys = {normalize_key(row) for row in rows}
    missing_keys = support_keys.difference(audit_keys)
    if missing_keys:
        sample = sorted(f"{table_name}:{field_name}" for table_name, field_name in missing_keys)[:10]
        raise ClinicalShortlistError(f"{label} contains field keys not present in the combined audit TSV: {sample}")

    return support_keys


def load_workflow_inputs(paths: WorkflowPaths) -> WorkflowInputs:
    core_audit_latest_pointer = load_json(paths.core_audit_latest_pointer)
    required_pointer_keys = {
        "audit_run_id",
        "parse_run_id",
        "source_run_id",
        "clinical_core_field_audit_tsv",
        "clinical_core_field_group_summary_tsv",
        "clinical_core_missingness_summary_tsv",
        "run_log_json",
    }
    missing_pointer_keys = required_pointer_keys.difference(core_audit_latest_pointer.keys())
    if missing_pointer_keys:
        raise ClinicalShortlistError(
            f"Clinical core audit latest pointer is missing required keys {sorted(missing_pointer_keys)}: "
            f"{paths.core_audit_latest_pointer}"
        )

    input_paths = {
        "core_audit_latest_json": paths.core_audit_latest_pointer,
        "clinical_core_field_audit_tsv": paths.repo_root / core_audit_latest_pointer["clinical_core_field_audit_tsv"],
        "clinical_core_field_group_summary_tsv": (
            paths.repo_root / core_audit_latest_pointer["clinical_core_field_group_summary_tsv"]
        ),
        "clinical_core_missingness_summary_tsv": (
            paths.repo_root / core_audit_latest_pointer["clinical_core_missingness_summary_tsv"]
        ),
        "core_audit_run_log_json": paths.repo_root / core_audit_latest_pointer["run_log_json"],
    }
    for support_label, filename in REVIEW_SUPPORT_TABLES.items():
        input_paths[support_label] = paths.results_root / filename

    shortlist_rows_raw = read_required_tsv(input_paths["clinical_core_field_audit_tsv"], "combined core field audit TSV")
    group_summary_rows = read_required_tsv(
        input_paths["clinical_core_field_group_summary_tsv"],
        "core field-group summary TSV",
    )
    missingness_rows = read_required_tsv(
        input_paths["clinical_core_missingness_summary_tsv"],
        "core missingness summary TSV",
    )
    core_audit_run_log = load_json(input_paths["core_audit_run_log_json"])
    if core_audit_run_log.get("status") != "completed":
        raise ClinicalShortlistError(
            f"Core audit run log is not marked completed: {input_paths['core_audit_run_log_json']}"
        )

    validate_unique_field_keys(shortlist_rows_raw, "combined core field audit TSV")
    validate_unique_field_keys(missingness_rows, "core missingness summary TSV")

    if len(shortlist_rows_raw) != len(missingness_rows):
        raise ClinicalShortlistError(
            "Combined core field audit row count does not match core missingness summary row count: "
            f"{len(shortlist_rows_raw)} vs {len(missingness_rows)}"
        )

    shortlist_rows: list[dict[str, Any]] = []
    for row in shortlist_rows_raw:
        shortlist_rows.append(
            {
                **row,
                "source_position": parse_int(row["source_position"], "source_position"),
                "missing_like_fraction": parse_float(row["missing_like_fraction"], "missing_like_fraction"),
                "non_missing_count": parse_int(row["non_missing_count"], "non_missing_count"),
                "distinct_non_missing_count": parse_int(
                    row["distinct_non_missing_count"],
                    "distinct_non_missing_count",
                ),
            }
        )

    audit_keys = {normalize_key(row) for row in shortlist_rows}
    support_key_sets = {
        support_label: build_support_key_set(
            path=input_paths[support_label],
            label=support_label,
            audit_keys=audit_keys,
        )
        for support_label in REVIEW_SUPPORT_TABLES
    }

    expected_field_count = parse_int(core_audit_latest_pointer.get("field_count"), "core audit pointer field_count")
    if expected_field_count != len(shortlist_rows):
        raise ClinicalShortlistError(
            f"Combined core audit TSV row count {len(shortlist_rows)} does not match pointer field_count {expected_field_count}."
        )

    return WorkflowInputs(
        core_audit_latest_pointer=core_audit_latest_pointer,
        core_audit_run_log=core_audit_run_log,
        shortlist_rows=shortlist_rows,
        group_summary_rows=group_summary_rows,
        missingness_rows=missingness_rows,
        support_key_sets=support_key_sets,
        input_paths=input_paths,
    )


def contains_any_pattern(value: str, patterns: tuple[str, ...]) -> bool:
    return any(pattern in value for pattern in patterns)


def is_response_like(field_name: str) -> bool:
    return contains_any_pattern(field_name, RESPONSE_SIGNAL_PATTERNS)


def is_endpoint_like(field_name: str) -> bool:
    return contains_any_pattern(field_name, ENDPOINT_SIGNAL_PATTERNS)


def is_treatment_like(field_name: str) -> bool:
    return contains_any_pattern(field_name, TREATMENT_SIGNAL_PATTERNS)


def determine_bucket_and_rule(
    *,
    row: dict[str, Any],
    support_flags: dict[str, bool],
) -> tuple[str, str]:
    field_name_lower = str(row["field_name"]).lower()
    probable_field_group = str(row["probable_field_group"])
    missing_like_fraction = float(row["missing_like_fraction"])
    non_missing_count = int(row["non_missing_count"])
    distinct_non_missing_count = int(row["distinct_non_missing_count"])
    table_name = str(row["table_name"])

    strong_enough_for_use = (
        missing_like_fraction <= 0.25
        and non_missing_count > 0
        and distinct_non_missing_count >= 2
        and not support_flags["support_likely_weak"]
    )
    field_is_identifier_admin = probable_field_group == "identifier / admin"
    field_is_other_unclear = probable_field_group == "other / unclear"
    field_is_nte_prefixed = field_name_lower.startswith("nte_")
    field_is_response_like = is_response_like(field_name_lower)
    field_is_endpoint_like = is_endpoint_like(field_name_lower)
    field_is_treatment_like = is_treatment_like(field_name_lower)
    field_is_direct_treatment_table = table_name in DIRECT_TREATMENT_TABLES
    field_is_baseline_group = probable_field_group in BASELINE_GROUPS

    if non_missing_count == 0:
        return "weak_or_unusable", "weak_all_missing"
    if missing_like_fraction >= 0.90:
        return "weak_or_unusable", "weak_extreme_missingness"
    if distinct_non_missing_count <= 1 or support_flags["support_likely_weak"]:
        return "weak_or_unusable", "weak_near_constant"

    if field_is_identifier_admin:
        return "unclear_manual_review", "manual_identifier_admin"
    if field_is_other_unclear:
        return "unclear_manual_review", "manual_other_unclear_group"
    if field_is_nte_prefixed:
        return "unclear_manual_review", "manual_nte_prefixed"
    if field_is_response_like:
        return "unclear_manual_review", "manual_response_like_ambiguous"

    raw_candidate_matches = 0
    if strong_enough_for_use and field_is_endpoint_like:
        raw_candidate_matches += 1
    if strong_enough_for_use and (field_is_direct_treatment_table or field_is_treatment_like):
        raw_candidate_matches += 1
    if strong_enough_for_use and field_is_baseline_group:
        raw_candidate_matches += 1
    if raw_candidate_matches > 1:
        return "unclear_manual_review", "manual_overlapping_candidate_signals"

    if strong_enough_for_use and field_is_endpoint_like:
        return "usable_endpoint_candidate", "endpoint_candidate_time_or_status_signal"

    if strong_enough_for_use and field_is_direct_treatment_table and not field_is_endpoint_like:
        return "usable_treatment_proxy", "treatment_proxy_direct_treatment_table"

    if strong_enough_for_use and field_is_treatment_like and not field_is_endpoint_like:
        return "usable_treatment_proxy", "treatment_proxy_treatment_signal"

    if (
        strong_enough_for_use
        and field_is_baseline_group
        and not field_is_endpoint_like
        and not field_is_treatment_like
        and not field_is_nte_prefixed
    ):
        return "usable_baseline", "baseline_candidate_structured_core_field"

    return "unclear_manual_review", "manual_default_fallback"


def assign_manual_review_priority(shortlist_bucket: str, missing_like_fraction: float) -> str:
    if shortlist_bucket in {"unclear_manual_review", "usable_endpoint_candidate"}:
        return "high"
    if shortlist_bucket == "usable_treatment_proxy":
        return "medium"
    if shortlist_bucket == "usable_baseline" and missing_like_fraction > 0.10:
        return "medium"
    return "low"


def build_shortlist_rows(
    *,
    shortlist_run_id: str,
    workflow_inputs: WorkflowInputs,
) -> list[dict[str, Any]]:
    core_audit_run_id = str(workflow_inputs.core_audit_latest_pointer["audit_run_id"])
    parse_run_id = str(workflow_inputs.core_audit_latest_pointer["parse_run_id"])
    source_run_id = str(workflow_inputs.core_audit_latest_pointer["source_run_id"])

    shortlist_rows: list[dict[str, Any]] = []
    for row in sorted(
        workflow_inputs.shortlist_rows,
        key=lambda item: (str(item["table_name"]), int(item["source_position"])),
    ):
        key = normalize_key(row)
        support_flags = {
            support_label: key in support_keys
            for support_label, support_keys in workflow_inputs.support_key_sets.items()
        }
        shortlist_bucket, shortlist_rule = determine_bucket_and_rule(row=row, support_flags=support_flags)
        if shortlist_rule not in SHORTLIST_RULES:
            raise ClinicalShortlistError(f"Unexpected shortlist rule produced for {key}: {shortlist_rule}")

        manual_review_priority = assign_manual_review_priority(
            shortlist_bucket=shortlist_bucket,
            missing_like_fraction=float(row["missing_like_fraction"]),
        )
        if manual_review_priority not in MANUAL_REVIEW_PRIORITIES:
            raise ClinicalShortlistError(
                f"Unexpected manual_review_priority produced for {key}: {manual_review_priority}"
            )

        shortlist_rows.append(
            {
                "shortlist_run_id": shortlist_run_id,
                "core_audit_run_id": core_audit_run_id,
                "parse_run_id": parse_run_id,
                "source_run_id": source_run_id,
                "table_name": row["table_name"],
                "field_name": row["field_name"],
                "source_position": row["source_position"],
                "alternate_column_name": row["alternate_column_name"],
                "probable_field_group": row["probable_field_group"],
                "missing_like_fraction": round(float(row["missing_like_fraction"]), 6),
                "non_missing_count": int(row["non_missing_count"]),
                "distinct_non_missing_count": int(row["distinct_non_missing_count"]),
                "example_values_small_sample": row["example_values_small_sample"],
                "support_high_completeness": support_flags["support_high_completeness"],
                "support_high_missingness": support_flags["support_high_missingness"],
                "support_likely_high_value": support_flags["support_likely_high_value"],
                "support_likely_weak": support_flags["support_likely_weak"],
                "shortlist_bucket": shortlist_bucket,
                "shortlist_rule": shortlist_rule,
                "manual_review_priority": manual_review_priority,
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )

    return shortlist_rows


def build_shortlist_summary_rows(
    *,
    shortlist_run_id: str,
    core_audit_run_id: str,
    parse_run_id: str,
    source_run_id: str,
    shortlist_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    bucket_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in shortlist_rows:
        bucket_rows[str(row["shortlist_bucket"])].append(row)

    total_fields = len(shortlist_rows)
    summary_rows: list[dict[str, Any]] = []
    for shortlist_bucket in SHORTLIST_BUCKETS:
        rows = sorted(
            bucket_rows.get(shortlist_bucket, []),
            key=lambda item: (str(item["table_name"]), int(item["source_position"])),
        )
        field_count = len(rows)
        summary_rows.append(
            {
                "shortlist_run_id": shortlist_run_id,
                "core_audit_run_id": core_audit_run_id,
                "parse_run_id": parse_run_id,
                "source_run_id": source_run_id,
                "shortlist_bucket": shortlist_bucket,
                "field_count": field_count,
                "field_fraction": round(field_count / total_fields, 6) if total_fields else 0.0,
                "high_priority_count": sum(row["manual_review_priority"] == "high" for row in rows),
                "medium_priority_count": sum(row["manual_review_priority"] == "medium" for row in rows),
                "low_priority_count": sum(row["manual_review_priority"] == "low" for row in rows),
                "field_names_json": json.dumps([row["field_name"] for row in rows], ensure_ascii=True),
            }
        )

    return summary_rows


def build_shortlist_by_table_rows(
    *,
    shortlist_run_id: str,
    core_audit_run_id: str,
    parse_run_id: str,
    source_run_id: str,
    shortlist_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows_by_table_bucket: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    table_field_counts: dict[str, int] = defaultdict(int)

    for row in shortlist_rows:
        table_name = str(row["table_name"])
        shortlist_bucket = str(row["shortlist_bucket"])
        rows_by_table_bucket[(table_name, shortlist_bucket)].append(row)
        table_field_counts[table_name] += 1

    shortlist_by_table_rows: list[dict[str, Any]] = []
    for table_name in sorted(table_field_counts):
        total_fields = table_field_counts[table_name]
        for shortlist_bucket in SHORTLIST_BUCKETS:
            rows = sorted(
                rows_by_table_bucket.get((table_name, shortlist_bucket), []),
                key=lambda item: int(item["source_position"]),
            )
            field_count = len(rows)
            shortlist_by_table_rows.append(
                {
                    "shortlist_run_id": shortlist_run_id,
                    "core_audit_run_id": core_audit_run_id,
                    "parse_run_id": parse_run_id,
                    "source_run_id": source_run_id,
                    "table_name": table_name,
                    "shortlist_bucket": shortlist_bucket,
                    "field_count": field_count,
                    "field_fraction_of_table": round(field_count / total_fields, 6) if total_fields else 0.0,
                    "high_priority_count": sum(row["manual_review_priority"] == "high" for row in rows),
                    "medium_priority_count": sum(row["manual_review_priority"] == "medium" for row in rows),
                    "low_priority_count": sum(row["manual_review_priority"] == "low" for row in rows),
                    "field_names_json": json.dumps([row["field_name"] for row in rows], ensure_ascii=True),
                }
            )

    return shortlist_by_table_rows


def write_failure_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    shortlist_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    shortlist_run_dir = paths.shortlist_runs_root / shortlist_run_id
    run_log_path = shortlist_run_dir / "run_log.json"

    try:
        trial_config = load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths)
        create_run_directory(shortlist_run_dir)

        shortlist_rows = build_shortlist_rows(
            shortlist_run_id=shortlist_run_id,
            workflow_inputs=workflow_inputs,
        )
        if len(shortlist_rows) != len(workflow_inputs.shortlist_rows):
            raise ClinicalShortlistError(
                f"Shortlist row count {len(shortlist_rows)} does not match core audit row count {len(workflow_inputs.shortlist_rows)}."
            )

        core_audit_run_id = str(workflow_inputs.core_audit_latest_pointer["audit_run_id"])
        parse_run_id = str(workflow_inputs.core_audit_latest_pointer["parse_run_id"])
        source_run_id = str(workflow_inputs.core_audit_latest_pointer["source_run_id"])

        shortlist_path = shortlist_run_dir / "clinical_shortlist.tsv"
        write_dict_rows_tsv(shortlist_path, SHORTLIST_FIELDNAMES, shortlist_rows)

        shortlist_summary_rows = build_shortlist_summary_rows(
            shortlist_run_id=shortlist_run_id,
            core_audit_run_id=core_audit_run_id,
            parse_run_id=parse_run_id,
            source_run_id=source_run_id,
            shortlist_rows=shortlist_rows,
        )
        shortlist_summary_path = shortlist_run_dir / "clinical_shortlist_summary.tsv"
        write_dict_rows_tsv(shortlist_summary_path, SHORTLIST_SUMMARY_FIELDNAMES, shortlist_summary_rows)

        shortlist_by_table_rows = build_shortlist_by_table_rows(
            shortlist_run_id=shortlist_run_id,
            core_audit_run_id=core_audit_run_id,
            parse_run_id=parse_run_id,
            source_run_id=source_run_id,
            shortlist_rows=shortlist_rows,
        )
        shortlist_by_table_path = shortlist_run_dir / "clinical_shortlist_by_table.tsv"
        write_dict_rows_tsv(shortlist_by_table_path, SHORTLIST_BY_TABLE_FIELDNAMES, shortlist_by_table_rows)

        bucket_validation = all(
            str(row["shortlist_bucket"]) in SHORTLIST_BUCKETS and str(row["shortlist_bucket"]).strip()
            for row in shortlist_rows
        )
        summary_total = sum(parse_int(row["field_count"], "summary field_count") for row in shortlist_summary_rows)
        by_table_totals_match = True
        counts_by_table_bucket: dict[tuple[str, str], int] = {
            (str(row["table_name"]), str(row["shortlist_bucket"])): parse_int(row["field_count"], "by-table field_count")
            for row in shortlist_by_table_rows
        }
        audit_table_counts = defaultdict(int)
        for row in shortlist_rows:
            audit_table_counts[str(row["table_name"])] += 1
        for table_name, total_count in audit_table_counts.items():
            counted = sum(counts_by_table_bucket[(table_name, bucket)] for bucket in SHORTLIST_BUCKETS)
            if counted != total_count:
                by_table_totals_match = False
                break

        latest_pointer_payload = {
            "updated_at_utc": format_utc_timestamp(utc_now()),
            "shortlist_run_id": shortlist_run_id,
            "core_audit_run_id": core_audit_run_id,
            "parse_run_id": parse_run_id,
            "source_run_id": source_run_id,
            "shortlist_profile": SHORTLIST_PROFILE,
            "shortlist_run_directory": repo_relative(shortlist_run_dir, paths.repo_root),
            "clinical_shortlist_tsv": repo_relative(shortlist_path, paths.repo_root),
            "clinical_shortlist_summary_tsv": repo_relative(shortlist_summary_path, paths.repo_root),
            "clinical_shortlist_by_table_tsv": repo_relative(shortlist_by_table_path, paths.repo_root),
            "run_log_json": repo_relative(run_log_path, paths.repo_root),
            "core_audit_latest_json": repo_relative(paths.core_audit_latest_pointer, paths.repo_root),
            "field_count": len(shortlist_rows),
        }

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
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
                "core_audit_latest_json": repo_relative(paths.core_audit_latest_pointer, paths.repo_root),
                "clinical_core_field_audit_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_core_field_audit_tsv"],
                    paths.repo_root,
                ),
                "clinical_core_field_group_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_core_field_group_summary_tsv"],
                    paths.repo_root,
                ),
                "clinical_core_missingness_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_core_missingness_summary_tsv"],
                    paths.repo_root,
                ),
                "core_audit_run_log_json": repo_relative(
                    workflow_inputs.input_paths["core_audit_run_log_json"],
                    paths.repo_root,
                ),
                "review_support_tables": {
                    support_label: repo_relative(path, paths.repo_root)
                    for support_label, path in workflow_inputs.input_paths.items()
                    if support_label in REVIEW_SUPPORT_TABLES
                },
            },
            "outputs": {
                "shortlist_run_directory": repo_relative(shortlist_run_dir, paths.repo_root),
                "clinical_shortlist_tsv": repo_relative(shortlist_path, paths.repo_root),
                "clinical_shortlist_summary_tsv": repo_relative(shortlist_summary_path, paths.repo_root),
                "clinical_shortlist_by_table_tsv": repo_relative(shortlist_by_table_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    len(shortlist_rows) == len(workflow_inputs.shortlist_rows)
                    and bucket_validation
                    and summary_total == len(shortlist_rows)
                    and by_table_totals_match
                ),
                "input_audit_tables_found": True,
                "input_group_summary_found": True,
                "input_missingness_summary_found": True,
                "input_review_tables_found": True,
                "output_row_count_matches_core_audit": len(shortlist_rows) == len(workflow_inputs.shortlist_rows),
                "every_field_has_exactly_one_shortlist_bucket": bucket_validation,
                "summary_bucket_counts_match_total": summary_total == len(shortlist_rows),
                "by_table_bucket_counts_match_table_totals": by_table_totals_match,
                "no_prior_run_overwrite": True,
            },
            "rules": {
                "shortlist_profile": SHORTLIST_PROFILE,
                "bucket_precedence": [
                    "weak_or_unusable",
                    "unclear_manual_review",
                    "usable_endpoint_candidate",
                    "usable_treatment_proxy",
                    "usable_baseline",
                    "unclear_manual_review",
                ],
                "strong_enough_for_use": {
                    "missing_like_fraction_max": 0.25,
                    "non_missing_count_min": 1,
                    "distinct_non_missing_count_min": 2,
                    "exclude_support_likely_weak": True,
                },
                "baseline_groups": sorted(BASELINE_GROUPS),
                "direct_treatment_tables": sorted(DIRECT_TREATMENT_TABLES),
                "endpoint_signal_patterns": list(ENDPOINT_SIGNAL_PATTERNS),
                "treatment_signal_patterns": list(TREATMENT_SIGNAL_PATTERNS),
                "response_signal_patterns": list(RESPONSE_SIGNAL_PATTERNS),
                "manual_review_priorities": list(MANUAL_REVIEW_PRIORITIES),
                "shortlist_rules": sorted(SHORTLIST_RULES),
                "notes_placeholder": NOTES_PLACEHOLDER,
            },
            "support_key_counts": {
                support_label: len(support_keys)
                for support_label, support_keys in workflow_inputs.support_key_sets.items()
            },
            "bucket_counts": {
                row["shortlist_bucket"]: parse_int(row["field_count"], "summary field_count")
                for row in shortlist_summary_rows
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "core_audit_latest_pointer": workflow_inputs.core_audit_latest_pointer,
                "core_audit_run_log": workflow_inputs.core_audit_run_log,
            },
        }

        write_json(run_log_path, run_log_payload)
        write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "shortlist_run_id": shortlist_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_clinical_shortlist",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA clinical shortlist audit complete.")
    print(f"Shortlist run ID: {run_log['shortlist_run_id']}")
    print(f"Core audit run ID: {run_log['core_audit_run_id']}")
    print(f"Parse run ID: {run_log['parse_run_id']}")
    print(f"Source run ID: {run_log['source_run_id']}")
    print(f"Shortlist output directory: {run_log['outputs']['shortlist_run_directory']}")
    print(f"Shortlist TSV: {run_log['outputs']['clinical_shortlist_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['clinical_shortlist_summary_tsv']}")
    print(f"By-table TSV: {run_log['outputs']['clinical_shortlist_by_table_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
