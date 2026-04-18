#!/usr/bin/env python
"""Build an auditable ambiguity-resolution layer from the saved TCGA-BRCA cohort blueprint."""

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

BLOCKER_CLASS_ORDER = (
    "blocking_now",
    "review_before_build",
    "non_blocking_for_minimal_build",
)
RESOLUTION_PRIORITY_ORDER = ("high", "medium", "low")
EVIDENCE_SUFFICIENCY_ORDER = (
    "enough_for_provisional_rule",
    "partial_needs_manual_review",
    "insufficient_needs_new_evidence",
)
RECOMMENDED_ACTION_TYPE_ORDER = (
    "provisional_rule",
    "manual_review",
    "later_xml_validation",
    "defer_from_minimal_build",
    "safe_to_ignore_for_now",
)
ACTION_ORDER = ("AR01", "AR02", "AR03", "AR04", "AR05", "AR06")
NOTES_PLACEHOLDER = "[fill in during blueprint ambiguity-resolution review]"
SUMMARY_FIELDNAMES = [
    "ambiguity_resolution_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
INVENTORY_FIELDNAMES = [
    "ambiguity_resolution_run_id",
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
    "evidence_still_missing",
    "blocker_class",
    "resolution_priority",
    "evidence_sufficiency",
    "recommended_next_action",
    "recommended_action_type",
    "linked_action_id",
    "notes_placeholder",
]
ACTIONS_FIELDNAMES = [
    "ambiguity_resolution_run_id",
    "action_id",
    "action_title",
    "action_scope",
    "linked_ambiguity_ids_json",
    "rationale",
    "provisional_rule_if_any",
    "evidence_basis",
    "risk_if_wrong",
    "needed_before_minimal_build",
    "notes_placeholder",
]
PRIORITY_FIELDNAMES = [
    "ambiguity_resolution_run_id",
    "priority_rank",
    "action_id",
    "action_title",
    "blocker_class",
    "resolution_priority",
    "recommended_action_type",
    "linked_ambiguity_count",
    "linked_high_severity_count",
    "linked_ambiguity_categories_json",
    "needed_before_minimal_build",
]


class BlueprintAmbiguityResolutionError(RuntimeError):
    """Raised when the blueprint ambiguity-resolution workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the ambiguity-resolution workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    blueprint_latest_pointer: Path
    clinical_shortlist_latest_pointer: Path
    endpoint_crosswalk_latest_pointer: Path
    biospecimen_crosswalk_latest_pointer: Path
    ambiguity_resolution_runs_root: Path
    latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved ambiguity-resolution inputs loaded from saved audit layers."""

    blueprint_latest_pointer: dict[str, Any]
    blueprint_run_log: dict[str, Any]
    clinical_shortlist_latest_pointer: dict[str, Any]
    clinical_shortlist_run_log: dict[str, Any]
    endpoint_crosswalk_latest_pointer: dict[str, Any]
    endpoint_crosswalk_run_log: dict[str, Any]
    biospecimen_crosswalk_latest_pointer: dict[str, Any]
    biospecimen_crosswalk_run_log: dict[str, Any]
    blueprint_required_rows: list[dict[str, str]]
    blueprint_optional_rows: list[dict[str, str]]
    blueprint_deferred_rows: list[dict[str, str]]
    blueprint_join_path_rows: list[dict[str, str]]
    blueprint_ambiguity_rows: list[dict[str, str]]
    blueprint_summary_rows: list[dict[str, str]]
    clinical_shortlist_rows: list[dict[str, str]]
    clinical_shortlist_summary_rows: list[dict[str, str]]
    clinical_shortlist_by_table_rows: list[dict[str, str]]
    endpoint_inventory_rows: list[dict[str, str]]
    endpoint_crosswalk_rows: list[dict[str, str]]
    endpoint_summary_rows: list[dict[str, str]]
    biospecimen_inventory_rows: list[dict[str, str]]
    biospecimen_crosswalk_rows: list[dict[str, str]]
    biospecimen_summary_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class ClassificationRule:
    """Deterministic rule applied to one blueprint ambiguity row."""

    blocker_class: str
    resolution_priority: str
    evidence_sufficiency: str
    recommended_action_type: str
    linked_action_id: str
    evidence_still_missing: str


@dataclass(frozen=True)
class ActionTemplate:
    """Static action-planning template used for grouped ambiguity actions."""

    action_title: str
    action_scope: str
    rationale: str
    provisional_rule_if_any: str
    risk_if_wrong: str
    needed_before_minimal_build: bool


CLASSIFICATION_RULES: dict[str, ClassificationRule] = {
    "A01": ClassificationRule(
        blocker_class="review_before_build",
        resolution_priority="medium",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="provisional_rule",
        linked_action_id="AR03",
        evidence_still_missing=(
            "A saved source-table precedence or row-level reconciliation rule between patient and follow-up vital_status is still missing."
        ),
    ),
    "A02": ClassificationRule(
        blocker_class="review_before_build",
        resolution_priority="medium",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="provisional_rule",
        linked_action_id="AR03",
        evidence_still_missing=(
            "A saved source-table precedence or row-level reconciliation rule between patient and follow-up last_contact_days_to is still missing."
        ),
    ),
    "A03": ClassificationRule(
        blocker_class="review_before_build",
        resolution_priority="medium",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="provisional_rule",
        linked_action_id="AR03",
        evidence_still_missing=(
            "A saved source-table precedence or row-level reconciliation rule between patient and follow-up tumor_status is still missing."
        ),
    ),
    "A04": ClassificationRule(
        blocker_class="review_before_build",
        resolution_priority="medium",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="provisional_rule",
        linked_action_id="AR03",
        evidence_still_missing=(
            "A saved source-table precedence or row-level reconciliation rule between patient and follow-up new_tumor_event_dx_indicator is still missing."
        ),
    ),
    "A05": ClassificationRule(
        blocker_class="review_before_build",
        resolution_priority="medium",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="defer_from_minimal_build",
        linked_action_id="AR04",
        evidence_still_missing=(
            "No sufficiently complete saved death timing source is available yet for endpoint freeze from disk-only audit outputs."
        ),
    ),
    "A06": ClassificationRule(
        blocker_class="review_before_build",
        resolution_priority="medium",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="defer_from_minimal_build",
        linked_action_id="AR04",
        evidence_still_missing=(
            "No non-missing saved progression timing fields are available yet for endpoint freeze from disk-only audit outputs."
        ),
    ),
    "A07": ClassificationRule(
        blocker_class="blocking_now",
        resolution_priority="high",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="provisional_rule",
        linked_action_id="AR01",
        evidence_still_missing=(
            "A saved cross-layer canonical patient identifier rule is still missing, so barcode and UUID evidence must stay parallel rather than collapsed."
        ),
    ),
    "A08": ClassificationRule(
        blocker_class="blocking_now",
        resolution_priority="high",
        evidence_sufficiency="partial_needs_manual_review",
        recommended_action_type="manual_review",
        linked_action_id="AR02",
        evidence_still_missing=(
            "A saved row-level explanation for the 1097 clinical / 1098 source / 1101 biospecimen patient-count gap is still missing."
        ),
    ),
    "A09": ClassificationRule(
        blocker_class="review_before_build",
        resolution_priority="medium",
        evidence_sufficiency="enough_for_provisional_rule",
        recommended_action_type="defer_from_minimal_build",
        linked_action_id="AR05",
        evidence_still_missing=(
            "A saved patient-level aggregation or prioritization rule for one-to-many treatment rows is still missing."
        ),
    ),
    "A10": ClassificationRule(
        blocker_class="non_blocking_for_minimal_build",
        resolution_priority="low",
        evidence_sufficiency="insufficient_needs_new_evidence",
        recommended_action_type="later_xml_validation",
        linked_action_id="AR06",
        evidence_still_missing=(
            "No saved direct analyte identifier appears in the slide layer, so any deeper child-layer validation would need later source expansion."
        ),
    ),
    "A11": ClassificationRule(
        blocker_class="non_blocking_for_minimal_build",
        resolution_priority="low",
        evidence_sufficiency="insufficient_needs_new_evidence",
        recommended_action_type="later_xml_validation",
        linked_action_id="AR06",
        evidence_still_missing=(
            "No saved direct analyte identifier appears in the aliquot layer, so any deeper child-layer validation would need later source expansion."
        ),
    ),
}
ACTION_TEMPLATES: dict[str, ActionTemplate] = {
    "AR01": ActionTemplate(
        action_title="Define provisional patient/case join rule",
        action_scope="clinical_patient + biospecimen sample anchor",
        rationale=(
            "The saved clinical and biospecimen layers expose strong patient identifiers, but they do not yet justify a single harmonized patient identifier."
        ),
        provisional_rule_if_any=(
            "Use clinical_patient as the provisional patient universe, preserve both bcr_patient_barcode and bcr_patient_uuid as parallel evidence-bearing join keys, and attach biospecimen only at the sample anchor without emitting a canonical patient identifier."
        ),
        risk_if_wrong=(
            "Collapsing barcode and UUID evidence too early could over-link or under-link clinical and biospecimen records."
        ),
        needed_before_minimal_build=True,
    ),
    "AR02": ActionTemplate(
        action_title="Document case-count discrepancy and containment rule",
        action_scope="clinical coverage vs source coverage vs biospecimen sample anchor",
        rationale=(
            "The current saved layers expose a patient-count mismatch that should remain explicit and auditable before any cohort size is treated as stable."
        ),
        provisional_rule_if_any="",
        risk_if_wrong=(
            "A minimal build could inherit an unexamined denominator mismatch and be misread as a frozen patient universe."
        ),
        needed_before_minimal_build=True,
    ),
    "AR03": ActionTemplate(
        action_title="Define provisional endpoint source-table coexistence rule",
        action_scope="clinical_patient + clinical_follow_up_v4_0 endpoint-like fields",
        rationale=(
            "Overlapping endpoint-like fields already have saved provenance and completeness evidence, but they should stay as separate candidate sources rather than be silently collapsed."
        ),
        provisional_rule_if_any=(
            "Carry overlapping patient and follow-up endpoint-like fields as separate provenance-bearing candidates and do not set source precedence or final endpoint semantics here."
        ),
        risk_if_wrong=(
            "Prematurely collapsing patient and follow-up signals could hard-code an unsupported event or censoring rule."
        ),
        needed_before_minimal_build=True,
    ),
    "AR04": ActionTemplate(
        action_title="Exclude sparse death/progression timing fields from minimal build",
        action_scope="death_days_to and progression timing candidates",
        rationale=(
            "Saved timing evidence is too sparse for endpoint freeze, but it can still be documented as deferred evidence."
        ),
        provisional_rule_if_any=(
            "Exclude sparse death and progression timing fields from a minimal build and keep them only as deferred endpoint-review evidence."
        ),
        risk_if_wrong=(
            "Sparse timing fields could be mistaken for validated endpoint inputs when they remain incomplete or all-missing."
        ),
        needed_before_minimal_build=True,
    ),
    "AR05": ActionTemplate(
        action_title="Exclude one-to-many treatment-detail fields from minimal build",
        action_scope="clinical_drug + clinical_radiation treatment detail",
        rationale=(
            "Treatment rows remain useful source-audit evidence, but later aggregation or prioritization rules are still undecided."
        ),
        provisional_rule_if_any=(
            "Exclude one-to-many treatment-detail fields from a minimal build and keep treatment as optional proxy-only evidence until later aggregation design."
        ),
        risk_if_wrong=(
            "A build could imply unsupported treatment harmonization or patient-level exposure summarization."
        ),
        needed_before_minimal_build=True,
    ),
    "AR06": ActionTemplate(
        action_title="Defer analyte child-layer expansion pending later validation",
        action_scope="analyte -> slide and analyte -> aliquot expansion",
        rationale=(
            "The saved crosswalk is sufficient for sample-level anchoring, but it does not expose a direct analyte child link for deeper expansion."
        ),
        provisional_rule_if_any=(
            "Keep analyte-to-slide and analyte-to-aliquot expansion out of the minimal build and revisit only if later child-layer work becomes in-scope."
        ),
        risk_if_wrong=(
            "Later child-layer joins could be treated as validated even though current saved evidence does not support them directly."
        ),
        needed_before_minimal_build=False,
    ),
}


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def detect_repo_root(start_path: Path) -> Path:
    for candidate in [start_path, *start_path.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise BlueprintAmbiguityResolutionError(
        "Unable to locate the repository root from the script path."
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
                raise BlueprintAmbiguityResolutionError(
                    f"Nested YAML content has no parent key: {source_path}:{line_number}"
                )
            if stripped.startswith("- "):
                if current_key not in parsed or parsed[current_key] == {}:
                    parsed[current_key] = []
                if not isinstance(parsed[current_key], list):
                    raise BlueprintAmbiguityResolutionError(
                        f"Cannot mix list and scalar values for key '{current_key}' in {source_path}:{line_number}"
                    )
                parsed[current_key].append(parse_yaml_scalar(stripped[2:]))
                continue

            if ":" not in stripped:
                raise BlueprintAmbiguityResolutionError(
                    f"Expected nested key/value pair in YAML: {source_path}:{line_number}"
                )
            child_key, child_value = stripped.split(":", 1)
            if current_key not in parsed:
                parsed[current_key] = {}
            if not isinstance(parsed[current_key], dict):
                raise BlueprintAmbiguityResolutionError(
                    f"Cannot mix mapping and scalar values for key '{current_key}' in {source_path}:{line_number}"
                )
            parsed[current_key][child_key.strip()] = parse_yaml_scalar(child_value)
            continue

        if ":" not in raw_line:
            raise BlueprintAmbiguityResolutionError(
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
        raise BlueprintAmbiguityResolutionError(f"Expected a mapping in YAML config: {path}")
    return data


def load_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise BlueprintAmbiguityResolutionError(f"Expected a JSON object in file: {path}")
    return data


def read_tsv_dict_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def write_json(path: Path, payload: Any, overwrite: bool = False) -> None:
    if path.exists() and not overwrite:
        raise BlueprintAmbiguityResolutionError(f"Refusing to overwrite existing file: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_dict_rows_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    if path.exists():
        raise BlueprintAmbiguityResolutionError(f"Refusing to overwrite existing TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def create_run_directory(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise BlueprintAmbiguityResolutionError(f"Run directory already exists: {path}")
    path.mkdir(parents=False, exist_ok=False)
    return path


def parse_int(value: Any, label: str) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError) as exc:
        raise BlueprintAmbiguityResolutionError(
            f"Expected integer-like value for {label}: {value!r}"
        ) from exc


def require_keys(
    payload: dict[str, Any],
    required_keys: set[str],
    label: str,
    source_path: Path,
) -> None:
    missing_keys = required_keys.difference(payload.keys())
    if missing_keys:
        raise BlueprintAmbiguityResolutionError(
            f"{label} is missing required keys {sorted(missing_keys)}: {source_path}"
        )


def resolve_existing_path(repo_root: Path, relative_path: str, label: str) -> Path:
    path = repo_root / relative_path
    if not path.exists():
        raise BlueprintAmbiguityResolutionError(f"Required {label} not found: {path}")
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
    cohort_root = audit_root / "tcga-brca" / "cohort"

    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        blueprint_latest_pointer=cohort_root / "tcga_brca_cohort_blueprint_latest.json",
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
        ambiguity_resolution_runs_root=cohort_root / "ambiguity_resolution_runs",
        latest_pointer=cohort_root / "tcga_brca_blueprint_ambiguity_resolution_latest.json",
    )


def load_workflow_inputs(paths: WorkflowPaths) -> WorkflowInputs:
    for pointer_path, label in [
        (paths.blueprint_latest_pointer, "cohort blueprint latest pointer"),
        (paths.clinical_shortlist_latest_pointer, "clinical shortlist latest pointer"),
        (paths.endpoint_crosswalk_latest_pointer, "endpoint crosswalk latest pointer"),
        (paths.biospecimen_crosswalk_latest_pointer, "biospecimen crosswalk latest pointer"),
    ]:
        if not pointer_path.exists():
            raise BlueprintAmbiguityResolutionError(f"Required {label} not found: {pointer_path}")

    blueprint_latest_pointer = load_json(paths.blueprint_latest_pointer)
    clinical_shortlist_latest_pointer = load_json(paths.clinical_shortlist_latest_pointer)
    endpoint_crosswalk_latest_pointer = load_json(paths.endpoint_crosswalk_latest_pointer)
    biospecimen_crosswalk_latest_pointer = load_json(paths.biospecimen_crosswalk_latest_pointer)

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
            "clinical_shortlist_latest_json",
            "endpoint_crosswalk_latest_json",
            "biospecimen_crosswalk_latest_json",
        },
        "Cohort blueprint latest pointer",
        paths.blueprint_latest_pointer,
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

    if str(blueprint_latest_pointer["shortlist_run_id"]) != str(
        clinical_shortlist_latest_pointer["shortlist_run_id"]
    ):
        raise BlueprintAmbiguityResolutionError(
            "Cohort blueprint latest pointer does not reference the current clinical shortlist run."
        )
    if str(blueprint_latest_pointer["endpoint_crosswalk_run_id"]) != str(
        endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise BlueprintAmbiguityResolutionError(
            "Cohort blueprint latest pointer does not reference the current endpoint crosswalk run."
        )
    if str(blueprint_latest_pointer["biospecimen_crosswalk_run_id"]) != str(
        biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise BlueprintAmbiguityResolutionError(
            "Cohort blueprint latest pointer does not reference the current biospecimen crosswalk run."
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
        "biospecimen_inventory_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_inventory_tsv"]),
            "biospecimen identifier inventory TSV",
        ),
        "biospecimen_crosswalk_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_crosswalk_tsv"]),
            "biospecimen identifier crosswalk TSV",
        ),
        "biospecimen_crosswalk_summary_tsv": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["biospecimen_identifier_crosswalk_summary_tsv"]),
            "biospecimen identifier crosswalk summary TSV",
        ),
        "biospecimen_crosswalk_run_log_json": resolve_existing_path(
            paths.repo_root,
            str(biospecimen_crosswalk_latest_pointer["run_log_json"]),
            "biospecimen identifier crosswalk run log",
        ),
    }

    blueprint_run_log = load_json(input_paths["blueprint_run_log_json"])
    clinical_shortlist_run_log = load_json(input_paths["clinical_shortlist_run_log_json"])
    endpoint_crosswalk_run_log = load_json(input_paths["endpoint_crosswalk_run_log_json"])
    biospecimen_crosswalk_run_log = load_json(input_paths["biospecimen_crosswalk_run_log_json"])
    for run_log, label in [
        (blueprint_run_log, "cohort blueprint"),
        (clinical_shortlist_run_log, "clinical shortlist"),
        (endpoint_crosswalk_run_log, "endpoint crosswalk"),
        (biospecimen_crosswalk_run_log, "biospecimen crosswalk"),
    ]:
        if run_log.get("status") != "completed":
            raise BlueprintAmbiguityResolutionError(f"Upstream {label} run log is not completed.")

    if str(blueprint_run_log.get("blueprint_run_id") or "") != str(
        blueprint_latest_pointer["blueprint_run_id"]
    ):
        raise BlueprintAmbiguityResolutionError(
            "Cohort blueprint run log ID does not match the cohort blueprint latest pointer."
        )
    if str(clinical_shortlist_run_log.get("shortlist_run_id") or "") != str(
        clinical_shortlist_latest_pointer["shortlist_run_id"]
    ):
        raise BlueprintAmbiguityResolutionError(
            "Clinical shortlist run log ID does not match the clinical shortlist latest pointer."
        )
    if str(endpoint_crosswalk_run_log.get("crosswalk_run_id") or "") != str(
        endpoint_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise BlueprintAmbiguityResolutionError(
            "Endpoint crosswalk run log ID does not match the endpoint crosswalk latest pointer."
        )
    if str(biospecimen_crosswalk_run_log.get("crosswalk_run_id") or "") != str(
        biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
    ):
        raise BlueprintAmbiguityResolutionError(
            "Biospecimen crosswalk run log ID does not match the biospecimen crosswalk latest pointer."
        )

    blueprint_required_rows = read_tsv_dict_rows(input_paths["blueprint_required_tsv"])
    blueprint_optional_rows = read_tsv_dict_rows(input_paths["blueprint_optional_tsv"])
    blueprint_deferred_rows = read_tsv_dict_rows(input_paths["blueprint_deferred_tsv"])
    blueprint_join_path_rows = read_tsv_dict_rows(input_paths["blueprint_join_path_tsv"])
    blueprint_ambiguity_rows = read_tsv_dict_rows(input_paths["blueprint_ambiguities_tsv"])
    blueprint_summary_rows = read_tsv_dict_rows(input_paths["blueprint_summary_tsv"])
    clinical_shortlist_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_tsv"])
    clinical_shortlist_summary_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_summary_tsv"])
    clinical_shortlist_by_table_rows = read_tsv_dict_rows(input_paths["clinical_shortlist_by_table_tsv"])
    endpoint_inventory_rows = read_tsv_dict_rows(input_paths["endpoint_candidate_inventory_tsv"])
    endpoint_crosswalk_rows = read_tsv_dict_rows(input_paths["endpoint_crosswalk_tsv"])
    endpoint_summary_rows = read_tsv_dict_rows(input_paths["endpoint_crosswalk_summary_tsv"])
    biospecimen_inventory_rows = read_tsv_dict_rows(input_paths["biospecimen_inventory_tsv"])
    biospecimen_crosswalk_rows = read_tsv_dict_rows(input_paths["biospecimen_crosswalk_tsv"])
    biospecimen_summary_rows = read_tsv_dict_rows(input_paths["biospecimen_crosswalk_summary_tsv"])

    for rows, label in [
        (blueprint_required_rows, "cohort blueprint required fields"),
        (blueprint_optional_rows, "cohort blueprint optional fields"),
        (blueprint_deferred_rows, "cohort blueprint deferred fields"),
        (blueprint_join_path_rows, "cohort blueprint join path"),
        (blueprint_ambiguity_rows, "cohort blueprint ambiguities"),
        (blueprint_summary_rows, "cohort blueprint summary"),
        (clinical_shortlist_rows, "clinical shortlist"),
        (clinical_shortlist_summary_rows, "clinical shortlist summary"),
        (clinical_shortlist_by_table_rows, "clinical shortlist by-table summary"),
        (endpoint_inventory_rows, "endpoint candidate inventory"),
        (endpoint_crosswalk_rows, "endpoint crosswalk"),
        (endpoint_summary_rows, "endpoint crosswalk summary"),
        (biospecimen_inventory_rows, "biospecimen inventory"),
        (biospecimen_crosswalk_rows, "biospecimen crosswalk"),
        (biospecimen_summary_rows, "biospecimen crosswalk summary"),
    ]:
        if not rows:
            raise BlueprintAmbiguityResolutionError(
                f"Expected saved {label} rows but found zero rows."
            )

    observed_ambiguity_ids = {str(row["ambiguity_id"]) for row in blueprint_ambiguity_rows}
    expected_ambiguity_ids = set(CLASSIFICATION_RULES)
    if observed_ambiguity_ids != expected_ambiguity_ids:
        raise BlueprintAmbiguityResolutionError(
            "Blueprint ambiguity IDs do not match the fixed ambiguity-resolution rules: "
            f"observed={sorted(observed_ambiguity_ids)} expected={sorted(expected_ambiguity_ids)}"
        )

    return WorkflowInputs(
        blueprint_latest_pointer=blueprint_latest_pointer,
        blueprint_run_log=blueprint_run_log,
        clinical_shortlist_latest_pointer=clinical_shortlist_latest_pointer,
        clinical_shortlist_run_log=clinical_shortlist_run_log,
        endpoint_crosswalk_latest_pointer=endpoint_crosswalk_latest_pointer,
        endpoint_crosswalk_run_log=endpoint_crosswalk_run_log,
        biospecimen_crosswalk_latest_pointer=biospecimen_crosswalk_latest_pointer,
        biospecimen_crosswalk_run_log=biospecimen_crosswalk_run_log,
        blueprint_required_rows=blueprint_required_rows,
        blueprint_optional_rows=blueprint_optional_rows,
        blueprint_deferred_rows=blueprint_deferred_rows,
        blueprint_join_path_rows=blueprint_join_path_rows,
        blueprint_ambiguity_rows=blueprint_ambiguity_rows,
        blueprint_summary_rows=blueprint_summary_rows,
        clinical_shortlist_rows=clinical_shortlist_rows,
        clinical_shortlist_summary_rows=clinical_shortlist_summary_rows,
        clinical_shortlist_by_table_rows=clinical_shortlist_by_table_rows,
        endpoint_inventory_rows=endpoint_inventory_rows,
        endpoint_crosswalk_rows=endpoint_crosswalk_rows,
        endpoint_summary_rows=endpoint_summary_rows,
        biospecimen_inventory_rows=biospecimen_inventory_rows,
        biospecimen_crosswalk_rows=biospecimen_crosswalk_rows,
        biospecimen_summary_rows=biospecimen_summary_rows,
        input_paths=input_paths,
    )


def find_single_row(
    rows: list[dict[str, str]],
    label: str,
    **conditions: str,
) -> dict[str, str]:
    matches = [
        row
        for row in rows
        if all(str(row.get(key, "")) == str(expected) for key, expected in conditions.items())
    ]
    if len(matches) != 1:
        raise BlueprintAmbiguityResolutionError(
            f"Expected exactly one {label} row matching {conditions}, found {len(matches)}."
        )
    return matches[0]


def build_summary_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        metric = str(row.get("summary_metric") or "")
        if metric in lookup:
            raise BlueprintAmbiguityResolutionError(
                f"Duplicate summary_metric detected in summary TSV: {metric}"
            )
        lookup[metric] = row
    return lookup


def build_shortlist_bucket_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        bucket = str(row.get("shortlist_bucket") or "")
        if bucket:
            lookup[bucket] = row
    return lookup


def json_list(values: list[str]) -> str:
    return json.dumps(values, ensure_ascii=True)


def summary_int(summary_lookup: dict[str, dict[str, str]], metric: str) -> int:
    if metric not in summary_lookup:
        raise BlueprintAmbiguityResolutionError(f"Required summary metric not found: {metric}")
    return parse_int(summary_lookup[metric]["summary_value"], metric)


def build_inventory_current_evidence(
    blueprint_row: dict[str, str],
    workflow_inputs: WorkflowInputs,
    blueprint_summary_lookup: dict[str, dict[str, str]],
    shortlist_bucket_lookup: dict[str, dict[str, str]],
) -> str:
    ambiguity_id = str(blueprint_row["ambiguity_id"])
    blueprint_evidence = str(blueprint_row["current_saved_evidence"])

    if ambiguity_id == "A01":
        patient_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="vital_status",
        )
        followup_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="vital_status",
        )
        return (
            f"{blueprint_evidence} Endpoint crosswalk roles: clinical_patient.vital_status="
            f"{patient_row['crosswalk_role']} with {patient_row['non_missing_count']} non-missing rows; "
            f"clinical_follow_up_v4_0.vital_status={followup_row['crosswalk_role']} with "
            f"{followup_row['non_missing_count']} non-missing rows."
        )

    if ambiguity_id == "A02":
        patient_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="last_contact_days_to",
        )
        followup_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="last_contact_days_to",
        )
        return (
            f"{blueprint_evidence} Endpoint crosswalk roles: clinical_follow_up_v4_0.last_contact_days_to="
            f"{followup_row['crosswalk_role']} with {followup_row['non_missing_count']} non-missing rows; "
            f"clinical_patient.last_contact_days_to={patient_row['crosswalk_role']} with "
            f"{patient_row['non_missing_count']} non-missing rows."
        )

    if ambiguity_id == "A03":
        patient_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="tumor_status",
        )
        followup_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="tumor_status",
        )
        return (
            f"{blueprint_evidence} Endpoint crosswalk roles: clinical_follow_up_v4_0.tumor_status="
            f"{followup_row['crosswalk_role']} with {followup_row['non_missing_count']} non-missing rows; "
            f"clinical_patient.tumor_status={patient_row['crosswalk_role']} with "
            f"{patient_row['non_missing_count']} non-missing rows."
        )

    if ambiguity_id == "A04":
        patient_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="new_tumor_event_dx_indicator",
        )
        followup_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="new_tumor_event_dx_indicator",
        )
        return (
            f"{blueprint_evidence} Endpoint crosswalk roles: clinical_follow_up_v4_0.new_tumor_event_dx_indicator="
            f"{followup_row['crosswalk_role']} with {followup_row['non_missing_count']} non-missing rows; "
            f"clinical_patient.new_tumor_event_dx_indicator={patient_row['crosswalk_role']} with "
            f"{patient_row['non_missing_count']} non-missing rows."
        )

    if ambiguity_id == "A05":
        patient_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="death_days_to",
        )
        followup_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="death_days_to",
        )
        return (
            f"{blueprint_evidence} Endpoint crosswalk roles: clinical_patient.death_days_to="
            f"{patient_row['crosswalk_role']} with {patient_row['non_missing_count']} non-missing rows; "
            f"clinical_follow_up_v4_0.death_days_to={followup_row['crosswalk_role']} with "
            f"{followup_row['non_missing_count']} non-missing rows."
        )

    if ambiguity_id == "A06":
        progression_free_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="days_to_patient_progression_free",
        )
        tumor_progression_row = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="days_to_tumor_progression",
        )
        return (
            f"{blueprint_evidence} Endpoint crosswalk marks both progression timing candidates as "
            f"{progression_free_row['crosswalk_role']} / {tumor_progression_row['crosswalk_role']} with "
            f"{progression_free_row['non_missing_count']} and {tumor_progression_row['non_missing_count']} non-missing rows."
        )

    if ambiguity_id == "A07":
        clinical_barcode_row = find_single_row(
            workflow_inputs.clinical_shortlist_rows,
            "clinical shortlist",
            table_name="clinical_patient",
            field_name="bcr_patient_barcode",
        )
        clinical_uuid_row = find_single_row(
            workflow_inputs.clinical_shortlist_rows,
            "clinical shortlist",
            table_name="clinical_patient",
            field_name="bcr_patient_uuid",
        )
        biospecimen_sample_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_sample",
            field_name="bcr_patient_uuid",
        )
        diagnostic_slide_barcode_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_diagnostic_slides",
            field_name="bcr_patient_barcode",
        )
        diagnostic_slide_uuid_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_diagnostic_slides",
            field_name="bcr_patient_uuid",
        )
        return (
            f"{blueprint_evidence} Clinical shortlist keeps clinical_patient.bcr_patient_barcode and "
            f"clinical_patient.bcr_patient_uuid with {clinical_barcode_row['non_missing_count']} and "
            f"{clinical_uuid_row['non_missing_count']} non-missing rows. Biospecimen crosswalk marks "
            f"biospecimen_sample.bcr_patient_uuid as {biospecimen_sample_row['crosswalk_role']} with "
            f"{biospecimen_sample_row['distinct_non_missing_count']} distinct patients, and "
            f"biospecimen_diagnostic_slides.bcr_patient_barcode / bcr_patient_uuid remain exact paired side evidence "
            f"at pattern pass fractions {diagnostic_slide_barcode_row['pattern_pass_fraction']} and "
            f"{diagnostic_slide_uuid_row['pattern_pass_fraction']}."
        )

    if ambiguity_id == "A08":
        return (
            f"{blueprint_evidence} Blueprint summary count signals: clinical_patient_case_count="
            f"{summary_int(blueprint_summary_lookup, 'clinical_patient_case_count')}, "
            f"clinical_source_case_submitter_coverage={summary_int(blueprint_summary_lookup, 'clinical_source_case_submitter_coverage')}, "
            f"biospecimen_patient_uuid_distinct_count={summary_int(blueprint_summary_lookup, 'biospecimen_patient_uuid_distinct_count')}, "
            f"biospecimen_sample_anchor_count={summary_int(blueprint_summary_lookup, 'biospecimen_sample_anchor_count')}."
        )

    if ambiguity_id == "A09":
        treatment_bucket = shortlist_bucket_lookup["usable_treatment_proxy"]
        return (
            f"{blueprint_evidence} Clinical shortlist currently carries "
            f"{treatment_bucket['field_count']} usable_treatment_proxy fields for treatment-proxy review only."
        )

    if ambiguity_id == "A10":
        analyte_barcode_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_analyte",
            field_name="bcr_analyte_barcode",
        )
        slide_sample_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_slide",
            field_name="bcr_sample_barcode",
        )
        slide_barcode_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_slide",
            field_name="bcr_slide_barcode",
        )
        return (
            f"{blueprint_evidence} Biospecimen crosswalk keeps analyte identifiers on biospecimen_analyte "
            f"({analyte_barcode_row['distinct_non_missing_count']} distinct analyte barcodes), while biospecimen_slide "
            f"carries slide ids and overlapping sample barcode evidence ({slide_barcode_row['distinct_non_missing_count']} "
            f"slide barcodes; {slide_sample_row['distinct_non_missing_count']} overlapping sample barcodes)."
        )

    if ambiguity_id == "A11":
        analyte_barcode_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_analyte",
            field_name="bcr_analyte_barcode",
        )
        aliquot_sample_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_aliquot",
            field_name="bcr_sample_barcode",
        )
        aliquot_barcode_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_aliquot",
            field_name="bcr_aliquot_barcode",
        )
        return (
            f"{blueprint_evidence} Biospecimen crosswalk keeps analyte identifiers on biospecimen_analyte "
            f"({analyte_barcode_row['distinct_non_missing_count']} distinct analyte barcodes), while "
            f"biospecimen_aliquot carries aliquot ids and overlapping sample barcode evidence "
            f"({aliquot_barcode_row['distinct_non_missing_count']} aliquot barcodes; "
            f"{aliquot_sample_row['distinct_non_missing_count']} overlapping sample barcodes)."
        )

    raise BlueprintAmbiguityResolutionError(f"Unhandled ambiguity_id for evidence build: {ambiguity_id}")


def build_inventory_rows(
    ambiguity_resolution_run_id: str,
    workflow_inputs: WorkflowInputs,
) -> list[dict[str, Any]]:
    blueprint_summary_lookup = build_summary_lookup(workflow_inputs.blueprint_summary_rows)
    shortlist_bucket_lookup = build_shortlist_bucket_lookup(workflow_inputs.clinical_shortlist_summary_rows)
    inventory_rows: list[dict[str, Any]] = []

    for blueprint_row in sorted(
        workflow_inputs.blueprint_ambiguity_rows,
        key=lambda row: str(row["ambiguity_id"]),
    ):
        ambiguity_id = str(blueprint_row["ambiguity_id"])
        rule = CLASSIFICATION_RULES[ambiguity_id]
        action_template = ACTION_TEMPLATES[rule.linked_action_id]
        inventory_rows.append(
            {
                "ambiguity_resolution_run_id": ambiguity_resolution_run_id,
                "blueprint_run_id": str(blueprint_row["blueprint_run_id"]),
                "ambiguity_id": ambiguity_id,
                "ambiguity_category": str(blueprint_row["ambiguity_category"]),
                "source_layer": str(blueprint_row["source_layer"]),
                "table_name": str(blueprint_row["table_name"]),
                "field_name_or_step": str(blueprint_row["field_name_or_step"]),
                "severity": str(blueprint_row["severity"]),
                "why_unresolved": str(blueprint_row["why_unresolved"]),
                "current_saved_evidence": build_inventory_current_evidence(
                    blueprint_row=blueprint_row,
                    workflow_inputs=workflow_inputs,
                    blueprint_summary_lookup=blueprint_summary_lookup,
                    shortlist_bucket_lookup=shortlist_bucket_lookup,
                ),
                "proposed_blueprint_handling": str(blueprint_row["proposed_blueprint_handling"]),
                "blocks_final_cohort_construction": str(
                    blueprint_row["blocks_final_cohort_construction"]
                ),
                "evidence_still_missing": rule.evidence_still_missing,
                "blocker_class": rule.blocker_class,
                "resolution_priority": rule.resolution_priority,
                "evidence_sufficiency": rule.evidence_sufficiency,
                "recommended_next_action": action_template.action_title,
                "recommended_action_type": rule.recommended_action_type,
                "linked_action_id": rule.linked_action_id,
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )

    return inventory_rows


def build_action_evidence_basis(
    action_id: str,
    inventory_rows: list[dict[str, Any]],
    workflow_inputs: WorkflowInputs,
) -> str:
    summary_lookup = build_summary_lookup(workflow_inputs.blueprint_summary_rows)

    if action_id == "AR01":
        sample_uuid_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_sample",
            field_name="bcr_patient_uuid",
        )
        slide_barcode_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_diagnostic_slides",
            field_name="bcr_patient_barcode",
        )
        return (
            "clinical_patient preserves both patient identifiers with 1097 non-missing rows each; "
            f"biospecimen_sample.bcr_patient_uuid carries {sample_uuid_row['distinct_non_missing_count']} distinct patients; "
            f"biospecimen_diagnostic_slides.bcr_patient_barcode pairs exactly with bcr_patient_uuid at "
            f"pattern pass fraction {slide_barcode_row['pattern_pass_fraction']}."
        )

    if action_id == "AR02":
        return (
            "Blueprint summary count signals currently disagree: "
            f"clinical_patient_case_count={summary_int(summary_lookup, 'clinical_patient_case_count')}, "
            f"clinical_source_case_submitter_coverage={summary_int(summary_lookup, 'clinical_source_case_submitter_coverage')}, "
            f"biospecimen_patient_uuid_distinct_count={summary_int(summary_lookup, 'biospecimen_patient_uuid_distinct_count')}."
        )

    if action_id == "AR03":
        vital_patient = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="vital_status",
        )
        vital_followup = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="vital_status",
        )
        last_contact_patient = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="last_contact_days_to",
        )
        last_contact_followup = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="last_contact_days_to",
        )
        return (
            "Endpoint crosswalk already distinguishes overlapping candidates instead of collapsing them: "
            f"vital_status patient/follow-up={vital_patient['crosswalk_role']}/{vital_followup['crosswalk_role']}; "
            f"last_contact_days_to patient/follow-up={last_contact_patient['crosswalk_role']}/{last_contact_followup['crosswalk_role']}."
        )

    if action_id == "AR04":
        death_patient = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="death_days_to",
        )
        death_followup = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_follow_up_v4_0",
            field_name="death_days_to",
        )
        progression_free = find_single_row(
            workflow_inputs.endpoint_crosswalk_rows,
            "endpoint crosswalk",
            table_name="clinical_patient",
            field_name="days_to_patient_progression_free",
        )
        return (
            "Saved timing evidence remains sparse: "
            f"clinical_patient.death_days_to={death_patient['non_missing_count']} non-missing, "
            f"clinical_follow_up_v4_0.death_days_to={death_followup['non_missing_count']} non-missing, "
            f"days_to_patient_progression_free={progression_free['non_missing_count']} non-missing."
        )

    if action_id == "AR05":
        treatment_rows = [
            row for row in inventory_rows if str(row["linked_action_id"]) == "AR05"
        ]
        if len(treatment_rows) != 1:
            raise BlueprintAmbiguityResolutionError("Expected exactly one treatment ambiguity row.")
        shortlist_bucket_lookup = build_shortlist_bucket_lookup(workflow_inputs.clinical_shortlist_summary_rows)
        treatment_bucket = shortlist_bucket_lookup["usable_treatment_proxy"]
        return (
            f"{treatment_rows[0]['current_saved_evidence']} Clinical shortlist keeps "
            f"{treatment_bucket['field_count']} treatment-proxy fields in a review-only bucket."
        )

    if action_id == "AR06":
        analyte_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_analyte",
            field_name="bcr_analyte_barcode",
        )
        slide_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_slide",
            field_name="bcr_sample_barcode",
        )
        aliquot_row = find_single_row(
            workflow_inputs.biospecimen_crosswalk_rows,
            "biospecimen crosswalk",
            table_name="biospecimen_aliquot",
            field_name="bcr_sample_barcode",
        )
        return (
            "Saved child-layer evidence stays indirect: "
            f"analyte identifiers remain on biospecimen_analyte ({analyte_row['distinct_non_missing_count']} distinct), "
            f"while slide and aliquot layers retain overlapping sample barcodes ({slide_row['distinct_non_missing_count']} and "
            f"{aliquot_row['distinct_non_missing_count']} distinct sample barcodes)."
        )

    raise BlueprintAmbiguityResolutionError(f"Unhandled action_id for evidence basis: {action_id}")


def build_action_rows(
    ambiguity_resolution_run_id: str,
    inventory_rows: list[dict[str, Any]],
    workflow_inputs: WorkflowInputs,
) -> list[dict[str, Any]]:
    action_rows: list[dict[str, Any]] = []
    for action_id in ACTION_ORDER:
        linked_rows = [
            row for row in inventory_rows if str(row["linked_action_id"]) == action_id
        ]
        if not linked_rows:
            raise BlueprintAmbiguityResolutionError(
                f"Action {action_id} has no linked ambiguity rows."
            )
        template = ACTION_TEMPLATES[action_id]
        action_rows.append(
            {
                "ambiguity_resolution_run_id": ambiguity_resolution_run_id,
                "action_id": action_id,
                "action_title": template.action_title,
                "action_scope": template.action_scope,
                "linked_ambiguity_ids_json": json_list(
                    sorted(str(row["ambiguity_id"]) for row in linked_rows)
                ),
                "rationale": template.rationale,
                "provisional_rule_if_any": template.provisional_rule_if_any,
                "evidence_basis": build_action_evidence_basis(
                    action_id=action_id,
                    inventory_rows=inventory_rows,
                    workflow_inputs=workflow_inputs,
                ),
                "risk_if_wrong": template.risk_if_wrong,
                "needed_before_minimal_build": template.needed_before_minimal_build,
                "notes_placeholder": NOTES_PLACEHOLDER,
            }
        )
    return action_rows


def build_priority_rows(
    ambiguity_resolution_run_id: str,
    inventory_rows: list[dict[str, Any]],
    action_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    action_row_by_id = {str(row["action_id"]): row for row in action_rows}
    unsorted_rows: list[dict[str, Any]] = []

    for action_id in ACTION_ORDER:
        linked_rows = [
            row for row in inventory_rows if str(row["linked_action_id"]) == action_id
        ]
        if not linked_rows:
            raise BlueprintAmbiguityResolutionError(
                f"Action {action_id} has no linked inventory rows for priority build."
            )
        blocker_class = min(
            (str(row["blocker_class"]) for row in linked_rows),
            key=BLOCKER_CLASS_ORDER.index,
        )
        resolution_priority = min(
            (str(row["resolution_priority"]) for row in linked_rows),
            key=RESOLUTION_PRIORITY_ORDER.index,
        )
        recommended_action_type = min(
            (str(row["recommended_action_type"]) for row in linked_rows),
            key=RECOMMENDED_ACTION_TYPE_ORDER.index,
        )
        categories = sorted({str(row["ambiguity_category"]) for row in linked_rows})
        action_row = action_row_by_id[action_id]
        unsorted_rows.append(
            {
                "ambiguity_resolution_run_id": ambiguity_resolution_run_id,
                "priority_rank": 0,
                "action_id": action_id,
                "action_title": str(action_row["action_title"]),
                "blocker_class": blocker_class,
                "resolution_priority": resolution_priority,
                "recommended_action_type": recommended_action_type,
                "linked_ambiguity_count": len(linked_rows),
                "linked_high_severity_count": sum(
                    str(row["severity"]) == "high" for row in linked_rows
                ),
                "linked_ambiguity_categories_json": json_list(categories),
                "needed_before_minimal_build": bool(action_row["needed_before_minimal_build"]),
            }
        )

    sorted_rows = sorted(
        unsorted_rows,
        key=lambda row: (
            0 if bool(row["needed_before_minimal_build"]) else 1,
            BLOCKER_CLASS_ORDER.index(str(row["blocker_class"])),
            RESOLUTION_PRIORITY_ORDER.index(str(row["resolution_priority"])),
            -parse_int(row["linked_ambiguity_count"], "linked_ambiguity_count"),
            str(row["action_id"]),
        ),
    )

    for index, row in enumerate(sorted_rows, start=1):
        row["priority_rank"] = index

    return sorted_rows


def determine_minimal_build_readiness(
    action_rows: list[dict[str, Any]],
    priority_rows: list[dict[str, Any]],
) -> tuple[str, str]:
    priority_by_action = {str(row["action_id"]): row for row in priority_rows}
    needed_rows = [row for row in action_rows if bool(row["needed_before_minimal_build"])]
    if not needed_rows:
        return (
            "ready_no_prebuild_actions_flagged",
            "No action is currently flagged as needed before a minimal build.",
        )

    if any(
        str(priority_by_action[str(row["action_id"])]["recommended_action_type"]) == "manual_review"
        for row in needed_rows
    ):
        return (
            "not_ready_pending_manual_review",
            "At least one action needed before a minimal build still requires manual review.",
        )

    return (
        "ready_after_documented_rules_applied",
        "Only provisional-rule or explicit deferral actions remain before a minimal build.",
    )


def build_summary_rows(
    ambiguity_resolution_run_id: str,
    workflow_inputs: WorkflowInputs,
    inventory_rows: list[dict[str, Any]],
    action_rows: list[dict[str, Any]],
    priority_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    readiness_signal, readiness_notes = determine_minimal_build_readiness(action_rows, priority_rows)

    def count_inventory(field_name: str, expected_value: str) -> int:
        return sum(str(row[field_name]) == expected_value for row in inventory_rows)

    def add_row(
        rows: list[dict[str, Any]],
        summary_section: str,
        summary_metric: str,
        summary_value: Any,
        notes: str,
    ) -> None:
        rows.append(
            {
                "ambiguity_resolution_run_id": ambiguity_resolution_run_id,
                "summary_section": summary_section,
                "summary_metric": summary_metric,
                "summary_value": summary_value,
                "notes": notes,
            }
        )

    summary_rows: list[dict[str, Any]] = []
    add_row(
        summary_rows,
        "inputs",
        "blueprint_run_id",
        workflow_inputs.blueprint_latest_pointer["blueprint_run_id"],
        "Blueprint run used as the direct ambiguity source.",
    )
    add_row(
        summary_rows,
        "inputs",
        "blueprint_ambiguity_count",
        len(workflow_inputs.blueprint_ambiguity_rows),
        "Count of ambiguity rows loaded from the current cohort blueprint.",
    )
    add_row(
        summary_rows,
        "inventory",
        "ambiguity_resolution_inventory_count",
        len(inventory_rows),
        "Count of ambiguity rows written into the ambiguity-resolution inventory.",
    )
    for blocker_class in BLOCKER_CLASS_ORDER:
        add_row(
            summary_rows,
            "inventory",
            f"{blocker_class}_count",
            count_inventory("blocker_class", blocker_class),
            "Inventory count by minimal-build blocker class.",
        )
    for evidence_sufficiency in EVIDENCE_SUFFICIENCY_ORDER:
        add_row(
            summary_rows,
            "inventory",
            f"{evidence_sufficiency}_count",
            count_inventory("evidence_sufficiency", evidence_sufficiency),
            "Inventory count by evidence sufficiency status.",
        )
    for action_type in RECOMMENDED_ACTION_TYPE_ORDER:
        add_row(
            summary_rows,
            "actions",
            f"{action_type}_count",
            count_inventory("recommended_action_type", action_type),
            "Inventory count by recommended action type.",
        )
    add_row(
        summary_rows,
        "actions",
        "action_count",
        len(action_rows),
        "Number of grouped ambiguity-resolution actions.",
    )
    add_row(
        summary_rows,
        "actions",
        "actions_needed_before_minimal_build_count",
        sum(bool(row["needed_before_minimal_build"]) for row in action_rows),
        "Number of grouped actions currently flagged as needed before a minimal build.",
    )
    add_row(
        summary_rows,
        "readiness",
        "minimal_build_readiness_signal",
        readiness_signal,
        readiness_notes,
    )
    add_row(
        summary_rows,
        "readiness",
        "top_priority_action_id",
        priority_rows[0]["action_id"],
        "Highest-priority grouped action after applying the deterministic ranking rules.",
    )
    add_row(
        summary_rows,
        "upstream_runs",
        "clinical_shortlist_run_id",
        workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"],
        "Current clinical shortlist run used by this ambiguity-resolution layer.",
    )
    add_row(
        summary_rows,
        "upstream_runs",
        "endpoint_crosswalk_run_id",
        workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"],
        "Current endpoint crosswalk run used by this ambiguity-resolution layer.",
    )
    add_row(
        summary_rows,
        "upstream_runs",
        "biospecimen_crosswalk_run_id",
        workflow_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"],
        "Current biospecimen crosswalk run used by this ambiguity-resolution layer.",
    )
    return summary_rows


def build_latest_pointer_payload(
    ambiguity_resolution_run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "ambiguity_resolution_run_id": ambiguity_resolution_run_id,
        "blueprint_run_id": str(workflow_inputs.blueprint_latest_pointer["blueprint_run_id"]),
        "shortlist_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
        "core_audit_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["core_audit_run_id"]),
        "clinical_parse_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["parse_run_id"]),
        "endpoint_crosswalk_run_id": str(workflow_inputs.endpoint_crosswalk_latest_pointer["crosswalk_run_id"]),
        "biospecimen_crosswalk_run_id": str(
            workflow_inputs.biospecimen_crosswalk_latest_pointer["crosswalk_run_id"]
        ),
        "biospecimen_parse_run_id": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["parse_run_id"]),
        "clinical_source_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["source_run_id"]),
        "biospecimen_source_run_id": str(workflow_inputs.biospecimen_crosswalk_latest_pointer["source_run_id"]),
        "ambiguity_resolution_run_directory": repo_relative(
            output_paths["ambiguity_resolution_run_directory"],
            paths.repo_root,
        ),
        "ambiguity_resolution_inventory_tsv": repo_relative(output_paths["inventory_tsv"], paths.repo_root),
        "ambiguity_resolution_priority_tsv": repo_relative(output_paths["priority_tsv"], paths.repo_root),
        "ambiguity_resolution_actions_tsv": repo_relative(output_paths["actions_tsv"], paths.repo_root),
        "ambiguity_resolution_summary_tsv": repo_relative(output_paths["summary_tsv"], paths.repo_root),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "cohort_blueprint_latest_json": repo_relative(paths.blueprint_latest_pointer, paths.repo_root),
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
    }


def write_failure_log(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    ambiguity_resolution_run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    paths = build_workflow_paths()
    ambiguity_resolution_run_dir = paths.ambiguity_resolution_runs_root / ambiguity_resolution_run_id
    run_log_path = ambiguity_resolution_run_dir / "run_log.json"

    try:
        trial_config = load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths)
        create_run_directory(ambiguity_resolution_run_dir)

        inventory_path = ambiguity_resolution_run_dir / "ambiguity_resolution_inventory.tsv"
        priority_path = ambiguity_resolution_run_dir / "ambiguity_resolution_priority.tsv"
        actions_path = ambiguity_resolution_run_dir / "ambiguity_resolution_actions.tsv"
        summary_path = ambiguity_resolution_run_dir / "ambiguity_resolution_summary.tsv"

        inventory_rows = build_inventory_rows(
            ambiguity_resolution_run_id=ambiguity_resolution_run_id,
            workflow_inputs=workflow_inputs,
        )
        if not inventory_rows:
            raise BlueprintAmbiguityResolutionError(
                "Ambiguity-resolution inventory rows were not generated."
            )
        action_rows = build_action_rows(
            ambiguity_resolution_run_id=ambiguity_resolution_run_id,
            inventory_rows=inventory_rows,
            workflow_inputs=workflow_inputs,
        )
        priority_rows = build_priority_rows(
            ambiguity_resolution_run_id=ambiguity_resolution_run_id,
            inventory_rows=inventory_rows,
            action_rows=action_rows,
        )
        summary_rows = build_summary_rows(
            ambiguity_resolution_run_id=ambiguity_resolution_run_id,
            workflow_inputs=workflow_inputs,
            inventory_rows=inventory_rows,
            action_rows=action_rows,
            priority_rows=priority_rows,
        )

        write_dict_rows_tsv(inventory_path, INVENTORY_FIELDNAMES, inventory_rows)
        write_dict_rows_tsv(priority_path, PRIORITY_FIELDNAMES, priority_rows)
        write_dict_rows_tsv(actions_path, ACTIONS_FIELDNAMES, action_rows)
        write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        inventory_row_count_matches_blueprint = len(inventory_rows) == len(
            workflow_inputs.blueprint_ambiguity_rows
        )
        blocker_class_values_valid = all(
            str(row["blocker_class"]) in BLOCKER_CLASS_ORDER for row in inventory_rows
        )
        resolution_priority_values_valid = all(
            str(row["resolution_priority"]) in RESOLUTION_PRIORITY_ORDER for row in inventory_rows
        )
        evidence_sufficiency_values_valid = all(
            str(row["evidence_sufficiency"]) in EVIDENCE_SUFFICIENCY_ORDER for row in inventory_rows
        )
        recommended_action_type_values_valid = all(
            str(row["recommended_action_type"]) in RECOMMENDED_ACTION_TYPE_ORDER for row in inventory_rows
        )
        every_ambiguity_has_one_linked_action = all(
            str(row["linked_action_id"]) in ACTION_TEMPLATES for row in inventory_rows
        )
        flattened_link_ids: list[str] = []
        for action_row in action_rows:
            flattened_link_ids.extend(json.loads(str(action_row["linked_ambiguity_ids_json"])))
        flattened_link_counts = {
            ambiguity_id: flattened_link_ids.count(ambiguity_id)
            for ambiguity_id in sorted({*flattened_link_ids})
        }
        inventory_ambiguity_ids = sorted(str(row["ambiguity_id"]) for row in inventory_rows)
        action_linkage_covers_each_ambiguity_once = (
            sorted(flattened_link_ids) == inventory_ambiguity_ids
            and all(count == 1 for count in flattened_link_counts.values())
        )
        priority_rows_match_action_rows = len(priority_rows) == len(action_rows)
        summary_lookup = build_summary_lookup(summary_rows)
        summary_counts_reconcile_to_inventory = (
            summary_int(summary_lookup, "ambiguity_resolution_inventory_count") == len(inventory_rows)
            and summary_int(summary_lookup, "blocking_now_count")
            == sum(str(row["blocker_class"]) == "blocking_now" for row in inventory_rows)
            and summary_int(summary_lookup, "review_before_build_count")
            == sum(str(row["blocker_class"]) == "review_before_build" for row in inventory_rows)
            and summary_int(summary_lookup, "non_blocking_for_minimal_build_count")
            == sum(
                str(row["blocker_class"]) == "non_blocking_for_minimal_build"
                for row in inventory_rows
            )
            and summary_int(summary_lookup, "action_count") == len(action_rows)
            and summary_int(summary_lookup, "actions_needed_before_minimal_build_count")
            == sum(bool(row["needed_before_minimal_build"]) for row in action_rows)
        )
        output_rows_positive = all(
            len(rows) > 0 for rows in (inventory_rows, priority_rows, action_rows, summary_rows)
        )

        output_paths = {
            "ambiguity_resolution_run_directory": ambiguity_resolution_run_dir,
            "inventory_tsv": inventory_path,
            "priority_tsv": priority_path,
            "actions_tsv": actions_path,
            "summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(
            ambiguity_resolution_run_id=ambiguity_resolution_run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            output_paths=output_paths,
        )

        completed_at = utc_now()
        run_log_payload = {
            "status": "completed",
            "ambiguity_resolution_run_id": ambiguity_resolution_run_id,
            "blueprint_run_id": str(workflow_inputs.blueprint_latest_pointer["blueprint_run_id"]),
            "shortlist_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["shortlist_run_id"]),
            "core_audit_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["core_audit_run_id"]),
            "clinical_parse_run_id": str(workflow_inputs.clinical_shortlist_latest_pointer["parse_run_id"]),
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
                "cohort_blueprint_latest_json": repo_relative(paths.blueprint_latest_pointer, paths.repo_root),
                "cohort_blueprint_required_fields_tsv": repo_relative(
                    workflow_inputs.input_paths["blueprint_required_tsv"],
                    paths.repo_root,
                ),
                "cohort_blueprint_optional_fields_tsv": repo_relative(
                    workflow_inputs.input_paths["blueprint_optional_tsv"],
                    paths.repo_root,
                ),
                "cohort_blueprint_deferred_fields_tsv": repo_relative(
                    workflow_inputs.input_paths["blueprint_deferred_tsv"],
                    paths.repo_root,
                ),
                "cohort_blueprint_join_path_tsv": repo_relative(
                    workflow_inputs.input_paths["blueprint_join_path_tsv"],
                    paths.repo_root,
                ),
                "cohort_blueprint_ambiguities_tsv": repo_relative(
                    workflow_inputs.input_paths["blueprint_ambiguities_tsv"],
                    paths.repo_root,
                ),
                "cohort_blueprint_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["blueprint_summary_tsv"],
                    paths.repo_root,
                ),
                "cohort_blueprint_run_log_json": repo_relative(
                    workflow_inputs.input_paths["blueprint_run_log_json"],
                    paths.repo_root,
                ),
                "clinical_shortlist_latest_json": repo_relative(
                    paths.clinical_shortlist_latest_pointer,
                    paths.repo_root,
                ),
                "clinical_shortlist_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_tsv"],
                    paths.repo_root,
                ),
                "clinical_shortlist_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_summary_tsv"],
                    paths.repo_root,
                ),
                "clinical_shortlist_by_table_tsv": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_by_table_tsv"],
                    paths.repo_root,
                ),
                "clinical_shortlist_run_log_json": repo_relative(
                    workflow_inputs.input_paths["clinical_shortlist_run_log_json"],
                    paths.repo_root,
                ),
                "endpoint_crosswalk_latest_json": repo_relative(
                    paths.endpoint_crosswalk_latest_pointer,
                    paths.repo_root,
                ),
                "endpoint_candidate_inventory_tsv": repo_relative(
                    workflow_inputs.input_paths["endpoint_candidate_inventory_tsv"],
                    paths.repo_root,
                ),
                "endpoint_crosswalk_tsv": repo_relative(
                    workflow_inputs.input_paths["endpoint_crosswalk_tsv"],
                    paths.repo_root,
                ),
                "endpoint_crosswalk_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["endpoint_crosswalk_summary_tsv"],
                    paths.repo_root,
                ),
                "endpoint_crosswalk_run_log_json": repo_relative(
                    workflow_inputs.input_paths["endpoint_crosswalk_run_log_json"],
                    paths.repo_root,
                ),
                "biospecimen_crosswalk_latest_json": repo_relative(
                    paths.biospecimen_crosswalk_latest_pointer,
                    paths.repo_root,
                ),
                "biospecimen_identifier_inventory_tsv": repo_relative(
                    workflow_inputs.input_paths["biospecimen_inventory_tsv"],
                    paths.repo_root,
                ),
                "biospecimen_identifier_crosswalk_tsv": repo_relative(
                    workflow_inputs.input_paths["biospecimen_crosswalk_tsv"],
                    paths.repo_root,
                ),
                "biospecimen_identifier_crosswalk_summary_tsv": repo_relative(
                    workflow_inputs.input_paths["biospecimen_crosswalk_summary_tsv"],
                    paths.repo_root,
                ),
                "biospecimen_identifier_crosswalk_run_log_json": repo_relative(
                    workflow_inputs.input_paths["biospecimen_crosswalk_run_log_json"],
                    paths.repo_root,
                ),
            },
            "outputs": {
                "ambiguity_resolution_run_directory": repo_relative(
                    ambiguity_resolution_run_dir,
                    paths.repo_root,
                ),
                "ambiguity_resolution_inventory_tsv": repo_relative(inventory_path, paths.repo_root),
                "ambiguity_resolution_priority_tsv": repo_relative(priority_path, paths.repo_root),
                "ambiguity_resolution_actions_tsv": repo_relative(actions_path, paths.repo_root),
                "ambiguity_resolution_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": {
                "passed": (
                    output_rows_positive
                    and inventory_row_count_matches_blueprint
                    and every_ambiguity_has_one_linked_action
                    and action_linkage_covers_each_ambiguity_once
                    and priority_rows_match_action_rows
                    and summary_counts_reconcile_to_inventory
                    and blocker_class_values_valid
                    and resolution_priority_values_valid
                    and evidence_sufficiency_values_valid
                    and recommended_action_type_values_valid
                ),
                "blueprint_latest_pointer_found": True,
                "blueprint_run_log_completed": True,
                "clinical_shortlist_latest_pointer_found": True,
                "clinical_shortlist_run_log_completed": True,
                "endpoint_crosswalk_latest_pointer_found": True,
                "endpoint_crosswalk_run_log_completed": True,
                "biospecimen_crosswalk_latest_pointer_found": True,
                "biospecimen_crosswalk_run_log_completed": True,
                "blueprint_required_row_count_positive": len(workflow_inputs.blueprint_required_rows) > 0,
                "blueprint_optional_row_count_positive": len(workflow_inputs.blueprint_optional_rows) > 0,
                "blueprint_deferred_row_count_positive": len(workflow_inputs.blueprint_deferred_rows) > 0,
                "blueprint_join_path_row_count_positive": len(workflow_inputs.blueprint_join_path_rows) > 0,
                "blueprint_ambiguity_row_count_positive": len(workflow_inputs.blueprint_ambiguity_rows) > 0,
                "blueprint_summary_row_count_positive": len(workflow_inputs.blueprint_summary_rows) > 0,
                "inventory_row_count_positive": len(inventory_rows) > 0,
                "priority_row_count_positive": len(priority_rows) > 0,
                "actions_row_count_positive": len(action_rows) > 0,
                "summary_row_count_positive": len(summary_rows) > 0,
                "inventory_row_count_matches_blueprint_ambiguities": inventory_row_count_matches_blueprint,
                "every_ambiguity_maps_to_exactly_one_linked_action": (
                    every_ambiguity_has_one_linked_action and action_linkage_covers_each_ambiguity_once
                ),
                "action_linkage_covers_each_ambiguity_once": action_linkage_covers_each_ambiguity_once,
                "priority_rows_match_action_rows": priority_rows_match_action_rows,
                "summary_counts_reconcile_to_inventory": summary_counts_reconcile_to_inventory,
                "blocker_class_values_valid": blocker_class_values_valid,
                "resolution_priority_values_valid": resolution_priority_values_valid,
                "evidence_sufficiency_values_valid": evidence_sufficiency_values_valid,
                "recommended_action_type_values_valid": recommended_action_type_values_valid,
                "no_prior_run_overwrite": True,
                "latest_pointer_written_after_success_only": True,
            },
            "rules": {
                "minimal_build_scope": (
                    "lean minimal build: patient/case unit, sample anchor only, no endpoint freeze, "
                    "and treatment/timing fields allowed to be deferred or excluded"
                ),
                "blocker_class_order": list(BLOCKER_CLASS_ORDER),
                "resolution_priority_order": list(RESOLUTION_PRIORITY_ORDER),
                "evidence_sufficiency_order": list(EVIDENCE_SUFFICIENCY_ORDER),
                "recommended_action_type_order": list(RECOMMENDED_ACTION_TYPE_ORDER),
                "action_order": list(ACTION_ORDER),
                "notes_placeholder": NOTES_PLACEHOLDER,
                "classification_rules": {
                    ambiguity_id: {
                        "blocker_class": rule.blocker_class,
                        "resolution_priority": rule.resolution_priority,
                        "evidence_sufficiency": rule.evidence_sufficiency,
                        "recommended_action_type": rule.recommended_action_type,
                        "linked_action_id": rule.linked_action_id,
                    }
                    for ambiguity_id, rule in CLASSIFICATION_RULES.items()
                },
            },
            "counts": {
                "blueprint_ambiguity_count": len(workflow_inputs.blueprint_ambiguity_rows),
                "ambiguity_resolution_inventory_count": len(inventory_rows),
                "ambiguity_resolution_priority_count": len(priority_rows),
                "ambiguity_resolution_actions_count": len(action_rows),
                "ambiguity_resolution_summary_count": len(summary_rows),
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "cohort_blueprint_latest_pointer": workflow_inputs.blueprint_latest_pointer,
                "cohort_blueprint_run_log": workflow_inputs.blueprint_run_log,
                "clinical_shortlist_latest_pointer": workflow_inputs.clinical_shortlist_latest_pointer,
                "clinical_shortlist_run_log": workflow_inputs.clinical_shortlist_run_log,
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
            "ambiguity_resolution_run_id": ambiguity_resolution_run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_blueprint_ambiguity_resolution",
        }
        write_failure_log(run_log_path, failure_payload)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    print("TCGA-BRCA blueprint ambiguity-resolution workflow complete.")
    print(f"Ambiguity-resolution run ID: {run_log['ambiguity_resolution_run_id']}")
    print(f"Blueprint run ID: {run_log['blueprint_run_id']}")
    print(
        "Ambiguity-resolution output directory: "
        f"{run_log['outputs']['ambiguity_resolution_run_directory']}"
    )
    print(f"Inventory TSV: {run_log['outputs']['ambiguity_resolution_inventory_tsv']}")
    print(f"Priority TSV: {run_log['outputs']['ambiguity_resolution_priority_tsv']}")
    print(f"Actions TSV: {run_log['outputs']['ambiguity_resolution_actions_tsv']}")
    print(f"Summary TSV: {run_log['outputs']['ambiguity_resolution_summary_tsv']}")
    print(f"Latest pointer: {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
