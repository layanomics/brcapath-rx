#!/usr/bin/env python
"""Audit treatment–OS overlap for TCGA-BRCA: verify treatment record overlap with the OS endpoint cohort."""

from __future__ import annotations

import importlib.util
import json
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

MISSING_LIKE_TOKENS = {
    "",
    "[not available]",
    "[not applicable]",
    "[unknown]",
    "[not evaluated]",
    "[discrepancy]",
    "na",
    "n/a",
    "null",
    "none",
    "nan",
}

OVERLAP_SUMMARY_FIELDNAMES = [
    "treatment_os_overlap_v1_run_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "provisional_patient_row_id",
    "os_event",
    "os_time_days",
    "os_requires_manual_review",
    "os_endpoint_inclusion_status",
    "has_drug_row",
    "drug_row_count",
    "has_radiation_row",
    "radiation_row_count",
    "drug_therapy_type_values_json",
    "drug_therapy_type_distinct_count",
    "drug_therapy_type_single_or_mixed",
    "dominant_therapy_type_if_any",
    "age_at_diagnosis",
    "er_status_by_ihc",
    "pr_status_by_ihc",
    "her2_status_by_ihc",
    "ajcc_pathologic_tumor_stage",
    "histological_type",
]

OVERLAP_COUNTS_FIELDNAMES = [
    "treatment_os_overlap_v1_run_id",
    "category",
    "subcategory",
    "count",
    "notes",
]

DRUG_THERAPY_TYPE_DISTRIBUTION_FIELDNAMES = [
    "treatment_os_overlap_v1_run_id",
    "bcr_patient_barcode",
    "drug_row_count",
    "non_missing_therapy_type_count",
    "distinct_non_missing_therapy_type_count",
    "therapy_type_values_json",
    "dominant_therapy_type_if_any",
    "drug_therapy_type_single_or_mixed",
]

DRUG_THERAPY_TYPE_PATIENT_SUMMARY_FIELDNAMES = [
    "treatment_os_overlap_v1_run_id",
    "group",
    "patient_count",
    "notes",
]

AUDIT_SUMMARY_FIELDNAMES = [
    "treatment_os_overlap_v1_run_id",
    "metric",
    "value",
    "notes",
]


class TreatmentOSOverlapV1Error(RuntimeError):
    """Raised when the treatment–OS overlap audit v1 workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the treatment–OS overlap audit workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    os_endpoint_v1_latest_pointer: Path
    clinical_biotabs_latest_pointer: Path
    baseline_model_input_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    os_endpoint_v1_latest_pointer: dict[str, Any]
    os_endpoint_v1_run_log: dict[str, Any]
    clinical_biotabs_latest_pointer: dict[str, Any]
    clinical_biotabs_run_log: dict[str, Any]
    baseline_model_input_latest_pointer: dict[str, Any]
    baseline_model_input_run_log: dict[str, Any]
    os_endpoint_rows: list[dict[str, str]]
    clinical_drug_rows: list[dict[str, str]]
    clinical_radiation_rows: list[dict[str, str]]
    clinical_patient_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class OverlapMetrics:
    """Scalar metrics derived from the patient-level overlap table."""

    total: int
    patients_with_drug: int
    patients_without_drug: int
    patients_with_radiation: int
    patients_without_radiation: int
    patients_single_therapy_type: int
    patients_mixed_therapy_type: int
    patients_missing_therapy_type_with_drug: int
    drug_row_count_total: int
    drug_row_counts_all: list[int]
    radiation_row_count_total: int
    radiation_row_counts_all: list[int]


# ---------------------------------------------------------------------------
# Local utility functions (mirrors script 18 pattern)
# ---------------------------------------------------------------------------

def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def json_list(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=True)


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def is_missing_like(value: str) -> bool:
    return value.strip().lower() in MISSING_LIKE_TOKENS


def therapy_type_slug(raw_value: str) -> str:
    """Produce a stable group-name slug from a raw therapy_type value."""
    return raw_value.strip().lower().replace(" ", "_").replace("-", "_")


# ---------------------------------------------------------------------------
# Helper module loading
# ---------------------------------------------------------------------------

def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise TreatmentOSOverlapV1Error(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise TreatmentOSOverlapV1Error(f"Unable to create an import spec for: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


# ---------------------------------------------------------------------------
# Path and input loading
# ---------------------------------------------------------------------------

def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_config = (
        repo_root / "09-trials" / "01-tcga-only-source-audited" / "04-config" / "trial_config.yaml"
    )
    if not trial_config.exists():
        raise TreatmentOSOverlapV1Error(f"Required trial config not found: {trial_config}")
    trial_config_data = helper_module.load_yaml(trial_config)
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    treatment_prep_root = audit_root / "tcga-brca" / "treatment-prep"
    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        audit_runs_root=treatment_prep_root / "treatment_os_overlap_v1_runs",
        latest_pointer=treatment_prep_root / "tcga_brca_treatment_os_overlap_v1_latest.json",
        os_endpoint_v1_latest_pointer=(
            audit_root / "tcga-brca" / "endpoint-prep" / "tcga_brca_os_endpoint_v1_latest.json"
        ),
        clinical_biotabs_latest_pointer=(
            audit_root / "tcga-brca" / "variables" / "tcga_brca_clinical_biotabs_latest.json"
        ),
        baseline_model_input_latest_pointer=(
            audit_root / "tcga-brca" / "model-input" / "tcga_brca_baseline_model_input_v1_latest.json"
        ),
    )


def _require_columns(rows: list[dict[str, str]], required: set[str], source: str) -> None:
    if not rows:
        raise TreatmentOSOverlapV1Error(f"{source} has no rows to validate columns against.")
    actual = set(rows[0].keys())
    missing = required - actual
    if missing:
        raise TreatmentOSOverlapV1Error(
            f"{source} is missing required columns: {sorted(missing)}"
        )


def _require_completed_run_log(
    repo_root: Path,
    run_log_relative: str,
    label: str,
    helper_module: Any,
) -> dict[str, Any]:
    run_log_path = helper_module.resolve_existing_path(repo_root, run_log_relative, f"{label} run log")
    run_log = helper_module.load_json(run_log_path)
    if run_log.get("status") != "completed":
        raise TreatmentOSOverlapV1Error(f"{label} run log is not completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise TreatmentOSOverlapV1Error(f"{label} run log does not report validation.passed == true.")
    return run_log


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    # Verify required pointer files exist
    for pointer_path, label in [
        (paths.os_endpoint_v1_latest_pointer, "OS endpoint v1 latest pointer"),
        (paths.clinical_biotabs_latest_pointer, "clinical biotabs latest pointer"),
        (paths.baseline_model_input_latest_pointer, "baseline model input v1 latest pointer"),
    ]:
        if not pointer_path.exists():
            raise TreatmentOSOverlapV1Error(f"Required {label} not found: {pointer_path}")

    # Load and validate pointer contents
    os_ep_pointer = helper_module.load_json(paths.os_endpoint_v1_latest_pointer)
    helper_module.require_keys(
        os_ep_pointer,
        {"os_endpoint_v1_run_id", "os_endpoint_v1_tsv", "run_log_json", "cohort_v1_build_id"},
        "OS endpoint v1 latest pointer",
        paths.os_endpoint_v1_latest_pointer,
    )
    biotabs_pointer = helper_module.load_json(paths.clinical_biotabs_latest_pointer)
    helper_module.require_keys(
        biotabs_pointer,
        {"parse_run_id", "source_run_id", "processed_run_directory", "run_log_json"},
        "clinical biotabs latest pointer",
        paths.clinical_biotabs_latest_pointer,
    )
    bmi_pointer = helper_module.load_json(paths.baseline_model_input_latest_pointer)
    helper_module.require_keys(
        bmi_pointer,
        {"baseline_model_input_v1_run_id", "run_log_json", "cohort_v1_build_id"},
        "baseline model input v1 latest pointer",
        paths.baseline_model_input_latest_pointer,
    )

    # Require completed run logs
    os_ep_run_log = _require_completed_run_log(
        paths.repo_root, str(os_ep_pointer["run_log_json"]), "OS endpoint v1", helper_module
    )
    biotabs_run_log = _require_completed_run_log(
        paths.repo_root, str(biotabs_pointer["run_log_json"]), "clinical biotabs", helper_module
    )
    bmi_run_log = _require_completed_run_log(
        paths.repo_root, str(bmi_pointer["run_log_json"]), "baseline model input v1", helper_module
    )

    # Resolve TSV paths
    os_endpoint_tsv = helper_module.resolve_existing_path(
        paths.repo_root, str(os_ep_pointer["os_endpoint_v1_tsv"]), "os_endpoint_v1.tsv"
    )
    biotabs_proc_dir = str(biotabs_pointer["processed_run_directory"])
    clinical_drug_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (Path(biotabs_proc_dir) / "clinical_drug.tsv").as_posix(),
        "clinical_drug.tsv",
    )
    clinical_radiation_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (Path(biotabs_proc_dir) / "clinical_radiation.tsv").as_posix(),
        "clinical_radiation.tsv",
    )
    clinical_patient_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (Path(biotabs_proc_dir) / "clinical_patient.tsv").as_posix(),
        "clinical_patient.tsv",
    )

    # Load rows
    os_endpoint_rows = helper_module.read_tsv_dict_rows(os_endpoint_tsv)
    clinical_drug_rows = helper_module.read_tsv_dict_rows(clinical_drug_tsv)
    clinical_radiation_rows = helper_module.read_tsv_dict_rows(clinical_radiation_tsv)
    clinical_patient_rows = helper_module.read_tsv_dict_rows(clinical_patient_tsv)

    # Validate columns
    _require_columns(
        os_endpoint_rows,
        {
            "bcr_patient_barcode", "bcr_patient_uuid", "provisional_patient_row_id",
            "os_event", "os_time_days", "os_requires_manual_review",
            "os_endpoint_inclusion_status",
        },
        "os_endpoint_v1.tsv",
    )
    _require_columns(
        clinical_drug_rows,
        {"bcr_patient_barcode", "pharmaceutical_therapy_type"},
        "clinical_drug.tsv",
    )
    _require_columns(clinical_radiation_rows, {"bcr_patient_barcode"}, "clinical_radiation.tsv")
    _require_columns(
        clinical_patient_rows,
        {
            "bcr_patient_barcode", "bcr_patient_uuid", "age_at_diagnosis",
            "er_status_by_ihc", "pr_status_by_ihc", "her2_status_by_ihc",
            "ajcc_pathologic_tumor_stage", "histological_type",
        },
        "clinical_patient.tsv",
    )

    if not os_endpoint_rows:
        raise TreatmentOSOverlapV1Error("os_endpoint_v1.tsv has no data rows.")
    if not clinical_drug_rows:
        raise TreatmentOSOverlapV1Error("clinical_drug.tsv has no data rows.")
    if not clinical_radiation_rows:
        raise TreatmentOSOverlapV1Error("clinical_radiation.tsv has no data rows.")

    return WorkflowInputs(
        os_endpoint_v1_latest_pointer=os_ep_pointer,
        os_endpoint_v1_run_log=os_ep_run_log,
        clinical_biotabs_latest_pointer=biotabs_pointer,
        clinical_biotabs_run_log=biotabs_run_log,
        baseline_model_input_latest_pointer=bmi_pointer,
        baseline_model_input_run_log=bmi_run_log,
        os_endpoint_rows=os_endpoint_rows,
        clinical_drug_rows=clinical_drug_rows,
        clinical_radiation_rows=clinical_radiation_rows,
        clinical_patient_rows=clinical_patient_rows,
        input_paths={
            "os_endpoint_v1_tsv": os_endpoint_tsv,
            "clinical_drug_tsv": clinical_drug_tsv,
            "clinical_radiation_tsv": clinical_radiation_tsv,
            "clinical_patient_tsv": clinical_patient_tsv,
        },
    )


# ---------------------------------------------------------------------------
# Output builders
# ---------------------------------------------------------------------------

def build_overlap_summary_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
) -> tuple[list[dict[str, Any]], OverlapMetrics]:
    """Build one overlap summary row per OS cohort patient."""
    drug_by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in workflow_inputs.clinical_drug_rows:
        barcode = row["bcr_patient_barcode"].strip().upper()
        if barcode:
            drug_by_patient[barcode].append(row)

    rad_by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in workflow_inputs.clinical_radiation_rows:
        barcode = row["bcr_patient_barcode"].strip().upper()
        if barcode:
            rad_by_patient[barcode].append(row)

    patient_by_barcode: dict[str, dict[str, str]] = {}
    for row in workflow_inputs.clinical_patient_rows:
        barcode = row["bcr_patient_barcode"].strip().upper()
        if barcode:
            patient_by_barcode[barcode] = row

    overlap_rows: list[dict[str, Any]] = []
    drug_counts_all: list[int] = []
    rad_counts_all: list[int] = []

    for os_row in workflow_inputs.os_endpoint_rows:
        barcode = os_row["bcr_patient_barcode"].strip().upper()
        drug_pt = drug_by_patient.get(barcode, [])
        rad_pt = rad_by_patient.get(barcode, [])
        clin_row = patient_by_barcode.get(barcode, {})

        drug_count = len(drug_pt)
        rad_count = len(rad_pt)
        drug_counts_all.append(drug_count)
        rad_counts_all.append(rad_count)

        raw_types = [r["pharmaceutical_therapy_type"] for r in drug_pt]
        non_missing = [t for t in raw_types if not is_missing_like(t)]
        distinct_nm = list(dict.fromkeys(t.strip() for t in non_missing))

        if not drug_pt:
            single_or_mixed = ""
            dominant = ""
        elif not non_missing:
            single_or_mixed = "missing"
            dominant = ""
        elif len(set(t.strip() for t in non_missing)) == 1:
            single_or_mixed = "single"
            dominant = non_missing[0].strip()
        else:
            single_or_mixed = "mixed"
            top = Counter(t.strip() for t in non_missing).most_common(2)
            dominant = top[0][0] if (len(top) == 1 or top[0][1] > top[1][1]) else ""

        overlap_rows.append({
            "treatment_os_overlap_v1_run_id": run_id,
            "bcr_patient_barcode": os_row["bcr_patient_barcode"],
            "bcr_patient_uuid": os_row["bcr_patient_uuid"],
            "provisional_patient_row_id": os_row["provisional_patient_row_id"],
            "os_event": os_row["os_event"],
            "os_time_days": os_row["os_time_days"],
            "os_requires_manual_review": os_row["os_requires_manual_review"],
            "os_endpoint_inclusion_status": os_row["os_endpoint_inclusion_status"],
            "has_drug_row": yes_no(drug_count > 0),
            "drug_row_count": str(drug_count),
            "has_radiation_row": yes_no(rad_count > 0),
            "radiation_row_count": str(rad_count),
            "drug_therapy_type_values_json": json_list(raw_types),
            "drug_therapy_type_distinct_count": str(len(distinct_nm)),
            "drug_therapy_type_single_or_mixed": single_or_mixed,
            "dominant_therapy_type_if_any": dominant,
            "age_at_diagnosis": clin_row.get("age_at_diagnosis", ""),
            "er_status_by_ihc": clin_row.get("er_status_by_ihc", ""),
            "pr_status_by_ihc": clin_row.get("pr_status_by_ihc", ""),
            "her2_status_by_ihc": clin_row.get("her2_status_by_ihc", ""),
            "ajcc_pathologic_tumor_stage": clin_row.get("ajcc_pathologic_tumor_stage", ""),
            "histological_type": clin_row.get("histological_type", ""),
        })

    total = len(overlap_rows)
    patients_with_drug = sum(1 for r in overlap_rows if r["has_drug_row"] == "yes")
    patients_with_rad = sum(1 for r in overlap_rows if r["has_radiation_row"] == "yes")
    patients_single = sum(1 for r in overlap_rows if r["drug_therapy_type_single_or_mixed"] == "single")
    patients_mixed = sum(1 for r in overlap_rows if r["drug_therapy_type_single_or_mixed"] == "mixed")
    patients_missing_type = sum(
        1 for r in overlap_rows if r["drug_therapy_type_single_or_mixed"] == "missing"
    )

    metrics = OverlapMetrics(
        total=total,
        patients_with_drug=patients_with_drug,
        patients_without_drug=total - patients_with_drug,
        patients_with_radiation=patients_with_rad,
        patients_without_radiation=total - patients_with_rad,
        patients_single_therapy_type=patients_single,
        patients_mixed_therapy_type=patients_mixed,
        patients_missing_therapy_type_with_drug=patients_missing_type,
        drug_row_count_total=sum(drug_counts_all),
        drug_row_counts_all=drug_counts_all,
        radiation_row_count_total=sum(rad_counts_all),
        radiation_row_counts_all=rad_counts_all,
    )
    return overlap_rows, metrics


def build_overlap_counts_rows(
    run_id: str,
    overlap_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build long-form count summary rows."""
    rows: list[dict[str, Any]] = []

    def add(category: str, subcategory: str, count: int, notes: str = "") -> None:
        rows.append({
            "treatment_os_overlap_v1_run_id": run_id,
            "category": category,
            "subcategory": subcategory,
            "count": str(count),
            "notes": notes,
        })

    # has_drug_row × os_event
    for drug_val in ("yes", "no"):
        for event_val in ("1", "0"):
            add(
                "has_drug_row_x_os_event",
                f"has_drug_row={drug_val}_os_event={event_val}",
                sum(1 for r in overlap_rows
                    if r["has_drug_row"] == drug_val and r["os_event"] == event_val),
            )

    # has_radiation_row × os_event
    for rad_val in ("yes", "no"):
        for event_val in ("1", "0"):
            add(
                "has_radiation_row_x_os_event",
                f"has_radiation_row={rad_val}_os_event={event_val}",
                sum(1 for r in overlap_rows
                    if r["has_radiation_row"] == rad_val and r["os_event"] == event_val),
            )

    # drug_row_count distribution
    drug_freq = Counter(int(r["drug_row_count"]) for r in overlap_rows)
    for count_val in sorted(drug_freq.keys()):
        add("drug_row_count_distribution", f"drug_row_count={count_val}", drug_freq[count_val])

    # radiation_row_count distribution
    rad_freq = Counter(int(r["radiation_row_count"]) for r in overlap_rows)
    for count_val in sorted(rad_freq.keys()):
        add("radiation_row_count_distribution", f"radiation_row_count={count_val}", rad_freq[count_val])

    # drug_therapy_type_single_or_mixed distribution
    som_freq: Counter[str] = Counter()
    for r in overlap_rows:
        val = r["drug_therapy_type_single_or_mixed"] or "not_applicable"
        som_freq[val] += 1
    for som_val in sorted(som_freq.keys()):
        add(
            "drug_therapy_type_single_or_mixed",
            f"drug_therapy_type_single_or_mixed={som_val}",
            som_freq[som_val],
        )

    return rows


def build_drug_therapy_type_distribution_rows(
    run_id: str,
    overlap_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """One row per patient who has at least one drug row."""
    rows: list[dict[str, Any]] = []
    for r in overlap_rows:
        if r["has_drug_row"] != "yes":
            continue
        raw_types_list: list[str] = json.loads(r["drug_therapy_type_values_json"])
        non_missing = [t for t in raw_types_list if not is_missing_like(t)]
        distinct_nm = list(dict.fromkeys(t.strip() for t in non_missing))
        rows.append({
            "treatment_os_overlap_v1_run_id": run_id,
            "bcr_patient_barcode": r["bcr_patient_barcode"],
            "drug_row_count": r["drug_row_count"],
            "non_missing_therapy_type_count": str(len(non_missing)),
            "distinct_non_missing_therapy_type_count": str(len(distinct_nm)),
            "therapy_type_values_json": r["drug_therapy_type_values_json"],
            "dominant_therapy_type_if_any": r["dominant_therapy_type_if_any"],
            "drug_therapy_type_single_or_mixed": r["drug_therapy_type_single_or_mixed"],
        })
    return rows


def build_drug_therapy_type_patient_summary_rows(
    run_id: str,
    overlap_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Aggregate patient counts by dominant therapy type group."""
    dominant_counts: Counter[str] = Counter()
    mixed_no_dominant = 0
    missing_type_count = 0
    no_drug_count = 0

    for r in overlap_rows:
        if r["has_drug_row"] != "yes":
            no_drug_count += 1
        elif r["drug_therapy_type_single_or_mixed"] == "missing":
            missing_type_count += 1
        elif r["drug_therapy_type_single_or_mixed"] == "mixed" and not r["dominant_therapy_type_if_any"]:
            mixed_no_dominant += 1
        else:
            dominant_counts[r["dominant_therapy_type_if_any"]] += 1

    rows: list[dict[str, Any]] = []

    def add(group: str, count: int, notes: str = "") -> None:
        rows.append({
            "treatment_os_overlap_v1_run_id": run_id,
            "group": group,
            "patient_count": str(count),
            "notes": notes,
        })

    for dom_val in sorted(dominant_counts.keys()):
        add(
            f"dominant_therapy_type__{therapy_type_slug(dom_val)}",
            dominant_counts[dom_val],
            dom_val,
        )
    add("mixed_no_dominant_therapy_type", mixed_no_dominant,
        "has drug rows; mixed therapy_type values; no single dominant")
    add("missing_therapy_type_despite_drug_rows", missing_type_count,
        "has drug rows but all pharmaceutical_therapy_type values are missing-like")
    add("no_drug_rows", no_drug_count,
        "no rows in clinical_drug.tsv for this patient")

    return rows


def build_audit_summary_rows(
    run_id: str,
    overlap_rows: list[dict[str, Any]],
    metrics: OverlapMetrics,
) -> tuple[list[dict[str, Any]], str]:
    """Build long-form audit summary rows. Returns (rows, feasibility_interpretation)."""
    rows: list[dict[str, Any]] = []

    def add(metric: str, value: Any, notes: str = "") -> None:
        rows.append({
            "treatment_os_overlap_v1_run_id": run_id,
            "metric": metric,
            "value": str(value),
            "notes": notes,
        })

    add("total_os_cohort_size", metrics.total)
    add("os_event_1_count", sum(1 for r in overlap_rows if r["os_event"] == "1"))
    add("os_event_0_count", sum(1 for r in overlap_rows if r["os_event"] == "0"))
    add("patients_with_any_drug_row", metrics.patients_with_drug)
    add("patients_with_no_drug_rows", metrics.patients_without_drug)
    pct_drug = 100.0 * metrics.patients_with_drug / metrics.total if metrics.total > 0 else 0.0
    add("pct_patients_with_drug_row", f"{pct_drug:.1f}%")
    add("patients_with_any_radiation_row", metrics.patients_with_radiation)
    add("patients_with_no_radiation_rows", metrics.patients_without_radiation)
    add("drug_row_count_total", metrics.drug_row_count_total)
    if metrics.drug_row_counts_all:
        drug_pos = [c for c in metrics.drug_row_counts_all if c > 0]
        add("drug_row_count_median_all_patients",
            f"{statistics.median(metrics.drug_row_counts_all):.1f}")
        add("drug_row_count_median_patients_with_drug",
            f"{statistics.median(drug_pos):.1f}" if drug_pos else "0.0")
        add("drug_row_count_max", max(metrics.drug_row_counts_all))
    add("radiation_row_count_total", metrics.radiation_row_count_total)
    if metrics.radiation_row_counts_all:
        add("radiation_row_count_max", max(metrics.radiation_row_counts_all))
    add("patients_with_single_therapy_type", metrics.patients_single_therapy_type)
    add("patients_with_mixed_therapy_type", metrics.patients_mixed_therapy_type)
    add("patients_with_missing_therapy_type_despite_drug_rows",
        metrics.patients_missing_therapy_type_with_drug)
    patients_with_dominant = sum(1 for r in overlap_rows if r["dominant_therapy_type_if_any"])
    add("patients_with_dominant_therapy_type_identifiable", patients_with_dominant)

    pct_identifiable = 100.0 * patients_with_dominant / metrics.total if metrics.total > 0 else 0.0
    add("pct_patients_with_identifiable_dominant_therapy_type", f"{pct_identifiable:.1f}%")

    # Feasibility interpretation
    frac_drug = metrics.patients_with_drug / metrics.total if metrics.total > 0 else 0.0
    if frac_drug >= 0.80:
        feasibility = "ready_for_coarse_exploratory_grouping"
        feasibility_notes = (
            f"{pct_drug:.0f}% of OS cohort have drug rows. "
            f"Coverage sufficient for coarse exploratory grouping."
        )
    elif frac_drug >= 0.60:
        feasibility = "ready_for_coarse_exploratory_grouping"
        feasibility_notes = (
            f"{pct_drug:.0f}% of OS cohort have drug rows. "
            f"Coarse grouping feasible with a 'treatment_unknown' bucket for patients without drug records."
        )
    else:
        feasibility = "not_ready"
        feasibility_notes = (
            f"Only {pct_drug:.0f}% of OS cohort have drug rows. "
            f"Drug coverage too low for meaningful coarse grouping."
        )

    add("feasibility_interpretation", feasibility, feasibility_notes)
    add("coarse_grouping_next_step",
        "patient_level_treatment_aggregation",
        "Aggregate drug rows per patient, normalize therapy_type, resolve mixed cases before freezing treatment arms.")
    add("remaining_blocker",
        "drug_name_not_normalized",
        "pharmaceutical_therapy_drug_name not yet audited or normalized; needed before arm-level grouping.")

    return rows, feasibility


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_outputs(
    workflow_inputs: WorkflowInputs,
    overlap_rows: list[dict[str, Any]],
    metrics: OverlapMetrics,
) -> dict[str, Any]:
    os_barcodes = [r["bcr_patient_barcode"].strip().upper() for r in workflow_inputs.os_endpoint_rows]
    ov_barcodes = [r["bcr_patient_barcode"].strip().upper() for r in overlap_rows]
    all_present = set(os_barcodes) == set(ov_barcodes)
    no_dupes = len(ov_barcodes) == len(set(ov_barcodes))
    counts_ok = (
        metrics.patients_with_drug + metrics.patients_without_drug == metrics.total
        and metrics.patients_with_radiation + metrics.patients_without_radiation == metrics.total
    )
    passed = (
        len(workflow_inputs.os_endpoint_rows) > 0
        and len(workflow_inputs.clinical_drug_rows) > 0
        and len(workflow_inputs.clinical_radiation_rows) > 0
        and len(overlap_rows) == len(workflow_inputs.os_endpoint_rows)
        and all_present
        and no_dupes
        and counts_ok
    )
    return {
        "passed": passed,
        "required_upstream_pointers_found": True,
        "required_source_tables_found": True,
        "os_endpoint_run_log_completed": True,
        "clinical_biotabs_run_log_completed": True,
        "baseline_model_input_run_log_completed": True,
        "os_cohort_row_count_positive": len(workflow_inputs.os_endpoint_rows) > 0,
        "drug_row_count_positive": len(workflow_inputs.clinical_drug_rows) > 0,
        "radiation_row_count_positive": len(workflow_inputs.clinical_radiation_rows) > 0,
        "overlap_summary_row_count_matches_os_cohort": (
            len(overlap_rows) == len(workflow_inputs.os_endpoint_rows)
        ),
        "all_os_patients_in_overlap_summary": all_present,
        "no_overlap_summary_duplicate_barcodes": no_dupes,
        "patient_counts_reconcile": counts_ok,
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
        "os_cohort_row_count": len(workflow_inputs.os_endpoint_rows),
        "overlap_summary_row_count": len(overlap_rows),
        "drug_source_row_count": len(workflow_inputs.clinical_drug_rows),
        "radiation_source_row_count": len(workflow_inputs.clinical_radiation_rows),
    }


# ---------------------------------------------------------------------------
# Latest pointer and failure log
# ---------------------------------------------------------------------------

def build_latest_pointer_payload(
    run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    audit_run_dir: Path,
    output_file_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "treatment_os_overlap_v1_run_id": run_id,
        "os_endpoint_v1_run_id": str(
            workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"]
        ),
        "clinical_biotab_parse_run_id": str(
            workflow_inputs.clinical_biotabs_latest_pointer["parse_run_id"]
        ),
        "baseline_model_input_v1_run_id": str(
            workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
        ),
        "cohort_v1_build_id": str(
            workflow_inputs.os_endpoint_v1_latest_pointer["cohort_v1_build_id"]
        ),
        "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
        "treatment_os_overlap_summary_tsv": repo_relative(
            output_file_paths["treatment_os_overlap_summary"], paths.repo_root
        ),
        "treatment_os_overlap_counts_tsv": repo_relative(
            output_file_paths["treatment_os_overlap_counts"], paths.repo_root
        ),
        "drug_therapy_type_distribution_tsv": repo_relative(
            output_file_paths["drug_therapy_type_distribution"], paths.repo_root
        ),
        "drug_therapy_type_patient_summary_tsv": repo_relative(
            output_file_paths["drug_therapy_type_patient_summary"], paths.repo_root
        ),
        "treatment_os_overlap_audit_summary_tsv": repo_relative(
            output_file_paths["treatment_os_overlap_audit_summary"], paths.repo_root
        ),
        "run_log_json": repo_relative(output_file_paths["run_log"], paths.repo_root),
        "os_endpoint_v1_latest_json": repo_relative(
            paths.os_endpoint_v1_latest_pointer, paths.repo_root
        ),
        "clinical_biotabs_latest_json": repo_relative(
            paths.clinical_biotabs_latest_pointer, paths.repo_root
        ),
        "baseline_model_input_v1_latest_json": repo_relative(
            paths.baseline_model_input_latest_pointer, paths.repo_root
        ),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    audit_run_dir = paths.audit_runs_root / run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)

        helper_module.create_run_directory(audit_run_dir)

        output_file_paths: dict[str, Path] = {
            "treatment_os_overlap_summary": audit_run_dir / "treatment_os_overlap_summary.tsv",
            "treatment_os_overlap_counts": audit_run_dir / "treatment_os_overlap_counts.tsv",
            "drug_therapy_type_distribution": audit_run_dir / "drug_therapy_type_distribution.tsv",
            "drug_therapy_type_patient_summary": audit_run_dir / "drug_therapy_type_patient_summary.tsv",
            "treatment_os_overlap_audit_summary": audit_run_dir / "treatment_os_overlap_audit_summary.tsv",
            "run_log": run_log_path,
        }

        # Build all output tables
        overlap_rows, metrics = build_overlap_summary_rows(run_id, workflow_inputs)
        counts_rows = build_overlap_counts_rows(run_id, overlap_rows)
        distribution_rows = build_drug_therapy_type_distribution_rows(run_id, overlap_rows)
        patient_summary_rows = build_drug_therapy_type_patient_summary_rows(run_id, overlap_rows)
        audit_summary_rows, feasibility = build_audit_summary_rows(run_id, overlap_rows, metrics)

        # Validate before writing
        validation = validate_outputs(workflow_inputs, overlap_rows, metrics)
        if not validation["passed"]:
            raise TreatmentOSOverlapV1Error(
                f"Output validation did not pass. Failed checks: "
                + str({k: v for k, v in validation.items() if v is False})
            )

        # Write TSVs
        helper_module.write_dict_rows_tsv(
            output_file_paths["treatment_os_overlap_summary"],
            OVERLAP_SUMMARY_FIELDNAMES, overlap_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_file_paths["treatment_os_overlap_counts"],
            OVERLAP_COUNTS_FIELDNAMES, counts_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_file_paths["drug_therapy_type_distribution"],
            DRUG_THERAPY_TYPE_DISTRIBUTION_FIELDNAMES, distribution_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_file_paths["drug_therapy_type_patient_summary"],
            DRUG_THERAPY_TYPE_PATIENT_SUMMARY_FIELDNAMES, patient_summary_rows,
        )
        helper_module.write_dict_rows_tsv(
            output_file_paths["treatment_os_overlap_audit_summary"],
            AUDIT_SUMMARY_FIELDNAMES, audit_summary_rows,
        )

        latest_pointer_payload = build_latest_pointer_payload(
            run_id=run_id,
            paths=paths,
            workflow_inputs=workflow_inputs,
            audit_run_dir=audit_run_dir,
            output_file_paths=output_file_paths,
        )

        completed_at = utc_now()
        run_log_payload: dict[str, Any] = {
            "status": "completed",
            "treatment_os_overlap_v1_run_id": run_id,
            "os_endpoint_v1_run_id": str(
                workflow_inputs.os_endpoint_v1_latest_pointer["os_endpoint_v1_run_id"]
            ),
            "clinical_biotab_parse_run_id": str(
                workflow_inputs.clinical_biotabs_latest_pointer["parse_run_id"]
            ),
            "baseline_model_input_v1_run_id": str(
                workflow_inputs.baseline_model_input_latest_pointer["baseline_model_input_v1_run_id"]
            ),
            "cohort_v1_build_id": str(
                workflow_inputs.os_endpoint_v1_latest_pointer["cohort_v1_build_id"]
            ),
            "clinical_source_run_id": str(
                workflow_inputs.clinical_biotabs_latest_pointer["source_run_id"]
            ),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "os_endpoint_v1_latest_json": repo_relative(
                    paths.os_endpoint_v1_latest_pointer, paths.repo_root
                ),
                "clinical_biotabs_latest_json": repo_relative(
                    paths.clinical_biotabs_latest_pointer, paths.repo_root
                ),
                "baseline_model_input_v1_latest_json": repo_relative(
                    paths.baseline_model_input_latest_pointer, paths.repo_root
                ),
                **{
                    key: repo_relative(path, paths.repo_root)
                    for key, path in workflow_inputs.input_paths.items()
                },
            },
            "outputs": {
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "treatment_os_overlap_summary_tsv": repo_relative(
                    output_file_paths["treatment_os_overlap_summary"], paths.repo_root
                ),
                "treatment_os_overlap_counts_tsv": repo_relative(
                    output_file_paths["treatment_os_overlap_counts"], paths.repo_root
                ),
                "drug_therapy_type_distribution_tsv": repo_relative(
                    output_file_paths["drug_therapy_type_distribution"], paths.repo_root
                ),
                "drug_therapy_type_patient_summary_tsv": repo_relative(
                    output_file_paths["drug_therapy_type_patient_summary"], paths.repo_root
                ),
                "treatment_os_overlap_audit_summary_tsv": repo_relative(
                    output_file_paths["treatment_os_overlap_audit_summary"], paths.repo_root
                ),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation,
            "rules": {
                "unit_of_analysis": "patient/case",
                "treatment_scope": "clinical_drug_and_radiation_biotabs_only",
                "no_drug_name_normalization": True,
                "no_treatment_arm_aggregation": True,
                "no_modeling": True,
                "no_causal_analysis": True,
                "audit_only": True,
                "os_cohort_source": "os_endpoint_v1",
                "confounder_snapshot_source": "clinical_patient_biotab",
                "missing_like_normalization": "value.strip().lower()",
                "missing_like_tokens": sorted(MISSING_LIKE_TOKENS),
            },
            "counts": {
                "os_endpoint_row_count": len(workflow_inputs.os_endpoint_rows),
                "clinical_drug_row_count": len(workflow_inputs.clinical_drug_rows),
                "clinical_radiation_row_count": len(workflow_inputs.clinical_radiation_rows),
                "clinical_patient_row_count": len(workflow_inputs.clinical_patient_rows),
                "overlap_summary_row_count": len(overlap_rows),
                "overlap_counts_row_count": len(counts_rows),
                "drug_therapy_type_distribution_row_count": len(distribution_rows),
                "drug_therapy_type_patient_summary_row_count": len(patient_summary_rows),
                "audit_summary_row_count": len(audit_summary_rows),
                "patients_with_drug": metrics.patients_with_drug,
                "patients_without_drug": metrics.patients_without_drug,
                "patients_with_radiation": metrics.patients_with_radiation,
                "patients_without_radiation": metrics.patients_without_radiation,
                "patients_single_therapy_type": metrics.patients_single_therapy_type,
                "patients_mixed_therapy_type": metrics.patients_mixed_therapy_type,
                "patients_missing_therapy_type_with_drug": metrics.patients_missing_therapy_type_with_drug,
                "drug_row_count_total": metrics.drug_row_count_total,
                "radiation_row_count_total": metrics.radiation_row_count_total,
                "feasibility_interpretation": feasibility,
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "os_endpoint_v1_latest_pointer": workflow_inputs.os_endpoint_v1_latest_pointer,
                "clinical_biotabs_latest_pointer": workflow_inputs.clinical_biotabs_latest_pointer,
                "baseline_model_input_v1_latest_pointer": (
                    workflow_inputs.baseline_model_input_latest_pointer
                ),
            },
        }

        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload

    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "treatment_os_overlap_v1_run_id": run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error": str(exc),
            "workflow": "tcga_brca_treatment_os_overlap_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    counts = run_log.get("counts", {})
    print("TCGA-BRCA treatment–OS overlap v1 audit complete.")
    print(f"  Run ID              : {run_log['treatment_os_overlap_v1_run_id']}")
    print(f"  OS cohort size      : {counts.get('os_endpoint_row_count', '?')}")
    print(f"  Patients with drug  : {counts.get('patients_with_drug', '?')}")
    print(f"  Patients with rad   : {counts.get('patients_with_radiation', '?')}")
    print(f"  Single therapy type : {counts.get('patients_single_therapy_type', '?')}")
    print(f"  Mixed therapy type  : {counts.get('patients_mixed_therapy_type', '?')}")
    print(f"  Feasibility         : {counts.get('feasibility_interpretation', '?')}")
    print(f"  Audit run dir       : {run_log['outputs']['audit_run_directory']}")
    print(f"  Latest pointer      : {run_log['outputs']['latest_pointer_json']}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
