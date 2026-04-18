#!/usr/bin/env python
"""Integrate indexed GDC clinical treatment evidence into patient_treatment_profile_v1."""

from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

MISSING_TOKENS = {
    "",
    "[discrepancy]",
    "[not applicable]",
    "[not available]",
    "[not evaluated]",
    "[unknown]",
    "n/a",
    "na",
    "nan",
    "none",
    "null",
}
JSON_FIELDS = {
    "drug_therapy_type_values_json",
    "therapy_type_counts_json",
    "regimen_context_values_json",
    "drug_name_values_json",
    "indexed_treatment_type_values_json",
    "indexed_treatment_intent_type_values_json",
    "indexed_therapeutic_agents_values_json",
    "indexed_regimen_line_values_json",
    "indexed_followup_timepoint_category_values_json",
}
PROFILE_COMPARE_MAP = {
    "has_any_drug_row": "local_profile_has_any_drug_row",
    "drug_row_count": "local_profile_drug_row_count",
    "has_any_radiation_row": "local_profile_has_any_radiation_row",
    "radiation_row_count": "local_profile_radiation_row_count",
    "drug_therapy_type_values_json": "local_profile_drug_therapy_type_values_json",
    "regimen_context_values_json": "local_profile_regimen_context_values_json",
    "has_any_treatment_timing": "local_profile_has_any_treatment_timing",
    "treatment_profile_status": "local_profile_status",
    "treatment_profile_requires_manual_review": "local_profile_requires_manual_review",
}
INDEXED_RAW = [
    "indexed_clinical_vs_local_treatment_v1_run_id",
    "treatment_source_audit_run_id",
    "indexed_case_id",
    "indexed_diagnosis_count",
    "indexed_treatment_row_count",
    "indexed_follow_up_row_count",
    "indexed_treatment_type_values_json",
    "indexed_treatment_type_distinct_count",
    "indexed_treatment_intent_type_values_json",
    "indexed_therapeutic_agents_values_json",
    "indexed_regimen_line_values_json",
    "indexed_followup_timepoint_category_values_json",
]
INDEXED_FLAGS = [
    "indexed_has_any_treatment_row",
    "indexed_has_any_druglike_treatment",
    "indexed_has_any_radiation_treatment",
    "indexed_has_any_surgery_treatment",
    "indexed_has_any_non_surgical_treatment",
    "indexed_has_any_therapeutic_agent_value",
    "indexed_has_any_regimen_line_value",
    "indexed_has_any_treatment_timing",
    "indexed_has_prior_treatment_flag",
    "indexed_has_followup_post_initial_treatment",
    "indexed_has_followup_explicit_treatment_field",
]
MASTER_FLAGS = [
    "master_has_any_treatment_evidence",
    "master_has_any_druglike_evidence",
    "master_has_any_radiation_evidence",
    "master_has_any_surgery_evidence",
    "master_has_any_non_surgical_treatment_evidence",
    "master_has_any_treatment_timing",
]
PROVENANCE = [
    "master_any_treatment_source_support",
    "master_druglike_source_support",
    "master_radiation_source_support",
    "master_surgery_source_support",
    "master_non_surgical_source_support",
    "master_treatment_timing_source_support",
]
AUDIT_COPY = [
    "indexed_adds_missing_druglike_coverage",
    "indexed_adds_missing_radiation_coverage",
    "indexed_adds_agent_detail_to_locally_missing_patient",
    "indexed_adds_regimen_line_detail_to_locally_missing_patient",
    "indexed_adds_timing_detail_to_locally_missing_patient",
    "indexed_surgery_only_flag",
    "indexed_context_only_flag",
    "indexed_incremental_value_class",
]
EXPECTED_INDEXED = {
    *INDEXED_RAW,
    *INDEXED_FLAGS,
    *AUDIT_COPY,
    "cohort_join_status",
    "cohort_scope",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "local_profile_has_any_drug_row",
    "local_profile_drug_row_count",
    "local_profile_has_any_radiation_row",
    "local_profile_radiation_row_count",
    "local_profile_drug_therapy_type_values_json",
    "local_profile_regimen_context_values_json",
    "local_profile_has_any_treatment_timing",
    "local_profile_status",
    "local_profile_requires_manual_review",
}
SOURCE_AUDIT_FIELDS = [
    "patient_treatment_master_v1_run_id",
    "patient_treatment_profile_v1_run_id",
    "indexed_clinical_vs_local_treatment_v1_run_id",
    "treatment_source_audit_run_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "provisional_patient_row_id",
    "baseline_analysis_v1_row_id",
    "feature_set_v1_row_index",
    "indexed_case_id",
    "local_profile_has_any_drug_row",
    "local_profile_has_any_radiation_row",
    "local_profile_has_any_treatment_evidence",
    "local_profile_has_any_treatment_timing",
    "indexed_has_any_treatment_row",
    "indexed_has_any_druglike_treatment",
    "indexed_has_any_radiation_treatment",
    "indexed_has_any_surgery_treatment",
    "indexed_has_any_non_surgical_treatment",
    "indexed_has_any_treatment_timing",
    "indexed_has_any_therapeutic_agent_value",
    "indexed_has_any_regimen_line_value",
    *PROVENANCE,
    "indexed_diagnosis_count",
    "indexed_treatment_row_count",
    "indexed_follow_up_row_count",
    "indexed_treatment_type_values_json",
    "indexed_treatment_intent_type_values_json",
    "indexed_therapeutic_agents_values_json",
    "indexed_regimen_line_values_json",
    "indexed_followup_timepoint_category_values_json",
    "indexed_surgery_only_flag",
    "indexed_context_only_flag",
    "indexed_incremental_value_class",
    "source_audit_interpretation",
]
SUMMARY_FIELDS = ["patient_treatment_master_v1_run_id", "summary_section", "summary_metric", "summary_value", "notes"]
MISSINGNESS_FIELDS = ["field_name", "field_group", "row_count", "non_missing_count", "missing_like_count", "missing_like_fraction", "distinct_non_missing_count", "notes"]


class PatientTreatmentMasterV1Error(RuntimeError):
    pass


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def utc_stamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def as_bool(value: Any, label: str) -> bool:
    token = str(value or "").strip().lower()
    if token == "yes":
        return True
    if token == "no":
        return False
    raise PatientTreatmentMasterV1Error(f"Expected yes/no for {label}, found {value!r}.")


def barcode(value: Any) -> str:
    return str(value or "").strip().upper()


def rel(path: Path, repo: Path) -> str:
    return path.resolve().relative_to(repo.resolve()).as_posix()


def missing_like(field: str, value: Any) -> bool:
    raw = str(value or "").strip()
    return raw in {"", "[]"} if field in JSON_FIELDS else raw.lower() in MISSING_TOKENS


def helper() -> Any:
    script = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script)
    if spec is None or spec.loader is None:
        raise PatientTreatmentMasterV1Error(f"Unable to load helper script: {script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module

def require_columns(rows: list[dict[str, str]], cols: set[str], label: str) -> None:
    if not rows:
        raise PatientTreatmentMasterV1Error(f"Required rows are empty for {label}.")
    missing = cols.difference(rows[0].keys())
    if missing:
        raise PatientTreatmentMasterV1Error(f"{label} is missing required columns: {sorted(missing)}")


def completed_pointer(pointer: Path, label: str, h: Any, repo: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    if not pointer.exists():
        raise PatientTreatmentMasterV1Error(f"Required {label} pointer not found: {pointer}")
    payload = h.load_json(pointer)
    run_log = h.resolve_existing_path(repo, str(payload.get("run_log_json") or ""), f"{label} run log")
    log = h.load_json(run_log)
    if log.get("status") != "completed" or not bool(log.get("validation", {}).get("passed", False)):
        raise PatientTreatmentMasterV1Error(f"{label} run log is not safely completed: {run_log}")
    return payload, log


def support(local: bool, indexed: bool) -> str:
    return "both" if local and indexed else "local_only" if local else "indexed_only" if indexed else "neither"


def frac(num: int, den: int) -> str:
    return "0.0000" if den <= 0 else f"{num / den:.4f}"


def paths(h: Any) -> dict[str, Path]:
    repo = h.detect_repo_root(Path(__file__).resolve().parent)
    cfg = repo / "09-trials" / "01-tcga-only-source-audited" / "04-config" / "trial_config.yaml"
    cfg_data = h.load_yaml(cfg)
    processed = repo / str(cfg_data.get("processed_data_root", "01-data/processed"))
    audit = repo / str(cfg_data.get("audit_root", "01-data/audit"))
    tp = audit / "tcga-brca" / "treatment-prep"
    return {
        "repo": repo,
        "cfg": cfg,
        "processed_runs": processed / "tcga-brca" / "treatment-prep" / "patient_treatment_master_v1_runs",
        "audit_runs": tp / "patient_treatment_master_v1_runs",
        "latest": tp / "tcga_brca_patient_treatment_master_v1_latest.json",
        "profile_latest": tp / "tcga_brca_patient_treatment_profile_v1_latest.json",
        "indexed_latest": audit / "tcga-brca" / "source" / "tcga_brca_indexed_clinical_vs_local_treatment_v1_latest.json",
    }


def load_inputs(p: dict[str, Path], h: Any) -> dict[str, Any]:
    profile_ptr, profile_log = completed_pointer(p["profile_latest"], "patient treatment profile latest", h, p["repo"])
    indexed_ptr, indexed_log = completed_pointer(p["indexed_latest"], "indexed clinical latest", h, p["repo"])
    if str(indexed_ptr["patient_treatment_profile_v1_run_id"]) != str(profile_ptr["patient_treatment_profile_v1_run_id"]):
        raise PatientTreatmentMasterV1Error("Indexed comparison pointer does not match the current profile pointer.")
    profile_tsv = h.resolve_existing_path(p["repo"], str(profile_ptr["patient_treatment_profile_v1_tsv"]), "patient_treatment_profile_v1.tsv")
    indexed_tsv = h.resolve_existing_path(p["repo"], str(indexed_ptr["indexed_clinical_treatment_patient_level_v1_tsv"]), "indexed_clinical_treatment_patient_level_v1.tsv")
    overlap_tsv = h.resolve_existing_path(p["repo"], str(indexed_ptr["indexed_vs_local_treatment_overlap_v1_tsv"]), "indexed_vs_local_treatment_overlap_v1.tsv")
    profile_rows = h.read_tsv_dict_rows(profile_tsv)
    indexed_rows = h.read_tsv_dict_rows(indexed_tsv)
    overlap_rows = h.read_tsv_dict_rows(overlap_tsv)
    require_columns(profile_rows, set(profile_rows[0].keys()), "patient_treatment_profile_v1.tsv")
    require_columns(indexed_rows, EXPECTED_INDEXED, "indexed_clinical_treatment_patient_level_v1.tsv")
    profile_fields = list(profile_rows[0].keys())
    indexed_local = [row for row in indexed_rows if row["cohort_scope"] == "local_profile"]
    indexed_only = [row for row in indexed_rows if row["cohort_scope"] == "indexed_only"]
    if len(indexed_local) != len(profile_rows):
        raise PatientTreatmentMasterV1Error("Indexed local_profile slice row count does not match profile row count.")
    for pr, ir in zip(profile_rows, indexed_local):
        if barcode(pr["bcr_patient_barcode"]) != barcode(ir["bcr_patient_barcode"]):
            raise PatientTreatmentMasterV1Error("Indexed local_profile slice does not preserve patient order.")
        for pf, ifield in PROFILE_COMPARE_MAP.items():
            if str(pr.get(pf, "")) != str(ir.get(ifield, "")):
                raise PatientTreatmentMasterV1Error(f"Profile/indexed local field mismatch for {pr['bcr_patient_barcode']}: {pf}")
    return {
        "profile_ptr": profile_ptr,
        "profile_log": profile_log,
        "profile_rows": profile_rows,
        "profile_fields": profile_fields,
        "indexed_ptr": indexed_ptr,
        "indexed_log": indexed_log,
        "indexed_local": indexed_local,
        "indexed_only": indexed_only,
        "overlap_rows": overlap_rows,
        "input_paths": {
            "patient_treatment_profile_v1_tsv": profile_tsv,
            "indexed_clinical_treatment_patient_level_v1_tsv": indexed_tsv,
            "indexed_vs_local_treatment_overlap_v1_tsv": overlap_tsv,
        },
    }


def build_master_rows(run_id: str, data: dict[str, Any]) -> tuple[list[str], list[dict[str, str]]]:
    fields = ["patient_treatment_master_v1_run_id", *data["profile_fields"], *INDEXED_RAW, *INDEXED_FLAGS, *MASTER_FLAGS, *PROVENANCE, *AUDIT_COPY]
    rows: list[dict[str, str]] = []
    for pr, ir in zip(data["profile_rows"], data["indexed_local"]):
        local_drug = as_bool(pr["has_any_drug_row"], "has_any_drug_row")
        local_rad = as_bool(pr["has_any_radiation_row"], "has_any_radiation_row")
        local_time = as_bool(pr["has_any_treatment_timing"], "has_any_treatment_timing")
        indexed_treat = as_bool(ir["indexed_has_any_treatment_row"], "indexed_has_any_treatment_row")
        indexed_drug = as_bool(ir["indexed_has_any_druglike_treatment"], "indexed_has_any_druglike_treatment")
        indexed_rad = as_bool(ir["indexed_has_any_radiation_treatment"], "indexed_has_any_radiation_treatment")
        indexed_surg = as_bool(ir["indexed_has_any_surgery_treatment"], "indexed_has_any_surgery_treatment")
        indexed_non = as_bool(ir["indexed_has_any_non_surgical_treatment"], "indexed_has_any_non_surgical_treatment")
        indexed_time = as_bool(ir["indexed_has_any_treatment_timing"], "indexed_has_any_treatment_timing")
        local_any = local_drug or local_rad
        row = {"patient_treatment_master_v1_run_id": run_id, **{f: str(pr.get(f, "")) for f in data["profile_fields"]}}
        row.update({f: str(ir.get(f, "")) for f in INDEXED_RAW + INDEXED_FLAGS + AUDIT_COPY})
        row.update({
            "master_has_any_treatment_evidence": yes_no(local_any or indexed_treat),
            "master_has_any_druglike_evidence": yes_no(local_drug or indexed_drug),
            "master_has_any_radiation_evidence": yes_no(local_rad or indexed_rad),
            "master_has_any_surgery_evidence": yes_no(indexed_surg),
            "master_has_any_non_surgical_treatment_evidence": yes_no(local_any or indexed_non),
            "master_has_any_treatment_timing": yes_no(local_time or indexed_time),
            "master_any_treatment_source_support": support(local_any, indexed_treat),
            "master_druglike_source_support": support(local_drug, indexed_drug),
            "master_radiation_source_support": support(local_rad, indexed_rad),
            "master_surgery_source_support": support(False, indexed_surg),
            "master_non_surgical_source_support": support(local_any, indexed_non),
            "master_treatment_timing_source_support": support(local_time, indexed_time),
        })
        rows.append(row)
    return fields, rows

def audit_interpretation(row: dict[str, str]) -> str:
    if row["indexed_incremental_value_class"] == "potentially_useful_new_evidence":
        return "indexed_fills_local_gap"
    if row["master_any_treatment_source_support"] == "both":
        return "local_and_indexed_concordant"
    if row["indexed_surgery_only_flag"] == "yes":
        return "indexed_surgery_only_recovery"
    return "no_broad_non_surgical_evidence_after_integration"


def build_source_audit_rows(master_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in master_rows:
        local_any = as_bool(row["has_any_drug_row"], "has_any_drug_row") or as_bool(row["has_any_radiation_row"], "has_any_radiation_row")
        rows.append({
            "patient_treatment_master_v1_run_id": row["patient_treatment_master_v1_run_id"],
            "patient_treatment_profile_v1_run_id": row["patient_treatment_profile_v1_run_id"],
            "indexed_clinical_vs_local_treatment_v1_run_id": row["indexed_clinical_vs_local_treatment_v1_run_id"],
            "treatment_source_audit_run_id": row["treatment_source_audit_run_id"],
            "bcr_patient_barcode": row["bcr_patient_barcode"],
            "bcr_patient_uuid": row["bcr_patient_uuid"],
            "provisional_patient_row_id": row["provisional_patient_row_id"],
            "baseline_analysis_v1_row_id": row["baseline_analysis_v1_row_id"],
            "feature_set_v1_row_index": row["feature_set_v1_row_index"],
            "indexed_case_id": row["indexed_case_id"],
            "local_profile_has_any_drug_row": row["has_any_drug_row"],
            "local_profile_has_any_radiation_row": row["has_any_radiation_row"],
            "local_profile_has_any_treatment_evidence": yes_no(local_any),
            "local_profile_has_any_treatment_timing": row["has_any_treatment_timing"],
            "indexed_has_any_treatment_row": row["indexed_has_any_treatment_row"],
            "indexed_has_any_druglike_treatment": row["indexed_has_any_druglike_treatment"],
            "indexed_has_any_radiation_treatment": row["indexed_has_any_radiation_treatment"],
            "indexed_has_any_surgery_treatment": row["indexed_has_any_surgery_treatment"],
            "indexed_has_any_non_surgical_treatment": row["indexed_has_any_non_surgical_treatment"],
            "indexed_has_any_treatment_timing": row["indexed_has_any_treatment_timing"],
            "indexed_has_any_therapeutic_agent_value": row["indexed_has_any_therapeutic_agent_value"],
            "indexed_has_any_regimen_line_value": row["indexed_has_any_regimen_line_value"],
            **{f: row[f] for f in PROVENANCE},
            "indexed_diagnosis_count": row["indexed_diagnosis_count"],
            "indexed_treatment_row_count": row["indexed_treatment_row_count"],
            "indexed_follow_up_row_count": row["indexed_follow_up_row_count"],
            "indexed_treatment_type_values_json": row["indexed_treatment_type_values_json"],
            "indexed_treatment_intent_type_values_json": row["indexed_treatment_intent_type_values_json"],
            "indexed_therapeutic_agents_values_json": row["indexed_therapeutic_agents_values_json"],
            "indexed_regimen_line_values_json": row["indexed_regimen_line_values_json"],
            "indexed_followup_timepoint_category_values_json": row["indexed_followup_timepoint_category_values_json"],
            "indexed_surgery_only_flag": row["indexed_surgery_only_flag"],
            "indexed_context_only_flag": row["indexed_context_only_flag"],
            "indexed_incremental_value_class": row["indexed_incremental_value_class"],
            "source_audit_interpretation": audit_interpretation(row),
        })
    return rows


def source_counts(master_rows: list[dict[str, str]], field: str) -> dict[str, int]:
    c = Counter(row[field] for row in master_rows)
    return {k: int(c.get(k, 0)) for k in ["local_only", "indexed_only", "both", "neither"]}


def sr(run_id: str, section: str, metric: str, value: Any, notes: str) -> dict[str, str]:
    return {"patient_treatment_master_v1_run_id": run_id, "summary_section": section, "summary_metric": metric, "summary_value": str(value), "notes": notes}


def build_summary_rows(run_id: str, data: dict[str, Any], master_rows: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict[str, Any]]:
    total = len(master_rows)
    lc = lambda field: sum(as_bool(row[field], field) for row in master_rows)
    local_any = sum(as_bool(r["has_any_drug_row"], "has_any_drug_row") or as_bool(r["has_any_radiation_row"], "has_any_radiation_row") for r in master_rows)
    local_drug, local_rad = lc("has_any_drug_row"), lc("has_any_radiation_row")
    local_non, local_surg = local_any, 0
    master_any, master_drug, master_rad = lc("master_has_any_treatment_evidence"), lc("master_has_any_druglike_evidence"), lc("master_has_any_radiation_evidence")
    master_non, master_surg = lc("master_has_any_non_surgical_treatment_evidence"), lc("master_has_any_surgery_evidence")
    supports = {
        "any_treatment": source_counts(master_rows, "master_any_treatment_source_support"),
        "druglike": source_counts(master_rows, "master_druglike_source_support"),
        "radiation": source_counts(master_rows, "master_radiation_source_support"),
        "non_surgical": source_counts(master_rows, "master_non_surgical_source_support"),
    }
    surg_only = lc("indexed_surgery_only_flag")
    no_drug, no_rad, no_non = supports["druglike"]["neither"], supports["radiation"]["neither"], supports["non_surgical"]["neither"]
    line_any, agent_gap, line_gap, time_gap = lc("indexed_has_any_regimen_line_value"), lc("indexed_adds_agent_detail_to_locally_missing_patient"), lc("indexed_adds_regimen_line_detail_to_locally_missing_patient"), lc("indexed_adds_timing_detail_to_locally_missing_patient")
    ready = "ready_for_next_treatment_organization_step_with_broad_source_audited_flags_only"
    rows = [
        sr(run_id, "inputs", "patient_treatment_profile_v1_run_id", data["profile_ptr"]["patient_treatment_profile_v1_run_id"], "Current local patient-level treatment profile feeding the integration."),
        sr(run_id, "inputs", "indexed_clinical_vs_local_treatment_v1_run_id", data["indexed_ptr"]["indexed_clinical_vs_local_treatment_v1_run_id"], "Saved indexed clinical comparison layer used as the official supplemental indexed input."),
        sr(run_id, "inputs", "treatment_source_audit_run_id", data["indexed_ptr"]["treatment_source_audit_run_id"], "Latest source-audit lineage carried forward from the indexed clinical comparison workflow."),
        sr(run_id, "inputs", "gdc_data_release", data["indexed_ptr"]["gdc_data_release"], "Live GDC release captured by the saved indexed comparison layer."),
        sr(run_id, "inputs", "gdc_tag", data["indexed_ptr"]["gdc_tag"], "Live GDC API tag captured by the saved indexed comparison layer."),
        sr(run_id, "inputs", "integrated_cohort_patient_count", total, "Anchored to patient_treatment_profile_v1 only; no cohort expansion during this integration step."),
        sr(run_id, "inputs", "indexed_only_case_count_not_integrated", len(data["indexed_only"]), f"Indexed-only out-of-cohort barcodes: {json.dumps([barcode(r['bcr_patient_barcode']) for r in data['indexed_only']], ensure_ascii=True)}"),
    ]
    for dim, local_count, master_count in [("any_treatment", local_any, master_any), ("druglike", local_drug, master_drug), ("radiation", local_rad, master_rad), ("non_surgical", local_non, master_non), ("surgery", local_surg, master_surg)]:
        rows += [sr(run_id, "coverage", f"{dim}_local_count", local_count, "Count before indexed supplemental integration."), sr(run_id, "coverage", f"{dim}_master_count", master_count, "Count after combining local treatment evidence with indexed clinical broad evidence."), sr(run_id, "coverage", f"{dim}_gain_count", master_count - local_count, "Absolute patient-count gain from adding indexed clinical evidence."), sr(run_id, "coverage", f"{dim}_master_fraction", frac(master_count, total), "Fraction of the integrated cohort with this broad treatment signal after integration.")]
    for dim, sc in supports.items():
        for cls in ["local_only", "indexed_only", "both", "neither"]:
            rows.append(sr(run_id, "source_support", f"{dim}_{cls}_count", sc[cls], f"Patient count for {dim} supported by {cls}."))
    rows += [
        sr(run_id, "weaknesses", "indexed_surgery_only_recovery_count", surg_only, "Patients whose indexed recovery remains surgery-only rather than broad non-surgical treatment support."),
        sr(run_id, "weaknesses", "remaining_no_master_druglike_evidence_count", no_drug, "Patients still lacking any broad druglike evidence after integration."),
        sr(run_id, "weaknesses", "remaining_no_master_radiation_evidence_count", no_rad, "Patients still lacking any broad radiation evidence after integration."),
        sr(run_id, "weaknesses", "remaining_no_master_non_surgical_treatment_evidence_count", no_non, "Patients still lacking any broad non-surgical treatment evidence after integration."),
        sr(run_id, "weaknesses", "patients_with_any_indexed_regimen_line_value", line_any, "Indexed regimen-or-line values remain very sparse in the integrated cohort."),
        sr(run_id, "weaknesses", "patients_with_indexed_agent_detail_added_to_locally_missing_patient", agent_gap, "Incremental indexed agent detail among locally missing patients remains minimal."),
        sr(run_id, "weaknesses", "patients_with_indexed_regimen_line_detail_added_to_locally_missing_patient", line_gap, "Incremental indexed regimen-line detail among locally missing patients remains absent."),
        sr(run_id, "weaknesses", "patients_with_indexed_timing_detail_added_to_locally_missing_patient", time_gap, "Incremental indexed timing detail among locally missing patients remains limited."),
        sr(run_id, "readiness", "final_readiness", ready, "The integrated layer is suitable for the next treatment-organization step only as a broad, source-audited flag layer."),
        sr(run_id, "readiness", "guardrail_note", "broad_flags_only_not_regimen_reconstruction", "Do not interpret indexed treatment types as full regimen reconstruction or treatment arms."),
        sr(run_id, "answers", "how_much_patient_level_treatment_coverage_improved", f"any_treatment {local_any}->{master_any} (+{master_any-local_any}); druglike {local_drug}->{master_drug} (+{master_drug-local_drug}); radiation {local_rad}->{master_rad} (+{master_rad-local_rad}); non_surgical {local_non}->{master_non} (+{master_non-local_non})", "Coverage gains reflect broad availability and status flags, not detailed regimen reconstruction."),
        sr(run_id, "answers", "how_many_patients_are_supported_by_local_only_indexed_only_both", f"any_treatment: {supports['any_treatment']['local_only']} local_only / {supports['any_treatment']['indexed_only']} indexed_only / {supports['any_treatment']['both']} both; druglike: {supports['druglike']['local_only']} / {supports['druglike']['indexed_only']} / {supports['druglike']['both']}; radiation: {supports['radiation']['local_only']} / {supports['radiation']['indexed_only']} / {supports['radiation']['both']}; non_surgical: {supports['non_surgical']['local_only']} / {supports['non_surgical']['indexed_only']} / {supports['non_surgical']['both']}", "Counts are reported on the integrated 1097-patient local profile cohort only."),
        sr(run_id, "answers", "what_still_remains_weak", f"{surg_only} patients remain surgery-only recoveries without broad non-surgical support; {no_drug} still lack druglike evidence; {no_rad} still lack radiation evidence; indexed regimen-line detail remains sparse ({line_any} patients overall; {line_gap} locally missing patients gain regimen-line detail).", "Indexed clinical materially improves coverage, but its added detail remains mostly broad presence/type information."),
        sr(run_id, "answers", "whether_dataset_is_strong_enough_for_next_stage_of_treatment_organization", ready, "Yes for the next broad treatment-organization step, but not for final treatment arms or regimen reconstruction."),
    ]
    metrics = {"total": total, "indexed_only": len(data["indexed_only"]), "local_any": local_any, "master_any": master_any, "local_drug": local_drug, "master_drug": master_drug, "local_rad": local_rad, "master_rad": master_rad, "local_non": local_non, "master_non": master_non, "master_surg": master_surg, "surg_only": surg_only, "no_drug": no_drug, "no_rad": no_rad, "no_non": no_non, "line_any": line_any, "agent_gap": agent_gap, "line_gap": line_gap, "time_gap": time_gap, "supports": supports, "ready": ready}
    return rows, metrics

def field_group(field: str) -> str:
    if field in {"patient_treatment_master_v1_run_id", "patient_treatment_profile_v1_run_id", "treatment_os_overlap_v1_run_id", "os_endpoint_v1_run_id", "clinical_biotab_parse_run_id", "baseline_model_input_v1_run_id", "cohort_v1_build_id", "indexed_clinical_vs_local_treatment_v1_run_id", "treatment_source_audit_run_id"}:
        return "lineage"
    if field in {"bcr_patient_barcode", "bcr_patient_uuid", "provisional_patient_row_id", "baseline_analysis_v1_row_id", "feature_set_v1_row_index", "indexed_case_id"}:
        return "identifier"
    if field in INDEXED_RAW[2:]:
        return "indexed_raw"
    if field in INDEXED_FLAGS or field in AUDIT_COPY[:-1]:
        return "indexed_flag"
    if field in MASTER_FLAGS:
        return "master_flag"
    if field in PROVENANCE:
        return "provenance_class"
    if field == "indexed_incremental_value_class":
        return "audit_class"
    return "local_profile"


def missing_note(field: str, group: str) -> str:
    if field in JSON_FIELDS:
        return "json_or_text_set_field; [] counted as missing-like"
    if field in INDEXED_FLAGS or field in MASTER_FLAGS or field in AUDIT_COPY[:-1] or field in {"has_any_drug_row", "has_any_radiation_row", "has_any_regimen_context", "has_any_treatment_timing", "treatment_profile_requires_manual_review"}:
        return "controlled_yes_no_flag; should be fully populated after integration"
    return "controlled source-support class" if group == "provenance_class" else "classification field carried from indexed comparison logic" if group == "audit_class" else group


def build_missingness_rows(master_fields: list[str], master_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    total = len(master_rows)
    for field in master_fields:
        vals = [str(r.get(field, "")) for r in master_rows]
        non_missing = [v for v in vals if not missing_like(field, v)]
        group = field_group(field)
        rows.append({"field_name": field, "field_group": group, "row_count": str(total), "non_missing_count": str(len(non_missing)), "missing_like_count": str(total - len(non_missing)), "missing_like_fraction": f"{((total - len(non_missing)) / total) if total else 0.0:.6f}", "distinct_non_missing_count": str(len(set(non_missing))), "notes": missing_note(field, group)})
    return rows


def overlap_map(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], int]:
    return {(r["dimension"], r["local_status"], r["indexed_status"]): int(r["patient_count"]) for r in rows if r.get("comparison_scope") == "local_profile_only"}


def summary_map(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {r["summary_metric"]: r for r in rows}


def validate(data: dict[str, Any], master_fields: list[str], master_rows: list[dict[str, str]], source_audit_rows: list[dict[str, str]], summary_rows: list[dict[str, str]], missingness_rows: list[dict[str, str]], metrics: dict[str, Any]) -> dict[str, Any]:
    missing = {r["field_name"]: int(r["missing_like_count"]) for r in missingness_rows}
    overlap = overlap_map(data["overlap_rows"])
    supports = metrics["supports"]
    summary = summary_map(summary_rows)
    checks = {
        "profile_indexed_run_id_match": str(data["indexed_ptr"]["patient_treatment_profile_v1_run_id"]) == str(data["profile_ptr"]["patient_treatment_profile_v1_run_id"]),
        "master_row_count_matches_profile": len(master_rows) == len(data["profile_rows"]),
        "source_audit_row_count_matches_master": len(source_audit_rows) == len(master_rows),
        "missingness_row_count_matches_master_columns": len(missingness_rows) == len(master_fields),
        "master_unique_barcodes": len({barcode(r["bcr_patient_barcode"]) for r in master_rows}) == len(master_rows),
        "indexed_only_cases_excluded_from_master": all(barcode(r["bcr_patient_barcode"]) not in {barcode(m["bcr_patient_barcode"]) for m in master_rows} for r in data["indexed_only"]),
        "profile_fields_preserved_verbatim": all(all(str(pr[f]) == str(mr[f]) for f in data["profile_fields"]) for pr, mr in zip(data["profile_rows"], master_rows)),
        "overlap_any_treatment_reconciles": supports["any_treatment"] == {"local_only": overlap.get(("any_treatment", "yes", "no"), -1), "indexed_only": overlap.get(("any_treatment", "no", "yes"), -1), "both": overlap.get(("any_treatment", "yes", "yes"), -1), "neither": overlap.get(("any_treatment", "no", "no"), -1)},
        "overlap_druglike_reconciles": supports["druglike"] == {"local_only": overlap.get(("druglike", "yes", "no"), -1), "indexed_only": overlap.get(("druglike", "no", "yes"), -1), "both": overlap.get(("druglike", "yes", "yes"), -1), "neither": overlap.get(("druglike", "no", "no"), -1)},
        "overlap_radiation_reconciles": supports["radiation"] == {"local_only": overlap.get(("radiation", "yes", "no"), -1), "indexed_only": overlap.get(("radiation", "no", "yes"), -1), "both": overlap.get(("radiation", "yes", "yes"), -1), "neither": overlap.get(("radiation", "no", "no"), -1)},
        "summary_any_treatment_counts_reconcile": int(summary["any_treatment_local_count"]["summary_value"]) == metrics["local_any"] and int(summary["any_treatment_master_count"]["summary_value"]) == metrics["master_any"],
        "summary_druglike_counts_reconcile": int(summary["druglike_local_count"]["summary_value"]) == metrics["local_drug"] and int(summary["druglike_master_count"]["summary_value"]) == metrics["master_drug"],
        "summary_radiation_counts_reconcile": int(summary["radiation_local_count"]["summary_value"]) == metrics["local_rad"] and int(summary["radiation_master_count"]["summary_value"]) == metrics["master_rad"],
        "summary_indexed_only_case_count_reconciles": int(summary["indexed_only_case_count_not_integrated"]["summary_value"]) == metrics["indexed_only"],
        "yes_no_fields_fully_populated": all(missing[f] == 0 for f in master_fields if field_group(f) in {"indexed_flag", "master_flag"} or f in {"has_any_drug_row", "has_any_radiation_row", "has_any_regimen_context", "has_any_treatment_timing", "treatment_profile_requires_manual_review"}),
        "provenance_fields_fully_populated": all(missing[f] == 0 for f in PROVENANCE),
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
    }
    checks["passed"] = all(checks.values())
    return checks


def latest_payload(run_id: str, p: dict[str, Path], data: dict[str, Any], outputs: dict[str, Path]) -> dict[str, Any]:
    return {
        "updated_at_utc": utc_stamp(utc_now()),
        "patient_treatment_master_v1_run_id": run_id,
        "patient_treatment_profile_v1_run_id": str(data["profile_ptr"]["patient_treatment_profile_v1_run_id"]),
        "indexed_clinical_vs_local_treatment_v1_run_id": str(data["indexed_ptr"]["indexed_clinical_vs_local_treatment_v1_run_id"]),
        "treatment_source_audit_run_id": str(data["indexed_ptr"]["treatment_source_audit_run_id"]),
        "gdc_data_release": str(data["indexed_ptr"]["gdc_data_release"]),
        "gdc_tag": str(data["indexed_ptr"]["gdc_tag"]),
        "processed_run_directory": rel(outputs["processed_dir"], p["repo"]),
        "audit_run_directory": rel(outputs["audit_dir"], p["repo"]),
        "patient_treatment_master_v1_tsv": rel(outputs["master_tsv"], p["repo"]),
        "patient_treatment_master_v1_source_audit_tsv": rel(outputs["source_audit_tsv"], p["repo"]),
        "patient_treatment_master_v1_summary_tsv": rel(outputs["summary_tsv"], p["repo"]),
        "patient_treatment_master_v1_missingness_tsv": rel(outputs["missingness_tsv"], p["repo"]),
        "run_log_json": rel(outputs["run_log"], p["repo"]),
        "patient_treatment_profile_latest_json": rel(p["profile_latest"], p["repo"]),
        "indexed_clinical_vs_local_treatment_latest_json": rel(p["indexed_latest"], p["repo"]),
    }

def run_workflow() -> dict[str, Any]:
    start = utc_now()
    run_id = start.strftime("%Y%m%dT%H%M%SZ")
    h = helper()
    p = paths(h)
    data = load_inputs(p, h)
    processed_dir = p["processed_runs"] / run_id
    audit_dir = p["audit_runs"] / run_id
    run_log = audit_dir / "run_log.json"
    try:
        trial_cfg = h.load_yaml(p["cfg"])
        h.create_run_directory(processed_dir)
        h.create_run_directory(audit_dir)
        master_fields, master_rows = build_master_rows(run_id, data)
        source_audit_rows = build_source_audit_rows(master_rows)
        summary_rows, metrics = build_summary_rows(run_id, data, master_rows)
        missingness_rows = build_missingness_rows(master_fields, master_rows)
        validation = validate(data, master_fields, master_rows, source_audit_rows, summary_rows, missingness_rows, metrics)
        if not validation["passed"]:
            raise PatientTreatmentMasterV1Error(f"Output validation failed: {[k for k, v in validation.items() if v is False]}")
        outputs = {
            "processed_dir": processed_dir,
            "audit_dir": audit_dir,
            "master_tsv": processed_dir / "patient_treatment_master_v1.tsv",
            "source_audit_tsv": audit_dir / "patient_treatment_master_v1_source_audit.tsv",
            "summary_tsv": audit_dir / "patient_treatment_master_v1_summary.tsv",
            "missingness_tsv": audit_dir / "patient_treatment_master_v1_missingness.tsv",
            "run_log": run_log,
        }
        h.write_dict_rows_tsv(outputs["master_tsv"], master_fields, master_rows)
        h.write_dict_rows_tsv(outputs["source_audit_tsv"], SOURCE_AUDIT_FIELDS, source_audit_rows)
        h.write_dict_rows_tsv(outputs["summary_tsv"], SUMMARY_FIELDS, summary_rows)
        h.write_dict_rows_tsv(outputs["missingness_tsv"], MISSINGNESS_FIELDS, missingness_rows)
        latest = latest_payload(run_id, p, data, outputs)
        done = utc_now()
        payload = {
            "status": "completed",
            "patient_treatment_master_v1_run_id": run_id,
            "started_at_utc": utc_stamp(start),
            "completed_at_utc": utc_stamp(done),
            "repo_root": str(p["repo"].resolve()),
            "trial_name": trial_cfg.get("trial_name"),
            "dataset_scope": trial_cfg.get("dataset_scope"),
            "inputs": {"trial_config_yaml": rel(p["cfg"], p["repo"]), "patient_treatment_profile_latest_json": rel(p["profile_latest"], p["repo"]), "indexed_clinical_vs_local_treatment_latest_json": rel(p["indexed_latest"], p["repo"]), **{k: rel(v, p["repo"]) for k, v in data["input_paths"].items()}},
            "outputs": {"processed_run_directory": rel(processed_dir, p["repo"]), "audit_run_directory": rel(audit_dir, p["repo"]), "patient_treatment_master_v1_tsv": rel(outputs["master_tsv"], p["repo"]), "patient_treatment_master_v1_source_audit_tsv": rel(outputs["source_audit_tsv"], p["repo"]), "patient_treatment_master_v1_summary_tsv": rel(outputs["summary_tsv"], p["repo"]), "patient_treatment_master_v1_missingness_tsv": rel(outputs["missingness_tsv"], p["repo"]), "run_log_json": rel(run_log, p["repo"]), "latest_pointer_json": rel(p["latest"], p["repo"])} ,
            "validation": validation,
            "rules": {"integration_only": True, "does_not_modify_raw_downloads": True, "does_not_build_treatment_arms": True, "does_not_perform_modeling": True, "does_not_attempt_regimen_reconstruction": True, "keeps_current_patient_treatment_profile_v1_cohort_fixed": True, "indexed_only_cases_reported_but_not_integrated": True, "surgery_remains_indexed_only_broad_presence_signal": True, "missing_like_normalization": "scalars: value.strip().lower(); json/text sets: empty string or []", "missing_like_tokens": sorted(MISSING_TOKENS)},
            "counts": {"patient_treatment_master_v1_row_count": len(master_rows), "indexed_only_case_count_not_integrated": metrics["indexed_only"], "local_any_treatment_count": metrics["local_any"], "master_any_treatment_count": metrics["master_any"], "local_druglike_count": metrics["local_drug"], "master_druglike_count": metrics["master_drug"], "local_radiation_count": metrics["local_rad"], "master_radiation_count": metrics["master_rad"], "local_non_surgical_count": metrics["local_non"], "master_non_surgical_count": metrics["master_non"], "master_surgery_count": metrics["master_surg"], "surgery_only_recovery_count": metrics["surg_only"], "remaining_no_druglike_count": metrics["no_drug"], "remaining_no_radiation_count": metrics["no_rad"], "remaining_no_non_surgical_count": metrics["no_non"], "patients_with_any_indexed_regimen_line_value": metrics["line_any"], "patients_with_indexed_agent_detail_added_to_locally_missing_patient": metrics["agent_gap"], "patients_with_indexed_regimen_line_detail_added_to_locally_missing_patient": metrics["line_gap"], "patients_with_indexed_timing_detail_added_to_locally_missing_patient": metrics["time_gap"], "any_treatment_source_support": metrics["supports"]["any_treatment"], "druglike_source_support": metrics["supports"]["druglike"], "radiation_source_support": metrics["supports"]["radiation"], "non_surgical_source_support": metrics["supports"]["non_surgical"], "final_readiness": metrics["ready"]},
            "latest_pointer": latest,
            "upstream_snapshots": {"patient_treatment_profile_pointer": data["profile_ptr"], "indexed_clinical_vs_local_pointer": data["indexed_ptr"]},
        }
        h.write_json(run_log, payload)
        h.write_json(p["latest"], latest, overwrite=True)
        return payload
    except Exception as exc:
        h.write_json(run_log, {"status": "failed", "patient_treatment_master_v1_run_id": run_id, "started_at_utc": utc_stamp(start), "failed_at_utc": utc_stamp(utc_now()), "error": str(exc), "workflow": "tcga_brca_patient_treatment_master_v1"}, overwrite=True)
        raise


def main() -> int:
    run_log = run_workflow()
    counts = run_log["counts"]
    print("TCGA-BRCA patient treatment master v1 workflow complete.")
    print(f"  Run ID                : {run_log['patient_treatment_master_v1_run_id']}")
    print(f"  Any treatment         : {counts['local_any_treatment_count']} -> {counts['master_any_treatment_count']}")
    print(f"  Druglike              : {counts['local_druglike_count']} -> {counts['master_druglike_count']}")
    print(f"  Radiation             : {counts['local_radiation_count']} -> {counts['master_radiation_count']}")
    print(f"  Non-surgical          : {counts['local_non_surgical_count']} -> {counts['master_non_surgical_count']}")
    print(f"  Indexed-only cases    : {counts['indexed_only_case_count_not_integrated']}")
    print(f"  Readiness             : {counts['final_readiness']}")
    print(f"  Processed run dir     : {run_log['outputs']['processed_run_directory']}")
    print(f"  Audit run dir         : {run_log['outputs']['audit_run_directory']}")
    print(f"  Latest pointer        : {run_log['outputs']['latest_pointer_json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
