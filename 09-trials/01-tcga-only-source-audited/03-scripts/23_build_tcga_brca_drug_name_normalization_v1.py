#!/usr/bin/env python
"""Build an auditable TCGA-BRCA drug-name normalization and provisional drug-class layer v1."""

from __future__ import annotations

import importlib.util
import json
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

MISSING_LIKE_NORMALIZATION = "trim -> lowercase -> punctuation collapsed to spaces"
MISSING_LIKE_MATCH_KEYS = {
    "",
    "na",
    "n a",
    "nan",
    "none",
    "not applicable",
    "not available",
    "not evaluated",
    "null",
    "unknown",
}
PROVISIONAL_CLASS_VOCABULARY = [
    "alkylating_agent",
    "anthracycline",
    "taxane",
    "platinum",
    "antimetabolite",
    "endocrine_serm",
    "endocrine_aromatase_inhibitor",
    "endocrine_other",
    "her2_targeted",
    "immunotherapy",
    "ancillary_supportive",
    "other_cytotoxic",
    "unknown_or_review_needed",
]
NORMALIZATION_CONFIDENCE_ALLOWED = {"confident", "review_needed", "unmapped"}
CLASS_MAPPING_CONFIDENCE_ALLOWED = {"confident", "review_needed"}
MAPPING_STATUS_ALLOWED = {"mapped_confidently", "mapped_but_review", "unmapped"}
SINGLE_OR_MULTI_ALLOWED = {"unknown_only", "single_class", "multi_class"}
READINESS_NOT_READY = "not_ready"
READINESS_DRUG_CLASS_REVIEW = "ready_for_drug_class_review"
READINESS_TREATMENT_ARM_FREEZE_REVIEW_NEXT = "ready_for_treatment_arm_freeze_review_next"
TREATMENT_ARM_FREEZE_REVIEW_BLOCKED = "blocked"
TREATMENT_ARM_FREEZE_REVIEW_PROVISIONAL_NEXT = "provisional_review_next_only"
TREATMENT_RECOMMENDATION_MODELING_STATUS = "out_of_scope"
INVENTORY_FIELDNAMES = [
    "drug_name_normalization_v1_run_id",
    "raw_drug_name",
    "normalized_raw_drug_name_for_matching",
    "drug_row_count",
    "distinct_patient_count",
    "therapy_type_values_json",
    "regimen_context_values_json",
    "example_patient_barcodes_json",
    "mapping_status",
    "notes",
]
NORMALIZATION_MAP_FIELDNAMES = [
    "drug_name_normalization_v1_run_id",
    "raw_drug_name",
    "normalized_raw_drug_name_for_matching",
    "provisional_normalized_drug_name",
    "normalization_rule",
    "normalization_confidence",
    "manual_review_needed",
    "notes",
]
CLASS_MAP_FIELDNAMES = [
    "drug_name_normalization_v1_run_id",
    "provisional_normalized_drug_name",
    "provisional_drug_class",
    "class_mapping_rule",
    "class_mapping_confidence",
    "manual_review_needed",
    "notes",
]
PATIENT_PROFILE_FIELDNAMES = [
    "drug_name_normalization_v1_run_id",
    "patient_treatment_profile_v1_run_id",
    "patient_treatment_grouping_v1_run_id",
    "treatment_os_overlap_v1_run_id",
    "os_endpoint_v1_run_id",
    "clinical_biotab_parse_run_id",
    "clinical_source_run_id",
    "baseline_model_input_v1_run_id",
    "cohort_v1_build_id",
    "bcr_patient_barcode",
    "bcr_patient_uuid",
    "provisional_patient_row_id",
    "baseline_analysis_v1_row_id",
    "feature_set_v1_row_index",
    "drug_row_count",
    "drug_therapy_type_values_json",
    "regimen_context_values_json",
    "raw_drug_name_values_json",
    "normalized_drug_name_values_json",
    "provisional_drug_class_values_json",
    "distinct_raw_drug_name_count",
    "distinct_normalized_drug_name_count",
    "distinct_provisional_drug_class_count",
    "single_or_multi_drug_class",
    "dominant_provisional_drug_class_if_any",
    "drug_name_normalization_requires_manual_review",
    "drug_class_mapping_requires_manual_review",
    "drug_name_profile_flags_json",
]
SUMMARY_FIELDNAMES = [
    "drug_name_normalization_v1_run_id",
    "summary_section",
    "summary_metric",
    "summary_value",
    "notes",
]
PATIENT_PROFILE_FLAG_ORDER = [
    "contains_unmapped_raw_drug_name",
    "contains_review_needed_name_normalization",
    "contains_unknown_or_review_needed_class",
    "contains_missing_like_raw_drug_name",
    "contains_compound_or_regimen_like_raw_drug_name",
    "contains_single_known_class_plus_unknown",
    "contains_multiple_known_drug_classes",
]
REGIMEN_OR_PLACEHOLDER_MATCH_KEYS = {
    "ac",
    "chemo nos",
    "hormone nos",
    "not otherwise specified",
    "taxane",
    "tc",
    "tch",
}
SPACE_ONLY_COMPOUND_MATCH_KEYS = {
    "adriamycin cyclophosphamid",
    "adriamycin cuclophosphamide",
    "adrimicin cyclophosphamide",
    "adrimycin cyclophosphamide",
    "cyclophosphamide methotrexatum fluorouracillum",
    "doxorubicin cyclophosphamide",
    "doxorubicine cyclophosphamide",
    "doxorubicine cyclophosphamide tamoxifen",
    "methotrexate 5 fluorouracil cyclophosphamide",
    "tamoxiphen anastrazolum",
    "tamoxiphene anastrozolum",
    "tamoxiphene leuporeline gosereline",
    "taxol adriamycin cyclophosphamide herceptin",
}


class DrugNameNormalizationV1Error(RuntimeError):
    """Raised when the drug-name normalization workflow cannot complete safely."""


@dataclass(frozen=True)
class WorkflowPaths:
    """Concrete repo paths used by the drug-name normalization workflow."""

    repo_root: Path
    trial_config: Path
    results_root: Path
    processed_runs_root: Path
    audit_runs_root: Path
    latest_pointer: Path
    patient_treatment_profile_latest_pointer: Path
    patient_treatment_grouping_latest_pointer: Path
    clinical_biotabs_latest_pointer: Path


@dataclass(frozen=True)
class WorkflowInputs:
    """Resolved workflow inputs loaded from saved audit layers."""

    patient_treatment_profile_latest_pointer: dict[str, Any]
    patient_treatment_profile_run_log: dict[str, Any]
    patient_treatment_grouping_latest_pointer: dict[str, Any]
    patient_treatment_grouping_run_log: dict[str, Any]
    clinical_biotabs_latest_pointer: dict[str, Any]
    clinical_biotabs_run_log: dict[str, Any]
    patient_treatment_profile_rows: list[dict[str, str]]
    patient_treatment_grouping_rows: list[dict[str, str]]
    clinical_drug_rows: list[dict[str, str]]
    input_paths: dict[str, Path]


@dataclass(frozen=True)
class NormalizationDecision:
    """Deterministic normalization decision for one distinct raw drug name."""

    raw_drug_name: str
    normalized_raw_drug_name_for_matching: str
    provisional_normalized_drug_name: str
    normalization_rule: str
    normalization_confidence: str
    manual_review_needed: str
    mapping_status: str
    notes: str
    is_missing_like: bool
    is_compound_or_regimen_like: bool


@dataclass(frozen=True)
class ClassDecision:
    """Deterministic provisional class decision for one normalized drug name."""

    provisional_normalized_drug_name: str
    provisional_drug_class: str
    class_mapping_rule: str
    class_mapping_confidence: str
    manual_review_needed: str
    notes: str


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def format_utc_timestamp(moment: datetime) -> str:
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_helper_module() -> Any:
    script_path = Path(__file__).resolve().with_name("11_build_tcga_brca_minimal_dry_run_cohort.py")
    if not script_path.exists():
        raise DrugNameNormalizationV1Error(f"Required helper script not found: {script_path}")
    spec = importlib.util.spec_from_file_location("tcga_brca_minimal_dry_run_cohort", script_path)
    if spec is None or spec.loader is None:
        raise DrugNameNormalizationV1Error(f"Unable to create an import spec for: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


def repo_relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def normalize_barcode(value: str) -> str:
    return value.strip().upper()


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def json_list(values: list[Any]) -> str:
    return json.dumps(values, ensure_ascii=True)


def ordered_unique(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def parse_int_or_none(value: str) -> int | None:
    stripped = value.strip()
    if not stripped:
        return None
    try:
        return int(stripped)
    except ValueError:
        return None


def parse_json_list(raw_value: str, field_name: str) -> list[str]:
    if not str(raw_value).strip():
        return []
    parsed = json.loads(raw_value)
    if not isinstance(parsed, list):
        raise DrugNameNormalizationV1Error(f"Expected a JSON list in {field_name}: {raw_value!r}")
    return [str(value) for value in parsed]


def require_columns(rows: list[dict[str, str]], required_columns: set[str], label: str) -> None:
    if not rows:
        raise DrugNameNormalizationV1Error(f"Required rows are empty for {label}.")
    missing = required_columns.difference(rows[0].keys())
    if missing:
        raise DrugNameNormalizationV1Error(f"{label} is missing required columns: {sorted(missing)}")


def require_completed_run_log(
    *,
    repo_root: Path,
    run_log_relative_path: str,
    label: str,
    helper_module: Any,
) -> dict[str, Any]:
    run_log_path = helper_module.resolve_existing_path(repo_root, run_log_relative_path, label)
    run_log = helper_module.load_json(run_log_path)
    if run_log.get("status") != "completed":
        raise DrugNameNormalizationV1Error(f"{label} does not report status == completed.")
    if not bool(run_log.get("validation", {}).get("passed", False)):
        raise DrugNameNormalizationV1Error(f"{label} does not report validation.passed == true.")
    return run_log


def normalize_raw_drug_name_for_matching(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", " ", value.strip().lower())
    return re.sub(r"\s+", " ", normalized).strip()


def is_missing_like_match_key(match_key: str) -> bool:
    return match_key in MISSING_LIKE_MATCH_KEYS


def ordered_non_missing_distinct_values(raw_values: list[str]) -> list[str]:
    ordered: list[str] = []
    seen: set[str] = set()
    for value in raw_values:
        stripped = value.strip()
        if is_missing_like_match_key(normalize_raw_drug_name_for_matching(stripped)):
            continue
        if stripped in seen:
            continue
        seen.add(stripped)
        ordered.append(stripped)
    return ordered


def build_unique_lookup_by_barcode(
    rows: list[dict[str, str]],
    *,
    barcode_field: str,
    label: str,
) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        barcode = normalize_barcode(str(row.get(barcode_field, "")))
        if not barcode:
            raise DrugNameNormalizationV1Error(f"{label} contains an empty barcode value.")
        if barcode in lookup:
            raise DrugNameNormalizationV1Error(f"{label} contains a duplicate barcode: {barcode}")
        lookup[barcode] = row
    return lookup


def group_rows_by_barcode(
    rows: list[dict[str, str]],
    *,
    barcode_field: str,
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        barcode = normalize_barcode(str(row.get(barcode_field, "")))
        if barcode:
            grouped[barcode].append(row)
    return grouped


def register_normalization_entries(
    target: dict[str, tuple[str, str, str]],
    *,
    keys: list[str],
    normalized_name: str,
    rule: str,
    notes: str,
) -> None:
    for key in keys:
        if key in target:
            raise DrugNameNormalizationV1Error(f"Duplicate normalization key registered: {key}")
        target[key] = (normalized_name, rule, notes)


def register_class_entries(
    target: dict[str, tuple[str, str, str]],
    *,
    names: list[str],
    provisional_drug_class: str,
    rule: str,
    notes: str,
) -> None:
    for name in names:
        if name in target:
            raise DrugNameNormalizationV1Error(f"Duplicate class mapping registered: {name}")
        target[name] = (provisional_drug_class, rule, notes)


def build_confident_normalization_map() -> dict[str, tuple[str, str, str]]:
    mapping: dict[str, tuple[str, str, str]] = {}
    register_normalization_entries(mapping, keys=["zoladex", "goserelin"], normalized_name="goserelin", rule="brand_generic_equivalence_or_case_normalization", notes="Case-normalized exact values and obvious brand/generic equivalents collapse to goserelin.")
    register_normalization_entries(mapping, keys=["poly e"], normalized_name="poly e", rule="case_punctuation_normalization_only", notes="Exact trial-agent label retained with case/punctuation normalization only.")
    register_normalization_entries(mapping, keys=["metformin"], normalized_name="metformin", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["trastuzumab", "herceptin"], normalized_name="trastuzumab", rule="brand_generic_equivalence_or_case_normalization", notes="Trastuzumab and Herceptin collapse to trastuzumab.")
    register_normalization_entries(mapping, keys=["doxorubicin", "adriamycin"], normalized_name="doxorubicin", rule="brand_generic_equivalence_or_case_normalization", notes="Doxorubicin and Adriamycin collapse to doxorubicin.")
    register_normalization_entries(mapping, keys=["cyclophosphamide", "cytoxan"], normalized_name="cyclophosphamide", rule="brand_generic_equivalence_or_case_normalization", notes="Cyclophosphamide and Cytoxan collapse to cyclophosphamide.")
    register_normalization_entries(mapping, keys=["docetaxel", "taxotere"], normalized_name="docetaxel", rule="brand_generic_equivalence_or_case_normalization", notes="Docetaxel and Taxotere collapse to docetaxel.")
    register_normalization_entries(mapping, keys=["carboplatin"], normalized_name="carboplatin", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["cisplatin"], normalized_name="cisplatin", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["5 fluorouracil", "5 fu", "fluorouracil"], normalized_name="5-fluorouracil", rule="generic_synonym_and_case_normalization", notes="5-FU and fluorouracil labels collapse to 5-fluorouracil.")
    register_normalization_entries(mapping, keys=["tamoxifen", "nolvadex", "tamoxifen novadex", "tamoxifen citrate"], normalized_name="tamoxifen", rule="brand_generic_equivalence_or_case_normalization", notes="Tamoxifen brand/salt variants collapse to tamoxifen.")
    register_normalization_entries(mapping, keys=["letrozole", "femara", "letrozole femara", "femara letrozole"], normalized_name="letrozole", rule="brand_generic_equivalence_or_parenthetical_equivalence", notes="Letrozole and Femara labels collapse to letrozole.")
    register_normalization_entries(mapping, keys=["anastrozole", "arimidex", "anastrozole arimidex", "arimidex anastrozole"], normalized_name="anastrozole", rule="brand_generic_equivalence_or_parenthetical_equivalence", notes="Anastrozole and Arimidex labels collapse to anastrozole.")
    register_normalization_entries(mapping, keys=["exemestane", "aromasin", "exemestane aromasin", "aromasin exemestane"], normalized_name="exemestane", rule="brand_generic_equivalence_or_parenthetical_equivalence", notes="Exemestane and Aromasin labels collapse to exemestane.")
    register_normalization_entries(mapping, keys=["epirubicin"], normalized_name="epirubicin", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["gemcitabine", "gemzar"], normalized_name="gemcitabine", rule="brand_generic_equivalence_or_case_normalization", notes="Gemcitabine and Gemzar collapse to gemcitabine.")
    register_normalization_entries(mapping, keys=["capecitabine", "xeloda", "xeloda capecitabine"], normalized_name="capecitabine", rule="brand_generic_equivalence_or_parenthetical_equivalence", notes="Capecitabine and Xeloda collapse to capecitabine.")
    register_normalization_entries(mapping, keys=["e 75"], normalized_name="e-75", rule="case_punctuation_normalization_only", notes="Exact trial-agent label retained with punctuation normalization only.")
    register_normalization_entries(mapping, keys=["ae 37"], normalized_name="ae-37", rule="case_punctuation_normalization_only", notes="Exact trial-agent label retained with punctuation normalization only.")
    register_normalization_entries(mapping, keys=["albumin bound paclitaxel", "paclitaxel protein bound", "abraxane"], normalized_name="albumin-bound paclitaxel", rule="formulation_or_brand_generic_equivalence", notes="Protein-bound paclitaxel and Abraxane collapse to albumin-bound paclitaxel.")
    register_normalization_entries(mapping, keys=["paclitaxel", "taxol"], normalized_name="paclitaxel", rule="brand_generic_equivalence_or_case_normalization", notes="Paclitaxel and Taxol collapse to paclitaxel.")
    register_normalization_entries(mapping, keys=["fulvestrant", "faslodex", "fulvestrant faslodex"], normalized_name="fulvestrant", rule="brand_generic_equivalence_or_parenthetical_equivalence", notes="Fulvestrant and Faslodex collapse to fulvestrant.")
    register_normalization_entries(mapping, keys=["denosumab", "xgeva"], normalized_name="denosumab", rule="brand_generic_equivalence_or_case_normalization", notes="Denosumab and Xgeva collapse to denosumab.")
    register_normalization_entries(mapping, keys=["doxorubicin liposome", "doxil"], normalized_name="liposomal doxorubicin", rule="formulation_specific_or_brand_generic_equivalence", notes="Liposomal doxorubicin labels collapse to liposomal doxorubicin.")
    register_normalization_entries(mapping, keys=["lupron", "leuprolide acetate lupron", "leuprolide"], normalized_name="leuprolide", rule="brand_generic_equivalence_or_parenthetical_equivalence", notes="Lupron and explicit leuprolide labels collapse to leuprolide.")
    register_normalization_entries(mapping, keys=["lapatinib"], normalized_name="lapatinib", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["bevacizumab", "avastin"], normalized_name="bevacizumab", rule="brand_generic_equivalence_or_case_normalization", notes="Bevacizumab and Avastin collapse to bevacizumab.")
    register_normalization_entries(mapping, keys=["vinorelbine", "navelbine"], normalized_name="vinorelbine", rule="brand_generic_equivalence_or_case_normalization", notes="Vinorelbine and Navelbine collapse to vinorelbine.")
    register_normalization_entries(mapping, keys=["methotrexate"], normalized_name="methotrexate", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["ibandronate"], normalized_name="ibandronate", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["pamidronate"], normalized_name="pamidronate", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["palonosetron", "aloxi"], normalized_name="palonosetron", rule="brand_generic_equivalence_or_case_normalization", notes="Aloxi and palonosetron collapse to palonosetron.")
    register_normalization_entries(mapping, keys=["tesetaxel"], normalized_name="tesetaxel", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["pemetrexed"], normalized_name="pemetrexed", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["rituximab"], normalized_name="rituximab", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["prednisone"], normalized_name="prednisone", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["triptorelin"], normalized_name="triptorelin", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["ifosfamide"], normalized_name="ifosfamide", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["mesna 1", "mesna 2"], normalized_name="mesna", rule="numbered_supportive_drug_label_normalization", notes="Number-suffixed Mesna labels collapse to mesna.")
    register_normalization_entries(mapping, keys=["vp 16"], normalized_name="etoposide", rule="alias_equivalence", notes="VP-16 collapses to etoposide.")
    register_normalization_entries(mapping, keys=["clodronate"], normalized_name="clodronate", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["fareston"], normalized_name="toremifene", rule="brand_generic_equivalence", notes="Fareston collapses to toremifene.")
    register_normalization_entries(mapping, keys=["everolimus"], normalized_name="everolimus", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["mitomycin"], normalized_name="mitomycin", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["mitoxantrone"], normalized_name="mitoxantrone", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["pamidronic acid"], normalized_name="pamidronate", rule="acid_salt_form_equivalence", notes="Pamidronic acid is normalized to pamidronate.")
    register_normalization_entries(mapping, keys=["leuprorelin"], normalized_name="leuprolide", rule="regional_name_equivalence", notes="Leuprorelin is normalized to leuprolide.")
    register_normalization_entries(mapping, keys=["vinblastine"], normalized_name="vinblastine", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["vincristine"], normalized_name="vincristine", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["megace"], normalized_name="megestrol acetate", rule="brand_generic_equivalence", notes="Megace collapses to megestrol acetate.")
    register_normalization_entries(mapping, keys=["yondelis"], normalized_name="trabectedin", rule="brand_generic_equivalence", notes="Yondelis collapses to trabectedin.")
    register_normalization_entries(mapping, keys=["ixabepilone"], normalized_name="ixabepilone", rule="case_punctuation_normalization_only", notes="Single-agent name retained with case normalization only.")
    register_normalization_entries(mapping, keys=["neulasta"], normalized_name="pegfilgrastim", rule="brand_generic_equivalence", notes="Neulasta collapses to pegfilgrastim.")
    register_normalization_entries(mapping, keys=["zoledronic acid", "zometa"], normalized_name="zoledronic acid", rule="brand_generic_equivalence_or_case_normalization", notes="Zoledronic acid and Zometa collapse to zoledronic acid.")
    return mapping


def build_review_normalization_map() -> dict[str, tuple[str, str, str]]:
    mapping: dict[str, tuple[str, str, str]] = {}
    register_normalization_entries(mapping, keys=["letrozol", "letrozolum"], normalized_name="letrozole", rule="obvious_spelling_variant_review_needed", notes="Letrozole spelling variants collapse provisionally to letrozole and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["anastrazole", "anastrozolum"], normalized_name="anastrozole", rule="obvious_spelling_variant_review_needed", notes="Anastrozole spelling variants collapse provisionally to anastrozole and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["capecetabine"], normalized_name="capecitabine", rule="obvious_spelling_variant_review_needed", notes="Capecitabine spelling variants collapse provisionally to capecitabine and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["cyclophasphamide", "cyclophospamide", "cyclophosphane", "cyclophosphamid", "cyclophosphamidum", "cyotxan", "cytoxen"], normalized_name="cyclophosphamide", rule="obvious_spelling_variant_review_needed", notes="Cyclophosphamide spelling variants collapse provisionally to cyclophosphamide and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["trustuzumab"], normalized_name="trastuzumab", rule="obvious_spelling_variant_review_needed", notes="Trastuzumab spelling variants collapse provisionally to trastuzumab and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["adriamyicin", "adrimycin", "adriamicin", "adrimicin", "doxorubicine", "doxorubicinum"], normalized_name="doxorubicin", rule="brand_or_generic_spelling_variant_review_needed", notes="Doxorubicin/Adriamycin spelling variants collapse provisionally to doxorubicin and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["epirubicoin"], normalized_name="epirubicin", rule="obvious_spelling_variant_review_needed", notes="Epirubicin spelling variants collapse provisionally to epirubicin and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["flourouracil", "5 flourouracil"], normalized_name="5-fluorouracil", rule="obvious_spelling_variant_review_needed", notes="Fluorouracil spelling variants collapse provisionally to 5-fluorouracil and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["metotreksat"], normalized_name="methotrexate", rule="obvious_spelling_variant_review_needed", notes="Methotrexate spelling variants collapse provisionally to methotrexate and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["paclitaxelum"], normalized_name="paclitaxel", rule="obvious_spelling_variant_review_needed", notes="Paclitaxel spelling variants collapse provisionally to paclitaxel and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["doxetaxel"], normalized_name="docetaxel", rule="obvious_spelling_variant_review_needed", notes="Docetaxel spelling variants collapse provisionally to docetaxel and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["clodronic acid"], normalized_name="clodronate", rule="acid_salt_or_spelling_variant_review_needed", notes="Clodronate spelling/acid-form variants collapse provisionally to clodronate and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["tamoxiphene", "tamoxiphen"], normalized_name="tamoxifen", rule="obvious_spelling_variant_review_needed", notes="Tamoxifen spelling variants collapse provisionally to tamoxifen and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["aromatase exemestane"], normalized_name="exemestane", rule="qualified_name_review_needed", notes="Qualified exemestane labels collapse provisionally to exemestane and remain flagged for manual review.")
    register_normalization_entries(mapping, keys=["doxorubicin hcl"], normalized_name="doxorubicin", rule="salt_form_review_needed", notes="Doxorubicin salt-form labels collapse provisionally to doxorubicin and remain flagged for manual review.")
    return mapping


def build_confident_class_map() -> dict[str, tuple[str, str, str]]:
    mapping: dict[str, tuple[str, str, str]] = {}
    register_class_entries(mapping, names=["cyclophosphamide", "ifosfamide"], provisional_drug_class="alkylating_agent", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to alkylating-agent class in v1.")
    register_class_entries(mapping, names=["doxorubicin", "epirubicin", "liposomal doxorubicin"], provisional_drug_class="anthracycline", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to anthracycline class in v1.")
    register_class_entries(mapping, names=["paclitaxel", "albumin-bound paclitaxel", "docetaxel", "tesetaxel"], provisional_drug_class="taxane", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to taxane class in v1.")
    register_class_entries(mapping, names=["carboplatin", "cisplatin"], provisional_drug_class="platinum", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to platinum class in v1.")
    register_class_entries(mapping, names=["5-fluorouracil", "capecitabine", "gemcitabine", "methotrexate", "pemetrexed"], provisional_drug_class="antimetabolite", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to antimetabolite class in v1.")
    register_class_entries(mapping, names=["tamoxifen", "toremifene"], provisional_drug_class="endocrine_serm", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to endocrine SERM class in v1.")
    register_class_entries(mapping, names=["anastrozole", "letrozole", "exemestane"], provisional_drug_class="endocrine_aromatase_inhibitor", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to aromatase-inhibitor class in v1.")
    register_class_entries(mapping, names=["fulvestrant", "goserelin", "leuprolide", "triptorelin", "megestrol acetate"], provisional_drug_class="endocrine_other", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to endocrine-other class in v1.")
    register_class_entries(mapping, names=["trastuzumab", "lapatinib"], provisional_drug_class="her2_targeted", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to HER2-targeted class in v1.")
    register_class_entries(mapping, names=["denosumab", "zoledronic acid", "ibandronate", "pamidronate", "clodronate", "mesna", "palonosetron", "pegfilgrastim", "prednisone"], provisional_drug_class="ancillary_supportive", rule="direct_supportive_agent_class_lookup", notes="These names are mapped directly to ancillary/supportive class in v1.")
    register_class_entries(mapping, names=["vinorelbine", "vinblastine", "vincristine", "mitomycin", "mitoxantrone", "ixabepilone", "etoposide", "trabectedin"], provisional_drug_class="other_cytotoxic", rule="direct_single_agent_class_lookup", notes="These names are mapped directly to other-cytotoxic class in v1.")
    return mapping


def build_review_class_map() -> dict[str, tuple[str, str, str]]:
    mapping: dict[str, tuple[str, str, str]] = {}
    register_class_entries(mapping, names=["poly e", "e-75", "ae-37", "rituximab"], provisional_drug_class="immunotherapy", rule="provisional_review_needed_class_lookup", notes="These names are provisionally treated as immunotherapy but remain flagged for manual review in v1.")
    register_class_entries(mapping, names=["bevacizumab", "everolimus", "metformin"], provisional_drug_class="unknown_or_review_needed", rule="outside_current_class_vocabulary_review_needed", notes="These normalized names are preserved, but v1 does not force them into a narrower provisional class.")
    return mapping


CONFIDENT_NORMALIZATION_MAP = build_confident_normalization_map()
REVIEW_NORMALIZATION_MAP = build_review_normalization_map()
CONFIDENT_CLASS_MAP = build_confident_class_map()
REVIEW_CLASS_MAP = build_review_class_map()


def raw_value_contains_compound_marker(raw_value: str) -> bool:
    lowered = raw_value.lower()
    return "+" in raw_value or "/" in raw_value or "," in raw_value or " and " in lowered or " or " in lowered


def determine_normalization_decision(raw_drug_name: str) -> NormalizationDecision:
    match_key = normalize_raw_drug_name_for_matching(raw_drug_name)
    is_missing_like = is_missing_like_match_key(match_key)
    if is_missing_like:
        return NormalizationDecision(raw_drug_name, match_key, raw_drug_name, "missing_like_value_preserved", "unmapped", "yes", "unmapped", "Missing-like raw values are preserved exactly and left unmapped in v1.", True, False)
    if match_key in REGIMEN_OR_PLACEHOLDER_MATCH_KEYS:
        return NormalizationDecision(raw_drug_name, match_key, raw_drug_name, "regimen_or_placeholder_label_preserved", "review_needed", "yes", "mapped_but_review", "Regimen-like or placeholder labels are preserved exactly and are not forced into a single normalized drug name in v1.", False, True)
    if "placebo" in match_key:
        return NormalizationDecision(raw_drug_name, match_key, raw_drug_name, "trial_or_placebo_language_preserved", "review_needed", "yes", "mapped_but_review", "Trial/placebo language is preserved exactly and remains for manual review in v1.", False, True)
    if match_key in SPACE_ONLY_COMPOUND_MATCH_KEYS or raw_value_contains_compound_marker(raw_drug_name):
        return NormalizationDecision(raw_drug_name, match_key, raw_drug_name, "compound_raw_name_preserved_without_split", "review_needed", "yes", "mapped_but_review", "Compound or regimen-like raw drug strings are preserved exactly and not split into component agents in v1.", False, True)
    confident_entry = CONFIDENT_NORMALIZATION_MAP.get(match_key)
    if confident_entry is not None:
        normalized_name, rule, notes = confident_entry
        return NormalizationDecision(raw_drug_name, match_key, normalized_name, rule, "confident", "no", "mapped_confidently", notes, False, False)
    review_entry = REVIEW_NORMALIZATION_MAP.get(match_key)
    if review_entry is not None:
        normalized_name, rule, notes = review_entry
        return NormalizationDecision(raw_drug_name, match_key, normalized_name, rule, "review_needed", "yes", "mapped_but_review", notes, False, False)
    return NormalizationDecision(raw_drug_name, match_key, raw_drug_name, "unrecognized_single_agent_label_preserved", "unmapped", "yes", "unmapped", "Raw value is preserved exactly because v1 does not assign a stronger normalized name without explicit audited support.", False, False)


def determine_class_decision(provisional_normalized_drug_name: str) -> ClassDecision:
    confident_entry = CONFIDENT_CLASS_MAP.get(provisional_normalized_drug_name)
    if confident_entry is not None:
        provisional_drug_class, rule, notes = confident_entry
        return ClassDecision(provisional_normalized_drug_name, provisional_drug_class, rule, "confident", "no", notes)
    review_entry = REVIEW_CLASS_MAP.get(provisional_normalized_drug_name)
    if review_entry is not None:
        provisional_drug_class, rule, notes = review_entry
        return ClassDecision(provisional_normalized_drug_name, provisional_drug_class, rule, "review_needed", "yes", notes)
    return ClassDecision(provisional_normalized_drug_name, "unknown_or_review_needed", "no_supported_v1_class_assignment", "review_needed", "yes", "V1 preserves this normalized name but does not force it into a narrower provisional class.")


def build_workflow_paths(helper_module: Any) -> WorkflowPaths:
    repo_root = helper_module.detect_repo_root(Path(__file__).resolve().parent)
    trial_root = repo_root / "09-trials" / "01-tcga-only-source-audited"
    trial_config = trial_root / "04-config" / "trial_config.yaml"
    if not trial_config.exists():
        raise DrugNameNormalizationV1Error(f"Required trial config not found: {trial_config}")
    trial_config_data = helper_module.load_yaml(trial_config)
    processed_root = repo_root / str(trial_config_data.get("processed_data_root", "01-data/processed"))
    audit_root = repo_root / str(trial_config_data.get("audit_root", "01-data/audit"))
    results_root = repo_root / str(
        trial_config_data.get("results_root", "09-trials/01-tcga-only-source-audited/05-results")
    )
    treatment_prep_root = audit_root / "tcga-brca" / "treatment-prep"
    variables_root = audit_root / "tcga-brca" / "variables"
    return WorkflowPaths(
        repo_root=repo_root,
        trial_config=trial_config,
        results_root=results_root,
        processed_runs_root=processed_root / "tcga-brca" / "treatment-prep" / "drug_name_normalization_v1_runs",
        audit_runs_root=treatment_prep_root / "drug_name_normalization_v1_runs",
        latest_pointer=treatment_prep_root / "tcga_brca_drug_name_normalization_v1_latest.json",
        patient_treatment_profile_latest_pointer=treatment_prep_root / "tcga_brca_patient_treatment_profile_v1_latest.json",
        patient_treatment_grouping_latest_pointer=treatment_prep_root / "tcga_brca_patient_treatment_grouping_v1_latest.json",
        clinical_biotabs_latest_pointer=variables_root / "tcga_brca_clinical_biotabs_latest.json",
    )


def load_workflow_inputs(paths: WorkflowPaths, helper_module: Any) -> WorkflowInputs:
    required_pointers = [
        (paths.patient_treatment_profile_latest_pointer, "patient treatment profile v1 latest pointer"),
        (paths.patient_treatment_grouping_latest_pointer, "patient treatment grouping v1 latest pointer"),
        (paths.clinical_biotabs_latest_pointer, "clinical biotabs latest pointer"),
    ]
    for pointer_path, label in required_pointers:
        if not pointer_path.exists():
            raise DrugNameNormalizationV1Error(f"Required {label} not found: {pointer_path}")

    profile_pointer = helper_module.load_json(paths.patient_treatment_profile_latest_pointer)
    helper_module.require_keys(
        profile_pointer,
        {
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "clinical_biotab_parse_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "patient_treatment_profile_v1_tsv",
            "run_log_json",
        },
        "patient treatment profile v1 latest pointer",
        paths.patient_treatment_profile_latest_pointer,
    )
    grouping_pointer = helper_module.load_json(paths.patient_treatment_grouping_latest_pointer)
    helper_module.require_keys(
        grouping_pointer,
        {
            "patient_treatment_grouping_v1_run_id",
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "patient_treatment_grouping_v1_tsv",
            "run_log_json",
        },
        "patient treatment grouping v1 latest pointer",
        paths.patient_treatment_grouping_latest_pointer,
    )
    biotabs_pointer = helper_module.load_json(paths.clinical_biotabs_latest_pointer)
    helper_module.require_keys(
        biotabs_pointer,
        {"parse_run_id", "source_run_id", "processed_run_directory", "run_log_json"},
        "clinical biotabs latest pointer",
        paths.clinical_biotabs_latest_pointer,
    )

    profile_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(profile_pointer["run_log_json"]),
        label="patient treatment profile v1 run log",
        helper_module=helper_module,
    )
    grouping_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(grouping_pointer["run_log_json"]),
        label="patient treatment grouping v1 run log",
        helper_module=helper_module,
    )
    biotabs_run_log = require_completed_run_log(
        repo_root=paths.repo_root,
        run_log_relative_path=str(biotabs_pointer["run_log_json"]),
        label="clinical biotabs run log",
        helper_module=helper_module,
    )

    for left_value, right_value, label in [
        (
            str(profile_pointer["patient_treatment_profile_v1_run_id"]),
            str(grouping_pointer["patient_treatment_profile_v1_run_id"]),
            "patient_treatment_profile_v1_run_id",
        ),
        (
            str(profile_pointer["treatment_os_overlap_v1_run_id"]),
            str(grouping_pointer["treatment_os_overlap_v1_run_id"]),
            "treatment_os_overlap_v1_run_id",
        ),
        (
            str(profile_pointer["os_endpoint_v1_run_id"]),
            str(grouping_pointer["os_endpoint_v1_run_id"]),
            "os_endpoint_v1_run_id",
        ),
        (
            str(profile_pointer["baseline_model_input_v1_run_id"]),
            str(grouping_pointer["baseline_model_input_v1_run_id"]),
            "baseline_model_input_v1_run_id",
        ),
        (
            str(profile_pointer["cohort_v1_build_id"]),
            str(grouping_pointer["cohort_v1_build_id"]),
            "cohort_v1_build_id",
        ),
        (
            str(profile_pointer["clinical_biotab_parse_run_id"]),
            str(biotabs_pointer["parse_run_id"]),
            "clinical_biotab_parse_run_id",
        ),
    ]:
        if left_value != right_value:
            raise DrugNameNormalizationV1Error(f"Mismatch for {label} across upstream pointers.")

    profile_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(profile_pointer["patient_treatment_profile_v1_tsv"]),
        "patient_treatment_profile_v1.tsv",
    )
    grouping_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        str(grouping_pointer["patient_treatment_grouping_v1_tsv"]),
        "patient_treatment_grouping_v1.tsv",
    )
    clinical_drug_tsv = helper_module.resolve_existing_path(
        paths.repo_root,
        (Path(str(biotabs_pointer["processed_run_directory"])) / "clinical_drug.tsv").as_posix(),
        "clinical_drug.tsv",
    )

    profile_rows = helper_module.read_tsv_dict_rows(profile_tsv)
    grouping_rows = helper_module.read_tsv_dict_rows(grouping_tsv)
    clinical_drug_rows = helper_module.read_tsv_dict_rows(clinical_drug_tsv)

    require_columns(
        profile_rows,
        {
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "clinical_biotab_parse_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
            "has_any_drug_row",
            "drug_row_count",
            "drug_therapy_type_values_json",
            "regimen_context_values_json",
            "drug_name_values_json",
        },
        "patient_treatment_profile_v1.tsv",
    )
    require_columns(
        grouping_rows,
        {
            "patient_treatment_grouping_v1_run_id",
            "patient_treatment_profile_v1_run_id",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
            "bcr_patient_barcode",
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
            "has_any_drug_row",
        },
        "patient_treatment_grouping_v1.tsv",
    )
    require_columns(
        clinical_drug_rows,
        {
            "bcr_patient_uuid",
            "bcr_patient_barcode",
            "pharmaceutical_therapy_drug_name",
            "pharmaceutical_therapy_type",
            "therapy_regimen",
            "pharm_regimen",
        },
        "clinical_drug.tsv",
    )

    return WorkflowInputs(
        patient_treatment_profile_latest_pointer=profile_pointer,
        patient_treatment_profile_run_log=profile_run_log,
        patient_treatment_grouping_latest_pointer=grouping_pointer,
        patient_treatment_grouping_run_log=grouping_run_log,
        clinical_biotabs_latest_pointer=biotabs_pointer,
        clinical_biotabs_run_log=biotabs_run_log,
        patient_treatment_profile_rows=profile_rows,
        patient_treatment_grouping_rows=grouping_rows,
        clinical_drug_rows=clinical_drug_rows,
        input_paths={
            "patient_treatment_profile_v1_tsv": profile_tsv,
            "patient_treatment_grouping_v1_tsv": grouping_tsv,
            "clinical_drug_tsv": clinical_drug_tsv,
        },
    )


def build_inventory_and_normalization_rows(
    run_id: str,
    clinical_drug_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, NormalizationDecision]]:
    inventory_by_raw: dict[str, dict[str, Any]] = {}
    raw_name_order: list[str] = []

    for row in clinical_drug_rows:
        raw_drug_name = str(row.get("pharmaceutical_therapy_drug_name", "")).strip()
        if raw_drug_name not in inventory_by_raw:
            raw_name_order.append(raw_drug_name)
            inventory_by_raw[raw_drug_name] = {
                "drug_row_count": 0,
                "patient_barcodes": [],
                "patient_barcode_seen": set(),
                "therapy_type_values": [],
                "therapy_type_seen": set(),
                "regimen_context_values": [],
                "regimen_context_seen": set(),
            }
        bucket = inventory_by_raw[raw_drug_name]
        bucket["drug_row_count"] += 1

        barcode = normalize_barcode(str(row.get("bcr_patient_barcode", "")))
        if barcode and barcode not in bucket["patient_barcode_seen"]:
            bucket["patient_barcode_seen"].add(barcode)
            bucket["patient_barcodes"].append(barcode)

        therapy_type = str(row.get("pharmaceutical_therapy_type", "")).strip()
        if (
            therapy_type
            and not is_missing_like_match_key(normalize_raw_drug_name_for_matching(therapy_type))
            and therapy_type not in bucket["therapy_type_seen"]
        ):
            bucket["therapy_type_seen"].add(therapy_type)
            bucket["therapy_type_values"].append(therapy_type)

        for field_name in ["therapy_regimen", "pharm_regimen"]:
            regimen_value = str(row.get(field_name, "")).strip()
            if (
                regimen_value
                and not is_missing_like_match_key(normalize_raw_drug_name_for_matching(regimen_value))
                and regimen_value not in bucket["regimen_context_seen"]
            ):
                bucket["regimen_context_seen"].add(regimen_value)
                bucket["regimen_context_values"].append(regimen_value)

    sorted_raw_names = sorted(
        raw_name_order,
        key=lambda raw_name: (-int(inventory_by_raw[raw_name]["drug_row_count"]), raw_name),
    )
    inventory_rows: list[dict[str, str]] = []
    normalization_rows: list[dict[str, str]] = []
    normalization_decisions: dict[str, NormalizationDecision] = {}

    for raw_drug_name in sorted_raw_names:
        decision = determine_normalization_decision(raw_drug_name)
        normalization_decisions[raw_drug_name] = decision
        bucket = inventory_by_raw[raw_drug_name]
        inventory_rows.append(
            {
                "drug_name_normalization_v1_run_id": run_id,
                "raw_drug_name": raw_drug_name,
                "normalized_raw_drug_name_for_matching": decision.normalized_raw_drug_name_for_matching,
                "drug_row_count": str(bucket["drug_row_count"]),
                "distinct_patient_count": str(len(bucket["patient_barcodes"])),
                "therapy_type_values_json": json_list(bucket["therapy_type_values"]),
                "regimen_context_values_json": json_list(bucket["regimen_context_values"]),
                "example_patient_barcodes_json": json_list(bucket["patient_barcodes"][:10]),
                "mapping_status": decision.mapping_status,
                "notes": decision.notes,
            }
        )
        normalization_rows.append(
            {
                "drug_name_normalization_v1_run_id": run_id,
                "raw_drug_name": raw_drug_name,
                "normalized_raw_drug_name_for_matching": decision.normalized_raw_drug_name_for_matching,
                "provisional_normalized_drug_name": decision.provisional_normalized_drug_name,
                "normalization_rule": decision.normalization_rule,
                "normalization_confidence": decision.normalization_confidence,
                "manual_review_needed": decision.manual_review_needed,
                "notes": decision.notes,
            }
        )

    return inventory_rows, normalization_rows, normalization_decisions


def build_class_rows(
    run_id: str,
    normalization_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], dict[str, ClassDecision]]:
    normalized_names = ordered_unique([str(row["provisional_normalized_drug_name"]) for row in normalization_rows])
    class_rows: list[dict[str, str]] = []
    class_decisions: dict[str, ClassDecision] = {}

    for provisional_normalized_drug_name in sorted(normalized_names):
        decision = determine_class_decision(provisional_normalized_drug_name)
        class_decisions[provisional_normalized_drug_name] = decision
        class_rows.append(
            {
                "drug_name_normalization_v1_run_id": run_id,
                "provisional_normalized_drug_name": provisional_normalized_drug_name,
                "provisional_drug_class": decision.provisional_drug_class,
                "class_mapping_rule": decision.class_mapping_rule,
                "class_mapping_confidence": decision.class_mapping_confidence,
                "manual_review_needed": decision.manual_review_needed,
                "notes": decision.notes,
            }
        )

    return class_rows, class_decisions


def build_patient_profile_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
    normalization_decisions: dict[str, NormalizationDecision],
    class_decisions: dict[str, ClassDecision],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    profile_rows = workflow_inputs.patient_treatment_profile_rows
    grouping_rows = workflow_inputs.patient_treatment_grouping_rows
    clinical_drug_rows = workflow_inputs.clinical_drug_rows
    grouping_by_barcode = build_unique_lookup_by_barcode(
        grouping_rows,
        barcode_field="bcr_patient_barcode",
        label="patient_treatment_grouping_v1.tsv",
    )
    drug_rows_by_barcode = group_rows_by_barcode(clinical_drug_rows, barcode_field="bcr_patient_barcode")

    treated_profile_rows = [row for row in profile_rows if str(row["has_any_drug_row"]).strip() == "yes"]
    treated_profile_barcodes = [normalize_barcode(str(row["bcr_patient_barcode"])) for row in treated_profile_rows]
    clinical_drug_barcodes = sorted(drug_rows_by_barcode.keys())
    if sorted(treated_profile_barcodes) != clinical_drug_barcodes:
        raise DrugNameNormalizationV1Error(
            "Treated-patient set from patient_treatment_profile_v1.tsv does not match the patient set in clinical_drug.tsv."
        )

    patient_rows: list[dict[str, str]] = []
    metrics = {
        "patients_with_any_normalized_name": 0,
        "patients_with_any_unknown_or_review_needed_mapping": 0,
        "patients_with_name_normalization_manual_review": 0,
        "patients_with_class_mapping_manual_review": 0,
        "patients_with_any_unmapped_raw_name": 0,
        "patients_with_any_review_needed_name_normalization": 0,
        "patients_with_single_drug_class": 0,
        "patients_with_multiple_drug_classes": 0,
        "patients_with_unknown_only_drug_class": 0,
        "patients_with_missing_like_only_raw_drug_names": 0,
        "patient_profile_row_order_preserved": True,
    }

    for expected_index, profile_row in enumerate(treated_profile_rows):
        barcode = normalize_barcode(str(profile_row["bcr_patient_barcode"]))
        grouping_row = grouping_by_barcode.get(barcode)
        if grouping_row is None:
            raise DrugNameNormalizationV1Error(
                f"Missing patient_treatment_grouping_v1.tsv row for treated patient barcode {barcode}."
            )
        if str(grouping_row["has_any_drug_row"]).strip() != "yes":
            raise DrugNameNormalizationV1Error(
                f"Grouping row does not preserve has_any_drug_row == yes for treated barcode {barcode}."
            )
        for field_name in [
            "bcr_patient_uuid",
            "provisional_patient_row_id",
            "baseline_analysis_v1_row_id",
            "feature_set_v1_row_index",
            "treatment_os_overlap_v1_run_id",
            "os_endpoint_v1_run_id",
            "baseline_model_input_v1_run_id",
            "cohort_v1_build_id",
        ]:
            if str(profile_row[field_name]).strip() != str(grouping_row[field_name]).strip():
                raise DrugNameNormalizationV1Error(
                    f"Profile/grouping lineage mismatch for field {field_name} at barcode {barcode}."
                )
        if str(grouping_row["patient_treatment_profile_v1_run_id"]).strip() != str(
            profile_row["patient_treatment_profile_v1_run_id"]
        ).strip():
            raise DrugNameNormalizationV1Error(
                f"Grouping/profile run-id continuity mismatch at barcode {barcode}."
            )

        patient_drug_rows = drug_rows_by_barcode.get(barcode, [])
        if not patient_drug_rows:
            raise DrugNameNormalizationV1Error(f"Missing clinical_drug.tsv rows for treated barcode {barcode}.")
        if parse_int_or_none(str(profile_row["drug_row_count"])) != len(patient_drug_rows):
            raise DrugNameNormalizationV1Error(
                f"drug_row_count mismatch between profile and clinical_drug rows for barcode {barcode}."
            )

        raw_values_all = [str(row.get("pharmaceutical_therapy_drug_name", "")).strip() for row in patient_drug_rows]
        raw_values_distinct = ordered_unique(raw_values_all)
        non_missing_raw_values_all = [
            raw_value
            for raw_value in raw_values_all
            if not is_missing_like_match_key(normalize_raw_drug_name_for_matching(raw_value))
        ]
        non_missing_raw_values_distinct = ordered_non_missing_distinct_values(raw_values_all)
        profile_raw_values = parse_json_list(str(profile_row["drug_name_values_json"]), "drug_name_values_json")
        allowed_profile_raw_values = [non_missing_raw_values_distinct, non_missing_raw_values_all, raw_values_distinct]
        if profile_raw_values not in allowed_profile_raw_values:
            raise DrugNameNormalizationV1Error(
                "drug_name_values_json continuity mismatch between patient_treatment_profile_v1.tsv "
                f"and clinical_drug.tsv for barcode {barcode}."
            )

        normalized_values_all = [
            normalization_decisions[raw_value].provisional_normalized_drug_name for raw_value in raw_values_all
        ]
        normalized_values_distinct = ordered_unique(normalized_values_all)
        class_values_all = [
            class_decisions[normalized_value].provisional_drug_class for normalized_value in normalized_values_all
        ]
        class_values_distinct = ordered_unique(class_values_all)
        known_class_values_distinct = [
            class_value for class_value in class_values_distinct if class_value != "unknown_or_review_needed"
        ]
        unknown_class_present = "unknown_or_review_needed" in class_values_distinct

        if any(not normalization_decisions[raw_value].is_missing_like for raw_value in raw_values_all):
            metrics["patients_with_any_normalized_name"] += 1
        else:
            metrics["patients_with_missing_like_only_raw_drug_names"] += 1

        name_review_needed = any(
            normalization_decisions[raw_value].manual_review_needed == "yes" for raw_value in raw_values_all
        )
        class_review_needed = any(
            class_decisions[normalized_value].manual_review_needed == "yes"
            for normalized_value in normalized_values_all
        )
        any_unmapped_raw_name = any(
            normalization_decisions[raw_value].normalization_confidence == "unmapped"
            for raw_value in raw_values_all
        )
        any_review_needed_name = any(
            normalization_decisions[raw_value].normalization_confidence == "review_needed"
            for raw_value in raw_values_all
        )
        any_compound_or_regimen_like_raw_name = any(
            normalization_decisions[raw_value].is_compound_or_regimen_like for raw_value in raw_values_all
        )
        any_missing_like_raw_name = any(
            normalization_decisions[raw_value].is_missing_like for raw_value in raw_values_all
        )
        any_unknown_or_review = name_review_needed or class_review_needed

        if name_review_needed:
            metrics["patients_with_name_normalization_manual_review"] += 1
        if class_review_needed:
            metrics["patients_with_class_mapping_manual_review"] += 1
        if any_unmapped_raw_name:
            metrics["patients_with_any_unmapped_raw_name"] += 1
        if any_review_needed_name:
            metrics["patients_with_any_review_needed_name_normalization"] += 1
        if any_unknown_or_review:
            metrics["patients_with_any_unknown_or_review_needed_mapping"] += 1

        known_class_count = len(known_class_values_distinct)
        if known_class_count == 0:
            single_or_multi_drug_class = "unknown_only"
            metrics["patients_with_unknown_only_drug_class"] += 1
        elif known_class_count == 1:
            single_or_multi_drug_class = "single_class"
            metrics["patients_with_single_drug_class"] += 1
        else:
            single_or_multi_drug_class = "multi_class"
            metrics["patients_with_multiple_drug_classes"] += 1

        dominant_provisional_drug_class_if_any = (
            known_class_values_distinct[0] if known_class_count == 1 and not unknown_class_present else ""
        )

        flags: list[str] = []
        if any_unmapped_raw_name:
            flags.append("contains_unmapped_raw_drug_name")
        if any_review_needed_name:
            flags.append("contains_review_needed_name_normalization")
        if any_unknown_or_review:
            flags.append("contains_unknown_or_review_needed_class")
        if any_missing_like_raw_name:
            flags.append("contains_missing_like_raw_drug_name")
        if any_compound_or_regimen_like_raw_name:
            flags.append("contains_compound_or_regimen_like_raw_drug_name")
        if known_class_count == 1 and unknown_class_present:
            flags.append("contains_single_known_class_plus_unknown")
        if known_class_count > 1:
            flags.append("contains_multiple_known_drug_classes")
        ordered_flags = [flag for flag in PATIENT_PROFILE_FLAG_ORDER if flag in flags]

        patient_rows.append(
            {
                "drug_name_normalization_v1_run_id": run_id,
                "patient_treatment_profile_v1_run_id": str(profile_row["patient_treatment_profile_v1_run_id"]),
                "patient_treatment_grouping_v1_run_id": str(grouping_row["patient_treatment_grouping_v1_run_id"]),
                "treatment_os_overlap_v1_run_id": str(profile_row["treatment_os_overlap_v1_run_id"]),
                "os_endpoint_v1_run_id": str(profile_row["os_endpoint_v1_run_id"]),
                "clinical_biotab_parse_run_id": str(profile_row["clinical_biotab_parse_run_id"]),
                "clinical_source_run_id": str(workflow_inputs.clinical_biotabs_latest_pointer["source_run_id"]),
                "baseline_model_input_v1_run_id": str(profile_row["baseline_model_input_v1_run_id"]),
                "cohort_v1_build_id": str(profile_row["cohort_v1_build_id"]),
                "bcr_patient_barcode": str(profile_row["bcr_patient_barcode"]),
                "bcr_patient_uuid": str(profile_row["bcr_patient_uuid"]),
                "provisional_patient_row_id": str(profile_row["provisional_patient_row_id"]),
                "baseline_analysis_v1_row_id": str(profile_row["baseline_analysis_v1_row_id"]),
                "feature_set_v1_row_index": str(profile_row["feature_set_v1_row_index"]),
                "drug_row_count": str(profile_row["drug_row_count"]),
                "drug_therapy_type_values_json": str(profile_row["drug_therapy_type_values_json"]),
                "regimen_context_values_json": str(profile_row["regimen_context_values_json"]),
                "raw_drug_name_values_json": json_list(raw_values_distinct),
                "normalized_drug_name_values_json": json_list(normalized_values_distinct),
                "provisional_drug_class_values_json": json_list(class_values_distinct),
                "distinct_raw_drug_name_count": str(len(raw_values_distinct)),
                "distinct_normalized_drug_name_count": str(len(normalized_values_distinct)),
                "distinct_provisional_drug_class_count": str(len(class_values_distinct)),
                "single_or_multi_drug_class": single_or_multi_drug_class,
                "dominant_provisional_drug_class_if_any": dominant_provisional_drug_class_if_any,
                "drug_name_normalization_requires_manual_review": yes_no(name_review_needed),
                "drug_class_mapping_requires_manual_review": yes_no(class_review_needed),
                "drug_name_profile_flags_json": json_list(ordered_flags),
            }
        )

        if normalize_barcode(patient_rows[expected_index]["bcr_patient_barcode"]) != barcode:
            metrics["patient_profile_row_order_preserved"] = False

    metrics["treated_patient_count"] = len(treated_profile_rows)
    return patient_rows, metrics


def build_summary_rows(
    run_id: str,
    workflow_inputs: WorkflowInputs,
    inventory_rows: list[dict[str, str]],
    normalization_rows: list[dict[str, str]],
    class_rows: list[dict[str, str]],
    patient_rows: list[dict[str, str]],
    patient_metrics: dict[str, Any],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    inventory_counts = Counter(str(row["mapping_status"]) for row in inventory_rows)

    def sum_inventory_rows(metric_name: str) -> int:
        return sum(int(row["drug_row_count"]) for row in inventory_rows if str(row["mapping_status"]) == metric_name)

    distinct_provisional_classes = len(ordered_unique([str(row["provisional_drug_class"]) for row in class_rows]))
    if not inventory_rows or (
        inventory_counts.get("mapped_confidently", 0) == 0 and inventory_counts.get("mapped_but_review", 0) == 0
    ):
        readiness = READINESS_NOT_READY
        treatment_arm_freeze_review_status = TREATMENT_ARM_FREEZE_REVIEW_BLOCKED
    elif inventory_counts.get("unmapped", 0) > 0 or patient_metrics["patients_with_any_unknown_or_review_needed_mapping"] > 0:
        readiness = READINESS_DRUG_CLASS_REVIEW
        treatment_arm_freeze_review_status = TREATMENT_ARM_FREEZE_REVIEW_BLOCKED
    else:
        readiness = READINESS_TREATMENT_ARM_FREEZE_REVIEW_NEXT
        treatment_arm_freeze_review_status = TREATMENT_ARM_FREEZE_REVIEW_PROVISIONAL_NEXT

    summary_rows = [
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "inputs", "summary_metric": "patient_treatment_profile_v1_run_id", "summary_value": str(workflow_inputs.patient_treatment_profile_latest_pointer["patient_treatment_profile_v1_run_id"]), "notes": ""},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "inputs", "summary_metric": "patient_treatment_grouping_v1_run_id", "summary_value": str(workflow_inputs.patient_treatment_grouping_latest_pointer["patient_treatment_grouping_v1_run_id"]), "notes": ""},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "inputs", "summary_metric": "clinical_biotab_parse_run_id", "summary_value": str(workflow_inputs.clinical_biotabs_latest_pointer["parse_run_id"]), "notes": ""},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "inputs", "summary_metric": "clinical_source_run_id", "summary_value": str(workflow_inputs.clinical_biotabs_latest_pointer["source_run_id"]), "notes": ""},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_clinical_drug_row_count", "summary_value": str(len(workflow_inputs.clinical_drug_rows)), "notes": "Total drug-detail rows reviewed from saved clinical_drug.tsv."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_treated_patient_count", "summary_value": str(patient_metrics["treated_patient_count"]), "notes": "Patients with has_any_drug_row == yes carried forward from patient_treatment_profile_v1.tsv."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_distinct_raw_drug_names", "summary_value": str(len(inventory_rows)), "notes": "One row per exact raw pharmaceutical_therapy_drug_name value."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_mapped_confidently", "summary_value": str(inventory_counts.get("mapped_confidently", 0)), "notes": "Distinct raw names with confident v1 normalization."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_mapped_but_review", "summary_value": str(inventory_counts.get("mapped_but_review", 0)), "notes": "Distinct raw names preserved or normalized provisionally but still requiring review."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_unmapped", "summary_value": str(inventory_counts.get("unmapped", 0)), "notes": "Distinct raw names left unmapped in v1."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "drug_rows_mapped_confidently", "summary_value": str(sum_inventory_rows("mapped_confidently")), "notes": "Clinical drug rows covered by confident v1 normalization."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "drug_rows_mapped_but_review", "summary_value": str(sum_inventory_rows("mapped_but_review")), "notes": "Clinical drug rows tied to review-needed normalization decisions."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "drug_rows_unmapped", "summary_value": str(sum_inventory_rows("unmapped")), "notes": "Clinical drug rows tied to unmapped raw-name decisions."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_distinct_normalized_drug_names", "summary_value": str(len(class_rows)), "notes": "Distinct provisional normalized names after v1 collapsing and preserved-review labels."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "coverage", "summary_metric": "total_distinct_provisional_drug_classes", "summary_value": str(distinct_provisional_classes), "notes": "Distinct provisional drug-class values observed in the v1 class map."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_any_normalized_name", "summary_value": str(patient_metrics["patients_with_any_normalized_name"]), "notes": "Patients with at least one non-missing raw drug label available for name normalization."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_any_unknown_or_review_needed_mapping", "summary_value": str(patient_metrics["patients_with_any_unknown_or_review_needed_mapping"]), "notes": "Patients carrying any review-needed name or class mapping burden."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_any_unmapped_raw_name", "summary_value": str(patient_metrics["patients_with_any_unmapped_raw_name"]), "notes": "Patients carrying at least one raw name left unmapped in v1."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_any_review_needed_name_normalization", "summary_value": str(patient_metrics["patients_with_any_review_needed_name_normalization"]), "notes": "Patients carrying review-needed provisional name normalization decisions."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_name_normalization_manual_review", "summary_value": str(patient_metrics["patients_with_name_normalization_manual_review"]), "notes": "Patients requiring name-level manual review in v1."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_class_mapping_manual_review", "summary_value": str(patient_metrics["patients_with_class_mapping_manual_review"]), "notes": "Patients requiring class-level manual review in v1."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_single_drug_class", "summary_value": str(patient_metrics["patients_with_single_drug_class"]), "notes": "Patient count based on known classes only."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_multiple_drug_classes", "summary_value": str(patient_metrics["patients_with_multiple_drug_classes"]), "notes": "Patient count based on more than one known provisional class."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_unknown_only_drug_class", "summary_value": str(patient_metrics["patients_with_unknown_only_drug_class"]), "notes": "Patients with no known provisional class after v1 mapping."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "patients", "summary_metric": "total_patients_with_missing_like_only_raw_drug_names", "summary_value": str(patient_metrics["patients_with_missing_like_only_raw_drug_names"]), "notes": "Treated patients whose raw drug-name values are missing-like only."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "readiness", "summary_metric": "provisional_readiness_interpretation", "summary_value": readiness, "notes": "Treatment-arm freeze remains provisional until the remaining manual-review burden is acceptable."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "readiness", "summary_metric": "treatment_arm_freeze_review_status", "summary_value": treatment_arm_freeze_review_status, "notes": "Drug-name normalization v1 supports review only; it does not itself freeze treatment arms."},
        {"drug_name_normalization_v1_run_id": run_id, "summary_section": "readiness", "summary_metric": "treatment_recommendation_modeling_status", "summary_value": TREATMENT_RECOMMENDATION_MODELING_STATUS, "notes": "Treatment recommendation modeling remains out of scope."},
    ]

    return summary_rows, {
        "inventory_counts": inventory_counts,
        "distinct_provisional_classes": distinct_provisional_classes,
        "readiness": readiness,
        "treatment_arm_freeze_review_status": treatment_arm_freeze_review_status,
        "drug_rows_mapped_confidently": sum_inventory_rows("mapped_confidently"),
        "drug_rows_mapped_but_review": sum_inventory_rows("mapped_but_review"),
        "drug_rows_unmapped": sum_inventory_rows("unmapped"),
    }


def validate_outputs(
    *,
    workflow_inputs: WorkflowInputs,
    inventory_rows: list[dict[str, str]],
    normalization_rows: list[dict[str, str]],
    class_rows: list[dict[str, str]],
    patient_rows: list[dict[str, str]],
    summary_rows: list[dict[str, str]],
    patient_metrics: dict[str, Any],
    summary_metrics: dict[str, Any],
) -> dict[str, Any]:
    treated_profile_rows = [row for row in workflow_inputs.patient_treatment_profile_rows if str(row["has_any_drug_row"]).strip() == "yes"]
    treated_barcodes = [normalize_barcode(str(row["bcr_patient_barcode"])) for row in treated_profile_rows]
    clinical_drug_barcodes = sorted(
        ordered_unique([normalize_barcode(str(row["bcr_patient_barcode"])) for row in workflow_inputs.clinical_drug_rows])
    )
    inventory_raw_names = [str(row["raw_drug_name"]) for row in inventory_rows]
    normalization_raw_names = [str(row["raw_drug_name"]) for row in normalization_rows]
    class_normalized_names = [str(row["provisional_normalized_drug_name"]) for row in class_rows]
    normalized_names_from_map = ordered_unique([str(row["provisional_normalized_drug_name"]) for row in normalization_rows])
    patient_barcodes = [normalize_barcode(str(row["bcr_patient_barcode"])) for row in patient_rows]
    summary_lookup = {str(row["summary_metric"]): str(row["summary_value"]) for row in summary_rows}

    def summary_int(metric_name: str) -> int | None:
        return parse_int_or_none(summary_lookup.get(metric_name, ""))

    summary_counts_reconcile = (
        summary_int("total_clinical_drug_row_count") == len(workflow_inputs.clinical_drug_rows)
        and summary_int("total_treated_patient_count") == len(treated_profile_rows)
        and summary_int("total_distinct_raw_drug_names") == len(inventory_rows)
        and summary_int("total_mapped_confidently") == summary_metrics["inventory_counts"].get("mapped_confidently", 0)
        and summary_int("total_mapped_but_review") == summary_metrics["inventory_counts"].get("mapped_but_review", 0)
        and summary_int("total_unmapped") == summary_metrics["inventory_counts"].get("unmapped", 0)
        and summary_int("drug_rows_mapped_confidently") == summary_metrics["drug_rows_mapped_confidently"]
        and summary_int("drug_rows_mapped_but_review") == summary_metrics["drug_rows_mapped_but_review"]
        and summary_int("drug_rows_unmapped") == summary_metrics["drug_rows_unmapped"]
        and summary_int("total_distinct_normalized_drug_names") == len(class_rows)
        and summary_int("total_distinct_provisional_drug_classes") == summary_metrics["distinct_provisional_classes"]
        and summary_int("total_patients_with_any_normalized_name") == patient_metrics["patients_with_any_normalized_name"]
        and summary_int("total_patients_with_any_unknown_or_review_needed_mapping")
        == patient_metrics["patients_with_any_unknown_or_review_needed_mapping"]
        and summary_int("total_patients_with_any_unmapped_raw_name")
        == patient_metrics["patients_with_any_unmapped_raw_name"]
        and summary_int("total_patients_with_any_review_needed_name_normalization")
        == patient_metrics["patients_with_any_review_needed_name_normalization"]
        and summary_int("total_patients_with_name_normalization_manual_review")
        == patient_metrics["patients_with_name_normalization_manual_review"]
        and summary_int("total_patients_with_class_mapping_manual_review")
        == patient_metrics["patients_with_class_mapping_manual_review"]
        and summary_int("total_patients_with_single_drug_class")
        == patient_metrics["patients_with_single_drug_class"]
        and summary_int("total_patients_with_multiple_drug_classes")
        == patient_metrics["patients_with_multiple_drug_classes"]
        and summary_int("total_patients_with_unknown_only_drug_class")
        == patient_metrics["patients_with_unknown_only_drug_class"]
        and summary_int("total_patients_with_missing_like_only_raw_drug_names")
        == patient_metrics["patients_with_missing_like_only_raw_drug_names"]
    )
    readiness_allowed = summary_lookup.get("provisional_readiness_interpretation", "") in {
        READINESS_NOT_READY,
        READINESS_DRUG_CLASS_REVIEW,
        READINESS_TREATMENT_ARM_FREEZE_REVIEW_NEXT,
    }
    class_vocab_allowed = all(str(row["provisional_drug_class"]) in PROVISIONAL_CLASS_VOCABULARY for row in class_rows)
    normalization_confidence_allowed = all(
        str(row["normalization_confidence"]) in NORMALIZATION_CONFIDENCE_ALLOWED for row in normalization_rows
    )
    class_confidence_allowed = all(
        str(row["class_mapping_confidence"]) in CLASS_MAPPING_CONFIDENCE_ALLOWED for row in class_rows
    )
    mapping_status_allowed = all(str(row["mapping_status"]) in MAPPING_STATUS_ALLOWED for row in inventory_rows)
    patient_single_or_multi_allowed = all(
        str(row["single_or_multi_drug_class"]) in SINGLE_OR_MULTI_ALLOWED for row in patient_rows
    )
    passed = all(
        [
            len(workflow_inputs.patient_treatment_profile_rows) > 0,
            len(workflow_inputs.patient_treatment_grouping_rows) > 0,
            len(workflow_inputs.clinical_drug_rows) > 0,
            len(treated_profile_rows) > 0,
            len(inventory_rows) > 0,
            len(normalization_rows) > 0,
            len(class_rows) > 0,
            len(patient_rows) > 0,
            len(inventory_rows) == len(set(inventory_raw_names)),
            len(normalization_rows) == len(set(normalization_raw_names)),
            sorted(inventory_raw_names) == sorted(normalization_raw_names),
            len(class_rows) == len(set(class_normalized_names)),
            sorted(class_normalized_names) == sorted(normalized_names_from_map),
            len(patient_rows) == len(treated_profile_rows),
            len(patient_barcodes) == len(set(patient_barcodes)),
            sorted(treated_barcodes) == sorted(patient_barcodes) == clinical_drug_barcodes,
            patient_metrics["patient_profile_row_order_preserved"],
            summary_counts_reconcile,
            readiness_allowed,
            class_vocab_allowed,
            normalization_confidence_allowed,
            class_confidence_allowed,
            mapping_status_allowed,
            patient_single_or_multi_allowed,
        ]
    )
    return {
        "passed": passed,
        "required_upstream_pointers_found": True,
        "required_source_tables_found": True,
        "patient_treatment_profile_run_log_completed": True,
        "patient_treatment_grouping_run_log_completed": True,
        "clinical_biotabs_run_log_completed": True,
        "clinical_drug_row_count_positive": len(workflow_inputs.clinical_drug_rows) > 0,
        "treated_patient_count_positive": len(treated_profile_rows) > 0,
        "output_rows_positive": all(len(rows) > 0 for rows in [inventory_rows, normalization_rows, class_rows, patient_rows, summary_rows]),
        "every_distinct_raw_drug_name_appears_once_in_inventory": len(inventory_rows) == len(set(inventory_raw_names)),
        "every_distinct_raw_drug_name_appears_once_in_normalization_map": len(normalization_rows) == len(set(normalization_raw_names)),
        "inventory_and_normalization_map_cover_same_raw_names": sorted(inventory_raw_names) == sorted(normalization_raw_names),
        "every_normalized_name_appears_once_in_class_map": len(class_rows) == len(set(class_normalized_names)),
        "class_map_covers_all_normalized_names": sorted(class_normalized_names) == sorted(normalized_names_from_map),
        "treated_patient_set_matches_clinical_drug_patient_set": sorted(treated_barcodes) == clinical_drug_barcodes,
        "patient_profile_row_count_matches_treated_patient_set": len(patient_rows) == len(treated_profile_rows),
        "no_patient_profile_duplicate_barcodes": len(patient_barcodes) == len(set(patient_barcodes)),
        "patient_profile_row_order_preserved_from_treated_profile": patient_metrics["patient_profile_row_order_preserved"],
        "summary_counts_reconcile_to_inventory_class_map_and_patient_profile": summary_counts_reconcile,
        "normalization_confidence_values_allowed": normalization_confidence_allowed,
        "inventory_mapping_status_values_allowed": mapping_status_allowed,
        "class_mapping_confidence_values_allowed": class_confidence_allowed,
        "class_values_restricted_to_allowed_vocabulary": class_vocab_allowed,
        "patient_single_or_multi_values_allowed": patient_single_or_multi_allowed,
        "summary_readiness_interpretation_allowed": readiness_allowed,
        "no_prior_run_overwrite": True,
        "latest_pointer_written_after_success_only": True,
        "clinical_drug_row_count": len(workflow_inputs.clinical_drug_rows),
        "treated_patient_count": len(treated_profile_rows),
        "inventory_row_count": len(inventory_rows),
        "normalization_map_row_count": len(normalization_rows),
        "class_map_row_count": len(class_rows),
        "patient_profile_row_count": len(patient_rows),
        "summary_row_count": len(summary_rows),
    }


def build_latest_pointer_payload(
    *,
    run_id: str,
    paths: WorkflowPaths,
    workflow_inputs: WorkflowInputs,
    output_paths: dict[str, Path],
) -> dict[str, Any]:
    return {
        "updated_at_utc": format_utc_timestamp(utc_now()),
        "drug_name_normalization_v1_run_id": run_id,
        "patient_treatment_profile_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["patient_treatment_profile_v1_run_id"]),
        "patient_treatment_grouping_v1_run_id": str(workflow_inputs.patient_treatment_grouping_latest_pointer["patient_treatment_grouping_v1_run_id"]),
        "treatment_os_overlap_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["treatment_os_overlap_v1_run_id"]),
        "os_endpoint_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["os_endpoint_v1_run_id"]),
        "clinical_biotab_parse_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["clinical_biotab_parse_run_id"]),
        "clinical_source_run_id": str(workflow_inputs.clinical_biotabs_latest_pointer["source_run_id"]),
        "baseline_model_input_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["baseline_model_input_v1_run_id"]),
        "cohort_v1_build_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["cohort_v1_build_id"]),
        "processed_run_directory": repo_relative(output_paths["processed_run_directory"], paths.repo_root),
        "audit_run_directory": repo_relative(output_paths["audit_run_directory"], paths.repo_root),
        "patient_drug_class_profile_v1_tsv": repo_relative(output_paths["patient_drug_class_profile_v1_tsv"], paths.repo_root),
        "drug_name_inventory_v1_tsv": repo_relative(output_paths["drug_name_inventory_v1_tsv"], paths.repo_root),
        "drug_name_normalization_map_v1_tsv": repo_relative(output_paths["drug_name_normalization_map_v1_tsv"], paths.repo_root),
        "drug_name_class_map_v1_tsv": repo_relative(output_paths["drug_name_class_map_v1_tsv"], paths.repo_root),
        "drug_name_normalization_v1_summary_tsv": repo_relative(output_paths["drug_name_normalization_v1_summary_tsv"], paths.repo_root),
        "run_log_json": repo_relative(output_paths["run_log_json"], paths.repo_root),
        "patient_treatment_profile_v1_latest_json": repo_relative(paths.patient_treatment_profile_latest_pointer, paths.repo_root),
        "patient_treatment_grouping_v1_latest_json": repo_relative(paths.patient_treatment_grouping_latest_pointer, paths.repo_root),
        "clinical_biotabs_latest_json": repo_relative(paths.clinical_biotabs_latest_pointer, paths.repo_root),
    }


def write_failure_log(path: Path, payload: dict[str, Any], helper_module: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    helper_module.write_json(path, payload, overwrite=True)


def run_workflow() -> dict[str, Any]:
    started_at = utc_now()
    run_id = started_at.strftime("%Y%m%dT%H%M%SZ")
    helper_module = load_helper_module()
    paths = build_workflow_paths(helper_module)
    processed_run_dir = paths.processed_runs_root / run_id
    audit_run_dir = paths.audit_runs_root / run_id
    run_log_path = audit_run_dir / "run_log.json"

    try:
        trial_config = helper_module.load_yaml(paths.trial_config)
        workflow_inputs = load_workflow_inputs(paths, helper_module)
        helper_module.create_run_directory(processed_run_dir)
        helper_module.create_run_directory(audit_run_dir)

        patient_profile_path = processed_run_dir / "patient_drug_class_profile_v1.tsv"
        inventory_path = audit_run_dir / "drug_name_inventory_v1.tsv"
        normalization_map_path = audit_run_dir / "drug_name_normalization_map_v1.tsv"
        class_map_path = audit_run_dir / "drug_name_class_map_v1.tsv"
        summary_path = audit_run_dir / "drug_name_normalization_v1_summary.tsv"

        inventory_rows, normalization_rows, normalization_decisions = build_inventory_and_normalization_rows(run_id, workflow_inputs.clinical_drug_rows)
        class_rows, class_decisions = build_class_rows(run_id, normalization_rows)
        patient_rows, patient_metrics = build_patient_profile_rows(run_id, workflow_inputs, normalization_decisions, class_decisions)
        summary_rows, summary_metrics = build_summary_rows(run_id, workflow_inputs, inventory_rows, normalization_rows, class_rows, patient_rows, patient_metrics)
        validation = validate_outputs(workflow_inputs=workflow_inputs, inventory_rows=inventory_rows, normalization_rows=normalization_rows, class_rows=class_rows, patient_rows=patient_rows, summary_rows=summary_rows, patient_metrics=patient_metrics, summary_metrics=summary_metrics)
        if not validation["passed"]:
            raise DrugNameNormalizationV1Error("Output validation did not pass. Failed checks: " + str({key: value for key, value in validation.items() if value is False}))

        helper_module.write_dict_rows_tsv(patient_profile_path, PATIENT_PROFILE_FIELDNAMES, patient_rows)
        helper_module.write_dict_rows_tsv(inventory_path, INVENTORY_FIELDNAMES, inventory_rows)
        helper_module.write_dict_rows_tsv(normalization_map_path, NORMALIZATION_MAP_FIELDNAMES, normalization_rows)
        helper_module.write_dict_rows_tsv(class_map_path, CLASS_MAP_FIELDNAMES, class_rows)
        helper_module.write_dict_rows_tsv(summary_path, SUMMARY_FIELDNAMES, summary_rows)

        output_paths = {
            "processed_run_directory": processed_run_dir,
            "audit_run_directory": audit_run_dir,
            "patient_drug_class_profile_v1_tsv": patient_profile_path,
            "drug_name_inventory_v1_tsv": inventory_path,
            "drug_name_normalization_map_v1_tsv": normalization_map_path,
            "drug_name_class_map_v1_tsv": class_map_path,
            "drug_name_normalization_v1_summary_tsv": summary_path,
            "run_log_json": run_log_path,
        }
        latest_pointer_payload = build_latest_pointer_payload(run_id=run_id, paths=paths, workflow_inputs=workflow_inputs, output_paths=output_paths)

        completed_at = utc_now()
        run_log_payload: dict[str, Any] = {
            "status": "completed",
            "drug_name_normalization_v1_run_id": run_id,
            "patient_treatment_profile_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["patient_treatment_profile_v1_run_id"]),
            "patient_treatment_grouping_v1_run_id": str(workflow_inputs.patient_treatment_grouping_latest_pointer["patient_treatment_grouping_v1_run_id"]),
            "treatment_os_overlap_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["treatment_os_overlap_v1_run_id"]),
            "os_endpoint_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["os_endpoint_v1_run_id"]),
            "clinical_biotab_parse_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["clinical_biotab_parse_run_id"]),
            "clinical_source_run_id": str(workflow_inputs.clinical_biotabs_latest_pointer["source_run_id"]),
            "baseline_model_input_v1_run_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["baseline_model_input_v1_run_id"]),
            "cohort_v1_build_id": str(workflow_inputs.patient_treatment_profile_latest_pointer["cohort_v1_build_id"]),
            "started_at_utc": format_utc_timestamp(started_at),
            "completed_at_utc": format_utc_timestamp(completed_at),
            "repo_root": str(paths.repo_root.resolve()),
            "trial_name": trial_config.get("trial_name"),
            "dataset_scope": trial_config.get("dataset_scope"),
            "inputs": {
                "trial_config_yaml": repo_relative(paths.trial_config, paths.repo_root),
                "patient_treatment_profile_v1_latest_json": repo_relative(paths.patient_treatment_profile_latest_pointer, paths.repo_root),
                "patient_treatment_grouping_v1_latest_json": repo_relative(paths.patient_treatment_grouping_latest_pointer, paths.repo_root),
                "clinical_biotabs_latest_json": repo_relative(paths.clinical_biotabs_latest_pointer, paths.repo_root),
                **{key: repo_relative(path, paths.repo_root) for key, path in workflow_inputs.input_paths.items()},
            },
            "outputs": {
                "processed_run_directory": repo_relative(processed_run_dir, paths.repo_root),
                "audit_run_directory": repo_relative(audit_run_dir, paths.repo_root),
                "patient_drug_class_profile_v1_tsv": repo_relative(patient_profile_path, paths.repo_root),
                "drug_name_inventory_v1_tsv": repo_relative(inventory_path, paths.repo_root),
                "drug_name_normalization_map_v1_tsv": repo_relative(normalization_map_path, paths.repo_root),
                "drug_name_class_map_v1_tsv": repo_relative(class_map_path, paths.repo_root),
                "drug_name_normalization_v1_summary_tsv": repo_relative(summary_path, paths.repo_root),
                "run_log_json": repo_relative(run_log_path, paths.repo_root),
                "latest_pointer_json": repo_relative(paths.latest_pointer, paths.repo_root),
            },
            "validation": validation,
            "rules": {
                "unit_of_analysis": "drug_name_inventory_and_treated_patient_profile",
                "input_layer": "saved_patient_treatment_profile_v1_plus_saved_patient_treatment_grouping_v1_plus_saved_clinical_drug_biotab",
                "output_layer": "drug_name_normalization_v1",
                "deterministic_mapping_only": True,
                "no_fuzzy_matching": True,
                "no_new_config_files": True,
                "raw_drug_rows_untouched": True,
                "compound_raw_names_not_split": True,
                "treatment_groups_not_redefined": True,
                "treatment_arms_not_frozen": True,
                "no_modeling": True,
                "no_causal_analysis": True,
                "unknown_preferred_over_wrong_certainty": True,
                "missing_like_normalization": MISSING_LIKE_NORMALIZATION,
                "allowed_normalization_confidence_values": sorted(NORMALIZATION_CONFIDENCE_ALLOWED),
                "allowed_class_mapping_confidence_values": sorted(CLASS_MAPPING_CONFIDENCE_ALLOWED),
                "allowed_mapping_status_values": sorted(MAPPING_STATUS_ALLOWED),
                "allowed_single_or_multi_drug_class_values": sorted(SINGLE_OR_MULTI_ALLOWED),
                "provisional_drug_class_vocabulary": PROVISIONAL_CLASS_VOCABULARY,
            },
            "counts": {
                "clinical_drug_row_count": len(workflow_inputs.clinical_drug_rows),
                "treated_patient_count": patient_metrics["treated_patient_count"],
                "inventory_row_count": len(inventory_rows),
                "normalization_map_row_count": len(normalization_rows),
                "class_map_row_count": len(class_rows),
                "patient_profile_row_count": len(patient_rows),
                "summary_row_count": len(summary_rows),
                "total_mapped_confidently": summary_metrics["inventory_counts"].get("mapped_confidently", 0),
                "total_mapped_but_review": summary_metrics["inventory_counts"].get("mapped_but_review", 0),
                "total_unmapped": summary_metrics["inventory_counts"].get("unmapped", 0),
                "drug_rows_mapped_confidently": summary_metrics["drug_rows_mapped_confidently"],
                "drug_rows_mapped_but_review": summary_metrics["drug_rows_mapped_but_review"],
                "drug_rows_unmapped": summary_metrics["drug_rows_unmapped"],
                "total_distinct_provisional_drug_classes": summary_metrics["distinct_provisional_classes"],
                "patients_with_any_unknown_or_review_needed_mapping": patient_metrics["patients_with_any_unknown_or_review_needed_mapping"],
                "patients_with_single_drug_class": patient_metrics["patients_with_single_drug_class"],
                "patients_with_multiple_drug_classes": patient_metrics["patients_with_multiple_drug_classes"],
                "patients_with_unknown_only_drug_class": patient_metrics["patients_with_unknown_only_drug_class"],
                "patients_with_missing_like_only_raw_drug_names": patient_metrics["patients_with_missing_like_only_raw_drug_names"],
                "provisional_readiness_interpretation": summary_metrics["readiness"],
            },
            "latest_pointer": latest_pointer_payload,
            "upstream_snapshots": {
                "patient_treatment_profile_v1_latest_pointer": workflow_inputs.patient_treatment_profile_latest_pointer,
                "patient_treatment_grouping_v1_latest_pointer": workflow_inputs.patient_treatment_grouping_latest_pointer,
                "clinical_biotabs_latest_pointer": workflow_inputs.clinical_biotabs_latest_pointer,
            },
        }
        helper_module.write_json(run_log_path, run_log_payload)
        helper_module.write_json(paths.latest_pointer, latest_pointer_payload, overwrite=True)
        return run_log_payload
    except Exception as exc:
        failure_payload = {
            "status": "failed",
            "drug_name_normalization_v1_run_id": run_id,
            "started_at_utc": format_utc_timestamp(started_at),
            "failed_at_utc": format_utc_timestamp(utc_now()),
            "error_type": type(exc).__name__,
            "error_message": str(exc),
            "workflow": "tcga_brca_drug_name_normalization_v1",
        }
        write_failure_log(run_log_path, failure_payload, helper_module)
        raise


def print_summary(run_log: dict[str, Any]) -> None:
    counts = run_log.get("counts", {})
    outputs = run_log.get("outputs", {})
    print("TCGA-BRCA drug-name normalization v1 completed")
    print(f"  Run ID                   : {run_log['drug_name_normalization_v1_run_id']}")
    print(f"  Profile run ID           : {run_log['patient_treatment_profile_v1_run_id']}")
    print(f"  Grouping run ID          : {run_log['patient_treatment_grouping_v1_run_id']}")
    print(f"  Clinical drug rows       : {counts.get('clinical_drug_row_count', 0)}")
    print(f"  Treated patients         : {counts.get('treated_patient_count', 0)}")
    print(f"  Distinct raw drug names  : {counts.get('inventory_row_count', 0)}")
    print(f"  Mapped confidently       : {counts.get('total_mapped_confidently', 0)}")
    print(f"  Mapped but review        : {counts.get('total_mapped_but_review', 0)}")
    print(f"  Unmapped                 : {counts.get('total_unmapped', 0)}")
    print(f"  Review-burden patients   : {counts.get('patients_with_any_unknown_or_review_needed_mapping', 0)}")
    print(f"  Readiness                : {counts.get('provisional_readiness_interpretation', '')}")
    print(f"  Processed run dir        : {outputs.get('processed_run_directory', '')}")
    print(f"  Audit run dir            : {outputs.get('audit_run_directory', '')}")
    print(f"  Latest pointer           : {outputs.get('latest_pointer_json', '')}")


def main() -> int:
    run_log = run_workflow()
    print_summary(run_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
