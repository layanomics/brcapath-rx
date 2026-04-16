"""
01_audit_phase1.py — Phase 1 Data Audit
========================================
Audits downloaded TCGA-BRCA XML files and METABRIC clinical TSV/JSON files.

TASK A: TCGA-BRCA XML Audit
TASK B: METABRIC Clinical Audit
"""

import json
import os
import re
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd

# ─── Paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]

TCGA_XML_DIR  = ROOT / "01-data" / "raw" / "tcga-brca" / "gdc" / "clinical"
TCGA_AUDIT_OUT = ROOT / "01-data" / "audit" / "tcga-brca" / "tcga_clinical_audit.json"

META_CLINICAL_DIR  = ROOT / "01-data" / "raw" / "metabric" / "clinical"
META_AUDIT_OUT     = ROOT / "01-data" / "audit" / "metabric" / "metabric_clinical_audit.json"

# ─── Treatment-relevant keyword pattern ───────────────────────────────────────
TREATMENT_KEYWORDS = re.compile(
    r"treatment|drug|therapy|radiation|chemo|chemotherapy|pharmaceutical|"
    r"agent|dose|response|outcome|survival|recurrence|follow|death|vital",
    re.IGNORECASE,
)


# ══════════════════════════════════════════════════════════════════════════════
# TASK A — TCGA-BRCA XML AUDIT
# ══════════════════════════════════════════════════════════════════════════════

def get_file_type(filename: str) -> str:
    """Derive a file-type key from everything before the first patient ID token."""
    # Filenames look like: nationwidechildrens.org_ssf.TCGA-XX-XXXX.xml
    # The patient ID segment starts with TCGA-
    parts = filename.split(".")
    type_parts = []
    for part in parts:
        if part.startswith("TCGA-") or part.lower() == "xml":
            break
        type_parts.append(part)
    return ".".join(type_parts) if type_parts else filename


def collect_xml_fields(path: Path) -> dict[str, set]:
    """
    Parse one XML file and return a dict mapping
    top-level element tag → set of child field names.
    """
    try:
        tree = ET.parse(path)
    except ET.ParseError as exc:
        print(f"  [WARN] Failed to parse {path.name}: {exc}")
        return {}

    root = tree.getroot()
    # Strip namespace from tag: {ns}tag → tag
    def strip_ns(tag: str) -> str:
        return re.sub(r"\{[^}]*\}", "", tag)

    structure: dict[str, set] = defaultdict(set)

    def walk(element, parent_tag):
        tag = strip_ns(element.tag)
        for child in element:
            child_tag = strip_ns(child.tag)
            structure[parent_tag].add(child_tag)
            # Recurse one more level so nested field names are captured too
            walk(child, child_tag)

    root_tag = strip_ns(root.tag)
    for top_elem in root:
        top_tag = strip_ns(top_elem.tag)
        structure[root_tag].add(top_tag)
        walk(top_elem, top_tag)

    return dict(structure)


def collect_all_leaf_fields(path: Path) -> set[str]:
    """
    Return every unique leaf-level field name in the XML
    (i.e., element tags that have text content or no children).
    """
    try:
        tree = ET.parse(path)
    except ET.ParseError:
        return set()

    def strip_ns(tag: str) -> str:
        return re.sub(r"\{[^}]*\}", "", tag)

    fields: set[str] = set()

    def walk(element):
        tag = strip_ns(element.tag)
        children = list(element)
        if not children:
            # Leaf node
            fields.add(tag)
        else:
            fields.add(tag)
            for child in children:
                walk(child)

    walk(tree.getroot())
    return fields


def run_tcga_audit():
    print("\n" + "=" * 60)
    print("TASK A — TCGA-BRCA XML AUDIT")
    print("=" * 60)

    # ── A1: Inventory ─────────────────────────────────────────────────────────
    print("\n[A1] Inventory")
    xml_files = sorted(TCGA_XML_DIR.glob("*.xml"))
    if not xml_files:
        print(f"  ERROR: No XML files found in {TCGA_XML_DIR}")
        return {}

    type_to_files: dict[str, list[Path]] = defaultdict(list)
    for f in xml_files:
        ftype = get_file_type(f.name)
        type_to_files[ftype].append(f)

    file_type_counts = {ft: len(files) for ft, files in sorted(type_to_files.items())}
    print(f"  Total XML files : {len(xml_files)}")
    print(f"  File types found: {len(file_type_counts)}")
    for ft, cnt in file_type_counts.items():
        print(f"    {ft:60s}  {cnt:4d} files")

    # ── A2: Parse a sample per type ───────────────────────────────────────────
    print("\n[A2] Sample parse — one file per type")
    sample_structures: dict[str, dict] = {}
    for ftype, files in sorted(type_to_files.items()):
        sample_path = files[0]
        structure = collect_xml_fields(sample_path)
        sample_structures[ftype] = {k: sorted(v) for k, v in structure.items()}
        print(f"\n  Type: {ftype}")
        print(f"  File: {sample_path.name}")
        for elem, children in sorted(structure.items()):
            print(f"    <{elem}>  ({len(children)} sub-elements)")
            for child in sorted(children)[:20]:
                print(f"      · {child}")
            if len(children) > 20:
                print(f"      … and {len(children) - 20} more")

    # ── A3: All fields across all files ───────────────────────────────────────
    print("\n[A3] Extracting all fields across all files …")
    all_fields_per_type: dict[str, dict] = {}

    for ftype, files in sorted(type_to_files.items()):
        n_files = len(files)
        field_file_counts: dict[str, int] = defaultdict(int)

        for fpath in files:
            fields = collect_all_leaf_fields(fpath)
            for field in fields:
                field_file_counts[field] += 1

        # Build coverage table
        coverage = {
            field: {
                "files_present": cnt,
                "coverage_pct": round(cnt / n_files * 100, 1),
            }
            for field, cnt in sorted(field_file_counts.items())
        }
        all_fields_per_type[ftype] = coverage
        print(f"  {ftype}: {len(coverage)} unique fields across {n_files} files")

    # ── A4: Treatment-relevant fields ─────────────────────────────────────────
    print("\n[A4] Treatment-relevant fields")
    treatment_relevant: dict[str, list[str]] = {}
    for ftype, coverage in all_fields_per_type.items():
        hits = [f for f in coverage if TREATMENT_KEYWORDS.search(f)]
        treatment_relevant[ftype] = sorted(hits)
        if hits:
            print(f"\n  [{ftype}]")
            for h in hits:
                pct = coverage[h]["coverage_pct"]
                print(f"    {h}  ({pct}% files)")

    # ── A5: Save TCGA audit JSON ───────────────────────────────────────────────
    TCGA_AUDIT_OUT.parent.mkdir(parents=True, exist_ok=True)
    tcga_audit = {
        "file_type_counts": file_type_counts,
        "all_fields_per_type": all_fields_per_type,
        "treatment_relevant_fields": treatment_relevant,
    }
    with open(TCGA_AUDIT_OUT, "w") as fh:
        json.dump(tcga_audit, fh, indent=2)
    print(f"\n[A5] Saved → {TCGA_AUDIT_OUT}")

    return tcga_audit


# ══════════════════════════════════════════════════════════════════════════════
# TASK B — METABRIC CLINICAL AUDIT
# ══════════════════════════════════════════════════════════════════════════════

def run_metabric_audit():
    print("\n" + "=" * 60)
    print("TASK B — METABRIC CLINICAL AUDIT")
    print("=" * 60)

    # ── B1: Wide TSV ──────────────────────────────────────────────────────────
    print("\n[B1] metabric_clinical_wide.tsv")
    wide_path = META_CLINICAL_DIR / "metabric_clinical_wide.tsv"
    df_wide = pd.read_csv(wide_path, sep="\t")
    shape = df_wide.shape
    print(f"  Shape  : {shape[0]} rows × {shape[1]} columns")
    print(f"  Columns: {list(df_wide.columns)}")
    print("\n  Dtypes and % missing:")
    missing_pct = (df_wide.isnull().mean() * 100).round(2)
    dtype_info = df_wide.dtypes.astype(str)
    col_report = pd.DataFrame({"dtype": dtype_info, "missing_%": missing_pct})
    print(col_report.to_string())

    wide_columns_set = set(df_wide.columns)

    # ── B2: Clinical attributes JSON ──────────────────────────────────────────
    print("\n[B2] metabric_clinical_attributes.json — all 36 attributes")
    attr_path = META_CLINICAL_DIR / "metabric_clinical_attributes.json"
    with open(attr_path) as fh:
        attributes: list[dict] = json.load(fh)

    print(f"  {'clinicalAttributeId':40s}  {'displayName':40s}  {'datatype':12s}  patientAttr")
    print("  " + "-" * 110)
    for attr in attributes:
        aid   = attr.get("clinicalAttributeId", "")
        dname = attr.get("displayName", "")
        dtype = attr.get("datatype", "")
        pat   = attr.get("patientAttribute", "")
        print(f"  {aid:40s}  {dname:40s}  {dtype:12s}  {pat}")

    # ── B3: Explain the 13 vs 36 gap ──────────────────────────────────────────
    print("\n[B3] Gap analysis: 36 attributes vs 13 TSV columns")

    # Load long-format clinical data to understand what's actually present per
    # attributeId and whether it is patient-level or sample-level
    clinical_data_path = META_CLINICAL_DIR / "metabric_clinical_data.json"
    with open(clinical_data_path) as fh:
        clinical_data: list[dict] = json.load(fh)

    # Build a set of attributeIds that appear in the long-format data
    attrs_in_long: set[str] = {rec.get("clinicalAttributeId", "") for rec in clinical_data}

    # Build presence table for all 36 attributes
    gap_table: list[dict] = []
    for attr in attributes:
        aid       = attr.get("clinicalAttributeId", "")
        dname     = attr.get("displayName", "")
        pat_attr  = attr.get("patientAttribute", False)  # True = patient-level
        dtype     = attr.get("datatype", "")

        in_wide = aid in wide_columns_set

        # Determine why it might be missing from the wide TSV
        if in_wide:
            reason = "present in wide TSV"
        elif pat_attr is False:
            # Sample-level attributes might not have been pivoted
            reason = "sample-level — not merged into patient-wide TSV"
        else:
            # Patient-level but absent: check if it exists in long data
            if aid in attrs_in_long:
                reason = "patient-level — exists in long-format data but not pivoted to wide TSV"
            else:
                reason = "patient-level — not found in downloaded long-format data"

        gap_table.append({
            "clinicalAttributeId": aid,
            "displayName": dname,
            "datatype": dtype,
            "patientAttribute": pat_attr,
            "in_wide_tsv": in_wide,
            "reason": reason,
        })

    # Print table
    present_count = sum(1 for r in gap_table if r["in_wide_tsv"])
    missing_count = len(gap_table) - present_count

    print(f"\n  Present in wide TSV : {present_count}")
    print(f"  Missing from wide TSV: {missing_count}")
    print()
    hdr = f"  {'clinicalAttributeId':35s}  {'patientAttr':11s}  {'inWide':6s}  reason"
    print(hdr)
    print("  " + "-" * 110)
    for row in gap_table:
        pid_str = str(row["patientAttribute"])
        in_w    = "YES" if row["in_wide_tsv"] else "NO"
        print(f"  {row['clinicalAttributeId']:35s}  {pid_str:11s}  {in_w:6s}  {row['reason']}")

    # Summarise the gap by category
    reasons: dict[str, int] = defaultdict(int)
    for row in gap_table:
        if not row["in_wide_tsv"]:
            reasons[row["reason"]] += 1
    print("\n  Gap breakdown:")
    for reason, cnt in reasons.items():
        print(f"    ({cnt}) {reason}")

    # ── B4: Treatment-relevant attributes ─────────────────────────────────────
    print("\n[B4] Treatment-relevant attributes")
    treatment_attrs: list[dict] = []
    for attr in attributes:
        aid   = attr.get("clinicalAttributeId", "")
        dname = attr.get("displayName", "")
        if TREATMENT_KEYWORDS.search(aid) or TREATMENT_KEYWORDS.search(dname):
            treatment_attrs.append({
                "clinicalAttributeId": aid,
                "displayName": dname,
                "patientAttribute": attr.get("patientAttribute", ""),
                "in_wide_tsv": aid in wide_columns_set,
            })

    if treatment_attrs:
        for a in treatment_attrs:
            flag = "✓ in wide TSV" if a["in_wide_tsv"] else "✗ missing from wide TSV"
            print(f"  {a['clinicalAttributeId']:35s}  {a['displayName']:40s}  {flag}")
    else:
        print("  None found.")

    # ── B5: Save METABRIC audit JSON ──────────────────────────────────────────
    META_AUDIT_OUT.parent.mkdir(parents=True, exist_ok=True)
    meta_audit = {
        "wide_tsv_shape": {"rows": shape[0], "cols": shape[1]},
        "all_36_attributes_with_presence": gap_table,
        "gap_explanation": {
            "total_attributes": len(attributes),
            "present_in_wide_tsv": present_count,
            "missing_from_wide_tsv": missing_count,
            "breakdown": dict(reasons),
        },
        "treatment_relevant_attributes": treatment_attrs,
    }
    with open(META_AUDIT_OUT, "w") as fh:
        json.dump(meta_audit, fh, indent=2)
    print(f"\n[B5] Saved → {META_AUDIT_OUT}")

    return meta_audit, df_wide, attributes


# ══════════════════════════════════════════════════════════════════════════════
# COMBINED SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

def print_summary(tcga_audit: dict, meta_audit: dict, df_wide: pd.DataFrame):
    print("\n" + "=" * 60)
    print("COMBINED SUMMARY")
    print("=" * 60)

    # TCGA
    print("\nTCGA-BRCA:")
    for ftype, cnt in tcga_audit["file_type_counts"].items():
        print(f"  File type '{ftype}': {cnt} files")

    total_unique_fields = sum(
        len(fields) for fields in tcga_audit["all_fields_per_type"].values()
    )
    print(f"  Total unique fields discovered: {total_unique_fields}")

    print("  Treatment-relevant fields found:")
    any_found = False
    for ftype, fields in tcga_audit["treatment_relevant_fields"].items():
        if fields:
            any_found = True
            for f in fields:
                print(f"    [{ftype}] {f}")
    if not any_found:
        print("    (none)")

    # METABRIC
    print("\nMETABRIC:")
    shape = meta_audit["wide_tsv_shape"]
    print(f"  Wide TSV shape: {shape['rows']} rows × {shape['cols']} cols")

    gap_info = meta_audit["gap_explanation"]
    print(f"  Attributes present in wide TSV : {gap_info['present_in_wide_tsv']} / {gap_info['total_attributes']}")
    print(f"  Attributes missing from wide TSV: {gap_info['missing_from_wide_tsv']} / {gap_info['total_attributes']}")
    for reason, cnt in gap_info["breakdown"].items():
        print(f"    ({cnt}) {reason}")

    print("  Treatment-relevant attributes found:")
    if meta_audit["treatment_relevant_attributes"]:
        for a in meta_audit["treatment_relevant_attributes"]:
            in_w = "in wide TSV" if a["in_wide_tsv"] else "NOT in wide TSV"
            print(f"    {a['clinicalAttributeId']} — {a['displayName']} [{in_w}]")
    else:
        print("    (none)")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    tcga_audit = run_tcga_audit()
    meta_audit, df_wide, attributes = run_metabric_audit()
    print_summary(tcga_audit, meta_audit, df_wide)
