"""
00_discover_sources.py

Queries GDC, UCSC Xena (GDC hub), and cBioPortal to produce a full
inventory of available data for TCGA-BRCA and METABRIC.
Downloads NO genomic data files.

Outputs:
  01-data/audit/source_discovery_report.json
  01-data/audit/source_discovery_report.md
"""

import json
import os
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
AUDIT_DIR = PROJECT_ROOT / "01-data" / "audit"
AUDIT_DIR.mkdir(parents=True, exist_ok=True)

JSON_OUT = AUDIT_DIR / "source_discovery_report.json"
MD_OUT = AUDIT_DIR / "source_discovery_report.md"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
TIMEOUT = 60  # seconds per request


def ts() -> str:
    return datetime.now(timezone.utc).isoformat()


def get_json(url: str, params: dict = None, method: str = "GET", body: dict = None) -> dict:
    """Thin wrapper — returns parsed JSON or raises with context."""
    try:
        if method.upper() == "POST":
            resp = requests.post(url, json=body, params=params, timeout=TIMEOUT)
        else:
            resp = requests.get(url, params=params, timeout=TIMEOUT)
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        raise RuntimeError(f"[{ts()}] HTTP error for {url}: {exc}") from exc


def safe(label: str, fn):
    """Call fn(); on exception return an error dict and print the traceback."""
    try:
        return fn()
    except Exception as exc:
        msg = f"[{ts()}] ERROR in {label}: {exc}"
        print(msg, file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return {"error": msg}


# ---------------------------------------------------------------------------
# TASK A  —  GDC API
# ---------------------------------------------------------------------------
GDC_BASE = "https://api.gdc.cancer.gov"

# Fields we consider "categorical" for faceting (exclude id/uuid/date fields)
_SKIP_SUFFIXES = (
    "_id", ".id", "uuid", "submitter_id", "file_name", "file_path",
    "created_datetime", "updated_datetime", "released", "state",
    "md5sum", "file_size",
)


def _is_categorical(field: str) -> bool:
    fl = field.lower()
    for s in _SKIP_SUFFIXES:
        if fl.endswith(s):
            return False
    # Skip pure numeric / byte-count fields
    for tok in ("size", "count", "number", "percent", "index", "position",
                "start", "end", "length", "score"):
        if fl.endswith(f".{tok}") or fl == tok:
            return False
    return True


def _extract_fields(mapping_response: dict) -> list:
    """Return flat list of field names from a GDC _mapping response.

    The GDC _mapping endpoint returns a JSON object with a top-level "fields"
    key whose value is already a flat list of all searchable field-name strings.
    Use that directly instead of walking the nested "_mapping" structure.
    """
    # Primary path: GDC includes a ready-made flat list under "fields"
    flat = mapping_response.get("fields")
    if isinstance(flat, list) and flat:
        print(f"  [_extract_fields] top-level 'fields' key found with "
              f"{len(flat)} entries (sample: {flat[:3]})")
        return sorted(set(flat))

    # Fallback: walk the nested _mapping structure
    print(f"  [_extract_fields] no 'fields' key — falling back to walking "
          f"'_mapping'. Top-level keys: {list(mapping_response.keys())[:10]}")
    fields = []

    def walk(obj, prefix=""):
        if isinstance(obj, dict):
            for k, v in obj.items():
                full = f"{prefix}.{k}" if prefix else k
                if isinstance(v, dict) and v.get("type") not in (None,):
                    # leaf node
                    fields.append(full)
                else:
                    walk(v, full)

    raw = mapping_response.get("_mapping") or mapping_response
    for top_key, top_val in raw.items():
        if isinstance(top_val, dict):
            props = top_val.get("properties") or top_val
            walk(props, top_key)
        else:
            fields.append(top_key)

    return sorted(set(fields))


def _find_project_id_field(fields: list) -> str:
    """Return the field whose name contains both 'project' and 'project_id'."""
    candidates = [f for f in fields
                  if "project" in f.lower() and "project_id" in f.lower()]
    if not candidates:
        raise ValueError(f"No project_id field found in: {fields[:30]}")
    # Prefer shortest (most direct) match
    return sorted(candidates, key=len)[0]


def _gdc_facets(endpoint: str, project_filter_field: str, facet_fields: list) -> dict:
    """Query a GDC endpoint with size=0 and facets; return bucket dict."""
    cat_fields = [f for f in facet_fields if _is_categorical(f)]
    # GDC accepts up to ~100 facets at once; chunk if needed
    CHUNK = 80
    all_buckets = {}
    for i in range(0, len(cat_fields), CHUNK):
        chunk = cat_fields[i: i + CHUNK]
        filter_payload = {
            "op": "=",
            "content": {"field": project_filter_field, "value": "TCGA-BRCA"},
        }
        params = {
            "filters": json.dumps(filter_payload),
            "facets": ",".join(chunk),
            "size": 0,
            "from": 0,
        }
        data = get_json(f"{GDC_BASE}/{endpoint}", params=params)
        aggs = (data.get("data") or {}).get("aggregations") or {}
        for fname, fbody in aggs.items():
            buckets = fbody.get("buckets") or []
            if buckets:
                all_buckets[fname] = buckets
    return all_buckets


def task_a_gdc() -> dict:
    result = {}

    # A1 — project summary
    def a1():
        data = get_json(f"{GDC_BASE}/projects", params={
            "filters": json.dumps({"op": "=", "content": {"field": "project_id", "value": "TCGA-BRCA"}}),
            "size": 1,
        })
        hits = (data.get("data") or {}).get("hits") or []
        return hits[0] if hits else data

    result["a1_project_summary"] = safe("A1 project summary", a1)

    # A2 — file-level fields
    def a2():
        data = get_json(f"{GDC_BASE}/files/_mapping")
        return _extract_fields(data)

    file_fields = safe("A2 file _mapping", a2)
    result["a2_file_fields"] = file_fields

    # A3 — identify project filter field for files
    def a3():
        fields = result["a2_file_fields"]
        if isinstance(fields, dict) and "error" in fields:
            return fields
        return _find_project_id_field(fields)

    file_project_field = safe("A3 file project_id field", a3)
    result["a3_file_project_filter_field"] = file_project_field

    # A4 — file facets
    def a4():
        fields = result["a2_file_fields"]
        ppf = result["a3_file_project_filter_field"]
        if isinstance(fields, dict) and "error" in fields:
            return fields
        if isinstance(ppf, dict) and "error" in ppf:
            return ppf
        return _gdc_facets("files", ppf, fields)

    result["a4_file_facets"] = safe("A4 file facets", a4)

    # A5 — case-level fields
    def a5():
        data = get_json(f"{GDC_BASE}/cases/_mapping")
        return _extract_fields(data)

    case_fields = safe("A5 cases _mapping", a5)
    result["a5_case_fields"] = case_fields

    # A6 — identify project filter field for cases
    def a6():
        fields = result["a5_case_fields"]
        if isinstance(fields, dict) and "error" in fields:
            return fields
        return _find_project_id_field(fields)

    case_project_field = safe("A6 case project_id field", a6)
    result["a6_case_project_filter_field"] = case_project_field

    # A7 — case facets
    def a7():
        fields = result["a5_case_fields"]
        ppf = result["a6_case_project_filter_field"]
        if isinstance(fields, dict) and "error" in fields:
            return fields
        if isinstance(ppf, dict) and "error" in ppf:
            return ppf
        return _gdc_facets("cases", ppf, fields)

    result["a7_case_facets"] = safe("A7 case facets", a7)

    return result


# ---------------------------------------------------------------------------
# TASK B  —  UCSC Xena GDC Hub
# ---------------------------------------------------------------------------
XENA_HUB = "https://gdc.xenahubs.net"


def task_b_xena(gdc_result: dict) -> dict:
    result = {}

    # B1 — list all datasets in TCGA-BRCA cohort
    def b1():
        try:
            import xenaPython as xena
        except ImportError:
            return {"error": "xenaPython not installed — run: pip install xenaPython"}

        # Diagnostic: show available functions so we can verify the API surface
        xena_api = [fn for fn in dir(xena) if not fn.startswith("_")]
        print(f"  [B1] xenaPython public API: {xena_api}")

        # Discover the right cohort name by listing all cohorts on the hub
        # all_cohorts(host, exclude) — pass empty exclude list
        try:
            all_cohorts = xena.all_cohorts(XENA_HUB, [])
            print(f"  [B1] All cohorts on {XENA_HUB}: {all_cohorts}")
        except Exception as exc:
            print(f"  [B1] Could not list cohorts: {exc}")
            all_cohorts = []

        # Try the expected cohort name; fall back to first available cohort
        preferred = "GDC TCGA Breast Cancer (BRCA)"
        if preferred in (all_cohorts or []):
            cohort = preferred
        elif all_cohorts:
            cohort = all_cohorts[0]
            print(f"  [B1] Preferred cohort not found; using first available: {cohort!r}")
        else:
            cohort = preferred  # last resort — let the call fail naturally

        # dataset_list(host, cohorts) — cohorts must be a LIST of strings;
        # returns list of dicts with keys: name, longtitle, count, type, etc.
        raw_datasets = xena.dataset_list(XENA_HUB, [cohort]) or []
        print(f"  [B1] dataset_list returned {len(raw_datasets)} entries for cohort {cohort!r}")
        # Extract dataset IDs (the 'name' field is the dataset identifier)
        if raw_datasets and isinstance(raw_datasets[0], dict):
            datasets = [d.get("name") or str(d) for d in raw_datasets]
        else:
            datasets = raw_datasets

        dataset_info = []
        for ds in datasets:
            try:
                meta = xena.dataset_field_metadata(XENA_HUB, ds) or {}
            except Exception:
                meta = {}
            try:
                samples = xena.dataset_samples(XENA_HUB, ds, None) or []
                sample_count = len(samples)
            except Exception:
                sample_count = None
            dataset_info.append({
                "dataset_id": ds,
                "sample_count": sample_count,
                "metadata": meta if isinstance(meta, dict) else {},
            })
        return {"cohort": cohort, "datasets": dataset_info}

    result["b1_xena_datasets"] = safe("B1 Xena datasets", b1)

    # B2 — gap analysis
    def b2():
        b1_data = result["b1_xena_datasets"]
        if isinstance(b1_data, dict) and "error" in b1_data:
            return b1_data

        # Collect GDC data categories from facets
        a4 = gdc_result.get("a4_file_facets") or {}
        gdc_categories = set()
        for fname, buckets in a4.items():
            if "data_category" in fname.lower():
                for b in buckets:
                    gdc_categories.add(b.get("key", ""))

        # Collect Xena dataset types (from dataset_id naming)
        xena_datasets = b1_data.get("datasets") or []
        xena_ids_lower = " ".join(d["dataset_id"].lower() for d in xena_datasets)

        # Simple heuristic: check if category keyword appears in any Xena id
        _KEYWORD_MAP = {
            "Transcriptome Profiling": ["rnaseq", "htseq", "fpkm", "expression"],
            "Simple Nucleotide Variation": ["mutation", "snv", "muse", "mutect", "somaticsniper", "varscan"],
            "Copy Number Variation": ["cnv", "copy_number", "gistic"],
            "DNA Methylation": ["methylation", "methyl"],
            "miRNA Expression Quantification": ["mirna"],
            "Clinical": ["clinical"],
            "Biospecimen": ["biospecimen"],
        }

        gaps = []
        covered = []
        for cat in sorted(gdc_categories):
            keywords = _KEYWORD_MAP.get(cat, [cat.lower().split()[0]])
            found = any(kw in xena_ids_lower for kw in keywords)
            if found:
                covered.append(cat)
            else:
                gaps.append(cat)

        return {
            "gdc_data_categories": sorted(gdc_categories),
            "xena_covered": covered,
            "gaps_in_xena": gaps,
        }

    result["b2_gap_analysis"] = safe("B2 gap analysis", b2)
    return result


# ---------------------------------------------------------------------------
# TASK C  —  cBioPortal
# ---------------------------------------------------------------------------
# Base URL must NOT end with /api — the OpenAPI spec paths already start with
# /api/, so appending them to a base that ends with /api would produce the
# double-prefix https://www.cbioportal.org/api/api/studies/...
CBIO_BASE = "https://www.cbioportal.org"
STUDY_ID = "brca_metabric"


def task_c_cbio() -> dict:
    result = {}

    # C1 — discover all per-study endpoints from OpenAPI spec
    def c1():
        spec = get_json("https://www.cbioportal.org/api/v3/api-docs")
        paths = spec.get("paths") or {}
        per_study = {}
        for path, methods in paths.items():
            if "{studyId}" in path:
                per_study[path] = methods
        result["_swagger_spec"] = spec  # cache for C3
        return {
            "total_paths": len(paths),
            "per_study_paths": list(per_study.keys()),
            "per_study_count": len(per_study),
        }

    result["c1_per_study_endpoints"] = safe("C1 cBioPortal endpoints", c1)

    # C2 — study metadata
    def c2():
        data = get_json(f"{CBIO_BASE}/api/studies/{STUDY_ID}")
        return data

    result["c2_study_metadata"] = safe("C2 study metadata", c2)

    # C3 — split & query
    def c3():
        c1_data = result["c1_per_study_endpoints"]
        if isinstance(c1_data, dict) and "error" in c1_data:
            return c1_data

        spec = result.get("_swagger_spec") or {}
        paths_spec = spec.get("paths") or {}
        per_study_paths = c1_data.get("per_study_paths") or []

        group1 = []  # satisfiable with studyId only
        group2 = []  # needs extra required params

        for path in per_study_paths:
            methods = paths_spec.get(path) or {}
            # Collect GET method params (most cBioPortal endpoints use GET)
            method_obj = methods.get("get") or methods.get("post") or {}
            parameters = method_obj.get("parameters") or []

            extra_required = []
            for p in parameters:
                pname = p.get("name") or ""
                prequired = p.get("required", False)
                pin = p.get("in") or ""
                if pname == "studyId":
                    continue
                if prequired and pin == "path":
                    # Required path param that is not studyId
                    ptype = (
                        (p.get("schema") or {}).get("type")
                        or p.get("type")
                        or "unknown"
                    )
                    extra_required.append({"name": pname, "type": ptype, "in": pin})

            if extra_required:
                group2.append({"path": path, "extra_required_params": extra_required})
            else:
                group1.append(path)

        # Query group 1
        group1_results = []
        for path in group1:
            url = f"{CBIO_BASE}{path.replace('{studyId}', STUDY_ID)}"
            try:
                data = get_json(url)
                count = len(data) if isinstance(data, list) else 1
                group1_results.append({
                    "path": path,
                    "url": url,
                    "item_count": count,
                    "response": data,
                    "status": "ok",
                })
            except Exception as exc:
                group1_results.append({
                    "path": path,
                    "url": url,
                    "item_count": 0,
                    "response": None,
                    "status": f"error: {exc}",
                })

        return {
            "group1_queried": group1_results,
            "group2_skipped": group2,
        }

    result["c3_endpoint_results"] = safe("C3 cBioPortal endpoint queries", c3)

    # Remove cached spec from result (too large for JSON output)
    result.pop("_swagger_spec", None)
    return result


# ---------------------------------------------------------------------------
# Report rendering
# ---------------------------------------------------------------------------

def _bucket_table(buckets: list) -> str:
    if not buckets:
        return "_none_\n"
    lines = ["| Value | Count |", "|-------|-------|"]
    for b in sorted(buckets, key=lambda x: -x.get("doc_count", 0)):
        lines.append(f"| {b.get('key', '')} | {b.get('doc_count', '')} |")
    return "\n".join(lines) + "\n"


def render_markdown(report: dict) -> str:
    gdc = report.get("gdc_tcga_brca") or {}
    xena = report.get("xena_gdc_tcga_brca") or {}
    cbio = report.get("cbio_metabric") or {}

    lines = [
        "# Source Discovery Report — BRCAPath-Rx",
        f"\n_Generated: {report.get('discovery_timestamp')}_\n",
    ]

    # ---- Section 1 --------------------------------------------------------
    lines.append("## 1. GDC TCGA-BRCA: Project Summary\n")
    proj = gdc.get("a1_project_summary") or {}
    if isinstance(proj, dict) and "error" not in proj:
        for k, v in proj.items():
            lines.append(f"- **{k}**: {v}")
    else:
        lines.append(f"_{proj}_")
    lines.append("")

    # ---- Section 2 --------------------------------------------------------
    lines.append("## 2. GDC TCGA-BRCA: File Fields Discovered via _mapping\n")
    file_fields = gdc.get("a2_file_fields") or []
    ppf = gdc.get("a3_file_project_filter_field") or "unknown"
    lines.append(f"**Project filter field identified:** `{ppf}`\n")
    lines.append(f"Total fields discovered: {len(file_fields)}\n")
    if isinstance(file_fields, list):
        lines.append("```")
        lines.extend(file_fields)
        lines.append("```")
    else:
        lines.append(f"_{file_fields}_")
    lines.append("")

    # ---- Section 3 --------------------------------------------------------
    lines.append("## 3. GDC TCGA-BRCA: File Inventory by Category\n")
    a4 = gdc.get("a4_file_facets") or {}
    if isinstance(a4, dict) and "error" not in a4:
        for fname, buckets in sorted(a4.items()):
            lines.append(f"### {fname}\n")
            lines.append(_bucket_table(buckets))
    else:
        lines.append(f"_{a4}_")
    lines.append("")

    # ---- Section 4 --------------------------------------------------------
    lines.append("## 4. GDC TCGA-BRCA: Case Inventory\n")
    a7 = gdc.get("a7_case_facets") or {}
    case_ppf = gdc.get("a6_case_project_filter_field") or "unknown"
    lines.append(f"**Case project filter field:** `{case_ppf}`\n")
    if isinstance(a7, dict) and "error" not in a7:
        for fname, buckets in sorted(a7.items()):
            lines.append(f"### {fname}\n")
            lines.append(_bucket_table(buckets))
    else:
        lines.append(f"_{a7}_")
    lines.append("")

    # ---- Section 5 --------------------------------------------------------
    lines.append("## 5. Xena GDC Hub: Available Datasets for TCGA-BRCA\n")
    b1 = xena.get("b1_xena_datasets") or {}
    if isinstance(b1, dict) and "error" not in b1:
        datasets = b1.get("datasets") or []
        lines.append(f"Cohort: **{b1.get('cohort')}**  |  Datasets found: **{len(datasets)}**\n")
        lines.append("| Dataset ID | Sample Count |")
        lines.append("|------------|-------------|")
        for d in datasets:
            lines.append(f"| {d['dataset_id']} | {d.get('sample_count', 'n/a')} |")
    else:
        lines.append(f"_{b1}_")
    lines.append("")

    # ---- Section 6 --------------------------------------------------------
    lines.append("## 6. Gap Analysis: In GDC but Not in Xena\n")
    b2 = xena.get("b2_gap_analysis") or {}
    if isinstance(b2, dict) and "error" not in b2:
        gaps = b2.get("gaps_in_xena") or []
        covered = b2.get("xena_covered") or []
        lines.append("### Covered by Xena")
        for c in covered:
            lines.append(f"- {c}")
        lines.append("\n### Gaps (GDC data categories with no Xena equivalent)")
        if gaps:
            for g in gaps:
                lines.append(f"- {g}")
        else:
            lines.append("_No gaps detected._")
    else:
        lines.append(f"_{b2}_")
    lines.append("")

    # ---- Section 7 --------------------------------------------------------
    lines.append("## 7. cBioPortal: All Per-Study Endpoints Discovered\n")
    c1 = cbio.get("c1_per_study_endpoints") or {}
    c3 = cbio.get("c3_endpoint_results") or {}
    lines.append(f"Total per-study endpoints: **{c1.get('per_study_count', 0)}**\n")

    group1 = (c3.get("group1_queried") or []) if isinstance(c3, dict) else []
    group2 = (c3.get("group2_skipped") or []) if isinstance(c3, dict) else []

    lines.append("### Table A — Endpoints Queried (studyId only required)\n")
    lines.append("| Endpoint Path | Item Count | Status |")
    lines.append("|---------------|------------|--------|")
    for ep in group1:
        lines.append(f"| `{ep['path']}` | {ep['item_count']} | {ep['status']} |")
    lines.append("")

    lines.append("### Table B — Endpoints Skipped (additional required parameters)\n")
    lines.append("| Endpoint Path | Required Additional Parameters | Parameter Types |")
    lines.append("|---------------|-------------------------------|-----------------|")
    for ep in group2:
        params = ep.get("extra_required_params") or []
        pnames = ", ".join(p["name"] for p in params)
        ptypes = ", ".join(p["type"] for p in params)
        lines.append(f"| `{ep['path']}` | {pnames} | {ptypes} |")
    lines.append("")

    # ---- Section 8 --------------------------------------------------------
    lines.append("## 8. METABRIC: Results from Each Queried Endpoint\n")
    for ep in group1:
        lines.append(f"### `{ep['path']}`\n")
        lines.append(f"- **Status:** {ep['status']}")
        lines.append(f"- **Item count:** {ep['item_count']}")
        resp = ep.get("response")
        if resp is not None:
            if isinstance(resp, list) and len(resp) > 0:
                lines.append(f"- **First item preview:** `{json.dumps(resp[0])[:300]}`")
            elif isinstance(resp, dict):
                lines.append(f"- **Response preview:** `{json.dumps(resp)[:300]}`")
        lines.append("")

    # ---- Section 9 --------------------------------------------------------
    lines.append("## 9. Recommended Download Strategy\n")
    lines.append("_[To be filled after reviewing this report]_\n")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"[{ts()}] Starting discovery...")

    report = {
        "discovery_timestamp": ts(),
        "gdc_tcga_brca": {},
        "xena_gdc_tcga_brca": {},
        "cbio_metabric": {},
    }

    print(f"[{ts()}] Task A: GDC TCGA-BRCA...")
    report["gdc_tcga_brca"] = task_a_gdc()

    print(f"[{ts()}] Task B: UCSC Xena GDC Hub...")
    report["xena_gdc_tcga_brca"] = task_b_xena(report["gdc_tcga_brca"])

    print(f"[{ts()}] Task C: cBioPortal METABRIC...")
    report["cbio_metabric"] = task_c_cbio()

    print(f"[{ts()}] Writing outputs...")

    # JSON
    with open(JSON_OUT, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, default=str)

    # Markdown
    md = render_markdown(report)
    with open(MD_OUT, "w", encoding="utf-8") as f:
        f.write(md)

    print("Discovery complete.")
    print(f"JSON: {JSON_OUT}")
    print(f"Markdown: {MD_OUT}")


if __name__ == "__main__":
    main()
