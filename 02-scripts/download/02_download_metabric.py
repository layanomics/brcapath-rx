"""
Download METABRIC clinical data for TCGA-BRCA Phase 1 via the cBioPortal REST API.

The S3 tarball URL is no longer available (403). This script uses the
cBioPortal API directly, which is confirmed working.

Endpoints used (all confirmed in discovery report):
  /studies/brca_metabric/clinical-attributes                          → 36 attribute definitions
  /studies/brca_metabric/patients                                     → 2,509 patient records
  /studies/brca_metabric/samples                                      → 2,509 sample records
  /studies/brca_metabric/clinical-data                                → sample-level values (long format)
  /studies/brca_metabric/clinical-data?clinicalDataType=PATIENT       → patient-level values (long format)

Outputs:
  01-data/raw/metabric/clinical/metabric_clinical_attributes.json
  01-data/raw/metabric/clinical/metabric_patients.json
  01-data/raw/metabric/clinical/metabric_samples.json
  01-data/raw/metabric/clinical/metabric_clinical_data.json          (sample-level)
  01-data/raw/metabric/clinical/metabric_clinical_data_patient.json  (patient-level)
  01-data/raw/metabric/clinical/metabric_clinical_wide.tsv           (merged wide, all 36 attrs)
  01-data/audit/metabric/metabric_manifest.json

Genomic endpoints (expression, mutations, CNA) are deferred to the genomic phase
and are NOT downloaded here.

Usage:
  python 02_download_metabric.py            # dry-run: show what will be called
  python 02_download_metabric.py --confirm  # run downloads and audit
"""

import argparse
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

# ── Configuration ─────────────────────────────────────────────────────────────
API_BASE   = "https://www.cbioportal.org/api"
STUDY_ID   = "brca_metabric"
CLINICAL_DIR  = Path("01-data/raw/metabric/clinical")
AUDIT_DIR     = Path("01-data/audit/metabric")
MANIFEST_PATH = AUDIT_DIR / "metabric_manifest.json"

PAGE_SIZE       = 10_000   # max records per page request
REQUEST_TIMEOUT = 60       # seconds
RETRY_WAIT      = 5        # seconds between retries
MAX_RETRIES     = 3

# ── Endpoints — Phase 1 (clinical only) ──────────────────────────────────────
ENDPOINTS: list[dict] = [
    {
        "key":       "clinical_attributes",
        "path":      f"/studies/{STUDY_ID}/clinical-attributes",
        "filename":  "metabric_clinical_attributes.json",
        "expected":  36,
        "paginated": False,
        "desc":      "Clinical attribute definitions (column metadata)",
    },
    {
        "key":       "patients",
        "path":      f"/studies/{STUDY_ID}/patients",
        "filename":  "metabric_patients.json",
        "expected":  2509,
        "paginated": True,
        "desc":      "Patient list",
    },
    {
        "key":       "samples",
        "path":      f"/studies/{STUDY_ID}/samples",
        "filename":  "metabric_samples.json",
        "expected":  2509,
        "paginated": True,
        "desc":      "Sample list",
    },
    {
        "key":       "clinical_data",
        "path":      f"/studies/{STUDY_ID}/clinical-data",
        "filename":  "metabric_clinical_data.json",
        "expected":  27868,
        "paginated": True,
        "desc":      "Clinical attribute values — long format, sample-level",
    },
    {
        "key":         "clinical_data_patient",
        "path":        f"/studies/{STUDY_ID}/clinical-data",
        "filename":    "metabric_clinical_data_patient.json",
        "expected":    60216,   # ~24 patient-level attrs × 2,509 patients
        "paginated":   True,
        "extra_params": {"clinicalDataType": "PATIENT"},
        "desc":        "Clinical attribute values — long format, patient-level only",
    },
]

DEFERRED_NOTE = (
    "Genomic endpoints (expression, mutations, CNA) are deferred to the genomic phase "
    "and are intentionally NOT downloaded here."
)

# ── Logging ───────────────────────────────────────────────────────────────────
AUDIT_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(AUDIT_DIR / "metabric_download.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)


# ── API helpers ───────────────────────────────────────────────────────────────
def api_get(path: str, params: dict | None = None, retries: int = MAX_RETRIES) -> list | dict:
    """GET from cBioPortal API with retry logic. Returns parsed JSON."""
    url = f"{API_BASE}{path}"
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT,
                                headers={"Accept": "application/json"})
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            log.warning("[attempt %d/%d] GET %s — %s", attempt, retries, url, exc)
            if attempt < retries:
                time.sleep(RETRY_WAIT * attempt)
    raise RuntimeError(f"All {retries} attempts failed for GET {url}")


def fetch_all_pages(
    path: str,
    expected: int | None = None,
    extra_params: dict | None = None,
) -> list:
    """
    Paginate a cBioPortal list endpoint using pageNumber / pageSize params.
    Returns the full list of records.
    """
    all_records: list = []
    page = 0
    while True:
        params = {"pageSize": PAGE_SIZE, "pageNumber": page}
        if extra_params:
            params.update(extra_params)
        batch = api_get(path, params=params)
        if not isinstance(batch, list):
            raise ValueError(f"Expected list response from {path}, got {type(batch)}")
        all_records.extend(batch)
        log.info("  page %d: %d records (total so far: %d)", page, len(batch), len(all_records))
        if len(batch) < PAGE_SIZE:
            break
        page += 1

    if expected is not None and len(all_records) != expected:
        log.warning(
            "Record count mismatch for %s: expected %d, got %d",
            path, expected, len(all_records),
        )
    return all_records


# ── Step 1: Download one endpoint ────────────────────────────────────────────
def download_endpoint(ep: dict, out_dir: Path) -> tuple[list, Path]:
    """
    Download one endpoint, save JSON, return (records, dest_path).
    Skips if the file already exists with content.
    """
    dest = out_dir / ep["filename"]

    if dest.exists() and dest.stat().st_size > 0:
        log.info("[skip-exists] %s  (%d bytes)", dest.name, dest.stat().st_size)
        with open(dest) as fh:
            return json.load(fh), dest

    log.info("=== Downloading: %s ===", ep["desc"])
    log.info("  endpoint: GET %s%s", API_BASE, ep["path"])

    if ep["paginated"]:
        records = fetch_all_pages(
            ep["path"],
            expected=ep["expected"],
            extra_params=ep.get("extra_params"),
        )
    else:
        records = api_get(ep["path"], params=ep.get("extra_params"))
        if not isinstance(records, list):
            records = [records]
        log.info("  fetched %d records", len(records))

    out_dir.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as fh:
        json.dump(records, fh, indent=2)
    log.info("  saved: %s  (%d records, %.1f KB)", dest.name, len(records), dest.stat().st_size / 1024)
    return records, dest


# ── Step 2: Pivot and merge clinical data long → wide ────────────────────────
def _pivot_to_wide(records: list, label: str) -> pd.DataFrame:
    """Pivot one set of long-format clinical records to wide (patientId × attrs)."""
    id_col   = "patientId"
    attr_col = "clinicalAttributeId"
    val_col  = "value"

    df_long = pd.DataFrame(records)
    log.info("  %s long-format shape: %s", label, df_long.shape)

    missing = {id_col, attr_col, val_col} - set(df_long.columns)
    if missing:
        raise ValueError(
            f"{label}: missing expected columns {missing}. Found: {list(df_long.columns)}"
        )

    df_wide = df_long.pivot_table(
        index=id_col,
        columns=attr_col,
        values=val_col,
        aggfunc="first",
    ).reset_index()
    df_wide.columns.name = None
    return df_wide


def build_wide_tsv(
    sample_data: list,
    patient_data: list,
    out_dir: Path,
) -> tuple[pd.DataFrame, Path]:
    """
    Pivot sample-level and patient-level clinical data to wide format separately,
    then merge on patientId.  Always rebuilds — never skips.
    Saves as metabric_clinical_wide.tsv.  Returns (DataFrame, dest_path).
    """
    log.info("=== Building merged wide TSV (sample + patient level) ===")
    dest = out_dir / "metabric_clinical_wide.tsv"

    df_sample  = _pivot_to_wide(sample_data,  "sample-level")
    df_patient = _pivot_to_wide(patient_data, "patient-level")

    # Merge: all patients from sample side; patient attrs joined in
    df_wide = df_sample.merge(df_patient, on="patientId", how="left", suffixes=("", "_pat"))

    df_wide.to_csv(dest, sep="\t", index=False)
    log.info(
        "  merged wide-format shape: %d patients × %d columns",
        df_wide.shape[0], df_wide.shape[1],
    )
    log.info("  saved: %s", dest.name)
    return df_wide, dest


# ── Step 3: Audit wide TSV ───────────────────────────────────────────────────
def audit_wide(df: pd.DataFrame) -> dict:
    missing_pct = (df.isnull().mean() * 100).round(2).to_dict()
    return {
        "rows":         int(df.shape[0]),
        "cols":         int(df.shape[1]),
        "columns":      list(df.columns),
        "dtypes":       {c: str(t) for c, t in df.dtypes.items()},
        "missing_pct":  missing_pct,
        "mean_missing": round(df.isnull().mean().mean() * 100, 2),
    }


# ── Step 4: Save manifest ─────────────────────────────────────────────────────
def save_manifest(
    download_ts: str,
    endpoint_results: list[dict],
    wide_audit: dict,
    wide_path: Path,
) -> None:
    manifest = {
        "download_timestamp": download_ts,
        "study_id":           STUDY_ID,
        "api_base":           API_BASE,
        "note":               DEFERRED_NOTE,
        "phase1_files": [
            {
                "endpoint":    f"{API_BASE}{ep['path']}",
                "filename":    ep["filename"],
                "description": ep["desc"],
                "expected":    ep["expected"],
                "actual":      result["actual"],
                "local_path":  result["local_path"],
                "size_bytes":  result["size_bytes"],
            }
            for ep, result in endpoint_results
        ],
        "wide_tsv": {
            "local_path": str(wide_path),
            "size_bytes": wide_path.stat().st_size if wide_path.exists() else None,
            **wide_audit,
        },
        "deferred_files": {
            "status": "deferred - genomic phase",
            "endpoints_not_downloaded": [
                f"{API_BASE}/studies/{STUDY_ID}/molecular-profiles",
                f"{API_BASE}/studies/{STUDY_ID}/expression (mrna)",
                f"{API_BASE}/studies/{STUDY_ID}/mutations",
                f"{API_BASE}/studies/{STUDY_ID}/copy-number-alterations",
            ],
        },
    }
    with open(MANIFEST_PATH, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log.info("Manifest saved: %s", MANIFEST_PATH)


# ── Step 5: Print summary ─────────────────────────────────────────────────────
def print_summary(endpoint_results: list, wide_audit: dict, wide_path: Path) -> None:
    print("\n" + "=" * 80)
    print("PHASE 1 — Clinical files downloaded and audited")
    print("-" * 80)
    print(f"  {'Endpoint':<45}  {'Expected':>8}  {'Actual':>8}  FILE")
    print("-" * 80)
    for ep, result in endpoint_results:
        match = "OK" if result["actual"] == ep["expected"] else "MISMATCH"
        print(f"  {ep['path']:<45}  {ep['expected']:>8}  {result['actual']:>8}  "
              f"{ep['filename']}  [{match}]")

    print()
    print(f"  Wide TSV: {wide_path.name}")
    print(f"    Patients : {wide_audit['rows']}")
    print(f"    Columns  : {wide_audit['cols']}")
    print(f"    Mean missing: {wide_audit['mean_missing']}%")

    print("\n" + "=" * 80)
    print("DEFERRED — Genomic endpoints (not downloaded in Phase 1)")
    print("-" * 80)
    print("  expression / mutations / CNA → deferred - genomic phase")
    print("=" * 80)


# ── Dry-run preview ───────────────────────────────────────────────────────────
def print_dry_run_preview() -> None:
    print("\n" + "=" * 80)
    print("DRY RUN — METABRIC cBioPortal API download plan")
    print("=" * 80)
    print(f"\n  API base: {API_BASE}")
    print(f"  Study ID: {STUDY_ID}")
    print()
    print(f"  {'ENDPOINT':<50}  {'EXPECTED':>8}  OUTPUT FILE")
    print("-" * 80)
    for ep in ENDPOINTS:
        print(f"  {ep['path']:<50}  {ep['expected']:>8}  {ep['filename']}")
    print()
    print("  Additionally produces:")
    print("    metabric_clinical_wide.tsv  — sample + patient data merged (patients × all 36 attrs)")
    print()
    print("  DEFERRED (not downloaded here):")
    print("    expression, mutations, CNA  → deferred - genomic phase")
    print()
    print(f"  NOTE: {DEFERRED_NOTE}")
    print()
    print("  Output directory : " + str(CLINICAL_DIR))
    print("  Manifest path    : " + str(MANIFEST_PATH))
    print("\n" + "=" * 80)
    print("Dry-run complete. Re-run with --confirm to start downloading.")


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download METABRIC clinical data via cBioPortal API."
    )
    parser.add_argument(
        "--confirm", action="store_true",
        help="Run the actual downloads and audit. Without this, only preview the plan.",
    )
    args = parser.parse_args()

    if not args.confirm:
        print_dry_run_preview()
        return

    download_ts = datetime.now(timezone.utc).isoformat()
    CLINICAL_DIR.mkdir(parents=True, exist_ok=True)

    try:
        # ── Step 1: Download all Phase 1 endpoints
        endpoint_results: list[tuple[dict, dict]] = []
        clinical_data: list        = []
        clinical_data_patient: list = []

        for ep in ENDPOINTS:
            records, dest = download_endpoint(ep, CLINICAL_DIR)
            result = {
                "actual":     len(records),
                "local_path": str(dest),
                "size_bytes": dest.stat().st_size if dest.exists() else None,
            }
            endpoint_results.append((ep, result))
            if ep["key"] == "clinical_data":
                clinical_data = records
            elif ep["key"] == "clinical_data_patient":
                clinical_data_patient = records

        # ── Step 2: Pivot + merge both datasets to wide TSV
        df_wide, wide_path = build_wide_tsv(clinical_data, clinical_data_patient, CLINICAL_DIR)

        # ── Step 3: Audit the wide TSV
        wide_audit = audit_wide(df_wide)

        # ── Step 4: Save manifest
        save_manifest(download_ts, endpoint_results, wide_audit, wide_path)

        # ── Step 5: Print summary
        print_summary(endpoint_results, wide_audit, wide_path)

    except Exception as exc:
        log.error("Fatal error: %s", exc, exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
