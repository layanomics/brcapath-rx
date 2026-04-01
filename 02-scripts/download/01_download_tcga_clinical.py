"""
Download Clinical Supplement and Biospecimen Supplement files for TCGA-BRCA from GDC.

Outputs:
  01-data/raw/tcga-brca/gdc/clinical/   — downloaded files
  01-data/audit/tcga-brca/clinical_download_manifest.json

Usage:
  python 01_download_tcga_clinical.py           # preview inventory, no download
  python 01_download_tcga_clinical.py --confirm  # run the actual download
"""

import argparse
import json
import hashlib
import logging
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import requests

# ── Configuration ─────────────────────────────────────────────────────────────
GDC_BASE       = "https://api.gdc.cancer.gov"
PROJECT_ID     = "TCGA-BRCA"
PROJECT_FIELD  = "cases.project.project_id"
OUT_DIR        = Path("01-data/raw/tcga-brca/gdc/clinical")
MANIFEST_DIR   = Path("01-data/audit/tcga-brca")
MANIFEST_PATH  = MANIFEST_DIR / "clinical_download_manifest.json"
TARGET_CATEGORIES  = {"clinical", "biospecimen"}
TARGET_DATA_TYPES  = {"Clinical Supplement", "Biospecimen Supplement"}
PAGE_SIZE      = 500
REQUEST_TIMEOUT = 120  # seconds per request
RETRY_WAIT     = 5     # seconds between retries
MAX_RETRIES    = 3

# ── Skip-reason lookup for known non-Phase-1 data types ──────────────────────
SKIP_REASONS: dict[str, str] = {
    "Slide Image":                   "imaging phase, not Phase 1",
    "Pathology Report":              "requires NLP processing, deferred",
    "Intermediate Analysis Archive": "pipeline artifacts, not needed",
}

# ── Logging ───────────────────────────────────────────────────────────────────
MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(MANIFEST_DIR / "tcga_clinical_download.log", mode="w"),
    ],
)
log = logging.getLogger(__name__)


# ── Helpers ───────────────────────────────────────────────────────────────────
def gdc_post(endpoint: str, payload: dict, retries: int = MAX_RETRIES) -> dict:
    """POST to GDC API with retry logic."""
    url = f"{GDC_BASE}/{endpoint}"
    for attempt in range(1, retries + 1):
        try:
            resp = requests.post(url, json=payload, timeout=REQUEST_TIMEOUT)
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            log.warning("[attempt %d/%d] %s — %s", attempt, retries, url, exc)
            if attempt < retries:
                time.sleep(RETRY_WAIT * attempt)
    raise RuntimeError(f"All {retries} attempts failed for {url}")


def gdc_get_download(file_id: str, dest: Path, retries: int = MAX_RETRIES) -> None:
    """Stream-download a single GDC file by file_id."""
    url = f"{GDC_BASE}/data/{file_id}"
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=REQUEST_TIMEOUT) as r:
                r.raise_for_status()
                dest.parent.mkdir(parents=True, exist_ok=True)
                with open(dest, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        fh.write(chunk)
            return
        except requests.RequestException as exc:
            log.warning("[attempt %d/%d] download %s — %s", attempt, retries, file_id, exc)
            if dest.exists():
                dest.unlink()
            if attempt < retries:
                time.sleep(RETRY_WAIT * attempt)
    raise RuntimeError(f"All {retries} download attempts failed for file_id={file_id}")


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ── Step 0 — Clean up any SVS files already downloaded ───────────────────────
def cleanup_svs_files() -> tuple[int, int]:
    """Delete any .svs files in OUT_DIR. Returns (count_deleted, bytes_freed)."""
    if not OUT_DIR.exists():
        return 0, 0
    svs_files = list(OUT_DIR.glob("*.svs"))
    total_bytes = 0
    for f in svs_files:
        size = f.stat().st_size
        total_bytes += size
        f.unlink()
        log.info("[cleanup-svs] Deleted: %s  (%.2f MB)", f.name, size / 1e6)
    return len(svs_files), total_bytes


# ── Step 1 — Discover available data_category / data_type values ──────────────
def discover_inventory() -> dict:
    """
    Query /files with facets on data_category and data_type.
    Returns a nested dict: {data_category: {data_type: count}}.
    """
    log.info("=== Step 1: Discovering file inventory via facets ===")
    payload = {
        "filters": {
            "op": "=",
            "content": {"field": PROJECT_FIELD, "value": PROJECT_ID},
        },
        "facets": "data_category,data_type,access",
        "size": 0,
        "format": "JSON",
    }
    data = gdc_post("files", payload)
    aggs = data.get("data", {}).get("aggregations", {})

    categories = {
        b["key"]: b["doc_count"]
        for b in aggs.get("data_category", {}).get("buckets", [])
    }
    types = {
        b["key"]: b["doc_count"]
        for b in aggs.get("data_type", {}).get("buckets", [])
    }
    access = {
        b["key"]: b["doc_count"]
        for b in aggs.get("access", {}).get("buckets", [])
    }

    log.info("data_category breakdown:")
    for cat, cnt in sorted(categories.items(), key=lambda x: -x[1]):
        marker = " ← TARGET" if cat.lower() in TARGET_CATEGORIES else ""
        log.info("  %-35s %5d%s", cat, cnt, marker)

    log.info("data_type breakdown:")
    for dt, cnt in sorted(types.items(), key=lambda x: -x[1]):
        log.info("  %-45s %5d", dt, cnt)

    log.info("access breakdown:")
    for ac, cnt in sorted(access.items(), key=lambda x: -x[1]):
        log.info("  %-15s %5d", ac, cnt)

    return {"categories": categories, "types": types, "access": access}


# ── Step 2/3 — Fetch file metadata for a given category ──────────────────────
def fetch_file_list(category: str) -> list[dict]:
    """
    Return metadata for all open-access files in the given data_category.
    Fields: file_id, file_name, data_category, data_type, file_size, md5sum, state.
    """
    log.info("Fetching file list for category: %s", category)
    fields = [
        "file_id", "file_name", "data_category", "data_type",
        "file_size", "md5sum", "state", "access",
    ]
    payload = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": PROJECT_FIELD, "value": PROJECT_ID}},
                {"op": "=", "content": {"field": "data_category", "value": category}},
                {"op": "=", "content": {"field": "access", "value": "open"}},
            ],
        },
        "fields": ",".join(fields),
        "size": PAGE_SIZE,
        "from": 0,
        "format": "JSON",
    }

    all_hits = []
    while True:
        data = gdc_post("files", payload)
        hits = data["data"]["hits"]
        all_hits.extend(hits)
        total = data["data"]["pagination"]["total"]
        fetched = payload["from"] + len(hits)
        log.info("  %d / %d", fetched, total)
        if fetched >= total:
            break
        payload["from"] = fetched

    log.info("  → %d open-access files in '%s'", len(all_hits), category)
    return all_hits


# ── Filter: keep only TARGET_DATA_TYPES, log everything else ─────────────────
def filter_and_log_skipped(all_files: list[dict]) -> list[dict]:
    """
    Return only files whose data_type is in TARGET_DATA_TYPES.
    Log each skipped data_type with file count and reason.
    """
    target: list[dict] = []
    skipped_list: list[dict] = []

    for f in all_files:
        dt = f.get("data_type", "unknown")
        if dt in TARGET_DATA_TYPES:
            target.append(f)
        else:
            skipped_list.append(f)

    if skipped_list:
        counts = Counter(f.get("data_type", "unknown") for f in skipped_list)
        log.info("=== Skipped data types (not in Phase 1 scope) ===")
        for dt, count in sorted(counts.items(), key=lambda x: -x[1]):
            reason = SKIP_REASONS.get(dt, "not in Phase 1 scope")
            log.info("  [SKIP] data_type=%-40s  files=%4d  reason=%s", dt, count, reason)
    else:
        log.info("No files skipped — all fetched files match target data types.")

    return target


# ── Inventory printout ────────────────────────────────────────────────────────
def print_inventory(target_files: list[dict]) -> None:
    counts = Counter(f.get("data_type", "unknown") for f in target_files)
    total_bytes = sum(f.get("file_size") or 0 for f in target_files)

    print("\n" + "=" * 60)
    print("PHASE 1 DOWNLOAD INVENTORY — TCGA-BRCA")
    print("=" * 60)
    print(f"  {'Data Type':<35}  {'Files':>6}")
    print("-" * 60)
    for dt in sorted(TARGET_DATA_TYPES):
        cnt = counts.get(dt, 0)
        marker = " ← target" if cnt > 0 else " ← NOT FOUND"
        print(f"  {dt:<35}  {cnt:>6}{marker}")
    print("-" * 60)
    print(f"  {'TOTAL':<35}  {len(target_files):>6}  ({round(total_bytes / 1e6, 1)} MB estimated)")
    print("=" * 60)


# ── Step 4 — Validate downloaded file ─────────────────────────────────────────
def validate_file(path: Path, expected_size: int | None, expected_md5: str | None) -> dict:
    result = {"size_ok": None, "md5_ok": None, "issues": []}
    if not path.exists():
        result["issues"].append("file missing")
        return result

    actual_size = path.stat().st_size
    if expected_size is not None:
        result["size_ok"] = actual_size == expected_size
        if not result["size_ok"]:
            result["issues"].append(
                f"size mismatch: expected {expected_size}, got {actual_size}"
            )

    if expected_md5:
        actual_md5 = md5_file(path)
        result["md5_ok"] = actual_md5 == expected_md5
        if not result["md5_ok"]:
            result["issues"].append(
                f"md5 mismatch: expected {expected_md5}, got {actual_md5}"
            )

    return result


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    parser = argparse.ArgumentParser(description="Download TCGA-BRCA clinical files from GDC.")
    parser.add_argument(
        "--confirm", action="store_true",
        help="Run the actual download. Without this flag, only preview the inventory.",
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Step 0: Clean up any SVS files already in the output directory
    svs_count, svs_bytes = cleanup_svs_files()
    if svs_count:
        log.info(
            "[cleanup] Removed %d SVS file(s) — freed %.2f MB",
            svs_count, svs_bytes / 1e6,
        )
    else:
        log.info("[cleanup] No SVS files found in %s", OUT_DIR)

    # ── Step 1: Discover inventory
    inventory = discover_inventory()

    available_categories = {k.lower() for k in inventory["categories"]}
    download_categories = TARGET_CATEGORIES & available_categories
    missing = TARGET_CATEGORIES - available_categories
    if missing:
        log.warning("Categories not found in GDC for this project: %s", missing)
    log.info("Will query categories: %s", sorted(download_categories))

    # ── Steps 2 & 3: Collect all file metadata for target categories
    all_files: list[dict] = []
    for cat in sorted(download_categories):
        all_files.extend(fetch_file_list(cat))

    log.info("Total open-access files in target categories: %d", len(all_files))

    # ── Filter to Phase 1 data types only
    target_files = filter_and_log_skipped(all_files)
    log.info("Files remaining after Phase 1 filter: %d", len(target_files))

    # ── Print inventory
    print_inventory(target_files)

    if not args.confirm:
        print("\nDry-run complete. Re-run with --confirm to start downloading.")
        return

    # ── Download + validate
    manifest_records = []
    attempted = succeeded = failed = 0
    total_bytes = 0

    for meta in target_files:
        file_id   = meta.get("file_id", "")
        file_name = meta.get("file_name", file_id)
        category  = meta.get("data_category", "unknown")
        data_type = meta.get("data_type", "unknown")
        file_size = meta.get("file_size")
        md5sum    = meta.get("md5sum", "")

        dest = OUT_DIR / file_name
        attempted += 1

        try:
            if dest.exists() and file_size and dest.stat().st_size == file_size:
                log.info("[skip-exists] %s", file_name)
            else:
                log.info("[download] %s  (%s bytes)", file_name, file_size)
                gdc_get_download(file_id, dest)

            # ── Step 4: Validate
            validation = validate_file(dest, file_size, md5sum)
            if validation["issues"]:
                for issue in validation["issues"]:
                    log.warning("VALIDATION ISSUE %s: %s", file_name, issue)
                failed += 1
            else:
                succeeded += 1
                total_bytes += dest.stat().st_size if dest.exists() else 0

        except Exception as exc:
            log.error("[FAILED] %s — %s", file_name, exc)
            failed += 1
            validation = {"issues": [str(exc)]}

        # ── Step 5: Build manifest record
        manifest_records.append({
            "file_id":            file_id,
            "file_name":          file_name,
            "data_category":      category,
            "data_type":          data_type,
            "file_size":          file_size,
            "md5sum":             md5sum,
            "download_timestamp": datetime.now(timezone.utc).isoformat(),
            "local_path":         str(dest),
            "validation":         validation,
        })

    # ── Step 5: Save manifest
    manifest = {
        "project":    PROJECT_ID,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "attempted":  attempted,
            "succeeded":  succeeded,
            "failed":     failed,
            "total_mb":   round(total_bytes / 1e6, 2),
        },
        "files": manifest_records,
    }
    with open(MANIFEST_PATH, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log.info("Manifest saved: %s", MANIFEST_PATH)

    # ── Final summary
    print("\n" + "=" * 60)
    print(f"files attempted : {attempted}")
    print(f"files succeeded : {succeeded}")
    print(f"files failed    : {failed}")
    print(f"total size MB   : {round(total_bytes / 1e6, 2)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
