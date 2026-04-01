from __future__ import annotations

import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests


ROOT = Path(__file__).resolve().parent
SOURCE_DIR = ROOT / "01-sources"
RAW_DIR = ROOT / "02-raw-downloads"

REQUEST_TIMEOUT = 120
PAGE_SIZE = 10_000

TCGA_XENA_FILES = [
    {
        "dataset_key": "tcga_brca_xena_gdc_clinical",
        "cohort": "TCGA-BRCA",
        "source_platform": "UCSC Xena GDC hub",
        "source_url": "https://gdc.xenahubs.net/download/TCGA-BRCA.clinical.tsv.gz",
        "local_name": "tcga_brca_xena_gdc_clinical.tsv.gz",
        "role": "primary",
        "used_for": "TCGA ready-made clinical descriptors and observational treatment availability",
        "notes": "Primary TCGA clinical table selected for tonight's Stage 1 non-XML workflow.",
    },
    {
        "dataset_key": "tcga_brca_xena_gdc_survival",
        "cohort": "TCGA-BRCA",
        "source_platform": "UCSC Xena GDC hub",
        "source_url": "https://gdc.xenahubs.net/download/TCGA-BRCA.survival.tsv.gz",
        "local_name": "tcga_brca_xena_gdc_survival.tsv.gz",
        "role": "primary",
        "used_for": "TCGA overall survival endpoint",
        "notes": "Primary TCGA survival table selected for tonight's Stage 1 workflow.",
    },
    {
        "dataset_key": "tcga_brca_xena_tcgahub_clinical_matrix",
        "cohort": "TCGA-BRCA",
        "source_platform": "UCSC Xena TCGA hub",
        "source_url": "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
        "local_name": "tcga_brca_xena_tcgahub_clinical_matrix.tsv",
        "role": "auxiliary",
        "used_for": "TCGA subtype and receptor augmentation only",
        "notes": "Used only because PAM50 and receptor fields are absent from the GDC-hub TCGA clinical table.",
    },
]

CBIO_BASE = "https://www.cbioportal.org/api"
METABRIC_ENDPOINTS = [
    {
        "dataset_key": "metabric_cbio_clinical_attributes",
        "cohort": "METABRIC",
        "source_platform": "cBioPortal REST API",
        "endpoint": "/studies/brca_metabric/clinical-attributes",
        "params": {},
        "local_name": "metabric_cbio_clinical_attributes.json",
        "role": "primary",
        "used_for": "Clinical field metadata",
        "notes": "Clinical attribute definitions for METABRIC.",
    },
    {
        "dataset_key": "metabric_cbio_clinical_data_sample",
        "cohort": "METABRIC",
        "source_platform": "cBioPortal REST API",
        "endpoint": "/studies/brca_metabric/clinical-data",
        "params": {},
        "local_name": "metabric_cbio_clinical_data_sample.json",
        "role": "primary",
        "used_for": "Sample-level clinical values including stage and receptor fields",
        "notes": "Sample-level clinical long table.",
    },
    {
        "dataset_key": "metabric_cbio_clinical_data_patient",
        "cohort": "METABRIC",
        "source_platform": "cBioPortal REST API",
        "endpoint": "/studies/brca_metabric/clinical-data",
        "params": {"clinicalDataType": "PATIENT"},
        "local_name": "metabric_cbio_clinical_data_patient.json",
        "role": "primary",
        "used_for": "Patient-level clinical values including survival, subtype, and treatment fields",
        "notes": "Patient-level clinical long table.",
    },
]


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def ensure_dirs() -> None:
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_delimited_rows(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        line_count = sum(1 for _ in handle)
    return max(line_count - 1, 0)


def download_binary(url: str, destination: Path) -> None:
    with requests.get(url, stream=True, timeout=REQUEST_TIMEOUT) as response:
        response.raise_for_status()
        with destination.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)


def fetch_json(endpoint: str, params: dict[str, Any]) -> list[dict[str, Any]]:
    response = requests.get(
        f"{CBIO_BASE}{endpoint}",
        params=params,
        timeout=REQUEST_TIMEOUT,
        headers={"Accept": "application/json"},
    )
    response.raise_for_status()
    payload = response.json()
    if isinstance(payload, list):
        return payload
    return [payload]


def fetch_paginated_json(endpoint: str, params: dict[str, Any]) -> list[dict[str, Any]]:
    all_rows: list[dict[str, Any]] = []
    page_number = 0
    while True:
        page_params = {"pageSize": PAGE_SIZE, "pageNumber": page_number}
        page_params.update(params)
        rows = fetch_json(endpoint, page_params)
        all_rows.extend(rows)
        if len(rows) < PAGE_SIZE:
            break
        page_number += 1
    return all_rows


def write_json(destination: Path, payload: list[dict[str, Any]]) -> None:
    with destination.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def write_manifest(rows: list[dict[str, Any]]) -> None:
    headers = [
        "dataset_key",
        "cohort",
        "source_platform",
        "role",
        "used_for",
        "source_url",
        "local_path",
        "retrieved_at_utc",
        "size_bytes",
        "sha256",
        "record_count",
        "notes",
    ]
    out_path = SOURCE_DIR / "dataset_manifest.tsv"
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(header, "")) for header in headers) + "\n")


def write_source_notes(manifest_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Source Notes",
        "",
        f"Generated: {utc_now_iso()}",
        "",
        "## Scope",
        "",
        "- Fast Stage 1 prognostic survival demo only.",
        "- No treatment-effect estimation.",
        "- No treatment recommendation claim.",
        "- No image or WSI data.",
        "- GDC XML is not the main path in this folder.",
        "",
        "## Primary Sources Used Tonight",
        "",
        "1. TCGA-BRCA from UCSC Xena GDC hub clinical and survival tables.",
        "2. METABRIC from cBioPortal clinical endpoints.",
        "",
        "## Auxiliary Source Used",
        "",
        "- UCSC Xena TCGA hub `TCGA.BRCA.sampleMap/BRCA_clinicalMatrix`.",
        "- Reason: the GDC-hub TCGA clinical table does not expose PAM50/receptor fields needed for subtype summaries.",
        "- Restriction: this auxiliary table is used only to augment TCGA subtype and receptor annotations.",
        "",
        "## Downloaded Files",
        "",
    ]

    for row in manifest_rows:
        lines.extend(
            [
                f"### {row['dataset_key']}",
                "",
                f"- Cohort: {row['cohort']}",
                f"- Platform: {row['source_platform']}",
                f"- Role: {row['role']}",
                f"- Purpose: {row['used_for']}",
                f"- URL: {row['source_url']}",
                f"- Local file: `{row['local_path']}`",
                f"- Retrieved at: {row['retrieved_at_utc']}",
                f"- Records: {row['record_count']}",
                f"- SHA256: `{row['sha256']}`",
                f"- Notes: {row['notes']}",
                "",
            ]
        )

    lines.extend(
        [
            "## Meeting-Safe Framing",
            "",
            "- Treatment fields in these public cohorts are observational descriptors only.",
            "- Treatment fields are retained for availability summaries and context, not for causal claims.",
            "- All prognostic modeling in this trial is cohort-specific and non-causal.",
            "",
        ]
    )

    (SOURCE_DIR / "source_notes.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    manifest_rows: list[dict[str, Any]] = []

    for dataset in TCGA_XENA_FILES:
        destination = RAW_DIR / dataset["local_name"]
        if not destination.exists():
            download_binary(dataset["source_url"], destination)

        manifest_rows.append(
            {
                "dataset_key": dataset["dataset_key"],
                "cohort": dataset["cohort"],
                "source_platform": dataset["source_platform"],
                "role": dataset["role"],
                "used_for": dataset["used_for"],
                "source_url": dataset["source_url"],
                "local_path": str(destination.relative_to(ROOT)).replace("\\", "/"),
                "retrieved_at_utc": utc_now_iso(),
                "size_bytes": destination.stat().st_size,
                "sha256": sha256_file(destination),
                "record_count": count_delimited_rows(destination),
                "notes": dataset["notes"],
            }
        )

    for dataset in METABRIC_ENDPOINTS:
        destination = RAW_DIR / dataset["local_name"]
        if destination.exists():
            with destination.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
        else:
            if dataset["dataset_key"] == "metabric_cbio_clinical_attributes":
                payload = fetch_json(dataset["endpoint"], dataset["params"])
            else:
                payload = fetch_paginated_json(dataset["endpoint"], dataset["params"])
            write_json(destination, payload)

        source_url = f"{CBIO_BASE}{dataset['endpoint']}"
        if dataset["params"]:
            query_string = "&".join(f"{key}={value}" for key, value in dataset["params"].items())
            source_url = f"{source_url}?{query_string}"

        manifest_rows.append(
            {
                "dataset_key": dataset["dataset_key"],
                "cohort": dataset["cohort"],
                "source_platform": dataset["source_platform"],
                "role": dataset["role"],
                "used_for": dataset["used_for"],
                "source_url": source_url,
                "local_path": str(destination.relative_to(ROOT)).replace("\\", "/"),
                "retrieved_at_utc": utc_now_iso(),
                "size_bytes": destination.stat().st_size,
                "sha256": sha256_file(destination),
                "record_count": len(payload),
                "notes": dataset["notes"],
            }
        )

    write_manifest(manifest_rows)
    write_source_notes(manifest_rows)

    print("Downloaded source files and wrote source metadata:")
    print(f"  - {SOURCE_DIR / 'dataset_manifest.tsv'}")
    print(f"  - {SOURCE_DIR / 'source_notes.md'}")


if __name__ == "__main__":
    main()
