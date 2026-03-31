# Discovery Scripts

00_discover_sources.py
  Queries GDC, UCSC Xena, and cBioPortal to produce a full
  inventory of what data exists. Downloads nothing.
  Run time: approximately 1-2 minutes.
  Output: 01-data/audit/source_discovery_report.md

  Note on output location: Reports are saved to 01-data/audit/
  which is excluded from Git. To preserve a report permanently,
  copy it here into 02-scripts/discovery/ after reviewing it.
