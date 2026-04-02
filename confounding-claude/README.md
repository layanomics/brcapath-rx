# BRCAPath-Rx · Stage 2 Confounding Demo (Claude version)

Interactive static website demonstrating confounding adjustment in breast cancer treatment comparisons.

## How to view

Open `index.html` directly in any modern browser — no build step required.

For GitHub Pages: push this folder to a repository and enable Pages from the root or `docs/` folder.

## Files

```
04-website-claude/
├── index.html                      ← Main website (single file)
├── confounding_demo_payload.json   ← All real analysis data as JSON
├── README.md                       ← This file
└── assets/
    ├── fig_naive_km.png            ← Unadjusted Kaplan-Meier curves
    ├── fig_balance_before.png      ← Covariate balance before adjustment
    ├── fig_balance_after.png       ← Covariate balance after overlap weighting
    ├── fig_adjusted_survival.png   ← Overlap-weighted survival curves
    └── fig_ps_overlap.png          ← Propensity score overlap density plot
```

## Six sections

| # | Section | Key interaction |
|---|---------|----------------|
| 1 | Quick Idea | Analogy cards |
| 2 | Naive Comparison | Confounding simulator with 3 toggles |
| 3 | Three Cases | Case A (real), B & C (future) tab selector + 4-view real case |
| 4 | Methods | OW vs PSM switcher |
| 5 | AI Model | Without/with adjustment toggle |
| 6 | Done vs Next | Progress tracker |

## What is real vs illustrative

| Content | Status |
|---------|--------|
| Case A — hormone therapy · ER+ METABRIC | **REAL completed analysis** |
| Naive HR = 1.50, Adjusted HR = 1.00 | **Real computed values** |
| 18 → 0 imbalances, max |SMD| 0.696 → 0.002 | **Real computed values** |
| 5 figures (balance, KM, PS overlap) | **Real generated figures** |
| Case B — chemotherapy | Future candidate only |
| Case C — radiotherapy | Future candidate only |
| Confounding simulator | Illustrative toy model |
| AI pipeline diagram | Conceptual illustration |

## Technical

- Static HTML/CSS/JS — no frameworks, no build step
- All interactivity via vanilla JavaScript
- GitHub Pages compatible — all paths are relative
- Figures referenced as `./assets/filename.png`
- Image fallback text shown if assets are missing
