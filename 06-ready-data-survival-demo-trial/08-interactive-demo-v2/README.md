# BRCAPath-Rx Stage 1 Interactive Demo — v2

Static, GitHub Pages-ready interactive survival demo built from the current verified Stage 1 package.

## Files

| File | Purpose |
|---|---|
| `index.html` | Single-file static web app |
| `demo_predictions.json` | Exported Cox model payload (copied from `08-interactive-demo/`) |
| `README.md` | This file |

## What is real

- All C-index values, coefficients, baseline survival curves, and risk thresholds come from the current Stage 1 package outputs
- TCGA subtype field: `PAM50Call_RNAseq`
- METABRIC subtype field: `PAM50_DERIVED_GENEFU_PAM50` (derived via genefu, not the original CLAUDIN_SUBTYPE)
- TCGA C-index: 0.729 · METABRIC C-index: 0.643

## How predictions work

No backend is required. The site:
1. Loads `demo_predictions.json` via `fetch()`
2. Reads the exported Cox coefficients, centering means, baseline survival, and risk thresholds
3. Computes the linear predictor and survival curve live in the browser using Chart.js

Per-profile predictions are **live in browser**, not precomputed.

## How to serve locally

```bash
cd D:\Projects\brcapath-rx\06-ready-data-survival-demo-trial\08-interactive-demo-v2
python -m http.server 8080
```

Then open: http://localhost:8080

The site will not work when opened as a local file (`file://`) because `fetch()` is blocked by browser CORS policy for local files.

## GitHub Pages deployment

Drop the contents of this folder into any GitHub Pages branch root or `/docs` folder.
No build step, no Node.js, no dependencies beyond the CDN links in index.html.

## Meeting-safe use

Say aloud during the meeting:

> "This is a cohort-specific Stage 1 prognostic demo. It uses a Cox proportional hazards model
> to rank relative survival risk within two independent breast cancer cohorts using age, clinical,
> and molecular subtype features. It does not estimate treatment effects, does not adjust for
> treatment confounding, and does not recommend treatment."

Additional points to be prepared to say:
- TCGA and METABRIC use different model inputs because the available data and population differ
- METABRIC subtype now uses derived PAM50 from genefu, not the original CLAUDIN_SUBTYPE field
- Treatment fields exist in both cohorts but are observational and incomplete — treatment confounding is a later-stage problem

## Differences from v1 (08-interactive-demo)

| Aspect | v1 (08-interactive-demo) | v2 (this folder) |
|---|---|---|
| CSS | Custom CSS variables | Tailwind CSS (CDN) |
| Chart library | Plotly.js 2.35 | Chart.js 4.4 |
| C-index tab | Brief text + two values | Visual concordant-pairs example + interpretation |
| How It Works tab | Two short paragraphs | 4-step flow + Stage 1 scope table |
| Data & Limits tab | JSON-loaded list items | Structured cohort cards + treatment confounding block |
| Subtype hints | None | Plain-language PAM50 descriptions per subtype |
| JSON payload | Same | Same (identical file, identical science) |
| Model math | Verified correct | Identical to v1 |
