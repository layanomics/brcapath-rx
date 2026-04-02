# Stage 2 Confounding Website V3

This is the shorter, more interactive confounding tutorial website.

## What is new in v3

- Direct static image references in HTML using `assets/...`
- Shorter sections with less wall-of-text
- Interactive toggles for:
  - wrong conclusion vs better interpretation
  - fictional breast-cancer-like scenarios
  - before adjustment vs after adjustment in the real case
- Interactive slider showing how imbalance alone can distort a crude treated-vs-untreated comparison
- A clearer split between:
  - the one real completed case study
  - candidate treatment questions for future work

## Files

- `index.html`
- `confounding_demo_payload.json`
- `assets/`

## Preview locally

From `07-stage2-confounding-prototype`:

```powershell
python -m http.server 8000
```

Then open:

`http://localhost:8000/04-website-v3/`

## GitHub Pages note

This folder is GitHub Pages compatible:

- no build step
- relative paths only
- image assets stored locally in `assets/`
- figures referenced directly from HTML instead of being injected later from payload paths
