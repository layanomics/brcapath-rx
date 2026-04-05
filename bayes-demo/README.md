# Bayes Theorem in Breast Cancer Subtype Prediction

**BRCAPath-Rx — Educational Demo**

A single-page static website explaining Bayesian inference through breast cancer subtype prediction, using real data from the METABRIC cohort.

## Files

| File | Purpose |
|------|---------|
| `index.html` | Single-page demo |
| `styles.css` | Academic dark-blue/pink styling |
| `script.js` | Interactive posterior calculator |
| `data.json` | Real-data-derived priors, likelihoods, precomputed posteriors |

## Data provenance

All values are **real-data-derived** from:
- **Source**: `01-data/processed/metabric/metabric_clinical_clean.tsv`
- **n_total**: 2,509 METABRIC patients
- **n_with_subtype**: 1,974 patients (CLAUDIN_SUBTYPE field)
- **Priors**: subtype frequency counts from METABRIC
- **Likelihoods**: per-subtype conditional rates for ER, PR, HER2, tumor grade
- **Posteriors**: naive Bayes computation over all 16 feature combinations

## Deploy to GitHub Pages

### Option A — Deploy from this subfolder

1. Create a new GitHub repo (e.g. `brcapath-bayes-demo`).
2. Copy the four files (`index.html`, `styles.css`, `script.js`, `data.json`) to the repo root.
3. Push to `main`.
4. In repo Settings → Pages → Source: `main` branch, `/ (root)`.
5. Your site will be live at `https://<username>.github.io/brcapath-bayes-demo/`.

### Option B — Deploy from a subfolder of the main repo

1. Commit this `08-bayes-demo/` folder to the `brcapath-rx` repo.
2. In repo Settings → Pages → Source: `main`, folder `/08-bayes-demo`.
3. Live at `https://<username>.github.io/brcapath-rx/08-bayes-demo/`.

### Local preview

Because `data.json` is loaded via `fetch()`, open with a local server:

```bash
# Python
cd 08-bayes-demo
python -m http.server 8080
# then open http://localhost:8080
```

Or use the VS Code Live Server extension.

## Methodological note

The interactive calculator uses the **naive Bayes assumption**: ER, PR, HER2, and tumor grade are treated as conditionally independent given subtype. This is a standard approximation for educational purposes. The demo is not a clinical decision tool.
