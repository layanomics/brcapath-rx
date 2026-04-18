# IPTW Confounding Demo — METABRIC BRCA

Interactive single-page demo showing the move from a raw survival paradox to a validated IPTW-corrected analysis pipeline.

---

## Local development

```bash
cd demos/iptw-confounding-demo
npm install
npm run dev
```

Opens at `http://localhost:5173`

## Build for deployment

```bash
npm run build
```

Output goes to `dist/`. Verify locally with:

```bash
npm run preview
```

---

## GitHub Pages deployment

### Option A — Deploy `dist/` folder to a dedicated repo

1. Build: `npm run build`
2. Create (or use) a GitHub Pages repo, e.g. `your-username/brcapath-rx-demos`
3. Copy the entire `dist/` folder contents into the repo root (or a subfolder like `iptw/`)
4. Commit and push
5. In repo Settings → Pages → set source branch to `main` and folder to `/` (or `/iptw`)

### Option B — Deploy from this repo using gh-pages

```bash
npm install --save-dev gh-pages
```

Add to `package.json` scripts:
```json
"deploy": "gh-pages -d dist"
```

Then:
```bash
npm run build
npm run deploy
```

### Important

- The `.nojekyll` file in the repo root tells GitHub Pages not to process the site with Jekyll.
  Copy it into your `dist/` output or the Pages repo root.
- `vite.config.js` uses `base: './'` — all asset paths are relative, so the site works
  in a subdirectory without any path changes.

---

## Section map → meeting story

| Section | Tells the audience |
|---|---|
| **Hero** | "We have a working IPTW pipeline — here's the update." |
| **01 The Paradox** | "The raw numbers looked broken — treated patients appeared to survive less." |
| **Tab: Raw Comparison** | "This is confounding by indication — the groups were not comparable to begin with." |
| **Tab: IPTW Workflow** | "Here's the correction method and the real numbers from our dataset." |
| **Tab: After Weighting** | "The weights are stable and the balance check passed." |
| **03 Balance Assessment** | "Toggle to see just how bad the imbalance was, and how completely it was resolved." |
| **04 Takeaway** | "Three sentences: the paradox is explained, IPTW is built in, the next output will be adjusted." |
| **Next Steps box** | "Here's what the next result will look like." |
