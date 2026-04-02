# Confounding Demo Website v4

This folder contains the `v4` static demo website for the Stage 2 confounding prototype.

## Purpose

`v4` is designed as a click-through tutorial for a nontechnical supervisor. It teaches:

- the basic idea of confounding
- breast-cancer-like treated-vs-untreated examples
- one real completed causal case study
- simple adjustment concepts
- how this supports the future BRCAPath-Rx treatment AI

## What Is Real vs Illustrative

Real completed case:

- ER-positive METABRIC
- hormone therapy YES vs NO
- overall survival
- naive comparison plus overlap-weighted adjustment

Illustrations only:

- the fictional breast-cancer-like cards
- the confounding simulator
- the AI workflow sketch
- the future candidate treatment panels

## File Layout

- `index.html`: static entry page
- `confounding_demo_payload.json`: content payload for text and labels
- `assets/`: real figure PNGs copied into the published site folder

## Serving Locally

From `07-stage2-confounding-prototype`:

```powershell
python -m http.server 8000
```

Then open:

```text
http://127.0.0.1:8000/04-website-v4/
```

## Publish Safety

- Uses relative paths only
- Figure PNGs are referenced directly from `assets/`
- No build step required
- GitHub Pages compatible as a static folder
