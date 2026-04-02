# Stage 2 Confounding Website V2

This is the tutorial-first version of the BRCAPath-Rx Stage 2 confounding website.

## What is different from the original website

- It teaches confounding in a story order rather than starting with technical project outputs.
- It begins with everyday examples, then a fictional patient story, and only then moves to the real METABRIC case study.
- It explains both overlap weighting and propensity score matching in plain language.
- It explicitly connects confounding control to the future BRCAPath-Rx treatment AI logic.

## Files

- `index.html`: the tutorial-first static website
- `confounding_demo_payload.json`: content and numbers used by the website
- `assets/`: copies of the real project figures used in the case-study section

## Preview locally

From `07-stage2-confounding-prototype`:

```powershell
python -m http.server 8000
```

Then open:

`http://localhost:8000/04-website-v2/`

## GitHub Pages

This folder is GitHub Pages-ready:

- no build step
- relative paths only
- static HTML/CSS/JS
- all real case-study figures stored locally in `assets/`
