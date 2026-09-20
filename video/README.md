# LinkedIn showcase video

A 44-second, 1080×1080 video for Biomodelling.jl 2.0, built with [Remotion](https://remotion.dev).
Every chart in it is rendered from the simulation outputs in `paper/output/`, so the
video shows the same numbers as the white paper.

## Build

```bash
npm install
python scripts/make_charts.py     # regenerate public/chart_*.png from paper/output
npm run studio                             # preview in the browser
npm run render                             # write out/biomodelling.mp4
```

Rendering needs a Chrome headless shell. On a machine without one, pass an existing
binary:

```bash
npx remotion render Showcase out/biomodelling-showcase.mp4 \
  --browser-executable=/path/to/chrome-headless-shell --codec=h264 --crf=20
```

## Layout

| File | What it holds |
|---|---|
| `src/Showcase.tsx` | the nine scenes and their timings |
| `src/theme.tsx` | palette, Inter faces, reveal/scene/chart primitives |
| `src/Cells.tsx` | the two SVG animations (drug response, cell cycle) |
| `public/chart_*.png` | charts rendered from real simulation output |
