# Tropism Autodecoder — web tool (`docs/`)

Client-side companion tool for the *Arabidopsis thaliana* Spaceflight Tropism Recognition System
(`Tropism_autodecoder_2026`). Upload an expression matrix (`.csv`/`.txt`), get a publication-ready
figure panel + auto-generated legend decoding four tropism signatures. **Runs entirely in the browser —
no server, no data leaves the visitor's machine.**

## Enable GitHub Pages

1. Copy this `docs/` folder (and `WEB_TOOL_PLAN.md`) into the root of the `Tropism_autodecoder_2026` repository and push.
2. GitHub → **Settings → Pages → Build and deployment → Deploy from a branch → `main` / `/docs` → Save**.
3. Live at: `https://dr-richard-barker.github.io/Tropism_autodecoder_2026/`
4. Add that URL to the manuscript ("Data & Code Availability") and the top-level `README.md`.

## Files

| File | Purpose |
|---|---|
| `index.html` | Landing page + embedded tool (single shareable URL) |
| `assets/style.css` | Accessible white/black theme, red–white–blue (RdBu) data scale |
| `assets/app.js` | Parse → singscore → composed SVG figure → legend → SVG/PNG/CSV export |
| `assets/signatures.js` | Curated tropism signature gene sets + symbol→AGI aliases (replaceable) |
| `assets/sample_data.csv` | **Synthetic** demo matrix for the "Load example" button (not experimental data) |
| `.nojekyll` | Serve `assets/` verbatim |

## Local preview

```bash
cd docs && python3 -m http.server 8791
# open http://localhost:8791
```
(Serve over HTTP, not `file://`, so the "Load example" fetch works.)

## Method (Phase 1)

Rank-based single-sample enrichment (**singscore**; Foroutan et al., *BMC Bioinformatics* 2018) of the uploaded
profile against curated, literature-verified tropism marker sets. Scale-invariant (counts, TPM/FPKM, or normalized
all work). These are marker-enrichment scores, **not** the manuscript's auto-decoder / deconvolution output —
that parity is Phase 2 (see `WEB_TOOL_PLAN.md`).

## Extending / replacing signatures

Edit `assets/signatures.js`. Each entry is `{ key, name, subtitle, genes:[AGI…], members }`. To swap in signatures
exported from the pipeline, replace the array — no changes to `app.js` needed.

## Verification

The scoring/figure core was validated by running the real `signatures.js` + `app.js` under
JavaScriptCore against the bundled demo data: correct dominant-tropism calls, transposed-matrix
auto-detection, gene-symbol / `.1`-suffix matching, and a well-formed multi-panel SVG. The
singscore math was independently cross-checked in Python (identical results).
