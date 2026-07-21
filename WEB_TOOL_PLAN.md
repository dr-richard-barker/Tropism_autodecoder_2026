# Tropism Autodecoder — Web Tool & GitHub Pages Plan

**Companion to:** `Tropism_autodecoder_2026` (Barker, 2026) — *Arabidopsis thaliana Spaceflight Tropism Recognition System*
**Purpose of this document:** define and record what we are building so the manuscript can cite a live, working URL where reviewers and readers upload their own *Arabidopsis* transcriptional data (`.csv`/`.txt`) and receive a publication-ready figure panel and matching figure legend that decodes tropism signatures.
**Status:** Phase 1 implemented in this repository's `docs/` folder (client-side, GitHub Pages). Phases 2–3 scoped below.
**Scope note:** This plan is *new* and complementary to the existing `PLAN.md`, which covers the compute pipeline only (it does not describe any web deployment). Nothing here changes the pipeline; it exposes a subset of its logic in the browser.

---

## 1. Goal & success criteria

A single public URL (GitHub Pages) that a peer reviewer can open with no install, upload an *Arabidopsis* expression matrix, and get back:

1. A **multi-panel, publication-ready figure** (SVG + PNG) decoding the four tropism signatures — gravitropism, phototropism, thigmotropism (touch), hydrotropism — as defined in the manuscript.
2. An **auto-generated figure legend** (Methods-style caption) that matches the figure and states exactly what was computed, how many genes matched, and the caveats.
3. **Downloadable results** (scores as CSV, figure as SVG/PNG, legend as Markdown).

**Acceptance criteria**
- [x] Works fully client-side (no server, no data leaves the browser) — required for GitHub Pages and for reviewer data privacy.
- [x] Accepts `.csv` and `.tsv/.txt`, auto-detects delimiter and matrix orientation, handles AGI locus IDs (`AT1G70940`) and common gene symbols.
- [x] "Load example" button produces a working figure in one click (bundled synthetic demo matrix, clearly labelled synthetic).
- [x] White background, black text, red–white–blue (RdBu) diverging data scale; WCAG-AA contrast; keyboard navigable.
- [x] Every number on screen is traceable to a transparent, documented method — no hidden or fabricated statistics.
- [ ] (Phase 2) Full parity with the manuscript auto-decoder via exported model artifacts.

---

## 2. Hard constraint that shapes the design: GitHub Pages is static

GitHub Pages serves static files only — there is no Python/R runtime. The manuscript pipeline (variational auto-decoder, DESeq2, NNLS deconvolution onto a 4000-HVG × 183-cluster atlas signature matrix, elastic-net classifier) cannot execute in a Pages site. Two honest options exist and we use both, staged:

- **Phase 1 (now):** a transparent, standard, training-free **rank-based single-sample enrichment** (singscore; Foroutan et al., 2018) of the uploaded profile against **curated, literature-referenced tropism signature gene sets**. This runs anywhere, on any single matrix, and every step is inspectable.
- **Phase 2 (upgrade path):** ship **exported model artifacts** from the real pipeline (signature matrix, stimulus codes, classifier coefficients) as static JSON/CSV and re-implement the pipeline's *inference* math (NNLS + linear projection + linear classifier) in JavaScript. This gives numerical parity with the manuscript without any server.

We do **not** claim Phase-1 marker scores are the auto-decoder output. The UI and legend state which method produced each number.

---

## 3. What the repository already provides (grounding)

Confirmed by inspecting the live repo (`main`):

| Asset | Use to the web tool |
|---|---|
| `README.md`, `Manuscript/manuscript.md` | Framing, four tropism definitions, headline metrics (AUC 0.919; 16,002 sig. genes; stem padj=0.003) shown on the landing page. |
| `Results/flight_vs_gc_feature_importance.csv` | Real elastic-net feature coefficients (cell-type clusters + `stim_dim_*`). Basis for the Phase-2 classifier and for the "what the model keys on" explainer. |
| `Results/celltype_flight_vs_ground.csv` | Real per-cell-type flight-vs-ground differences (stem_9, seed_425d_7, silique_20 significant). Basis for Phase-2 tissue projection. |
| `Results/meta_analysis_results.tsv` (26,402 genes) | Source of the real **flight-response gene signature** to add as a fifth panel in Phase 2 (export top BH-significant genes as a static list). |
| `metadata.json` | DataCite citation block (title, creator Richard Barker, CC-BY-4.0, DOI placeholder). |

The web tool reuses these rather than inventing content, consistent with the project's "ground all docs in real content" rule.

---

## 4. Phase 1 — the tool we ship now (client-side)

### 4.1 Input contract
- File: `.csv`, `.tsv`, or `.txt`. Delimiter auto-detected (`,`, `\t`, or `;`).
- Shape: genes × samples (default) **or** samples × genes (auto-detected by which axis matches known gene IDs). First column/row = gene IDs; remaining = numeric expression (raw counts, TPM, FPKM, or normalized — the rank-based method is scale-invariant, so no pre-normalization is required).
- Gene IDs: AGI locus IDs (`AT5G57090`, case-insensitive, `.1` transcript suffix stripped) or symbols (`PIN2`) resolved via a bundled alias table.

### 4.2 Method (transparent)
1. Parse → numeric matrix + gene IDs + sample names.
2. Per sample, rank all genes by expression (average ranks for ties).
3. For each tropism signature set S, compute **singscore**: mean normalized rank of matched genes, rescaled against the theoretical min/max mean rank for |S| genes, then centered to **[−1, +1]** (positive = coordinately high expression of the signature; negative = low).
4. Report a samples × 4-tropism score table; per sample, call the **dominant tropism** (argmax) with a confidence = gap to the runner-up.
5. Report signature **coverage** (matched/total genes) per set; warn when coverage is low (e.g., hydrotropism is intentionally a small set — see §4.4).

### 4.3 Output figure (single composed SVG, RdBu)
- **Panel A — Heatmap:** samples (rows) × tropisms (columns), RdBu diverging, colorbar, value labels for accessibility.
- **Panel B — Bar chart:** mean signature score per tropism across all samples, with per-sample spread.
- **Panel C — Coverage / call strip:** genes matched per signature and the dominant-tropism call per sample.
- Composed as one downloadable figure (SVG vector + PNG raster). RdBu red = high, blue = low, white ≈ 0.

### 4.4 Curated signature sets (Phase 1) — all loci literature-verified
Marker genes are canonical tropism regulators; AGI loci verified against TAIR/primary literature before hard-coding (no invented identifiers).

- **Gravitropism:** PIN2/EIR1 (AT5G57090), PIN3 (AT1G70940), PIN7 (AT1G23080), AUX1 (AT2G38120), LAZY1 (AT5G14090).
- **Phototropism:** PHOT1/NPH1 (AT3G45780), PHOT2 (AT5G58140), NPH3 (AT5G64330), RPT2 (AT2G30520), PKS1 (AT2G02950).
- **Thigmotropism (touch):** TCH1/CAM2 (AT5G37780), TCH2/CML24 (AT5G37770), TCH3/CML12 (AT2G41100), TCH4/XTH22 (AT5G57560), MSL10 (AT5G12080).
- **Hydrotropism:** MIZ1 (AT2G41660), GNOM/MIZ2 (AT1G13980). *Deliberately small — mirrors the manuscript's stated hydrotropism data limitation; flagged low-coverage in the UI.*

Signatures live in `docs/assets/signatures.js` with a documented schema so they can be edited/extended or **replaced wholesale by exported pipeline signatures** without touching the app code.

### 4.5 Auto figure legend
Generated from the run: figure title, per-panel description, exact signature membership, genes matched (X/Y per set), method sentence (singscore + citation), color-scale statement, tool version + date, and caveats (rank-based; not the auto-decoder; low-coverage sets). Downloadable as Markdown.

### 4.6 Non-goals for Phase 1
No cell-type deconvolution, no 32-dim stimulus projection, no flight/ground classifier (all require exported artifacts — Phase 2). No cross-species input (Phase 3). Stated plainly in-app.

---

## 5. Phase 2 — parity with the manuscript pipeline (artifact export)

To reproduce the manuscript's *inference* in-browser, export these static artifacts from the pipeline into `docs/assets/model/`:

| Artifact | Shape | Produced by | Enables |
|---|---|---|---|
| `signature_matrix.csv` | 4000 HVG × 183 clusters | atlas preprocessing | NNLS cell-type deconvolution of an uploaded bulk profile |
| `stimulus_codes.csv` | 183 clusters × 32 dims | `train_autodecoder.py` | stimulus-activation projection (weighted avg of per-cluster codes) |
| `classifier_coefficients.json` | features → weights + intercept | `meta_classifier.py` | Flight-vs-GC score (AUC 0.919 model) reproduced as a dot product |
| `flight_signature_genes.json` | top BH-sig genes + log2FC | `meta_classifier.py` meta-analysis | fifth "spaceflight response" panel |
| `gene_index.json` | HVG AGI order | atlas preprocessing | align uploaded genes to the signature matrix |

Browser re-implements: NNLS (active-set, small solver in JS) → fractions → stimulus projection → linear classifier. Pure linear algebra; fully static. The app already isolates the scoring module so this drops in behind a "Method: curated markers ▸ full auto-decoder" toggle.

**Reviewer-facing benefit:** the tool then returns the *same* cell-type fractions, stimulus scores, and flight-vs-ground probability that the paper reports, on the reviewer's own data.

---

## 6. Phase 3 — orthology network (future, per the request)

Allow non-*Arabidopsis* transcriptomes by mapping input genes to *Arabidopsis* AGI loci through an ortholog table (e.g., PLAZA / OrthoFinder / Ensembl Plants). Ship `docs/assets/ortho/<species>.json` (source gene → AGI, with confidence), let the user pick species, map, then run the same scoring. Report mapping coverage and flag one-to-many/low-confidence orthologs. This is the "eventually add an orthology network" item; deferred, schema sketched, not implemented in Phase 1.

---

## 7. Accessibility & theme spec

- Background `#ffffff`, body text `#000000`; links/accents from the RdBu ends — blue `#2166ac`, red `#b2182b`.
- Data color scale: ColorBrewer **RdBu** diverging (red high / white mid / blue low) — a standard, reasonably color-vision-deficiency-tolerant diverging palette; reinforced with on-cell numeric labels so meaning never depends on color alone.
- WCAG-AA contrast; visible focus rings; full keyboard operation; ARIA labels on controls and the figure; responsive down to mobile; system font stack; no external fonts/CDNs (works offline and within any CSP).

---

## 8. File layout (this deliverable)

```
docs/                      # GitHub Pages source (set Pages → Deploy from branch → /docs)
├── index.html             # Landing + embedded tool (one shareable URL)
├── assets/
│   ├── style.css          # Accessible white/black + RdBu theme
│   ├── app.js             # Parse → singscore → SVG panel → legend → export
│   ├── signatures.js      # Curated tropism signatures + alias table (replaceable)
│   └── sample_data.csv    # SYNTHETIC demo matrix (clearly labelled)
├── .nojekyll              # Serve assets verbatim
└── README.md              # How to enable Pages + how to extend
WEB_TOOL_PLAN.md           # This document
```

## 9. Deployment

1. Copy `docs/` and `WEB_TOOL_PLAN.md` into the repository root of `Tropism_autodecoder_2026`.
2. GitHub → **Settings → Pages → Build and deployment → Deploy from a branch → `main` / `/docs`**.
3. Live URL: `https://dr-richard-barker.github.io/Tropism_autodecoder_2026/`.
4. Add that URL to the manuscript ("Data and Code Availability") and to `README.md`.

## 10. Testing & acceptance

- Load-example path renders a full figure + legend with zero uploads.
- Round-trip a real matrix (genes × samples and the transposed orientation) → both parse.
- AGI-with-version (`AT5G57090.1`) and symbol (`PIN2`) inputs both match.
- Empty/low-coverage signature → warning, not a crash.
- Exports (SVG/PNG/CSV/MD) open correctly and the legend matches the figure.
- Keyboard-only walkthrough; contrast check; mobile width.

## 11. Honesty & limitations (carried into the UI and legend)

- Phase-1 scores are **rank-based marker enrichment**, not the auto-decoder; stated in-app and in every legend.
- Curated signatures are small and literature-based; hydrotropism especially is low-coverage (consistent with the manuscript's own limitation).
- The tool decodes **relative signature activity within the uploaded samples**; it is not a diagnostic and does not assign spaceflight status until Phase 2 ships the real classifier.
- Synthetic example data is labelled synthetic and must not be presented as experimental results.
