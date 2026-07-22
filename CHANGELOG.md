# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/); this project uses semantic versioning.

## [1.1.0] — 2026-07-21

Tooling release: adds an interactive web tool and a GPU / full-atlas training workflow.
**The published auto-decoder model is unchanged** — the full 432,919-nucleus retrain is prepared
but has not yet been run, so all v1.0.0 results (including Flight-vs-GC AUC = 0.919) still stand.

### Added
- **Interactive web tool** (`docs/`, GitHub Pages) — upload an *Arabidopsis* expression matrix
  (`.csv`/`.txt`) and get a publication-ready figure panel + auto-generated legend decoding four
  tropism signatures (gravitropism, phototropism, thigmotropism, hydrotropism). Fully client-side;
  uploaded data never leaves the browser. Accessible white/black theme with a red–white–blue
  (RdBu) data scale. Phase-1 method: rank-based single-sample enrichment (singscore) against
  curated, literature-verified marker sets.
- **`WEB_TOOL_PLAN.md`** — full design, accessibility spec, and roadmap (Phase 2 pipeline parity;
  Phase 3 orthology network for non-*Arabidopsis* input).
- **GPU / full-atlas training workflow:**
  - `Code/train_autodecoder_gpu.py` — GPU-enabled, memory-scalable retrain. Same model and loss as
    `train_autodecoder.py`; adds device/precision auto-selection (bf16-mixed), memmap streaming of
    the full atlas, a startup cardinality guard (183 clusters / 12 organs / 10 stages), and a CLI.
  - `Code/export_atlas.R` — atlas-extraction step (full atlas or reproducible stratified subsample)
    emitting the exact binary contract the trainer consumes; reconstructed from the manuscript
    Methods and byte-verified against the reader.
  - `environment.gpu.yml`, `Dockerfile.gpu` — CUDA 12.1 environment and container.
  - `GPU_TRAINING_PLAN.md` — plan, compute/memory estimates, run recipe, and acceptance criteria.
- **`CHANGELOG.md`** — this file.
- **Web tool Phase 2 plumbing** — `docs/assets/model.js` reproduces the pipeline's inference in-browser
  (NNLS deconvolution → stimulus projection → elastic-net Flight-vs-GC classifier); `Code/export_web_artifacts.py`
  exports the model bundle to `docs/assets/model/`; `meta_classifier.py` now also dumps `classifier_params.json`
  (coefficients + intercept + scaler). A method toggle activates Phase 2 automatically once artifacts are present;
  it ships gated (`manifest.phase2_ready: false`) so the live tool stays on Phase 1 until they're exported.

### Changed
- **README**: added a "Web Tool" section; listed the new scripts and files in the repository
  structure; reworded Limitation #1 to describe the v1.0.0 CPU / 60,792-cell stratified subsample
  and the now-available (planned) full-atlas GPU retraining path.
- **`metadata.json`**: `version` 1.0.0 → 1.1.0; added an `Updated` date (2026-07-21).

### Notes
- Limitation #1 remains open until the full-atlas GPU retrain is executed; update it (and confirm
  AUC ≥ 0.919) once results exist.

## [1.0.0] — 2026-07-20
- Initial release: end-to-end FAIR pipeline integrating 1,337 bulk transcriptomic samples with the
  Salk *Arabidopsis* Developmental Atlas foundation model; auto-decoder deconvolution, DESeq2
  differential expression, DerSimonian-Laird meta-analysis, elastic-net classification
  (Flight-vs-GC AUC = 0.919), and ggPlantMap / KEGG visualizations.
