# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/); this project uses semantic versioning.

## [2.0.0] — 2026-11-02

Major release: completes Phase 3 Cross-Species Functional Conservation across 5 kingdoms/models,
overhauls the web portal to flagship CoSE styling with an embedded PDF reader, and completely
revamps Figure 8 into an integrated spatial-molecular systems biology architecture.

### Added
- **Figure 8 Spatial-Molecular Systems Architecture**:
  - Re-engineered `Figures/kegg_pathway_network.*` into a 4-panel systems biology model:
    - **Panel A**: Subcellular mechanotransduction & biochemical reaction cascade across cell wall (pH ~4.8 acid growth via `AHA1/2`, `XTH22`, `EXPA1`), plasma membrane (`MSL10`/`MCA1`, `PIN2/3/7`, `AUX1`), statoliths (`CPN60A`), cytosol (`TCH1/CAM2`, `RABA1A`, `ARF1A1C`), nucleus, mitochondria (`ATPA`), and ER/26S proteasome (`HSP70-1`, `HSP90.1`, `PAA1`).
    - **Panel B**: Single-cell tissue activation landscape mapping pathways against the 183 single-cell clusters, highlighting peak cell-wall loosening in inflorescence stem vascular cambium (`stem_9`, $p_{\mathrm{adj}} = 0.003$) and $\text{Ca}^{2+}$ signaling in root tip columella statocytes.
    - **Panel C**: Curated pathway cross-talk network replacing cryptic 5-digit KEGG IDs with 6 annotated biological hubs.
    - **Panel D**: Chronological systems model tracing the microgravity perception cascade.
  - Python generator script: `Code/generate_figure8_revamped.py` producing 300 DPI publication PNGs and valid vector SVGs in both print (light) and web (obsidian dark) palettes.
- **Phase 3 Cross-Species Functional Conservation**:
  - Multi-species orthology mapping engine supporting *Arabidopsis*, *Homo sapiens*, *Mus musculus*, *Drosophila melanogaster*, and *Caenorhabditis elegans*.
  - Precomputed static data bundles in `docs/assets/ortho/` (`orthomap.json`, `conserved_degs.json`, `function_enrichment.json`, `demo_cross_species.json`).
  - Interactive client-side conservation dashboard in Section 05 featuring an UpSet conservation breakdown, 5-species Pfam domain/function enrichment heatmap, and searchable gene-domain drill-down explorer.
- **Flagship CoSE Theme & Portal Upgrades**:
  - Obsidian dark palette with smooth OS auto-detection and manual toggle.
  - Interactive embedded PDF reader in Section 02 with fullscreen mode and direct draft download.
  - Figure gallery lightbox with full-resolution image links.
  - Left navigation document map and ecosystem integration with companion CoSE decoders.
- **Updated Manuscript**:
  - Recompiled Springer Nature `sn-jnl` manuscript (`Manuscript/latex/main.pdf`, 13 pages) and synchronized with `docs/Tropism_autodecoder_2026.pdf`.

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
- **Phase 2 activated (full-atlas):** the auto-decoder was retrained on all 432,919 nuclei (GSE226097, Apple MPS),
  bulk OSDR RNA-seq (14 studies) was re-acquired + re-run through DESeq2 → NNLS deconvolution → elastic-net classifier
  (AUC 0.915 ± 0.043), and the exported bundle in `docs/assets/model/` (`manifest.phase2_ready: true`) now powers the
  live Phase 2 method — real cell-type deconvolution, stimulus projection, and Flight-vs-GC probability in the browser.

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
