# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/); this project uses semantic versioning.

## [Unreleased]

Corrections from an evidence check of the repository (28 Oct 2026). Re-running the pipeline is still to do.

### Fixed
- **Meta-analysis τ²** (`Code/meta_classifier.py`): the DerSimonian–Laird between-study variance was divided by
  Σw too many times, so it came out too small by a factor of Σw and random-effects p-values were far too small.
  The fixed function matches `metafor::rma(method = "DL")` on test inputs. `Results/meta_analysis_results.tsv`,
  the volcano plot and the "16,002 significant genes" still come from the old code (README, Limitation 7).
- **GSE8300 removed** from `Data/harmonized_metadata.tsv`: GEO lists it as a *Haemophilus influenzae* genomic
  hybridization series, not Arabidopsis. The dataset is now 1,253 samples (63 from 5 GEO series); counts in the
  README, manuscripts and `metadata.json` follow. It was not among the 398 deconvolved samples.
- **Citations** (`Manuscript/latex/references.bib`, `Manuscript/main.tex`, `Manuscript/manuscript.md`), checked
  against CrossRef: Toyota & Gilroy is 2013, *Am J Bot* 100:111–125; CIBERSORTx is Newman et al. 2019,
  *Nat Biotechnol* 37:773–782; PhysioSpace is Lenz et al. 2013, *PLoS ONE* 8:e77627 (the "Schiller 2021"
  entry matched no paper); the ggPlantmap title and pages are corrected (*J Exp Bot* 75:5366–5376).
- **`Manuscript/main.tex`**: GSE111555 (a human exercise study) replaced by GSE115554; spaceflight / ground
  control counts and the top classifier features now match the data.
- **`metadata.json`**: the identifier `10.5281/zenodo.arabidopsis-tropism` resolved nowhere; it is now the
  placeholder `10.5281/zenodo.XXXXXXX` until a Zenodo DOI is minted.
- `seed_425d_7` flight difference is +0.04 % (`Results/celltype_flight_vs_ground.csv`), not +0.05 %.
- `Figures/dataset_summary` subtitle is computed from the metadata instead of typed in.

### Marked TODO
- AUC 0.919 ± 0.047, 85 % accuracy, tropism F1 1.000 and the auto-decoder metrics have no committed output
  (README, Limitation 8). The ggPlantMap maps add random jitter per region (README, Limitation 9).

## [1.1.0] — 2026-07-21

Tooling release: adds an interactive web tool and a GPU / full-atlas training workflow.
**The published auto-decoder model is unchanged.** At release the full 432,919-nucleus retrain was prepared
but not run. It was run by 2026-07-23 (commits 4d85728 and 7dc04df; "Phase 2 activated" below): the re-fit Flight-vs-GC classifier gives
AUC 0.915 ± 0.043 on 124 OSDR RNA-seq samples (`docs/assets/model/classifier.json`), against 0.919 ± 0.047 for
v1.0.0. The manuscript, the `Results/` tables and the figures still report v1.0.0.

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
- Limitation #1 was updated after the full-atlas GPU retrain ran (2026-07-23). The re-fit AUC, 0.915 ± 0.043,
  is below the ≥ 0.919 target in `GPU_TRAINING_PLAN.md` §5.

## [1.0.0] — 2026-07-20
- Initial release: end-to-end FAIR pipeline integrating 1,337 bulk transcriptomic samples with the
  Salk *Arabidopsis* Developmental Atlas foundation model; auto-decoder deconvolution, DESeq2
  differential expression, DerSimonian-Laird meta-analysis, elastic-net classification
  (Flight-vs-GC AUC = 0.919), and ggPlantMap / KEGG visualizations.
