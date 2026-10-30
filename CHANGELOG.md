# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/); this project uses semantic versioning.

## [Unreleased]

Corrections from an evidence check of the repository (28 Oct 2026). `Results/PROVENANCE.md` says where each
result file now comes from.

### Fixed
- **Meta-analysis** (`Code/meta_classifier.py`): the DerSimonian–Laird τ² was too small by a factor of Σw, and
  p-values were computed as `1 − cdf`, which underflows to 0. Both are fixed; the function matches
  `metafor::rma(method = "DL")`. Re-run on DESeq2 results for the same 11 studies from OSDR counts re-acquired in
  July 2026 (the v1.0.0 per-study tables were not retained): 26,527 genes tested, **248 significant** (was 16,002).
  The volcano, forest and KEGG figures, the web tool's `flight_signature_genes.json` (now padj < 0.05 only, 215
  genes from the 14-study retrain inputs) and all manuscripts follow.
- **Classifier metrics recorded**: re-running `Code/meta_classifier.py` on the committed v1.0.0 projection tables
  reproduces AUC 0.919 ± 0.047, accuracy 0.85 and tropism F1 1.000; the script now writes
  `classifier_metrics.json`, and `tropism_confusion_matrix.csv` (80 phototropism, 44 gravitropism) is committed.
  The cross-validation is stratified 5-fold with fixed hyperparameters, not nested; the text and the web
  tool's figure captions now say so, and the web page labels its 0.919 headline as the v1.0.0 model's.
- **Full-atlas training log** committed (`Results/autodecoder_fullatlas_training_log.csv`).
- **ggPlantMap maps** redrawn without the random per-region jitter: every region of a map has its organ's mean
  difference, on one colour scale (±6.2 × 10⁻⁴) shared by all five maps, and each subtitle states that value.
  (A colour scale fitted to a map's single value would paint every region the midpoint colour.) The text
  reports organ-level differences only; the map titles and the Discussion no longer claim tropism signalling
  or say where in the plant it is perceived.
- **Organ-level results**: the 0.003 / 0.027 / 0.030 reported as organ "mean padj" are each organ's lowest
  cluster padj (`Code/visualization.R` computed `min(padj)` under the name `mean_padj`, now `min_padj`); the
  text now gives organ means. The bar-plot captions list the six organs the plot shows, ordered by magnitude,
  with no error bars (the plot has none). The abstracts no longer say the deconvolution recovers
  "tissue-level tropism signatures".
- **Discussion** (LaTeX and md): no cluster in the repository is annotated as vascular, so the classifier's top
  features and the significant shifts no longer "suggest vascular tissue differentiation". The three
  significant shifts are stem_9, seed_425d_7 and silique_20 (seed, not "late-seedling" as written).
- **Deconvolution results** (all three manuscripts, README): "all 398 samples successfully deconvolved" counted
  a parsing artefact (GSE143760 `gene.1`, every fraction zero); 397 are samples. The atlas has no root clusters,
  and root samples map mostly to silique (36 %), seedling (32 %) and seed (21 %) clusters, so "root samples mapped
  predominantly to seedling and root clusters, as expected" is withdrawn and a limitation added. The standalone
  `Manuscript/main.tex` also named "columella clusters" (there are none), silique samples (none were
  deconvolved) and fraction ranges that do not recompute; that paragraph is rewritten from
  `Results/cell_type_fractions_all.csv`.
- **GSE8300 removed** from `Data/harmonized_metadata.tsv`: GEO lists it as a *Haemophilus influenzae* genomic
  hybridization series, leaving 1,253 metadata rows (63 from 5 GEO series).
- **Duplicate metadata rows**: 48 rows of `Data/harmonized_metadata.tsv` were exact copies (OSD-120 36, OSD-121
  12) and are removed, leaving 1,205 rows. 42 samples are listed under two OSDR studies (OSD-218/OSD-219,
  OSD-147/OSD-205), so the dataset is **1,163 distinct samples** (was reported as 1,253): OSDR 1,100, GEO 63;
  gravitropism 901, phototropism 198, both 46 (10 samples are annotated gravitropism in OSD-147 and
  phototropism in OSD-205), mechanotropism 12, hydrotropism 6; 787 flight, 376 ground; 1,004 RNA-seq, 159
  microarray. DE, the classifier and deconvolution are unaffected: they read count-matrix columns or match
  by GSM. `Figures/dataset_summary` counts distinct samples.
- **Citations** checked against CrossRef: Toyota & Gilroy 2013; CIBERSORTx = Newman et al. 2019 *Nat Biotechnol*;
  PhysioSpace = Lenz et al. 2013 (the "Schiller 2021" entry matched no paper); ggPlantmap title and pages.
- **`Manuscript/main.tex`**: GSE111555 (a human exercise study) replaced by GSE115554; counts, classifier
  features, tropism class sizes and the classifier method (scikit-learn, not glmnet; no inner CV) now match the
  code and data.
- **`metadata.json`**: the identifier `10.5281/zenodo.arabidopsis-tropism` resolved nowhere; it is now the
  placeholder `10.5281/zenodo.XXXXXXX` until a Zenodo DOI is minted.
- `seed_425d_7` flight difference is +0.04 %; figure subtitles are computed, not typed.

### Still TODO
- v1.0.0 auto-decoder metrics (val_loss 0.4300 etc.) have no surviving log (README, Limitation 8).
- Deposit the re-acquired DESeq2 tables on Zenodo with the other retrain artefacts.

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
