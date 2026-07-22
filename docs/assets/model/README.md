# Phase 2 model artifacts

This directory holds the exported model artifacts that let the web tool reproduce the manuscript
pipeline's inference in the browser (NNLS cell-type deconvolution → stimulus projection →
Flight-vs-Ground-Control elastic-net classifier).

**It ships empty** except for a stub `manifest.json` with `phase2_ready: false`, so the live tool
stays on Phase 1 (curated markers) and shows "Phase 2 unavailable" until real artifacts are added.

## To activate Phase 2

Run the pipeline (or the full-atlas GPU retrain), then:

```bash
python Code/export_web_artifacts.py \
  --sig-matrix  /path/to/cell_type_signatures.csv \
  --stim-codes  /path/to/cluster_stimulus_codes.json \
  --classifier  /path/to/classifier_params.json \
  --meta-analysis /path/to/meta_analysis_results.tsv \
  --out docs/assets/model
```

This writes `signature_matrix.bin`, `signature_index.json`, `stimulus_codes.json`,
`classifier.json`, `flight_signature_genes.json`, and a `manifest.json` with
`phase2_ready: true`. Commit the directory and redeploy Pages — the tool enables the
"Full auto-decoder (Phase 2)" method automatically.

`classifier_params.json` is produced by `Code/meta_classifier.py` (it now dumps the fitted
coefficients, intercept, and StandardScaler parameters — the feature-importance CSV alone is
insufficient because it omits the intercept and scaler).
