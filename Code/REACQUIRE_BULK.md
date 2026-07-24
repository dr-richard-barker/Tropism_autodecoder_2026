# Reproducing the bulk deconvolution + Flight-vs-GC classifier from public data

The original count-acquisition and intermediate workspace were not retained, so this documents how
the bulk pipeline (and the web tool's Phase 2 classifier) is reproduced entirely from public sources:
NASA **OSDR** for the bulk RNA-seq counts and the **GSE226097** Salk atlas for the reference.

All steps below were used to build the shipped `docs/assets/model/` bundle (`version 2.0.0-fullatlas`,
Flight-vs-GC AUC ≈ 0.915). Paths are examples — adjust to your layout.

## 0. Atlas reference (once)
Download `GSE226097_global_integration_221009.rds` from GEO, then:
```bash
# 4000 HVGs x 183-cluster signature matrix (see the manuscript Methods)
Rscript Code/build_signature_matrix.R      # -> cell_type_signatures.csv   (adapt the paths at top)
# full-atlas extraction + GPU/MPS retrain (produces cluster_stimulus_codes.json, embeddings)
Rscript Code/export_atlas.R                 # -> autodecoder/{expr_sub.bin,meta_sub.csv,dims.txt}
python  Code/train_autodecoder_gpu.py --accelerator mps --in-memory ...
```
(`build_signature_matrix.R` is the atlas-derived signature builder; on Apple Silicon use
`environment.mps.yml` and `--accelerator mps`.)

## 1. Re-acquire bulk counts from OSDR
```bash
python Code/fetch_osdr_counts.py --out bulk/counts
```
Finds each study's *unnormalized* count matrix (RSEM preferred, then STAR) via the live OSDR file
API and downloads `bulk/counts/OSD-<id>_raw.csv`. Default study set = the 15 RNA-seq studies in
`Data/harmonized_metadata.tsv` that have both Space Flight and Ground Control samples.

## 2. Build DESeq2 inputs (conditions from the committed metadata)
```bash
python Code/build_de_inputs.py --counts bulk/counts \
       --meta Data/harmonized_metadata.tsv --de-dir bulk/de_results
```
Resolves each count column to a condition by exact `sample_id`, then GSM, then a GC/FLT name token,
writing `OSD-<id>_counts_for_de.csv` + `OSD-<id>_conditions.csv`.

## 3. Differential expression
```bash
Rscript Code/run_de.R bulk/de_results      # or set the DE_DIR env var
```
Per-study DESeq2 (Flight vs Ground Control) → `*_de_results.tsv` + `*_normalized_counts.csv`.

## 4. Deconvolution onto the retrained atlas
```bash
python Code/project_bulk.py --atlas-dir autodecoder --proc-dir . \
       --de-dir bulk/de_results --out-dir bulk/projection
```
NNLS of each sample onto `cell_type_signatures.csv` (in `--proc-dir`) → 183 cell-type fractions,
then `fractions · cluster_stimulus_codes` → 32 stimulus scores. (`--geo-dir` optionally adds the GEO
count matrices.)

## 5. Classifier + meta-analysis
```bash
python Code/meta_classifier.py --meta Data/harmonized_metadata.tsv \
       --de-dir bulk/de_results --proj-dir bulk/projection --out-dir bulk/classifier
```
Elastic-net logistic on [183 fractions + 32 stim], nested 5-fold CV → `classifier_params.json`.

## 6. Export the web-tool Phase 2 bundle
```bash
python Code/export_web_artifacts.py \
  --sig-matrix cell_type_signatures.csv \
  --stim-codes autodecoder/cluster_stimulus_codes.json \
  --classifier classifier_params.json \
  --meta-analysis meta_analysis_results.tsv \
  --out docs/assets/model
```
Writes `manifest.phase2_ready: true`; commit `docs/assets/model/` and redeploy Pages — the tool's
Full-auto-decoder (Phase 2) method activates automatically.

## Scope notes
- Only **OSDR RNA-seq** studies are included in the classifier re-fit (matches the original 124
  Flight/GC samples matched by GSM). GEO series and microarray studies are not yet folded in.
- `run_de.R`, `project_bulk.py`, and `meta_classifier.py` take their paths from CLI flags / env vars
  (shown above); no in-file editing needed. `build_signature_matrix.R` and `export_atlas.R` still use
  top-of-file path constants — set those two before running.
