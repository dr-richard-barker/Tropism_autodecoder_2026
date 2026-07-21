# Full-Atlas GPU Retraining Plan — Stimulus Auto-Decoder

**Repo:** `Tropism_autodecoder_2026` (Barker, 2026) · **Subsystem 3** (auto-decoder training)
**Goal:** retrain the conditional-VAE stimulus auto-decoder on the **complete Salk atlas (432,919 nuclei)** using **GPU acceleration**, replacing the first draft's **~60k-cell CPU subsample**.
**Why it matters:** this is **Limitation #1** in the repository README ("No GPU: auto-decoder trained on CPU with 60k-cell subsample"). Closing it is a direct, reviewer-visible improvement and yields better per-cluster stimulus codes for all 183 clusters.

---

## 1. Where the 60k limit actually lives (important)

I read `Code/train_autodecoder.py`. It does **not** subsample — it consumes a pre-extracted binary:

```
<atlas_dir>/expr_sub.bin   float64, genes × cells, column-major (Fortran)
<atlas_dir>/meta_sub.csv   orig.cluster, orig.ident, dataset  (one row per cell)
<atlas_dir>/dims.txt        "<n_genes> <n_cells>"
```

The `.bin` is built upstream by an extraction step (referenced in the script as `fix_signatures.R`) that is **not committed to the repo's `Code/`**. **That is where the 60k subsample is applied.** So a true full-atlas run requires changes in **two** places:

| # | Change | File |
|---|--------|------|
| A | Export **all 432,919 nuclei** (no subsample) in the existing binary format | `Code/export_atlas.R` *(reconstructed from the Methods + trainer contract; `N_CELLS=0` = full, or set to `60792` to reproduce the draft's stratified subsample)* |
| B | Train on **GPU**, streaming the large matrix so it fits in memory | `Code/train_autodecoder_gpu.py` *(provided)* |

The trainer itself was also CPU-locked: `accelerator="cpu"` is hardcoded. That is fixed in the GPU version.

---

## 2. What the GPU trainer changes (and what it deliberately keeps)

`Code/train_autodecoder_gpu.py` — **model architecture and loss are byte-for-byte the same** as the draft (encoder 4000→512→256, latent 32, decoder 96→512→512→4000, cluster/organ/stage embeddings, MSE + 0.5·KLD + 0.1·CE), so results stay comparable. Changes are confined to device, data path, and configuration:

1. **Device / precision:** `accelerator=auto` (GPU when present), `precision=bf16-mixed` on Ampere+ with automatic fallback to `16-mixed`, then `32` on CPU. `torch.set_float32_matmul_precision("high")`.
2. **Scale to 432k:** the expression binary is **memory-mapped** and cells are read on demand in the `DataLoader` (no full-matrix RAM load). `--in-memory` optionally loads the whole atlas as float32 for speed when RAM allows.
3. **Reproducible log-transform:** `--log-transform auto|on|off` (auto samples ≤20k cells to decide, instead of the draft's global-max check on an in-RAM matrix).
4. **Scheduler fix:** `CosineAnnealingLR(T_max=epochs)` now tracks `--epochs` (was hardcoded 40).
5. **Batched GPU embedding pass** writes the same outputs the draft did — `atlas_latent_embeddings.npy`, `atlas_cls_predictions.npy`, and `cluster_stimulus_codes.json` (the 183×32 per-cluster stimulus codes).
6. **Fully argparse-driven**; defaults reproduce the manuscript methodology on GPU + full atlas.
7. **Cardinality guard at startup:** errors out if `meta_sub.csv` is missing a required column, or if the distinct cluster/organ/stage counts differ from the expected **183 / 12 / 10** (the sizes the embeddings assume). Blocks by default so a wrong-object or partial export can't reach training unnoticed; tune with `--expect-clusters/-organs/-stages` (0 disables a check) or bypass with `--allow-cardinality-mismatch`.

**Format contract verified:** the R chunked column-major writer and the trainer's `np.fromfile(...).reshape(n_hvgs, n_cells, order="F")` were confirmed (round-trip test) to reconstruct the identical matrix, including per-cell columns.

---

## 3. Compute & memory estimate (full atlas)

Model is small (~4.6M parameters); the workload is **memory/IO-bound, not compute-bound**.

| Resource | Full-atlas figure | Note |
|---|---|---|
| Cells | 432,919 | vs ~60,000 draft (≈7×) |
| Steps/epoch @ batch 512 | ~846 | |
| Expression file (float64) | ~13.9 GB on disk | `--input-dtype float64` |
| …as float32 | ~6.9 GB | pre-convert to halve IO (recommended) |
| Host RAM for `--in-memory` (float32) | ~6.9 GB | recommend ≥16 GB host |
| GPU VRAM | < 2 GB needed | any ≥8 GB card is ample |
| **Wall-clock, 40 epochs** | **≈10–40 min on one A100 / RTX 4090-class GPU** | IO-dominated; `--in-memory` + float32 lands at the low end |

Even a mid-range consumer GPU finishes well inside an hour. Early stopping (patience 10) may end sooner (the draft's best was epoch 5).

---

## 4. Run recipe

```bash
# 0) one-time env
mamba env create -f environment.gpu.yml && conda activate arabidopsis-tropism-gpu

# 1) export the FULL atlas (edit paths at the top of the script first; N_CELLS=0 = full)
Rscript Code/export_atlas.R               # -> expr_sub.bin (all 432,919 cells), meta_sub.csv, dims.txt

# 2) GPU training on the full atlas
python Code/train_autodecoder_gpu.py \
    --atlas-dir /mnt/shared-workspace/autodecoder \
    --sig-matrix /mnt/shared-workspace/processed/cell_type_signatures.csv \
    --out-dir   /mnt/shared-workspace/autodecoder \
    --in-memory --input-dtype float64 --epochs 40 --batch-size 512 --precision bf16-mixed

# or containerized:
docker build -f Dockerfile.gpu -t arabidopsis-tropism-gpu .
docker run --gpus all -v /path/to/workspace:/mnt/shared-workspace arabidopsis-tropism-gpu \
    --atlas-dir /mnt/shared-workspace/autodecoder --in-memory
```

Downstream (unchanged) after retraining: re-run `project_bulk.py`, `meta_classifier.py`, `visualization.R`, `viz_ggplantmap.R` so deconvolution, classifier, and figures use the new full-atlas codes.

---

## 5. Acceptance criteria

- [ ] `dims.txt` reports **~432,919 cells** and all **183 clusters** appear in `meta_sub.csv` (full atlas ⇒ every cluster covered by real cells, not the subsample approximation).
- [ ] Training completes on GPU; `val_recon` ≤ the draft's **0.140** (more data should match or improve reconstruction); auxiliary cluster-classification loss improves.
- [ ] `cluster_stimulus_codes.json` contains **183** entries × **32** dims; `atlas_latent_embeddings.npy` is `(432919, 32)`.
- [ ] Re-run downstream: **Flight-vs-GC AUC ≥ 0.919** (should hold or improve); regenerate figures.
- [ ] Update README Limitation #1 to reflect full-atlas GPU training; bump `metadata.json` version.

## 6. Reproducibility & caveats

- `pl.seed_everything(seed, workers=True)` is set. GPU kernels are not bit-deterministic by default; for strict determinism add `torch.use_deterministic_algorithms(True)` and `CUBLAS_WORKSPACE_CONFIG=:4096:8` (small speed cost). Record the GPU model, driver, CUDA, and library versions alongside results.
- `export_atlas.R` is a **faithful reconstruction** of the original (uncommitted) `fix_signatures.R`, rebuilt from the manuscript Methods (Seurat v5, log-normalized `RNA/data`, 4,000 HVGs, `orig.cluster`, organ=12, stage=10) and the trainer's exact binary contract. It is **not** the original bytes — it runs cardinality checks (183 clusters / 12 organs / 10 stages) at startup and errors out if a column name doesn't match, so a wrong-object/column mistake fails loudly rather than silently. Confirm the three column mappings once against your atlas object.
- float16 input (`--input-dtype float16`) quarters IO; acceptable here because data are log-normalized and clipped to [0,10], but validate `val_recon` if used.

## 7. Connection to the web tool

The retrained `cluster_stimulus_codes.json` (183×32) and the signature matrix are exactly the artifacts the **web tool's Phase 2** needs to reproduce the auto-decoder in-browser (see `WEB_TOOL_PLAN.md` §5). Doing the full-atlas retrain first means Phase 2 ships the definitive codes rather than draft ones.
