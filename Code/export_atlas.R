#!/usr/bin/env Rscript
# =============================================================================
# export_atlas.R
#
# Reconstruction of the atlas-extraction step that feeds the stimulus auto-decoder
# (the original `fix_signatures.R` was not committed to the repo — it ran on the
# /mnt/shared-workspace compute box). This is rebuilt from two authoritative sources
# so it matches the real pipeline:
#   * Manuscript Methods §"Foundation Model" — Seurat v5 object, log-normalized data,
#     4,000 HVGs x 183 clusters keyed on `orig.cluster`, conditioning on organ (12) and
#     developmental stage (10); draft used a 60,792-cell stratified subsample over all 183 clusters.
#   * Code/train_autodecoder.py — the exact binary I/O contract it consumes (below).
#
# MODE:
#   N_CELLS = 0            -> export the FULL atlas (all 432,919 nuclei)  [GPU run — default]
#   N_CELLS = 60792 (etc.) -> stratified subsample covering all 183 clusters (reproduce the draft)
#
# Output contract (must match the trainer, byte-verified):
#   <OUT_DIR>/expr_sub.bin : float64, genes x cells, COLUMN-MAJOR (R default = NumPy order='F').
#                            Rows = HVGs in signature-matrix order; columns = cells.
#   <OUT_DIR>/meta_sub.csv : one row per cell (same order as columns of expr_sub.bin),
#                            columns: orig.cluster, orig.ident, dataset
#   <OUT_DIR>/dims.txt      : "<n_genes> <n_cells>"
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

# ---- configuration (EDIT paths / column names to your object) ---------------
ATLAS_RDS  <- "/mnt/shared-workspace/processed/atlas_compact_reference.rds"
SIG_MATRIX <- "/mnt/shared-workspace/processed/cell_type_signatures.csv"   # defines HVG set + order
OUT_DIR    <- "/mnt/shared-workspace/autodecoder"

ASSAY   <- "RNA"          # assay holding log-normalized data
LAYER   <- "data"         # Seurat v5 layer (log-normalized); on v4 this maps to slot="data"

# Metadata columns the trainer expects, mapped from your object's columns.
# Per Methods, cluster == "orig.cluster". Confirm organ/stage names on your object.
COL_CLUSTER <- "orig.cluster"   # 183 cell-type clusters
COL_ORGAN   <- "orig.ident"     # organ / sample  (expected 12 categories -> organ embedding 12->16)
COL_STAGE   <- "dataset"        # developmental stage (expected 10 -> stage embedding 10->16)

# 0 = full atlas (GPU). Set to 60792 to reproduce the draft's stratified subsample.
N_CELLS <- 0
MIN_PER_CLUSTER <- 20     # coverage floor when subsampling, so every cluster is represented
SEED    <- 42
# -----------------------------------------------------------------------------

set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- HVG order comes from the signature matrix so rows align with the trainer/deconvolution ---
sig  <- read.csv(SIG_MATRIX, row.names = 1, check.names = FALSE)
hvgs <- rownames(sig)
cat(sprintf("Signature HVGs: %d (expected 4000)\n", length(hvgs)))

cat("Loading Seurat v5 atlas ...\n")
atlas <- readRDS(ATLAS_RDS)
DefaultAssay(atlas) <- ASSAY

# --- validate metadata columns + cardinality (catches wrong object / wrong column names) ---
md_all <- atlas@meta.data
for (col in c(COL_CLUSTER, COL_ORGAN, COL_STAGE)) {
  if (!col %in% colnames(md_all))
    stop(sprintf("Metadata column '%s' not found. Columns present: %s",
                 col, paste(colnames(md_all), collapse = ", ")))
}
n_clu <- length(unique(md_all[[COL_CLUSTER]]))
n_org <- length(unique(md_all[[COL_ORGAN]]))
n_stg <- length(unique(md_all[[COL_STAGE]]))
cat(sprintf("Cardinality  clusters=%d (exp 183)  organs=%d (exp 12)  stages=%d (exp 10)\n",
            n_clu, n_org, n_stg))
if (n_clu != 183) cat("  WARNING: cluster count != 183 — check COL_CLUSTER.\n")
if (n_org != 12)  cat("  WARNING: organ count != 12 — check COL_ORGAN (embedding is 12->16).\n")
if (n_stg != 10)  cat("  WARNING: stage count != 10 — check COL_STAGE (embedding is 10->16).\n")

# --- choose cells: full atlas, or stratified subsample covering all clusters ---
all_cells <- colnames(atlas)
if (N_CELLS <= 0 || N_CELLS >= length(all_cells)) {
  sel_cells <- all_cells
  cat(sprintf("MODE: FULL atlas — exporting all %d nuclei.\n", length(sel_cells)))
} else {
  cat(sprintf("MODE: stratified subsample — target %d cells, >=%d per cluster.\n",
              N_CELLS, MIN_PER_CLUSTER))
  clu <- as.character(md_all[[COL_CLUSTER]]); names(clu) <- all_cells
  by_cluster <- split(all_cells, clu)
  # proportional allocation with a per-cluster floor, capped by availability
  sizes <- sapply(by_cluster, length)
  alloc <- pmax(MIN_PER_CLUSTER, round(N_CELLS * sizes / sum(sizes)))
  alloc <- pmin(alloc, sizes)
  sel_cells <- unlist(mapply(function(cells, k) sample(cells, k),
                             by_cluster, alloc[names(by_cluster)], SIMPLIFY = FALSE),
                      use.names = FALSE)
  sel_cells <- sample(sel_cells)  # shuffle so train/val split isn't cluster-ordered
  cat(sprintf("  selected %d cells across %d clusters.\n", length(sel_cells), length(by_cluster)))
}

# --- expression: genes x cells, HVG rows in signature order, log-normalized ---
expr <- tryCatch(
  GetAssayData(atlas, assay = ASSAY, layer = LAYER),          # Seurat v5
  error = function(e) GetAssayData(atlas, assay = ASSAY, slot = LAYER))  # v4 fallback

# Subset to the HVG ROWS FIRST (e.g. 27522 -> 4000) *before* any column op. A
# character-indexed column reorder on the full-gene matrix can transiently blow
# past available RAM (hit R's 24 GB vector cap on the 432k atlas); doing rows
# first keeps every subsequent op on the small 4000-row matrix.
missing <- setdiff(hvgs, rownames(expr))
if (length(missing) > 0) {
  cat(sprintf("WARNING: %d signature HVGs absent from atlas; filling with 0.\n", length(missing)))
  pad <- Matrix(0, nrow = length(missing), ncol = ncol(expr), sparse = TRUE)
  rownames(pad) <- missing
  expr <- rbind(expr, pad)
}
expr <- expr[hvgs, , drop = FALSE]                            # 4000 x n_cells sparse
# Select cells: a no-op reorder in full mode (skip it), a real subset otherwise.
if (!identical(as.character(sel_cells), colnames(expr)))
  expr <- expr[, sel_cells, drop = FALSE]
n_genes <- length(hvgs); n_cells <- ncol(expr)
cat(sprintf("Exporting %d genes x %d cells (%.1f GB float64); densifying per chunk to bound memory.\n",
            n_genes, n_cells, n_genes * n_cells * 8 / 1e9))

# --- metadata aligned to expr columns ---
md <- atlas@meta.data[sel_cells, , drop = FALSE]
write.csv(data.frame(
  orig.cluster = as.character(md[[COL_CLUSTER]]),
  orig.ident   = as.character(md[[COL_ORGAN]]),
  dataset      = as.character(md[[COL_STAGE]]),
  stringsAsFactors = FALSE),
  file.path(OUT_DIR, "meta_sub.csv"), row.names = FALSE)

# --- column-major binary, chunked over cells to bound memory on the full atlas ---
con <- file(file.path(OUT_DIR, "expr_sub.bin"), "wb")
chunk <- 5000L
for (start in seq(1L, n_cells, by = chunk)) {
  end <- min(start + chunk - 1L, n_cells)
  block <- as.matrix(expr[, start:end, drop = FALSE])         # densify ONLY this chunk (~160 MB)
  writeBin(as.double(block), con)                             # column-major within the chunk
  cat(sprintf("  wrote cells %d-%d\r", start, end))
}
close(con); cat("\n")
writeLines(paste(n_genes, n_cells), file.path(OUT_DIR, "dims.txt"))

cat("DONE. Next: python Code/train_autodecoder_gpu.py --in-memory  (see GPU_TRAINING_PLAN.md)\n")
cat("     (trainer reads dtype float64 by default; pass --input-dtype float64)\n")
