#!/usr/bin/env Rscript
# Build cell_type_signatures.csv (4000 HVGs x 183 clusters) from the GSE226097 global
# integration Seurat object — reproduces the pipeline's signature matrix from public data.
# Per the manuscript Methods: 4,000 HVGs (vst), per-cluster (orig.cluster) mean of the
# log-normalized (RNA/data) expression.  Edit the paths below for your setup.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

ATLAS   <- "GSE226097_global_integration_221009.rds"   # download from GEO GSE226097 supplementary
OUT     <- "cell_type_signatures.csv"
ASSAY   <- "RNA"; LAYER <- "data"; CLUSTER <- "orig.cluster"; NHVG <- 4000

cat("loading atlas ...\n"); t0 <- Sys.time()
obj <- readRDS(ATLAS)
DefaultAssay(obj) <- ASSAY
cat(sprintf("loaded %d genes x %d cells in %.0fs\n", nrow(obj), ncol(obj),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

obj  <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = NHVG, verbose = FALSE)
hvgs <- VariableFeatures(obj)
cat(sprintf("selected %d HVGs\n", length(hvgs)))

expr <- GetAssayData(obj, assay = ASSAY, layer = LAYER)[hvgs, , drop = FALSE]  # sparse, genes x cells
clu  <- as.character(obj@meta.data[[CLUSTER]])
clusters <- sort(unique(clu))

# per-cluster mean via a sparse cells x clusters indicator matrix-multiply (memory-frugal)
ind <- sparse.model.matrix(~ 0 + factor(clu, levels = clusters))
colnames(ind) <- clusters
sig <- as.matrix(expr %*% ind)                        # genes x clusters (summed)
sig <- sweep(sig, 2, Matrix::colSums(ind), "/")       # -> per-cluster means
rownames(sig) <- hvgs; colnames(sig) <- clusters

cat(sprintf("signature matrix: %d genes x %d clusters -> %s\n", nrow(sig), ncol(sig), OUT))
write.csv(sig, OUT)
