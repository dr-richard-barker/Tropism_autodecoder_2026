#!/usr/bin/env Rscript
# Subsystem 4: Run DESeq2 on each OSDR RNA-seq study
.libPaths(c("/workspace/.Rlib", .libPaths()))
suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
})

DE_DIR <- "/mnt/shared-workspace/processed/de_results"
dir.create(DE_DIR, showWarnings = FALSE, recursive = TRUE)

count_files <- list.files(DE_DIR, pattern="_counts_for_de.csv$", full.names=TRUE)
cat(sprintf("Found %d studies for DE analysis\n", length(count_files)))

run_de <- function(cf) {
  sid <- sub("_counts_for_de.csv$", "", basename(cf))
  cond_file <- file.path(DE_DIR, paste0(sid, "_conditions.csv"))
  if (!file.exists(cond_file)) return(data.frame(sid=sid, status="no_conditions"))

  cat(sprintf("\n=== %s ===\n", sid))
  counts <- fread(cf, data.table=FALSE)
  gene_col <- colnames(counts)[1]
  rownames(counts) <- counts[[gene_col]]
  counts[[gene_col]] <- NULL
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  counts <- counts[rowSums(counts) > 0, ]

  conds <- fread(cond_file, data.table=FALSE)
  conds <- conds[match(colnames(counts), conds$sample), ]
  keep <- !is.na(conds$condition)
  counts <- counts[, keep]
  conds <- conds[keep, ]
  conds$condition <- factor(conds$condition, levels=c("Ground Control","Space Flight"))

  cat(sprintf("  %d genes x %d samples (FLT=%d, GC=%d)\n",
              nrow(counts), ncol(counts), sum(conds$condition=="Space Flight"), sum(conds$condition=="Ground Control")))

  if (sum(conds$condition=="Space Flight") < 2 || sum(conds$condition=="Ground Control") < 2) {
    cat("  SKIP: need >=2 per group\n")
    return(data.frame(sid=sid, status="skipped_few_samples"))
  }

  dds <- DESeqDataSetFromMatrix(countData=counts, colData=conds, design=~condition)
  tryCatch({
    dds <- DESeq(dds, quiet=TRUE)
    res <- results(dds, contrast=c("condition","Space Flight","Ground Control"), alpha=0.05)
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df <- res_df[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    colnames(res_df) <- c("gene","baseMean","log2FC","lfcSE","stat","pvalue","padj")

    tryCatch({
      res_shrunk <- lfcShrink(dds, contrast=c("condition","Space Flight","Ground Control"), type="apeglm", quiet=TRUE)
      res_df$log2FC_shrunk <- res_shrunk$log2FoldChange
    }, error=function(e) res_df$log2FC_shrunk <<- res_df$log2FC)

    out_file <- file.path(DE_DIR, paste0(sid, "_de_results.tsv"))
    fwrite(res_df, out_file, sep="\t")

    n_sig <- sum(res_df$padj < 0.05, na.rm=TRUE)
    n_up <- sum(res_df$padj < 0.05 & res_df$log2FC > 0, na.rm=TRUE)
    n_down <- sum(res_df$padj < 0.05 & res_df$log2FC < 0, na.rm=TRUE)
    cat(sprintf("  DE: %d significant (padj<0.05): %d up, %d down\n", n_sig, n_up, n_down))

    norm_counts <- counts(dds, normalized=TRUE)
    norm_out <- file.path(DE_DIR, paste0(sid, "_normalized_counts.csv"))
    fwrite(as.data.frame(norm_counts, keep.rownames="gene"), norm_out)

    return(data.frame(sid=sid, status="done", n_genes=nrow(counts), n_sig=n_sig, n_up=n_up, n_down=n_down,
                      n_flight=sum(conds$condition=="Space Flight"), n_control=sum(conds$condition=="Ground Control")))
  }, error=function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    return(data.frame(sid=sid, status="error", message=conditionMessage(e)))
  })
}

results_list <- lapply(count_files, run_de)
results_df <- do.call(rbind, results_list)
summary_out <- file.path(DE_DIR, "de_analysis_summary.tsv")
fwrite(results_df, summary_out, sep="\t")
cat("\n=== SUMMARY ===\n")
print(results_df)
cat("\nDONE_DE\n")
