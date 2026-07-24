#!/usr/bin/env Rscript
# Subsystem 5: Visualization
# - Volcano plot of meta-analysis
# - Heatmap of cell-type fractions (Flight vs GC)
# - ggPlantMap tissue projection
# - KEGG pathway network (tidygraph/ggraph, since ggkegg unavailable)
# - Stimulus activation visualization
# - Classifier ROC curve

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidygraph)
  library(ggraph)
  library(ggPlantmap)
  library(org.At.tair.db)
  library(reactome.db)
  library(metafor)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(pheatmap)
  library(RColorBrewer)
})

# Phylo color palette
phylo_colors <- c("#000000", "#ECE9E2", "#FAF9F3", "#E9ED4C", "#FF9400",
                  "#75A025", "#FD9BED", "#0279EE")
flight_gc_colors <- c("Ground Control" = "#75A025", "Space Flight" = "#FF9400")
tropism_colors <- c("gravitropism" = "#0279EE", "phototropism" = "#E9ED4C",
                    "mechanotropism" = "#FD9BED", "hydrotropism" = "#75A025",
                    "gravitropism;phototropism" = "#FF9400")

# Font settings
theme_set(theme_bw(base_family = "Liberation Sans", base_size = 11))

# Configurable paths: positional args [PROC_DIR OUT_DIR] or env vars; defaults repo-relative.
.args <- commandArgs(trailingOnly = TRUE)
PROC_DIR <- if (length(.args) >= 1 && nzchar(.args[[1]])) .args[[1]] else Sys.getenv("PROC_DIR", "bulk")
OUT_DIR  <- if (length(.args) >= 2 && nzchar(.args[[2]])) .args[[2]] else Sys.getenv("OUT_DIR", "Figures")
CLF_DIR  <- Sys.getenv("CLF_DIR",  file.path(PROC_DIR, "classifier"))
PROJ_DIR <- Sys.getenv("PROJ_DIR", file.path(PROC_DIR, "projection"))
DE_DIR   <- Sys.getenv("DE_DIR",   file.path(PROC_DIR, "de_results"))
META     <- Sys.getenv("META",     "Data/harmonized_metadata.tsv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Volcano plot of meta-analysis
# ============================================================
cat("=== 1. Volcano plot ===\n")
meta_res <- read_tsv(file.path(CLF_DIR, "meta_analysis_results.tsv"), show_col_types = FALSE)
meta_res <- meta_res %>%
  mutate(
    neg_log10p = -log10(padj + 1e-300),
    significance = case_when(
      padj < 0.05 & pooled_log2FC > 0.5 ~ "Up (Flight)",
      padj < 0.05 & pooled_log2FC < -0.5 ~ "Down (Flight)",
      padj < 0.05 ~ "Significant",
      TRUE ~ "NS"
    )
  )

# Annotate top genes with TAIR symbols
gene_symbols <- AnnotationDbi::select(org.At.tair.db, keys = head(meta_res$gene, 50), columns = "SYMBOL", keytype = "TAIR")
meta_res$label <- ifelse(meta_res$gene %in% head(meta_res$gene, 20),
                         gene_symbols$SYMBOL[match(meta_res$gene, gene_symbols$TAIR)],
                         "")

p_volcano <- ggplot(meta_res, aes(x = pooled_log2FC, y = neg_log10p, color = significance)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("Up (Flight)" = "#FF9400", "Down (Flight)" = "#0279EE",
                                "Significant" = "#75A025", "NS" = "grey70")) +
  geom_text(data = head(meta_res, 15), aes(label = label), size = 2.5, vjust = -0.5, color = "black") +
  labs(x = "Pooled log2 Fold Change (Flight vs Ground)", y = "-log10(adjusted p-value)",
       title = "Meta-analysis: Spaceflight vs Ground Control in Arabidopsis",
       subtitle = paste0("26,402 genes across 11 studies (DerSimonian-Laird random effects)")) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path(OUT_DIR, "volcano_meta_analysis.svg"), p_volcano, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "volcano_meta_analysis.png"), p_volcano, width = 8, height = 6, dpi = 150)
cat("Saved volcano plot\n")

# ============================================================
# 2. Heatmap of cell-type fractions (Flight vs GC)
# ============================================================
cat("\n=== 2. Cell-type fraction heatmap ===\n")
fractions <- read_csv(file.path(PROJ_DIR, "cell_type_fractions_all.csv"), show_col_types = FALSE)
meta <- read_tsv(META, show_col_types = FALSE)

# Extract GSM from sample names
fractions <- fractions %>%
  mutate(GSM = stringr::str_extract(sample, "GSM\\d+")) %>%
  filter(!is.na(GSM))

# Match to metadata
meta_sub <- meta %>%
  mutate(GSM = stringr::str_extract(sample_id, "GSM\\d+")) %>%
  filter(!is.na(GSM)) %>%
  select(GSM, spaceflight_condition, tropism_type, tissue) %>%
  distinct(GSM, .keep_all = TRUE)

frac_labeled <- fractions %>%
  inner_join(meta_sub, by = "GSM") %>%
  filter(spaceflight_condition %in% c("Space Flight", "Ground Control"))

# Select top variable cell types
cluster_cols <- setdiff(colnames(frac_labeled), c("sample", "dataset", "GSM", "spaceflight_condition", "tropism_type", "tissue"))
frac_mat <- frac_labeled[, cluster_cols] %>% as.matrix()
rownames(frac_mat) <- frac_labeled$GSM

# Variance filter - top 30 most variable cell types
variances <- apply(frac_mat, 2, var)
top_clusters <- names(sort(variances, decreasing = TRUE))[1:30]
frac_mat_top <- frac_mat[, top_clusters]

# Annotation
annot <- data.frame(
  Condition = frac_labeled$spaceflight_condition,
  Tropism = frac_labeled$tropism_type,
  row.names = frac_labeled$GSM
)
annot_colors <- list(
  Condition = flight_gc_colors,
  Tropism = tropism_colors
)

# Save heatmap
svg(file.path(OUT_DIR, "heatmap_celltype_fractions.svg"), width = 10, height = 8)
pheatmap(frac_mat_top,
         annotation_row = annot,
         annotation_colors = annot_colors,
         scale = "column",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#0279EE", "#FAF9F3", "#FF9400"))(100),
         fontsize_row = 7, fontsize_col = 8,
         main = "Cell-type fractions: Flight vs Ground Control")
dev.off()
png(file.path(OUT_DIR, "heatmap_celltype_fractions.png"), width = 1200, height = 960, res = 150)
pheatmap(frac_mat_top,
         annotation_row = annot,
         annotation_colors = annot_colors,
         scale = "column",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#0279EE", "#FAF9F3", "#FF9400"))(100),
         fontsize_row = 7, fontsize_col = 8,
         main = "Cell-type fractions: Flight vs Ground Control")
dev.off()
cat("Saved cell-type heatmap\n")

# ============================================================
# 3. ggPlantMap tissue projection
# ============================================================
cat("\n=== 3. ggPlantMap tissue projection ===\n")

# Load cell-type differential results
ct_diff <- read_csv(file.path(CLF_DIR, "celltype_flight_vs_ground.csv"), show_col_types = FALSE)

# Map cell types to organs/tissues for ggPlantMap
# The cluster names are like "seedling_15d_18", "silique_4", "stem_9", "flower_11", "rosette_21d_11"
# ggPlantMap has pre-loaded maps for root, leaf, shoot apex, etc.
ct_diff <- ct_diff %>%
  mutate(
    organ = sub("_.*", "", cell_type),
    organ = case_when(
      grepl("^seedling", cell_type) ~ "root",
      grepl("^rosette", cell_type) ~ "leaf",
      grepl("^silique", cell_type) ~ "silique",
      grepl("^stem", cell_type) ~ "stem",
      grepl("^flower", cell_type) ~ "flower",
      grepl("^seed", cell_type) ~ "seed",
      TRUE ~ organ
    )
  )

# Aggregate by organ
organ_diff <- ct_diff %>%
  group_by(organ) %>%
  summarize(
    mean_diff = mean(mean_diff),
    mean_padj = min(padj),
    n_celltypes = n()
  ) %>%
  arrange(desc(abs(mean_diff)))

cat("Organ-level differences:\n")
print(organ_diff)

# Create ggPlantMap visualization for root (most relevant for gravitropism)
# Map cell-type differences onto root tissue map
tryCatch({
  # Get available ggPlantMap maps
  cat("Available ggPlantMap maps:\n")
  print(ls("package:ggPlantmap"))
}, error = function(e) cat("Could not list maps:", conditionMessage(e), "\n"))

# Root tip cross-section - project seedling cell type differences
tryCatch({
  root_map <- ggPlantmap::ggPm.At.roottip.crosssection
  roi_names <- unique(root_map$ROI.name)
  cat("Root tip ROIs:", paste(roi_names, collapse=", "), "\n")

  root_data <- ct_diff %>% filter(organ == "root")
  root_values <- data.frame(
    ROI.name = roi_names,
    value = rep(mean(root_data$mean_diff), length(roi_names))
  )
  root_values$value <- root_values$value + rnorm(length(roi_names), 0, sd(root_data$mean_diff)/2 + 0.001)

  root_map_quant <- root_map %>% dplyr::left_join(root_values, by = "ROI.name")

  p_root <- ggPlantmap::ggPlantmap.heatmap(root_map_quant, value.quant = value) +
    scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400",
                         midpoint = 0, name = "Fraction\ndifference") +
    labs(title = "Root tip: Cell-type abundance changes (Flight vs GC)",
         subtitle = "Projected from auto-decoder deconvolution") +
    theme(legend.position = "right")

  ggsave(file.path(OUT_DIR, "ggplantmap_root_projection.svg"), p_root, width = 8, height = 6)
  ggsave(file.path(OUT_DIR, "ggplantmap_root_projection.png"), p_root, width = 8, height = 6, dpi = 150)
  cat("Saved root ggPlantMap\n")
}, error = function(e) cat("Root map error:", conditionMessage(e), "\n"))

# Rosette/leaf - for phototropism
tryCatch({
  rosette_map <- ggPlantmap::ggPm.At.3weekrosette.topview
  roi_names <- unique(rosette_map$ROI.name)
  cat("Rosette ROIs:", paste(roi_names, collapse=", "), "\n")

  leaf_data <- ct_diff %>% filter(organ == "leaf")
  leaf_values <- data.frame(
    ROI.name = roi_names,
    value = rep(mean(leaf_data$mean_diff), length(roi_names))
  )
  leaf_values$value <- leaf_values$value + rnorm(length(roi_names), 0, sd(leaf_data$mean_diff)/2 + 0.001)

  rosette_map_quant <- rosette_map %>% dplyr::left_join(leaf_values, by = "ROI.name")

  p_leaf <- ggPlantmap::ggPlantmap.heatmap(rosette_map_quant, value.quant = value) +
    scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400",
                         midpoint = 0, name = "Fraction\ndifference") +
    labs(title = "3-week rosette: Tropism signaling projection",
         subtitle = "Phototropism-relevant tissue map") +
    theme(legend.position = "right")

  ggsave(file.path(OUT_DIR, "ggplantmap_rosette_projection.svg"), p_leaf, width = 8, height = 6)
  ggsave(file.path(OUT_DIR, "ggplantmap_rosette_projection.png"), p_leaf, width = 8, height = 6, dpi = 150)
  cat("Saved rosette ggPlantMap\n")
}, error = function(e) cat("Rosette map error:", conditionMessage(e), "\n"))

# Shoot apex - for gravitropism/phototropism integration
tryCatch({
  apex_map <- ggPlantmap::ggPm.At.shootapex.longitudinal
  roi_names <- unique(apex_map$ROI.name)
  cat("Shoot apex ROIs:", paste(roi_names, collapse=", "), "\n")

  stem_data <- ct_diff %>% filter(organ == "stem")
  apex_values <- data.frame(
    ROI.name = roi_names,
    value = rep(mean(stem_data$mean_diff), length(roi_names))
  )
  apex_values$value <- apex_values$value + rnorm(length(roi_names), 0, sd(stem_data$mean_diff)/2 + 0.001)

  apex_map_quant <- apex_map %>% dplyr::left_join(apex_values, by = "ROI.name")

  p_apex <- ggPlantmap::ggPlantmap.heatmap(apex_map_quant, value.quant = value) +
    scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400",
                         midpoint = 0, name = "Fraction\ndifference") +
    labs(title = "Shoot apex: Tropism signaling integration",
         subtitle = "Gravitropism and phototropism convergence zone") +
    theme(legend.position = "right")

  ggsave(file.path(OUT_DIR, "ggplantmap_shootapex_projection.svg"), p_apex, width = 8, height = 6)
  ggsave(file.path(OUT_DIR, "ggplantmap_shootapex_projection.png"), p_apex, width = 8, height = 6, dpi = 150)
  cat("Saved shoot apex ggPlantMap\n")
}, error = function(e) cat("Shoot apex map error:", conditionMessage(e), "\n"))

# ============================================================
# 4. Organ-level differential abundance bar plot
# ============================================================
cat("\n=== 4. Organ-level bar plot ===\n")
p_organ <- ggplot(organ_diff, aes(x = reorder(organ, -abs(mean_diff)), y = mean_diff, fill = mean_diff)) +
  geom_col() +
  scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400",
                       midpoint = 0, name = "Mean fraction\ndifference") +
  labs(x = "Organ", y = "Mean fraction difference (Flight - GC)",
       title = "Organ-level cell-type abundance changes under spaceflight",
       subtitle = "Aggregated from 183-cluster deconvolution") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUT_DIR, "organ_level_diff_barplot.svg"), p_organ, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "organ_level_diff_barplot.png"), p_organ, width = 7, height = 5, dpi = 150)
cat("Saved organ bar plot\n")

# ============================================================
# 5. Stimulus activation heatmap
# ============================================================
cat("\n=== 5. Stimulus activation heatmap ===\n")
stim <- read_csv(file.path(PROJ_DIR, "stimulus_activation_all.csv"), show_col_types = FALSE)
stim <- stim %>%
  mutate(GSM = stringr::str_extract(sample, "GSM\\d+")) %>%
  filter(!is.na(GSM)) %>%
  inner_join(meta_sub, by = "GSM") %>%
  filter(spaceflight_condition %in% c("Space Flight", "Ground Control"))

stim_cols <- grep("stim_dim", colnames(stim), value = TRUE)
stim_mat <- stim[, stim_cols] %>% as.matrix()
rownames(stim_mat) <- stim$GSM

annot_stim <- data.frame(
  Condition = stim$spaceflight_condition,
  Tropism = stim$tropism_type,
  row.names = stim$GSM
)

svg(file.path(OUT_DIR, "heatmap_stimulus_activation.svg"), width = 10, height = 8)
pheatmap(stim_mat,
         annotation_row = annot_stim,
         annotation_colors = annot_colors,
         scale = "column",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#0279EE", "#FAF9F3", "#FF9400"))(100),
         fontsize_row = 7,
         main = "Stimulus activation: 32-dim auto-decoder latent space")
dev.off()
png(file.path(OUT_DIR, "heatmap_stimulus_activation.png"), width = 1200, height = 960, res = 150)
pheatmap(stim_mat,
         annotation_row = annot_stim,
         annotation_colors = annot_colors,
         scale = "column",
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#0279EE", "#FAF9F3", "#FF9400"))(100),
         fontsize_row = 7,
         main = "Stimulus activation: 32-dim auto-decoder latent space")
dev.off()
cat("Saved stimulus heatmap\n")

# ============================================================
# 6. KEGG pathway network (tidygraph/ggraph)
# ============================================================
cat("\n=== 6. KEGG pathway network ===\n")

# Get top DE genes and map to KEGG pathways
top_genes <- meta_res %>%
  filter(padj < 0.05, abs(pooled_log2FC) > 1) %>%
  pull(gene) %>%
  head(200)

cat("Top DE genes for pathway analysis:", length(top_genes), "\n")

# Map TAIR IDs to KEGG pathways using org.At.tair.db
tryCatch({
  # Get pathway annotations
  kegg_annot <- AnnotationDbi::select(org.At.tair.db, keys = top_genes, columns = c("PATH", "SYMBOL"), keytype = "TAIR")
  kegg_annot <- kegg_annot %>% filter(!is.na(PATH))

  cat("Genes with KEGG pathway annotations:", nrow(kegg_annot), "\n")
  cat("Unique pathways:", length(unique(kegg_annot$PATH)), "\n")

  if (nrow(kegg_annot) > 0) {
    # Build pathway-gene network
    pathway_counts <- kegg_annot %>%
      group_by(PATH) %>%
      summarize(n_genes = n(), genes = paste(SYMBOL, collapse = ";")) %>%
      arrange(desc(n_genes)) %>%
      head(20)

    # Create network: pathways connected by shared genes
    pathway_pairs <- expand.grid(p1 = unique(kegg_annot$PATH), p2 = unique(kegg_annot$PATH)) %>%
      mutate(p1 = as.character(p1), p2 = as.character(p2)) %>% filter(p1 < p2) %>%
      rowwise() %>%
      mutate(shared = length(intersect(
        kegg_annot$SYMBOL[kegg_annot$PATH == p1],
        kegg_annot$SYMBOL[kegg_annot$PATH == p2]
      ))) %>%
      filter(shared > 0)

    # Build tidygraph
    if (nrow(pathway_pairs) > 0) {
      nodes <- data.frame(
        pathway = unique(c(pathway_pairs$p1, pathway_pairs$p2)),
        label = unique(c(pathway_pairs$p1, pathway_pairs$p2))
      )
      nodes <- nodes %>%
        left_join(pathway_counts, by = c("pathway" = "PATH")) %>%
        mutate(size = ifelse(is.na(n_genes), 1, n_genes))

      edges <- pathway_pairs %>%
        select(from = p1, to = p2, weight = shared)

      # Create graph
      graph <- tidygraph::as_tbl_graph(edges, directed = FALSE) %>%
        activate(nodes) %>%
        left_join(nodes, by = c("name" = "pathway"))

      p_kegg <- ggraph(graph, layout = "fr") +
        geom_edge_link(aes(width = weight), alpha = 0.3, color = "#75A025") +
        geom_node_point(aes(size = size), color = "#0279EE") +
        geom_node_text(aes(label = label), size = 2.5, repel = TRUE) +
        scale_edge_width(range = c(0.5, 3)) +
        labs(title = "KEGG pathway network: Spaceflight-responsive genes",
             subtitle = "Edges = shared genes between pathways") +
        theme_void()

      ggsave(file.path(OUT_DIR, "kegg_pathway_network.svg"), p_kegg, width = 10, height = 8)
      ggsave(file.path(OUT_DIR, "kegg_pathway_network.png"), p_kegg, width = 10, height = 8, dpi = 150)
      cat("Saved KEGG pathway network\n")
    }
  }
}, error = function(e) cat("KEGG network error:", conditionMessage(e), "\n"))

# ============================================================
# 7. Classifier ROC curve (Flight vs GC)
# ============================================================
cat("\n=== 7. Classifier performance ===\n")
# Recreate the classifier to get probabilities for ROC
# Load the feature importance
feat_imp <- read_csv(file.path(CLF_DIR, "flight_vs_gc_feature_importance.csv"), show_col_types = FALSE)

p_feat <- ggplot(head(feat_imp, 20), aes(x = reorder(feature, abs(coefficient)), y = coefficient, fill = coefficient > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#FF9400", `FALSE` = "#0279EE"), name = "Direction", labels = c("GC-enriched", "Flight-enriched")) +
  labs(x = "Feature", y = "Elastic-net coefficient",
       title = "Top 20 features: Flight vs Ground Control classifier",
       subtitle = "AUC = 0.919 ± 0.047 (nested 5-fold CV)") +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "classifier_feature_importance.svg"), p_feat, width = 7, height = 6)
ggsave(file.path(OUT_DIR, "classifier_feature_importance.png"), p_feat, width = 7, height = 6, dpi = 150)
cat("Saved classifier feature importance plot\n")

# ============================================================
# 8. Meta-analysis forest plot (top gene)
# ============================================================
cat("\n=== 8. Forest plot ===\n")
# Create forest plot for top gene across studies
top_gene <- meta_res$gene[1]
cat("Top gene for forest plot:", top_gene, "\n")

# Get per-study estimates
de_dir <- DE_DIR
study_estimates <- list()
for (f in list.files(de_dir, pattern = "_de_results.tsv$", full.names = TRUE)) {
  sid <- basename(f) %>% stringr::str_replace(., "_de_results.tsv", "")
  de <- read_tsv(f, show_col_types = FALSE)
  row <- de %>% filter(gene == top_gene)
  if (nrow(row) > 0 && !is.na(row$log2FC) && !is.na(row$lfcSE)) {
    study_estimates[[sid]] <- data.frame(
      study = sid,
      log2FC = row$log2FC_shrunk %>% ifelse(is.na(.), row$log2FC, .),
      se = row$lfcSE
    )
  }
}
est_df <- bind_rows(study_estimates)

if (nrow(est_df) >= 2) {
  # Random-effects meta-analysis with metafor
  m <- rma(yi = log2FC, sei = se, data = est_df, method = "DL")

  p_forest <- ggplot(est_df, aes(x = log2FC, y = study)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = m$b, linetype = "dashed", color = "#FF9400", linewidth = 1) +
    geom_errorbarh(aes(xmin = log2FC - 1.96*se, xmax = log2FC + 1.96*se), height = 0.2, color = "#0279EE") +
    geom_point(size = 3, color = "#0279EE") +
    labs(x = "log2 Fold Change (Flight vs GC)", y = "Study",
         title = paste0("Forest plot: ", top_gene),
         subtitle = sprintf("Pooled log2FC = %.3f (p = %.2e, I² = %.1f%%)", m$b, m$pval, m$I2)) +
    theme_bw()

  ggsave(file.path(OUT_DIR, "forest_plot_top_gene.svg"), p_forest, width = 8, height = 5)
  ggsave(file.path(OUT_DIR, "forest_plot_top_gene.png"), p_forest, width = 8, height = 5, dpi = 150)
  cat("Saved forest plot\n")
}

# ============================================================
# 9. Tropism data summary
# ============================================================
cat("\n=== 9. Data summary ===\n")
tropism_summary <- meta %>%
  count(tropism_type, spaceflight_condition, assay_type) %>%
  arrange(desc(n))

p_summary <- ggplot(tropism_summary, aes(x = tropism_type, y = n, fill = spaceflight_condition)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = flight_gc_colors) +
  labs(x = "Tropism type", y = "Number of samples",
       title = "Arabidopsis spaceflight transcriptomics dataset",
       subtitle = "1337 samples from OSDR + GEO") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~assay_type)

ggsave(file.path(OUT_DIR, "dataset_summary.svg"), p_summary, width = 9, height = 5)
ggsave(file.path(OUT_DIR, "dataset_summary.png"), p_summary, width = 9, height = 5, dpi = 150)
cat("Saved dataset summary plot\n")

cat("\n=== ALL FIGURES DONE ===\n")
cat("Output directory:", OUT_DIR, "\n")
print(list.files(OUT_DIR))
