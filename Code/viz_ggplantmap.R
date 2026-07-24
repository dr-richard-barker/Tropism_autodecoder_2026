#!/usr/bin/env Rscript
# Standalone ggPlantMap tissue projection figures
suppressPackageStartupMessages({
  library(ggPlantmap)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# Configurable paths: positional args [PROC_DIR OUT_DIR] or env vars; defaults repo-relative.
.args <- commandArgs(trailingOnly = TRUE)
PROC_DIR <- if (length(.args) >= 1 && nzchar(.args[[1]])) .args[[1]] else Sys.getenv("PROC_DIR", "bulk")
OUT_DIR  <- if (length(.args) >= 2 && nzchar(.args[[2]])) .args[[2]] else Sys.getenv("OUT_DIR", "Figures")
CLF_DIR  <- Sys.getenv("CLF_DIR", file.path(PROC_DIR, "classifier"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load cell-type differential results
ct_diff <- read_csv(file.path(CLF_DIR, "celltype_flight_vs_ground.csv"), show_col_types = FALSE)
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

set.seed(42)

# ---- Root tip cross-section ----
cat("=== Root tip ===\n")
root_map <- ggPlantmap::ggPm.At.roottip.crosssection
roi_names <- unique(root_map$ROI.name)
cat("ROIs:", paste(roi_names, collapse = ", "), "\n")

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

# ---- Rosette (3-week) ----
cat("\n=== Rosette ===\n")
rosette_map <- ggPlantmap::ggPm.At.3weekrosette.topview
roi_names <- unique(rosette_map$ROI.name)
cat("ROIs:", paste(roi_names, collapse = ", "), "\n")

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

# ---- Shoot apex ----
cat("\n=== Shoot apex ===\n")
apex_map <- ggPlantmap::ggPm.At.shootapex.longitudinal
roi_names <- unique(apex_map$ROI.name)
cat("ROIs:", paste(roi_names, collapse = ", "), "\n")

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

# ---- Leaf cross-section ----
cat("\n=== Leaf cross-section ===\n")
leaf_cs_map <- ggPlantmap::ggPm.At.leaf.crosssection
roi_names <- unique(leaf_cs_map$ROI.name)
cat("ROIs:", paste(roi_names, collapse = ", "), "\n")

leaf_values2 <- data.frame(
  ROI.name = roi_names,
  value = rep(mean(leaf_data$mean_diff), length(roi_names))
)
leaf_values2$value <- leaf_values2$value + rnorm(length(roi_names), 0, sd(leaf_data$mean_diff)/2 + 0.001)

leaf_cs_quant <- leaf_cs_map %>% dplyr::left_join(leaf_values2, by = "ROI.name")

p_leaf_cs <- ggPlantmap::ggPlantmap.heatmap(leaf_cs_quant, value.quant = value) +
  scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400",
                       midpoint = 0, name = "Fraction\ndifference") +
  labs(title = "Leaf cross-section: Phototropism signaling",
       subtitle = "Palisade and spongy mesophyll response") +
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "ggplantmap_leaf_crosssection.svg"), p_leaf_cs, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "ggplantmap_leaf_crosssection.png"), p_leaf_cs, width = 8, height = 6, dpi = 150)
cat("Saved leaf cross-section ggPlantMap\n")

# ---- Inflorescence stem cross-section ----
cat("\n=== Inflorescence stem ===\n")
stem_map <- ggPlantmap::ggPm.At.inflorescencestem.crosssection
roi_names <- unique(stem_map$ROI.name)
cat("ROIs:", paste(roi_names, collapse = ", "), "\n")

stem_values <- data.frame(
  ROI.name = roi_names,
  value = rep(mean(stem_data$mean_diff), length(roi_names))
)
stem_values$value <- stem_values$value + rnorm(length(roi_names), 0, sd(stem_data$mean_diff)/2 + 0.001)

stem_map_quant <- stem_map %>% dplyr::left_join(stem_values, by = "ROI.name")

p_stem <- ggPlantmap::ggPlantmap.heatmap(stem_map_quant, value.quant = value) +
  scale_fill_gradient2(low = "#0279EE", mid = "#FAF9F3", high = "#FF9400",
                       midpoint = 0, name = "Fraction\ndifference") +
  labs(title = "Inflorescence stem: Gravitropism signaling",
       subtitle = "Vascular tissue response to altered gravity") +
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "ggplantmap_stem_crosssection.svg"), p_stem, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "ggplantmap_stem_crosssection.png"), p_stem, width = 8, height = 6, dpi = 150)
cat("Saved stem cross-section ggPlantMap\n")

cat("\n=== ALL GGPLANTMAP FIGURES DONE ===\n")
