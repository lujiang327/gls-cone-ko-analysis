#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
})

project_dir <- normalizePath(getwd())
out_dir <- file.path(project_dir, "results", "umap_review")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

condition_from_sample <- function(sample) {
  dplyr::recode(
    sample,
    "15dayS1" = "Ctrl",
    "30dayS1" = "Ctrl",
    "15dayS2" = "KO",
    "30dayS2" = "KO",
    .default = as.character(sample)
  )
}

age_from_sample <- function(sample) {
  dplyr::case_when(
    sample %in% c("15dayS1", "15dayS2") ~ "P15",
    sample %in% c("30dayS1", "30dayS2") ~ "P35",
    TRUE ~ NA_character_
  )
}

save_plot_pair <- function(plot, stem, width = 10, height = 7) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300)
}

message("Loading full annotated object...")
full <- readRDS(file.path(project_dir, "scRNA_4_thomas", "photoreceptors_annotated.rds"))
DefaultAssay(full) <- "RNA"
full$cell_type <- as.character(Idents(full))
full$condition <- condition_from_sample(full$sample)
full$age <- age_from_sample(full$sample)

full_counts <- full@meta.data %>%
  count(age, condition, sample, cell_type, name = "n_cells") %>%
  arrange(age, condition, sample, cell_type)
write_csv(full_counts, file.path(out_dir, "full_object_cell_counts.csv"))

message("Writing full-object UMAPs...")
p_full_cell_type <- DimPlot(
  full,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2
) +
  ggtitle("Full photoreceptor object: annotated cell types") +
  theme_classic(base_size = 12)

p_full_sample <- DimPlot(
  full,
  reduction = "umap",
  group.by = "sample",
  pt.size = 0.2
) +
  ggtitle("Full photoreceptor object: sample") +
  theme_classic(base_size = 12)

p_full_condition <- DimPlot(
  full,
  reduction = "umap",
  group.by = "condition",
  pt.size = 0.2
) +
  ggtitle("Full photoreceptor object: Ctrl vs KO") +
  theme_classic(base_size = 12)

p_full_age <- DimPlot(
  full,
  reduction = "umap",
  group.by = "age",
  pt.size = 0.2
) +
  ggtitle("Full photoreceptor object: age") +
  theme_classic(base_size = 12)

save_plot_pair(p_full_cell_type, "full_umap_cell_type", width = 10, height = 7)
save_plot_pair(p_full_sample, "full_umap_sample", width = 10, height = 7)
save_plot_pair(p_full_condition, "full_umap_condition", width = 8, height = 6)
save_plot_pair(p_full_age, "full_umap_age", width = 8, height = 6)

pdf(file.path(out_dir, "full_umap_review_multipage.pdf"), width = 11, height = 8)
print(p_full_cell_type)
print(p_full_sample)
print(p_full_condition)
print(p_full_age)
if (all(c("Gls", "Gls2") %in% rownames(full))) {
  print(FeaturePlot(full, features = c("Gls", "Gls2"), reduction = "umap", order = TRUE, pt.size = 0.2) +
    plot_annotation(title = "Full photoreceptor object: Gls / Gls2"))
}
dev.off()

message("Loading cone object...")
cones <- readRDS(file.path(project_dir, "scRNA_4_thomas", "photoreceptors_reClusteredCones.rds"))
DefaultAssay(cones) <- "RNA"
cones$condition <- condition_from_sample(cones$sample)
cones$age <- age_from_sample(cones$sample)
cones$sample_label <- as.character(cones$sample)

cone_counts <- cones@meta.data %>%
  count(age, condition, sample, name = "n_cells") %>%
  arrange(age, condition, sample)
write_csv(cone_counts, file.path(out_dir, "cone_object_cell_counts.csv"))

message("Writing cone-object UMAPs...")
p_cones_sample <- DimPlot(
  cones,
  reduction = "umap",
  group.by = "sample",
  pt.size = 0.8
) +
  ggtitle("Reclustered cones: sample") +
  theme_classic(base_size = 12)

p_cones_condition <- DimPlot(
  cones,
  reduction = "umap",
  group.by = "condition",
  pt.size = 0.8
) +
  ggtitle("Reclustered cones: Ctrl vs KO") +
  theme_classic(base_size = 12)

p_cones_age <- DimPlot(
  cones,
  reduction = "umap",
  group.by = "age",
  pt.size = 0.8
) +
  ggtitle("Reclustered cones: age") +
  theme_classic(base_size = 12)

save_plot_pair(p_cones_sample, "cones_umap_sample", width = 8, height = 6)
save_plot_pair(p_cones_condition, "cones_umap_condition", width = 8, height = 6)
save_plot_pair(p_cones_age, "cones_umap_age", width = 8, height = 6)

pdf(file.path(out_dir, "cones_umap_review_multipage.pdf"), width = 10, height = 7)
print(p_cones_sample)
print(p_cones_condition)
print(p_cones_age)
feature_genes <- intersect(c("Gls", "Gls2", "Arr3", "Gnat2", "Opn1mw", "Opn1sw"), rownames(cones))
if (length(feature_genes) > 0) {
  print(FeaturePlot(cones, features = feature_genes, reduction = "umap", order = TRUE, pt.size = 0.8) +
    plot_annotation(title = "Reclustered cones: marker/features"))
}
dev.off()

message("Done. Outputs written to: ", out_dir)
