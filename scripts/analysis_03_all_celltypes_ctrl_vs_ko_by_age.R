#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

project_dir <- normalizePath(getwd())
out_dir <- file.path(project_dir, "results", "analysis_03_all_celltypes_ctrl_vs_ko_by_age")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_rds <- file.path(project_dir, "scRNA_4_thomas", "photoreceptors_annotated.rds")
genes <- c("Gls", "Gls2")
cell_type_order <- c("Cones", "Rod", "MG", "AC", "BC", "HC", "RPE", "Vasculature cells")

sample_info <- tibble(
  sample = c("15dayS1", "15dayS2", "30dayS1", "30dayS2"),
  age = c("P15", "P15", "P35", "P35"),
  condition = c("Ctrl", "KO", "Ctrl", "KO")
)

message("Loading full annotated object: ", input_rds)
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

obj$sample_original <- as.character(obj$sample)
obj$cell_type <- as.character(Idents(obj))
obj@meta.data <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(sample_info, by = c("sample_original" = "sample")) %>%
  column_to_rownames("cell")
obj$condition <- factor(obj$condition, levels = c("Ctrl", "KO"))
obj$age <- factor(obj$age, levels = c("P15", "P35"))
obj$cell_type <- factor(obj$cell_type, levels = cell_type_order)
obj$age_condition <- factor(
  paste(obj$age, obj$condition, sep = "_"),
  levels = c("P15_Ctrl", "P15_KO", "P35_Ctrl", "P35_KO")
)

genes_present <- intersect(genes, rownames(obj))
if (length(genes_present) == 0) {
  stop("None of the requested genes were found: ", paste(genes, collapse = ", "))
}

message("Summarizing cell composition...")
cell_counts <- obj@meta.data %>%
  count(age, condition, sample_original, cell_type, name = "n_cells") %>%
  group_by(age, condition, sample_original) %>%
  mutate(total_cells = sum(n_cells), pct_of_sample = n_cells / total_cells * 100) %>%
  ungroup() %>%
  arrange(age, condition, cell_type)
write_csv(cell_counts, file.path(out_dir, "cell_counts_and_percent_by_age_condition_cell_type.csv"))

composition_plot <- ggplot(cell_counts, aes(x = condition, y = pct_of_sample, fill = cell_type)) +
  geom_col(width = 0.75) +
  facet_wrap(~ age, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(
    x = NULL,
    y = "% of cells",
    fill = "Cell type",
    title = "Cell type composition by age and condition"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "cell_type_composition_by_age_condition.pdf"), composition_plot, width = 8, height = 5)
ggsave(file.path(out_dir, "cell_type_composition_by_age_condition.png"), composition_plot, width = 8, height = 5, dpi = 300)

message("Extracting Gls/Gls2 expression...")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")[genes_present, , drop = FALSE]
counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes_present, , drop = FALSE]
meta <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, sample_original, age, condition, age_condition, cell_type)

expr_long <- as.data.frame(t(as.matrix(data_mat))) %>%
  rownames_to_column("cell") %>%
  pivot_longer(all_of(genes_present), names_to = "gene", values_to = "lognorm_expr") %>%
  left_join(meta, by = "cell")

counts_long <- as.data.frame(t(as.matrix(counts_mat))) %>%
  rownames_to_column("cell") %>%
  pivot_longer(all_of(genes_present), names_to = "gene", values_to = "raw_count") %>%
  select(cell, gene, raw_count)

expr_long <- expr_long %>%
  left_join(counts_long, by = c("cell", "gene")) %>%
  mutate(detected = raw_count > 0)

summary_by_group <- expr_long %>%
  group_by(gene, age, condition, sample_original, cell_type) %>%
  summarise(
    n_cells = n(),
    pct_detected = mean(detected) * 100,
    avg_lognorm_expr = mean(lognorm_expr),
    median_lognorm_expr = median(lognorm_expr),
    .groups = "drop"
  ) %>%
  arrange(gene, age, cell_type, condition)

write_csv(expr_long, file.path(out_dir, "all_celltypes_gls_gls2_single_cell_expression.csv"))
write_csv(summary_by_group, file.path(out_dir, "all_celltypes_gls_gls2_summary_by_age_condition_cell_type.csv"))

message("Running exploratory cell-level Ctrl vs KO tests within age and cell type...")
stats <- expr_long %>%
  group_by(gene, age, cell_type) %>%
  group_modify(function(.x, .y) {
    ctrl <- .x %>% filter(condition == "Ctrl") %>% pull(lognorm_expr)
    ko <- .x %>% filter(condition == "KO") %>% pull(lognorm_expr)
    if (length(ctrl) == 0 || length(ko) == 0) {
      return(tibble())
    }
    if (all(ctrl == ko[1]) && length(unique(c(ctrl, ko))) == 1) {
      p_value <- NA_real_
    } else {
      p_value <- suppressWarnings(wilcox.test(ctrl, ko)$p.value)
    }
    tibble(
      ctrl_n = length(ctrl),
      ko_n = length(ko),
      ctrl_avg_lognorm_expr = mean(ctrl),
      ko_avg_lognorm_expr = mean(ko),
      difference_ko_minus_ctrl = mean(ko) - mean(ctrl),
      ctrl_pct_detected = mean(ctrl > 0) * 100,
      ko_pct_detected = mean(ko > 0) * 100,
      p_value = p_value
    )
  }) %>%
  ungroup() %>%
  group_by(gene, age) %>%
  mutate(p_adj_bh_within_age_gene = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(gene, age, cell_type)

write_csv(stats, file.path(out_dir, "all_celltypes_gls_gls2_ctrl_vs_ko_by_age_cell_type_wilcox_exploratory.csv"))

dot_plot <- ggplot(summary_by_group, aes(x = condition, y = cell_type)) +
  geom_point(aes(size = pct_detected, color = avg_lognorm_expr)) +
  facet_grid(gene ~ age) +
  scale_color_viridis_c(option = "magma", name = "Avg log-normalized\nexpression") +
  scale_size_continuous(name = "% cells detected", range = c(1.2, 8), limits = c(0, 100)) +
  labs(
    x = NULL,
    y = NULL,
    title = "Gls/Gls2 expression across cell types: Ctrl vs KO by age",
    subtitle = "Full annotated object: photoreceptors_annotated.rds"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "all_celltypes_gls_gls2_dotplot_ctrl_vs_ko_by_age_cell_type.pdf"), dot_plot, width = 9, height = 7)
ggsave(file.path(out_dir, "all_celltypes_gls_gls2_dotplot_ctrl_vs_ko_by_age_cell_type.png"), dot_plot, width = 9, height = 7, dpi = 300)

gls_summary <- summary_by_group %>% filter(gene == "Gls")
gls_heatmap <- ggplot(gls_summary, aes(x = condition, y = cell_type, fill = avg_lognorm_expr)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f%%", pct_detected)), size = 3) +
  facet_wrap(~ age, nrow = 1) +
  scale_fill_viridis_c(option = "magma", name = "Avg log-normalized\nexpression") +
  labs(
    x = NULL,
    y = NULL,
    title = "Gls expression across cell types",
    subtitle = "Tile text shows % cells detected"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "all_celltypes_gls_heatmap_ctrl_vs_ko_by_age_cell_type.pdf"), gls_heatmap, width = 8, height = 5.5)
ggsave(file.path(out_dir, "all_celltypes_gls_heatmap_ctrl_vs_ko_by_age_cell_type.png"), gls_heatmap, width = 8, height = 5.5, dpi = 300)

feature_plot <- FeaturePlot(
  obj,
  features = genes_present,
  split.by = "age_condition",
  reduction = "umap",
  order = TRUE,
  pt.size = 0.15,
  ncol = 4
) +
  plot_annotation(title = "Full object: Gls/Gls2 on UMAP by age and condition")

ggsave(file.path(out_dir, "all_celltypes_gls_gls2_featureplot_by_age_condition.pdf"), feature_plot, width = 15, height = 8)
ggsave(file.path(out_dir, "all_celltypes_gls_gls2_featureplot_by_age_condition.png"), feature_plot, width = 15, height = 8, dpi = 300)

message("Done. Outputs written to: ", out_dir)
