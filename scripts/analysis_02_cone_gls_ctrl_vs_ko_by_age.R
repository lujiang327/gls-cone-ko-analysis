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
out_dir <- file.path(project_dir, "results", "analysis_02_cone_gls_ctrl_vs_ko_by_age")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

primary_rds <- file.path(project_dir, "scRNA_4_thomas", "photoreceptors_reClusteredCones.rds")
sensitivity_rds <- file.path(project_dir, "scRNA_4_thomas", "photoreceptors_SubsetCones.rds")
genes <- c("Gls", "Gls2")
condition_colors <- c("Ctrl" = "#2C7BB6", "KO" = "#D7191C")

sample_info <- tibble(
  sample = c("15dayS1", "15dayS2", "30dayS1", "30dayS2"),
  age = c("P15", "P15", "P35", "P35"),
  condition = c("Ctrl", "KO", "Ctrl", "KO")
)

add_design <- function(obj) {
  obj$sample_original <- as.character(obj$sample)
  obj@meta.data <- obj@meta.data %>%
    rownames_to_column("cell") %>%
    left_join(sample_info, by = c("sample_original" = "sample")) %>%
    column_to_rownames("cell")
  obj$age_condition <- factor(
    paste(obj$age, obj$condition, sep = "_"),
    levels = c("P15_Ctrl", "P15_KO", "P35_Ctrl", "P35_KO")
  )
  obj$condition <- factor(obj$condition, levels = c("Ctrl", "KO"))
  obj$age <- factor(obj$age, levels = c("P15", "P35"))
  obj
}

summarize_gene_expression <- function(obj, object_label) {
  DefaultAssay(obj) <- "RNA"
  genes_present <- intersect(genes, rownames(obj))
  if (length(genes_present) == 0) {
    stop("No requested genes found in ", object_label)
  }

  data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")[genes_present, , drop = FALSE]
  counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes_present, , drop = FALSE]
  meta <- obj@meta.data %>%
    rownames_to_column("cell") %>%
    select(cell, sample_original, age, condition, age_condition)

  expr_long <- as.data.frame(t(as.matrix(data_mat))) %>%
    rownames_to_column("cell") %>%
    pivot_longer(all_of(genes_present), names_to = "gene", values_to = "lognorm_expr") %>%
    left_join(meta, by = "cell")

  counts_long <- as.data.frame(t(as.matrix(counts_mat))) %>%
    rownames_to_column("cell") %>%
    pivot_longer(all_of(genes_present), names_to = "gene", values_to = "raw_count") %>%
    select(cell, gene, raw_count)

  expr_long %>%
    left_join(counts_long, by = c("cell", "gene")) %>%
    mutate(
      detected = raw_count > 0,
      object = object_label
    )
}

age_wilcox <- function(expr_long) {
  expr_long %>%
    group_by(object, gene, age) %>%
    group_modify(function(.x, .y) {
      ctrl <- .x %>% filter(condition == "Ctrl") %>% pull(lognorm_expr)
      ko <- .x %>% filter(condition == "KO") %>% pull(lognorm_expr)
      if (length(ctrl) == 0 || length(ko) == 0) {
        return(tibble())
      }
      test <- suppressWarnings(wilcox.test(ctrl, ko))
      tibble(
        ctrl_n = length(ctrl),
        ko_n = length(ko),
        ctrl_avg_lognorm_expr = mean(ctrl),
        ko_avg_lognorm_expr = mean(ko),
        difference_ko_minus_ctrl = mean(ko) - mean(ctrl),
        ctrl_pct_detected = mean(ctrl > 0) * 100,
        ko_pct_detected = mean(ko > 0) * 100,
        p_value = test$p.value
      )
    }) %>%
    ungroup() %>%
    group_by(object, gene) %>%
    mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
    ungroup()
}

format_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x))
  )
}

message("Loading primary cone object: ", primary_rds)
cones <- readRDS(primary_rds) %>% add_design()
DefaultAssay(cones) <- "RNA"

message("Loading sensitivity cone object: ", sensitivity_rds)
cones_all <- readRDS(sensitivity_rds) %>% add_design()
DefaultAssay(cones_all) <- "RNA"

cell_counts <- bind_rows(
  cones@meta.data %>%
    count(age, condition, sample_original, name = "n_cells") %>%
    mutate(object = "reclustered_cones"),
  cones_all@meta.data %>%
    count(age, condition, sample_original, name = "n_cells") %>%
    mutate(object = "all_initial_cones")
) %>%
  select(object, age, condition, sample_original, n_cells) %>%
  arrange(object, age, condition)
write_csv(cell_counts, file.path(out_dir, "cone_cell_counts_by_object_age_condition.csv"))

primary_expr <- summarize_gene_expression(cones, "reclustered_cones")
sensitivity_expr <- summarize_gene_expression(cones_all, "all_initial_cones")
combined_expr <- bind_rows(primary_expr, sensitivity_expr)

summary_by_group <- combined_expr %>%
  group_by(object, gene, age, condition, sample_original) %>%
  summarise(
    n_cells = n(),
    pct_detected = mean(detected) * 100,
    avg_lognorm_expr = mean(lognorm_expr),
    median_lognorm_expr = median(lognorm_expr),
    .groups = "drop"
  ) %>%
  arrange(object, gene, age, condition)

stats <- age_wilcox(combined_expr)

write_csv(primary_expr, file.path(out_dir, "reclustered_cones_gls_gls2_single_cell_expression.csv"))
write_csv(summary_by_group, file.path(out_dir, "cone_gls_gls2_summary_by_age_condition.csv"))
write_csv(stats, file.path(out_dir, "cone_gls_gls2_ctrl_vs_ko_by_age_wilcox_exploratory.csv"))

primary_summary <- summary_by_group %>%
  filter(object == "reclustered_cones")

dot_plot <- ggplot(primary_summary, aes(x = condition, y = gene)) +
  geom_point(aes(size = pct_detected, color = avg_lognorm_expr)) +
  facet_wrap(~ age, nrow = 1) +
  scale_color_viridis_c(option = "magma", name = "Avg log-normalized\nexpression") +
  scale_size_continuous(name = "% cells detected", range = c(1.5, 9), limits = c(0, 100)) +
  labs(
    x = NULL,
    y = NULL,
    title = "Cone Gls/Gls2 expression: Ctrl vs KO by age",
    subtitle = "Primary object: photoreceptors_reClusteredCones.rds"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "cone_gls_gls2_dotplot_ctrl_vs_ko_by_age.pdf"), dot_plot, width = 7, height = 4.5)
ggsave(file.path(out_dir, "cone_gls_gls2_dotplot_ctrl_vs_ko_by_age.png"), dot_plot, width = 7, height = 4.5, dpi = 300)

violin_plot <- primary_expr %>%
  mutate(condition = factor(condition, levels = c("Ctrl", "KO"))) %>%
  ggplot(aes(x = condition, y = lognorm_expr, fill = condition)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.65) +
  facet_grid(gene ~ age, scales = "free_y") +
  scale_fill_manual(values = c("Ctrl" = "#2C7BB6", "KO" = "#D7191C")) +
  labs(
    x = NULL,
    y = "Log-normalized expression",
    title = "Cone Gls/Gls2 expression: Ctrl vs KO by age",
    subtitle = "Cell-level distributions; p-values are exploratory because each age has one Ctrl and one KO sample"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "cone_gls_gls2_violin_ctrl_vs_ko_by_age.pdf"), violin_plot, width = 8, height = 6.5)
ggsave(file.path(out_dir, "cone_gls_gls2_violin_ctrl_vs_ko_by_age.png"), violin_plot, width = 8, height = 6.5, dpi = 300)

feature_genes <- intersect(genes, rownames(cones))
feature_plot <- FeaturePlot(
  cones,
  features = feature_genes,
  split.by = "age_condition",
  reduction = "umap",
  order = TRUE,
  pt.size = 0.8,
  ncol = 4
) +
  plot_annotation(title = "Reclustered cones: Gls/Gls2 on UMAP by age and condition")

ggsave(file.path(out_dir, "cone_gls_gls2_featureplot_by_age_condition.pdf"), feature_plot, width = 14, height = 7)
ggsave(file.path(out_dir, "cone_gls_gls2_featureplot_by_age_condition.png"), feature_plot, width = 14, height = 7, dpi = 300)

sample_umap <- DimPlot(
  cones,
  reduction = "umap",
  group.by = "age_condition",
  pt.size = 0.8
) +
  ggtitle("Reclustered cones: age and condition") +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "cone_umap_age_condition.pdf"), sample_umap, width = 8, height = 6)
ggsave(file.path(out_dir, "cone_umap_age_condition.png"), sample_umap, width = 8, height = 6, dpi = 300)

primary_stats <- stats %>%
  filter(object == "reclustered_cones") %>%
  mutate(
    age = factor(age, levels = c("P15", "P35")),
    stat_label = paste0("p=", format_p(p_value), "\nFDR=", format_p(p_adj_bh))
  )

for (gene_i in genes) {
  gene_expr <- primary_expr %>%
    filter(gene == gene_i) %>%
    mutate(
      age = factor(age, levels = c("P15", "P35")),
      condition = factor(condition, levels = c("Ctrl", "KO"))
    )
  gene_stats <- primary_stats %>% filter(gene == gene_i)
  y_max <- gene_expr %>%
    group_by(age) %>%
    summarise(y = max(lognorm_expr, na.rm = TRUE) * 1.08 + 0.05, .groups = "drop")
  gene_stats <- gene_stats %>% left_join(y_max, by = "age")

  gene_violin <- ggplot(gene_expr, aes(x = condition, y = lognorm_expr, fill = condition)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.75) +
    geom_boxplot(width = 0.14, outlier.size = 0.25, alpha = 0.85) +
    geom_text(
      data = gene_stats,
      aes(x = 1.5, y = y, label = stat_label),
      inherit.aes = FALSE,
      size = 5.2,
      lineheight = 0.95
    ) +
    facet_wrap(~ age, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = condition_colors) +
    labs(
      x = NULL,
      y = "Log-normalized expression",
      title = paste0("Cones: ", gene_i, " Ctrl vs KO by age"),
      subtitle = "Primary reclustered cone object; exploratory cell-level Wilcoxon tests"
    ) +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(out_dir, paste0("moloy_review_cones_", tolower(gene_i), "_violin_ctrl_vs_ko_by_age_with_stats.pdf")),
    gene_violin,
    width = 10,
    height = 5.8
  )
  ggsave(
    file.path(out_dir, paste0("moloy_review_cones_", tolower(gene_i), "_violin_ctrl_vs_ko_by_age_with_stats.png")),
    gene_violin,
    width = 10,
    height = 5.8,
    dpi = 300
  )
}

count_labels <- cell_counts %>%
  filter(object == "reclustered_cones") %>%
  mutate(
    age = factor(age, levels = c("P15", "P35")),
    condition = factor(condition, levels = c("Ctrl", "KO")),
    label = paste0("n=", n_cells)
  )

count_plot <- ggplot(count_labels, aes(x = condition, y = n_cells, fill = condition)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_text(aes(label = label), vjust = -0.35, size = 6) +
  facet_wrap(~ age, nrow = 1) +
  scale_fill_manual(values = condition_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = NULL,
    y = "Number of reclustered cone cells",
    title = "Recovered cone cells used for cone-only analysis"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "moloy_review_cone_cell_counts_ctrl_vs_ko_by_age.pdf"), count_plot, width = 9.5, height = 5.2)
ggsave(file.path(out_dir, "moloy_review_cone_cell_counts_ctrl_vs_ko_by_age.png"), count_plot, width = 9.5, height = 5.2, dpi = 300)

message("Done. Outputs written to: ", out_dir)
