#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
})

project_dir <- normalizePath(getwd())
input_rds <- file.path(project_dir, "scRNA_4_thomas", "photoreceptors_annotated.rds")
out_dir <- file.path(project_dir, "results", "analysis_01_wt_gls_expression")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes_requested <- c("Gls", "Gls2")
ctrl_samples <- c("15dayS1", "30dayS1")

message("Loading: ", input_rds)
obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

obj$cell_type <- as.character(Idents(obj))
obj$condition <- dplyr::recode(
  obj$sample,
  "15dayS1" = "Ctrl",
  "30dayS1" = "Ctrl",
  "15dayS2" = "KO",
  "30dayS2" = "KO",
  .default = as.character(obj$sample)
)
obj$age <- dplyr::case_when(
  obj$sample %in% c("15dayS1", "15dayS2") ~ "P15",
  obj$sample %in% c("30dayS1", "30dayS2") ~ "P35",
  TRUE ~ NA_character_
)

genes_present <- intersect(genes_requested, rownames(obj))
genes_missing <- setdiff(genes_requested, genes_present)
if (length(genes_present) == 0) {
  stop("None of the requested genes were found: ", paste(genes_requested, collapse = ", "))
}
if (length(genes_missing) > 0) {
  warning("Missing requested genes: ", paste(genes_missing, collapse = ", "))
}

wt <- subset(obj, subset = sample %in% ctrl_samples)
cell_type_order <- c("Cones", "Rod", "MG", "AC", "BC", "HC", "RPE", "Vasculature cells")
wt$cell_type <- factor(wt$cell_type, levels = cell_type_order)

cell_counts <- wt@meta.data %>%
  count(age, sample, cell_type, name = "n_cells") %>%
  arrange(age, sample, cell_type)
write_csv(cell_counts, file.path(out_dir, "wt_cell_counts_by_age_cell_type.csv"))

data_mat <- GetAssayData(wt, assay = "RNA", layer = "data")[genes_present, , drop = FALSE]
counts_mat <- GetAssayData(wt, assay = "RNA", layer = "counts")[genes_present, , drop = FALSE]
meta <- wt@meta.data %>%
  mutate(cell = rownames(wt@meta.data)) %>%
  select(cell, sample, age, condition, cell_type)

expr_long <- as.data.frame(t(as.matrix(data_mat))) %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(all_of(genes_present), names_to = "gene", values_to = "lognorm_expr") %>%
  left_join(meta, by = "cell")

counts_long <- as.data.frame(t(as.matrix(counts_mat))) %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(all_of(genes_present), names_to = "gene", values_to = "raw_count") %>%
  select(cell, gene, raw_count)

expr_long <- expr_long %>%
  left_join(counts_long, by = c("cell", "gene")) %>%
  mutate(detected = raw_count > 0)

summary_by_cell_type <- expr_long %>%
  group_by(gene, cell_type) %>%
  summarise(
    n_cells = n(),
    pct_detected = mean(detected) * 100,
    avg_lognorm_expr = mean(lognorm_expr),
    median_lognorm_expr = median(lognorm_expr),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(avg_lognorm_expr))

summary_by_age_cell_type <- expr_long %>%
  group_by(gene, age, sample, cell_type) %>%
  summarise(
    n_cells = n(),
    pct_detected = mean(detected) * 100,
    avg_lognorm_expr = mean(lognorm_expr),
    median_lognorm_expr = median(lognorm_expr),
    .groups = "drop"
  ) %>%
  arrange(gene, age, desc(avg_lognorm_expr))

write_csv(summary_by_cell_type, file.path(out_dir, "wt_gls_gls2_expression_by_cell_type.csv"))
write_csv(summary_by_age_cell_type, file.path(out_dir, "wt_gls_gls2_expression_by_age_cell_type.csv"))
write_csv(expr_long, file.path(out_dir, "wt_gls_gls2_single_cell_expression.csv"))

wilcox_results <- lapply(genes_present, function(gene_i) {
  gene_df <- expr_long %>% filter(gene == gene_i)
  cone_values <- gene_df %>% filter(cell_type == "Cones") %>% pull(lognorm_expr)
  lapply(setdiff(levels(wt$cell_type), "Cones"), function(cell_type_i) {
    other_values <- gene_df %>% filter(cell_type == cell_type_i) %>% pull(lognorm_expr)
    if (length(cone_values) == 0 || length(other_values) == 0) {
      return(NULL)
    }
    test <- suppressWarnings(wilcox.test(cone_values, other_values))
    data.frame(
      gene = gene_i,
      comparison = paste("Cones vs", cell_type_i),
      cones_n = length(cone_values),
      other_cell_type = cell_type_i,
      other_n = length(other_values),
      cones_mean = mean(cone_values),
      other_mean = mean(other_values),
      mean_difference = mean(cone_values) - mean(other_values),
      p_value = test$p.value
    )
  }) %>% bind_rows()
}) %>%
  bind_rows() %>%
  group_by(gene) %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write_csv(wilcox_results, file.path(out_dir, "wt_gls_gls2_cones_vs_other_cell_types_wilcox_exploratory.csv"))

format_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), sprintf("%.3f", x))
  )
}

dot_df <- summary_by_cell_type %>%
  mutate(cell_type = factor(cell_type, levels = levels(wt$cell_type)))

dot_plot <- ggplot(dot_df, aes(x = cell_type, y = gene)) +
  geom_point(aes(size = pct_detected, color = avg_lognorm_expr)) +
  scale_color_viridis_c(option = "magma", name = "Avg log-normalized\nexpression") +
  scale_size_continuous(name = "% cells detected", range = c(1.5, 9), limits = c(0, 100)) +
  labs(x = "WT/Ctrl cell type", y = NULL, title = "WT expression of Gls and Gls2 across retinal cell types") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "wt_gls_gls2_dotplot_by_cell_type.pdf"), dot_plot, width = 8.5, height = 4.5)
ggsave(file.path(out_dir, "wt_gls_gls2_dotplot_by_cell_type.png"), dot_plot, width = 8.5, height = 4.5, dpi = 300)

age_dot_df <- summary_by_age_cell_type %>%
  mutate(cell_type = factor(cell_type, levels = levels(wt$cell_type)))

age_dot_plot <- ggplot(age_dot_df, aes(x = cell_type, y = gene)) +
  geom_point(aes(size = pct_detected, color = avg_lognorm_expr)) +
  facet_wrap(~ age, ncol = 1) +
  scale_color_viridis_c(option = "magma", name = "Avg log-normalized\nexpression") +
  scale_size_continuous(name = "% cells detected", range = c(1.5, 9), limits = c(0, 100)) +
  labs(x = "WT/Ctrl cell type", y = NULL, title = "WT expression split by age/sample") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "wt_gls_gls2_dotplot_by_age_cell_type.pdf"), age_dot_plot, width = 8.5, height = 6.5)
ggsave(file.path(out_dir, "wt_gls_gls2_dotplot_by_age_cell_type.png"), age_dot_plot, width = 8.5, height = 6.5, dpi = 300)

violin_plot <- ggplot(expr_long, aes(x = cell_type, y = lognorm_expr, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.1, alpha = 0.6) +
  facet_wrap(~ gene, ncol = 1, scales = "free_y") +
  labs(x = "WT/Ctrl cell type", y = "Log-normalized expression", title = "WT single-cell expression") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, "wt_gls_gls2_violin_by_cell_type.pdf"), violin_plot, width = 8.5, height = 7)
ggsave(file.path(out_dir, "wt_gls_gls2_violin_by_cell_type.png"), violin_plot, width = 8.5, height = 7, dpi = 300)

for (gene_i in genes_present) {
  gene_summary <- summary_by_age_cell_type %>%
    filter(gene == gene_i) %>%
    mutate(cell_type = factor(cell_type, levels = rev(cell_type_order)))

  gene_plot <- ggplot(gene_summary, aes(x = age, y = cell_type)) +
    geom_point(
      aes(size = pct_detected, color = avg_lognorm_expr),
      stroke = 0.35
    ) +
    scale_color_viridis_c(option = "viridis", name = "Avg log-normalized\nexpression") +
    scale_size_continuous(name = "% cells detected", range = c(3, 12), limits = c(0, 100)) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("WT/Ctrl ", gene_i, " expression by retinal cell type"),
      subtitle = "Cones are shown at the top; dot size is detection rate"
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )

  ggsave(
    file.path(out_dir, paste0("moloy_review_wt_", tolower(gene_i), "_dotplot_by_age_cell_type.pdf")),
    gene_plot,
    width = 8.5,
    height = 6.2
  )
  ggsave(
    file.path(out_dir, paste0("moloy_review_wt_", tolower(gene_i), "_dotplot_by_age_cell_type.png")),
    gene_plot,
    width = 8.5,
    height = 6.2,
    dpi = 300
  )
}

gls_cone_stats <- wilcox_results %>%
  filter(gene == "Gls") %>%
  mutate(
    other_cell_type = factor(other_cell_type, levels = rev(setdiff(cell_type_order, "Cones"))),
    fdr_label = paste0("FDR=", format_p(p_adj_bh))
  )

if (nrow(gls_cone_stats) > 0) {
  gls_stats_plot <- ggplot(gls_cone_stats, aes(x = mean_difference, y = other_cell_type)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_col(fill = "#2C7BB6", width = 0.65) +
    geom_text(aes(label = fdr_label), hjust = -0.05, size = 5) +
    labs(
      x = "Mean log-normalized expression difference (Cones - other cell type)",
      y = NULL,
      title = "WT/Ctrl Gls: cones compared with other retinal cell types",
      subtitle = "Exploratory Wilcoxon tests; labels show BH FDR"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 18) +
    theme(plot.margin = margin(5.5, 85, 5.5, 5.5))

  ggsave(file.path(out_dir, "moloy_review_wt_gls_cones_vs_other_fdr_barplot.pdf"), gls_stats_plot, width = 10, height = 5.8)
  ggsave(file.path(out_dir, "moloy_review_wt_gls_cones_vs_other_fdr_barplot.png"), gls_stats_plot, width = 10, height = 5.8, dpi = 300)
}

message("Done. Outputs written to: ", out_dir)
