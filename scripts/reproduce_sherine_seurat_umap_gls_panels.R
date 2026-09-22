#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

project_dir <- normalizePath(getwd(), mustWork = TRUE)

in_rds <- file.path(project_dir, "scRNA_4_thomas", "photoreceptors_annotated.rds")
out_dir <- file.path(project_dir, "results", "sherine_style_umap_gls_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_map <- data.frame(
  sample = c("15dayS1", "15dayS2", "30dayS1", "30dayS2"),
  age = c("P15", "P15", "P35", "P35"),
  condition = c("WT", "cKO", "WT", "cKO"),
  stringsAsFactors = FALSE
)

celltype_order <- c(
  "Rod", "BC", "Cones", "AC", "MG", "HC", "Vasculature cells", "RPE"
)
celltype_labels <- c(
  Rod = "Rod",
  BC = "Bipolar",
  Cones = "Cone",
  AC = "Amacrine",
  MG = "Muller glia",
  HC = "Horizontal",
  `Vasculature cells` = "V/E cells",
  RPE = "RPE"
)
celltype_colors <- c(
  Rod = "#F8766D",
  BC = "#E69F00",
  Cones = "#7CAE00",
  AC = "#00BA38",
  MG = "#00BFC4",
  HC = "#C77CFF",
  `Vasculature cells` = "#619CFF",
  RPE = "#F564E3"
)

message("Reading ", in_rds)
obj <- readRDS(in_rds)
DefaultAssay(obj) <- "RNA"

umap <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap) <- c("UMAP_1", "UMAP_2")

expr <- FetchData(obj, vars = "Gls", layer = "data")[, 1]
df <- cbind(
  umap,
  data.frame(
    sample = as.character(obj$sample),
    celltype = as.character(Idents(obj)),
    Gls = as.numeric(expr),
    stringsAsFactors = FALSE
  )
)
df <- merge(df, sample_map, by = "sample", all.x = TRUE, sort = FALSE)
df$celltype <- factor(df$celltype, levels = celltype_order)

cell_counts <- as.data.frame(table(df$age, df$condition, df$celltype), stringsAsFactors = FALSE)
colnames(cell_counts) <- c("age", "condition", "celltype", "n_cells")
write.csv(cell_counts, file.path(out_dir, "seurat_cells_used_for_sherine_style_panels.csv"), row.names = FALSE)

cone_box <- function(plot_df, pad = 0.25) {
  cones <- plot_df[plot_df$celltype == "Cones", ]
  if (nrow(cones) == 0) return(NULL)
  xr <- as.numeric(quantile(cones$UMAP_1, c(0.05, 0.95), na.rm = TRUE))
  yr <- as.numeric(quantile(cones$UMAP_2, c(0.05, 0.95), na.rm = TRUE))
  data.frame(
    xmin = xr[1] - pad,
    xmax = xr[2] + pad,
    ymin = yr[1] - pad,
    ymax = yr[2] + pad
  )
}

base_umap_theme <- function() {
  theme_classic(base_size = 13) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 11),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

plot_celltypes <- function(plot_df) {
  labels <- aggregate(cbind(UMAP_1, UMAP_2) ~ celltype, data = plot_df, FUN = median)
  labels$label <- unname(celltype_labels[as.character(labels$celltype)])
  box <- cone_box(plot_df)

  ggplot(plot_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = celltype), size = 0.25, alpha = 0.9, stroke = 0) +
    geom_text(
      data = labels,
      aes(x = UMAP_1, y = UMAP_2, label = label),
      color = "black",
      size = 3.1,
      inherit.aes = FALSE
    ) +
    geom_rect(
      data = box,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = NA, color = "#0072B2", linewidth = 0.8
    ) +
    annotate("text", x = mean(c(box$xmin, box$xmax)), y = mean(c(box$ymin, box$ymax)),
             label = "Cone", size = 3.3, fontface = "bold") +
    scale_color_manual(
      values = celltype_colors,
      breaks = names(celltype_colors),
      labels = unname(celltype_labels[names(celltype_colors)]),
      drop = TRUE
    ) +
    coord_fixed() +
    labs(x = "UMAP_1", y = "UMAP_2", color = NULL) +
    base_umap_theme() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.key.height = unit(0.38, "cm"),
      legend.key.width = unit(0.38, "cm")
    )
}

plot_gls <- function(plot_df, condition_label) {
  sub <- plot_df[plot_df$condition == condition_label, ]
  box <- cone_box(sub)
  vmax <- max(quantile(plot_df$Gls, 0.995, na.rm = TRUE), 1e-6)

  ggplot(sub, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "#D0D0D0", size = 0.25, alpha = 0.75, stroke = 0) +
    geom_point(
      data = sub[sub$Gls > 0, ],
      aes(color = Gls),
      size = 0.32,
      alpha = 0.9,
      stroke = 0
    ) +
    geom_rect(
      data = box,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = NA, color = "#0072B2", linewidth = 0.8
    ) +
    annotate("text", x = median(sub$UMAP_1, na.rm = TRUE), y = median(sub$UMAP_2, na.rm = TRUE),
             label = condition_label, size = 4.2, fontface = "bold") +
    scale_color_gradient(low = "#FEE5D9", high = "#A50F15", limits = c(0, vmax), guide = "none") +
    coord_fixed() +
    labs(x = "UMAP_1", y = "UMAP_2") +
    base_umap_theme()
}

make_panel <- function(age_label) {
  age_df <- df[df$age == age_label, ]
  panel_a <- plot_celltypes(age_df) +
    labs(tag = "A") +
    theme(plot.tag = element_text(size = 24, face = "bold"))
  panel_b <- (
    (plot_gls(age_df, "WT") +
       labs(tag = "B") +
       theme(plot.tag = element_text(size = 24, face = "bold"))) |
      (plot_gls(age_df, "cKO") +
         ggtitle("Gls Gene plot") +
         theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  )

  final <- panel_a + panel_b +
    plot_layout(widths = c(1.12, 1.55))

  pdf_file <- file.path(out_dir, paste0("sherine_seurat_style_umap_gls_", age_label, ".pdf"))
  png_file <- file.path(out_dir, paste0("sherine_seurat_style_umap_gls_", age_label, ".png"))
  ggsave(pdf_file, final, width = 12.8, height = 5.0, units = "in", device = cairo_pdf)
  ggsave(png_file, final, width = 12.8, height = 5.0, units = "in", dpi = 300)
}

make_panel("P15")
make_panel("P35")

message("Done. Outputs written to: ", out_dir)
