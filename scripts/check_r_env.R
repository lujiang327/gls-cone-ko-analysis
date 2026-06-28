required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Signac",
  "GenomeInfoDb",
  "EnsDb.Mmusculus.v79",
  "BSgenome.Mmusculus.UCSC.mm10",
  "future",
  "ggplot2",
  "patchwork",
  "Matrix",
  "dplyr",
  "tidyr",
  "tibble",
  "scales",
  "RColorBrewer",
  "viridis",
  "gplots",
  "heatmaply",
  "hdf5r"
)

cat("R version:", R.version.string, "\n\n")

status <- data.frame(
  package = required_packages,
  installed = vapply(required_packages, requireNamespace, logical(1), quietly = TRUE),
  version = NA_character_,
  stringsAsFactors = FALSE
)

for (pkg in status$package[status$installed]) {
  status$version[status$package == pkg] <- as.character(utils::packageVersion(pkg))
}

print(status, row.names = FALSE)

missing <- status$package[!status$installed]
if (length(missing) > 0) {
  stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

cat("\nAll required packages are installed.\n")
