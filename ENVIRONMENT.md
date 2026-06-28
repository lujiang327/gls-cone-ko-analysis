# R Environment

This project continues Sherine's `scRNA_4_thomas` analysis. The existing scripts are
R/Seurat scripts that load Seurat, Signac, Bioconductor mouse genome annotation
packages, and plotting packages.

## Detected Requirements

The active scripts require:

- R with Seurat/SeuratObject and Signac
- Bioconductor annotation packages: `GenomeInfoDb`, `EnsDb.Mmusculus.v79`,
  `BSgenome.Mmusculus.UCSC.mm10`
- Plotting/data packages: `ggplot2`, `patchwork`, `Matrix`, `dplyr`, `tidyr`,
  `tibble`, `scales`, `RColorBrewer`, `viridis`, `gplots`, `heatmaply`
- `hdf5r` for reading 10x `.h5` matrices through Seurat
- optional `snakemake-minimal` for the existing `Snakefile`

The original Great Lakes scripts activated a conda environment named `archr`,
but the current workflow is Seurat/Signac-based.

## Conda/Mamba Setup

This repository now includes a self-contained micromamba binary under
`tools/micromamba/` and a local environment under `.micromamba/`.

Create the environment:

```bash
mamba env create -f envs/r-seurat.yml
```

or, if mamba is unavailable:

```bash
conda env create -f envs/r-seurat.yml
```

Activate it:

```bash
conda activate gls-cone-ko-r
```

Validate it:

```bash
Rscript scripts/check_r_env.R
```

With the self-contained project install, validate using:

```bash
.micromamba/envs/gls-cone-ko-r/bin/Rscript scripts/check_r_env.R
```

To run Sherine's scripts from the project environment, use the same `Rscript`
path, for example:

```bash
cd scRNA_4_thomas
../.micromamba/envs/gls-cone-ko-r/bin/Rscript plot_figures.R photoreceptors
```

## Notes

This local machine initially did not have `R`, `Rscript`, `conda`, `mamba`, or
`brew` available in `PATH`, so a self-contained micromamba install was used.
`micromamba run` may try to write locks under `~/.cache/mamba`; direct calls to
`.micromamba/envs/gls-cone-ko-r/bin/Rscript` avoid that sandbox/home-cache issue.

## Scanpy Environment

Sherine's Great Lakes Scanpy scripts activated a conda environment named
`archr`, and the Slurm log points to Python 3.7 under that environment. For this
Mac/local continuation, the project uses a compatible modern Scanpy environment:

```bash
tools/micromamba/bin/micromamba create -y -r .micromamba -f envs/scanpy.yml
```

Run the Scanpy-specific Analysis 03 with:

```bash
MPLCONFIGDIR=.matplotlib-cache \
  .micromamba/envs/gls-cone-ko-scanpy/bin/python \
  scripts/analysis_03_scanpy_specific_celltypes_ctrl_vs_ko_by_age.py
```

This script uses `scRNA_4_thomas/scanpy/photoreceptors.h5ad` for raw counts and
`scRNA_4_thomas/scanpy/annotated_photoreceptors.h5ad` for the detailed cell type
labels and UMAP coordinates. The sample labels follow Moloy's confirmation:
`15dayS1=P15 Ctrl`, `15dayS2=P15 KO`, `35dayS1=P35 Ctrl`, and `35dayS2=P35 KO`.

Validated Scanpy package versions:

| Package | Version |
| --- | --- |
| Python | 3.10.20 |
| scanpy | 1.11.5 |
| anndata | 0.11.4 |
| pandas | 2.3.3 |
| numpy | 2.2.6 |
| scipy | 1.15.2 |
| h5py | 3.16.0 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |

Validated package versions:

| Package | Version |
| --- | --- |
| R | 4.4.3 |
| Seurat | 5.5.0 |
| SeuratObject | 5.4.0 |
| Signac | 1.16.0 |
| GenomeInfoDb | 1.42.0 |
| EnsDb.Mmusculus.v79 | 2.99.0 |
| BSgenome.Mmusculus.UCSC.mm10 | 1.4.3 |
| ggplot2 | 4.0.3 |
| dplyr | 1.2.1 |
| Matrix | 1.7.5 |
