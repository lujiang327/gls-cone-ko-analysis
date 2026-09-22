#!/usr/bin/env python3

from pathlib import Path
import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / ".matplotlib-cache"),
)

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCANPY_DIR = PROJECT_DIR / "scRNA_4_thomas" / "scanpy"
COUNTS_H5AD = SCANPY_DIR / "photoreceptors.h5ad"
ANNOTATED_H5AD = SCANPY_DIR / "annotated_photoreceptors.h5ad"
OUT_DIR = PROJECT_DIR / "results" / "sherine_style_umap_gls_panels"

SAMPLE_MAP = {
    "15dayS1": {"age": "P15", "condition": "WT"},
    "15dayS2": {"age": "P15", "condition": "cKO"},
    "35dayS1": {"age": "P35", "condition": "WT"},
    "35dayS2": {"age": "P35", "condition": "cKO"},
}

CELLTYPE_ORDER = [
    "Rod",
    "BC",
    "Cones",
    "AC",
    "MG",
    "HC",
    "Endothelial",
    "Microglia",
    "RGC",
    "MG_AC",
]

CELLTYPE_LABELS = {
    "Rod": "Rod",
    "BC": "Bipolar",
    "Cones": "Cone",
    "AC": "Amacrine",
    "MG": "Muller glia",
    "HC": "Horizontal",
    "Endothelial": "V/E cells",
    "Microglia": "Microglia",
    "RGC": "RGC",
    "MG_AC": "MG/AC",
}

COLORS = {
    "Rod": "#F8766D",
    "BC": "#E69F00",
    "Cones": "#7CAE00",
    "AC": "#00BA38",
    "MG": "#00BFC4",
    "HC": "#C77CFF",
    "Endothelial": "#619CFF",
    "Microglia": "#F564E3",
    "RGC": "#B79F00",
    "MG_AC": "#FF61C3",
}


def dense_vector(x):
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def read_data():
    counts = sc.read_h5ad(COUNTS_H5AD)
    ann = sc.read_h5ad(ANNOTATED_H5AD, backed="r")
    if counts.n_obs != ann.n_obs:
        raise ValueError("Raw-count and annotated h5ad objects have different n_obs.")
    if not np.array_equal(counts.obs["sample"].astype(str).values, ann.obs["sample"].astype(str).values):
        raise ValueError("Raw-count and annotated h5ad objects are not aligned by sample order.")

    meta = ann.obs[["sample", "celltype"]].copy()
    meta["sample"] = meta["sample"].astype(str)
    meta["celltype"] = meta["celltype"].astype(str)
    meta["age"] = meta["sample"].map(lambda s: SAMPLE_MAP[s]["age"])
    meta["condition"] = meta["sample"].map(lambda s: SAMPLE_MAP[s]["condition"])
    meta["umap_1"] = ann.obsm["X_umap"][:, 0]
    meta["umap_2"] = ann.obsm["X_umap"][:, 1]

    gene_idx = counts.var_names.get_loc("Gls")
    raw = dense_vector(counts[:, gene_idx].X).astype(float)
    library_size = counts.obs["total_counts"].to_numpy(dtype=float)
    meta["Gls_lognorm"] = np.log1p(raw / library_size * 1e4)
    meta["Gls_detected"] = raw > 0
    ann.file.close()
    return meta


def cone_box(df, pad=0.35):
    cones = df[df["celltype"] == "Cones"]
    if cones.empty:
        return None
    x1, x2 = np.percentile(cones["umap_1"], [1, 99])
    y1, y2 = np.percentile(cones["umap_2"], [1, 99])
    return (x1 - pad, y1 - pad, (x2 - x1) + 2 * pad, (y2 - y1) + 2 * pad)


def add_cone_box(ax, df):
    box = cone_box(df)
    if box is None:
        return
    ax.add_patch(Rectangle((box[0], box[1]), box[2], box[3], fill=False, lw=1.8, ec="#0072B2"))


def style_umap_axis(ax):
    ax.set_xlabel("UMAP_1", fontsize=13)
    ax.set_ylabel("UMAP_2", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=10)
    ax.set_aspect("equal", adjustable="box")


def plot_celltype_panel(ax, df):
    for ct in CELLTYPE_ORDER:
        sub = df[df["celltype"] == ct]
        if sub.empty:
            continue
        ax.scatter(sub["umap_1"], sub["umap_2"], s=2.0, c=COLORS[ct], linewidth=0, alpha=0.9)

    for ct in CELLTYPE_ORDER:
        sub = df[df["celltype"] == ct]
        if sub.empty:
            continue
        x, y = sub[["umap_1", "umap_2"]].median()
        ax.text(x, y, CELLTYPE_LABELS.get(ct, ct), fontsize=7.5, ha="center", va="center")

    add_cone_box(ax, df)
    ax.text(0.30, 0.58, "Cone", transform=ax.transAxes, fontsize=10, fontweight="bold")
    style_umap_axis(ax)

    handles = [
        Line2D([0], [0], marker="o", linestyle="", label=CELLTYPE_LABELS.get(ct, ct),
               markerfacecolor=COLORS[ct], markeredgecolor=COLORS[ct], markersize=7)
        for ct in CELLTYPE_ORDER
        if ct in set(df["celltype"])
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=9)


def plot_gls_panel(ax, df, condition, title=None):
    sub = df[df["condition"] == condition]
    ax.scatter(sub["umap_1"], sub["umap_2"], s=2.0, c="#D0D0D0", linewidth=0, alpha=0.75)
    detected = sub[sub["Gls_detected"]]
    vmax = max(df["Gls_lognorm"].quantile(0.995), 1e-6)
    ax.scatter(
        detected["umap_1"],
        detected["umap_2"],
        s=2.5,
        c=detected["Gls_lognorm"],
        cmap="Reds",
        vmin=0,
        vmax=vmax,
        linewidth=0,
        alpha=0.92,
    )
    add_cone_box(ax, sub)
    ax.text(0.47, 0.52, condition, transform=ax.transAxes, fontsize=11, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=13, fontweight="bold")
    style_umap_axis(ax)


def composite_for_age(df, age):
    age_df = df[df["age"] == age].copy()
    fig = plt.figure(figsize=(12.5, 5.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1], wspace=0.42)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b1 = fig.add_subplot(gs[0, 1])
    ax_b2 = fig.add_subplot(gs[0, 2])

    plot_celltype_panel(ax_a, age_df)
    plot_gls_panel(ax_b1, age_df, "WT")
    plot_gls_panel(ax_b2, age_df, "cKO", title="Gls Gene plot")

    ax_a.text(-0.28, 1.05, "A", transform=ax_a.transAxes, fontsize=24, fontweight="bold")
    ax_b1.text(-0.42, 1.05, "B", transform=ax_b1.transAxes, fontsize=24, fontweight="bold")
    fig.suptitle(f"Sherine-style UMAP / Gls feature panels ({age})", fontsize=15, y=1.02)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / f"sherine_style_umap_gls_{age}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"sherine_style_umap_gls_{age}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    df = read_data()
    for age in ["P15", "P35"]:
        composite_for_age(df, age)
    df.groupby(["age", "condition", "celltype"], observed=True).size().reset_index(name="n_cells").to_csv(
        OUT_DIR / "cells_used_for_sherine_style_panels.csv", index=False
    )
    print(f"Done. Outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
