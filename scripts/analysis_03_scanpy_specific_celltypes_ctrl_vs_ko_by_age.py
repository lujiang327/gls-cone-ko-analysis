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
from scipy import sparse, stats
import matplotlib.pyplot as plt
import seaborn as sns


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCANPY_DIR = PROJECT_DIR / "scRNA_4_thomas" / "scanpy"
COUNTS_H5AD = SCANPY_DIR / "photoreceptors.h5ad"
ANNOTATED_H5AD = SCANPY_DIR / "annotated_photoreceptors.h5ad"
OUT_DIR = PROJECT_DIR / "results" / "analysis_03_scanpy_specific_celltypes_ctrl_vs_ko_by_age"

GENES = ["Gls", "Gls2"]
SAMPLE_MAP = {
    "15dayS1": {"age": "P15", "condition": "Ctrl"},
    "15dayS2": {"age": "P15", "condition": "KO"},
    "35dayS1": {"age": "P35", "condition": "Ctrl"},
    "35dayS2": {"age": "P35", "condition": "KO"},
}
AGE_CONDITION_ORDER = ["P15_Ctrl", "P15_KO", "P35_Ctrl", "P35_KO"]
CELLTYPE_ORDER = [
    "Cones",
    "Rod",
    "BC",
    "AC",
    "RGC",
    "HC",
    "MG",
    "MG_AC",
    "Endothelial",
    "Microglia",
]


def dense_vector(matrix):
    if sparse.issparse(matrix):
        return np.asarray(matrix.toarray()).ravel()
    return np.asarray(matrix).ravel()


def savefig(fig, stem, width=None, height=None):
    if width and height:
        fig.set_size_inches(width, height)
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def bh_adjust(p_values):
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    out[valid] = restored
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    print(f"Loading raw-count object: {COUNTS_H5AD}")
    counts_ad = sc.read_h5ad(COUNTS_H5AD)
    print(f"Loading annotated object: {ANNOTATED_H5AD}")
    ann_ad = sc.read_h5ad(ANNOTATED_H5AD)

    if counts_ad.n_obs != ann_ad.n_obs:
        raise ValueError("Raw-count and annotated objects have different cell counts.")
    if not np.array_equal(counts_ad.obs["sample"].astype(str).values, ann_ad.obs["sample"].astype(str).values):
        raise ValueError("Raw-count and annotated objects are not aligned by sample order.")
    missing_genes = [gene for gene in GENES if gene not in counts_ad.var_names]
    if missing_genes:
        raise ValueError(f"Missing genes in raw-count object: {missing_genes}")
    if "celltype" not in ann_ad.obs:
        raise ValueError("Annotated object does not contain obs['celltype'].")
    if "X_umap" not in ann_ad.obsm:
        raise ValueError("Annotated object does not contain obsm['X_umap'].")

    meta = ann_ad.obs[["sample", "celltype"]].copy()
    meta["sample"] = meta["sample"].astype(str)
    meta["celltype"] = meta["celltype"].astype(str)
    meta["age"] = meta["sample"].map(lambda x: SAMPLE_MAP[x]["age"])
    meta["condition"] = meta["sample"].map(lambda x: SAMPLE_MAP[x]["condition"])
    meta["age_condition"] = meta["age"] + "_" + meta["condition"]
    meta["cell"] = [f"cell_{i + 1}" for i in range(meta.shape[0])]
    meta["umap_1"] = ann_ad.obsm["X_umap"][:, 0]
    meta["umap_2"] = ann_ad.obsm["X_umap"][:, 1]
    meta["age"] = pd.Categorical(meta["age"], categories=["P15", "P35"], ordered=True)
    meta["condition"] = pd.Categorical(meta["condition"], categories=["Ctrl", "KO"], ordered=True)
    meta["age_condition"] = pd.Categorical(meta["age_condition"], categories=AGE_CONDITION_ORDER, ordered=True)
    present_celltypes = [ct for ct in CELLTYPE_ORDER if ct in set(meta["celltype"])]
    meta["celltype"] = pd.Categorical(meta["celltype"], categories=present_celltypes, ordered=True)

    cell_counts = (
        meta.groupby(["age", "condition", "sample", "celltype"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    totals = cell_counts.groupby(["age", "condition", "sample"], observed=True)["n_cells"].transform("sum")
    cell_counts["total_cells"] = totals
    cell_counts["pct_of_sample"] = cell_counts["n_cells"] / cell_counts["total_cells"] * 100
    cell_counts.to_csv(OUT_DIR / "scanpy_cell_counts_and_percent_by_age_condition_celltype.csv", index=False)

    expr_frames = []
    library_size = counts_ad.obs["total_counts"].to_numpy(dtype=float)
    for gene in GENES:
        gene_idx = counts_ad.var_names.get_loc(gene)
        raw = dense_vector(counts_ad[:, gene_idx].X).astype(float)
        lognorm = np.log1p(raw / library_size * 1e4)
        frame = meta[["cell", "sample", "age", "condition", "age_condition", "celltype", "umap_1", "umap_2"]].copy()
        frame["gene"] = gene
        frame["raw_count"] = raw
        frame["detected"] = raw > 0
        frame["lognorm_expr"] = lognorm
        expr_frames.append(frame)

    expr_long = pd.concat(expr_frames, ignore_index=True)
    expr_long.to_csv(OUT_DIR / "scanpy_gls_gls2_single_cell_expression.csv", index=False)

    summary = (
        expr_long.groupby(["gene", "age", "condition", "sample", "age_condition", "celltype"], observed=True)
        .agg(
            n_cells=("cell", "size"),
            pct_detected=("detected", lambda x: float(np.mean(x) * 100)),
            avg_lognorm_expr=("lognorm_expr", "mean"),
            median_lognorm_expr=("lognorm_expr", "median"),
            avg_raw_count=("raw_count", "mean"),
        )
        .reset_index()
        .sort_values(["gene", "age", "celltype", "condition"])
    )
    summary.to_csv(OUT_DIR / "scanpy_gls_gls2_summary_by_age_condition_celltype.csv", index=False)

    stat_rows = []
    for (gene, age, celltype), sub in expr_long.groupby(["gene", "age", "celltype"], observed=True):
        ctrl = sub.loc[sub["condition"] == "Ctrl", "lognorm_expr"].to_numpy()
        ko = sub.loc[sub["condition"] == "KO", "lognorm_expr"].to_numpy()
        if len(ctrl) == 0 or len(ko) == 0:
            continue
        if np.unique(np.concatenate([ctrl, ko])).size == 1:
            p_value = np.nan
        else:
            p_value = stats.mannwhitneyu(ctrl, ko, alternative="two-sided").pvalue
        stat_rows.append(
            {
                "gene": gene,
                "age": age,
                "celltype": celltype,
                "ctrl_n": len(ctrl),
                "ko_n": len(ko),
                "ctrl_avg_lognorm_expr": float(np.mean(ctrl)),
                "ko_avg_lognorm_expr": float(np.mean(ko)),
                "difference_ko_minus_ctrl": float(np.mean(ko) - np.mean(ctrl)),
                "ctrl_pct_detected": float(np.mean(ctrl > 0) * 100),
                "ko_pct_detected": float(np.mean(ko > 0) * 100),
                "p_value": p_value,
            }
        )
    stats_df = pd.DataFrame(stat_rows)
    if not stats_df.empty:
        stats_df["p_adj_bh_within_age_gene"] = np.nan
        for (gene, age), idx in stats_df.groupby(["gene", "age"], observed=True).groups.items():
            stats_df.loc[idx, "p_adj_bh_within_age_gene"] = bh_adjust(stats_df.loc[idx, "p_value"])
        stats_df = stats_df.sort_values(["gene", "age", "celltype"])
    stats_df.to_csv(OUT_DIR / "scanpy_gls_gls2_ctrl_vs_ko_by_age_celltype_wilcox_exploratory.csv", index=False)

    composition_plot(cell_counts)
    dot_plot(summary)
    gls_heatmap(summary)
    umap_overview(meta)
    feature_umaps(expr_long)

    run_info = pd.DataFrame(
        [
            {"item": "raw_count_input", "value": str(COUNTS_H5AD.relative_to(PROJECT_DIR))},
            {"item": "annotation_umap_input", "value": str(ANNOTATED_H5AD.relative_to(PROJECT_DIR))},
            {"item": "normalization", "value": "log1p(raw_count / cell_total_counts * 10000)"},
            {"item": "sample_mapping", "value": "15dayS1=P15 Ctrl; 15dayS2=P15 KO; 35dayS1=P35 Ctrl; 35dayS2=P35 KO"},
        ]
    )
    run_info.to_csv(OUT_DIR / "scanpy_analysis_03_run_info.csv", index=False)
    print(f"Done. Outputs written to: {OUT_DIR}")


def composition_plot(cell_counts):
    fig, axes = plt.subplots(1, 2, figsize=(9, 5), sharey=True)
    palette = sns.color_palette("tab10", n_colors=len(CELLTYPE_ORDER))
    for ax, age in zip(axes, ["P15", "P35"]):
        sub = cell_counts[cell_counts["age"].astype(str) == age]
        pivot = sub.pivot_table(
            index="condition",
            columns="celltype",
            values="pct_of_sample",
            aggfunc="sum",
            observed=True,
        ).reindex(["Ctrl", "KO"])
        pivot = pivot[[ct for ct in CELLTYPE_ORDER if ct in pivot.columns]]
        pivot.plot(kind="bar", stacked=True, ax=ax, color=palette[: pivot.shape[1]], width=0.75)
        ax.set_title(age)
        ax.set_xlabel("")
        ax.set_ylabel("% of cells")
        ax.tick_params(axis="x", rotation=0)
        ax.legend_.remove()
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, title="Cell type", bbox_to_anchor=(1.02, 0.5), loc="center left")
    savefig(fig, "scanpy_celltype_composition_by_age_condition", 11, 5)


def dot_plot(summary):
    fig, axes = plt.subplots(1, len(GENES), figsize=(10, 5.5), sharey=True)
    if len(GENES) == 1:
        axes = [axes]
    vmax = max(summary["avg_lognorm_expr"].max(), 1e-6)
    for ax, gene in zip(axes, GENES):
        sub = summary[summary["gene"] == gene].copy()
        x_map = {label: i for i, label in enumerate(AGE_CONDITION_ORDER)}
        y_order = [ct for ct in CELLTYPE_ORDER if ct in set(sub["celltype"].astype(str))]
        y_map = {label: i for i, label in enumerate(y_order)}
        sca = ax.scatter(
            sub["age_condition"].astype(str).map(x_map),
            sub["celltype"].astype(str).map(y_map),
            s=np.clip(sub["pct_detected"], 0, 100) * 3 + 8,
            c=sub["avg_lognorm_expr"],
            cmap="viridis",
            vmin=0,
            vmax=vmax,
            edgecolor="black",
            linewidth=0.3,
        )
        ax.set_title(gene)
        ax.set_xticks(range(len(AGE_CONDITION_ORDER)))
        ax.set_xticklabels(AGE_CONDITION_ORDER, rotation=45, ha="right")
        ax.set_yticks(range(len(y_order)))
        ax.set_yticklabels(y_order)
        ax.set_xlabel("")
        ax.grid(False)
    cbar = fig.colorbar(sca, ax=axes, fraction=0.025, pad=0.03)
    cbar.set_label("Avg log-normalized expression")
    for pct, size in [(25, 25 * 3 + 8), (50, 50 * 3 + 8), (75, 75 * 3 + 8)]:
        axes[-1].scatter([], [], s=size, c="lightgray", edgecolor="black", linewidth=0.3, label=f"{pct}%")
    axes[-1].legend(title="% detected", bbox_to_anchor=(1.32, 0.3), loc="center left", frameon=True)
    savefig(fig, "scanpy_gls_gls2_dotplot_ctrl_vs_ko_by_age_celltype", 12, 5.5)


def gls_heatmap(summary):
    gls = summary[summary["gene"] == "Gls"].copy()
    gls["group"] = gls["age_condition"].astype(str)
    avg = gls.pivot_table(index="celltype", columns="group", values="avg_lognorm_expr", observed=True)
    pct = gls.pivot_table(index="celltype", columns="group", values="pct_detected", observed=True)
    y_order = [ct for ct in CELLTYPE_ORDER if ct in avg.index]
    avg = avg.reindex(index=y_order, columns=AGE_CONDITION_ORDER)
    pct = pct.reindex(index=y_order, columns=AGE_CONDITION_ORDER)
    labels = pct.map(lambda x: "" if pd.isna(x) else f"{x:.1f}%")
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    sns.heatmap(
        avg,
        annot=labels,
        fmt="",
        cmap="viridis",
        cbar_kws={"label": "Avg log-normalized expression"},
        linewidths=0.5,
        linecolor="white",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Scanpy annotated object: Gls expression")
    savefig(fig, "scanpy_gls_heatmap_ctrl_vs_ko_by_age_celltype", 7.5, 5.5)


def umap_overview(meta):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    sns.scatterplot(
        data=meta,
        x="umap_1",
        y="umap_2",
        hue="celltype",
        s=4,
        linewidth=0,
        ax=axes[0],
        palette="tab10",
    )
    axes[0].set_title("Cell type")
    axes[0].legend(markerscale=4, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    sns.scatterplot(
        data=meta,
        x="umap_1",
        y="umap_2",
        hue="age_condition",
        s=4,
        linewidth=0,
        ax=axes[1],
        palette="Set2",
    )
    axes[1].set_title("Age and condition")
    axes[1].legend(markerscale=4, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True)
    for ax in axes:
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_aspect("equal", adjustable="box")
    savefig(fig, "scanpy_umap_celltype_and_age_condition", 13, 5)


def feature_umaps(expr_long):
    for gene in GENES:
        sub_gene = expr_long[expr_long["gene"] == gene]
        vmax = np.percentile(sub_gene["lognorm_expr"], 99.5)
        fig, axes = plt.subplots(1, 4, figsize=(14, 3.5), sharex=True, sharey=True)
        for ax, group in zip(axes, AGE_CONDITION_ORDER):
            sub = sub_gene[sub_gene["age_condition"].astype(str) == group]
            ax.scatter(sub["umap_1"], sub["umap_2"], s=3, c="lightgray", linewidth=0, alpha=0.35)
            sca = ax.scatter(
                sub["umap_1"],
                sub["umap_2"],
                s=4,
                c=sub["lognorm_expr"],
                cmap="magma",
                vmin=0,
                vmax=vmax,
                linewidth=0,
            )
            ax.set_title(group)
            ax.set_xlabel("UMAP 1")
            ax.set_ylabel("UMAP 2")
            ax.set_aspect("equal", adjustable="box")
        cbar = fig.colorbar(sca, ax=axes, fraction=0.02, pad=0.02)
        cbar.set_label(f"{gene} log-normalized expression")
        savefig(fig, f"scanpy_{gene.lower()}_feature_umap_by_age_condition", 14, 3.5)


if __name__ == "__main__":
    main()
