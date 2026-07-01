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
    sns.set_theme(style="whitegrid", context="talk")

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
    cone_abundance_plot(cell_counts)
    dot_plot(summary)
    for gene in GENES:
        gene_dot_plot(summary, gene)
        gene_heatmap(summary, gene)
        effect_plot(stats_df, gene)
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


def cone_abundance_plot(cell_counts):
    cone_counts = cell_counts[cell_counts["celltype"].astype(str) == "Cones"].copy()
    if cone_counts.empty:
        return
    cone_counts["age"] = pd.Categorical(cone_counts["age"].astype(str), categories=["P15", "P35"], ordered=True)
    cone_counts["condition"] = pd.Categorical(cone_counts["condition"].astype(str), categories=["Ctrl", "KO"], ordered=True)
    cone_counts["label"] = cone_counts.apply(
        lambda row: f"n={int(row['n_cells'])}\n{row['pct_of_sample']:.1f}%",
        axis=1,
    )
    fig, axes = plt.subplots(1, 2, figsize=(9, 5.5), sharey=True)
    colors = {"Ctrl": "#2C7BB6", "KO": "#D7191C"}
    for ax, age in zip(axes, ["P15", "P35"]):
        sub = cone_counts[cone_counts["age"].astype(str) == age]
        ax.bar(
            sub["condition"].astype(str),
            sub["pct_of_sample"],
            color=[colors[c] for c in sub["condition"].astype(str)],
            width=0.65,
            alpha=0.9,
        )
        for _, row in sub.iterrows():
            ax.text(
                str(row["condition"]),
                row["pct_of_sample"] + 0.12,
                row["label"],
                ha="center",
                va="bottom",
                fontsize=15,
            )
        ax.set_title(age)
        ax.set_xlabel("")
        ax.set_ylim(0, max(cone_counts["pct_of_sample"].max() * 1.35, 1))
        ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("Cones (% of total annotated cells)")
    fig.suptitle("Cone abundance by age and condition", fontweight="bold")
    savefig(fig, "moloy_review_scanpy_cone_abundance_percent_and_counts_by_age_condition", 9, 5.5)


def dot_plot(summary):
    fig, axes = plt.subplots(1, len(GENES), figsize=(11, 6.5), sharey=True)
    if len(GENES) == 1:
        axes = [axes]
    vmax = max(summary["avg_lognorm_expr"].max(), 1e-6)
    for ax, gene in zip(axes, GENES):
        sub = summary[summary["gene"] == gene].copy()
        x_map = {label: i for i, label in enumerate(AGE_CONDITION_ORDER)}
        y_order = list(reversed([ct for ct in CELLTYPE_ORDER if ct in set(sub["celltype"].astype(str))]))
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
        cone_sub = sub[sub["celltype"].astype(str) == "Cones"]
        ax.scatter(
            cone_sub["age_condition"].astype(str).map(x_map),
            cone_sub["celltype"].astype(str).map(y_map),
            s=np.clip(cone_sub["pct_detected"], 0, 100) * 3 + 55,
            facecolors="none",
            edgecolors="#FFD700",
            linewidth=2.6,
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
    savefig(fig, "scanpy_gls_gls2_dotplot_ctrl_vs_ko_by_age_celltype", 13, 6.5)


def format_p_value(value):
    if pd.isna(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def gene_dot_plot(summary, gene):
    sub = summary[summary["gene"] == gene].copy()
    fig, ax = plt.subplots(figsize=(8.5, 7))
    x_map = {label: i for i, label in enumerate(AGE_CONDITION_ORDER)}
    y_order = list(reversed([ct for ct in CELLTYPE_ORDER if ct in set(sub["celltype"].astype(str))]))
    y_map = {label: i for i, label in enumerate(y_order)}
    vmax = max(sub["avg_lognorm_expr"].max(), 1e-6)
    sca = ax.scatter(
        sub["age_condition"].astype(str).map(x_map),
        sub["celltype"].astype(str).map(y_map),
        s=np.clip(sub["pct_detected"], 0, 100) * 4 + 18,
        c=sub["avg_lognorm_expr"],
        cmap="viridis",
        vmin=0,
        vmax=vmax,
        edgecolor="black",
        linewidth=0.4,
    )
    cone_sub = sub[sub["celltype"].astype(str) == "Cones"]
    ax.scatter(
        cone_sub["age_condition"].astype(str).map(x_map),
        cone_sub["celltype"].astype(str).map(y_map),
        s=np.clip(cone_sub["pct_detected"], 0, 100) * 4 + 80,
        facecolors="none",
        edgecolors="#FFD700",
        linewidth=3,
    )
    ax.set_xticks(range(len(AGE_CONDITION_ORDER)))
    ax.set_xticklabels(AGE_CONDITION_ORDER, rotation=45, ha="right")
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels(y_order)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(f"{gene} expression across retinal cell types")
    ax.grid(False)
    cbar = fig.colorbar(sca, ax=ax, fraction=0.045, pad=0.04)
    cbar.set_label("Avg log-normalized expression")
    for pct, size in [(25, 25 * 4 + 18), (50, 50 * 4 + 18), (75, 75 * 4 + 18)]:
        ax.scatter([], [], s=size, c="lightgray", edgecolor="black", linewidth=0.4, label=f"{pct}%")
    ax.legend(title="% detected", bbox_to_anchor=(1.28, 0.35), loc="center left", frameon=True)
    savefig(fig, f"moloy_review_scanpy_{gene.lower()}_dotplot_ctrl_vs_ko_by_age_celltype", 10.5, 7)


def gene_heatmap(summary, gene):
    gene_summary = summary[summary["gene"] == gene].copy()
    gene_summary["group"] = gene_summary["age_condition"].astype(str)
    avg = gene_summary.pivot_table(index="celltype", columns="group", values="avg_lognorm_expr", observed=True)
    pct = gene_summary.pivot_table(index="celltype", columns="group", values="pct_detected", observed=True)
    y_order = [ct for ct in CELLTYPE_ORDER if ct in avg.index]
    avg = avg.reindex(index=y_order, columns=AGE_CONDITION_ORDER)
    pct = pct.reindex(index=y_order, columns=AGE_CONDITION_ORDER)
    labels = pct.map(lambda x: "" if pd.isna(x) else f"{x:.0f}%")
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sns.heatmap(
        avg,
        annot=labels,
        fmt="",
        cmap="viridis",
        cbar_kws={"label": "Avg log-normalized expression"},
        linewidths=0.7,
        linecolor="white",
        annot_kws={"fontsize": 15, "fontweight": "bold"},
        ax=ax,
    )
    ax.add_patch(plt.Rectangle((0, 0), len(AGE_CONDITION_ORDER), 1, fill=False, edgecolor="#FFD700", linewidth=4))
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(f"{gene} expression; labels show % detected")
    savefig(fig, f"moloy_review_scanpy_{gene.lower()}_heatmap_ctrl_vs_ko_by_age_celltype", 8, 6.5)


def effect_plot(stats_df, gene):
    sub = stats_df[stats_df["gene"] == gene].copy()
    if sub.empty:
        return
    sub["age"] = pd.Categorical(sub["age"].astype(str), categories=["P15", "P35"], ordered=True)
    sub["celltype"] = pd.Categorical(sub["celltype"].astype(str), categories=CELLTYPE_ORDER, ordered=True)
    sub["label"] = sub.apply(
        lambda row: f"p={format_p_value(row['p_value'])}\nFDR={format_p_value(row['p_adj_bh_within_age_gene'])}",
        axis=1,
    )
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
    y_order = list(reversed([ct for ct in CELLTYPE_ORDER if ct in set(sub["celltype"].astype(str))]))
    y_map = {label: i for i, label in enumerate(y_order)}
    max_abs = max(abs(sub["difference_ko_minus_ctrl"]).max(), 0.01)
    for ax, age in zip(axes, ["P15", "P35"]):
        age_sub = sub[sub["age"].astype(str) == age].copy()
        colors = ["#F28E2B" if ct == "Cones" else "#4E79A7" for ct in age_sub["celltype"].astype(str)]
        ax.barh(
            age_sub["celltype"].astype(str).map(y_map),
            age_sub["difference_ko_minus_ctrl"],
            color=colors,
            edgecolor=["#FFD700" if ct == "Cones" else "white" for ct in age_sub["celltype"].astype(str)],
            linewidth=[3 if ct == "Cones" else 0.5 for ct in age_sub["celltype"].astype(str)],
            height=0.7,
        )
        for _, row in age_sub.iterrows():
            x = row["difference_ko_minus_ctrl"]
            ha = "right" if x < 0 else "left"
            offset = -0.03 * max_abs if x < 0 else 0.03 * max_abs
            ax.text(x + offset, y_map[str(row["celltype"])], row["label"], va="center", ha=ha, fontsize=10)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlim(-max_abs * 1.65, max_abs * 1.65)
        ax.set_yticks(range(len(y_order)))
        ax.set_yticklabels(y_order)
        ax.set_title(age)
        ax.set_xlabel("KO - Ctrl mean log-normalized expression")
        ax.grid(axis="y", visible=False)
    axes[0].set_ylabel("")
    fig.suptitle(f"{gene} Ctrl-vs-KO effect by cell type; labels show p value and BH FDR", fontweight="bold")
    savefig(fig, f"moloy_review_scanpy_{gene.lower()}_ko_minus_ctrl_effect_with_p_fdr", 14, 7)


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
