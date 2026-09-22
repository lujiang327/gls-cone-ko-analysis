#!/usr/bin/env python3

from pathlib import Path
import os
import warnings

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
import seaborn as sns

try:
    import gseapy as gp
except ImportError:
    gp = None


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCANPY_DIR = PROJECT_DIR / "scRNA_4_thomas" / "scanpy"
COUNTS_H5AD = SCANPY_DIR / "photoreceptors.h5ad"
ANNOTATED_H5AD = SCANPY_DIR / "annotated_photoreceptors.h5ad"
OUT_DIR = PROJECT_DIR / "results" / "analysis_04_degs_and_pathways"
TABLE_DIR = OUT_DIR / "tables"
FIG_DIR = OUT_DIR / "figures"
PATHWAY_DIR = OUT_DIR / "pathway_enrichment"

SAMPLE_MAP = {
    "15dayS1": {"age": "P15", "condition": "Ctrl"},
    "15dayS2": {"age": "P15", "condition": "KO"},
    "35dayS1": {"age": "P35", "condition": "Ctrl"},
    "35dayS2": {"age": "P35", "condition": "KO"},
}
AGE_ORDER = ["P15", "P35"]
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

MIN_CELLS_PER_GROUP = 10
MIN_PCT_EXPRESSED = 0.10
FDR_CUTOFF = 0.05
LOGFC_CUTOFF = 0.25
PATHWAY_TOP_N = 15


def dense_vector(x):
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def savefig(fig, filename, width=None, height=None):
    if width and height:
        fig.set_size_inches(width, height)
    fig.savefig(FIG_DIR / f"{filename}.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{filename}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def read_and_prepare():
    print(f"Loading raw-count object: {COUNTS_H5AD}")
    adata = sc.read_h5ad(COUNTS_H5AD)
    print(f"Loading annotation object: {ANNOTATED_H5AD}")
    ann = sc.read_h5ad(ANNOTATED_H5AD, backed="r")

    if adata.n_obs != ann.n_obs:
        raise ValueError("Raw-count and annotated objects have different numbers of cells.")
    if not np.array_equal(adata.obs["sample"].astype(str).values, ann.obs["sample"].astype(str).values):
        raise ValueError("Raw-count and annotated objects are not aligned by sample order.")

    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata.obs["celltype"] = ann.obs["celltype"].astype(str).values
    adata.obs["age"] = adata.obs["sample"].map(lambda x: SAMPLE_MAP[x]["age"])
    adata.obs["condition"] = adata.obs["sample"].map(lambda x: SAMPLE_MAP[x]["condition"])
    adata.obs["age_condition"] = adata.obs["age"] + "_" + adata.obs["condition"]
    adata.obsm["X_umap"] = ann.obsm["X_umap"][:]
    ann.file.close()

    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def matrix_col_mean(matrix):
    if sparse.issparse(matrix):
        return np.asarray(matrix.mean(axis=0)).ravel()
    return np.asarray(matrix).mean(axis=0)


def expression_summaries(sub):
    frames = []
    for condition in ["Ctrl", "KO"]:
        cells = sub.obs["condition"].astype(str).values == condition
        if cells.sum() == 0:
            continue
        expr = sub[cells, :].X
        counts = sub[cells, :].layers["counts"]
        avg_expr = matrix_col_mean(expr)
        pct_detected = matrix_col_mean(counts > 0) * 100
        frames.append(
            pd.DataFrame(
                {
                    "gene": sub.var_names,
                    "condition": condition,
                    "n_cells": int(cells.sum()),
                    "avg_lognorm_expr": avg_expr,
                    "pct_detected": pct_detected,
                }
            )
        )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def run_de_for_group(adata, age, celltype):
    mask = (adata.obs["age"].astype(str) == age) & (adata.obs["celltype"].astype(str) == celltype)
    sub = adata[mask].copy()
    ctrl_n = int((sub.obs["condition"].astype(str) == "Ctrl").sum())
    ko_n = int((sub.obs["condition"].astype(str) == "KO").sum())
    status = {
        "age": age,
        "celltype": celltype,
        "ctrl_n": ctrl_n,
        "ko_n": ko_n,
        "tested_genes": 0,
        "status": "ok",
    }
    if ctrl_n < MIN_CELLS_PER_GROUP or ko_n < MIN_CELLS_PER_GROUP:
        status["status"] = f"skipped_min_cells_{MIN_CELLS_PER_GROUP}"
        return pd.DataFrame(), status

    counts = sub.layers["counts"]
    ctrl_mask = sub.obs["condition"].astype(str).values == "Ctrl"
    ko_mask = sub.obs["condition"].astype(str).values == "KO"
    ctrl_pct = np.asarray((counts[ctrl_mask, :] > 0).mean(axis=0)).ravel()
    ko_pct = np.asarray((counts[ko_mask, :] > 0).mean(axis=0)).ravel()
    keep = (ctrl_pct >= MIN_PCT_EXPRESSED) | (ko_pct >= MIN_PCT_EXPRESSED)
    sub = sub[:, keep].copy()
    status["tested_genes"] = int(sub.n_vars)
    if sub.n_vars == 0:
        status["status"] = "skipped_no_detected_genes"
        return pd.DataFrame(), status

    sub.obs["condition"] = pd.Categorical(sub.obs["condition"], categories=["Ctrl", "KO"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            sub,
            groupby="condition",
            groups=["KO"],
            reference="Ctrl",
            method="wilcoxon",
            tie_correct=True,
            pts=True,
        )
    de = sc.get.rank_genes_groups_df(sub, group="KO")
    de = de.rename(
        columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "log2fc_ko_vs_ctrl",
            "pvals": "p_value",
            "pvals_adj": "p_adj_bh",
            "pct_nz_group": "ko_pct_detected_rank_genes_groups",
            "pct_nz_reference": "ctrl_pct_detected_rank_genes_groups",
        }
    )
    summaries = expression_summaries(sub)
    ctrl = summaries[summaries["condition"] == "Ctrl"].drop(columns=["condition"]).rename(
        columns={
            "n_cells": "ctrl_n_cells",
            "avg_lognorm_expr": "ctrl_avg_lognorm_expr",
            "pct_detected": "ctrl_pct_detected",
        }
    )
    ko = summaries[summaries["condition"] == "KO"].drop(columns=["condition"]).rename(
        columns={
            "n_cells": "ko_n_cells",
            "avg_lognorm_expr": "ko_avg_lognorm_expr",
            "pct_detected": "ko_pct_detected",
        }
    )
    de = de.merge(ctrl, on="gene", how="left").merge(ko, on="gene", how="left")
    de.insert(0, "celltype", celltype)
    de.insert(0, "age", age)
    de["direction"] = "not_significant"
    de.loc[(de["p_adj_bh"] < FDR_CUTOFF) & (de["log2fc_ko_vs_ctrl"] >= LOGFC_CUTOFF), "direction"] = "up_in_KO"
    de.loc[(de["p_adj_bh"] < FDR_CUTOFF) & (de["log2fc_ko_vs_ctrl"] <= -LOGFC_CUTOFF), "direction"] = "down_in_KO"
    return de, status


def plot_volcano(de, age, celltype):
    fig, ax = plt.subplots(figsize=(7, 6))
    de = de.copy()
    de["neg_log10_fdr"] = -np.log10(de["p_adj_bh"].clip(lower=1e-300))
    palette = {
        "up_in_KO": "#D7191C",
        "down_in_KO": "#2C7BB6",
        "not_significant": "lightgray",
    }
    for direction, group in de.groupby("direction"):
        ax.scatter(
            group["log2fc_ko_vs_ctrl"],
            group["neg_log10_fdr"],
            s=9 if direction == "not_significant" else 16,
            c=palette.get(direction, "lightgray"),
            alpha=0.65 if direction == "not_significant" else 0.9,
            linewidth=0,
            label=direction.replace("_", " "),
        )
    ax.axvline(LOGFC_CUTOFF, color="black", linewidth=0.8, linestyle="--")
    ax.axvline(-LOGFC_CUTOFF, color="black", linewidth=0.8, linestyle="--")
    ax.axhline(-np.log10(FDR_CUTOFF), color="black", linewidth=0.8, linestyle="--")
    top = de[de["direction"] != "not_significant"].nsmallest(10, "p_adj_bh")
    for _, row in top.iterrows():
        ax.text(row["log2fc_ko_vs_ctrl"], row["neg_log10_fdr"], row["gene"], fontsize=8)
    ax.set_xlabel("log2FC (KO vs Ctrl)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(f"{celltype} DEGs at {age}")
    ax.legend(frameon=True, fontsize=9, loc="upper right")
    savefig(fig, f"volcano_{age}_{celltype}_ko_vs_ctrl")


def plot_top_genes_heatmap(cone_de, age):
    sig = cone_de[cone_de["direction"] != "not_significant"].copy()
    if sig.empty:
        return
    top_down = sig[sig["direction"] == "down_in_KO"].nsmallest(12, "p_adj_bh")
    top_up = sig[sig["direction"] == "up_in_KO"].nsmallest(12, "p_adj_bh")
    top = pd.concat([top_down, top_up], ignore_index=True)
    if top.empty:
        return
    plot_df = top[["gene", "ctrl_avg_lognorm_expr", "ko_avg_lognorm_expr", "log2fc_ko_vs_ctrl", "direction"]].copy()
    plot_df = plot_df.drop_duplicates("gene")
    mat = plot_df.set_index("gene")[["ctrl_avg_lognorm_expr", "ko_avg_lognorm_expr"]]
    fig, ax = plt.subplots(figsize=(5.8, max(4, 0.28 * len(mat) + 1.5)))
    sns.heatmap(
        mat,
        cmap="viridis",
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Avg log-normalized expression"},
        ax=ax,
    )
    ax.set_title(f"Top cone DEGs at {age}")
    ax.set_xlabel("")
    ax.set_ylabel("")
    savefig(fig, f"top_cone_degs_heatmap_{age}", 5.8, max(4, 0.28 * len(mat) + 1.5))


def plot_deg_counts(counts):
    plot_df = counts[counts["status"] == "ok"].copy()
    if plot_df.empty:
        return
    long = plot_df.melt(
        id_vars=["age", "celltype"],
        value_vars=["n_up_in_KO", "n_down_in_KO"],
        var_name="direction",
        value_name="n_degs",
    )
    long["direction"] = long["direction"].replace({"n_up_in_KO": "Up in KO", "n_down_in_KO": "Down in KO"})
    long["celltype"] = pd.Categorical(long["celltype"], categories=CELLTYPE_ORDER, ordered=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    for ax, age in zip(axes, AGE_ORDER):
        sub = long[long["age"] == age]
        sns.barplot(
            data=sub,
            x="n_degs",
            y="celltype",
            hue="direction",
            palette={"Up in KO": "#D7191C", "Down in KO": "#2C7BB6"},
            ax=ax,
        )
        ax.set_title(age)
        ax.set_xlabel("Number of DEGs")
        ax.set_ylabel("")
        ax.grid(axis="y", visible=False)
    axes[-1].legend(title="")
    fig.suptitle("DEG counts by cell type and age", fontweight="bold")
    savefig(fig, "deg_counts_by_celltype_age_direction", 13, 6)


def plot_pathway_results(enrich_df, stem, title):
    if enrich_df.empty:
        return
    df = enrich_df.copy().head(PATHWAY_TOP_N)
    p_col = "Adjusted P-value"
    term_col = "Term"
    df["neg_log10_fdr"] = -np.log10(df[p_col].clip(lower=1e-300))
    df[term_col] = df[term_col].str.replace(r" \\(GO:\\d+\\)", "", regex=True)
    df = df.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, max(4, 0.35 * len(df) + 1.5)))
    ax.barh(df[term_col], df["neg_log10_fdr"], color="#4E79A7")
    ax.set_xlabel("-log10 adjusted p value")
    ax.set_ylabel("")
    ax.set_title(title)
    savefig(fig, stem, 8.5, max(4, 0.35 * len(df) + 1.5))


def run_enrichment(cone_de):
    if os.environ.get("RUN_ENRICHR", "0") != "1":
        status = pd.DataFrame(
            [
                {
                    "age": "all",
                    "direction": "all",
                    "n_genes": np.nan,
                    "status": "skipped_online_enrichr_not_authorized",
                }
            ]
        )
        status.to_csv(PATHWAY_DIR / "pathway_status.csv", index=False)
        return status.to_dict("records")

    if gp is None:
        pd.DataFrame([{"status": "gseapy_not_installed"}]).to_csv(PATHWAY_DIR / "pathway_status.csv", index=False)
        return []

    statuses = []
    libraries = ["GO_Biological_Process_2023"]
    for age in AGE_ORDER:
        age_de = cone_de[cone_de["age"] == age].copy()
        for direction in ["up_in_KO", "down_in_KO"]:
            genes = age_de.loc[age_de["direction"] == direction, "gene"].dropna().astype(str).unique().tolist()
            status = {"age": age, "direction": direction, "n_genes": len(genes), "status": "not_run"}
            if len(genes) < 5:
                status["status"] = "skipped_fewer_than_5_genes"
                statuses.append(status)
                continue
            try:
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=libraries,
                    organism="Mouse",
                    outdir=None,
                    cutoff=1.0,
                    verbose=False,
                )
                res = enr.results.copy()
                res.insert(0, "direction", direction)
                res.insert(0, "age", age)
                out_csv = PATHWAY_DIR / f"cones_{age}_{direction}_GO_BP_2023_enrichr.csv"
                res.to_csv(out_csv, index=False)
                sig = res[res["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")
                plot_pathway_results(
                    sig,
                    f"pathway_cones_{age}_{direction}_GO_BP_2023",
                    f"Cones {age}: {direction.replace('_', ' ')} GO BP enrichment",
                )
                status["status"] = "ok"
                status["n_enriched_fdr_0_05"] = int(sig.shape[0])
            except Exception as exc:
                status["status"] = f"failed: {exc}"
            statuses.append(status)
    pd.DataFrame(statuses).to_csv(PATHWAY_DIR / "pathway_status.csv", index=False)
    return statuses


def write_summary(deg_counts, pathway_status):
    lines = [
        "Analysis 04: DEGs and pathway enrichment",
        "",
        "Design:",
        "- Main input: Scanpy raw-count object photoreceptors.h5ad plus annotations from annotated_photoreceptors.h5ad.",
        "- Comparisons: KO vs Ctrl within each age and cell type.",
        "- Sample mapping: 15dayS1=P15 Ctrl, 15dayS2=P15 KO, 35dayS1=P35 Ctrl, 35dayS2=P35 KO.",
        "",
        "Important caveat:",
        "- Each age/condition has one sample, so these are exploratory cell-level Wilcoxon DE results, not replicate-aware sample-level DE tests.",
        "",
        f"DEG threshold: FDR < {FDR_CUTOFF} and absolute log2FC >= {LOGFC_CUTOFF}.",
        f"Gene filter before DE: detected in at least {MIN_PCT_EXPRESSED * 100:.0f}% of cells in Ctrl or KO.",
        "",
        "Top non-cone DEG burdens:",
    ]
    ok = deg_counts[deg_counts["status"] == "ok"].copy()
    ok_noncone = ok[ok["celltype"] != "Cones"].sort_values(["age", "n_total_degs"], ascending=[True, False])
    for age in AGE_ORDER:
        top = ok_noncone[ok_noncone["age"] == age].head(5)
        lines.append(f"- {age}: " + "; ".join([f"{r.celltype}={int(r.n_total_degs)}" for r in top.itertuples()]))
    if pathway_status:
        lines.append("")
        lines.append("Pathway status:")
        for st in pathway_status:
            n_genes = st.get("n_genes")
            gene_text = "not run" if pd.isna(n_genes) else f"{int(n_genes)} genes"
            lines.append(f"- {st['age']} {st['direction']}: {st['status']} ({gene_text})")
    (OUT_DIR / "analysis_04_summary.txt").write_text("\n".join(lines) + "\n")


def main():
    for directory in [OUT_DIR, TABLE_DIR, FIG_DIR, PATHWAY_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")

    adata = read_and_prepare()
    all_de = []
    status_rows = []
    for age in AGE_ORDER:
        for celltype in CELLTYPE_ORDER:
            print(f"DE: {age} {celltype}")
            de, status = run_de_for_group(adata, age, celltype)
            status_rows.append(status)
            if not de.empty:
                out_name = f"{age}_{celltype}_KO_vs_Ctrl_degs.csv".replace("/", "_")
                de.to_csv(TABLE_DIR / out_name, index=False)
                all_de.append(de)
                if celltype == "Cones":
                    plot_volcano(de, age, celltype)
                    plot_top_genes_heatmap(de, age)

    de_all = pd.concat(all_de, ignore_index=True) if all_de else pd.DataFrame()
    de_all.to_csv(TABLE_DIR / "all_celltypes_all_de_results.csv", index=False)
    status_df = pd.DataFrame(status_rows)
    if not de_all.empty:
        counts = (
            de_all.groupby(["age", "celltype"], observed=True)
            .agg(
                n_up_in_KO=("direction", lambda x: int((x == "up_in_KO").sum())),
                n_down_in_KO=("direction", lambda x: int((x == "down_in_KO").sum())),
                n_total_degs=("direction", lambda x: int((x != "not_significant").sum())),
                n_tested_genes=("gene", "size"),
                ctrl_n=("ctrl_n_cells", "max"),
                ko_n=("ko_n_cells", "max"),
            )
            .reset_index()
        )
        counts = counts.merge(
            status_df[["age", "celltype", "ctrl_n", "ko_n", "status"]],
            on=["age", "celltype"],
            how="outer",
            suffixes=("", "_status"),
        )
        counts["ctrl_n"] = counts["ctrl_n"].fillna(counts["ctrl_n_status"])
        counts["ko_n"] = counts["ko_n"].fillna(counts["ko_n_status"])
        counts = counts.drop(columns=["ctrl_n_status", "ko_n_status"])
    else:
        counts = status_df.copy()
        counts["n_total_degs"] = 0
        counts["n_up_in_KO"] = 0
        counts["n_down_in_KO"] = 0
    counts.to_csv(TABLE_DIR / "all_celltypes_deg_counts_by_age_celltype.csv", index=False)
    status_df.to_csv(TABLE_DIR / "de_status_by_age_celltype.csv", index=False)
    plot_deg_counts(counts)

    cone_de = de_all[de_all["celltype"] == "Cones"].copy() if not de_all.empty else pd.DataFrame()
    cone_de.to_csv(TABLE_DIR / "cones_all_de_results.csv", index=False)
    cone_sig = cone_de[cone_de["direction"] != "not_significant"].copy()
    cone_sig.to_csv(TABLE_DIR / "cones_significant_degs_fdr0.05_abslogfc0.25.csv", index=False)
    gene_list_dir = PATHWAY_DIR / "input_gene_lists"
    gene_list_dir.mkdir(parents=True, exist_ok=True)
    for age in AGE_ORDER:
        for direction in ["up_in_KO", "down_in_KO"]:
            genes = cone_sig.loc[
                (cone_sig["age"] == age) & (cone_sig["direction"] == direction),
                "gene",
            ].dropna().astype(str).drop_duplicates().tolist()
            (gene_list_dir / f"cones_{age}_{direction}_genes.txt").write_text("\n".join(genes) + ("\n" if genes else ""))
    cone_counts = counts[counts["celltype"] == "Cones"].copy()
    cone_counts.to_csv(TABLE_DIR / "cones_deg_counts_by_age_direction.csv", index=False)

    pathway_status = run_enrichment(cone_de)
    write_summary(counts, pathway_status)
    print(f"Done. Outputs written to: {OUT_DIR}")


if __name__ == "__main__":
    main()
