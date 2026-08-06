#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
scCODA analysis for differential cell-type composition
Extended version with publication-style plots

Usage example:
    python run_sccoda_from_metadata_with_figures.py \
        --meta cell_metadata.csv \
        --sample_col orig.ident \
        --group_col group \
        --celltype_col cell.cluster \
        --control_level homing \
        --case_level tumbler \
        --outdir sccoda_cerebellum \
        --reference automatic \
        --plot

Input metadata table must contain at least:
    sample_col, group_col, celltype_col
"""

import os
import re
import sys
import argparse
import warnings
import numpy as np
import pandas as pd

# silence some tensorflow warnings if present
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
warnings.filterwarnings("ignore")

try:
    import anndata as ad
except ImportError:
    raise ImportError("Please install anndata first: pip install anndata")

try:
    from sccoda.util import comp_ana as mod
except ImportError:
    raise ImportError(
        "Cannot import scCODA. Please install it first: pip install sccoda"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run scCODA on sample-level cell type count table from metadata."
    )
    parser.add_argument("--meta", required=True, help="Metadata CSV/TSV file")
    parser.add_argument("--sample_col", required=True, help="Sample column name, e.g. orig.ident")
    parser.add_argument("--group_col", required=True, help="Group column name, e.g. group")
    parser.add_argument("--celltype_col", required=True, help="Cell type column name, e.g. cell.cluster or cell_type")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--control_level", required=True, help="Control level, e.g. homing")
    parser.add_argument("--case_level", required=True, help="Case level, e.g. tumbler")
    parser.add_argument("--reference", default="automatic",
                        help='Reference cell type name, index, or "automatic"')
    parser.add_argument("--min_total_cells", type=int, default=50,
                        help="Minimum total cells per sample to keep")
    parser.add_argument("--min_celltype_total", type=int, default=10,
                        help="Minimum total counts across all samples for a cell type to keep")
    parser.add_argument("--min_presence", type=int, default=2,
                        help="Minimum number of samples in which a cell type must appear (>0)")
    parser.add_argument("--covariates", default="",
                        help="Additional covariates separated by commas, e.g. batch,sex")
    parser.add_argument("--delimiter", default=",",
                        help='Input delimiter: "," for CSV, "\\t" for TSV')
    parser.add_argument("--seed", type=int, default=1234, help="Random seed")
    parser.add_argument("--draws", type=int, default=20000,
                        help="MCMC draws (posterior samples)")
    parser.add_argument("--warmup", type=int, default=5000,
                        help="MCMC burnin / warmup")
    parser.add_argument("--plot", action="store_true",
                        help="Save publication-style plots")
    return parser.parse_args()


def read_table(path, delimiter=","):
    if delimiter == r"\t":
        delimiter = "\t"
    df = pd.read_csv(path, sep=delimiter)
    return df


def build_count_table(df, sample_col, group_col, celltype_col):
    sample_meta = (
        df[[sample_col, group_col]]
        .drop_duplicates()
        .copy()
    )

    count_tab = (
        df.groupby([sample_col, celltype_col])
          .size()
          .unstack(fill_value=0)
          .copy()
    )

    sample_meta = sample_meta.set_index(sample_col).loc[count_tab.index].reset_index()
    return count_tab, sample_meta


def filter_samples_and_celltypes(count_tab, sample_meta,
                                 min_total_cells=50,
                                 min_celltype_total=10,
                                 min_presence=2):
    sample_totals = count_tab.sum(axis=1)
    keep_samples = sample_totals >= min_total_cells
    count_tab = count_tab.loc[keep_samples].copy()
    sample_meta = sample_meta.set_index(sample_meta.columns[0]).loc[count_tab.index].reset_index()

    celltype_totals = count_tab.sum(axis=0)
    celltype_presence = (count_tab > 0).sum(axis=0)
    keep_celltypes = (celltype_totals >= min_celltype_total) & (celltype_presence >= min_presence)
    count_tab = count_tab.loc[:, keep_celltypes].copy()

    return count_tab, sample_meta


def add_additional_covariates(df, sample_meta, sample_col, covariates):
    if not covariates:
        return sample_meta

    cov_list = [x.strip() for x in covariates.split(",") if x.strip()]
    if len(cov_list) == 0:
        return sample_meta

    extra = df[[sample_col] + cov_list].drop_duplicates().copy()

    for c in cov_list:
        n_dup = extra.groupby(sample_col)[c].nunique()
        bad = n_dup[n_dup > 1]
        if len(bad) > 0:
            raise ValueError(
                f"Covariate '{c}' is not unique within samples. "
                f"Please ensure one value per sample."
            )

    sample_meta = sample_meta.merge(extra, on=sample_col, how="left")
    return sample_meta


def make_anndata(count_tab, sample_meta, sample_col):
    adata = ad.AnnData(X=count_tab.values.astype(np.int64))
    adata.obs_names = count_tab.index.astype(str)
    adata.var_names = count_tab.columns.astype(str)

    obs = sample_meta.set_index(sample_col).loc[count_tab.index].copy()
    adata.obs = obs
    return adata


def make_formula(group_col, additional_covariates=""):
    cov_list = [x.strip() for x in additional_covariates.split(",") if x.strip()]
    if len(cov_list) == 0:
        return group_col
    return " + ".join([group_col] + cov_list)


def safe_name(x):
    return str(x).replace("/", "_").replace(" ", "_")


def save_text(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(str(text))


def summarize_proportions(prop_long, group_col, sample_col, control_level, case_level):
    def sem(x):
        x = pd.Series(x).dropna()
        if len(x) <= 1:
            return np.nan
        return x.std(ddof=1) / np.sqrt(len(x))

    summary = (
        prop_long.groupby([group_col, "cell_type"])["proportion"]
        .agg(["mean", "median", "std", "min", "max", sem, "count"])
        .reset_index()
    )

    summary_wide = summary.pivot(index="cell_type", columns=group_col)
    summary_wide.columns = [f"{a}_{b}" for a, b in summary_wide.columns]
    summary_wide = summary_wide.reset_index()

    mean_ctrl = f"mean_{control_level}"
    mean_case = f"mean_{case_level}"
    med_ctrl = f"median_{control_level}"
    med_case = f"median_{case_level}"

    eps = 1e-9
    if mean_ctrl in summary_wide.columns and mean_case in summary_wide.columns:
        summary_wide["log2FC_case_vs_control"] = np.log2(
            (summary_wide[mean_case] + eps) / (summary_wide[mean_ctrl] + eps)
        )

    if med_ctrl in summary_wide.columns and med_case in summary_wide.columns:
        summary_wide["median_diff_case_minus_control"] = (
            summary_wide[med_case] - summary_wide[med_ctrl]
        )

    return summary_wide


def parse_sccoda_summary_to_tables(summary_text):
    """
    Parse scCODA text summary like:

    Intercepts:
                                    Final Parameter  Expected Sample
    Cell Type
    astrocytes                                3.039      1345.595792

    Effects:
                                                     Final Parameter  ...  log2-fold change
    Covariate        Cell Type
    group[T.tumbler] astrocytes                         0.0          ... 0.0
    """
    lines = str(summary_text).splitlines()
    lines = [x.rstrip("\n") for x in lines]

    intercept_rows = []
    effect_rows = []

    section = None
    for line in lines:
        s = line.strip()
        if not s:
            continue
        if s.startswith("Intercepts:"):
            section = "intercepts"
            continue
        if s.startswith("Effects:"):
            section = "effects"
            continue

        if section == "intercepts":
            if s.startswith("Cell Type") or s.startswith("Final Parameter") or s.startswith("Expected Sample"):
                continue
            m = re.match(r"^(.*?)\s+([-+]?\d*\.?\d+)\s+([-+]?\d*\.?\d+)$", s)
            if m:
                cell_type = m.group(1).strip()
                final_parameter = float(m.group(2))
                expected_sample = float(m.group(3))
                intercept_rows.append({
                    "cell_type": cell_type,
                    "final_parameter": final_parameter,
                    "expected_sample": expected_sample
                })

        elif section == "effects":
            if s.startswith("Covariate") or s.startswith("Final Parameter") or "log2-fold change" in s:
                continue

            # try to parse the last numeric value as log2-fold change
            # and the preceding numeric as final parameter
            parts = re.split(r"\s{2,}", s)
            parts = [p for p in parts if p != ""]

            if len(parts) >= 3:
                # e.g. group[T.tumbler], astrocytes, 0.0, ..., 0.0
                covariate = parts[0].strip()
                cell_type = parts[1].strip()

                nums = []
                for p in parts[2:]:
                    try:
                        nums.append(float(p))
                    except Exception:
                        pass

                if len(nums) >= 1:
                    row = {
                        "covariate": covariate,
                        "cell_type": cell_type,
                        "final_parameter": nums[0]
                    }
                    if len(nums) >= 2:
                        row["log2_fold_change"] = nums[-1]
                    effect_rows.append(row)

    intercept_df = pd.DataFrame(intercept_rows)
    effect_df = pd.DataFrame(effect_rows)

    return intercept_df, effect_df


def try_extract_sccoda_tables(result, summary_text, outdir):
    intercept_df = pd.DataFrame()
    effect_df = pd.DataFrame()

    if hasattr(result, "summary_prepare"):
        try:
            prep = result.summary_prepare()
            if isinstance(prep, pd.DataFrame):
                prep.to_csv(os.path.join(outdir, "summary_prepare.csv"), index=True)
        except Exception:
            pass

    if hasattr(result, "credible_effects"):
        try:
            cred = result.credible_effects()
            if isinstance(cred, pd.DataFrame):
                cred.to_csv(os.path.join(outdir, "credible_effects.csv"), index=True)
            else:
                pd.DataFrame(cred).to_csv(os.path.join(outdir, "credible_effects.csv"), index=True)
        except Exception:
            pass

    try:
        intercept_df, effect_df = parse_sccoda_summary_to_tables(summary_text)
    except Exception:
        pass

    if not intercept_df.empty:
        intercept_df.to_csv(os.path.join(outdir, "sccoda_intercepts_parsed.csv"), index=False)
    if not effect_df.empty:
        effect_df.to_csv(os.path.join(outdir, "sccoda_effects_parsed.csv"), index=False)

    return intercept_df, effect_df


def make_color_map(celltypes):
    import matplotlib.pyplot as plt

    n = len(celltypes)
    cmap1 = plt.cm.get_cmap("tab20", 20)
    cmap2 = plt.cm.get_cmap("tab20b", 20)
    cmap3 = plt.cm.get_cmap("tab20c", 20)

    colors = []
    for cmap in [cmap1, cmap2, cmap3]:
        for i in range(cmap.N):
            colors.append(cmap(i))
    colors = colors[:n]

    return dict(zip(celltypes, colors))


def plot_stacked_bar(prop_tab, sample_meta, sample_col, group_col,
                     control_level, case_level, plot_dir):
    import matplotlib.pyplot as plt

    meta = sample_meta.copy()
    meta[group_col] = pd.Categorical(
        meta[group_col],
        categories=[control_level, case_level],
        ordered=True
    )
    meta = meta.sort_values([group_col, sample_col])

    prop_tab = prop_tab.loc[meta[sample_col].values]
    celltypes = list(prop_tab.columns)
    color_map = make_color_map(celltypes)

    fig, ax = plt.subplots(figsize=(max(8, 0.6 * prop_tab.shape[0] + 3), 5.5))

    bottoms = np.zeros(prop_tab.shape[0])
    x = np.arange(prop_tab.shape[0])

    for ct in celltypes:
        vals = prop_tab[ct].values
        ax.bar(x, vals, bottom=bottoms, label=ct, color=color_map[ct], width=0.85)
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(prop_tab.index, rotation=90, fontsize=8)
    ax.set_ylabel("Cell proportion")
    ax.set_xlabel("Sample")
    ax.set_ylim(0, 1.0)
    ax.set_title("Sample-level cell-type composition")

    group_vals = meta[group_col].astype(str).values
    boundary = None
    for i in range(1, len(group_vals)):
        if group_vals[i] != group_vals[i-1]:
            boundary = i - 0.5
            break
    if boundary is not None:
        ax.axvline(boundary, color="black", linestyle="--", linewidth=1)

    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "stacked_barplot_sample_composition.pdf"))
    plt.savefig(os.path.join(plot_dir, "stacked_barplot_sample_composition.png"), dpi=300)
    plt.close(fig)


def plot_boxplot_facets(prop_long, group_col, control_level, case_level, plot_dir):
    import matplotlib.pyplot as plt

    celltypes = list(prop_long["cell_type"].drop_duplicates())
    n = len(celltypes)
    ncol = 4
    nrow = int(np.ceil(n / ncol))

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 3.6 * nrow), squeeze=False)
    axes = axes.flatten()

    for i, ct in enumerate(celltypes):
        ax = axes[i]
        sub = prop_long[prop_long["cell_type"] == ct].copy()

        groups = [control_level, case_level]
        data = [
            sub.loc[sub[group_col] == g, "proportion"].values
            for g in groups
        ]

        ax.boxplot(
            data,
            labels=groups,
            widths=0.55,
            patch_artist=True,
            boxprops=dict(facecolor="white", edgecolor="black"),
            medianprops=dict(color="black"),
            whiskerprops=dict(color="black"),
            capprops=dict(color="black")
        )

        for j, g in enumerate(groups, start=1):
            y = sub.loc[sub[group_col] == g, "proportion"].values
            if len(y) > 0:
                x = np.random.normal(j, 0.05, size=len(y))
                ax.scatter(x, y, s=28, alpha=0.85)

        ax.set_title(ct, fontsize=10)
        ax.set_ylabel("Proportion")
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Cell-type proportions by group", fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "proportion_boxplot_by_celltype.pdf"))
    plt.savefig(os.path.join(plot_dir, "proportion_boxplot_by_celltype.png"), dpi=300)
    plt.close(fig)


def plot_heatmap(prop_tab, sample_meta, sample_col, group_col,
                 control_level, case_level, plot_dir):
    import matplotlib.pyplot as plt

    meta = sample_meta.copy()
    meta[group_col] = pd.Categorical(
        meta[group_col],
        categories=[control_level, case_level],
        ordered=True
    )
    meta = meta.sort_values([group_col, sample_col])

    mat = prop_tab.loc[meta[sample_col].values].copy()
    mat_z = mat.copy()

    for col in mat_z.columns:
        mu = mat_z[col].mean()
        sd = mat_z[col].std(ddof=0)
        if sd == 0:
            mat_z[col] = 0.0
        else:
            mat_z[col] = (mat_z[col] - mu) / sd

    fig, ax = plt.subplots(figsize=(max(8, 0.45 * mat_z.shape[1] + 2), max(5, 0.4 * mat_z.shape[0] + 2)))
    im = ax.imshow(mat_z.values, aspect="auto", interpolation="nearest")

    ax.set_xticks(np.arange(mat_z.shape[1]))
    ax.set_xticklabels(mat_z.columns, rotation=90)
    ax.set_yticks(np.arange(mat_z.shape[0]))
    ax.set_yticklabels(mat_z.index)
    ax.set_title("Composition heatmap (column z-score)")
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Sample")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Z-score")

    group_vals = meta[group_col].astype(str).values
    boundary = None
    for i in range(1, len(group_vals)):
        if group_vals[i] != group_vals[i-1]:
            boundary = i - 0.5
            break
    if boundary is not None:
        ax.axhline(boundary, color="white", linestyle="--", linewidth=1.5)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "composition_heatmap_zscore.pdf"))
    plt.savefig(os.path.join(plot_dir, "composition_heatmap_zscore.png"), dpi=300)
    plt.close(fig)


def plot_sccoda_effects(effect_df, plot_dir, case_level):
    import matplotlib.pyplot as plt

    if effect_df is None or effect_df.empty:
        print("No parsed scCODA effect table available; skip effect plot.")
        return

    df = effect_df.copy()

    if "log2_fold_change" not in df.columns:
        df["log2_fold_change"] = df["final_parameter"]

    df = df.sort_values("log2_fold_change")
    y = np.arange(df.shape[0])

    fig, ax = plt.subplots(figsize=(7, max(4, 0.5 * df.shape[0] + 1.5)))
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.scatter(df["log2_fold_change"], y, s=55)

    for i, (_, row) in enumerate(df.iterrows()):
        ax.plot([0, row["log2_fold_change"]], [i, i], linewidth=1)

    ax.set_yticks(y)
    ax.set_yticklabels(df["cell_type"])
    ax.set_xlabel(f"scCODA effect (approx. log2 fold change)\n{case_level} vs control")
    ax.set_ylabel("Cell type")
    ax.set_title("scCODA effect plot")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "sccoda_effect_plot.pdf"))
    plt.savefig(os.path.join(plot_dir, "sccoda_effect_plot.png"), dpi=300)
    plt.close(fig)


def main():
    args = parse_args()
    np.random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    print(">>> Reading metadata ...")
    df = read_table(args.meta, args.delimiter)

    required_cols = [args.sample_col, args.group_col, args.celltype_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in metadata: {missing}")

    df = df.copy()
    df[args.sample_col] = df[args.sample_col].astype(str).str.strip()
    df[args.group_col] = df[args.group_col].astype(str).str.strip()
    df[args.celltype_col] = df[args.celltype_col].astype(str).str.strip()

    df = df[df[args.group_col].isin([args.control_level, args.case_level])].copy()
    if df.empty:
        raise ValueError("No cells left after filtering control/case groups.")

    print(">>> Building count table ...")
    count_tab, sample_meta = build_count_table(
        df=df,
        sample_col=args.sample_col,
        group_col=args.group_col,
        celltype_col=args.celltype_col
    )

    print(">>> Adding extra covariates (if any) ...")
    sample_meta = add_additional_covariates(
        df=df,
        sample_meta=sample_meta,
        sample_col=args.sample_col,
        covariates=args.covariates
    )

    print(">>> Filtering low-depth samples and rare cell types ...")
    count_tab, sample_meta = filter_samples_and_celltypes(
        count_tab=count_tab,
        sample_meta=sample_meta,
        min_total_cells=args.min_total_cells,
        min_celltype_total=args.min_celltype_total,
        min_presence=args.min_presence
    )

    if count_tab.shape[0] < 4:
        raise ValueError(
            f"Too few samples after filtering: {count_tab.shape[0]}. "
            f"scCODA needs biological replicates."
        )
    if count_tab.shape[1] < 2:
        raise ValueError(
            f"Too few cell types after filtering: {count_tab.shape[1]}."
        )

    sample_meta[args.group_col] = pd.Categorical(
        sample_meta[args.group_col],
        categories=[args.control_level, args.case_level],
        ordered=True
    )

    count_tab.to_csv(os.path.join(args.outdir, "cell_count_matrix.csv"))
    sample_meta.to_csv(os.path.join(args.outdir, "sample_metadata.csv"), index=False)

    print(">>> Creating AnnData ...")
    adata = make_anndata(count_tab, sample_meta, args.sample_col)
    adata.write(os.path.join(args.outdir, "sccoda_input.h5ad"))

    formula = make_formula(args.group_col, args.covariates)
    print(f">>> Formula: {formula}")
    print(f">>> Reference: {args.reference}")

    print(">>> Running scCODA model ...")
    model = mod.CompositionalAnalysis(
        adata,
        formula=formula,
        reference_cell_type=args.reference
    )

    result = model.sample_hmc(
        num_results=args.draws,
        num_burnin=args.warmup
    )

    print(">>> Saving summaries ...")
    summary_text = result.summary()
    save_text(os.path.join(args.outdir, "sccoda_summary.txt"), summary_text)

    intercept_df, effect_df = try_extract_sccoda_tables(
        result=result,
        summary_text=summary_text,
        outdir=args.outdir
    )

    prop_tab = count_tab.div(count_tab.sum(axis=1), axis=0)
    prop_long = (
        prop_tab.reset_index()
               .melt(id_vars=count_tab.index.name or "index",
                     var_name="cell_type",
                     value_name="proportion")
    )
    prop_long.rename(columns={prop_long.columns[0]: args.sample_col}, inplace=True)
    prop_long = prop_long.merge(sample_meta, on=args.sample_col, how="left")
    prop_long.to_csv(os.path.join(args.outdir, "sample_level_proportions.csv"), index=False)

    summary_prop = summarize_proportions(
        prop_long=prop_long,
        group_col=args.group_col,
        sample_col=args.sample_col,
        control_level=args.control_level,
        case_level=args.case_level
    )
    summary_prop.to_csv(os.path.join(args.outdir, "celltype_proportion_summary.csv"), index=False)

    if args.plot:
        print(">>> Generating plots ...")
        plot_dir = os.path.join(args.outdir, "plots")
        os.makedirs(plot_dir, exist_ok=True)

        plot_stacked_bar(
            prop_tab=prop_tab,
            sample_meta=sample_meta,
            sample_col=args.sample_col,
            group_col=args.group_col,
            control_level=args.control_level,
            case_level=args.case_level,
            plot_dir=plot_dir
        )

        plot_boxplot_facets(
            prop_long=prop_long,
            group_col=args.group_col,
            control_level=args.control_level,
            case_level=args.case_level,
            plot_dir=plot_dir
        )

        plot_heatmap(
            prop_tab=prop_tab,
            sample_meta=sample_meta,
            sample_col=args.sample_col,
            group_col=args.group_col,
            control_level=args.control_level,
            case_level=args.case_level,
            plot_dir=plot_dir
        )

        plot_sccoda_effects(
            effect_df=effect_df,
            plot_dir=plot_dir,
            case_level=args.case_level
        )

    print(">>> Done.")
    print(f">>> Outputs saved in: {args.outdir}")


if __name__ == "__main__":
    main()
