import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_and_prepare_lsbl(file_path):
    """
    读取 LSBL 文件并计算 Z-score
    预期至少包含列: CHR, START, END, LSBL
    """
    print(f"🚀 正在加载: {file_path}")
    df = pd.read_csv(file_path, sep=r"\s+")

    required_cols = ["CHR", "START", "END", "LSBL"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"文件缺少必要列: {missing}")

    df = df[required_cols].copy()

    df["LSBL"] = pd.to_numeric(df["LSBL"], errors="coerce")
    df["START"] = pd.to_numeric(df["START"], errors="coerce")
    df["END"] = pd.to_numeric(df["END"], errors="coerce")
    df = df.dropna(subset=["LSBL", "START", "END"])

    df["START"] = df["START"].astype(int)
    df["END"] = df["END"].astype(int)
    df["CHR"] = df["CHR"].astype(str)

    # 按全基因组标准化；若你希望按染色体标准化，可改这里
    lsbl_std = df["LSBL"].std(ddof=1)
    if lsbl_std == 0 or pd.isna(lsbl_std):
        raise ValueError(f"{file_path} 的 LSBL 标准差为 0，无法计算 Z-score。")

    df["z"] = (df["LSBL"] - df["LSBL"].mean()) / lsbl_std
    df = df.sort_values(["CHR", "START", "END"]).reset_index(drop=True)

    return df


def align_datasets(df1, df2, name1="xj", name2="pr"):
    """
    按 CHR/START/END 对齐两个数据集
    """
    merged = pd.merge(
        df1[["CHR", "START", "END", "z"]],
        df2[["CHR", "START", "END", "z"]],
        on=["CHR", "START", "END"],
        how="inner",
        suffixes=(f"_{name1}", f"_{name2}")
    )

    merged = merged.sort_values(["CHR", "START", "END"]).reset_index(drop=True)
    if merged.empty:
        raise ValueError("两个数据集对齐后为空，请检查窗口定义是否一致。")

    return merged


def mask_to_peaks(chr_df, sig_mask):
    """
    将单条染色体上的显著窗口布尔向量转换为 peak 区间
    连续 True 窗口会被合并为一个 peak
    """
    if len(chr_df) != len(sig_mask):
        raise ValueError("chr_df 与 sig_mask 长度不一致。")

    peaks = []
    if len(sig_mask) == 0:
        return peaks

    starts = chr_df["START"].to_numpy()
    ends = chr_df["END"].to_numpy()
    chrom = chr_df["CHR"].iloc[0]

    in_peak = False
    peak_start = None
    peak_end = None

    for i, is_sig in enumerate(sig_mask):
        if is_sig:
            if not in_peak:
                in_peak = True
                peak_start = starts[i]
                peak_end = ends[i]
            else:
                # 连续显著窗口合并；直接延伸终点
                peak_end = max(peak_end, ends[i])
        else:
            if in_peak:
                peaks.append((chrom, peak_start, peak_end))
                in_peak = False
                peak_start = None
                peak_end = None

    if in_peak:
        peaks.append((chrom, peak_start, peak_end))

    return peaks


def build_peak_table(merged_df, sig_col):
    """
    从全基因组显著性布尔列构建 peak 表
    返回 DataFrame: CHR, PEAK_START, PEAK_END
    """
    peak_records = []

    for chrom, chr_df in merged_df.groupby("CHR", sort=False):
        sig_mask = chr_df[sig_col].to_numpy(dtype=bool)
        peaks = mask_to_peaks(chr_df, sig_mask)
        peak_records.extend(peaks)

    peak_df = pd.DataFrame(peak_records, columns=["CHR", "PEAK_START", "PEAK_END"])
    return peak_df


def count_overlapping_peaks(peak_a, peak_b):
    """
    统计 peak_a 中有多少个 peak 与 peak_b 至少发生一次重叠
    按染色体分别用双指针计算
    """
    if peak_a.empty or peak_b.empty:
        return 0

    total_overlap = 0

    chroms = sorted(set(peak_a["CHR"]).intersection(set(peak_b["CHR"])))
    for chrom in chroms:
        a = peak_a.loc[peak_a["CHR"] == chrom, ["PEAK_START", "PEAK_END"]].sort_values("PEAK_START").to_numpy()
        b = peak_b.loc[peak_b["CHR"] == chrom, ["PEAK_START", "PEAK_END"]].sort_values("PEAK_START").to_numpy()

        j = 0
        for a_start, a_end in a:
            while j < len(b) and b[j][1] < a_start:
                j += 1

            k = j
            overlapped = False
            while k < len(b) and b[k][0] <= a_end:
                if b[k][1] >= a_start:
                    overlapped = True
                    break
                k += 1

            if overlapped:
                total_overlap += 1

    return total_overlap


def circular_shift_boolean(arr, shift):
    """
    对布尔数组做循环位移
    """
    if len(arr) == 0:
        return arr
    shift = shift % len(arr)
    if shift == 0:
        return arr.copy()
    return np.roll(arr, shift)


def permuted_pr_peak_table_by_chr(merged_df, pr_sig_col, rng):
    """
    对每条染色体内 PR 显著窗口布尔向量做 circular permutation，
    再构建置换后的 peak 表
    """
    permuted_peak_records = []

    for chrom, chr_df in merged_df.groupby("CHR", sort=False):
        sig = chr_df[pr_sig_col].to_numpy(dtype=bool)
        n = len(sig)

        if n == 0:
            continue
        elif n == 1:
            shifted = sig.copy()
        else:
            shift = rng.integers(0, n)
            shifted = circular_shift_boolean(sig, shift)

        peaks = mask_to_peaks(chr_df, shifted)
        permuted_peak_records.extend(peaks)

    perm_peak_df = pd.DataFrame(permuted_peak_records, columns=["CHR", "PEAK_START", "PEAK_END"])
    return perm_peak_df


def summarize_peak_lengths(peak_df):
    if peak_df.empty:
        return {
            "n_peaks": 0,
            "median_peak_size": 0,
            "mean_peak_size": 0
        }

    lengths = peak_df["PEAK_END"] - peak_df["PEAK_START"] + 1
    return {
        "n_peaks": len(peak_df),
        "median_peak_size": float(np.median(lengths)),
        "mean_peak_size": float(np.mean(lengths))
    }


def main():
    # =========================
    # 1. 输入文件与参数
    # =========================
    file_xj = "../Target_XJtumbler38_fele20otherTumbler6.lsbl.W10s5.filterSNP40"
    file_pr = "../Target_ParlorRoller22_fele20otherTumbler6.lsbl.W10s5.filterSNP40"

    z_thresh = 4.0
    n_permutations = 100000
    random_seed = 20260323

    # 输出前缀
    out_prefix = "LSBL_peak_circular_permutation"

    # =========================
    # 2. 读取并对齐数据
    # =========================
    df_xj = load_and_prepare_lsbl(file_xj)
    df_pr = load_and_prepare_lsbl(file_pr)

    merged = align_datasets(df_xj, df_pr, name1="xj", name2="pr")
    print(f"\n📌 对齐后的总窗口数: {len(merged)}")

    merged["xj_sig"] = merged["z_xj"] >= z_thresh
    merged["pr_sig"] = merged["z_pr"] >= z_thresh

    print(f"   - XJ 显著窗口数 (Z >= {z_thresh}): {merged['xj_sig'].sum()}")
    print(f"   - PR 显著窗口数 (Z >= {z_thresh}): {merged['pr_sig'].sum()}")

    # =========================
    # 3. 真实观察值：peak-level overlap
    # =========================
    xj_peak_df = build_peak_table(merged, "xj_sig")
    pr_peak_df = build_peak_table(merged, "pr_sig")

    observed_overlap = count_overlapping_peaks(xj_peak_df, pr_peak_df)

    xj_peak_summary = summarize_peak_lengths(xj_peak_df)
    pr_peak_summary = summarize_peak_lengths(pr_peak_df)

    print("\n📊 Peak 统计:")
    print(f"   - XJ peak 数: {xj_peak_summary['n_peaks']}")
    print(f"   - PR peak 数: {pr_peak_summary['n_peaks']}")
    print(f"   - XJ peak 中位长度: {xj_peak_summary['median_peak_size']:.1f} bp")
    print(f"   - PR peak 中位长度: {pr_peak_summary['median_peak_size']:.1f} bp")
    print(f"   🌟 真实观察到的重叠 peak 数（XJ peaks overlapped by PR peaks）: {observed_overlap}")

    # 保存观察到的 peaks
    xj_peak_df.to_csv(f"{out_prefix}.observed_XJ_peaks.tsv", sep="\t", index=False)
    pr_peak_df.to_csv(f"{out_prefix}.observed_PR_peaks.tsv", sep="\t", index=False)

    # =========================
    # 4. Circular permutation
    # =========================
    rng = np.random.default_rng(random_seed)
    null_overlaps = np.empty(n_permutations, dtype=int)

    print(f"\n🎲 正在执行 {n_permutations:,} 次按染色体 circular permutation ...")
    for i in range(n_permutations):
        perm_pr_peak_df = permuted_pr_peak_table_by_chr(merged, "pr_sig", rng)
        null_overlaps[i] = count_overlapping_peaks(xj_peak_df, perm_pr_peak_df)

        if (i + 1) % 10000 == 0:
            print(f"   已完成 {i + 1:,} / {n_permutations:,}")

    # =========================
    # 5. 经验 P 值
    # =========================
    p_value = (np.sum(null_overlaps >= observed_overlap) + 1) / (n_permutations + 1)

    null_mean = float(np.mean(null_overlaps))
    null_std = float(np.std(null_overlaps, ddof=1))
    null_median = float(np.median(null_overlaps))

    print("\n✅ 分析完成！")
    print(f"   - Observed overlap: {observed_overlap}")
    print(f"   - Null mean: {null_mean:.4f}")
    print(f"   - Null median: {null_median:.4f}")
    print(f"   - Null SD: {null_std:.4f}")
    print(f"   - Empirical P-value: {p_value:.6g}")

    # =========================
    # 6. 保存统计摘要
    # =========================
    summary = pd.DataFrame({
        "metric": [
            "aligned_windows",
            "z_threshold",
            "n_permutations",
            "xj_sig_windows",
            "pr_sig_windows",
            "xj_peaks",
            "pr_peaks",
            "observed_overlap_peaks",
            "null_mean_overlap",
            "null_median_overlap",
            "null_sd_overlap",
            "empirical_p_value"
        ],
        "value": [
            len(merged),
            z_thresh,
            n_permutations,
            int(merged["xj_sig"].sum()),
            int(merged["pr_sig"].sum()),
            len(xj_peak_df),
            len(pr_peak_df),
            int(observed_overlap),
            null_mean,
            null_median,
            null_std,
            p_value
        ]
    })
    summary.to_csv(f"{out_prefix}.summary.tsv", sep="\t", index=False)

    pd.DataFrame({"null_overlap_peaks": null_overlaps}).to_csv(
        f"{out_prefix}.null_distribution.tsv",
        sep="\t",
        index=False
    )

    # =========================
    # 7. 绘图
    # =========================
    plt.figure(figsize=(8, 6))

    min_val = int(null_overlaps.min())
    max_val = int(null_overlaps.max())
    bins = np.arange(min_val - 0.5, max_val + 1.5, 1)

    plt.hist(null_overlaps, bins=bins, alpha=0.75)
    plt.axvline(
        observed_overlap,
        linestyle="--",
        linewidth=2,
        label=f"Observed overlap = {observed_overlap}"
    )

    plt.title(
        "Circular permutation test for LSBL peak overlap\n"
        f"Empirical P = {p_value:.6g}",
        fontsize=13,
        fontweight="bold"
    )
    plt.xlabel("Number of overlapping XJ peaks")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()

    fig_png = f"{out_prefix}.png"
    fig_pdf = f"{out_prefix}.pdf"
    plt.savefig(fig_png, dpi=300)
    plt.savefig(fig_pdf)
    plt.close()

    print(f"\n📁 输出文件:")
    print(f"   - {out_prefix}.observed_XJ_peaks.tsv")
    print(f"   - {out_prefix}.observed_PR_peaks.tsv")
    print(f"   - {out_prefix}.null_distribution.tsv")
    print(f"   - {out_prefix}.summary.tsv")
    print(f"   - {fig_png}")
    print(f"   - {fig_pdf}")

    if p_value < 0.05:
        print("\n🧬 结论：两群体的 LSBL peak overlap 显著高于按染色体 circular permutation 的随机期望，")
        print("      支持存在非随机共享选择信号（consistent with convergent/shared selection signals）。")
    else:
        print("\n🧬 结论：观察到的 LSBL peak overlap 未显著高于按染色体 circular permutation 的随机期望。")


if __name__ == "__main__":
    main()
