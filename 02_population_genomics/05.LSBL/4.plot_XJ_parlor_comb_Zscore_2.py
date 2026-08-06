import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
import io
import sys

# --- 1. 环境与样式设置 ---
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

def load_lsbl(file_path):
    print(f"🚀 正在加载: {file_path}")
    try:
        df = pd.read_csv(file_path, sep=r'\s+', header=None, comment='#', low_memory=False)
        cols = ['scaffold', 'start', 'end', 'sum_lsbl', 'count', 'lsbl']
        df.columns = cols[:df.shape[1]]
        for col in ['start', 'end', 'lsbl']:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        return df.dropna(subset=['lsbl'])
    except Exception as e:
        print(f"❌ 无法读取 {file_path}: {e}")
        return pd.DataFrame()

# --- 2. 加载数据 ---
file_xj = 'Target_XJtumbler38_fele20otherTumbler6.lsbl.W10s5.filterSNP40'
file_pr = 'Target_ParlorRoller22_fele20otherTumbler6.lsbl.W10s5.filterSNP40'
file_comb = 'Target_combineXJtumblerParlorRoller60.lsbl.W10s5.filterSNP40'

df_xj = load_lsbl(file_xj)
df_pr = load_lsbl(file_pr)
df_comb = load_lsbl(file_comb)

if df_xj.empty or df_pr.empty or df_comb.empty:
    sys.exit("❌ 数据文件缺失，请检查当前路径。")

# --- 3. 合并与 Z-score 计算 ---
merged = pd.merge(df_xj, df_pr, on=['scaffold', 'start', 'end'], suffixes=('_xj', '_pr'))
merged = pd.merge(merged, df_comb, on=['scaffold', 'start', 'end'])
merged.rename(columns={'lsbl': 'lsbl_comb'}, inplace=True)

for suffix in ['xj', 'pr', 'comb']:
    col = f'lsbl_{suffix}'
    merged[f'z_{suffix}'] = (merged[col] - merged[col].mean()) / merged[col].std()

Z_THRESHOLD = 4

# --- 3.5 提取并保存 Combine 组显著位点 (Z >= 4) ---
# 筛选单文件中 Z-score >= 4 的所有位点
combine_sig_all = merged[merged['z_comb'] >= Z_THRESHOLD].copy()
# 保存为 BED 格式 (Scaffold, Start, End, Z-score) 方便 IGV 查看
combine_sig_bed = combine_sig_all[['scaffold', 'start', 'end', 'z_comb']]
combine_sig_bed.to_csv('Combine60_Significant_Z4.bed', sep='\t', header=False, index=False)

print(f"📌 单独提取 Combine 组显著位点完成:")
print(f"   - 发现 Z_comb >= {Z_THRESHOLD} 的位点共 {len(combine_sig_all)} 个")

# --- 4. 基因列表处理 ---
raw_gene_txt = """
scaffold_1 144607553 144626551 4.466 6.33 144607554 144643305 SYT10
scaffold_1 144631551 144643305 4.301 4.314 144607554 144643305 SYT10
scaffold_1 153591551 153611551 5.478 5.765 153510088 153637910 TRHDE
scaffold_15 13248228 13257190 4.628 4.821 13248229 13292253 FAT2
scaffold_15 13237190 13241017 4.628 4.821 13223659 13241017 SLC36A1
scaffold_15 13272190 13292253 4.731 4.95 13248229 13292253 FAT2
scaffold_24 8084773 8085072 4.83 5.403 8084774 8101689 PJA2
scaffold_6 39091476 39121476 5.038 4.792 39072714 39125740 ERCC6L2
scaffold_6 40740210 40751476 4.159 5.118 40740211 40831232 NTRK2
scaffold_7 19651822 19663329 5.246 4.17 19651823 19712065 LOC102092085
scaffold_7 19651822 19663329 5.246 4.17 19651823 19712065 SCN3A
scaffold_9 17091397 17097528 6.511 4.855 17091398 17097528 TROVE2
scaffold_9 17105095 17106413 6.511 4.855 17105096 17106413 GLRX2
"""
df_genes = pd.read_csv(io.StringIO(raw_gene_txt.strip()), sep=r'\s+', header=None)
df_genes.columns = ['scaffold', 's', 'e', 'z1', 'z2', 'gs', 'ge', 'gene_name']

# --- 5. 绘图 ---
fig, ax = plt.subplots(figsize=(11, 8.5))

# A. 背景点 (Z_comb < 4)
non_sig = merged[merged['z_comb'] < Z_THRESHOLD]
ax.scatter(non_sig['z_xj'], non_sig['z_pr'], 
           c=non_sig['z_comb'], cmap='Blues', 
           s=25, alpha=0.3, zorder=1, vmin=0, vmax=Z_THRESHOLD)

# B. 显著点 (Z_comb >= 4) - 全量着色
sig_comb = merged[merged['z_comb'] >= Z_THRESHOLD]
vmax = max(sig_comb['z_comb'].max(), 6)
sc = ax.scatter(sig_comb['z_xj'], sig_comb['z_pr'], 
                c=sig_comb['z_comb'], cmap='YlOrRd', 
                s=45, alpha=0.9, zorder=3, vmin=Z_THRESHOLD, vmax=vmax)

# C. 核心位点判定 (三方显著)
triple_sig = merged[(merged['z_xj'] >= Z_THRESHOLD) & 
                    (merged['z_pr'] >= Z_THRESHOLD) & 
                    (merged['z_comb'] >= Z_THRESHOLD)].copy()

texts = []
for idx, row in triple_sig.iterrows():
    # 匹配基因
    mask = (df_genes['scaffold'] == row['scaffold']) & \
           (df_genes['s'] <= row['end'] + 5000) & \
           (df_genes['e'] >= row['start'] - 5000)
    matched = sorted(df_genes[mask]['gene_name'].unique())
    
    # 只有当匹配到基因时才绘制黑圈并添加标注
    if matched:
        label = ",".join(matched)
        # 绘制黑圈
        ax.scatter(row['z_xj'], row['z_pr'], s=65, facecolors='none', 
                   edgecolors='black', linewidths=1.0, zorder=5)
        # 准备标注
        t = ax.text(row['z_xj'], row['z_pr'], label, fontsize=8, fontweight='bold', style='italic',
                    bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.8, ec='black', lw=0.5))
        texts.append(t)

# 自动文字避让
if texts:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=0.5, alpha=0.5))

# --- 6. 细节修饰 ---
# 坐标轴从 0 开始
ax.set_xlim(0, max(merged['z_xj'].max() + 0.5, 10))
ax.set_ylim(0, max(merged['z_pr'].max() + 0.5, 8))

# 辅助线
ax.axvline(Z_THRESHOLD, color='gray', linestyle='--', alpha=0.6)
ax.axhline(Z_THRESHOLD, color='gray', linestyle='--', alpha=0.6)

# 颜色条
cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
cbar.set_label(f'Z-score Combine (Significance ≥ {Z_THRESHOLD})', fontweight='bold')

ax.set_xlabel('Z-score (Xinjiang Tumbler)', fontsize=12)
ax.set_ylabel('Z-score (Parlor Roller)', fontsize=12)
ax.set_title('Screened Selection Signals (Black Circle: Annotated Triple-Significant Only)', fontsize=13)

plt.tight_layout()
output_base = 'LSBL_Triple_Comaprison_Result2'
plt.savefig(f'{output_base}.png', dpi=300)
plt.savefig(f'{output_base}.pdf')
print("✅ 绘图完成！图片已保存。")
