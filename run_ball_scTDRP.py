#!/usr/bin/env python3
"""
scTDRP B-ALL 验证分析脚本 (v3 - single-cell TDI)

数据来源: Witkowski et al. 2020, Cancer Cell (GSE130116)
正常B细胞发育轨迹从blood_map中提取:
  Pro-B VDJ -> Large Pre-B -> Immature B -> Mature B
"""

import os
import sys
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from scTDRP import TDRPAnalyzer

warnings.filterwarnings("ignore")

# ========================== 配置 ==========================
BLOOD_MAP_PATH = "../../1.data/processed/scRNAseq/merged/blood_map.h5ad"
BALL_PREPROCESSED_PATH = "../../1.data/raw/scRNAseq/Witkowski2019_BALL/GSE130116_qc_filtered_loose.h5ad"
MODULES_PATH = "./modules.json"
OUT_DIR = "./results_ball"
FIG_DIR = "./figures_ball"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

STAGE_ORDER = ["Pro-B VDJ", "Large Pre-B", "Immature B", "Mature B"]
TERMINAL_STAGE = "Mature B"
STAGE_KEY = "AuthorCellType"

OT_METRIC = "sqeuclidean"
OT_REG = 0.05
USE_REP = "X_pca"
N_TOP_GENES = 2000
N_PCS = 50

B_CELL_MARKERS = ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5', 'VPREB1', 'RAG1', 'RAG2']
T_CELL_MARKERS = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'TRAC']
MYELOID_MARKERS = ['CD14', 'LYZ', 'S100A8', 'S100A9', 'CD33']

# ========================== 步骤1: 正常B细胞轨迹 ==========================
print("="*60)
print("Step 1: Normal B-cell trajectory from blood_map")
print("="*60)

blood_map = sc.read_h5ad(BLOOD_MAP_PATH)
adata_b_normal = blood_map[blood_map.obs[STAGE_KEY].isin(STAGE_ORDER)].copy()
print(f"  B-cell stages: {adata_b_normal.shape}")
print(adata_b_normal.obs[STAGE_KEY].value_counts().to_string())

if 'gene_symbols' in adata_b_normal.var.columns:
    adata_b_normal.var['gene_symbols'] = adata_b_normal.var['gene_symbols'].astype(str)
    dup_mask = adata_b_normal.var['gene_symbols'].duplicated(keep=False)
    if dup_mask.sum() > 0:
        mean_expr = np.array(adata_b_normal.X.mean(axis=0)).flatten()
        if hasattr(mean_expr, 'toarray'):
            mean_expr = mean_expr.toarray().flatten()
        adata_b_normal.var['_mean_expr'] = mean_expr
        keep_idx = []
        for sym, group in adata_b_normal.var.groupby('gene_symbols'):
            if len(group) > 1:
                keep_idx.append(group['_mean_expr'].idxmax())
            else:
                keep_idx.append(group.index[0])
        adata_b_normal = adata_b_normal[:, keep_idx].copy()
    adata_b_normal.var_names = adata_b_normal.var['gene_symbols'].values
    adata_b_normal.var_names_make_unique()

print(f"  After gene mapping: {adata_b_normal.shape}")

# ========================== 步骤2: B-ALL数据 ==========================
print("\n" + "="*60)
print("Step 2: Loading QC-filtered B-ALL data")
print("="*60)

adata_ball = sc.read_h5ad(BALL_PREPROCESSED_PATH)
print(f"  QC-filtered: {adata_ball.shape}")

# 只保留诊断样本
adata_ball = adata_ball[adata_ball.obs['condition'] == 'B-ALL_Diagnosis'].copy()
print(f"  Diagnosis only: {adata_ball.shape}")

# ========================== 步骤3: B细胞筛选 ==========================
print("\n" + "="*60)
print("Step 3: B-cell selection")
print("="*60)

sc.pp.normalize_total(adata_ball, target_sum=1e4)
sc.pp.log1p(adata_ball)

available_b = [g for g in B_CELL_MARKERS if g in adata_ball.var_names]
available_t = [g for g in T_CELL_MARKERS if g in adata_ball.var_names]
available_m = [g for g in MYELOID_MARKERS if g in adata_ball.var_names]

if available_b:
    sc.tl.score_genes(adata_ball, gene_list=available_b, score_name='b_cell_score')
if available_t:
    sc.tl.score_genes(adata_ball, gene_list=available_t, score_name='t_cell_score')
if available_m:
    sc.tl.score_genes(adata_ball, gene_list=available_m, score_name='myeloid_score')

# 宽松筛选：B score > median 即可
b_mask = adata_ball.obs['b_cell_score'] > adata_ball.obs['b_cell_score'].median()
adata_ball_b = adata_ball[b_mask].copy()
print(f"  B-cell selected: {adata_ball_b.shape[0]} / {adata_ball.shape[0]}")
print(adata_ball_b.obs['sample_id'].value_counts().to_string())

# 正常B细胞标准化
sc.pp.normalize_total(adata_b_normal, target_sum=1e4)
sc.pp.log1p(adata_b_normal)

# ========================== 步骤4: 联合PCA ==========================
print("\n" + "="*60)
print("Step 4: Preprocessing")
print("="*60)

common_genes = list(set(adata_b_normal.var_names) & set(adata_ball_b.var_names))
print(f"  Common genes: {len(common_genes)}")

adata_b_normal = adata_b_normal[:, common_genes].copy()
adata_ball_b = adata_ball_b[:, common_genes].copy()

adata_combined = sc.concat(
    [adata_b_normal, adata_ball_b],
    label="_source",
    keys=["normal", "disease"],
    index_unique="-",
    join="outer"
)
sc.pp.highly_variable_genes(adata_combined, n_top_genes=N_TOP_GENES, flavor="seurat_v3")
sc.pp.scale(adata_combined, max_value=10)
sc.pp.pca(adata_combined, n_comps=N_PCS, use_highly_variable=True)

normal_mask = adata_combined.obs["_source"] == "normal"
adata_b_normal = adata_combined[normal_mask].copy()
adata_ball_b = adata_combined[~normal_mask].copy()

print(f"  Normal B cells: {adata_b_normal.shape}")
print(f"  B-ALL B cells: {adata_ball_b.shape}")

# ========================== 步骤5: scTDRP (single-cell TDI) ==========================
print("\n" + "="*60)
print("Step 5: Running scTDRP (single-cell level)")
print("="*60)

analyzer = TDRPAnalyzer(
    normal_adata=adata_b_normal,
    stage_key=STAGE_KEY,
    terminal_stage=TERMINAL_STAGE,
    stage_order=STAGE_ORDER,
    use_rep=USE_REP,
    n_top_genes=N_TOP_GENES,
)
analyzer.prepare_data()
analyzer.build_ot_cost_map(metric=OT_METRIC, reg=OT_REG)

# B-ALL单个细胞TDI
print("\n  Computing TDI for B-ALL cells...")
analyzer.compute_tdi(adata_ball_b, metric=OT_METRIC, reg=OT_REG)
tdi_ball = analyzer.tdi_results.copy()
tdi_ball['condition'] = 'B-ALL'
tdi_ball['true_stage'] = 'B-ALL'
tdi_ball['sample_id'] = adata_ball_b.obs['sample_id'].values

# 正常单个细胞TDI
print("\n  Computing TDI for Normal B cells...")
analyzer.compute_tdi(adata_b_normal, metric=OT_METRIC, reg=OT_REG)
tdi_normal = analyzer.tdi_results.copy()
tdi_normal['condition'] = 'Normal'
tdi_normal['true_stage'] = adata_b_normal.obs[STAGE_KEY].values

# 合并
tdi_all = pd.concat([tdi_normal, tdi_ball], ignore_index=True)

# ========================== 步骤6: 可视化 ==========================
print("\n" + "="*60)
print("Step 6: Visualization")
print("="*60)

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Panel A: TDI by condition
ax = axes[0, 0]
sns.boxplot(data=tdi_all, x='condition', y='TDI', palette=['#2ecc71', '#e74c3c'], ax=ax)
ax.set_title("TDI: Normal B cells vs B-ALL", fontsize=13, fontweight='bold')
ax.set_ylabel("TDI")

# Panel B: B-ALL Best Match Stage
ax = axes[0, 1]
ball_stages = tdi_ball['Best_Match_Stage'].value_counts(normalize=True)
ball_stages = ball_stages.reindex(STAGE_ORDER, fill_value=0)
bars = ax.bar(ball_stages.index, ball_stages.values, color='#e74c3c')
ax.set_title("B-ALL Best Match Stage", fontsize=13, fontweight='bold')
ax.set_ylabel("Proportion")
ax.tick_params(axis='x', rotation=45)
# Add text labels
for bar, val in zip(bars, ball_stages.values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, 
            f"{val*100:.1f}%", ha='center', va='bottom', fontsize=10)

# Panel C: Normal B cells TDI by true stage
ax = axes[1, 0]
sns.boxplot(data=tdi_normal, x='true_stage', y='TDI', order=STAGE_ORDER, palette='Greens', ax=ax)
ax.set_title("Normal B cells TDI by Stage", fontsize=13, fontweight='bold')
ax.set_ylabel("TDI")
ax.tick_params(axis='x', rotation=45)

# Panel D: B-ALL vs matched Normal Large Pre-B
ax = axes[1, 1]
matched_normal = tdi_normal[tdi_normal['true_stage'] == 'Large Pre-B']['TDI']
ball_tdi = tdi_ball['TDI']
sns.kdeplot(matched_normal, color='#2ecc71', fill=True, label=f'Normal Large Pre-B (n={len(matched_normal)})', ax=ax)
sns.kdeplot(ball_tdi, color='#e74c3c', fill=True, label=f'B-ALL (n={len(ball_tdi)})', ax=ax)
ax.set_title("TDI Distribution: B-ALL vs Normal Large Pre-B", fontsize=13, fontweight='bold')
ax.set_xlabel("TDI")
ax.set_ylabel("Density")
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Figure4_BALL_scTDRP.pdf"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, "Figure4_BALL_scTDRP.png"), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: Figure4_BALL_scTDRP")

# ========================== 步骤7: 统计与摘要 ==========================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)

print(f"\n1. Normal B cells: {len(tdi_normal)}")
for stage in STAGE_ORDER:
    n = (tdi_normal['true_stage'] == stage).sum()
    mean_tdi = tdi_normal[tdi_normal['true_stage'] == stage]['TDI'].mean()
    print(f"   {stage}: n={n}, mean TDI={mean_tdi:.4f}")

print(f"\n2. B-ALL B cells: n={len(tdi_ball)}")
print(f"   Mean TDI: {tdi_ball['TDI'].mean():.4f}")
print(f"   Median TDI: {tdi_ball['TDI'].median():.4f}")

print(f"\n3. B-ALL Best Match Stage:")
for stage, prop in ball_stages.items():
    print(f"   {stage}: {prop*100:.1f}%")

# 统计检验
ball_tdi = tdi_ball['TDI']
norm_tdi = tdi_normal['TDI']
_, p_all = stats.mannwhitneyu(ball_tdi, norm_tdi, alternative='two-sided')
print(f"\n4. Mann-Whitney U (B-ALL vs All Normal): p={p_all:.2e}")

matched_normal = tdi_normal[tdi_normal['true_stage'] == 'Large Pre-B']['TDI']
_, p_matched = stats.mannwhitneyu(ball_tdi, matched_normal, alternative='two-sided')
print(f"   Mann-Whitney U (B-ALL vs Normal Large Pre-B): p={p_matched:.2e}")

# Effect size (Cohen's d)
def cohens_d(x, y):
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx-1)*x.var() + (ny-1)*y.var()) / dof)
    return (x.mean() - y.mean()) / pooled_std

d_all = cohens_d(ball_tdi, norm_tdi)
d_matched = cohens_d(ball_tdi, matched_normal)
print(f"   Cohen's d (B-ALL vs All Normal): {d_all:.3f}")
print(f"   Cohen's d (B-ALL vs Large Pre-B): {d_matched:.3f}")

# 保存结果
tdi_all.to_csv(os.path.join(OUT_DIR, "tdi_ball_singlecell.csv"), index=False)
stage_props = ball_stages.to_frame(name='proportion').reset_index()
stage_props.to_csv(os.path.join(OUT_DIR, "ball_best_match_stage.csv"), index=False)

print("\n" + "="*60)
print(f"Done! Results -> {OUT_DIR}")
print(f"Figures -> {FIG_DIR}")
print("="*60)
