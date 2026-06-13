#!/usr/bin/env python3
"""
B-ALL 遗传亚型比较分析：ETV-RUNX1 vs Ph+

数据来源: Witkowski et al. 2019, Cancer Cell (GSE130116)
"""

import os
import sys
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

# ========================== 步骤1: 正常B细胞轨迹 ==========================
print("="*60)
print("Step 1: Normal B-cell trajectory from blood_map")
print("="*60)

blood_map = sc.read_h5ad(BLOOD_MAP_PATH)
adata_b_normal = blood_map[blood_map.obs[STAGE_KEY].isin(STAGE_ORDER)].copy()
print(f"  B-cell stages: {adata_b_normal.shape}")

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

# ========================== 步骤2: B-ALL恶性B细胞按亚型分组 ==========================
print("\n" + "="*60)
print("Step 2: Loading B-ALL data and subtype classification")
print("="*60)

adata_ball = sc.read_h5ad(BALL_PREPROCESSED_PATH)
adata_ball = adata_ball[adata_ball.obs['condition'] == 'B-ALL_Diagnosis'].copy()
print(f"  Diagnosis cells: {adata_ball.shape}")

sc.pp.normalize_total(adata_ball, target_sum=1e4)
sc.pp.log1p(adata_ball)

# B-cell selection
b_markers = ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5', 'VPREB1', 'RAG1', 'RAG2']
avail_b = [g for g in b_markers if g in adata_ball.var_names]
sc.tl.score_genes(adata_ball, gene_list=avail_b, score_name='b_cell_score')

b_mask = adata_ball.obs['b_cell_score'] > adata_ball.obs['b_cell_score'].median()
adata_ball_b = adata_ball[b_mask].copy()
print(f"  B-cell selected: {adata_ball_b.shape[0]}")

# 亚型分组
etv_samples = [s for s in adata_ball_b.obs['sample_id'].unique() if 'ETV' in s]
ph_samples = [s for s in adata_ball_b.obs['sample_id'].unique() if 'PH' in s]

adata_ball_b.obs['subtype'] = 'Unknown'
adata_ball_b.obs.loc[adata_ball_b.obs['sample_id'].isin(etv_samples), 'subtype'] = 'ETV-RUNX1'
adata_ball_b.obs.loc[adata_ball_b.obs['sample_id'].isin(ph_samples), 'subtype'] = 'Ph+'

print(f"\n  ETV-RUNX1 samples: {etv_samples}")
print(f"  Ph+ samples: {ph_samples}")
print(f"\n  Subtype distribution:")
print(adata_ball_b.obs['subtype'].value_counts().to_string())
print(f"\n  By sample:")
print(adata_ball_b.obs.groupby('sample_id').size().to_string())

# ========================== 步骤3: 联合PCA ==========================
print("\n" + "="*60)
print("Step 3: Preprocessing")
print("="*60)

sc.pp.normalize_total(adata_b_normal, target_sum=1e4)
sc.pp.log1p(adata_b_normal)

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

# ========================== 步骤4: scTDRP ==========================
print("\n" + "="*60)
print("Step 4: Running scTDRP")
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

# 全量B-ALL单个细胞TDI
print("\n  Computing TDI for all B-ALL cells...")
analyzer.compute_tdi(adata_ball_b, metric=OT_METRIC, reg=OT_REG)
tdi_all = analyzer.tdi_results.copy()
tdi_all['subtype'] = adata_ball_b.obs['subtype'].values
tdi_all['sample_id'] = adata_ball_b.obs['sample_id'].values

# ========================== 步骤5: 亚型比较可视化 ==========================
print("\n" + "="*60)
print("Step 5: Visualization")
print("="*60)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel A: TDI by subtype
ax = axes[0]
sns.boxplot(data=tdi_all, x='subtype', y='TDI', palette=['#3498db', '#e74c3c'], ax=ax)
ax.set_title("TDI by Genetic Subtype", fontsize=13, fontweight='bold')
ax.set_ylabel("TDI")

# Panel B: TDI density by subtype
ax = axes[1]
for subtype, color in [('ETV-RUNX1', '#3498db'), ('Ph+', '#e74c3c')]:
    subset = tdi_all[tdi_all['subtype'] == subtype]['TDI']
    sns.kdeplot(subset, color=color, fill=True, label=f'{subtype} (n={len(subset)})', ax=ax)
ax.set_title("TDI Distribution by Subtype", fontsize=13, fontweight='bold')
ax.set_xlabel("TDI")
ax.set_ylabel("Density")
ax.legend()

# Panel C: Best match stage by subtype
ax = axes[2]
stage_props = pd.crosstab(tdi_all['subtype'], tdi_all['Best_Match_Stage'], normalize='index')
stage_props = stage_props.reindex(columns=STAGE_ORDER, fill_value=0)
stage_props.plot(kind='bar', ax=ax, color=['#2ecc71', '#f39c12', '#9b59b6', '#e74c3c'])
ax.set_title("Best Match Stage by Subtype", fontsize=13, fontweight='bold')
ax.set_ylabel("Proportion")
ax.tick_params(axis='x', rotation=45)
ax.legend(title='Stage', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Figure4S_BALL_Subtype_Comparison.pdf"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, "Figure4S_BALL_Subtype_Comparison.png"), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: Figure4S_BALL_Subtype_Comparison")

# ========================== 步骤6: 统计摘要 ==========================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)

for subtype in ['ETV-RUNX1', 'Ph+']:
    subset = tdi_all[tdi_all['subtype'] == subtype]
    print(f"\n{subtype}: n={len(subset)}")
    print(f"  Mean TDI: {subset['TDI'].mean():.4f}")
    print(f"  Median TDI: {subset['TDI'].median():.4f}")
    print(f"  Best match stage:")
    for stage, prop in subset['Best_Match_Stage'].value_counts(normalize=True).items():
        print(f"    {stage}: {prop*100:.1f}%")

# 统计检验
etv_tdi = tdi_all[tdi_all['subtype'] == 'ETV-RUNX1']['TDI']
ph_tdi = tdi_all[tdi_all['subtype'] == 'Ph+']['TDI']
_, p_val = stats.mannwhitneyu(etv_tdi, ph_tdi, alternative='two-sided')
print(f"\nMann-Whitney U (ETV vs Ph+): p={p_val:.2e}")

# Effect size
def cohens_d(x, y):
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx-1)*x.var() + (ny-1)*y.var()) / dof)
    return (x.mean() - y.mean()) / pooled_std

d = cohens_d(ph_tdi, etv_tdi)
print(f"Cohen's d (Ph+ vs ETV): {d:.3f}")

# 保存结果
tdi_all.to_csv(os.path.join(OUT_DIR, "tdi_ball_subtype.csv"), index=False)

print("\n" + "="*60)
print(f"Done! Results -> {OUT_DIR}")
print(f"Figures -> {FIG_DIR}")
print("="*60)
