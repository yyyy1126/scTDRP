#!/usr/bin/env python3
"""
B-ALL 肿瘤微环境 scTDRP 分析

分析B-ALL样本中的non-B细胞（T细胞、髓系细胞）与正常对应细胞类型的偏离
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

BLOOD_MAP_PATH = "../../1.data/processed/scRNAseq/merged/blood_map.h5ad"
BALL_PREPROCESSED_PATH = "../../1.data/raw/scRNAseq/Witkowski2019_BALL/GSE130116_qc_filtered_loose.h5ad"
OUT_DIR = "./results_ball"
FIG_DIR = "./figures_ball"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

OT_METRIC = "sqeuclidean"
OT_REG = 0.05
USE_REP = "X_pca"
N_TOP_GENES = 2000
N_PCS = 50

# ========================== 辅助函数 ==========================
def prepare_normal_reference(blood_map, stage_key, stage_order, terminal_stage):
    """从blood_map提取正常参考并做基因映射"""
    adata = blood_map[blood_map.obs[stage_key].isin(stage_order)].copy()
    if 'gene_symbols' in adata.var.columns:
        adata.var['gene_symbols'] = adata.var['gene_symbols'].astype(str)
        dup_mask = adata.var['gene_symbols'].duplicated(keep=False)
        if dup_mask.sum() > 0:
            mean_expr = np.array(adata.X.mean(axis=0)).flatten()
            if hasattr(mean_expr, 'toarray'):
                mean_expr = mean_expr.toarray().flatten()
            adata.var['_mean_expr'] = mean_expr
            keep_idx = []
            for sym, group in adata.var.groupby('gene_symbols'):
                if len(group) > 1:
                    keep_idx.append(group['_mean_expr'].idxmax())
                else:
                    keep_idx.append(group.index[0])
            adata = adata[:, keep_idx].copy()
        adata.var_names = adata.var['gene_symbols'].values
        adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

# ========================== 1. T细胞微环境分析 ==========================
print("="*60)
print("Part 1: T-cell Microenvironment Analysis")
print("="*60)

blood_map = sc.read_h5ad(BLOOD_MAP_PATH)

# T细胞正常参考
t_stages = ['CD4 Naive', 'CD4 Central Memory', 'CD4 Effector Memory',
            'CD8 Naive', 'CD8 Central Memory', 'CD8 Effector Memory 1', 'CD8 Effector Memory 2']
adata_t_normal = prepare_normal_reference(blood_map, 'AuthorCellType', t_stages, 'CD4 Effector Memory')
print(f"  Normal T cells: {adata_t_normal.shape}")
print(adata_t_normal.obs['AuthorCellType'].value_counts().to_string())

# B-ALL T细胞
adata_ball = sc.read_h5ad(BALL_PREPROCESSED_PATH)
adata_ball = adata_ball[adata_ball.obs['condition'] == 'B-ALL_Diagnosis'].copy()
sc.pp.normalize_total(adata_ball, target_sum=1e4)
sc.pp.log1p(adata_ball)

b_markers = ['CD19', 'MS4A1', 'CD79A', 'CD79B', 'PAX5', 'VPREB1', 'RAG1', 'RAG2']
t_markers = ['CD3D', 'CD3E', 'CD4', 'CD8A', 'TRAC']
sc.tl.score_genes(adata_ball, gene_list=[g for g in b_markers if g in adata_ball.var_names], score_name='b_cell_score')
sc.tl.score_genes(adata_ball, gene_list=[g for g in t_markers if g in adata_ball.var_names], score_name='t_cell_score')

b_mask = adata_ball.obs['b_cell_score'] > adata_ball.obs['b_cell_score'].median()
t_mask = (~b_mask) & (adata_ball.obs['t_cell_score'] > adata_ball.obs['t_cell_score'].median())
adata_t_disease = adata_ball[t_mask].copy()
print(f"\n  B-ALL T cells: {adata_t_disease.shape}")

# 联合PCA
common_genes = list(set(adata_t_normal.var_names) & set(adata_t_disease.var_names))
adata_t_normal = adata_t_normal[:, common_genes].copy()
adata_t_disease = adata_t_disease[:, common_genes].copy()

adata_combined = sc.concat([adata_t_normal, adata_t_disease], label="_source",
                           keys=["normal", "disease"], index_unique="-", join="outer")
sc.pp.highly_variable_genes(adata_combined, n_top_genes=N_TOP_GENES, flavor="seurat_v3")
sc.pp.scale(adata_combined, max_value=10)
sc.pp.pca(adata_combined, n_comps=N_PCS, use_highly_variable=True)

normal_mask = adata_combined.obs["_source"] == "normal"
adata_t_normal = adata_combined[normal_mask].copy()
adata_t_disease = adata_combined[~normal_mask].copy()

# scTDRP
print("\n  Running scTDRP for T cells...")
analyzer_t = TDRPAnalyzer(
    normal_adata=adata_t_normal,
    stage_key='AuthorCellType',
    terminal_stage='CD4 Effector Memory',
    stage_order=t_stages,
    use_rep=USE_REP,
    n_top_genes=N_TOP_GENES,
)
analyzer_t.prepare_data()
analyzer_t.build_ot_cost_map(metric=OT_METRIC, reg=OT_REG)
analyzer_t.compute_tdi(adata_t_disease, metric=OT_METRIC, reg=OT_REG)
tdi_t = analyzer_t.tdi_results.copy()
tdi_t['cell_type'] = 'T_cell'

print(f"\n  T-cell TDI: n={len(tdi_t)}, mean={tdi_t['TDI'].mean():.4f}, median={tdi_t['TDI'].median():.4f}")
print(f"  Best match stage: {tdi_t['Best_Match_Stage'].value_counts().to_string()}")

# ========================== 2. 髓系微环境分析 ==========================
print("\n" + "="*60)
print("Part 2: Myeloid Microenvironment Analysis")
print("="*60)

myeloid_stages = ['CD14 Mono', 'CD16 Mono', 'Early ProMono', 'Late ProMono']
adata_my_normal = prepare_normal_reference(blood_map, 'AuthorCellType', myeloid_stages, 'CD14 Mono')
print(f"  Normal myeloid: {adata_my_normal.shape}")
print(adata_my_normal.obs['AuthorCellType'].value_counts().to_string())

# B-ALL髓系细胞
myeloid_markers = ['CD14', 'LYZ', 'S100A8', 'S100A9', 'CD33']
sc.tl.score_genes(adata_ball, gene_list=[g for g in myeloid_markers if g in adata_ball.var_names], score_name='myeloid_score')

myeloid_mask = (~b_mask) & (~t_mask) & (adata_ball.obs['myeloid_score'] > adata_ball.obs['myeloid_score'].median())
adata_my_disease = adata_ball[myeloid_mask].copy()
print(f"\n  B-ALL myeloid cells: {adata_my_disease.shape}")

if adata_my_disease.shape[0] >= 10:
    common_genes = list(set(adata_my_normal.var_names) & set(adata_my_disease.var_names))
    adata_my_normal = adata_my_normal[:, common_genes].copy()
    adata_my_disease = adata_my_disease[:, common_genes].copy()
    
    adata_combined = sc.concat([adata_my_normal, adata_my_disease], label="_source",
                               keys=["normal", "disease"], index_unique="-", join="outer")
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=N_TOP_GENES, flavor="seurat_v3")
    sc.pp.scale(adata_combined, max_value=10)
    sc.pp.pca(adata_combined, n_comps=N_PCS, use_highly_variable=True)
    
    normal_mask = adata_combined.obs["_source"] == "normal"
    adata_my_normal = adata_combined[normal_mask].copy()
    adata_my_disease = adata_combined[~normal_mask].copy()
    
    print("\n  Running scTDRP for myeloid cells...")
    analyzer_my = TDRPAnalyzer(
        normal_adata=adata_my_normal,
        stage_key='AuthorCellType',
        terminal_stage='CD14 Mono',
        stage_order=myeloid_stages,
        use_rep=USE_REP,
        n_top_genes=N_TOP_GENES,
    )
    analyzer_my.prepare_data()
    analyzer_my.build_ot_cost_map(metric=OT_METRIC, reg=OT_REG)
    analyzer_my.compute_tdi(adata_my_disease, metric=OT_METRIC, reg=OT_REG)
    tdi_my = analyzer_my.tdi_results.copy()
    tdi_my['cell_type'] = 'Myeloid'
    
    print(f"\n  Myeloid TDI: n={len(tdi_my)}, mean={tdi_my['TDI'].mean():.4f}, median={tdi_my['TDI'].median():.4f}")
    print(f"  Best match stage: {tdi_my['Best_Match_Stage'].value_counts().to_string()}")
else:
    print("  Too few myeloid cells (<10), skipping myeloid analysis")
    tdi_my = pd.DataFrame()

# ========================== 3. 正常对应细胞TDI（对照）======================
print("\n" + "="*60)
print("Part 3: Normal counterpart TDI (controls)")
print("="*60)

# T细胞正常对照
print("  Computing TDI for normal T cells...")
analyzer_t.compute_tdi(adata_t_normal, metric=OT_METRIC, reg=OT_REG)
tdi_t_normal = analyzer_t.tdi_results.copy()
tdi_t_normal['cell_type'] = 'T_cell_normal'
print(f"  Normal T-cell TDI: n={len(tdi_t_normal)}, mean={tdi_t_normal['TDI'].mean():.4f}, median={tdi_t_normal['TDI'].median():.4f}")

if adata_my_disease.shape[0] >= 10:
    print("  Computing TDI for normal myeloid cells...")
    analyzer_my.compute_tdi(adata_my_normal, metric=OT_METRIC, reg=OT_REG)
    tdi_my_normal = analyzer_my.tdi_results.copy()
    tdi_my_normal['cell_type'] = 'Myeloid_normal'
    print(f"  Normal myeloid TDI: n={len(tdi_my_normal)}, mean={tdi_my_normal['TDI'].mean():.4f}, median={tdi_my_normal['TDI'].median():.4f}")
else:
    tdi_my_normal = pd.DataFrame()

# ========================== 4. 可视化 ==========================
print("\n" + "="*60)
print("Part 4: Visualization")
print("="*60)

tdi_all = pd.concat([tdi_t, tdi_t_normal, tdi_my, tdi_my_normal], ignore_index=True)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# T细胞
ax = axes[0]
tdi_t_plot = tdi_all[tdi_all['cell_type'].isin(['T_cell', 'T_cell_normal'])]
sns.boxplot(data=tdi_t_plot, x='cell_type', y='TDI', 
            palette=['#e74c3c', '#2ecc71'], order=['T_cell_normal', 'T_cell'], ax=ax)
ax.set_title("T-cell Microenvironment TDI", fontsize=13, fontweight='bold')
ax.set_ylabel("TDI")
ax.set_xticklabels(['Normal T', 'B-ALL TME T'])

# 髓系
ax = axes[1]
if len(tdi_my) > 0:
    tdi_my_plot = tdi_all[tdi_all['cell_type'].isin(['Myeloid', 'Myeloid_normal'])]
    sns.boxplot(data=tdi_my_plot, x='cell_type', y='TDI',
                palette=['#e74c3c', '#2ecc71'], order=['Myeloid_normal', 'Myeloid'], ax=ax)
    ax.set_title("Myeloid Microenvironment TDI", fontsize=13, fontweight='bold')
    ax.set_ylabel("TDI")
    ax.set_xticklabels(['Normal Myeloid', 'B-ALL TME Myeloid'])
else:
    ax.text(0.5, 0.5, 'Too few myeloid cells\n(n<10)', ha='center', va='center', transform=ax.transAxes)
    ax.set_title("Myeloid Microenvironment TDI", fontsize=13, fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Figure4S_BALL_Microenvironment.pdf"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, "Figure4S_BALL_Microenvironment.png"), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: Figure4S_BALL_Microenvironment")

# ========================== 5. 统计摘要 ==========================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)

# T细胞
print(f"\nT-cell microenvironment:")
print(f"  B-ALL TME T cells: n={len(tdi_t)}, mean TDI={tdi_t['TDI'].mean():.4f}")
print(f"  Normal T cells: n={len(tdi_t_normal)}, mean TDI={tdi_t_normal['TDI'].mean():.4f}")
_, p_t = stats.mannwhitneyu(tdi_t['TDI'], tdi_t_normal['TDI'], alternative='two-sided')
print(f"  Mann-Whitney U (B-ALL vs Normal): p={p_t:.2e}")

def cohens_d(x, y):
    nx, ny = len(x), len(y)
    dof = nx + ny - 2
    pooled_std = np.sqrt(((nx-1)*x.var() + (ny-1)*y.var()) / dof)
    return (x.mean() - y.mean()) / pooled_std

d_t = cohens_d(tdi_t['TDI'], tdi_t_normal['TDI'])
print(f"  Cohen's d: {d_t:.3f}")

# 髓系
if len(tdi_my) > 0:
    print(f"\nMyeloid microenvironment:")
    print(f"  B-ALL TME myeloid: n={len(tdi_my)}, mean TDI={tdi_my['TDI'].mean():.4f}")
    print(f"  Normal myeloid: n={len(tdi_my_normal)}, mean TDI={tdi_my_normal['TDI'].mean():.4f}")
    _, p_my = stats.mannwhitneyu(tdi_my['TDI'], tdi_my_normal['TDI'], alternative='two-sided')
    print(f"  Mann-Whitney U (B-ALL vs Normal): p={p_my:.2e}")
    d_my = cohens_d(tdi_my['TDI'], tdi_my_normal['TDI'])
    print(f"  Cohen's d: {d_my:.3f}")

# 保存
tdi_all.to_csv(os.path.join(OUT_DIR, "tdi_ball_microenvironment.csv"), index=False)

print("\n" + "="*60)
print(f"Done! Results -> {OUT_DIR}")
print(f"Figures -> {FIG_DIR}")
print("="*60)
