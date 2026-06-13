#!/usr/bin/env python3
"""
AML5 (AEL) 内部深挖分析

1. 恶性 vs 残留正常 vs 其他细胞的TDI比较
2. 不同细胞周期阶段的TDI
3. CNV score亚克隆的TDI差异
4. Seurat cluster与TDI的关系
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from scTDRP import TDRPAnalyzer
from scTDRP.distance import compute_wasserstein_distance

# ========================== 配置 ==========================
NORMAL_PATH = "../../1.data/processed/erythroid_lineage_from_MEP.h5ad"
DISEASE_PATH = "../../4.results/4.infercnv/infercnv_aml5/aml5_annotated_with_cnv_cycle.h5ad"
MODULES_PATH = "./modules.json"
OUT_DIR = "./results_aml5_deep"
FIG_DIR = "./figures_aml5_deep"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

STAGE_ORDER = ["MEP","BFU-E","CFU-E","Pro-Erythroblast","Basophilic Erythroblast","Polychromatic Erythroblast","Orthochromatic Erythroblast"]
TERMINAL_STAGE = "Orthochromatic Erythroblast"

OT_METRIC = "sqeuclidean"
OT_REG = 0.05
USE_REP = "X_pca"
N_TOP_GENES = 2000

METACELL_RESOLUTION = 1.0

# ========================== 加载数据 ==========================
print("="*60)
print("Loading data")
print("="*60)

adata_normal = sc.read_h5ad(NORMAL_PATH)
adata_aml5 = sc.read_h5ad(DISEASE_PATH)

print(f"Normal: {adata_normal.shape}")
print(f"AML5: {adata_aml5.shape}")

# ========================== 基因名统一 ==========================
if 'gene_symbols' in adata_normal.var.columns:
    adata_normal.var['gene_symbols'] = adata_normal.var['gene_symbols'].astype(str)
    dup_mask = adata_normal.var['gene_symbols'].duplicated(keep=False)
    if dup_mask.sum() > 0:
        mean_expr = np.array(adata_normal.X.mean(axis=0)).flatten()
        if hasattr(mean_expr, 'toarray'):
            mean_expr = mean_expr.toarray().flatten()
        adata_normal.var['_mean_expr'] = mean_expr
        keep_idx = []
        for sym, group in adata_normal.var.groupby('gene_symbols'):
            if len(group) > 1:
                keep_idx.append(group['_mean_expr'].idxmax())
            else:
                keep_idx.append(group.index[0])
        adata_normal = adata_normal[:, keep_idx].copy()
    adata_normal.var_names = adata_normal.var['gene_symbols'].values
    adata_normal.var_names_make_unique()

# ========================== 共同基因 + PCA ==========================
common_genes = list(set(adata_normal.var_names) & set(adata_aml5.var_names))
print(f"Common genes: {len(common_genes)}")

adata_normal = adata_normal[:, common_genes].copy()
adata_aml5 = adata_aml5[:, common_genes].copy()

adata_combined = sc.concat([adata_normal, adata_aml5], label="_source",
                             keys=["normal","disease"], index_unique="-", join="outer")
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)
sc.pp.highly_variable_genes(adata_combined, n_top_genes=N_TOP_GENES, flavor="seurat_v3")
sc.pp.scale(adata_combined)
sc.pp.pca(adata_combined, n_comps=50, use_highly_variable=True)

normal_mask = adata_combined.obs["_source"] == "normal"
adata_normal = adata_combined[normal_mask].copy()
adata_aml5 = adata_combined[~normal_mask].copy()

# ========================== 构建正常元细胞 ==========================
print("\n" + "="*60)
print("Building normal metacells")
print("="*60)

def build_metacells(adata, use_rep="X_pca", resolution=1.0, label=None):
    ad = adata.copy()
    rep_key = use_rep if use_rep in ad.obsm else "X_pca"
    sc.pp.neighbors(ad, use_rep=rep_key, n_neighbors=15)
    sc.tl.leiden(ad, resolution=resolution)
    n_meta = ad.obs['leiden'].nunique()
    clusters = ad.obs['leiden'].unique()
    X_list, obs_list = [], []
    for cl in clusters:
        mask = ad.obs['leiden'] == cl
        cells = ad[mask]
        X_cl = cells.X.mean(axis=0)
        if hasattr(X_cl, 'toarray'):
            X_cl = X_cl.toarray().flatten()
        X_list.append(np.asarray(X_cl).flatten())
        obs_list.append({'metacell_id': f"{label}_MC{cl}", 'n_cells': int(mask.sum()), 'leiden': cl})
    meta_X = np.vstack(X_list)
    meta_obs = pd.DataFrame(obs_list)
    meta_var = ad.var.copy()
    meta_ad = sc.AnnData(X=meta_X, obs=meta_obs, var=meta_var)
    if rep_key in ad.obsm:
        meta_rep = []
        for cl in clusters:
            mask = ad.obs['leiden'] == cl
            meta_rep.append(ad.obsm[rep_key][mask].mean(axis=0))
        meta_ad.obsm[rep_key] = np.vstack(meta_rep)
    return meta_ad

import anndata
meta_normal_list = []
for stage in STAGE_ORDER:
    if stage not in adata_normal.obs['AuthorCellType'].values:
        continue
    stage_ad = adata_normal[adata_normal.obs['AuthorCellType'] == stage].copy()
    n_cells = stage_ad.shape[0]
    if n_cells < 10:
        continue
    res = max(0.5, min(3.0, 30 * METACELL_RESOLUTION / max(1, n_cells / 50)))
    meta_stage = build_metacells(stage_ad, use_rep=USE_REP, resolution=res, label=stage)
    meta_stage.obs['AuthorCellType'] = stage
    meta_normal_list.append(meta_stage)
    print(f"  {stage}: {n_cells} cells -> {meta_stage.shape[0]} metacells")

adata_meta_normal = anndata.concat(meta_normal_list, label="_stage_batch", index_unique="-")
adata_meta_normal.obs['AuthorCellType'] = pd.Categorical(
    [s for ad in meta_normal_list for s in [ad.obs['AuthorCellType'].iloc[0]] * ad.shape[0]],
    categories=STAGE_ORDER, ordered=True
)
print(f"Total normal metacells: {adata_meta_normal.shape[0]}")

# ========================== scTDRP分析器 ==========================
analyzer = TDRPAnalyzer(
    normal_adata=adata_meta_normal,
    stage_key='AuthorCellType',
    terminal_stage=TERMINAL_STAGE,
    stage_order=STAGE_ORDER,
    use_rep=USE_REP,
    n_top_genes=N_TOP_GENES,
)
analyzer.build_ot_cost_map(metric=OT_METRIC, reg=OT_REG)

# ========================== 分析1: 恶性 vs 正常残留 vs 其他 ==========================
print("\n" + "="*60)
print("Analysis 1: Malignancy subgroups")
print("="*60)

malignancy_groups = {
    'Malignant Erythroid': adata_aml5[adata_aml5.obs['malignancy'] == 'Malignant Erythroid'].copy(),
    'Normal Residual': adata_aml5[adata_aml5.obs['malignancy'] == 'Normal Residual'].copy(),
    'Other': adata_aml5[adata_aml5.obs['malignancy'] == 'Other'].copy(),
}

tdi_malignancy = []
for name, adata_sub in malignancy_groups.items():
    if adata_sub.shape[0] == 0:
        continue
    analyzer.compute_tdi(adata_sub, metric=OT_METRIC, reg=OT_REG)
    tdi_df = analyzer.tdi_results.copy()
    tdi_df['malignancy'] = name
    tdi_df['n_cells'] = adata_sub.shape[0]
    tdi_malignancy.append(tdi_df)
    print(f"  {name}: n={adata_sub.shape[0]}, mean_TDI={tdi_df['TDI'].mean():.4f}, median={tdi_df['TDI'].median():.4f}")

tdi_malignancy_df = pd.concat(tdi_malignancy, ignore_index=True)
tdi_malignancy_df.to_csv(os.path.join(OUT_DIR, "tdi_by_malignancy.csv"), index=False)

# 统计检验
mal_tdi = tdi_malignancy_df[tdi_malignancy_df['malignancy'] == 'Malignant Erythroid']['TDI']
norm_tdi = tdi_malignancy_df[tdi_malignancy_df['malignancy'] == 'Normal Residual']['TDI']
other_tdi = tdi_malignancy_df[tdi_malignancy_df['malignancy'] == 'Other']['TDI']

_, p_mal_norm = stats.mannwhitneyu(mal_tdi, norm_tdi, alternative='two-sided')
_, p_mal_other = stats.mannwhitneyu(mal_tdi, other_tdi, alternative='two-sided')
print(f"\n  Mann-Whitney U:")
print(f"    Malignant vs Normal Residual: p={p_mal_norm:.2e}")
print(f"    Malignant vs Other: p={p_mal_other:.2e}")

# ========================== 分析2: 细胞周期 ==========================
print("\n" + "="*60)
print("Analysis 2: Cell cycle phases")
print("="*60)

# 只分析恶性细胞
data_malignant = adata_aml5[adata_aml5.obs['malignancy'] == 'Malignant Erythroid'].copy()

cycle_groups = {
    'G1': data_malignant[data_malignant.obs['phase'] == 'G1'].copy(),
    'S': data_malignant[data_malignant.obs['phase'] == 'S'].copy(),
    'G2M': data_malignant[data_malignant.obs['phase'] == 'G2M'].copy(),
}

tdi_cycle = []
for name, adata_sub in cycle_groups.items():
    if adata_sub.shape[0] == 0:
        continue
    analyzer.compute_tdi(adata_sub, metric=OT_METRIC, reg=OT_REG)
    tdi_df = analyzer.tdi_results.copy()
    tdi_df['phase'] = name
    tdi_df['n_cells'] = adata_sub.shape[0]
    tdi_cycle.append(tdi_df)
    print(f"  {name}: n={adata_sub.shape[0]}, mean_TDI={tdi_df['TDI'].mean():.4f}, median={tdi_df['TDI'].median():.4f}")

tdi_cycle_df = pd.concat(tdi_cycle, ignore_index=True)
tdi_cycle_df.to_csv(os.path.join(OUT_DIR, "tdi_by_cell_cycle.csv"), index=False)

# 统计检验
g1_tdi = tdi_cycle_df[tdi_cycle_df['phase'] == 'G1']['TDI']
s_tdi = tdi_cycle_df[tdi_cycle_df['phase'] == 'S']['TDI']
g2m_tdi = tdi_cycle_df[tdi_cycle_df['phase'] == 'G2M']['TDI']

_, p_g1_s = stats.mannwhitneyu(g1_tdi, s_tdi, alternative='two-sided')
_, p_g1_g2m = stats.mannwhitneyu(g1_tdi, g2m_tdi, alternative='two-sided')
_, p_s_g2m = stats.mannwhitneyu(s_tdi, g2m_tdi, alternative='two-sided')
print(f"\n  Mann-Whitney U (Malignant cells):")
print(f"    G1 vs S: p={p_g1_s:.2e}")
print(f"    G1 vs G2M: p={p_g1_g2m:.2e}")
print(f"    S vs G2M: p={p_s_g2m:.2e}")

# ========================== 分析3: CNV Score亚克隆 ==========================
print("\n" + "="*60)
print("Analysis 3: CNV score subclones")
print("="*60)

# 看CNV score的分布，按四分位数分组
cnv_scores = data_malignant.obs['cnv_score'].astype(float)
print(f"  CNV score range: {cnv_scores.min():.4f} - {cnv_scores.max():.4f}")
print(f"  Unique scores: {cnv_scores.nunique()}")
print(cnv_scores.value_counts().sort_index().to_string())

# 按CNV score中位数二分
cnv_median = cnv_scores.median()
data_malignant.obs['cnv_group'] = (cnv_scores > cnv_median).map({True: 'High_CNV', False: 'Low_CNV'})

cnv_groups = {
    'Low_CNV': data_malignant[data_malignant.obs['cnv_group'] == 'Low_CNV'].copy(),
    'High_CNV': data_malignant[data_malignant.obs['cnv_group'] == 'High_CNV'].copy(),
}

tdi_cnv = []
for name, adata_sub in cnv_groups.items():
    if adata_sub.shape[0] == 0:
        continue
    analyzer.compute_tdi(adata_sub, metric=OT_METRIC, reg=OT_REG)
    tdi_df = analyzer.tdi_results.copy()
    tdi_df['cnv_group'] = name
    tdi_df['n_cells'] = adata_sub.shape[0]
    tdi_cnv.append(tdi_df)
    print(f"  {name}: n={adata_sub.shape[0]}, mean_TDI={tdi_df['TDI'].mean():.4f}, median={tdi_df['TDI'].median():.4f}")

tdi_cnv_df = pd.concat(tdi_cnv, ignore_index=True)
tdi_cnv_df.to_csv(os.path.join(OUT_DIR, "tdi_by_cnv_group.csv"), index=False)

low_tdi = tdi_cnv_df[tdi_cnv_df['cnv_group'] == 'Low_CNV']['TDI']
high_tdi = tdi_cnv_df[tdi_cnv_df['cnv_group'] == 'High_CNV']['TDI']
_, p_cnv = stats.mannwhitneyu(low_tdi, high_tdi, alternative='two-sided')
print(f"\n  Low vs High CNV: p={p_cnv:.2e}")

# ========================== 分析4: Seurat Cluster与TDI ==========================
print("\n" + "="*60)
print("Analysis 4: Seurat clusters")
print("="*60)

# 只看恶性细胞中的cluster分布
cluster_groups = {}
for cl in data_malignant.obs['seurat_clusters'].unique():
    sub = data_malignant[data_malignant.obs['seurat_clusters'] == cl].copy()
    if sub.shape[0] >= 20:  # 至少20个细胞
        cluster_groups[f'Cluster_{cl}'] = sub

tdi_cluster = []
for name, adata_sub in cluster_groups.items():
    analyzer.compute_tdi(adata_sub, metric=OT_METRIC, reg=OT_REG)
    tdi_df = analyzer.tdi_results.copy()
    tdi_df['cluster'] = name
    tdi_df['n_cells'] = adata_sub.shape[0]
    tdi_cluster.append(tdi_df)
    print(f"  {name}: n={adata_sub.shape[0]}, mean_TDI={tdi_df['TDI'].mean():.4f}, median={tdi_df['TDI'].median():.4f}")

tdi_cluster_df = pd.concat(tdi_cluster, ignore_index=True)
tdi_cluster_df.to_csv(os.path.join(OUT_DIR, "tdi_by_cluster.csv"), index=False)

# ========================== 可视化 ==========================
print("\n" + "="*60)
print("Visualization")
print("="*60)

fig = plt.figure(figsize=(16, 12))

# Panel 1: Malignancy
ax1 = fig.add_subplot(2, 2, 1)
sns.boxplot(data=tdi_malignancy_df, x='malignancy', y='TDI', 
            palette=['#e74c3c', '#3498db', '#95a5a6'], ax=ax1)
ax1.set_title(f'A. TDI by Malignancy Status\n(Mal vs Norm: p={p_mal_norm:.2e})', fontsize=12, fontweight='bold')
ax1.set_xlabel('')
ax1.tick_params(axis='x', rotation=15)

# Panel 2: Cell Cycle
ax2 = fig.add_subplot(2, 2, 2)
sns.boxplot(data=tdi_cycle_df, x='phase', y='TDI', 
            palette=['#2ecc71', '#f39c12', '#9b59b6'], ax=ax2)
ax2.set_title(f'B. TDI by Cell Cycle (Malignant Cells)\n(G1 vs G2M: p={p_g1_g2m:.2e})', fontsize=12, fontweight='bold')
ax2.set_xlabel('')

# Panel 3: CNV Group
ax3 = fig.add_subplot(2, 2, 3)
sns.boxplot(data=tdi_cnv_df, x='cnv_group', y='TDI', 
            palette=['#3498db', '#e74c3c'], ax=ax3)
ax3.set_title(f'C. TDI by CNV Load (Malignant Cells)\n(p={p_cnv:.2e})', fontsize=12, fontweight='bold')
ax3.set_xlabel('')

# Panel 4: Cluster
ax4 = fig.add_subplot(2, 2, 4)
cluster_means = tdi_cluster_df.groupby('cluster')['TDI'].mean().sort_values(ascending=False)
sns.barplot(x=cluster_means.values, y=cluster_means.index, palette='RdYlBu_r', ax=ax4)
ax4.set_title('D. Mean TDI by Seurat Cluster\n(Malignant Cells)', fontsize=12, fontweight='bold')
ax4.set_xlabel('Mean TDI')

plt.tight_layout(pad=2.0)
plt.savefig(os.path.join(FIG_DIR, "aml5_deep_dive.pdf"), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, "aml5_deep_dive.png"), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: aml5_deep_dive")

# ========================== 摘要 ==========================
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"\n1. Malignancy:")
for name in ['Malignant Erythroid', 'Normal Residual', 'Other']:
    subset = tdi_malignancy_df[tdi_malignancy_df['malignancy'] == name]
    if len(subset) > 0:
        print(f"   {name}: n={len(subset)}, mean={subset['TDI'].mean():.4f}, median={subset['TDI'].median():.4f}")

print(f"\n2. Cell Cycle (Malignant):")
for name in ['G1', 'S', 'G2M']:
    subset = tdi_cycle_df[tdi_cycle_df['phase'] == name]
    if len(subset) > 0:
        print(f"   {name}: n={len(subset)}, mean={subset['TDI'].mean():.4f}, median={subset['TDI'].median():.4f}")

print(f"\n3. CNV Load (Malignant):")
for name in ['Low_CNV', 'High_CNV']:
    subset = tdi_cnv_df[tdi_cnv_df['cnv_group'] == name]
    if len(subset) > 0:
        print(f"   {name}: n={len(subset)}, mean={subset['TDI'].mean():.4f}, median={subset['TDI'].median():.4f}")

print(f"\n4. Top clusters by TDI:")
for cl, tdi in cluster_means.head(5).items():
    print(f"   {cl}: mean_TDI={tdi:.4f}")

print("\n" + "="*60)
print(f"Done! Results -> {OUT_DIR}")
print(f"Figures -> {FIG_DIR}")
print("="*60)
