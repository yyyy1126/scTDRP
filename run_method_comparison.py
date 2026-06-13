#!/usr/bin/env python3
"""
Method Comparison: scTDRP TDI vs Baseline Methods on AEL Dataset
Corrected alignment: TDI CSV stores cell_idx per malignancy group, not globally.
"""

import os
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import cdist
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

AEL_PATH = "/data1/yja/zhongzhuan/4.results/4.infercnv/infercnv_aml5/aml5_annotated_with_cnv_cycle.h5ad"
NORMAL_PATH = "/data1/yja/zhongzhuan/1.data/processed/erythroid_lineage_from_MEP.h5ad"
MODULES_PATH = "/data1/yja/zhongzhuan/5.external/scTDRP/modules.json"
TDI_PATH = "/data1/yja/zhongzhuan/5.external/scTDRP/results_aml5_deep/tdi_by_malignancy.csv"
OUT_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/results_method_comparison"
FIG_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/figures_method_comparison"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

np.random.seed(42)

print("=" * 60)
print("Method Comparison: scTDRP vs Baselines on AEL")
print("=" * 60)

# ========================== 1. Load Data ==========================
print("\n[1] Loading data...")
adata_ael = sc.read_h5ad(AEL_PATH)
adata_normal = sc.read_h5ad(NORMAL_PATH)
print(f"  AEL: {adata_ael.shape}")
print(f"  Normal: {adata_normal.shape}")

# Save normal stage labels before manipulation
normal_stage_labels = adata_normal.obs['AuthorCellType'].copy()

# ========================== 2. Gene Mapping ==========================
print("\n[2] Gene mapping...")
adata_normal.var['gene_symbol'] = adata_normal.var['gene_symbols'].astype(str)
mean_expr = np.array(adata_normal.X.mean(axis=0)).flatten()
if hasattr(mean_expr, 'toarray'):
    mean_expr = mean_expr.toarray().flatten()
adata_normal.var['_mean_expr'] = mean_expr

keep_idx = []
for sym, group in adata_normal.var.groupby('gene_symbol'):
    if len(group) > 1:
        keep_idx.append(group['_mean_expr'].idxmax())
    else:
        keep_idx.append(group.index[0])
adata_normal = adata_normal[:, keep_idx].copy()
adata_normal.var_names = adata_normal.var['gene_symbol'].values
adata_normal.var_names_make_unique()

common_genes = list(set(adata_ael.var_names) & set(adata_normal.var_names))
print(f"  Common genes: {len(common_genes)}")
adata_ael = adata_ael[:, common_genes].copy()
adata_normal = adata_normal[:, common_genes].copy()

# ========================== 3. Preprocessing ==========================
print("\n[3] Preprocessing...")
adata_ael.X = adata_ael.layers['logcounts']

adata_combined = sc.concat([adata_normal, adata_ael], label='_source',
                            keys=['normal', 'disease'], index_unique='-')
sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, flavor='seurat')
sc.pp.scale(adata_combined, max_value=10)
sc.pp.pca(adata_combined, n_comps=30, use_highly_variable=True)

mask_normal = adata_combined.obs['_source'] == 'normal'
mask_disease = adata_combined.obs['_source'] == 'disease'
normal_pca = adata_combined.obsm['X_pca'][mask_normal]
disease_pca = adata_combined.obsm['X_pca'][mask_disease]
print(f"  PCA: normal={normal_pca.shape}, disease={disease_pca.shape}")

# ========================== 4. Load & Align TDI ==========================
print("\n[4] Loading and aligning TDI...")
tdi_df = pd.read_csv(TDI_PATH)

# TDI CSV stores cell_idx per malignancy group. We need to align by group order.
# The order in TDI CSV for each group matches the order of cells in adata_ael for that group.
ael_tdi = np.zeros(adata_ael.shape[0])
for mal_group in tdi_df['malignancy'].unique():
    group_tdi = tdi_df[tdi_df['malignancy'] == mal_group].sort_values('cell_idx').reset_index(drop=True)
    # Find cells in adata_ael belonging to this group, in original order
    mask = adata_ael.obs['malignancy'] == mal_group
    n_cells = mask.sum()
    n_tdi = len(group_tdi)
    if n_cells != n_tdi:
        print(f"  WARNING: {mal_group} - adata has {n_cells} cells, TDI has {n_tdi} rows")
        n_match = min(n_cells, n_tdi)
    else:
        n_match = n_cells
    ael_tdi[mask.values] = group_tdi['TDI'].values[:n_match]
    print(f"  {mal_group}: aligned {n_match} cells, mean TDI = {group_tdi['TDI'].mean():.3f}")

print(f"  Overall TDI mean: {ael_tdi.mean():.3f}")

# ========================== 5. Baseline 1: PCA Euclidean Distance ==========================
print("\n[5] Baseline 1: PCA Euclidean Distance to nearest normal stage centroid...")
stage_order = ['MEP', 'BFU-E', 'CFU-E', 'Pro-Erythroblast', 'Basophilic Erythroblast',
               'Polychromatic Erythroblast', 'Orthochromatic Erythroblast']
stage_centroids = {}
for stage in stage_order:
    mask = normal_stage_labels.values == stage
    if mask.sum() > 0:
        stage_centroids[stage] = normal_pca[mask].mean(axis=0)

centroid_matrix = np.array(list(stage_centroids.values()))
distances = cdist(disease_pca, centroid_matrix, metric='euclidean')
pca_euc_dist = distances.min(axis=1)
print(f"  PCA EucDist: mean={pca_euc_dist.mean():.3f}")

# ========================== 6. Baseline 2: Gene-set Differentiation Score ==========================
print("\n[6] Baseline 2: Gene-set differentiation score...")
with open(MODULES_PATH) as f:
    modules = json.load(f)
p7_genes = [g for g in modules['P7_TerminalPrep'] if g in adata_ael.var_names]
p8_genes = [g for g in modules['P8_ExecutionPrecursor'] if g in adata_ael.var_names]
print(f"  P7: {len(p7_genes)}/100, P8: {len(p8_genes)}/100")

sc.tl.score_genes(adata_ael, gene_list=p7_genes, score_name='P7_score')
sc.tl.score_genes(adata_ael, gene_list=p8_genes, score_name='P8_score')

def minmax_scale(x):
    return (x - x.min()) / (x.max() - x.min() + 1e-10)

p7 = adata_ael.obs['P7_score'].values
p8 = adata_ael.obs['P8_score'].values
geneset_diff_score = minmax_scale(minmax_scale(p7) - minmax_scale(p8))
print(f"  DiffScore: mean={geneset_diff_score.mean():.3f}")

# ========================== 7. Baseline 3: Pseudotime Deviation ==========================
print("\n[7] Baseline 3: Pseudotime deviation...")
normal_pc1 = normal_pca[:, 0]
stage_pc1_means = []
for stage in stage_order:
    mask = normal_stage_labels.values == stage
    stage_pc1_means.append(normal_pc1[mask].mean() if mask.sum() > 0 else np.nan)

stage_corr = stats.spearmanr(range(len(stage_order)), stage_pc1_means)[0]
print(f"  PC1 vs stage order: r={stage_corr:.3f}")

if stage_corr < 0:
    ael_pc1 = -disease_pca[:, 0]
    normal_pc1 = -normal_pc1
else:
    ael_pc1 = disease_pca[:, 0]

pt_min, pt_max = normal_pc1.min(), normal_pc1.max()
ael_pseudotime = np.clip((ael_pc1 - pt_min) / (pt_max - pt_min + 1e-10), 0, 1)
pseudotime_dev = 1.0 - ael_pseudotime
print(f"  PseudotimeDev: mean={pseudotime_dev.mean():.3f}")

# ========================== 8. Comparison ==========================
print("\n[8] Computing metrics...")

malignancy = adata_ael.obs['malignancy'].values
phase = adata_ael.obs['phase'].values
cnv_score = adata_ael.obs['cnv_score'].values
cnv_median = np.median(cnv_score[malignancy == 'Malignant Erythroid'])
cnv_group = np.array(['High_CNV' if s > cnv_median else 'Low_CNV' for s in cnv_score])
cycle_score = np.array([{'G1': 0.0, 'S': 1.0, 'G2M': 2.0}.get(p, np.nan) for p in phase])

comp_df = pd.DataFrame({
    'malignancy': malignancy, 'phase': phase, 'cnv_score': cnv_score,
    'cnv_group': cnv_group, 'cycle_score': cycle_score,
    'TDI': ael_tdi, 'PCA_EucDist': pca_euc_dist,
    'GeneSet_DiffScore': geneset_diff_score, 'Pseudotime_Dev': pseudotime_dev,
})

methods = ['TDI', 'PCA_EucDist', 'GeneSet_DiffScore', 'Pseudotime_Dev']
results = []

# AUC: Malignant vs Normal
mask_malnorm = comp_df['malignancy'].isin(['Malignant Erythroid', 'Normal Residual'])
y_true_malnorm = (comp_df.loc[mask_malnorm, 'malignancy'] == 'Malignant Erythroid').astype(int)
print("\n  AUC (Malignant vs Normal Residual):")
for method in methods:
    auc = roc_auc_score(y_true_malnorm, comp_df.loc[mask_malnorm, method])
    results.append({'Metric': 'AUC_Malignant_vs_Normal', 'Method': method, 'Value': auc})
    print(f"    {method}: {auc:.3f}")

# AUC: High vs Low CNV
mask_mal = comp_df['malignancy'] == 'Malignant Erythroid'
y_true_cnv = (comp_df.loc[mask_mal, 'cnv_group'] == 'High_CNV').astype(int)
print("\n  AUC (High-CNV vs Low-CNV, malignant):")
for method in methods:
    auc = roc_auc_score(y_true_cnv, comp_df.loc[mask_mal, method])
    results.append({'Metric': 'AUC_HighCNV_vs_LowCNV', 'Method': method, 'Value': auc})
    print(f"    {method}: {auc:.3f}")

# Spearman: Cycle
mask_valid = mask_mal & (~np.isnan(cycle_score))
print("\n  Spearman r (Cell Cycle, malignant):")
for method in methods:
    r, p = stats.spearmanr(comp_df.loc[mask_valid, method], comp_df.loc[mask_valid, 'cycle_score'])
    results.append({'Metric': 'Spearman_CycleScore', 'Method': method, 'Value': r, 'Pvalue': p})
    print(f"    {method}: r={r:+.3f}, p={p:.2e}")

# Spearman: CNV
print("\n  Spearman r (CNV Score, malignant):")
for method in methods:
    r, p = stats.spearmanr(comp_df.loc[mask_mal, method], comp_df.loc[mask_mal, 'cnv_score'])
    results.append({'Metric': 'Spearman_CNVScore', 'Method': method, 'Value': r, 'Pvalue': p})
    print(f"    {method}: r={r:+.3f}, p={p:.2e}")

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(OUT_DIR, 'method_comparison_metrics.csv'), index=False)
comp_df.to_csv(os.path.join(OUT_DIR, 'method_comparison_scores.csv'), index=False)

# ========================== 9. Visualization ==========================
print("\n[9] Plotting...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
auc_mal = results_df[results_df['Metric'] == 'AUC_Malignant_vs_Normal']
sns.barplot(data=auc_mal, x='Method', y='Value', palette=['#e74c3c', '#3498db', '#2ecc71', '#f39c12'], ax=axes[0])
axes[0].set_ylim(0.5, 1.0)
axes[0].axhline(0.5, color='gray', linestyle='--', alpha=0.5)
axes[0].set_title("AUC: Malignant vs Normal", fontsize=13, fontweight='bold')
for i, (_, row) in enumerate(auc_mal.iterrows()):
    axes[0].text(i, row['Value'] + 0.01, f"{row['Value']:.3f}", ha='center', fontsize=10, fontweight='bold')

auc_cnv = results_df[results_df['Metric'] == 'AUC_HighCNV_vs_LowCNV']
sns.barplot(data=auc_cnv, x='Method', y='Value', palette=['#e74c3c', '#3498db', '#2ecc71', '#f39c12'], ax=axes[1])
axes[1].set_ylim(0.5, 1.0)
axes[1].axhline(0.5, color='gray', linestyle='--', alpha=0.5)
axes[1].set_title("AUC: High-CNV vs Low-CNV", fontsize=13, fontweight='bold')
for i, (_, row) in enumerate(auc_cnv.iterrows()):
    axes[1].text(i, row['Value'] + 0.01, f"{row['Value']:.3f}", ha='center', fontsize=10, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'FigureS2_MethodComparison_AUC.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, 'FigureS2_MethodComparison_AUC.png'), dpi=300, bbox_inches='tight')
plt.close()

corr_data = results_df[results_df['Metric'].isin(['Spearman_CycleScore', 'Spearman_CNVScore'])].copy()
corr_pivot = corr_data.pivot(index='Method', columns='Metric', values='Value')
fig, ax = plt.subplots(figsize=(8, 5))
sns.heatmap(corr_pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0, vmin=-1, vmax=1, ax=ax, linewidths=0.5)
ax.set_title("Spearman Correlation", fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'FigureS2_MethodComparison_Correlation.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, 'FigureS2_MethodComparison_Correlation.png'), dpi=300, bbox_inches='tight')
plt.close()

fig, axes = plt.subplots(1, 4, figsize=(18, 4))
for ax, method in zip(axes, methods):
    sns.boxplot(data=comp_df[mask_malnorm], x='malignancy', y=method, palette=['#2ecc71', '#e74c3c'], ax=ax)
    ax.set_title(method, fontsize=11, fontweight='bold')
    ax.set_xticklabels(['Normal', 'Malignant'], rotation=15)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'FigureS2_MethodComparison_Boxplots.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, 'FigureS2_MethodComparison_Boxplots.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"\nDone! -> {OUT_DIR}")
print("=" * 60)
