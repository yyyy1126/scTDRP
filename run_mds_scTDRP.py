#!/usr/bin/env python3
"""
scTDRP MDS-5q 无效红细胞生成分析脚本

数据来源: Doty et al. 2023, Blood Advances (GSE222368)
- 正常对照: N1-N7 (骨髓细胞红系扩增培养 Day0/3/6)
- MDS-5q: M1-M2 (5q-综合征患者)
- DBA: D1-D2 (Diamond-Blackfan贫血, 可选对比)

分析目标: 验证MDS-5q红系细胞的分化阻滞和TDI升高
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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from scTDRP import TDRPAnalyzer
from scTDRP.distance import compute_wasserstein_distance
from scTDRP.repair import compute_gene_repair_delta
from scTDRP.module import module_repair_strategy
from scTDRP.io import save_results

warnings.filterwarnings("ignore")

# ========================== 配置 ==========================
DATA_PATH = "../../1.data/raw/scRNAseq/Doty2023_MDS/GSE222368_heme_velo_aggr_5exp_211012.h5ad"
MODULES_PATH = "./modules.json"
OUT_DIR = "./results_mds"
FIG_DIR = "./figures_mds"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

STAGE_ORDER = ["BFU-E", "CFU-E", "Pro-Erythroblast", "Basophilic Erythroblast",
               "Polychromatic Erythroblast", "Orthochromatic Erythroblast"]
TERMINAL_STAGE = "Orthochromatic Erythroblast"

OT_METRIC = "sqeuclidean"
OT_REG = 0.05
TOP_N_REPAIR = 100
USE_REP = "X_pca"
N_TOP_GENES = 2000
N_PCS = 50

METACELL_RESOLUTION = 1.0

# CITE-seq阈值
ERY_THRESHOLD_CD71 = 10
ERY_THRESHOLD_CD235A = 5
ERY_THRESHOLD_CD36 = 10
NON_ERY_THRESHOLD_CD3 = 1
NON_ERY_THRESHOLD_CD19 = 1
NON_ERY_THRESHOLD_CD11B = 5
NON_ERY_THRESHOLD_CD14 = 1

# ========================== 辅助函数 ==========================
def build_metacells(adata, use_rep="X_pca", resolution=1.0, label=None):
    """用Leiden聚类聚合元细胞"""
    ad = adata.copy()
    rep_key = use_rep if use_rep in ad.obsm else "X_pca"
    sc.pp.neighbors(ad, use_rep=rep_key, n_neighbors=15)
    sc.tl.leiden(ad, resolution=resolution)
    
    n_meta = ad.obs['leiden'].nunique()
    print(f"    {label}: {ad.shape[0]} cells -> {n_meta} metacells")
    
    clusters = ad.obs['leiden'].unique()
    X_list, obs_list = [], []
    for cl in clusters:
        mask = ad.obs['leiden'] == cl
        cells = ad[mask]
        X_cl = cells.X.mean(axis=0)
        if hasattr(X_cl, 'toarray'):
            X_cl = X_cl.toarray().flatten()
        X_list.append(np.asarray(X_cl).flatten())
        obs_list.append({
            'metacell_id': f"{label}_MC{cl}",
            'n_cells': int(mask.sum()),
            'leiden': cl,
        })
    
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
    
    return meta_ad, ad.obs['leiden'].to_dict()


def assign_stages_by_maturity(adata, stage_key='stage', inplace=True):
    """
    基于CITE-seq表面蛋白标记计算成熟度分数，分6个阶段
    
    maturity_score = log1p(CD235a) + log1p(CD71) - log1p(CD117)
    分数越高 = 越成熟
    """
    ad = adata if inplace else adata.copy()
    
    # 计算成熟度分数
    cd235a = np.log1p(ad.obs['tsb-CD235a'].astype(float).values)
    cd71 = np.log1p(ad.obs['tsb-CD71'].astype(float).values)
    cd117 = np.log1p(ad.obs['tsb-CD117'].astype(float).values)
    
    ad.obs['maturity_score'] = cd235a + cd71 - cd117
    
    # 分6个quantile bin
    ad.obs[stage_key] = pd.qcut(ad.obs['maturity_score'], q=6, labels=STAGE_ORDER)
    
    # 统计
    print(f"\n  Stage distribution:")
    print(ad.obs[stage_key].value_counts().sort_index().to_string())
    
    return ad


# ========================== 步骤1: 加载数据 ==========================
print("="*60)
print("Step 1: Loading GSE222368 data")
print("="*60)

adata = sc.read_h5ad(DATA_PATH)
print(f"  Loaded: {adata.shape}")
print(f"  obs: {list(adata.obs.columns)}")

# ========================== 步骤2: 基础质控和过滤 ==========================
print("\n" + "="*60)
print("Step 2: QC and filtering")
print("="*60)

# 排除未分配/双细胞
adata = adata[~adata.obs['id'].isin(['unassigned', 'doublet', 'unknown'])].copy()
print(f"  After removing unassigned/doublet/unknown: {adata.shape}")

# ID映射到condition
id_to_condition = {
    'N1': 'Normal', 'N2': 'Normal', 'N3': 'Normal', 'N4': 'Normal',
    'N5': 'Normal', 'N6': 'Normal', 'N7': 'Normal',
    'D1': 'DBA', 'D2': 'DBA',
    'M1': 'MDS-5q', 'M2': 'MDS-5q'
}
adata.obs['condition'] = adata.obs['id'].map(id_to_condition)
print(f"\n  Condition distribution:")
print(adata.obs['condition'].value_counts().to_string())

# ========================== 步骤3: 红系细胞筛选 ==========================
print("\n" + "="*60)
print("Step 3: Erythroid cell selection (CITE-seq based)")
print("="*60)

# 提取CITE-seq标记
cd71 = adata.obs['tsb-CD71'].astype(float).values
cd235a = adata.obs['tsb-CD235a'].astype(float).values
cd36 = adata.obs['tsb-CD36'].astype(float).values
cd117 = adata.obs['tsb-CD117'].astype(float).values

cd3 = adata.obs['tsb-CD3'].astype(float).values
cd19 = adata.obs['tsb-CD19'].astype(float).values
cd11b = adata.obs['tsb-CD11b'].astype(float).values
cd14 = adata.obs['tsb-CD14'].astype(float).values

# 红系阳性: CD71+ 或 CD235a+ 或 CD36+
ery_mask = (
    (cd71 > ERY_THRESHOLD_CD71) |
    (cd235a > ERY_THRESHOLD_CD235A) |
    (cd36 > ERY_THRESHOLD_CD36)
)

# 排除明显的非红系: CD3+, CD19+, CD11b+, CD14+
non_ery_mask = (
    (cd3 > NON_ERY_THRESHOLD_CD3) |
    (cd19 > NON_ERY_THRESHOLD_CD19) |
    (cd11b > NON_ERY_THRESHOLD_CD11B) |
    (cd14 > NON_ERY_THRESHOLD_CD14)
)

final_mask = ery_mask & (~non_ery_mask)
adata_ery = adata[final_mask].copy()
print(f"  Erythroid cells selected: {adata_ery.shape[0]} / {adata.shape[0]}")

# 各condition的红系细胞数
print(f"\n  By condition:")
print(adata_ery.obs['condition'].value_counts().to_string())

# ========================== 步骤4: 阶段定义 ==========================
print("\n" + "="*60)
print("Step 4: Stage assignment by CITE-seq maturity score")
print("="*60)

# 仅用Normal样本定义阶段边界
adata_normal = adata_ery[adata_ery.obs['condition'] == 'Normal'].copy()
adata_normal = assign_stages_by_maturity(adata_normal, stage_key='stage', inplace=True)

# 获取quantile边界（用于给疾病细胞分阶段）
maturity_q = adata_normal.obs['maturity_score'].quantile([0, 1/6, 2/6, 3/6, 4/6, 5/6, 1]).values
print(f"\n  Maturity score quantiles (from Normal): {maturity_q}")

# 给所有细胞分阶段
def assign_stage_from_quantiles(score, quantiles, labels):
    """根据Normal样本的quantile边界给所有细胞分阶段"""
    stages = []
    for s in score:
        for i in range(len(quantiles)-1):
            if quantiles[i] <= s < quantiles[i+1]:
                stages.append(labels[i])
                break
        else:
            stages.append(labels[-1])
    return stages

adata_ery.obs['maturity_score'] = (
    np.log1p(adata_ery.obs['tsb-CD235a'].astype(float).values) +
    np.log1p(adata_ery.obs['tsb-CD71'].astype(float).values) -
    np.log1p(adata_ery.obs['tsb-CD117'].astype(float).values)
)
adata_ery.obs['stage'] = assign_stage_from_quantiles(
    adata_ery.obs['maturity_score'].values, maturity_q, STAGE_ORDER
)
adata_ery.obs['stage'] = pd.Categorical(adata_ery.obs['stage'], categories=STAGE_ORDER, ordered=True)

print(f"\n  All cells stage distribution:")
print(adata_ery.obs['stage'].value_counts().sort_index().to_string())

# ========================== 步骤5: 数据预处理 ==========================
print("\n" + "="*60)
print("Step 5: Preprocessing")
print("="*60)

# 使用spliced counts作为基因表达（最可靠）
if 'spliced' in adata_ery.layers:
    adata_ery.X = adata_ery.layers['spliced'].copy()
    print(f"  Using spliced counts as X")

# 标准化
sc.pp.normalize_total(adata_ery, target_sum=1e4)
sc.pp.log1p(adata_ery)

# 高变基因
sc.pp.highly_variable_genes(adata_ery, n_top_genes=N_TOP_GENES, flavor="seurat_v3")

# 对全体做PCA（用于OT计算）
sc.pp.scale(adata_ery, max_value=10)
sc.pp.pca(adata_ery, n_comps=N_PCS, use_highly_variable=True)
print(f"  PCA done: {adata_ery.obsm['X_pca'].shape}")

# ========================== 步骤6: 分离Normal和Disease ==========================
print("\n" + "="*60)
print("Step 6: Split Normal and Disease")
print("="*60)

adata_normal = adata_ery[adata_ery.obs['condition'] == 'Normal'].copy()
adata_mds = adata_ery[adata_ery.obs['condition'] == 'MDS-5q'].copy()
adata_dba = adata_ery[adata_ery.obs['condition'] == 'DBA'].copy()

print(f"  Normal: {adata_normal.shape}")
print(f"  MDS-5q: {adata_mds.shape}")
print(f"  DBA: {adata_dba.shape}")

# 模块基因
with open(MODULES_PATH, "r") as f:
    modules = json.load(f)
module_genes = list(set(modules["P7_TerminalPrep"] + modules["P8_ExecutionPrecursor"]))

# ========================== 步骤7: 构建元细胞 ==========================
print("\n" + "="*60)
print("Step 7: Building metacells")
print("="*60)

meta_normal_list = []
for stage in STAGE_ORDER:
    if stage not in adata_normal.obs['stage'].values:
        continue
    stage_ad = adata_normal[adata_normal.obs['stage'] == stage].copy()
    n_cells = stage_ad.shape[0]
    if n_cells < 10:
        print(f"    Skipping {stage}: only {n_cells} cells")
        continue
    res = max(0.5, min(3.0, 30 * METACELL_RESOLUTION / max(1, n_cells / 50)))
    meta_stage, _ = build_metacells(stage_ad, use_rep=USE_REP, resolution=res, label=stage)
    meta_stage.obs['stage'] = stage
    meta_normal_list.append(meta_stage)
    print(f"      {stage}: {n_cells} cells -> {meta_stage.shape[0]} metacells")

adata_meta_normal = anndata.concat(meta_normal_list, label="_stage_batch", index_unique="-")
adata_meta_normal.obs['stage'] = pd.Categorical(
    [s for ad in meta_normal_list for s in [ad.obs['stage'].iloc[0]] * ad.shape[0]],
    categories=STAGE_ORDER, ordered=True
)
print(f"  Total normal metacells: {adata_meta_normal.shape[0]}")

# MDS元细胞
meta_mds, _ = build_metacells(adata_mds, use_rep=USE_REP,
                               resolution=METACELL_RESOLUTION, label="MDS")
print(f"  Total MDS metacells: {meta_mds.shape[0]}")

# DBA元细胞（可选对比）
meta_dba, _ = build_metacells(adata_dba, use_rep=USE_REP,
                               resolution=METACELL_RESOLUTION, label="DBA")
print(f"  Total DBA metacells: {meta_dba.shape[0]}")

# ========================== 步骤8: scTDRP分析 ==========================
print("\n" + "="*60)
print("Step 8: Running scTDRP")
print("="*60)

analyzer = TDRPAnalyzer(
    normal_adata=adata_meta_normal,
    stage_key='stage',
    terminal_stage=TERMINAL_STAGE,
    stage_order=STAGE_ORDER,
    use_rep=USE_REP,
    n_top_genes=N_TOP_GENES,
)
analyzer.build_ot_cost_map(metric=OT_METRIC, reg=OT_REG)

# 计算MDS的TDI
print("\n  Computing TDI for MDS...")
analyzer.compute_tdi(meta_mds, metric=OT_METRIC, reg=OT_REG)
tdi_mds = analyzer.tdi_results.copy()
tdi_mds['condition'] = 'MDS-5q'

# 计算DBA的TDI
print("\n  Computing TDI for DBA...")
analyzer.compute_tdi(meta_dba, metric=OT_METRIC, reg=OT_REG)
tdi_dba = analyzer.tdi_results.copy()
tdi_dba['condition'] = 'DBA'

# 计算Normal的TDI（作为对照）
print("\n  Computing TDI for Normal...")
# 对Normal的元细胞也计算TDI
analyzer.compute_tdi(adata_meta_normal, metric=OT_METRIC, reg=OT_REG)
tdi_normal = analyzer.tdi_results.copy()
tdi_normal['condition'] = 'Normal'

# 合并TDI结果
tdi_all = pd.concat([tdi_normal, tdi_mds, tdi_dba], ignore_index=True)

# ========================== 步骤9: 基因级修复路径 ==========================
print("\n" + "="*60)
print("Step 9: Gene-level repair pathway (MDS -> Terminal)")
print("="*60)

common_genes = list(set(adata_meta_normal.var_names) & set(meta_mds.var_names))
covered_module = [g for g in module_genes if g in common_genes]
print(f"  Common genes: {len(common_genes)}, Module covered: {len(covered_module)}")

terminal_meta = adata_meta_normal[adata_meta_normal.obs['stage'] == TERMINAL_STAGE]

X_disease = meta_mds[:, covered_module].X
X_terminal = terminal_meta[:, covered_module].X
if hasattr(X_disease, "toarray"):
    X_disease = X_disease.toarray()
if hasattr(X_terminal, "toarray"):
    X_terminal = X_terminal.toarray()
X_disease = np.asarray(X_disease, dtype=np.float64)
X_terminal = np.asarray(X_terminal, dtype=np.float64)

w_disease = meta_mds.obs['n_cells'].values.astype(np.float64)
w_disease = w_disease / w_disease.sum()
w_terminal = terminal_meta.obs['n_cells'].values.astype(np.float64)
w_terminal = w_terminal / w_terminal.sum()

dist, transport_plan = compute_wasserstein_distance(
    X_disease, X_terminal, w_disease, w_terminal,
    metric=OT_METRIC, reg=OT_REG, method="sinkhorn"
)
print(f"  Wasserstein Distance (MDS -> Terminal): {dist:.4f}")

repair_deltas = compute_gene_repair_delta(
    X_disease, X_terminal, transport_plan, gene_names=covered_module
)
sorted_repair = sorted(repair_deltas.items(), key=lambda x: x[1], reverse=True)
top_up = [g for g, d in sorted_repair[:TOP_N_REPAIR] if d > 0]
top_down = [g for g, d in sorted_repair[-TOP_N_REPAIR:] if d < 0]

analyzer.repair_results = {
    "wasserstein_distance": dist,
    "transport_plan": transport_plan,
    "repair_deltas": repair_deltas,
    "top_up_targets": top_up,
    "top_down_targets": top_down,
}

# ========================== 步骤10: 模块策略 ==========================
print("\n" + "="*60)
print("Step 10: Module repair strategy")
print("="*60)

analyzer.module_strategy = module_repair_strategy(
    repair_deltas=repair_deltas, module_dict=modules, threshold=0.0, method="mean"
)
print(analyzer.module_strategy.to_string(index=False))

# ========================== 步骤11: 可视化 ==========================
print("\n" + "="*60)
print("Step 11: Visualization")
print("="*60)

# 1. OT Cost Map
fig, ax = plt.subplots(figsize=(10, 8))
cost_matrix = np.zeros((len(STAGE_ORDER), len(STAGE_ORDER)))
for i, s1 in enumerate(STAGE_ORDER):
    for j, s2 in enumerate(STAGE_ORDER):
        if (s1, s2) in analyzer.cost_map:
            cost_matrix[i, j] = analyzer.cost_map[(s1, s2)]

sns.heatmap(cost_matrix, xticklabels=STAGE_ORDER, yticklabels=STAGE_ORDER,
            annot=True, fmt=".3f", cmap="YlOrRd", ax=ax, square=True)
ax.set_title("OT Cost Map (Normal Erythroid Metacells)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "01_ot_cost_map.pdf"))
plt.savefig(os.path.join(FIG_DIR, "01_ot_cost_map.png"), dpi=300)
plt.close()
print("  Saved: 01_ot_cost_map")

# 2. TDI分布比较
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# TDI violin
sns.violinplot(data=tdi_all, x='condition', y='TDI', palette=['#2ecc71', '#e74c3c', '#3498db'], ax=axes[0])
axes[0].set_title("TDI Distribution by Condition")
axes[0].set_ylabel("TDI (Wasserstein Distance)")

# TDI箱线图
sns.boxplot(data=tdi_all, x='condition', y='TDI', palette=['#2ecc71', '#e74c3c', '#3498db'], ax=axes[1])
axes[1].set_title("TDI Boxplot by Condition")
axes[1].set_ylabel("TDI (Wasserstein Distance)")

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "02_tdi_distribution.pdf"))
plt.savefig(os.path.join(FIG_DIR, "02_tdi_distribution.png"), dpi=300)
plt.close()
print("  Saved: 02_tdi_distribution")

# 3. Best Match Stage
fig, ax = plt.subplots(figsize=(10, 6))
stage_props = tdi_all.groupby('condition')['Best_Match_Stage'].value_counts(normalize=True).unstack(fill_value=0)
stage_props = stage_props.reindex(columns=STAGE_ORDER, fill_value=0)
stage_props.plot(kind='bar', stacked=True, ax=ax, colormap='RdYlGn_r')
ax.set_title("Best Match Stage Distribution")
ax.set_ylabel("Proportion")
ax.legend(title="Stage", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "03_best_match_stage.pdf"))
plt.savefig(os.path.join(FIG_DIR, "03_best_match_stage.png"), dpi=300)
plt.close()
print("  Saved: 03_best_match_stage")

# 4. MDS trajectory B vs A 的TDI比较
if 'trajectoryb' in adata_mds.obs.columns:
    # 给MDS元细胞附加trajectoryb信息
    # 需要重新计算：对每个metacell，看有多少原始细胞是trajectoryb=True
    print("\n  Analyzing Trajectory B effect in MDS...")
    
    # 简化：直接用细胞级别的TDI（如果内存允许）
    # 或者用元细胞的加权平均
    
    # 这里我们在原始adata_mds上计算cell-level TDI
    print("  Computing cell-level TDI for MDS cells...")
    analyzer.compute_tdi(adata_mds, metric=OT_METRIC, reg=OT_REG)
    tdi_mds_cells = analyzer.tdi_results.copy()
    tdi_mds_cells['trajectoryb'] = adata_mds.obs['trajectoryb'].values
    tdi_mds_cells['hap'] = adata_mds.obs['hap'].values
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    sns.boxplot(data=tdi_mds_cells, x='trajectoryb', y='TDI', palette=['#3498db', '#e74c3c'], ax=axes[0])
    axes[0].set_title("MDS: TDI by Trajectory (A vs B)")
    axes[0].set_xlabel("Trajectory B (Ineffective)")
    
    sns.boxplot(data=tdi_mds_cells, x='hap', y='TDI', palette=['#3498db', '#e74c3c'], ax=axes[1])
    axes[1].set_title("MDS: TDI by 5q- Status")
    axes[1].set_xlabel("5q- Haploinsufficiency")
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "04_mds_trajectory_hap.pdf"))
    plt.savefig(os.path.join(FIG_DIR, "04_mds_trajectory_hap.png"), dpi=300)
    plt.close()
    print("  Saved: 04_mds_trajectory_hap")
    
    # 统计检验
    from scipy import stats
    tdi_b_false = tdi_mds_cells[tdi_mds_cells['trajectoryb'] == False]['TDI']
    tdi_b_true = tdi_mds_cells[tdi_mds_cells['trajectoryb'] == True]['TDI']
    t_stat, p_val = stats.mannwhitneyu(tdi_b_false, tdi_b_true, alternative='two-sided')
    print(f"\n  Trajectory A vs B TDI: p={p_val:.2e} (Mann-Whitney U)")
    print(f"    Trajectory A (n={len(tdi_b_false)}): median={tdi_b_false.median():.4f}, mean={tdi_b_false.mean():.4f}")
    print(f"    Trajectory B (n={len(tdi_b_true)}): median={tdi_b_true.median():.4f}, mean={tdi_b_true.mean():.4f}")

# 5. 模块策略图
fig, ax = plt.subplots(figsize=(8, 4))
strategy_data = analyzer.module_strategy.copy()
colors = ['#e74c3c' if s == 'Up-regulate' else '#3498db' if s == 'Down-regulate' else '#95a5a6'
          for s in strategy_data['Strategy']]
ax.barh(strategy_data['Module'], strategy_data['Repair_Score'], color=colors)
ax.axvline(0, color='black', linewidth=0.5)
ax.set_xlabel("Repair Score")
ax.set_title("Module Repair Strategy (MDS-5q)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "05_module_strategy.pdf"))
plt.savefig(os.path.join(FIG_DIR, "05_module_strategy.png"), dpi=300)
plt.close()
print("  Saved: 05_module_strategy")

# ========================== 步骤12: 保存结果 ==========================
print("\n" + "="*60)
print("Step 12: Saving results")
print("="*60)

results = {
    "cost_map": analyzer.cost_map,
    "tdi_all": tdi_all,
    "repair": analyzer.repair_results,
    "module_strategy": analyzer.module_strategy,
}
save_results(results, OUT_DIR)

# 保存TDI详细结果
tdi_all.to_csv(os.path.join(OUT_DIR, "tdi_all_conditions.csv"), index=False)
if 'tdi_mds_cells' in dir():
    tdi_mds_cells.to_csv(os.path.join(OUT_DIR, "tdi_mds_cell_level.csv"), index=False)

# ========================== 步骤13: 摘要 ==========================
print("\n" + "="*60)
print("RESULTS SUMMARY")
print("="*60)

print(f"\n1. OT Cost Map (Normal Metacells):")
for (s1, s2), d in analyzer.cost_map.items():
    print(f"   {s1} -> {s2}: {d:.4f}")

print(f"\n2. TDI by Condition:")
for cond in ['Normal', 'MDS-5q', 'DBA']:
    subset = tdi_all[tdi_all['condition'] == cond]
    print(f"   {cond}: n={len(subset)}, median={subset['TDI'].median():.4f}, mean={subset['TDI'].mean():.4f}")

print(f"\n3. Best Match Stage (MDS-5q):")
mds_stages = tdi_all[tdi_all['condition'] == 'MDS-5q']['Best_Match_Stage'].value_counts(normalize=True)
for stage, prop in mds_stages.head(5).items():
    print(f"   {stage}: {prop*100:.1f}%")

print(f"\n4. Repair: W={analyzer.repair_results['wasserstein_distance']:.4f}")
print(f"   Top Up: {analyzer.repair_results['top_up_targets'][:5]}")
print(f"   Top Down: {analyzer.repair_results['top_down_targets'][:5]}")

print(f"\n5. Module Strategy:")
for _, row in analyzer.module_strategy.iterrows():
    print(f"   {row['Module']}: {row['Strategy']} ({row['Repair_Score']:.4f})")

print("\n" + "="*60)
print(f"Done! Results -> {OUT_DIR}")
print(f"Figures -> {FIG_DIR}")
print("="*60)
