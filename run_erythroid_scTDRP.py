#!/usr/bin/env python3
"""
scTDRP 红系脱核分析脚本（元细胞版）

用 Leiden 聚类将大量细胞聚合为元细胞，再跑 OT，速度提升 10-100 倍。
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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from scTDRP import TDRPAnalyzer
from scTDRP.distance import compute_wasserstein_distance
from scTDRP.repair import compute_gene_repair_delta
from scTDRP.module import module_repair_strategy
from scTDRP.io import save_results

# ========================== 配置 ==========================
NORMAL_PATH = "../../1.data/processed/erythroid_lineage_from_MEP.h5ad"
DISEASE_PATH = "../../4.results/4.infercnv/infercnv_aml5/aml5_annotated_with_cnv_cycle.h5ad"
MODULES_PATH = "./modules.json"
OUT_DIR = "./results"
FIG_DIR = "./figures"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

STAGE_KEY = "AuthorCellType"
STAGE_ORDER = ["MEP","BFU-E","CFU-E","Pro-Erythroblast","Basophilic Erythroblast","Polychromatic Erythroblast","Orthochromatic Erythroblast"]
TERMINAL_STAGE = "Orthochromatic Erythroblast"
DISEASE_CONDITION_KEY = "malignancy"
DISEASE_CONDITION_VALUE = "Malignant Erythroid"

OT_METRIC = "sqeuclidean"
OT_REG = 0.01
TOP_N_REPAIR = 100
USE_REP = "X_pca"
N_TOP_GENES = 2000

# 元细胞参数
METACELL_RESOLUTION = 1.0   # Leiden分辨率，越大元细胞越多
TARGET_METACELLS_PER_STAGE = 30  # 每个阶段目标元细胞数

# ========================== 辅助函数：构建元细胞 ==========================
def build_metacells(adata, use_rep="X_pca", resolution=1.0, label=None):
    """
    用 Leiden 聚类聚合元细胞
    
    Returns:
        metacell_adata: AnnData，每个 obs 是一个元细胞
        cluster_map: dict，{metacell_id: [cell1, cell2, ...]}
    """
    ad = adata.copy()
    # 用现有表示做邻居图
    rep_key = use_rep if use_rep in ad.obsm else "X_pca"
    sc.pp.neighbors(ad, use_rep=rep_key, n_neighbors=15)
    sc.tl.leiden(ad, resolution=resolution)
    
    n_meta = ad.obs['leiden'].nunique()
    print(f"    {label}: {ad.shape[0]} cells -> {n_meta} metacells (resolution={resolution})")
    
    # 计算每个 cluster 的平均表达
    clusters = ad.obs['leiden'].unique()
    X_list = []
    obs_list = []
    for cl in clusters:
        mask = ad.obs['leiden'] == cl
        cells = ad[mask]
        # 平均表达
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
    # 保留表示层（用于后续OT）
    if rep_key in ad.obsm:
        meta_rep = []
        for cl in clusters:
            mask = ad.obs['leiden'] == cl
            meta_rep.append(ad.obsm[rep_key][mask].mean(axis=0))
        meta_ad.obsm[rep_key] = np.vstack(meta_rep)
    
    return meta_ad, ad.obs['leiden'].to_dict()

# ========================== 步骤1: 加载数据 ==========================
print("="*60)
print("Step 1: Loading data")
print("="*60)

adata_normal = sc.read_h5ad(NORMAL_PATH)
print(f"  Normal: {adata_normal.shape}")

adata_aml5 = sc.read_h5ad(DISEASE_PATH)
print(f"  AML5: {adata_aml5.shape}")

# ========================== 步骤2: 基因名统一 ==========================
print("\n" + "="*60)
print("Step 2: Mapping gene names to HGNC symbols")
print("="*60)

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
    print(f"  Normal mapped: {adata_normal.shape}")

# ========================== 步骤3: 提取恶性细胞 ==========================
print("\n" + "="*60)
print("Step 3: Extracting malignant erythroid cells")
print("="*60)

if DISEASE_CONDITION_KEY in adata_aml5.obs.columns and DISEASE_CONDITION_VALUE in adata_aml5.obs[DISEASE_CONDITION_KEY].values:
    adata_malignant = adata_aml5[adata_aml5.obs[DISEASE_CONDITION_KEY] == DISEASE_CONDITION_VALUE].copy()
    print(f"  Selected {adata_malignant.shape[0]} malignant cells")
else:
    adata_malignant = adata_aml5.copy()
    print(f"  Using all {adata_malignant.shape[0]} cells")

# ========================== 步骤4: 共同基因 + 联合PCA ==========================
print("\n" + "="*60)
print("Step 4: Joint PCA on common genes")
print("="*60)

with open(MODULES_PATH, "r") as f:
    modules = json.load(f)
module_genes = list(set(modules["P7_TerminalPrep"] + modules["P8_ExecutionPrecursor"]))

common_genes = list(set(adata_normal.var_names) & set(adata_malignant.var_names))
print(f"  Common genes: {len(common_genes)}")
covered_module = [g for g in module_genes if g in common_genes]
print(f"  Module genes covered: {len(covered_module)}")

adata_normal = adata_normal[:, common_genes].copy()
adata_malignant = adata_malignant[:, common_genes].copy()

adata_combined = sc.concat([adata_normal, adata_malignant], label="_source",
                             keys=["normal","disease"], index_unique="-", join="outer")
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)
sc.pp.highly_variable_genes(adata_combined, n_top_genes=N_TOP_GENES, flavor="seurat_v3")
sc.pp.scale(adata_combined)
sc.pp.pca(adata_combined, n_comps=50, use_highly_variable=True)

normal_mask = adata_combined.obs["_source"] == "normal"
adata_normal = adata_combined[normal_mask].copy()
adata_malignant = adata_combined[~normal_mask].copy()

# ========================== 步骤5: 构建元细胞 ==========================
print("\n" + "="*60)
print("Step 5: Building metacells")
print("="*60)

# 正常数据：按阶段分别聚类，保证每个阶段都有元细胞
meta_normal_list = []
for stage in STAGE_ORDER:
    if stage not in adata_normal.obs[STAGE_KEY].values:
        continue
    stage_ad = adata_normal[adata_normal.obs[STAGE_KEY] == stage].copy()
    n_cells = stage_ad.shape[0]
    if n_cells == 0:
        continue
    # 根据细胞数调整分辨率，目标每个阶段 ~30 个元细胞
    res = max(0.5, min(3.0, 30 * METACELL_RESOLUTION / max(1, n_cells / 50)))
    meta_stage, _ = build_metacells(stage_ad, use_rep=USE_REP, resolution=res, label=stage)
    meta_stage.obs[STAGE_KEY] = stage
    meta_normal_list.append(meta_stage)
    print(f"      {stage}: {n_cells} cells -> {meta_stage.shape[0]} metacells")

# 合并所有正常元细胞
import anndata
adata_meta_normal = anndata.concat(meta_normal_list, label="_stage_batch", index_unique="-")
adata_meta_normal.obs[STAGE_KEY] = pd.Categorical(
    [s for ad in meta_normal_list for s in [ad.obs[STAGE_KEY].iloc[0]] * ad.shape[0]],
    categories=STAGE_ORDER,
    ordered=True
)
print(f"  Total normal metacells: {adata_meta_normal.shape[0]}")

# 疾病数据：整体聚类
meta_disease, _ = build_metacells(adata_malignant, use_rep=USE_REP,
                                    resolution=METACELL_RESOLUTION, label="Malignant")
print(f"  Total disease metacells: {meta_disease.shape[0]}")

# ========================== 步骤6: scTDRP 分析 ==========================
print("\n" + "="*60)
print("Step 6: Running scTDRP on metacells")
print("="*60)

analyzer = TDRPAnalyzer(
    normal_adata=adata_meta_normal,
    stage_key=STAGE_KEY,
    terminal_stage=TERMINAL_STAGE,
    stage_order=STAGE_ORDER,
    use_rep=USE_REP,
    n_top_genes=N_TOP_GENES,
)
analyzer.build_ot_cost_map(metric=OT_METRIC, reg=OT_REG)
analyzer.compute_tdi(meta_disease, metric=OT_METRIC, reg=OT_REG)

# ========================== 步骤7: 基因级修复路径 ==========================
print("\n" + "="*60)
print("Step 7: Gene-level repair pathway")
print("="*60)

# 提取模块基因表达（用 log1p 后的数据）
terminal_meta = adata_meta_normal[adata_meta_normal.obs[STAGE_KEY] == TERMINAL_STAGE]

X_disease = meta_disease[:, covered_module].X
X_terminal = terminal_meta[:, covered_module].X
if hasattr(X_disease, "toarray"):
    X_disease = X_disease.toarray()
if hasattr(X_terminal, "toarray"):
    X_terminal = X_terminal.toarray()
X_disease = np.asarray(X_disease, dtype=np.float64)
X_terminal = np.asarray(X_terminal, dtype=np.float64)

# 权重 = 元细胞包含的细胞数
w_disease = meta_disease.obs['n_cells'].values.astype(np.float64)
w_disease = w_disease / w_disease.sum()
w_terminal = terminal_meta.obs['n_cells'].values.astype(np.float64)
w_terminal = w_terminal / w_terminal.sum()

print(f"  Disease metacells: {X_disease.shape}, weights sum={w_disease.sum():.3f}")
print(f"  Terminal metacells: {X_terminal.shape}, weights sum={w_terminal.sum():.3f}")

dist, transport_plan = compute_wasserstein_distance(
    X_disease, X_terminal, w_disease, w_terminal,
    metric=OT_METRIC, reg=OT_REG, method="sinkhorn"
)
print(f"  Wasserstein Distance: {dist:.4f}")

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
print(f"  Top Up: {len(top_up)}, Top Down: {len(top_down)}")

# ========================== 步骤8: 模块策略 ==========================
print("\n" + "="*60)
print("Step 8: Module repair strategy")
print("="*60)

analyzer.module_strategy = module_repair_strategy(
    repair_deltas=repair_deltas, module_dict=modules, threshold=0.0, method="mean"
)
print(analyzer.module_strategy.to_string(index=False))

# ========================== 步骤9: 保存 + 可视化 ==========================
print("\n" + "="*60)
print("Step 9: Saving results and figures")
print("="*60)

results = {
    "cost_map": analyzer.cost_map,
    "tdi": analyzer.tdi_results,
    "repair": analyzer.repair_results,
    "module_strategy": analyzer.module_strategy,
}
save_results(results, OUT_DIR)

analyzer.plot_ot_cost_map(save_path=os.path.join(FIG_DIR, "01_ot_cost_map.pdf"))
analyzer.plot_tdi_distribution(save_path=os.path.join(FIG_DIR, "02_tdi_distribution.pdf"))
analyzer.plot_repair_heatmap(top_n=50, save_path=os.path.join(FIG_DIR, "03_repair_heatmap.pdf"))
analyzer.plot_module_strategy(save_path=os.path.join(FIG_DIR, "04_module_strategy.pdf"))

# ========================== 步骤10: 摘要 ==========================
print("\n" + "="*60)
print("Results Summary")
print("="*60)

print(f"\n1. OT Cost Map (Normal Metacells):")
for (s1, s2), d in analyzer.cost_map.items():
    print(f"   {s1} -> {s2}: {d:.4f}")

print(f"\n2. AML5 Malignant Metacell TDI:")
print(f"   Mean: {results['tdi']['TDI'].mean():.4f}")
print(f"   Median: {results['tdi']['TDI'].median():.4f}")

stage_dist = results['tdi']['Best_Match_Stage'].value_counts(normalize=True)
print(f"\n3. Best Match Stage:")
for stage, prop in stage_dist.head(5).items():
    print(f"   {stage}: {prop*100:.1f}%")

print(f"\n4. Repair: W={results['repair']['wasserstein_distance']:.4f}")
print(f"   Top Up: {results['repair']['top_up_targets'][:5]}")
print(f"   Top Down: {results['repair']['top_down_targets'][:5]}")

print(f"\n5. Module Strategy:")
for _, row in results['module_strategy'].iterrows():
    print(f"   {row['Module']}: {row['Strategy']} ({row['Repair_Score']:.4f})")

print("\n" + "="*60)
print(f"Done! Results -> {OUT_DIR}")
print(f"Figures -> {FIG_DIR}")
print("="*60)
