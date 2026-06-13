#!/usr/bin/env python3
"""Debug: compare metacell data between original and benchmark pipelines."""
import os, sys, json, numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from scTDRP import TDRPAnalyzer

NORMAL_PATH = "/data1/yja/zhongzhuan/1.data/processed/erythroid_lineage_from_MEP.h5ad"
DISEASE_PATH = "/data1/yja/zhongzhuan/4.results/4.infercnv/infercnv_aml5/aml5_annotated_with_cnv_cycle.h5ad"
MODULES_PATH = "/data1/yja/zhongzhuan/5.external/scTDRP/modules.json"
STAGE_ORDER = ["MEP","BFU-E","CFU-E","Pro-Erythroblast","Basophilic Erythroblast","Polychromatic Erythroblast","Orthochromatic Erythroblast"]
TERMINAL_STAGE = "Orthochromatic Erythroblast"
STAGE_KEY = "AuthorCellType"

# ---- Load raw data ----
adata_normal_raw = sc.read_h5ad(NORMAL_PATH)
adata_aml5_raw = sc.read_h5ad(DISEASE_PATH)

# Gene mapping (same as both scripts)
adata_normal_raw.var['gene_symbol'] = adata_normal_raw.var['gene_symbols'].astype(str)
mean_expr = np.array(adata_normal_raw.X.mean(axis=0)).flatten()
if hasattr(mean_expr, 'toarray'): mean_expr = mean_expr.toarray().flatten()
adata_normal_raw.var['_mean_expr'] = mean_expr
keep_idx = []
for sym, group in adata_normal_raw.var.groupby('gene_symbol'):
    if len(group) > 1:
        keep_idx.append(group['_mean_expr'].idxmax())
    else:
        keep_idx.append(group.index[0])
adata_normal_raw = adata_normal_raw[:, keep_idx].copy()
adata_normal_raw.var_names = adata_normal_raw.var['gene_symbol'].values
adata_normal_raw.var_names_make_unique()

common_genes = list(set(adata_aml5_raw.var_names) & set(adata_normal_raw.var_names))
print(f"Common genes: {len(common_genes)}")
adata_aml5_raw = adata_aml5_raw[:, common_genes].copy()
adata_normal_raw = adata_normal_raw[:, common_genes].copy()

# Subset malignant
adata_malignant_raw = adata_aml5_raw[adata_aml5_raw.obs['malignancy'] == 'Malignant Erythroid'].copy()
print(f"Malignant cells: {adata_malignant_raw.shape[0]}")

# ---- Pipeline A: Original (run_erythroid_scTDRP.py) ----
print("\n=== Pipeline A: Original ===")
adata_normal_A = adata_normal_raw.copy()
adata_malignant_A = adata_malignant_raw.copy()
adata_combined_A = sc.concat([adata_normal_A, adata_malignant_A], label="_source",
                               keys=["normal","disease"], index_unique="-", join="outer")
sc.pp.normalize_total(adata_combined_A, target_sum=1e4)
sc.pp.log1p(adata_combined_A)
sc.pp.highly_variable_genes(adata_combined_A, n_top_genes=2000, flavor="seurat_v3")
sc.pp.scale(adata_combined_A)
sc.pp.pca(adata_combined_A, n_comps=50, use_highly_variable=True)

normal_mask = adata_combined_A.obs["_source"] == "normal"
adata_normal_A = adata_combined_A[normal_mask].copy()
adata_malignant_A = adata_combined_A[~normal_mask].copy()

# Build metacells (same function as original)
def build_metacells(adata, use_rep="X_pca", resolution=1.0, label=None):
    ad = adata.copy()
    rep_key = use_rep if use_rep in ad.obsm else "X_pca"
    sc.pp.neighbors(ad, use_rep=rep_key, n_neighbors=15)
    sc.tl.leiden(ad, resolution=resolution)
    clusters = ad.obs['leiden'].unique()
    X_list = []
    obs_list = []
    meta_rep = []
    for cl in clusters:
        mask = ad.obs['leiden'] == cl
        cells = ad[mask]
        X_cl = cells.X.mean(axis=0)
        if hasattr(X_cl, 'toarray'): X_cl = X_cl.toarray().flatten()
        X_list.append(np.asarray(X_cl).flatten())
        obs_list.append({'metacell_id': f"{label}_MC{cl}", 'n_cells': int(mask.sum()), 'leiden': cl})
        if rep_key in ad.obsm:
            meta_rep.append(ad.obsm[rep_key][mask].mean(axis=0))
    meta_X = np.vstack(X_list)
    meta_obs = pd.DataFrame(obs_list)
    meta_var = ad.var.copy()
    meta_ad = sc.AnnData(X=meta_X, obs=meta_obs, var=meta_var)
    meta_ad.obs_names = meta_obs['metacell_id'].values.astype(str)
    if meta_rep:
        meta_ad.obsm[rep_key] = np.vstack(meta_rep)
    return meta_ad

meta_normal_list_A = []
for stage in STAGE_ORDER:
    if stage not in adata_normal_A.obs[STAGE_KEY].values: continue
    stage_ad = adata_normal_A[adata_normal_A.obs[STAGE_KEY] == stage].copy()
    n_cells = stage_ad.shape[0]
    if n_cells == 0: continue
    res = max(0.5, min(3.0, 30 * 1.0 / max(1, n_cells / 50)))
    meta_stage = build_metacells(stage_ad, use_rep="X_pca", resolution=res, label=stage)
    meta_stage.obs[STAGE_KEY] = stage
    meta_normal_list_A.append(meta_stage)

adata_meta_normal_A = sc.concat(meta_normal_list_A, label="_stage_batch", index_unique="-")
adata_meta_normal_A.obs[STAGE_KEY] = pd.Categorical(
    [s for ad in meta_normal_list_A for s in [ad.obs[STAGE_KEY].iloc[0]] * ad.shape[0]],
    categories=STAGE_ORDER, ordered=True
)
meta_disease_A = build_metacells(adata_malignant_A, use_rep="X_pca", resolution=1.0, label="Malignant")
print(f"Normal metacells A: {adata_meta_normal_A.shape[0]}")
print(f"Disease metacells A: {meta_disease_A.shape[0]}")

# Run TDI
analyzer_A = TDRPAnalyzer(
    normal_adata=adata_meta_normal_A,
    stage_key=STAGE_KEY,
    terminal_stage=TERMINAL_STAGE,
    stage_order=STAGE_ORDER,
    use_rep="X_pca",
    n_top_genes=2000,
)
analyzer_A.build_ot_cost_map(metric="sqeuclidean", reg=0.01)
tdi_A = analyzer_A.compute_tdi(meta_disease_A, metric="sqeuclidean", reg=0.01)
print(f"A stage dist: {tdi_A['Best_Match_Stage'].value_counts().to_dict()}")

# ---- Pipeline B: Benchmark (run_comprehensive_benchmark.py) ----
print("\n=== Pipeline B: Benchmark ===")
adata_normal_B = adata_normal_raw.copy()
adata_ael_B = adata_aml5_raw.copy()
adata_ael_B.X = adata_ael_B.layers['logcounts']

adata_combined_B = sc.concat(
    [adata_normal_B, adata_ael_B], label='_source',
    keys=['normal', 'disease'], index_unique='-', join='outer'
)
sc.pp.normalize_total(adata_combined_B, target_sum=1e4)
sc.pp.log1p(adata_combined_B)
sc.pp.highly_variable_genes(adata_combined_B, n_top_genes=2000, flavor='seurat_v3')
sc.pp.scale(adata_combined_B)
sc.pp.pca(adata_combined_B, n_comps=50, use_highly_variable=True)

mask_normal = adata_combined_B.obs['_source'] == 'normal'
mask_disease = adata_combined_B.obs['_source'] == 'disease'
normal_pca = adata_combined_B.obsm['X_pca'][mask_normal]
disease_pca = adata_combined_B.obsm['X_pca'][mask_disease]
adata_normal_B.obsm['X_pca'] = normal_pca
adata_ael_B.obsm['X_pca'] = disease_pca
adata_normal_B.var['highly_variable'] = adata_combined_B.var['highly_variable'].values

if 'AuthorCellType' in adata_normal_B.obs.columns:
    adata_normal_B.obs['stage'] = adata_normal_B.obs['AuthorCellType'].values
else:
    adata_normal_B.obs['stage'] = pd.Series(adata_normal_B.obs['AuthorCellType'].values, index=adata_normal_B.obs_names).reindex(adata_normal_B.obs_names).values

adata_malignant_B = adata_ael_B[adata_ael_B.obs['malignancy'] == 'Malignant Erythroid'].copy()
print(f"Malignant cells B: {adata_malignant_B.shape[0]}")

meta_normal_list_B = []
for stage in STAGE_ORDER:
    if stage not in adata_normal_B.obs['stage'].values: continue
    stage_ad = adata_normal_B[adata_normal_B.obs['stage'] == stage].copy()
    if stage_ad.shape[0] == 0: continue
    res = max(0.5, min(3.0, 30 * 1.0 / max(1, stage_ad.shape[0] / 50)))
    meta_stage = build_metacells(stage_ad, use_rep='X_pca', resolution=res, label=stage)
    meta_stage.obs['stage'] = stage
    meta_normal_list_B.append(meta_stage)

adata_meta_normal_B = sc.concat(meta_normal_list_B, label="_stage_batch", index_unique="-")
adata_meta_normal_B.obs['stage'] = pd.Categorical(
    [s for ad in meta_normal_list_B for s in [ad.obs['stage'].iloc[0]] * ad.shape[0]],
    categories=STAGE_ORDER, ordered=True
)
meta_disease_B = build_metacells(adata_malignant_B, use_rep='X_pca', resolution=1.0, label='Malignant')
print(f"Normal metacells B: {adata_meta_normal_B.shape[0]}")
print(f"Disease metacells B: {meta_disease_B.shape[0]}")

analyzer_B = TDRPAnalyzer(
    normal_adata=adata_meta_normal_B,
    stage_key='stage',
    terminal_stage=TERMINAL_STAGE,
    stage_order=STAGE_ORDER,
    use_rep='X_pca',
    n_top_genes=2000,
)
analyzer_B.prepare_data(flavor='seurat_v3')
analyzer_B.build_ot_cost_map(metric='sqeuclidean', reg=0.01)
tdi_B = analyzer_B.compute_tdi(meta_disease_B, metric='sqeuclidean', reg=0.01)
print(f"B stage dist: {tdi_B['Best_Match_Stage'].value_counts().to_dict()}")

# ---- Compare metacell PCA directly ----
print("\n=== Comparison ===")
print(f"A disease metacell X_pca shape: {meta_disease_A.obsm['X_pca'].shape}")
print(f"B disease metacell X_pca shape: {meta_disease_B.obsm['X_pca'].shape}")
print(f"A normal metacell X_pca shape: {adata_meta_normal_A.obsm['X_pca'].shape}")
print(f"B normal metacell X_pca shape: {adata_meta_normal_B.obsm['X_pca'].shape}")

# Compute PCA difference
if meta_disease_A.obsm['X_pca'].shape == meta_disease_B.obsm['X_pca'].shape:
    diff = np.abs(meta_disease_A.obsm['X_pca'] - meta_disease_B.obsm['X_pca']).mean()
    print(f"Mean abs diff in disease X_pca: {diff:.6f}")
else:
    print("Shapes differ, cannot compare directly")

# Compare normal metacell stage centroids
for stage in ['Polychromatic Erythroblast', 'CFU-E', 'Orthochromatic Erythroblast']:
    mask_A = adata_meta_normal_A.obs[STAGE_KEY] == stage
    mask_B = adata_meta_normal_B.obs['stage'] == stage
    if mask_A.sum() > 0 and mask_B.sum() > 0:
        centroid_A = adata_meta_normal_A.obsm['X_pca'][mask_A].mean(axis=0)
        centroid_B = adata_meta_normal_B.obsm['X_pca'][mask_B].mean(axis=0)
        diff = np.abs(centroid_A - centroid_B).mean()
        print(f"Stage {stage}: centroid diff = {diff:.6f}")
    else:
        print(f"Stage {stage}: A={mask_A.sum()}, B={mask_B.sum()}")

print("\nDone!")
