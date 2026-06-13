#!/usr/bin/env python3
"""
Comprehensive Benchmark: scTDRP vs Baselines and State-of-the-Art Methods

This script systematically compares scTDRP against:
  1. Naive baselines: PCA Euclidean distance, Gene-set score, Pseudotime deviation
  2. Ablation: scTDRP without OT (nearest centroid Euclidean)
  3. DGE baseline: Differential expression feature score
  4. Waddington-OT (WOT): Adapted for disease-normal mapping
  5. MOSCOT (TemporalProblem): Adapted for disease-normal alignment

Evaluation metrics:
  - AUC: Malignant vs Normal Residual
  - AUC: High-CNV vs Low-CNV (malignant only)
  - Spearman r: Cell cycle score
  - Spearman r: CNV score
  - Stage concordance: Best-match stage agreement with expected biology
  - Computation time

Output:
  - results_benchmark/method_comparison_metrics.csv
  - results_benchmark/method_comparison_scores.csv
  - figures_benchmark/Figure_Comparison_*.pdf/png
"""

import os
import sys
import time
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import cdist
from sklearn.metrics import roc_auc_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import anndata

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# ============================================================================
# Configuration
# ============================================================================
AEL_PATH = "/data1/yja/zhongzhuan/4.results/4.infercnv/infercnv_aml5/aml5_annotated_with_cnv_cycle.h5ad"
NORMAL_PATH = "/data1/yja/zhongzhuan/1.data/processed/erythroid_lineage_from_MEP.h5ad"
MODULES_PATH = "/data1/yja/zhongzhuan/5.external/scTDRP/modules.json"
OUT_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/results_benchmark"
FIG_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/figures_benchmark"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

STAGE_ORDER = [
    'MEP', 'BFU-E', 'CFU-E', 'Pro-Erythroblast',
    'Basophilic Erythroblast', 'Polychromatic Erythroblast',
    'Orthochromatic Erythroblast'
]
TERMINAL_STAGE = 'Orthochromatic Erythroblast'

np.random.seed(42)

# ============================================================================
# Optional imports
# ============================================================================
try:
    import wot
    WOT_AVAILABLE = True
except Exception as e:
    WOT_AVAILABLE = False
    print(f"[INFO] Waddington-OT not available: {e}")

try:
    from moscot.problems.time import TemporalProblem
    MOSCOT_AVAILABLE = True
except Exception as e:
    MOSCOT_AVAILABLE = False
    print(f"[INFO] MOSCOT not available: {e}")

try:
    from scTDRP import TDRPAnalyzer
    SCTDRP_AVAILABLE = True
except Exception as e:
    SCTDRP_AVAILABLE = False
    print(f"[WARNING] scTDRP not available: {e}")
    sys.path.insert(0, "/data1/yja/zhongzhuan/5.external/scTDRP/src")
    from scTDRP import TDRPAnalyzer


# ============================================================================
# Data Loading & Preprocessing
# ============================================================================
def load_and_preprocess():
    """Load and preprocess normal and disease data with unified gene space and PCA."""
    print("=" * 70)
    print("[1] Loading data...")
    print("=" * 70)
    adata_ael = sc.read_h5ad(AEL_PATH)
    adata_normal = sc.read_h5ad(NORMAL_PATH)
    print(f"  AEL raw: {adata_ael.shape}")
    print(f"  Normal raw: {adata_normal.shape}")

    # Save normal stage labels
    normal_stage_labels = adata_normal.obs['AuthorCellType'].copy()

    # Gene mapping: normal uses Ensembl -> map to HGNC symbol
    print("\n[2] Gene mapping (Ensembl -> HGNC)...")
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

    # Preprocessing (mirror run_erythroid_scTDRP.py exactly)
    print("\n[3] Preprocessing (normalize -> log1p -> HVG -> scale -> PCA)...")
    adata_ael.X = adata_ael.layers['logcounts']

    adata_combined = sc.concat(
        [adata_normal, adata_ael], label='_source',
        keys=['normal', 'disease'], index_unique='-', join='outer'
    )
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, flavor='seurat_v3')
    sc.pp.scale(adata_combined)
    sc.pp.pca(adata_combined, n_comps=50, use_highly_variable=True)

    mask_normal = adata_combined.obs['_source'] == 'normal'
    mask_disease = adata_combined.obs['_source'] == 'disease'
    normal_pca = adata_combined.obsm['X_pca'][mask_normal]
    disease_pca = adata_combined.obsm['X_pca'][mask_disease]

    # Attach PCA back to individual objects
    adata_normal.obsm['X_pca'] = normal_pca
    adata_ael.obsm['X_pca'] = disease_pca

    # Ensure HVG info on normal data for scTDRP
    adata_normal.var['highly_variable'] = adata_combined.var['highly_variable'].values

    print(f"  PCA: normal={normal_pca.shape}, disease={disease_pca.shape}")

    # Subset malignant cells for disease analysis
    mask_malignant = adata_ael.obs['malignancy'] == 'Malignant Erythroid'
    adata_malignant = adata_ael[mask_malignant].copy()
    print(f"  Malignant cells: {adata_malignant.shape[0]}")

    # Prepare normal data with proper stage key
    if 'AuthorCellType' in adata_normal.obs.columns:
        adata_normal.obs['stage'] = adata_normal.obs['AuthorCellType'].values
    else:
        adata_normal.obs['stage'] = normal_stage_labels.reindex(adata_normal.obs_names).values

    # CRITICAL: Re-run joint preprocessing on normal + malignant ONLY
    # This ensures PCA/HVG are computed on the exact same cells as run_erythroid_scTDRP.py
    print("\n[3b] Re-processing normal + malignant for scTDRP consistency...")
    adata_joint = sc.concat(
        [adata_normal, adata_malignant], label='_source',
        keys=['normal', 'disease'], index_unique='-', join='outer'
    )
    sc.pp.normalize_total(adata_joint, target_sum=1e4)
    sc.pp.log1p(adata_joint)
    sc.pp.highly_variable_genes(adata_joint, n_top_genes=2000, flavor='seurat_v3')
    sc.pp.scale(adata_joint)
    sc.pp.pca(adata_joint, n_comps=50, use_highly_variable=True)

    jmask_normal = adata_joint.obs['_source'] == 'normal'
    jmask_disease = adata_joint.obs['_source'] == 'disease'
    adata_normal.obsm['X_pca'] = adata_joint.obsm['X_pca'][jmask_normal]
    adata_malignant.obsm['X_pca'] = adata_joint.obsm['X_pca'][jmask_disease]
    adata_normal.var['highly_variable'] = adata_joint.var['highly_variable'].values
    adata_malignant.var['highly_variable'] = adata_joint.var['highly_variable'].values

    return adata_normal, adata_ael, adata_malignant, normal_stage_labels


# ============================================================================
# Metacell builder (mirrors run_erythroid_scTDRP.py)
# ============================================================================
def build_metacells(adata, use_rep="X_pca", resolution=1.0, label=None):
    """Aggregate cells into metacells via Leiden clustering."""
    ad = adata.copy()
    rep_key = use_rep if use_rep in ad.obsm else "X_pca"
    sc.pp.neighbors(ad, use_rep=rep_key, n_neighbors=15)
    sc.tl.leiden(ad, resolution=resolution)
    n_meta = ad.obs['leiden'].nunique()
    print(f"    {label}: {ad.shape[0]} cells -> {n_meta} metacells")
    clusters = ad.obs['leiden'].unique()
    X_list = []
    obs_list = []
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
    meta_ad.obs_names = meta_obs['metacell_id'].values.astype(str)
    if rep_key in ad.obsm:
        meta_rep = []
        for cl in clusters:
            mask = ad.obs['leiden'] == cl
            meta_rep.append(ad.obsm[rep_key][mask].mean(axis=0))
        meta_ad.obsm[rep_key] = np.vstack(meta_rep)
    # Carry over cell-to-metacell mapping
    meta_ad.uns['cell_map'] = {cell: f"{label}_MC{cl}" for cell, cl in ad.obs['leiden'].items()}
    return meta_ad


# ============================================================================
# Method 1: scTDRP (Full)
# ============================================================================
def run_scTDRP_full(adata_meta_normal, meta_disease, adata_malignant, adata_normal):
    """Run full scTDRP pipeline on pre-built metacells."""
    print("\n" + "=" * 70)
    print("[Method 1] scTDRP (Full with Metacells)")
    print("=" * 70)
    t0 = time.time()

    # Run scTDRP on metacells
    analyzer = TDRPAnalyzer(
        normal_adata=adata_meta_normal,
        stage_key='stage',
        terminal_stage=TERMINAL_STAGE,
        stage_order=STAGE_ORDER,
        use_rep='X_pca',
        n_top_genes=2000,
    )
    analyzer.prepare_data(flavor='seurat_v3')
    analyzer.build_ot_cost_map(metric='sqeuclidean', reg=0.01)
    tdi_df = analyzer.compute_tdi(meta_disease, metric='sqeuclidean', reg=0.01)

    print(f"  Metacell stage distribution:")
    for stage, count in tdi_df['Best_Match_Stage'].value_counts().items():
        print(f"    {stage}: {count} ({count/len(tdi_df)*100:.1f}%)")

    # Map metacell TDI back to original cells using saved mapping
    cell_map = meta_disease.uns.get('cell_map', {})
    meta_tdi = dict(zip(meta_disease.obs_names, tdi_df['TDI'].values))
    meta_stage = dict(zip(meta_disease.obs_names, tdi_df['Best_Match_Stage'].values))

    cell_tdi = []
    cell_stage = []
    for cell in adata_malignant.obs_names:
        meta_name = cell_map.get(cell)
        if meta_name is None:
            cell_tdi.append(np.nan)
            cell_stage.append('Unknown')
        else:
            cell_tdi.append(meta_tdi.get(meta_name, np.nan))
            cell_stage.append(meta_stage.get(meta_name, 'Unknown'))

    # For repair pathway, use module genes on metacells
    with open(MODULES_PATH) as f:
        modules = json.load(f)
    module_genes = set()
    for glist in modules.values():
        module_genes.update(glist)
    module_genes = [g for g in module_genes if g in adata_normal.var_names]
    print(f"  Module genes for repair: {len(module_genes)}")

    analyzer.gene_list = module_genes
    analyzer.use_rep = 'X'
    analyzer.infer_repair_pathway(meta_disease, metric='sqeuclidean', reg=0.01, top_n=100)

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")

    scores = pd.DataFrame({
        'cell_idx': range(adata_malignant.shape[0]),
        'scTDRP_TDI': cell_tdi,
        'scTDRP_BestStage': cell_stage,
    })
    return scores, elapsed


# ============================================================================
# Method 2: scTDRP-NoOT (Ablation)
# ============================================================================
def run_scTDRP_noOT(adata_meta_normal, meta_disease, adata_malignant):
    """Ablation: replace Wasserstein with Euclidean distance to stage centroids on metacells."""
    print("\n" + "=" * 70)
    print("[Method 2] scTDRP-NoOT (Euclidean on Metacells)")
    print("=" * 70)
    t0 = time.time()

    # Compute stage centroids from normal metacells
    stage_centroids = {}
    for stage in STAGE_ORDER:
        mask = adata_meta_normal.obs['stage'] == stage
        if mask.sum() > 0:
            stage_centroids[stage] = adata_meta_normal.obsm['X_pca'][mask].mean(axis=0)

    centroid_matrix = np.array(list(stage_centroids.values()))
    stage_names = list(stage_centroids.keys())

    # Distance from each disease metacell to each centroid
    distances = cdist(meta_disease.obsm['X_pca'], centroid_matrix, metric='euclidean')
    tdi_noOT = distances.min(axis=1)
    best_stage_idx = distances.argmin(axis=1)
    best_stages = [stage_names[i] for i in best_stage_idx]

    print(f"  NoOT metacell stage distribution:")
    from collections import Counter
    for stage, count in Counter(best_stages).most_common():
        print(f"    {stage}: {count} ({count/len(best_stages)*100:.1f}%)")

    # Map back to original cells using saved mapping
    cell_map = meta_disease.uns.get('cell_map', {})
    meta_tdi = dict(zip(meta_disease.obs_names, tdi_noOT))
    meta_stage = dict(zip(meta_disease.obs_names, best_stages))

    cell_tdi = []
    cell_stage = []
    for cell in adata_malignant.obs_names:
        meta_name = cell_map.get(cell)
        if meta_name is None:
            cell_tdi.append(np.nan)
            cell_stage.append('Unknown')
        else:
            cell_tdi.append(meta_tdi.get(meta_name, np.nan))
            cell_stage.append(meta_stage.get(meta_name, 'Unknown'))

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")

    scores = pd.DataFrame({
        'cell_idx': range(adata_malignant.shape[0]),
        'NoOT_TDI': cell_tdi,
        'NoOT_BestStage': cell_stage,
    })
    return scores, elapsed


# ============================================================================
# Method 3: Naive Baselines
# ============================================================================
def run_naive_baselines(adata_normal, adata_ael, adata_malignant, normal_stage_labels):
    """Run naive baseline methods."""
    print("\n" + "=" * 70)
    print("[Method 3] Naive Baselines")
    print("=" * 70)
    t0_total = time.time()

    normal_pca = adata_normal.obsm['X_pca']
    disease_pca = adata_malignant.obsm['X_pca']

    # --- 3a: PCA Euclidean distance to nearest stage centroid ---
    t0 = time.time()
    stage_centroids = {}
    for stage in STAGE_ORDER:
        mask = normal_stage_labels.values == stage
        if mask.sum() > 0:
            stage_centroids[stage] = normal_pca[mask].mean(axis=0)
    centroid_matrix = np.array(list(stage_centroids.values()))
    distances = cdist(disease_pca, centroid_matrix, metric='euclidean')
    pca_euc_dist = distances.min(axis=1)
    t_euc = time.time() - t0
    print(f"  3a PCA_EucDist: {t_euc:.2f}s")

    # --- 3b: Gene-set differentiation score ---
    t0 = time.time()
    with open(MODULES_PATH) as f:
        modules = json.load(f)
    p7_genes = [g for g in modules['P7_TerminalPrep'] if g in adata_ael.var_names]
    p8_genes = [g for g in modules['P8_ExecutionPrecursor'] if g in adata_ael.var_names]

    sc.tl.score_genes(adata_malignant, gene_list=p7_genes, score_name='P7_score')
    sc.tl.score_genes(adata_malignant, gene_list=p8_genes, score_name='P8_score')

    def minmax_scale(x):
        return (x - x.min()) / (x.max() - x.min() + 1e-10)

    p7 = adata_malignant.obs['P7_score'].values
    p8 = adata_malignant.obs['P8_score'].values
    geneset_diff_score = minmax_scale(minmax_scale(p7) - minmax_scale(p8))
    t_geneset = time.time() - t0
    print(f"  3b GeneSet_DiffScore: {t_geneset:.2f}s")

    # --- 3c: Pseudotime deviation ---
    t0 = time.time()
    normal_pc1 = normal_pca[:, 0]
    stage_pc1_means = []
    for stage in STAGE_ORDER:
        mask = normal_stage_labels.values == stage
        stage_pc1_means.append(normal_pc1[mask].mean() if mask.sum() > 0 else np.nan)

    stage_corr = stats.spearmanr(range(len(STAGE_ORDER)), stage_pc1_means)[0]
    print(f"  PC1 vs stage order: r={stage_corr:.3f}")

    if stage_corr < 0:
        ael_pc1 = -disease_pca[:, 0]
        normal_pc1 = -normal_pc1
    else:
        ael_pc1 = disease_pca[:, 0]

    pt_min, pt_max = normal_pc1.min(), normal_pc1.max()
    ael_pseudotime = np.clip((ael_pc1 - pt_min) / (pt_max - pt_min + 1e-10), 0, 1)
    pseudotime_dev = 1.0 - ael_pseudotime
    t_pseudo = time.time() - t0
    print(f"  3c Pseudotime_Dev: {t_pseudo:.2f}s")

    elapsed_total = time.time() - t0_total

    scores = pd.DataFrame({
        'cell_idx': range(adata_malignant.shape[0]),
        'PCA_EucDist': pca_euc_dist,
        'GeneSet_DiffScore': geneset_diff_score,
        'Pseudotime_Dev': pseudotime_dev,
    })
    return scores, elapsed_total


# ============================================================================
# Method 4: DGE + Module Score Baseline
# ============================================================================
def run_dge_baseline(adata_normal, adata_ael, adata_malignant):
    """
    Differential expression baseline:
    1. Compute malignant vs normal DE genes
    2. For each malignant cell, compute a 'disease feature score'
       = mean(up-regulated genes) - mean(down-regulated genes)
    This mimics what a typical DGE+pathway analysis would output
    if converted to a per-cell metric.
    """
    print("\n" + "=" * 70)
    print("[Method 4] DGE + Disease Feature Score")
    print("=" * 70)
    t0 = time.time()

    # Use the full AEL (malignant + normal residual) vs normal reference
    # For simplicity, compare malignant vs normal reference
    adata_de = sc.concat(
        [adata_normal, adata_malignant], label='group',
        keys=['normal', 'malignant'], index_unique='-'
    )
    adata_de.obs['group'] = adata_de.obs['group'].astype('category')

    sc.tl.rank_genes_groups(
        adata_de, groupby='group', groups=['malignant'],
        reference='normal', method='wilcoxon', n_genes=200
    )

    # Extract top up/down genes
    up_genes = []
    down_genes = []
    try:
        de_df = sc.get.rank_genes_groups_df(adata_de, group='malignant')
        up_genes = de_df[(de_df['logfoldchanges'] > 0) & (de_df['pvals_adj'] < 0.05)]['names'].head(200).tolist()
        down_genes = de_df[(de_df['logfoldchanges'] < 0) & (de_df['pvals_adj'] < 0.05)]['names'].head(200).tolist()
    except Exception as e:
        print(f"  Warning: DE extraction failed ({e}), using fallback")
        up_genes = de_df.nlargest(200, 'logfoldchanges')['names'].tolist()
        down_genes = de_df.nsmallest(200, 'logfoldchanges')['names'].tolist()

    print(f"  DE up genes: {len(up_genes)}, down genes: {len(down_genes)}")

    # Compute per-cell disease feature score
    X_mal = adata_malignant.X
    if hasattr(X_mal, 'toarray'):
        X_mal = X_mal.toarray()

    gene_names = adata_malignant.var_names.tolist()
    up_idx = [gene_names.index(g) for g in up_genes if g in gene_names]
    down_idx = [gene_names.index(g) for g in down_genes if g in gene_names]

    up_mean = X_mal[:, up_idx].mean(axis=1) if len(up_idx) > 0 else np.zeros(X_mal.shape[0])
    down_mean = X_mal[:, down_idx].mean(axis=1) if len(down_idx) > 0 else np.zeros(X_mal.shape[0])

    disease_score = np.asarray(up_mean).flatten() - np.asarray(down_mean).flatten()
    # Normalize to [0,1] for fair AUC comparison
    disease_score = (disease_score - disease_score.min()) / (disease_score.max() - disease_score.min() + 1e-10)

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.1f}s")

    scores = pd.DataFrame({
        'cell_idx': range(adata_malignant.shape[0]),
        'DGE_FeatureScore': disease_score,
    })
    return scores, elapsed


# ============================================================================
# Method 5: Waddington-OT (Adapted)
# ============================================================================
def run_wot_baseline(adata_meta_normal, meta_disease, adata_malignant):
    """
    Adapt Waddington-OT for disease-normal comparison on metacells.
    """
    print("\n" + "=" * 70)
    print("[Method 5] Waddington-OT (Adapted, Metacells)")
    print("=" * 70)

    if not WOT_AVAILABLE:
        print("  [SKIP] Waddington-OT not available")
        return None, 0.0

    t0 = time.time()

    adata_normal_copy = adata_meta_normal.copy()
    adata_mal_copy = meta_disease.copy()

    day_map = {stage: float(i) for i, stage in enumerate(STAGE_ORDER)}
    adata_normal_copy.obs['day'] = adata_normal_copy.obs['stage'].map(day_map).astype(float)
    adata_mal_copy.obs['day'] = float(len(STAGE_ORDER))

    adata_wot = sc.concat([adata_normal_copy, adata_mal_copy], label='source',
                          keys=['normal', 'disease'], index_unique='-')

    print(f"  WOT metacell data shape: {adata_wot.shape}")

    try:
        ot_model = wot.ot.OTModel(adata_wot, day_field='day',
                                   epsilon=0.05, lambda1=1, lambda2=50)
        t_normal_max = float(len(STAGE_ORDER) - 1)
        t_disease = float(len(STAGE_ORDER))
        print(f"  Computing transport map: day {t_normal_max:.0f} -> day {t_disease:.0f}")
        tmap = ot_model.compute_transport_map(t_normal_max, t_disease)

        if hasattr(tmap, 'X'):
            tmap_matrix = tmap.X
            if hasattr(tmap_matrix, 'toarray'):
                tmap_matrix = tmap_matrix.toarray()
            tmap_matrix = np.asarray(tmap_matrix)
            tmap_rows = tmap.obs_names.tolist()
            tmap_cols = tmap.var_names.tolist()
        else:
            tmap_matrix = np.asarray(tmap)
            tmap_rows = adata_normal_copy[adata_normal_copy.obs['day'] == t_normal_max].obs_names.tolist()
            tmap_cols = adata_mal_copy.obs_names.tolist()

        disease_cell_names = adata_mal_copy.obs_names.tolist()
        normal_terminal_names = adata_normal_copy[adata_normal_copy.obs['day'] == t_normal_max].obs_names.tolist()

        row_idx = [tmap_rows.index(n) for n in normal_terminal_names if n in tmap_rows]
        col_idx = [tmap_cols.index(n) for n in disease_cell_names if n in tmap_cols]

        if len(row_idx) == 0 or len(col_idx) == 0:
            print("  [WARNING] WOT tmap name mismatch, falling back to all cells")
            col_sums = tmap_matrix.sum(axis=0)
        else:
            tmap_sub = tmap_matrix[np.ix_(row_idx, col_idx)]
            col_sums = tmap_sub.sum(axis=0)

        col_sums = np.asarray(col_sums).flatten()
        wot_dev = 1.0 - (col_sums - col_sums.min()) / (col_sums.max() - col_sums.min() + 1e-10)

        # Map back to original cells
        disease_cell_map = {}
        sc.pp.neighbors(adata_malignant, use_rep='X_pca', n_neighbors=15)
        sc.tl.leiden(adata_malignant, resolution=1.0)
        for cell, cluster in adata_malignant.obs['leiden'].items():
            disease_cell_map[cell] = int(cluster)

        meta_dev = dict(zip(meta_disease.obs_names, wot_dev))
        cell_dev = []
        for cell in adata_malignant.obs_names:
            cl = str(disease_cell_map[cell])
            meta_name = f"Malignant_MC{cl}"
            cell_dev.append(meta_dev.get(meta_name, np.nan))

        elapsed = time.time() - t0
        print(f"  Elapsed: {elapsed:.1f}s")

        scores = pd.DataFrame({
            'cell_idx': range(adata_malignant.shape[0]),
            'WOT_Deviation': cell_dev,
        })
        return scores, elapsed

    except Exception as e:
        print(f"  [ERROR] WOT failed: {e}")
        import traceback
        traceback.print_exc()
        return None, 0.0


# ============================================================================
# Method 6: MOSCOT TemporalProblem (Adapted)
# ============================================================================
def run_moscot_baseline(adata_meta_normal, meta_disease):
    """
    Adapt MOSCOT TemporalProblem for disease-normal comparison on metacells.
    """
    print("\n" + "=" * 70)
    print("[Method 6] MOSCOT TemporalProblem (Adapted, Metacells)")
    print("=" * 70)

    if not MOSCOT_AVAILABLE:
        print("  [SKIP] MOSCOT not available")
        return None, 0.0

    t0 = time.time()

    adata_normal_copy = adata_meta_normal.copy()
    adata_mal_copy = meta_disease.copy()

    day_map = {stage: i for i, stage in enumerate(STAGE_ORDER)}
    adata_normal_copy.obs['day'] = adata_normal_copy.obs['stage'].map(day_map).astype(int)
    adata_mal_copy.obs['day'] = len(STAGE_ORDER)

    adata_moscot = sc.concat([adata_normal_copy, adata_mal_copy], label='source',
                              keys=['normal', 'disease'], index_unique='-')

    print(f"  MOSCOT metacell data shape: {adata_moscot.shape}")

    try:
        tp = TemporalProblem(adata_moscot)
        # MOSCOT's prepare may internally call PCA with default n_comps=30,
        # which can exceed small time-point sample sizes. Pre-compute PCA
        # and attempt to disable callbacks; if it still fails, gracefully skip.
        n_comps = min(30, max(2, adata_moscot.shape[0] - 1))
        sc.pp.pca(adata_moscot, n_comps=n_comps)
        tp = tp.prepare(time_key="day", xy_callback=None)
        tp = tp.solve(epsilon=0.01)

        t_normal_max = len(STAGE_ORDER) - 1
        t_disease = len(STAGE_ORDER)

        print(f"  Computing coupling: day {t_normal_max} -> day {t_disease}")
        coupling = tp[(t_normal_max, t_disease)].solution.transport_matrix

        col_sums = np.array(coupling.sum(axis=0)).flatten()

        # Deviation score
        moscot_dev = 1.0 - (col_sums - col_sums.min()) / (col_sums.max() - col_sums.min() + 1e-10)

        # Map metacell deviation back to original cells using saved mapping
        cell_map = meta_disease.uns.get('cell_map', {})
        meta_dev_map = dict(zip(meta_disease.obs_names, moscot_dev))

        cell_dev = []
        for cell in adata_malignant.obs_names:
            meta_name = cell_map.get(cell)
            if meta_name is None:
                cell_dev.append(np.nan)
            else:
                cell_dev.append(meta_dev_map.get(meta_name, np.nan))

        elapsed = time.time() - t0
        print(f"  Elapsed: {elapsed:.1f}s")

        scores = pd.DataFrame({
            'cell_idx': range(meta_disease.shape[0]),
            'MOSCOT_Deviation': cell_dev,
        })
        return scores, elapsed

    except Exception as e:
        print(f"  [SKIP] MOSCOT failed on metacell data (sample size too small for internal PCA): {e}")
        return None, 0.0


# ============================================================================
# Evaluation
# ============================================================================
def evaluate_all_methods(comp_df, results_dict):
    """Compute evaluation metrics for all methods."""
    print("\n" + "=" * 70)
    print("[Evaluation] Computing metrics...")
    print("=" * 70)

    malignancy = comp_df['malignancy'].values
    phase = comp_df['phase'].values
    cnv_score = comp_df['cnv_score'].values
    cnv_median = np.median(cnv_score[malignancy == 'Malignant Erythroid'])
    cnv_group = np.array(['High_CNV' if s > cnv_median else 'Low_CNV' for s in cnv_score])
    cycle_score = np.array([{'G1': 0.0, 'S': 1.0, 'G2M': 2.0}.get(p, np.nan) for p in phase])

    methods = []
    for key in results_dict:
        if key.endswith('_scores') and results_dict[key] is not None:
            method_name = key.replace('_scores', '')
            # Only include numeric score columns (exclude stage labels, indices, etc.)
            score_cols = [c for c in results_dict[key].columns
                          if c not in ['cell_idx']
                          and 'Stage' not in c
                          and pd.api.types.is_numeric_dtype(results_dict[key][c])]
            for col in score_cols:
                methods.append((method_name, col))

    records = []

    # AUC: Malignant vs Normal Residual
    mask_malnorm = comp_df['malignancy'].isin(['Malignant Erythroid', 'Normal Residual'])
    y_true_malnorm = (comp_df.loc[mask_malnorm, 'malignancy'] == 'Malignant Erythroid').astype(int)
    print("\n  AUC (Malignant vs Normal Residual):")
    for method_name, col in methods:
        vals = comp_df.loc[mask_malnorm, col].values
        valid = ~pd.isna(vals)
        if valid.sum() < 2 or len(np.unique(y_true_malnorm[valid])) < 2:
            auc = np.nan
        else:
            auc = roc_auc_score(y_true_malnorm[valid], vals[valid])
        records.append({'Metric': 'AUC_Malignant_vs_Normal', 'Method': f"{method_name}:{col}", 'Value': auc})
        print(f"    {method_name}:{col}: {auc:.3f}")

    # AUC: High-CNV vs Low-CNV (malignant only)
    mask_mal = comp_df['malignancy'] == 'Malignant Erythroid'
    cnv_score = comp_df.loc[mask_mal, 'cnv_score'].values
    cnv_median = np.median(cnv_score)
    cnv_group = np.array(['High_CNV' if s > cnv_median else 'Low_CNV' for s in cnv_score])
    y_true_cnv = (cnv_group == 'High_CNV').astype(int)
    print("\n  AUC (High-CNV vs Low-CNV, malignant):")
    for method_name, col in methods:
        vals = comp_df.loc[mask_mal, col].values
        valid = ~pd.isna(vals)
        if valid.sum() < 2 or len(np.unique(y_true_cnv[valid])) < 2:
            auc = np.nan
        else:
            auc = roc_auc_score(y_true_cnv[valid], vals[valid])
        records.append({'Metric': 'AUC_HighCNV_vs_LowCNV', 'Method': f"{method_name}:{col}", 'Value': auc})
        print(f"    {method_name}:{col}: {auc:.3f}")

    # Spearman: Cell Cycle
    mask_valid = mask_mal & (~np.isnan(cycle_score))
    print("\n  Spearman r (Cell Cycle, malignant):")
    for method_name, col in methods:
        vals = comp_df.loc[mask_valid, col].values
        valid = ~pd.isna(vals)
        if valid.sum() < 3:
            r, p = np.nan, np.nan
        else:
            r, p = stats.spearmanr(vals[valid], cycle_score[mask_valid][valid])
        records.append({'Metric': 'Spearman_CycleScore', 'Method': f"{method_name}:{col}", 'Value': r, 'Pvalue': p})
        print(f"    {method_name}:{col}: r={r:+.3f}, p={p:.2e}")

    # Spearman: CNV
    print("\n  Spearman r (CNV Score, malignant):")
    for method_name, col in methods:
        vals = comp_df.loc[mask_mal, col].values
        valid = ~pd.isna(vals)
        if valid.sum() < 3:
            r, p = np.nan, np.nan
        else:
            r, p = stats.spearmanr(vals[valid], cnv_score[mask_mal][valid])
        records.append({'Metric': 'Spearman_CNVScore', 'Method': f"{method_name}:{col}", 'Value': r, 'Pvalue': p})
        print(f"    {method_name}:{col}: r={r:+.3f}, p={p:.2e}")

    # Stage concordance (for methods that output best stage)
    print("\n  Stage concordance with expected (Polychromatic Erythroblast dominant):")
    for key in results_dict:
        if key.endswith('_scores'):
            df = results_dict[key]
            best_stage_col = [c for c in df.columns if 'BestStage' in c or 'Best_Stage' in c]
            if best_stage_col:
                method_name = key.replace('_scores', '')
                stages = df[best_stage_col[0]].values
                # Exclude Unknown/NaN stages
                valid_stages = [s for s in stages if pd.notna(s) and s != 'Unknown']
                if len(valid_stages) > 0:
                    poly_frac = (np.array(valid_stages) == 'Polychromatic Erythroblast').mean()
                else:
                    poly_frac = np.nan
                records.append({'Metric': 'Frac_Polychromatic', 'Method': method_name, 'Value': poly_frac})
                print(f"    {method_name}: {poly_frac:.1%} assigned to Polychromatic Erythroblast (n_valid={len(valid_stages)})")

    results_df = pd.DataFrame(records)
    results_df.to_csv(os.path.join(OUT_DIR, 'method_comparison_metrics.csv'), index=False)
    comp_df.to_csv(os.path.join(OUT_DIR, 'method_comparison_scores.csv'), index=False)

    return results_df


# ============================================================================
# Visualization
# ============================================================================
def plot_results(comp_df, results_df, results_dict):
    """Generate comparison figures."""
    print("\n" + "=" * 70)
    print("[Visualization] Generating figures...")
    print("=" * 70)

    # Simplify method names for plotting
    method_rename = {
        'scTDRP:scTDRP_TDI': 'scTDRP',
        'NoOT:NoOT_TDI': 'scTDRP-NoOT',
        'Naive:PCA_EucDist': 'PCA EucDist',
        'Naive:GeneSet_DiffScore': 'GeneSet Score',
        'Naive:Pseudotime_Dev': 'Pseudotime Dev',
        'DGE:DGE_FeatureScore': 'DGE Feature',
        'WOT:WOT_Deviation': 'Waddington-OT',
        'MOSCOT:MOSCOT_Deviation': 'MOSCOT',
    }
    results_df['Method_Short'] = results_df['Method'].map(method_rename).fillna(results_df['Method'])

    # Color palette
    palette_map = {
        'scTDRP': '#e74c3c',
        'scTDRP-NoOT': '#c0392b',
        'PCA EucDist': '#3498db',
        'GeneSet Score': '#2ecc71',
        'Pseudotime Dev': '#f39c12',
        'DGE Feature': '#9b59b6',
        'Waddington-OT': '#1abc9c',
        'MOSCOT': '#34495e',
    }

    # --- Figure 1: AUC comparison ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    auc_mal = results_df[results_df['Metric'] == 'AUC_Malignant_vs_Normal'].copy()
    if auc_mal['Value'].isna().all():
        axes[0].text(0.5, 0.5, 'N/A\n(requires scores for normal residual cells)',
                     ha='center', va='center', transform=axes[0].transAxes,
                     fontsize=12, style='italic', color='gray')
        axes[0].set_title("AUC: Malignant vs Normal", fontsize=13, fontweight='bold')
        axes[0].set_xticks([])
        axes[0].set_yticks([])
    else:
        order_mal = auc_mal.sort_values('Value', ascending=False)['Method_Short'].tolist()
        sns.barplot(data=auc_mal, x='Method_Short', y='Value',
                    order=order_mal,
                    palette=[palette_map.get(m, '#95a5a6') for m in order_mal],
                    ax=axes[0])
        axes[0].set_ylim(0.5, 1.0)
        axes[0].axhline(0.5, color='gray', linestyle='--', alpha=0.5)
        axes[0].set_title("AUC: Malignant vs Normal", fontsize=13, fontweight='bold')
        axes[0].set_xlabel("")
        axes[0].set_ylabel("AUC")
        for i, row in enumerate(auc_mal.set_index('Method_Short').reindex(order_mal).itertuples()):
            if not pd.isna(row.Value):
                axes[0].text(i, row.Value + 0.01, f"{row.Value:.3f}", ha='center', fontsize=9, fontweight='bold')
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')

    auc_cnv = results_df[results_df['Metric'] == 'AUC_HighCNV_vs_LowCNV'].copy()
    order_cnv = auc_cnv.sort_values('Value', ascending=False)['Method_Short'].tolist()
    sns.barplot(data=auc_cnv, x='Method_Short', y='Value',
                order=order_cnv,
                palette=[palette_map.get(m, '#95a5a6') for m in order_cnv],
                ax=axes[1])
    axes[1].set_ylim(0.5, 1.0)
    axes[1].axhline(0.5, color='gray', linestyle='--', alpha=0.5)
    axes[1].set_title("AUC: High-CNV vs Low-CNV", fontsize=13, fontweight='bold')
    axes[1].set_xlabel("")
    axes[1].set_ylabel("AUC")
    for i, row in enumerate(auc_cnv.set_index('Method_Short').reindex(order_cnv).itertuples()):
        if not pd.isna(row.Value):
            axes[1].text(i, row.Value + 0.01, f"{row.Value:.3f}", ha='center', fontsize=9, fontweight='bold')
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=30, ha='right')
    # Add annotation explaining that TDI measures differentiation deviation, not CNV burden
    axes[1].annotate('Note: TDI measures trajectory deviation,\nnot CNV burden per se',
                     xy=(0.98, 0.02), xycoords='axes fraction',
                     ha='right', va='bottom', fontsize=8, style='italic',
                     color='dimgray',
                     bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_AUC.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_AUC.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIG_DIR}/Figure_Comparison_AUC.pdf")

    # --- Figure 2: Correlation heatmap ---
    corr_data = results_df[results_df['Metric'].isin(['Spearman_CycleScore', 'Spearman_CNVScore'])].copy()
    corr_pivot = corr_data.pivot(index='Method_Short', columns='Metric', values='Value')
    fig, ax = plt.subplots(figsize=(8, max(4, len(corr_pivot) * 0.6)))
    sns.heatmap(corr_pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0,
                vmin=-1, vmax=1, ax=ax, linewidths=0.5)
    ax.set_title("Spearman Correlation with Biological Annotations", fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Correlation.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Correlation.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIG_DIR}/Figure_Comparison_Correlation.pdf")

    # --- Figure 3: Boxplots ---
    mask_malnorm = comp_df['malignancy'].isin(['Malignant Erythroid', 'Normal Residual'])
    score_cols = [c for c in comp_df.columns if c not in [
        'malignancy', 'phase', 'cnv_score', 'cnv_group', 'cycle_score',
        'cell_idx', 'scTDRP_BestStage', 'NoOT_BestStage'
    ]]
    n_cols = len(score_cols)
    fig, axes = plt.subplots(1, n_cols, figsize=(4 * n_cols, 4))
    if n_cols == 1:
        axes = [axes]
    for ax, col in zip(axes, score_cols):
        plot_df = comp_df[mask_malnorm]
        sns.boxplot(data=plot_df, x='malignancy', y=col,
                    palette=['#2ecc71', '#e74c3c'], ax=ax)
        ax.set_title(col, fontsize=11, fontweight='bold')
        unique_labels = plot_df['malignancy'].unique()
        if len(unique_labels) == 1:
            label = unique_labels[0].replace(' Erythroid', '').replace(' Residual', '')
            ax.set_xticklabels([label], rotation=15)
        else:
            ax.set_xticklabels(['Normal', 'Malignant'], rotation=15)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Boxplots.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Boxplots.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIG_DIR}/Figure_Comparison_Boxplots.pdf")

    # --- Figure 4: Stage assignment comparison (for methods with stage output) ---
    stage_cols = [c for c in comp_df.columns if 'BestStage' in c or 'Best_Stage' in c]
    if len(stage_cols) > 0:
        fig, axes = plt.subplots(1, len(stage_cols), figsize=(5 * len(stage_cols), 4))
        if len(stage_cols) == 1:
            axes = [axes]
        for ax, col in zip(axes, stage_cols):
            stage_counts = comp_df.loc[comp_df['malignancy'] == 'Malignant Erythroid', col].value_counts()
            stage_counts = stage_counts.reindex(STAGE_ORDER, fill_value=0)
            bars = stage_counts.plot(kind='barh', ax=ax, color='steelblue')
            ax.set_title(col.replace('_', ' '), fontsize=11, fontweight='bold')
            ax.set_xlabel("Cell Count")
            # Add percentage labels
            total = stage_counts.sum()
            for i, (stage, count) in enumerate(stage_counts.items()):
                if count > 0:
                    pct = count / total * 100
                    ax.text(count + total * 0.01, i, f"{pct:.1f}%",
                            va='center', fontsize=9, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Stages.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Stages.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {FIG_DIR}/Figure_Comparison_Stages.pdf")

    # --- Figure 5: Runtime comparison ---
    runtime_data = []
    for key, val in results_dict.items():
        if key.endswith('_time'):
            method_name = key.replace('_time', '')
            runtime_data.append({'Method': method_rename.get(method_name, method_name), 'Time_s': val})
    if runtime_data:
        rt_df = pd.DataFrame(runtime_data)
        fig, ax = plt.subplots(figsize=(8, 4))
        sns.barplot(data=rt_df, x='Method', y='Time_s', palette='viridis', ax=ax)
        ax.set_title("Computation Time Comparison", fontsize=13, fontweight='bold')
        ax.set_ylabel("Time (seconds)")
        ax.set_xlabel("")
        for i, row in rt_df.iterrows():
            ax.text(i, row['Time_s'] + max(rt_df['Time_s']) * 0.02,
                    f"{row['Time_s']:.1f}s", ha='center', fontsize=9)
        plt.xticks(rotation=30, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Runtime.pdf'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(FIG_DIR, 'Figure_Comparison_Runtime.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {FIG_DIR}/Figure_Comparison_Runtime.pdf")


# ============================================================================
# Main
# ============================================================================
def main():
    print("\n" + "=" * 70)
    print("COMPREHENSIVE BENCHMARK: scTDRP vs State-of-the-Art Methods")
    print("=" * 70)

    # Load data
    adata_normal, adata_ael, adata_malignant, normal_stage_labels = load_and_preprocess()

    # Build metacells for OT-based methods
    print("\n" + "=" * 70)
    print("[Preprocessing] Building metacells for OT-based methods")
    print("=" * 70)
    meta_normal_list = []
    for stage in STAGE_ORDER:
        if stage not in adata_normal.obs['stage'].values:
            continue
        stage_ad = adata_normal[adata_normal.obs['stage'] == stage].copy()
        if stage_ad.shape[0] == 0:
            continue
        res = max(0.5, min(3.0, 30 * 1.0 / max(1, stage_ad.shape[0] / 50)))
        meta_stage = build_metacells(stage_ad, use_rep='X_pca', resolution=res, label=stage)
        meta_stage.obs['stage'] = stage
        meta_normal_list.append(meta_stage)

    adata_meta_normal = anndata.concat(meta_normal_list, label="_stage_batch", index_unique="-")
    adata_meta_normal.obs['stage'] = pd.Categorical(
        [s for ad in meta_normal_list for s in [ad.obs['stage'].iloc[0]] * ad.shape[0]],
        categories=STAGE_ORDER, ordered=True
    )
    print(f"  Total normal metacells: {adata_meta_normal.shape[0]}")

    meta_disease = build_metacells(adata_malignant, use_rep='X_pca', resolution=1.0, label='Malignant')
    print(f"  Total disease metacells: {meta_disease.shape[0]}")

    results_dict = {}

    # Method 1: scTDRP Full
    scores, elapsed = run_scTDRP_full(adata_meta_normal, meta_disease, adata_malignant, adata_normal)
    results_dict['scTDRP_scores'] = scores
    results_dict['scTDRP_time'] = elapsed

    # Method 2: scTDRP-NoOT
    scores, elapsed = run_scTDRP_noOT(adata_meta_normal, meta_disease, adata_malignant)
    results_dict['NoOT_scores'] = scores
    results_dict['NoOT_time'] = elapsed

    # Method 3: Naive baselines (on original cells)
    scores, elapsed = run_naive_baselines(adata_normal, adata_ael, adata_malignant, normal_stage_labels)
    results_dict['Naive_scores'] = scores
    results_dict['Naive_time'] = elapsed

    # Method 4: DGE baseline (on original cells)
    scores, elapsed = run_dge_baseline(adata_normal, adata_ael, adata_malignant)
    results_dict['DGE_scores'] = scores
    results_dict['DGE_time'] = elapsed

    # Method 5: WOT (on metacells)
    scores, elapsed = run_wot_baseline(adata_meta_normal, meta_disease, adata_malignant)
    if scores is not None:
        results_dict['WOT_scores'] = scores
        results_dict['WOT_time'] = elapsed

    # Method 6: MOSCOT (on metacells)
    scores, elapsed = run_moscot_baseline(adata_meta_normal, meta_disease)
    if scores is not None:
        results_dict['MOSCOT_scores'] = scores
        results_dict['MOSCOT_time'] = elapsed

    # Build comparison dataframe
    print("\n[Merge] Building comparison dataframe...")
    comp_df = pd.DataFrame({
        'malignancy': adata_malignant.obs['malignancy'].values,
        'phase': adata_malignant.obs['phase'].values,
        'cnv_score': adata_malignant.obs['cnv_score'].values,
    })

    for key, df in results_dict.items():
        if key.endswith('_scores') and df is not None:
            for col in df.columns:
                if col != 'cell_idx':
                    comp_df[col] = df[col].values

    # Evaluate
    results_df = evaluate_all_methods(comp_df, results_dict)

    # Plot
    plot_results(comp_df, results_df, results_dict)

    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print(f"Results: {OUT_DIR}")
    print(f"Figures: {FIG_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
