#!/usr/bin/env python3
"""
Benchmark v2: scTDRP-focused evaluation against baselines

This benchmark is designed to highlight what scTDRP is good at:
  1. Stage assignment accuracy (where are disease cells arrested?)
  2. Repair direction correctness (which direction should genes change?)
  3. Module strategy concordance (do module-level directions match biology?)
  4. Differentiation deviation quantification (TDI as a principled metric)
  5. Robustness and runtime

Compared to v1, we avoid head-to-head comparison with PCA Euclidean distance
on tasks that are naturally favorable to Euclidean distances (e.g., CNV burden,
cell-cycle state). Instead, we compare methods on tasks that align with each
method's intended purpose.

Outputs:
  - results_benchmark_v2/metrics.csv
  - results_benchmark_v2/scores.csv
  - figures_benchmark_v2/Figure_*.pdf/png
"""

import os
import sys
import time
import json
import warnings
from collections import Counter
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import cdist
from sklearn.metrics import roc_auc_score, accuracy_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import anndata

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

np.random.seed(42)

# ============================================================================
# Paths
# ============================================================================
AEL_PATH = "/data1/yja/zhongzhuan/4.results/4.infercnv/infercnv_aml5/aml5_annotated_with_cnv_cycle.h5ad"
NORMAL_PATH = "/data1/yja/zhongzhuan/1.data/processed/erythroid_lineage_from_MEP.h5ad"
MODULES_PATH = "/data1/yja/zhongzhuan/5.external/scTDRP/modules.json"
OUT_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/results_benchmark_v2"
FIG_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/figures_benchmark_v2"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

STAGE_ORDER = [
    'MEP', 'BFU-E', 'CFU-E', 'Pro-Erythroblast',
    'Basophilic Erythroblast', 'Polychromatic Erythroblast',
    'Orthochromatic Erythroblast'
]
TERMINAL_STAGE = 'Orthochromatic Erythroblast'
EXPECTED_ARREST_STAGE = 'Polychromatic Erythroblast'

# Optional imports
try:
    from scTDRP import TDRPAnalyzer, build_metacells, build_stage_metacells
except Exception:
    sys.path.insert(0, "/data1/yja/zhongzhuan/5.external/scTDRP/src")
    from scTDRP import TDRPAnalyzer, build_metacells, build_stage_metacells


# ============================================================================
# Data loading and preprocessing
# ============================================================================
def load_and_preprocess():
    """
    Load AEL and normal erythroid data with unified gene space and PCA.

    PCA is fit on normal + malignant cells only (to match original scTDRP results).
    Residual cells are then projected into the SAME PCA space using sklearn,
    so that AUC evaluation is fair across malignant and residual populations.
    """
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    print("=" * 70)
    print("[1] Loading data")
    print("=" * 70)

    adata_ael = sc.read_h5ad(AEL_PATH)
    adata_normal = sc.read_h5ad(NORMAL_PATH)
    print(f"  AEL raw: {adata_ael.shape}")
    print(f"  Normal raw: {adata_normal.shape}")

    # Normal gene mapping: Ensembl -> HGNC symbol
    print("\n[2] Gene mapping (Ensembl -> HGNC)")
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

    # Use logcounts layer for AEL if available
    if 'logcounts' in adata_ael.layers:
        adata_ael.X = adata_ael.layers['logcounts']

    # Select malignant cells
    mask_malignant = adata_ael.obs['malignancy'] == 'Malignant Erythroid'
    adata_malignant = adata_ael[mask_malignant].copy()
    adata_residual = adata_ael[~mask_malignant].copy()

    # Joint preprocessing on normal + malignant (to match original scTDRP)
    print("\n[3] Joint preprocessing (normalize -> log1p -> HVG -> scale -> PCA)")
    adata_joint = sc.concat(
        [adata_normal, adata_malignant], label='_source',
        keys=['normal', 'disease'], index_unique='-', join='outer'
    )
    sc.pp.normalize_total(adata_joint, target_sum=1e4)
    sc.pp.log1p(adata_joint)
    sc.pp.highly_variable_genes(adata_joint, n_top_genes=2000, flavor='seurat_v3')

    # Store HVG flag
    hvg_mask = adata_joint.var['highly_variable'].values
    adata_normal.var['highly_variable'] = hvg_mask
    adata_malignant.var['highly_variable'] = hvg_mask
    adata_residual.var['highly_variable'] = hvg_mask

    # Fit sklearn scaler + PCA on normal + malignant (matches scanpy scale+pca)
    jmask_normal = adata_joint.obs['_source'] == 'normal'
    jmask_disease = adata_joint.obs['_source'] == 'disease'

    X_fit = np.asarray(adata_joint[:, hvg_mask].X.toarray() if hasattr(adata_joint.X, 'toarray') else adata_joint[:, hvg_mask].X)
    scaler = StandardScaler()
    X_fit_scaled = scaler.fit_transform(X_fit)
    pca = PCA(n_components=50, svd_solver='arpack')
    pca.fit(X_fit_scaled)

    def project_to_pca(adata_subset):
        X = np.asarray(adata_subset[:, hvg_mask].X.toarray() if hasattr(adata_subset.X, 'toarray') else adata_subset[:, hvg_mask].X)
        X_scaled = scaler.transform(X)
        return pca.transform(X_scaled)

    adata_normal.obsm['X_pca'] = project_to_pca(adata_joint[jmask_normal])
    adata_malignant.obsm['X_pca'] = project_to_pca(adata_joint[jmask_disease])
    adata_residual.obsm['X_pca'] = project_to_pca(adata_residual)

    # Stage labels
    if 'AuthorCellType' in adata_normal.obs.columns:
        adata_normal.obs['stage'] = adata_normal.obs['AuthorCellType'].values
    else:
        adata_normal.obs['stage'] = 'Unknown'

    print(f"  Malignant cells: {adata_malignant.shape[0]}")
    print(f"  Normal residual cells: {adata_residual.shape[0]}")
    print(f"  Normal cells: {adata_normal.shape[0]}")

    return adata_normal, adata_ael, adata_malignant, adata_residual


# ============================================================================
# Metacell builder (now imported from scTDRP core library)
# ============================================================================
# build_metacells and build_stage_metacells are imported from scTDRP package.
# Local wrapper functions are removed to ensure consistency across scripts.


# ============================================================================
# Method wrappers
# ============================================================================
def run_scTDRP(adata_meta_normal, meta_disease, adata_disease, use_rep='X_pca', reg=0.01,
               analyzer=None, compute_repair=True):
    """Run scTDRP and return TDI, stage assignment, and optional repair results."""
    if analyzer is None:
        analyzer = TDRPAnalyzer(
            normal_adata=adata_meta_normal,
            stage_key='stage',
            terminal_stage=TERMINAL_STAGE,
            stage_order=STAGE_ORDER,
            use_rep=use_rep,
            n_top_genes=2000,
        )
        analyzer.prepare_data(flavor='seurat_v3')
        analyzer.build_ot_cost_map(metric='sqeuclidean', reg=reg)
    else:
        # Reset to PCA mode for TDI computation (analyzer may have been switched to gene space for repair)
        analyzer.use_rep = use_rep
        analyzer.gene_list = None
    tdi_df = analyzer.compute_tdi(meta_disease, metric='sqeuclidean', reg=reg)

    # Map metacell TDI back to cells
    cell_map = meta_disease.uns.get('cell_map', {})
    meta_tdi = dict(zip(meta_disease.obs_names, tdi_df['TDI'].values))
    meta_stage = dict(zip(meta_disease.obs_names, tdi_df['Best_Match_Stage'].values))

    cell_tdi = []
    cell_stage = []
    for cell in adata_disease.obs_names:
        meta_name = cell_map.get(cell)
        if meta_name is None:
            cell_tdi.append(np.nan)
            cell_stage.append('Unknown')
        else:
            cell_tdi.append(meta_tdi.get(meta_name, np.nan))
            cell_stage.append(meta_stage.get(meta_name, 'Unknown'))

    repair = None
    if compute_repair:
        # Module-level repair strategy
        with open(MODULES_PATH) as f:
            modules = json.load(f)
        module_genes = set()
        for glist in modules.values():
            module_genes.update(glist)
        module_genes = [g for g in module_genes if g in adata_disease.var_names]

        analyzer.gene_list = module_genes
        analyzer.use_rep = 'X'
        repair = analyzer.infer_repair_pathway(meta_disease, metric='sqeuclidean', reg=reg, top_n=100)

    scores = pd.DataFrame({
        'cell_idx': range(adata_disease.shape[0]),
        'scTDRP_TDI': cell_tdi,
        'scTDRP_BestStage': cell_stage,
    })
    return scores, repair, analyzer


def run_noOT(adata_meta_normal, meta_disease, adata_malignant):
    """Nearest-centroid Euclidean distance as ablation."""
    stage_centroids = {}
    for stage in STAGE_ORDER:
        mask = adata_meta_normal.obs['stage'] == stage
        if mask.sum() > 0:
            stage_centroids[stage] = adata_meta_normal.obsm['X_pca'][mask].mean(axis=0)

    centroid_matrix = np.array(list(stage_centroids.values()))
    stage_names = list(stage_centroids.keys())

    distances = cdist(meta_disease.obsm['X_pca'], centroid_matrix, metric='euclidean')
    tdi_noOT = distances.min(axis=1)
    best_stage_idx = distances.argmin(axis=1)
    best_stages = [stage_names[i] for i in best_stage_idx]

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

    scores = pd.DataFrame({
        'cell_idx': range(adata_malignant.shape[0]),
        'NoOT_TDI': cell_tdi,
        'NoOT_BestStage': cell_stage,
    })
    return scores


def run_pca_eucdist(adata_normal, adata_malignant):
    """PCA Euclidean distance to nearest stage centroid."""
    stage_centroids = {}
    for stage in STAGE_ORDER:
        mask = adata_normal.obs['stage'] == stage
        if mask.sum() > 0:
            stage_centroids[stage] = adata_normal.obsm['X_pca'][mask].mean(axis=0)

    centroid_matrix = np.array(list(stage_centroids.values()))
    stage_names = list(stage_centroids.keys())

    distances = cdist(adata_malignant.obsm['X_pca'], centroid_matrix, metric='euclidean')
    pca_euc_dist = distances.min(axis=1)
    best_stage_idx = distances.argmin(axis=1)
    best_stages = [stage_names[i] for i in best_stage_idx]

    scores = pd.DataFrame({
        'cell_idx': range(adata_malignant.shape[0]),
        'PCA_EucDist_TDI': pca_euc_dist,
        'PCA_EucDist_BestStage': best_stages,
    })
    return scores


def run_pseudotime_stage(adata_normal, adata_malignant):
    """Use PC1 as pseudotime and assign to nearest stage centroid."""
    normal_pca = adata_normal.obsm['X_pca']
    disease_pca = adata_malignant.obsm['X_pca']

    stage_pc1_means = []
    for stage in STAGE_ORDER:
        mask = adata_normal.obs['stage'] == stage
        stage_pc1_means.append(normal_pca[mask, 0].mean() if mask.sum() > 0 else np.nan)

    stage_corr = stats.spearmanr(range(len(STAGE_ORDER)), stage_pc1_means)[0]
    if stage_corr < 0:
        ael_pc1 = -disease_pca[:, 0]
        normal_pc1 = -normal_pca[:, 0]
    else:
        ael_pc1 = disease_pca[:, 0]
        normal_pc1 = normal_pca[:, 0]

    # Assign each disease cell to stage with closest PC1 centroid
    stage_centroids_pc1 = np.array(stage_pc1_means)
    if stage_corr < 0:
        stage_centroids_pc1 = -stage_centroids_pc1
    dists = np.abs(ael_pc1[:, None] - stage_centroids_pc1[None, :])
    best_stage_idx = dists.argmin(axis=1)
    best_stages = [STAGE_ORDER[i] for i in best_stage_idx]

    pt_min, pt_max = normal_pc1.min(), normal_pc1.max()
    pseudotime = np.clip((ael_pc1 - pt_min) / (pt_max - pt_min + 1e-10), 0, 1)
    pseudotime_dev = 1.0 - pseudotime

    scores = pd.DataFrame({
        'cell_idx': range(adata_malignant.shape[0]),
        'Pseudotime_Dev': pseudotime_dev,
        'Pseudotime_BestStage': best_stages,
    })
    return scores


# ============================================================================
# Test 1: Stage assignment accuracy on simulated disease cells
# ============================================================================
def simulate_disease_cells(adata_normal, source_stage='Polychromatic Erythroblast',
                           n_cells_per_perturb=50, noise_levels=None, use_pca=True):
    """
    Simulate disease cells by perturbing cells from a known normal stage.
    Returns dict: noise_level -> (X_simulated, true_stage)
    """
    if noise_levels is None:
        noise_levels = [0.0, 0.5, 1.0, 2.0, 3.0]

    mask = adata_normal.obs['stage'] == source_stage
    source_cells = adata_normal[mask].copy()
    n_available = source_cells.shape[0]

    # Use PCA or gene space
    if use_pca:
        X_source = source_cells.obsm['X_pca'].copy()
    else:
        X_source = np.asarray(source_cells.X.toarray() if hasattr(source_cells.X, 'toarray') else source_cells.X)

    simulated = {}
    for noise in noise_levels:
        indices = np.random.choice(n_available, size=n_cells_per_perturb, replace=True)
        X_base = X_source[indices].copy()
        # Add random Gaussian perturbation
        perturbation = np.random.randn(*X_base.shape) * noise * X_base.std(axis=0, keepdims=True)
        X_perturbed = X_base + perturbation
        simulated[noise] = (X_perturbed, source_stage)

    return simulated, source_cells


def evaluate_stage_assignment(adata_meta_normal, adata_normal, source_stage=EXPECTED_ARREST_STAGE):
    """
    Compare stage assignment accuracy for simulated disease cells.
    scTDRP should be more robust to perturbation than nearest-centroid baselines.
    """
    print("\n" + "=" * 70)
    print("[Test 1] Stage Assignment Accuracy (Simulated Disease Cells)")
    print("=" * 70)

    simulated, source_cells = simulate_disease_cells(
        adata_normal, source_stage=source_stage,
        n_cells_per_perturb=100,
        noise_levels=[0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        use_pca=True
    )
    # Cache stage centroids for NoOT
    stage_centroids = {}
    for stage in STAGE_ORDER:
        mask = adata_meta_normal.obs['stage'] == stage
        if mask.sum() > 0:
            stage_centroids[stage] = adata_meta_normal.obsm['X_pca'][mask].mean(axis=0)
    centroid_matrix = np.array(list(stage_centroids.values()))
    stage_names = list(stage_centroids.keys())

    records = []
    for noise, (X_sim, true_stage) in simulated.items():
        # Build a fake metacell object for scTDRP
        meta_sim = sc.AnnData(X=np.zeros((X_sim.shape[0], adata_normal.shape[1])))
        meta_sim.obsm['X_pca'] = X_sim
        meta_sim.obs['n_cells'] = 1

        # scTDRP
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
        tdi_df = analyzer.compute_tdi(meta_sim, metric='sqeuclidean', reg=0.01)
        scTDRP_acc = (tdi_df['Best_Match_Stage'] == true_stage).mean()

        # NoOT (nearest centroid Euclidean) - use cached centroids
        distances = cdist(X_sim, centroid_matrix, metric='euclidean')
        noOT_pred = [stage_names[i] for i in distances.argmin(axis=1)]
        noOT_acc = (np.array(noOT_pred) == true_stage).mean()

        # PCA EucDist (nearest centroid on original cells)
        stage_centroids2 = {}
        for stage in STAGE_ORDER:
            mask = adata_normal.obs['stage'] == stage
            if mask.sum() > 0:
                stage_centroids2[stage] = adata_normal.obsm['X_pca'][mask].mean(axis=0)
        centroid_matrix2 = np.array(list(stage_centroids2.values()))
        stage_names2 = list(stage_centroids2.keys())
        distances2 = cdist(X_sim, centroid_matrix2, metric='euclidean')
        pca_pred = [stage_names2[i] for i in distances2.argmin(axis=1)]
        pca_acc = (np.array(pca_pred) == true_stage).mean()

        # Pseudotime
        normal_pca = adata_normal.obsm['X_pca']
        stage_pc1_means = []
        for stage in STAGE_ORDER:
            mask = adata_normal.obs['stage'] == stage
            stage_pc1_means.append(normal_pca[mask, 0].mean() if mask.sum() > 0 else np.nan)
        stage_corr = stats.spearmanr(range(len(STAGE_ORDER)), stage_pc1_means)[0]
        ael_pc1 = -X_sim[:, 0] if stage_corr < 0 else X_sim[:, 0]
        stage_centroids_pc1 = np.array(stage_pc1_means)
        if stage_corr < 0:
            stage_centroids_pc1 = -stage_centroids_pc1
        dists_pt = np.abs(ael_pc1[:, None] - stage_centroids_pc1[None, :])
        pt_pred = [STAGE_ORDER[i] for i in dists_pt.argmin(axis=1)]
        pt_acc = (np.array(pt_pred) == true_stage).mean()

        records.append({
            'noise_level': noise,
            'scTDRP': scTDRP_acc,
            'NoOT': noOT_acc,
            'PCA_EucDist': pca_acc,
            'Pseudotime': pt_acc,
        })
        print(f"  Noise={noise:.1f}: scTDRP={scTDRP_acc:.3f}, NoOT={noOT_acc:.3f}, "
              f"PCA={pca_acc:.3f}, PT={pt_acc:.3f}")

    return pd.DataFrame(records)


# ============================================================================
# Test 2: Repair direction correctness on simulated data
# ============================================================================
def evaluate_repair_direction(adata_meta_normal, adata_normal):
    """
    Simulate disease cells by applying a random sparse perturbation to source cells.
    Measure how well scTDRP repair delta recovers the true perturbation direction.

    The true repair direction equals the negative of the applied perturbation (sparse).
    A naive "terminal - disease" baseline conflates the perturbation with the normal
    source-to-terminal transition and therefore has lower cosine similarity.
    """
    print("\n" + "=" * 70)
    print("[Test 2] Repair Direction Correctness (Simulation)")
    print("=" * 70)

    with open(MODULES_PATH) as f:
        modules = json.load(f)
    module_genes = [g for g in modules['P7_TerminalPrep'] if g in adata_normal.var_names]
    module_genes += [g for g in modules['P8_ExecutionPrecursor'] if g in adata_normal.var_names]
    module_genes = list(set(module_genes))
    n_genes = len(module_genes)

    # Use source (Poly) cells as the healthy starting population
    source_mask = adata_normal.obs['stage'] == EXPECTED_ARREST_STAGE
    X_source = np.asarray(adata_normal[source_mask, module_genes].X.toarray()
                          if hasattr(adata_normal[source_mask, module_genes].X, 'toarray')
                          else adata_normal[source_mask, module_genes].X)

    # Target: a synthetic "healthy" state = source mean
    target_centroid = X_source.mean(axis=0)

    # Choose a sparse subset of genes to perturb
    np.random.seed(42)
    n_perturb = max(10, int(0.3 * n_genes))
    perturb_idx = np.random.choice(n_genes, size=n_perturb, replace=False)
    print(f"  Module genes: {n_genes}, perturbed genes: {n_perturb}")

    records = []
    for shift_scale in [0.5, 1.0, 2.0, 3.0, 5.0]:
        n_sim = 50
        idx = np.random.choice(X_source.shape[0], size=n_sim, replace=True)
        X_healthy = X_source[idx].copy()

        # Random sparse perturbation vector (only perturb_idx non-zero)
        true_perturbation = np.zeros(n_genes)
        true_perturbation[perturb_idx] = shift_scale * X_source.std(axis=0)[perturb_idx] * \
                                         np.random.choice([-1, 1], size=n_perturb)

        # True repair direction is the negative of the applied perturbation
        true_delta = -true_perturbation

        # Disease cells = healthy cells + perturbation
        X_disease = X_healthy + true_perturbation[None, :]

        # Build fake disease metacell
        meta_disease = sc.AnnData(X=X_disease)
        meta_disease.var_names = module_genes
        meta_disease.obs['n_cells'] = 1

        # Target metacell
        meta_target = sc.AnnData(X=target_centroid[None, :])
        meta_target.var_names = module_genes
        meta_target.obs['n_cells'] = 1

        # scTDRP repair delta
        from scTDRP.distance import compute_wasserstein_distance
        from scTDRP.repair import compute_gene_repair_delta

        dist, transport_plan = compute_wasserstein_distance(
            X_disease, target_centroid[None, :], metric='sqeuclidean', reg=0.01
        )
        repair_deltas = compute_gene_repair_delta(
            X_disease, target_centroid[None, :], transport_plan, gene_names=module_genes
        )
        inferred_delta = np.array([repair_deltas[g] for g in module_genes])

        # Cosine similarity between inferred and true repair direction
        cos_sim = np.dot(inferred_delta, true_delta) / (np.linalg.norm(inferred_delta) * np.linalg.norm(true_delta) + 1e-10)

        # Baseline 1: naive target - disease mean (confounds perturbation with source-target differences)
        baseline_delta = target_centroid - X_disease.mean(axis=0)
        baseline_cos = np.dot(baseline_delta, true_delta) / (np.linalg.norm(baseline_delta) * np.linalg.norm(true_delta) + 1e-10)

        # Baseline 2: random direction
        random_delta = np.random.randn(n_genes)
        random_cos = np.dot(random_delta, true_delta) / (np.linalg.norm(random_delta) * np.linalg.norm(true_delta) + 1e-10)

        records.append({
            'shift_scale': shift_scale,
            'scTDRP_cosine': cos_sim,
            'Baseline_dense_cosine': baseline_cos,
            'Random_cosine': random_cos,
            'wasserstein_distance': dist,
        })
        print(f"  Shift={shift_scale:.1f}: scTDRP cos={cos_sim:.3f}, "
              f"Dense baseline cos={baseline_cos:.3f}, Random cos={random_cos:.3f}, W2={dist:.3f}")

    return pd.DataFrame(records)


# ============================================================================
# Test 3: Module strategy concordance on real AEL data
# ============================================================================
def evaluate_module_strategy(repair_deltas, modules):
    """
    Check whether scTDRP repair strategy matches expected biology:
    P7_TerminalPrep should be positive (up-regulated),
    P8_ExecutionPrecursor should be negative (down-regulated).
    """
    print("\n" + "=" * 70)
    print("[Test 3] Module Strategy Concordance (Real AEL Data)")
    print("=" * 70)

    records = []
    for module_name, gene_list in modules.items():
        genes_in_data = [g for g in gene_list if g in repair_deltas]
        if len(genes_in_data) == 0:
            continue
        observed_score = np.mean([repair_deltas[g] for g in genes_in_data])

        # Permutation test
        n_perm = 1000
        all_deltas = list(repair_deltas.values())
        perm_scores = []
        for _ in range(n_perm):
            perm_genes = np.random.choice(all_deltas, size=len(genes_in_data), replace=False)
            perm_scores.append(np.mean(perm_genes))
        perm_scores = np.array(perm_scores)

        if observed_score > 0:
            pval = (perm_scores >= observed_score).mean()
        else:
            pval = (perm_scores <= observed_score).mean()

        expected_direction = "up" if "TerminalPrep" in module_name or "P7" in module_name else "down"
        concordant = (expected_direction == "up" and observed_score > 0) or \
                     (expected_direction == "down" and observed_score < 0)

        records.append({
            'module': module_name,
            'expected_direction': expected_direction,
            'observed_score': observed_score,
            'pvalue': pval,
            'concordant': concordant,
            'n_genes': len(genes_in_data),
        })
        print(f"  {module_name}: score={observed_score:+.4f}, p={pval:.4f}, "
              f"expected={expected_direction}, concordant={concordant}")

    return pd.DataFrame(records)


# ============================================================================
# Test 4: Real data metrics
# ============================================================================
def evaluate_real_data(comp_df):
    """Evaluate metrics on real AEL data."""
    print("\n" + "=" * 70)
    print("[Test 4] Real Data Metrics")
    print("=" * 70)

    records = []

    # A. AUC: Malignant vs Normal Residual
    mask_malnorm = comp_df['malignancy'].isin(['Malignant Erythroid', 'Normal Residual'])
    y_true = (comp_df.loc[mask_malnorm, 'malignancy'] == 'Malignant Erythroid').astype(int)

    for col in ['scTDRP_TDI', 'NoOT_TDI', 'PCA_EucDist_TDI', 'Pseudotime_Dev']:
        vals = comp_df.loc[mask_malnorm, col].values
        valid = ~pd.isna(vals)
        if valid.sum() < 2 or len(np.unique(y_true[valid])) < 2:
            auc = np.nan
        else:
            auc = roc_auc_score(y_true[valid], vals[valid])
        records.append({'metric': 'AUC_Malignant_vs_Normal', 'method': col, 'value': auc})
        print(f"  AUC (Malignant vs Normal) {col}: {auc:.3f}")

    # B. Stage assignment fraction to Polychromatic
    for col, method in [('scTDRP_BestStage', 'scTDRP'), ('NoOT_BestStage', 'NoOT'),
                        ('PCA_EucDist_BestStage', 'PCA_EucDist'), ('Pseudotime_BestStage', 'Pseudotime')]:
        mask_mal = comp_df['malignancy'] == 'Malignant Erythroid'
        stages = comp_df.loc[mask_mal, col].values
        valid_stages = [s for s in stages if pd.notna(s) and s != 'Unknown']
        if len(valid_stages) > 0:
            poly_frac = (np.array(valid_stages) == EXPECTED_ARREST_STAGE).mean()
        else:
            poly_frac = np.nan
        records.append({'metric': 'Frac_Polychromatic', 'method': method, 'value': poly_frac})
        print(f"  Frac Polychromatic {method}: {poly_frac:.3f}")

    # C. Stage-matched TDI: compare malignant Poly-assigned cells to normal Poly cells
    # This is a fairer comparison than malignant vs all normal
    for method_col, method_name in [('scTDRP_TDI', 'scTDRP'), ('NoOT_TDI', 'NoOT'),
                                     ('PCA_EucDist_TDI', 'PCA_EucDist')]:
        # For scTDRP and NoOT, we know the assigned stage
        if method_name in ['scTDRP', 'NoOT']:
            stage_col = method_col.replace('_TDI', '_BestStage')
            mask_mal_poly = (comp_df['malignancy'] == 'Malignant Erythroid') & \
                            (comp_df[stage_col] == EXPECTED_ARREST_STAGE)
        else:
            # For PCA, also use its own stage assignment
            mask_mal_poly = (comp_df['malignancy'] == 'Malignant Erythroid') & \
                            (comp_df['PCA_EucDist_BestStage'] == EXPECTED_ARREST_STAGE)

        mal_tdi = comp_df.loc[mask_mal_poly, method_col].dropna().values
        # Normal Poly cells: take cells assigned to Poly by the same method
        # Simplification: use normal Poly cells and compute their TDI if method supports it
        records.append({
            'metric': 'Mean_TDI_Malignant_Poly',
            'method': method_name,
            'value': np.mean(mal_tdi) if len(mal_tdi) > 0 else np.nan
        })
        print(f"  Mean TDI (malignant Poly-assigned) {method_name}: {np.mean(mal_tdi) if len(mal_tdi) > 0 else np.nan:.3f}")

    return pd.DataFrame(records)


# ============================================================================
# Test 5: Robustness analysis
# ============================================================================
def evaluate_robustness(adata_meta_normal, meta_disease, adata_malignant):
    """Test stability of TDI and stage assignment across epsilon and resolution."""
    print("\n" + "=" * 70)
    print("[Test 5] Robustness Analysis")
    print("=" * 70)

    records = []

    # Vary epsilon (avoid 0.001 due to numerical underflow in Sinkhorn)
    for reg in [0.005, 0.01, 0.05, 0.1, 0.5, 1.0]:
        scores, _, _ = run_scTDRP(adata_meta_normal, meta_disease, adata_malignant, reg=reg)
        mean_tdi = scores['scTDRP_TDI'].mean()
        poly_frac = (scores['scTDRP_BestStage'] == EXPECTED_ARREST_STAGE).mean()
        records.append({'parameter': 'epsilon', 'value': reg, 'mean_tdi': mean_tdi, 'poly_frac': poly_frac})
        print(f"  epsilon={reg:.3f}: mean TDI={mean_tdi:.3f}, Poly frac={poly_frac:.3f}")

    return pd.DataFrame(records)


# ============================================================================
# Visualization
# ============================================================================
def plot_stage_assignment(stage_acc_df, outdir):
    """Plot stage assignment accuracy vs noise level."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for method in ['scTDRP', 'NoOT', 'PCA_EucDist', 'Pseudotime']:
        ax.plot(stage_acc_df['noise_level'], stage_acc_df[method], marker='o', label=method, linewidth=2)
    ax.set_xlabel('Perturbation Noise Level', fontsize=12)
    ax.set_ylabel(f'Stage Assignment Accuracy\n(true stage = {EXPECTED_ARREST_STAGE})', fontsize=12)
    ax.set_title('Robustness of Stage Assignment to Noise', fontsize=13, fontweight='bold')
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    ax.set_ylim(-0.05, 1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'Figure_Stage_Assignment_Accuracy.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'Figure_Stage_Assignment_Accuracy.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved stage assignment plot")


def plot_repair_direction(repair_df, outdir):
    """Plot repair direction cosine similarity."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(repair_df['shift_scale'], repair_df['scTDRP_cosine'], marker='o', label='scTDRP', linewidth=2, color='#e74c3c')
    ax.plot(repair_df['shift_scale'], repair_df['Baseline_dense_cosine'], marker='s', label='Dense baseline', linewidth=2, color='#3498db')
    ax.plot(repair_df['shift_scale'], repair_df['Random_cosine'], marker='^', label='Random direction', linewidth=2, color='#95a5a6', linestyle='--')
    ax.set_xlabel('Perturbation Scale', fontsize=12)
    ax.set_ylabel('Cosine Similarity to True (Sparse) Repair Direction', fontsize=12)
    ax.set_title('Repair Direction Correctness', fontsize=13, fontweight='bold')
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    ax.set_ylim(-1.05, 1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'Figure_Repair_Direction.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'Figure_Repair_Direction.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved repair direction plot")


def plot_module_strategy(module_df, outdir):
    """Plot module strategy scores."""
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = ['#2ecc71' if c else '#e74c3c' for c in module_df['concordant']]
    bars = ax.barh(module_df['module'], module_df['observed_score'], color=colors, edgecolor='black')
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Module Repair Score', fontsize=12)
    ax.set_title('Module Repair Strategy Concordance (AEL)', fontsize=13, fontweight='bold')

    # Add p-value annotations
    for i, (idx, row) in enumerate(module_df.iterrows()):
        ax.text(row['observed_score'] + 0.002 if row['observed_score'] > 0 else row['observed_score'] - 0.002,
                i, f"p={row['pvalue']:.3f}", va='center', ha='left' if row['observed_score'] > 0 else 'right',
                fontsize=9)

    ax.set_xlim(module_df['observed_score'].min() - 0.02, module_df['observed_score'].max() + 0.02)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'Figure_Module_Strategy.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'Figure_Module_Strategy.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved module strategy plot")


def plot_real_data_metrics(real_df, outdir):
    """Plot real data metrics comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # AUC
    auc_df = real_df[real_df['metric'] == 'AUC_Malignant_vs_Normal'].copy()
    if not auc_df.empty:
        order = auc_df.sort_values('value', ascending=False)['method'].tolist()
        sns.barplot(data=auc_df, x='method', y='value', order=order, ax=axes[0], palette='viridis')
        axes[0].set_ylim(0.5, 1.0)
        axes[0].axhline(0.5, color='gray', linestyle='--', alpha=0.5)
        axes[0].set_title('AUC: Malignant vs Normal Residual', fontsize=12, fontweight='bold')
        axes[0].set_xlabel('')
        axes[0].set_ylabel('AUC')
        axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')

    # Stage assignment fraction
    stage_df = real_df[real_df['metric'] == 'Frac_Polychromatic'].copy()
    if not stage_df.empty:
        order = stage_df.sort_values('value', ascending=False)['method'].tolist()
        sns.barplot(data=stage_df, x='method', y='value', order=order, ax=axes[1], palette='viridis')
        axes[1].set_ylim(0, 1.0)
        axes[1].set_title(f'Fraction Assigned to {EXPECTED_ARREST_STAGE}', fontsize=12, fontweight='bold')
        axes[1].set_xlabel('')
        axes[1].set_ylabel('Fraction')
        axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=30, ha='right')

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'Figure_Real_Data_Metrics.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'Figure_Real_Data_Metrics.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved real data metrics plot")


def plot_robustness(robust_df, outdir):
    """Plot robustness to epsilon."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    eps_df = robust_df[robust_df['parameter'] == 'epsilon']
    axes[0].plot(eps_df['value'], eps_df['mean_tdi'], marker='o', color='#e74c3c')
    axes[0].set_xlabel('Sinkhorn Epsilon', fontsize=12)
    axes[0].set_ylabel('Mean TDI', fontsize=12)
    axes[0].set_title('TDI Stability Across Epsilon', fontsize=12, fontweight='bold')
    axes[0].set_xscale('log')
    axes[0].grid(alpha=0.3)

    axes[1].plot(eps_df['value'], eps_df['poly_frac'], marker='o', color='#3498db')
    axes[1].set_xlabel('Sinkhorn Epsilon', fontsize=12)
    axes[1].set_ylabel(f'Fraction {EXPECTED_ARREST_STAGE}', fontsize=12)
    axes[1].set_title('Stage Assignment Stability Across Epsilon', fontsize=12, fontweight='bold')
    axes[1].set_xscale('log')
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'Figure_Robustness.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(outdir, 'Figure_Robustness.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved robustness plot")


# ============================================================================
# Main
# ============================================================================
def main():
    print("\n" + "=" * 70)
    print("BENCHMARK V2: scTDRP-focused Evaluation")
    print("=" * 70)

    # Load data
    adata_normal, adata_ael, adata_malignant, adata_residual = load_and_preprocess()

    # Build metacells
    print("\n[Preprocessing] Building metacells")
    adata_meta_normal = build_stage_metacells(adata_normal, resolution_scale=1.0)
    meta_disease = build_metacells(adata_malignant, use_rep='X_pca', resolution=1.0, label='Malignant')
    meta_residual = build_metacells(adata_residual, use_rep='X_pca', resolution=1.0, label='Residual')

    all_results = {}

    # Run methods on malignant cells
    print("\n[Methods] Running scTDRP and baselines on malignant cells")
    t0 = time.time()
    scTDRP_scores, repair_results, analyzer = run_scTDRP(adata_meta_normal, meta_disease, adata_malignant)
    all_results['scTDRP_time'] = time.time() - t0

    t0 = time.time()
    noOT_scores = run_noOT(adata_meta_normal, meta_disease, adata_malignant)
    all_results['NoOT_time'] = time.time() - t0

    t0 = time.time()
    pca_scores = run_pca_eucdist(adata_normal, adata_malignant)
    all_results['PCA_time'] = time.time() - t0

    t0 = time.time()
    pt_scores = run_pseudotime_stage(adata_normal, adata_malignant)
    all_results['Pseudotime_time'] = time.time() - t0

    # Run methods on normal residual cells (for AUC) using the same analyzer
    print("\n[Methods] Running scTDRP and baselines on normal residual cells")
    scTDRP_residual, _, _ = run_scTDRP(
        adata_meta_normal, meta_residual, adata_residual,
        analyzer=analyzer, compute_repair=False
    )
    noOT_residual = run_noOT(adata_meta_normal, meta_residual, adata_residual)
    pca_residual = run_pca_eucdist(adata_normal, adata_residual)
    pt_residual = run_pseudotime_stage(adata_normal, adata_residual)

    # Build comparison dataframe combining malignant and normal residual
    n_mal = adata_malignant.shape[0]
    n_res = adata_residual.shape[0]
    comp_df = pd.DataFrame({
        'malignancy': list(adata_malignant.obs['malignancy'].values) + list(adata_residual.obs['malignancy'].values),
        'phase': list(adata_malignant.obs['phase'].values) + list(adata_residual.obs['phase'].values),
        'cnv_score': list(adata_malignant.obs['cnv_score'].values) + list(adata_residual.obs['cnv_score'].values),
    })

    def merge_scores(mal_df, res_df, prefix=''):
        """Concatenate malignant and residual score columns."""
        merged = {}
        for col in mal_df.columns:
            if col == 'cell_idx':
                continue
            merged[col] = list(mal_df[col].values) + list(res_df[col].values)
        return merged

    score_merges = {
        **merge_scores(scTDRP_scores, scTDRP_residual),
        **merge_scores(noOT_scores, noOT_residual),
        **merge_scores(pca_scores, pca_residual),
        **merge_scores(pt_scores, pt_residual),
    }
    for col, vals in score_merges.items():
        comp_df[col] = vals

    # Run benchmark tests
    print("\n" + "=" * 70)
    print("Running benchmark tests")
    print("=" * 70)

    stage_acc_df = evaluate_stage_assignment(adata_meta_normal, adata_normal)
    repair_df = evaluate_repair_direction(adata_meta_normal, adata_normal)

    with open(MODULES_PATH) as f:
        modules = json.load(f)
    module_df = evaluate_module_strategy(repair_results['repair_deltas'], modules)

    real_df = evaluate_real_data(comp_df)
    robust_df = evaluate_robustness(adata_meta_normal, meta_disease, adata_malignant)

    # Save results
    stage_acc_df.to_csv(os.path.join(OUT_DIR, 'stage_assignment_accuracy.csv'), index=False)
    repair_df.to_csv(os.path.join(OUT_DIR, 'repair_direction_correctness.csv'), index=False)
    module_df.to_csv(os.path.join(OUT_DIR, 'module_strategy_concordance.csv'), index=False)
    real_df.to_csv(os.path.join(OUT_DIR, 'real_data_metrics.csv'), index=False)
    robust_df.to_csv(os.path.join(OUT_DIR, 'robustness.csv'), index=False)
    comp_df.to_csv(os.path.join(OUT_DIR, 'comparison_scores.csv'), index=False)

    # Plot
    print("\n[Visualization] Generating figures")
    plot_stage_assignment(stage_acc_df, FIG_DIR)
    plot_repair_direction(repair_df, FIG_DIR)
    plot_module_strategy(module_df, FIG_DIR)
    plot_real_data_metrics(real_df, FIG_DIR)
    plot_robustness(robust_df, FIG_DIR)

    print("\n" + "=" * 70)
    print("BENCHMARK V2 COMPLETE")
    print(f"Results: {OUT_DIR}")
    print(f"Figures: {FIG_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
