#!/usr/bin/env python3
"""
Simulation Benchmark: Quantitative evaluation of scTDRP vs baselines on controlled data.

This script creates a synthetic dataset with known ground truth:
  - 5-stage normal differentiation (Stage_0 -> Stage_4)
  - Disease cells arrested at Stage_2 with a known perturbation offset
  - Ground truth: each disease cell needs a specific transcriptional shift to reach Stage_4

Methods compared:
  1. scTDRP (Full)
  2. scTDRP-NoOT (Ablation: Euclidean instead of Wasserstein)
  3. Waddington-OT (if available)
  4. MOSCOT (if available)
  5. Random baseline

Evaluation metrics:
  1. Arrest Stage Accuracy: fraction of disease cells correctly assigned to Stage_2
  2. Repair Direction Cosine Similarity: cosine similarity between inferred and true repair vector
  3. Top-K Target Recovery: overlap of top up/down targets with ground truth
  4. Computation time
"""

import os
import sys
import time
import json
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

OUT_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/results_simulation_benchmark"
FIG_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/figures_simulation_benchmark"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

np.random.seed(42)

# Optional imports
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
except Exception:
    sys.path.insert(0, "/data1/yja/zhongzhuan/5.external/scTDRP/src")
    from scTDRP import TDRPAnalyzer


# =============================================================================
# Simple metacell builder
# =============================================================================
def build_metacells_simple(adata, use_rep="X", n_metacells=10, label=None):
    """Simple k-means based metacell builder for small datasets."""
    from sklearn.cluster import KMeans
    X = adata.obsm[use_rep] if use_rep in adata.obsm else adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    n_meta = min(n_metacells, X.shape[0])
    kmeans = KMeans(n_clusters=n_meta, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X)
    X_meta = kmeans.cluster_centers_
    obs_meta = pd.DataFrame({
        'metacell_id': [f"{label}_MC{i}" for i in range(n_meta)],
        'n_cells': [int((labels == i).sum()) for i in range(n_meta)],
    })
    meta_ad = sc.AnnData(X=X_meta, obs=obs_meta, var=adata.var.copy())
    meta_ad.obs_names = obs_meta['metacell_id'].values.astype(str)
    meta_ad.uns['cell_labels'] = labels
    return meta_ad


# =============================================================================
# Simulation
# =============================================================================
def simulate_data(
    n_stages=5,
    n_cells_per_stage=200,
    n_genes=100,
    n_disease_cells=100,
    arrest_stage=2,
    terminal_stage=4,
    perturbation_vector=None,
    noise_std=0.5,
):
    """
    Simulate a linear differentiation trajectory with disease arrest.

    Ground truth:
      - Normal cells progress along a straight line in gene space
      - Disease cells are arrested at arrest_stage with an additional perturbation
      - True repair = (terminal_mean - disease_mean) for each disease cell
    """
    print("=" * 60)
    print("[Simulation] Generating synthetic data...")
    print("=" * 60)

    # Define stage means along a line in gene space
    # Use first 2 genes as "main axis", rest as noise
    stage_means = []
    for i in range(n_stages):
        mean = np.zeros(n_genes)
        mean[0] = i * 5.0  # X-axis: progression
        mean[1] = i * 3.0  # Y-axis: correlated progression
        stage_means.append(mean)

    # Generate normal cells
    normal_X = []
    normal_labels = []
    for i, mean in enumerate(stage_means):
        cells = mean + np.random.normal(0, noise_std, size=(n_cells_per_stage, n_genes))
        normal_X.append(cells)
        normal_labels.extend([f"Stage_{i}"] * n_cells_per_stage)

    normal_X = np.vstack(normal_X)

    # Generate disease cells: arrested at arrest_stage + perturbation
    disease_mean = stage_means[arrest_stage].copy()
    if perturbation_vector is None:
        # Default: push away from terminal in a specific direction
        perturbation_vector = np.zeros(n_genes)
        perturbation_vector[0] = -2.0  # Backward on main axis
        perturbation_vector[1] = 3.0   # Sideways offset
    disease_mean = disease_mean + perturbation_vector

    disease_X = disease_mean + np.random.normal(0, noise_std, size=(n_disease_cells, n_genes))
    disease_labels = ['Disease'] * n_disease_cells

    # Ground truth repair: for each disease cell, vector to terminal stage mean
    terminal_mean = stage_means[terminal_stage]
    true_repair = terminal_mean - disease_X  # (n_disease, n_genes)

    # Gene names
    gene_names = [f"gene_{i:03d}" for i in range(n_genes)]

    # Create AnnData objects
    adata_normal = sc.AnnData(X=normal_X)
    adata_normal.obs_names = [f"normal_{i}" for i in range(normal_X.shape[0])]
    adata_normal.var_names = gene_names
    adata_normal.obs['stage'] = normal_labels

    adata_disease = sc.AnnData(X=disease_X)
    adata_disease.obs_names = [f"disease_{i}" for i in range(disease_X.shape[0])]
    adata_disease.var_names = gene_names
    adata_disease.obs['stage'] = disease_labels

    print(f"  Normal: {adata_normal.shape}")
    print(f"  Disease: {adata_disease.shape}")
    print(f"  Arrest stage: Stage_{arrest_stage}")
    print(f"  Terminal stage: Stage_{terminal_stage}")
    print(f"  Perturbation: X={perturbation_vector[0]:.2f}, Y={perturbation_vector[1]:.2f}")
    print(f"  True repair mean: X={true_repair[:,0].mean():.2f}, Y={true_repair[:,1].mean():.2f}")

    return adata_normal, adata_disease, true_repair, stage_means, perturbation_vector


# =============================================================================
# Method 1: scTDRP Full
# =============================================================================
def run_scTDRP_sim(adata_normal, adata_disease, stage_order, terminal_stage):
    print("\n" + "=" * 60)
    print("[Method 1] scTDRP (Full with Metacells)")
    print("=" * 60)
    t0 = time.time()

    # Build normal metacells per stage
    meta_normal_list = []
    for stage in stage_order:
        mask = adata_normal.obs['stage'] == stage
        stage_ad = adata_normal[mask].copy()
        n_meta = min(10, stage_ad.shape[0])
        meta_stage = build_metacells_simple(stage_ad, use_rep='X', n_metacells=n_meta, label=stage)
        meta_stage.obs['stage'] = stage
        meta_normal_list.append(meta_stage)

    adata_meta_normal = sc.concat(meta_normal_list, label="_batch")
    adata_meta_normal.obs['stage'] = [s for ad in meta_normal_list for s in [ad.obs['stage'].iloc[0]] * ad.shape[0]]

    # Build disease metacells
    n_disease_meta = min(10, adata_disease.shape[0])
    meta_disease = build_metacells_simple(adata_disease, use_rep='X', n_metacells=n_disease_meta, label='Disease')

    print(f"  Normal metacells: {adata_meta_normal.shape[0]}, Disease metacells: {meta_disease.shape[0]}")

    analyzer = TDRPAnalyzer(
        normal_adata=adata_meta_normal,
        stage_key='stage',
        terminal_stage=terminal_stage,
        stage_order=stage_order,
        use_rep='X',
        n_top_genes=None,
    )
    analyzer.prepare_data()
    analyzer.build_ot_cost_map(metric='sqeuclidean', reg=0.01)
    tdi_df = analyzer.compute_tdi(meta_disease, metric='sqeuclidean', reg=0.01)
    repair = analyzer.infer_repair_pathway(meta_disease, metric='sqeuclidean', reg=0.01, top_n=20)

    # Map metacell results back to original cells
    disease_labels = meta_disease.uns['cell_labels']
    meta_tdi = dict(zip(meta_disease.obs_names, tdi_df['TDI'].values))
    meta_stage = dict(zip(meta_disease.obs_names, tdi_df['Best_Match_Stage'].values))

    cell_tdi = []
    cell_stage = []
    for i, cell in enumerate(adata_disease.obs_names):
        cl = disease_labels[i]
        meta_name = f"Disease_MC{cl}"
        cell_tdi.append(meta_tdi.get(meta_name, np.nan))
        cell_stage.append(meta_stage.get(meta_name, 'Unknown'))

    # Repair deltas (population level)
    repair_deltas = repair['repair_deltas']
    gene_names = adata_disease.var_names.tolist()
    repair_matrix = np.array([repair_deltas.get(g, 0.0) for g in gene_names])
    repair_matrix_per_cell = np.tile(repair_matrix, (adata_disease.shape[0], 1))

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.2f}s")

    return {
        'best_stages': np.array(cell_stage),
        'repair_matrix': repair_matrix_per_cell,
        'tdi': np.array(cell_tdi),
        'time': elapsed,
    }


# =============================================================================
# Method 2: scTDRP-NoOT (Ablation)
# =============================================================================
def run_noOT_sim(adata_normal, adata_disease, stage_order):
    print("\n" + "=" * 60)
    print("[Method 2] scTDRP-NoOT (Euclidean centroids on metacells)")
    print("=" * 60)
    t0 = time.time()

    # Use same metacells as scTDRP for fair comparison
    meta_normal_list = []
    for stage in stage_order:
        mask = adata_normal.obs['stage'] == stage
        stage_ad = adata_normal[mask].copy()
        n_meta = min(10, stage_ad.shape[0])
        meta_stage = build_metacells_simple(stage_ad, use_rep='X', n_metacells=n_meta, label=stage)
        meta_stage.obs['stage'] = stage
        meta_normal_list.append(meta_stage)

    adata_meta_normal = sc.concat(meta_normal_list, label="_batch")
    adata_meta_normal.obs['stage'] = [s for ad in meta_normal_list for s in [ad.obs['stage'].iloc[0]] * ad.shape[0]]

    n_disease_meta = min(10, adata_disease.shape[0])
    meta_disease = build_metacells_simple(adata_disease, use_rep='X', n_metacells=n_disease_meta, label='Disease')

    # Compute stage centroids from normal metacells
    centroids = {}
    for stage in stage_order:
        mask = adata_meta_normal.obs['stage'] == stage
        centroids[stage] = adata_meta_normal.X[mask].mean(axis=0)

    centroid_matrix = np.array([centroids[s] for s in stage_order])
    distances = cdist(meta_disease.X, centroid_matrix, metric='euclidean')
    best_idx = distances.argmin(axis=1)
    best_stages = [stage_order[i] for i in best_idx]

    # Map back to original cells
    disease_labels = meta_disease.uns['cell_labels']
    meta_stage_map = dict(zip(meta_disease.obs_names, best_stages))
    cell_stage = []
    for i in range(adata_disease.shape[0]):
        cl = disease_labels[i]
        meta_name = f"Disease_MC{cl}"
        cell_stage.append(meta_stage_map.get(meta_name, 'Unknown'))

    # Repair = terminal centroid - disease cell (population level)
    terminal_idx = len(stage_order) - 1
    repair_matrix = centroid_matrix[terminal_idx] - adata_disease.X

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.2f}s")

    return {
        'best_stages': np.array(cell_stage),
        'repair_matrix': repair_matrix,
        'tdi': distances.min(axis=1),  # metacell-level, approximate
        'time': elapsed,
    }


# =============================================================================
# Method 3: Waddington-OT
# =============================================================================
def run_wot_sim(adata_normal, adata_disease, stage_order):
    print("\n" + "=" * 60)
    print("[Method 3] Waddington-OT")
    print("=" * 60)

    if not WOT_AVAILABLE:
        print("  [SKIP] WOT not available")
        return None

    t0 = time.time()

    # Assign days
    adata_normal_copy = adata_normal.copy()
    adata_disease_copy = adata_disease.copy()
    day_map = {stage: float(i) for i, stage in enumerate(stage_order)}
    adata_normal_copy.obs['day'] = adata_normal_copy.obs['stage'].map(day_map).astype(float)
    adata_disease_copy.obs['day'] = float(len(stage_order))

    adata_combined = sc.concat([adata_normal_copy, adata_disease_copy], label='source',
                                keys=['normal', 'disease'], index_unique='-')

    try:
        ot_model = wot.ot.OTModel(adata_combined, day_field='day',
                                   epsilon=0.05, lambda1=1, lambda2=50)

        t_last_normal = float(len(stage_order) - 1)
        t_disease = float(len(stage_order))
        tmap = ot_model.compute_transport_map(t_last_normal, t_disease)

        # WOT returns AnnData in newer versions
        if hasattr(tmap, 'X'):
            tmap_matrix = tmap.X
            if hasattr(tmap_matrix, 'toarray'):
                tmap_matrix = tmap_matrix.toarray()
            tmap_matrix = np.asarray(tmap_matrix)
            tmap_rows = tmap.obs_names.tolist()
            tmap_cols = tmap.var_names.tolist()
        else:
            tmap_matrix = np.asarray(tmap)
            tmap_rows = adata_normal_copy[adata_normal_copy.obs['day'] == t_last_normal].obs_names.tolist()
            tmap_cols = adata_disease_copy.obs_names.tolist()

        disease_names = adata_disease_copy.obs_names.tolist()
        normal_terminal_names = adata_normal_copy[adata_normal_copy.obs['day'] == t_last_normal].obs_names.tolist()

        row_idx = [tmap_rows.index(n) for n in normal_terminal_names if n in tmap_rows]
        col_idx = [tmap_cols.index(n) for n in disease_names if n in tmap_cols]

        if len(row_idx) == 0 or len(col_idx) == 0:
            print("  [WARNING] WOT tmap name mismatch")
            col_sums = tmap_matrix.sum(axis=0)
        else:
            tmap_sub = tmap_matrix[np.ix_(row_idx, col_idx)]
            col_sums = tmap_sub.sum(axis=0)

        col_sums = np.asarray(col_sums).flatten()
        wot_dev = 1.0 - (col_sums - col_sums.min()) / (col_sums.max() - col_sums.min() + 1e-10)

        best_stages = np.array([stage_order[-1]] * adata_disease.shape[0])

        terminal_mean = adata_normal_copy[adata_normal_copy.obs['day'] == t_last_normal].X.mean(axis=0)
        repair_matrix = terminal_mean - adata_disease_copy.X

        elapsed = time.time() - t0
        print(f"  Elapsed: {elapsed:.2f}s")

        return {
            'best_stages': best_stages,
            'repair_matrix': repair_matrix,
            'tdi': wot_dev,
            'time': elapsed,
        }
    except Exception as e:
        print(f"  [ERROR] WOT failed: {e}")
        import traceback
        traceback.print_exc()
        return None


# =============================================================================
# Method 4: MOSCOT
# =============================================================================
def run_moscot_sim(adata_normal, adata_disease, stage_order):
    print("\n" + "=" * 60)
    print("[Method 4] MOSCOT")
    print("=" * 60)

    if not MOSCOT_AVAILABLE:
        print("  [SKIP] MOSCOT not available")
        return None

    t0 = time.time()

    adata_normal_copy = adata_normal.copy()
    adata_disease_copy = adata_disease.copy()
    day_map = {stage: i for i, stage in enumerate(stage_order)}
    adata_normal_copy.obs['day'] = adata_normal_copy.obs['stage'].map(day_map).astype(int)
    adata_disease_copy.obs['day'] = len(stage_order)

    adata_combined = sc.concat([adata_normal_copy, adata_disease_copy], label='source',
                                keys=['normal', 'disease'], index_unique='-')

    try:
        tp = TemporalProblem(adata_combined)
        tp = tp.prepare(time_key="day")
        tp = tp.solve(epsilon=0.01)

        t_last = len(stage_order) - 1
        t_disease = len(stage_order)
        coupling = tp[(t_last, t_disease)].solution.transport_matrix
        col_sums = np.array(coupling.sum(axis=0)).flatten()

        moscot_dev = 1.0 - (col_sums - col_sums.min()) / (col_sums.max() - col_sums.min() + 1e-10)

        terminal_mean = adata_normal_copy[adata_normal_copy.obs['day'] == t_last].X.mean(axis=0)
        repair_matrix = terminal_mean - adata_disease_copy.X

        best_stages = np.array([stage_order[-1]] * adata_disease.shape[0])

        elapsed = time.time() - t0
        print(f"  Elapsed: {elapsed:.2f}s")

        return {
            'best_stages': best_stages,
            'repair_matrix': repair_matrix,
            'tdi': moscot_dev,
            'time': elapsed,
        }
    except Exception as e:
        print(f"  [ERROR] MOSCOT failed: {e}")
        import traceback
        traceback.print_exc()
        return None


# =============================================================================
# Method 5: Random Baseline
# =============================================================================
def run_random_sim(adata_disease, stage_order):
    print("\n" + "=" * 60)
    print("[Method 5] Random Baseline")
    print("=" * 60)
    t0 = time.time()

    n = adata_disease.shape[0]
    best_stages = np.random.choice(stage_order, size=n)
    repair_matrix = np.random.normal(0, 1, size=adata_disease.X.shape)
    tdi = np.random.uniform(0, 1, size=n)

    elapsed = time.time() - t0
    print(f"  Elapsed: {elapsed:.4f}s")

    return {
        'best_stages': best_stages,
        'repair_matrix': repair_matrix,
        'tdi': tdi,
        'time': elapsed,
    }


# =============================================================================
# Evaluation
# =============================================================================
def evaluate_simulation(results, true_repair, arrest_stage, stage_order):
    print("\n" + "=" * 60)
    print("[Evaluation] Simulation Benchmark Results")
    print("=" * 60)

    records = []
    n_disease = true_repair.shape[0]
    true_stage = f"Stage_{arrest_stage}"

    for method_name, res in results.items():
        if res is None:
            continue

        # 1. Arrest stage accuracy
        best_stages = res['best_stages']
        if best_stages is not None:
            acc = (best_stages == true_stage).mean()
        else:
            acc = np.nan

        # 2. Repair direction cosine similarity
        repair_matrix = res['repair_matrix']
        if repair_matrix is not None and repair_matrix.shape == true_repair.shape:
            # Cosine similarity per cell, then average
            cos_sims = []
            for i in range(n_disease):
                a = true_repair[i:i+1, :]
                b = repair_matrix[i:i+1, :]
                if np.linalg.norm(b) > 1e-10:
                    sim = cosine_similarity(a, b)[0, 0]
                    cos_sims.append(sim)
            mean_cos = np.mean(cos_sims) if cos_sims else np.nan
        else:
            mean_cos = np.nan

        # 3. Top-10 target recovery (directional)
        # For population-level methods, compare population mean repair vector
        pop_true = true_repair.mean(axis=0)
        pop_pred = repair_matrix.mean(axis=0)

        top10_true_up = np.argsort(pop_true)[-10:]
        top10_pred_up = np.argsort(pop_pred)[-10:]
        up_recovery = len(set(top10_true_up) & set(top10_pred_up)) / 10.0

        top10_true_down = np.argsort(pop_true)[:10]
        top10_pred_down = np.argsort(pop_pred)[:10]
        down_recovery = len(set(top10_true_down) & set(top10_pred_down)) / 10.0

        records.append({
            'Method': method_name,
            'Stage_Accuracy': acc,
            'Repair_CosineSim': mean_cos,
            'Top10_Up_Recovery': up_recovery,
            'Top10_Down_Recovery': down_recovery,
            'Time_s': res['time'],
        })

        print(f"\n  {method_name}:")
        print(f"    Stage Accuracy: {acc:.1%}")
        print(f"    Repair Cosine Sim: {mean_cos:.3f}")
        print(f"    Top-10 Up Recovery: {up_recovery:.1%}")
        print(f"    Top-10 Down Recovery: {down_recovery:.1%}")
        print(f"    Time: {res['time']:.2f}s")

    df = pd.DataFrame(records)
    df.to_csv(os.path.join(OUT_DIR, 'simulation_benchmark_metrics.csv'), index=False)
    return df


# =============================================================================
# Visualization
# =============================================================================
def plot_simulation_results(eval_df, stage_order, results, true_repair, stage_means, perturbation_vector):
    print("\n" + "=" * 60)
    print("[Visualization] Generating simulation figures...")
    print("=" * 60)

    # --- Figure 1: Metrics barplot ---
    metrics = ['Stage_Accuracy', 'Repair_CosineSim', 'Top10_Up_Recovery', 'Top10_Down_Recovery']
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]

    colors = plt.cm.tab10(np.linspace(0, 1, len(eval_df)))
    for ax, metric in zip(axes, metrics):
        bars = ax.bar(eval_df['Method'], eval_df[metric], color=colors)
        ax.set_ylim(0, 1.05)
        ax.set_title(metric.replace('_', ' '), fontsize=12, fontweight='bold')
        ax.set_ylabel("Score")
        ax.set_xticklabels(eval_df['Method'], rotation=30, ha='right')
        for bar, val in zip(bars, eval_df[metric]):
            if not pd.isna(val):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                        f"{val:.2f}", ha='center', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Sim_Metrics.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Sim_Metrics.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIG_DIR}/Figure_Sim_Metrics.pdf")

    # --- Figure 2: 2D trajectory visualization ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot normal trajectory
    ax = axes[0]
    for i, stage in enumerate(stage_order):
        # We don't have the original normal X here, reconstruct roughly
        pass

    # Actually, let's just plot the stage means and disease cells
    ax = axes[0]
    for i, mean in enumerate(stage_means):
        ax.scatter(mean[0], mean[1], s=200, c='blue', marker='o', zorder=5)
        ax.text(mean[0] + 0.1, mean[1] + 0.1, stage_order[i], fontsize=10, fontweight='bold')
    # Disease mean
    disease_mean = stage_means[2] + perturbation_vector
    ax.scatter(disease_mean[0], disease_mean[1], s=300, c='red', marker='X', zorder=5)
    ax.text(disease_mean[0] + 0.1, disease_mean[1] + 0.1, 'Disease', fontsize=11, color='red', fontweight='bold')
    # Terminal
    ax.scatter(stage_means[-1][0], stage_means[-1][1], s=300, c='green', marker='*', zorder=5)
    ax.text(stage_means[-1][0] + 0.1, stage_means[-1][1] + 0.1, 'Terminal', fontsize=11, color='green', fontweight='bold')
    # Arrows: true repair
    ax.arrow(disease_mean[0], disease_mean[1],
             stage_means[-1][0] - disease_mean[0],
             stage_means[-1][1] - disease_mean[1],
             head_width=0.3, head_length=0.3, fc='black', ec='black', linewidth=2, zorder=4)
    ax.set_title("Ground Truth: Disease -> Terminal", fontsize=13, fontweight='bold')
    ax.set_xlabel("Gene 1 (Main axis)")
    ax.set_ylabel("Gene 2 (Correlated axis)")
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Plot repair vectors comparison
    ax = axes[1]
    pop_true = true_repair.mean(axis=0)
    methods_to_plot = [m for m in results if results[m] is not None]
    n_show = 10
    x_genes = np.arange(n_show)
    width = 0.15
    n_methods = len(methods_to_plot)

    for i, method in enumerate(methods_to_plot):
        pop_pred = results[method]['repair_matrix'].mean(axis=0)
        offset = (i - n_methods / 2) * width
        ax.bar(x_genes + offset, pop_pred[:n_show], width, label=method, alpha=0.8)
    # Ground truth
    offset = (-n_methods / 2 - 1) * width
    ax.bar(x_genes + offset, pop_true[:n_show], width, label='Ground Truth', color='black', alpha=0.5)
    ax.set_title("Repair Vectors (First 10 Genes)", fontsize=13, fontweight='bold')
    ax.set_xlabel("Gene Index")
    ax.set_ylabel("Repair Delta")
    ax.set_xticks(x_genes)
    ax.set_xticklabels([f"G{i}" for i in range(n_show)])
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Sim_Trajectory.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Sim_Trajectory.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIG_DIR}/Figure_Sim_Trajectory.pdf")

    # --- Figure 3: Runtime ---
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.barplot(data=eval_df, x='Method', y='Time_s', palette='viridis', ax=ax)
    ax.set_title("Computation Time", fontsize=13, fontweight='bold')
    ax.set_ylabel("Time (seconds)")
    for i, row in eval_df.iterrows():
        ax.text(i, row['Time_s'] + eval_df['Time_s'].max() * 0.02,
                f"{row['Time_s']:.2f}s", ha='center', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Sim_Runtime.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'Figure_Sim_Runtime.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {FIG_DIR}/Figure_Sim_Runtime.pdf")


# =============================================================================
# Main
# =============================================================================
def main():
    print("\n" + "=" * 70)
    print("SIMULATION BENCHMARK")
    print("=" * 70)

    stage_order = [f"Stage_{i}" for i in range(5)]
    terminal_stage = "Stage_4"
    arrest_stage = 2

    adata_normal, adata_disease, true_repair, stage_means, perturbation_vector = simulate_data(
        n_stages=5,
        n_cells_per_stage=500,
        n_genes=100,
        n_disease_cells=200,
        arrest_stage=arrest_stage,
        terminal_stage=4,
        perturbation_vector=None,
        noise_std=0.3,
    )

    results = {}
    results['scTDRP'] = run_scTDRP_sim(adata_normal, adata_disease, stage_order, terminal_stage)
    results['scTDRP-NoOT'] = run_noOT_sim(adata_normal, adata_disease, stage_order)
    results['WOT'] = run_wot_sim(adata_normal, adata_disease, stage_order)
    results['MOSCOT'] = run_moscot_sim(adata_normal, adata_disease, stage_order)
    results['Random'] = run_random_sim(adata_disease, stage_order)

    eval_df = evaluate_simulation(results, true_repair, arrest_stage, stage_order)
    plot_simulation_results(eval_df, stage_order, results, true_repair, stage_means, perturbation_vector)

    print("\n" + "=" * 70)
    print("SIMULATION BENCHMARK COMPLETE")
    print(f"Results: {OUT_DIR}")
    print(f"Figures: {FIG_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
