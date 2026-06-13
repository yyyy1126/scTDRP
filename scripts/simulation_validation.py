#!/usr/bin/env python3
"""
scTDRP Concept Validation Simulation
====================================
Generate synthetic 5-stage differentiation + arrested disease cells
Validate: (1) OT cost map monotonicity, (2) TDI accuracy, (3) repair directionality
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/data1/yja/zhongzhuan/5.external/scTDRP/src')

from scTDRP.distance import build_stage_distributions, compute_wasserstein_distance, compute_ot_cost_map
from scTDRP.repair import compute_gene_repair_delta
from scTDRP.utils import compute_expression_distribution

# ------------------------------------------------------------------
# 1. Generate synthetic data
# ------------------------------------------------------------------
np.random.seed(42)

N_STAGES = 5          # normal developmental stages
N_CELLS_PER_STAGE = 100
N_GENES = 200         # gene expression space
NOISE_LEVEL = 0.3     # intra-stage biological noise
OFFSET_MAG = 2.0      # disease offset magnitude

# Differentiation axis: cells move along a primary direction in gene space
diff_axis = np.zeros(N_GENES)
diff_axis[:20] = 1.0  # first 20 genes drive differentiation
diff_axis = diff_axis / np.linalg.norm(diff_axis)

# Disease perturbation: specific genes (21-40) are dysregulated
perturb_axis = np.zeros(N_GENES)
perturb_axis[20:40] = 1.0
perturb_axis = perturb_axis / np.linalg.norm(perturb_axis)

normal_data = []
stage_labels = []

for s in range(N_STAGES):
    # Stage centers move along differentiation axis
    center = s * 2.0 * diff_axis
    # Intra-stage noise (Gaussian)
    cells = center + np.random.normal(0, NOISE_LEVEL, size=(N_CELLS_PER_STAGE, N_GENES))
    normal_data.append(cells)
    stage_labels.extend([f"Stage_{s}"] * N_CELLS_PER_STAGE)

normal_data = np.vstack(normal_data)  # (500, 200)

# Disease cells: arrest at Stage_3 (index 3) + perturbation
# They retain the Stage_3 differentiation position but add disease-specific offset
disease_center = 3 * 2.0 * diff_axis + OFFSET_MAG * perturb_axis
disease_data = disease_center + np.random.normal(0, NOISE_LEVEL, size=(N_CELLS_PER_STAGE, N_GENES))
disease_labels = ["Disease"] * N_CELLS_PER_STAGE

# ------------------------------------------------------------------
# 2. Build AnnData objects
# ------------------------------------------------------------------
all_data = np.vstack([normal_data, disease_data])
all_labels = stage_labels + disease_labels

gene_names = [f"gene_{i:03d}" for i in range(N_GENES)]
obs = pd.DataFrame({"cell_type": all_labels})
var = pd.DataFrame(index=gene_names)

adata = sc.AnnData(X=all_data, obs=obs, var=var)
adata.layers["counts"] = adata.X.copy()

# Split into normal and disease
adata_normal = adata[adata.obs["cell_type"].str.startswith("Stage_")].copy()
adata_disease = adata[adata.obs["cell_type"] == "Disease"].copy()

# PCA on normal data
sc.pp.pca(adata_normal, n_comps=20)
# Project disease onto same PCA
adata_disease.obsm["X_pca"] = adata_disease.X @ adata_normal.varm["PCs"]

# ------------------------------------------------------------------
# 3. Build stage distributions (using gene space for repair, PCA for TDI)
# ------------------------------------------------------------------
stage_order = [f"Stage_{s}" for s in range(N_STAGES)]

# For TDI and cost map: use PCA space
stage_dist_pca = build_stage_distributions(
    adata_normal, stage_key="cell_type", stage_order=stage_order, use_rep="X_pca"
)

# For repair: use original gene space
stage_dist_gene = build_stage_distributions(
    adata_normal, stage_key="cell_type", stage_order=stage_order, use_rep="X"
)

# Disease distributions
disease_pca = compute_expression_distribution(adata_disease, use_rep="X_pca")
disease_gene = compute_expression_distribution(adata_disease, use_rep="X")

# ------------------------------------------------------------------
# 4. Validate OT Cost Map Monotonicity
# ------------------------------------------------------------------
print("=" * 60)
print("TEST 1: OT Cost Map Monotonicity")
print("=" * 60)

cost_map = compute_ot_cost_map(stage_dist_pca, reg=0.05)

print("\nAdjacent stage W2 distances:")
distances = []
for i in range(N_STAGES - 1):
    pair = (stage_order[i], stage_order[i+1])
    d = cost_map[pair]
    distances.append(d)
    print(f"  {pair[0]} -> {pair[1]}: W2 = {d:.4f}")

# Check monotonicity: later stages should have accumulated larger distance from Stage_0
print("\nCumulative distances from Stage_0:")
cum_dists = []
for i in range(1, N_STAGES):
    X0, w0 = stage_dist_pca["Stage_0"]
    Xi, wi = stage_dist_pca[f"Stage_{i}"]
    d, _ = compute_wasserstein_distance(X0, Xi, w0, wi, reg=0.05)
    cum_dists.append(d)
    print(f"  Stage_0 -> Stage_{i}: W2 = {d:.4f}")

monotonic = all(cum_dists[i] <= cum_dists[i+1] for i in range(len(cum_dists)-1))
print(f"\n✓ Monotonicity check: {'PASS' if monotonic else 'FAIL'}")

# ------------------------------------------------------------------
# 5. Validate TDI Accuracy
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 2: TDI Accuracy (Disease cells should map to Stage_3)")
print("=" * 60)

# Compute TDI for each disease cell
tdi_results = []
for idx in range(min(20, disease_pca.shape[0])):  # sample 20 cells
    X_cell = disease_pca[idx:idx+1, :]
    w_cell = np.ones(1)
    
    dists = {}
    for stage, (X_stage, w_stage) in stage_dist_pca.items():
        d, _ = compute_wasserstein_distance(X_cell, X_stage, w_cell, w_stage, reg=0.05)
        dists[stage] = d
    
    min_stage = min(dists, key=dists.get)
    tdi_results.append(min_stage)

from collections import Counter
stage_counts = Counter(tdi_results)
print(f"\nTDI attribution for 20 disease cells:")
for stage in stage_order:
    print(f"  {stage}: {stage_counts[stage]} cells")

pred_stage = stage_counts.most_common(1)[0][0]
correct = (pred_stage == "Stage_3")
print(f"\n✓ Ground truth arrest: Stage_3")
print(f"✓ Predicted arrest: {pred_stage}")
print(f"✓ TDI accuracy: {'PASS' if correct else 'FAIL'}")

# ------------------------------------------------------------------
# 6. Validate Repair Pathway Directionality
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 3: Repair Pathway Directionality")
print("=" * 60)

# Disease -> Terminal (Stage_4)
X_disease = disease_gene
X_terminal, w_terminal = stage_dist_gene["Stage_4"]

# Uniform weights for disease cells
w_disease = np.ones(X_disease.shape[0]) / X_disease.shape[0]

# Compute transport plan
dist, transport_plan = compute_wasserstein_distance(
    X_disease, X_terminal, w_disease, w_terminal, reg=0.05
)

# Compute gene repair deltas
repair_deltas = compute_gene_repair_delta(
    X_disease, X_terminal, transport_plan, gene_names=gene_names
)

# Genes 0-19: differentiation driver (should be upregulated to reach terminal)
# Genes 20-39: disease perturbation (should be downregulated to remove offset)
# Genes 40+: neutral (should have delta ~ 0)

diff_genes = gene_names[:20]
perturb_genes = gene_names[20:40]
neutral_genes = gene_names[40:]

diff_mean = np.mean([repair_deltas[g] for g in diff_genes])
perturb_mean = np.mean([repair_deltas[g] for g in perturb_genes])
neutral_mean = np.mean([repair_deltas[g] for g in neutral_genes])

print(f"\nMean repair delta by gene category:")
print(f"  Differentiation drivers (genes 0-19):  {diff_mean:+.4f}  (expected: positive)")
print(f"  Disease perturbation (genes 20-39):    {perturb_mean:+.4f}  (expected: negative)")
print(f"  Neutral genes (genes 40+):             {neutral_mean:+.4f}  (expected: ~0)")

direction_ok = (diff_mean > 0) and (perturb_mean < 0) and (abs(neutral_mean) < 0.5)
print(f"\n✓ Repair directionality: {'PASS' if direction_ok else 'FAIL'}")

# Top up/down targets
top_up = sorted(repair_deltas.items(), key=lambda x: x[1], reverse=True)[:10]
top_down = sorted(repair_deltas.items(), key=lambda x: x[1])[:10]

print(f"\nTop 10 up-regulation targets:")
for g, d in top_up:
    cat = "diff" if g in diff_genes else ("perturb" if g in perturb_genes else "neutral")
    print(f"  {g}: {d:+.4f} [{cat}]")

print(f"\nTop 10 down-regulation targets:")
for g, d in top_down:
    cat = "diff" if g in diff_genes else ("perturb" if g in perturb_genes else "neutral")
    print(f"  {g}: {d:+.4f} [{cat}]")

# ------------------------------------------------------------------
# 7. Parameter Sensitivity (epsilon)
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 4: Parameter Sensitivity (regularization epsilon)")
print("=" * 60)

epsilons = [0.001, 0.01, 0.05, 0.1, 0.5]
tdi_eps_results = []

for eps in epsilons:
    dists = {}
    for stage, (X_stage, w_stage) in stage_dist_pca.items():
        d, _ = compute_wasserstein_distance(
            disease_pca[:1], X_stage, np.ones(1), w_stage, reg=eps
        )
        dists[stage] = d
    min_stage = min(dists, key=dists.get)
    tdi_eps_results.append(min_stage)
    print(f"  eps={eps:5.3f}: TDI = {min_stage}, distances = {[f'{d:.4f}' for d in dists.values()]}")

stable = all(s == "Stage_3" for s in tdi_eps_results)
print(f"\n✓ TDI stability across epsilon: {'PASS' if stable else 'FAIL'}")

# ------------------------------------------------------------------
# 8. Visualization
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("Generating figures...")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: PCA visualization
ax = axes[0, 0]
for s in range(N_STAGES):
    mask = adata_normal.obs["cell_type"] == f"Stage_{s}"
    ax.scatter(
        adata_normal.obsm["X_pca"][mask, 0],
        adata_normal.obsm["X_pca"][mask, 1],
        label=f"Stage {s}", alpha=0.5, s=10
    )
ax.scatter(
    adata_disease.obsm["X_pca"][:, 0],
    adata_disease.obsm["X_pca"][:, 1],
    label="Disease", c="red", marker="x", s=30
)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("A. Simulated Developmental Landscape")
ax.legend(loc="upper left", fontsize=7)

# Panel B: OT Cost Map
ax = axes[0, 1]
pairs = [f"S{i}\u2192S{i+1}" for i in range(N_STAGES - 1)]
vals = distances
bars = ax.bar(pairs, vals, color="steelblue", edgecolor="black")
ax.set_ylabel("Wasserstein-2 Distance")
ax.set_title("B. OT Cost Map (Adjacent Stages)")
ax.axhline(y=0, color="black", linewidth=0.5)
for bar, val in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
            f"{val:.3f}", ha="center", va="bottom", fontsize=8)

# Panel C: TDI distribution
ax = axes[1, 0]
all_tdi_dists = []
for idx in range(disease_pca.shape[0]):
    X_cell = disease_pca[idx:idx+1, :]
    w_cell = np.ones(1)
    dists = []
    for stage in stage_order:
        X_stage, w_stage = stage_dist_pca[stage]
        d, _ = compute_wasserstein_distance(X_cell, X_stage, w_cell, w_stage, reg=0.05)
        dists.append(d)
    all_tdi_dists.append(dists)

all_tdi_dists = np.array(all_tdi_dists)
mean_dists = all_tdi_dists.mean(axis=0)
std_dists = all_tdi_dists.std(axis=0)

x_pos = np.arange(N_STAGES)
ax.bar(x_pos, mean_dists, yerr=std_dists, capsize=5, color="coral", edgecolor="black")
ax.set_xticks(x_pos)
ax.set_xticklabels(stage_order, rotation=45, ha="right")
ax.set_ylabel("Mean W2 Distance")
ax.set_title("C. TDI: Disease-to-Normal-Stage Distances")
ax.axvline(x=3, color="green", linestyle="--", linewidth=2, label="Ground truth (Stage_3)")
ax.legend()

# Panel D: Repair deltas by gene category
ax = axes[1, 1]
categories = ["Diff. drivers\n(genes 0-19)", "Disease perturb.\n(genes 20-39)", "Neutral\n(genes 40+)"]
means = [diff_mean, perturb_mean, neutral_mean]
colors = ["#2ecc71", "#e74c3c", "#95a5a6"]
bars = ax.bar(categories, means, color=colors, edgecolor="black")
ax.axhline(y=0, color="black", linewidth=0.5)
ax.set_ylabel("Mean Repair Delta")
ax.set_title("D. Repair Directionality by Gene Category")
for bar, val in zip(bars, means):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05 * np.sign(val),
            f"{val:+.3f}", ha="center", va="bottom" if val > 0 else "top", fontsize=9)

plt.tight_layout()
out_path = "/data1/yja/zhongzhuan/5.external/scTDRP/scripts/simulation_validation_results.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"\nFigure saved to: {out_path}")

# ------------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("VALIDATION SUMMARY")
print("=" * 60)
print(f"  [1] OT Cost Map Monotonicity:  {'PASS' if monotonic else 'FAIL'}")
print(f"  [2] TDI Accuracy:               {'PASS' if correct else 'FAIL'}")
print(f"  [3] Repair Directionality:      {'PASS' if direction_ok else 'FAIL'}")
print(f"  [4] Parameter Stability:        {'PASS' if stable else 'FAIL'}")
print("=" * 60)
