#!/usr/bin/env python3
"""
scTDRP Concept Validation Simulation v2
========================================
Fixed: use EMD (no normalization) for TDI to ensure comparability
"""

import numpy as np
import pandas as pd
import scanpy as sc
import ot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/data1/yja/zhongzhuan/5.external/scTDRP/src')

from scTDRP.distance import compute_wasserstein_distance, compute_ot_cost_map
from scTDRP.repair import compute_gene_repair_delta
from scTDRP.utils import compute_expression_distribution

# ------------------------------------------------------------------
# 1. Generate synthetic data (higher dim, more cells)
# ------------------------------------------------------------------
np.random.seed(42)

N_STAGES = 5
N_CELLS_PER_STAGE = 200
N_GENES = 500         # increased for more stable OT
NOISE_LEVEL = 0.4
OFFSET_MAG = 3.0      # larger offset for clearer signal

# Differentiation axis: first 30 genes
diff_axis = np.zeros(N_GENES)
diff_axis[:30] = np.random.randn(30)
diff_axis = diff_axis / np.linalg.norm(diff_axis)

# Disease perturbation: next 30 genes, orthogonal to diff_axis
perturb_axis = np.random.randn(N_GENES)
# Gram-Schmidt orthogonalize
perturb_axis = perturb_axis - np.dot(perturb_axis, diff_axis) * diff_axis
perturb_axis = perturb_axis / np.linalg.norm(perturb_axis)

normal_data = []
stage_labels = []

for s in range(N_STAGES):
    center = s * 2.5 * diff_axis
    cells = center + np.random.normal(0, NOISE_LEVEL, size=(N_CELLS_PER_STAGE, N_GENES))
    normal_data.append(cells)
    stage_labels.extend([f"Stage_{s}"] * N_CELLS_PER_STAGE)

normal_data = np.vstack(normal_data)

# Disease cells: arrest at Stage_3 + orthogonal perturbation
disease_center = 3 * 2.5 * diff_axis + OFFSET_MAG * perturb_axis
disease_data = disease_center + np.random.normal(0, NOISE_LEVEL, size=(N_CELLS_PER_STAGE, N_GENES))

# ------------------------------------------------------------------
# 2. Build AnnData
# ------------------------------------------------------------------
all_data = np.vstack([normal_data, disease_data])
all_labels = stage_labels + ["Disease"] * N_CELLS_PER_STAGE

gene_names = [f"gene_{i:03d}" for i in range(N_GENES)]
obs = pd.DataFrame({"cell_type": all_labels})
var = pd.DataFrame(index=gene_names)

adata = sc.AnnData(X=all_data, obs=obs, var=var)
adata_normal = adata[adata.obs["cell_type"].str.startswith("Stage_")].copy()
adata_disease = adata[adata.obs["cell_type"] == "Disease"].copy()

# PCA
sc.pp.pca(adata_normal, n_comps=20)
adata_disease.obsm["X_pca"] = adata_disease.X @ adata_normal.varm["PCs"]

# ------------------------------------------------------------------
# Helper: EMD-based W2 (no normalization, comparable across stages)
# ------------------------------------------------------------------
def emd_w2(X_a, w_a, X_b, w_b):
    M = ot.dist(X_a, X_b, metric="sqeuclidean")
    plan = ot.emd(w_a, w_b, M)
    return float(np.sum(plan * M)), plan

# ------------------------------------------------------------------
# 3. Stage distributions
# ------------------------------------------------------------------
stage_order = [f"Stage_{s}" for s in range(N_STAGES)]

stage_dist_pca = {}
stage_dist_gene = {}
for stage in stage_order:
    mask = adata_normal.obs["cell_type"] == stage
    X_pca = adata_normal.obsm["X_pca"][mask]
    X_gene = adata_normal.X[mask]
    w = np.ones(X_pca.shape[0]) / X_pca.shape[0]
    stage_dist_pca[stage] = (X_pca, w)
    stage_dist_gene[stage] = (X_gene, w)

disease_pca = adata_disease.obsm["X_pca"]
disease_gene = adata_disease.X

# Disease metacell (average)
disease_pca_mean = disease_pca.mean(axis=0).reshape(1, -1)
disease_gene_mean = disease_gene.mean(axis=0).reshape(1, -1)

# ------------------------------------------------------------------
# 4. Test 1: OT Cost Map Monotonicity
# ------------------------------------------------------------------
print("=" * 60)
print("TEST 1: OT Cost Map Monotonicity")
print("=" * 60)

distances = []
for i in range(N_STAGES - 1):
    X1, w1 = stage_dist_pca[stage_order[i]]
    X2, w2 = stage_dist_pca[stage_order[i+1]]
    d, _ = compute_wasserstein_distance(X1, X2, w1, w2, reg=0.05)
    distances.append(d)
    print(f"  {stage_order[i]} -> {stage_order[i+1]}: W2 = {d:.4f}")

# Cumulative from Stage_0
cum_dists = []
for i in range(1, N_STAGES):
    X0, w0 = stage_dist_pca["Stage_0"]
    Xi, wi = stage_dist_pca[f"Stage_{i}"]
    d, _ = emd_w2(X0, w0, Xi, wi)
    cum_dists.append(d)

monotonic = all(cum_dists[i] <= cum_dists[i+1] for i in range(len(cum_dists)-1))
print(f"\n✓ Monotonicity: {'PASS' if monotonic else 'FAIL'}")
for i, d in enumerate(cum_dists, 1):
    print(f"  Stage_0 -> Stage_{i}: {d:.4f}")

# ------------------------------------------------------------------
# 5. Test 2: TDI Accuracy (using EMD for fair comparison)
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 2: TDI Accuracy")
print("=" * 60)

w_disease = np.ones(1)

tdi_dists = {}
for stage in stage_order:
    X_stage, w_stage = stage_dist_pca[stage]
    d, _ = emd_w2(disease_pca_mean, w_disease, X_stage, w_stage)
    tdi_dists[stage] = d

print("\nDisease metacell to normal stages (EMD W2):")
for stage in stage_order:
    marker = " <-- GROUND TRUTH" if stage == "Stage_3" else ""
    print(f"  {stage}: {tdi_dists[stage]:.4f}{marker}")

pred_stage = min(tdi_dists, key=tdi_dists.get)
correct = (pred_stage == "Stage_3")
print(f"\n✓ Predicted arrest: {pred_stage}")
print(f"✓ TDI accuracy: {'PASS' if correct else 'FAIL'}")

# ------------------------------------------------------------------
# 6. Test 3: Repair Directionality
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 3: Repair Pathway Directionality")
print("=" * 60)

X_disease = disease_gene_mean
X_terminal, w_terminal = stage_dist_gene["Stage_4"]

# Sinkhorn for repair (as in real analysis)
dist, transport_plan = compute_wasserstein_distance(
    X_disease, X_terminal, np.ones(1), w_terminal, reg=0.05
)

repair_deltas = compute_gene_repair_delta(
    X_disease, X_terminal, transport_plan, gene_names=gene_names
)

diff_genes = gene_names[:30]
perturb_genes = gene_names[30:60]
neutral_genes = gene_names[60:]

diff_mean = np.mean([repair_deltas[g] for g in diff_genes])
perturb_mean = np.mean([repair_deltas[g] for g in perturb_genes])
neutral_mean = np.mean([repair_deltas[g] for g in neutral_genes])

print(f"\nMean repair delta:")
print(f"  Diff. drivers (genes 0-29):   {diff_mean:+.4f}  (expected: positive)")
print(f"  Disease perturb. (genes 30-59): {perturb_mean:+.4f}  (expected: negative)")
print(f"  Neutral (genes 60+):          {neutral_mean:+.4f}  (expected: ~0)")

direction_ok = (diff_mean > 0.1) and (perturb_mean < -0.1) and (abs(neutral_mean) < 0.1)
print(f"\n✓ Repair directionality: {'PASS' if direction_ok else 'FAIL'}")

# Top targets
top_up = sorted(repair_deltas.items(), key=lambda x: x[1], reverse=True)[:10]
top_down = sorted(repair_deltas.items(), key=lambda x: x[1])[:10]

print(f"\nTop 10 UP targets:")
for g, d in top_up:
    cat = "diff" if g in diff_genes else ("perturb" if g in perturb_genes else "neutral")
    print(f"  {g}: {d:+.4f} [{cat}]")

print(f"\nTop 10 DOWN targets:")
for g, d in top_down:
    cat = "diff" if g in diff_genes else ("perturb" if g in perturb_genes else "neutral")
    print(f"  {g}: {d:+.4f} [{cat}]")

# ------------------------------------------------------------------
# 7. Test 4: Epsilon Sensitivity (repair stability)
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 4: Repair Sensitivity to epsilon")
print("=" * 60)

epsilons = [0.001, 0.01, 0.05, 0.1, 0.5]
repair_scores = []

for eps in epsilons:
    _, tp = compute_wasserstein_distance(
        X_disease, X_terminal, np.ones(1), w_terminal, reg=eps
    )
    rd = compute_gene_repair_delta(X_disease, X_terminal, tp, gene_names=gene_names)
    repair_scores.append({
        "eps": eps,
        "diff": np.mean([rd[g] for g in diff_genes]),
        "perturb": np.mean([rd[g] for g in perturb_genes]),
        "neutral": np.mean([rd[g] for g in neutral_genes]),
    })

print(f"{'eps':>8} {'diff':>10} {'perturb':>10} {'neutral':>10}")
print("-" * 42)
for rs in repair_scores:
    print(f"{rs['eps']:8.3f} {rs['diff']:10.4f} {rs['perturb']:10.4f} {rs['neutral']:10.4f}")

# Check sign consistency
signs_ok = all(
    (rs["diff"] > 0) and (rs["perturb"] < 0) and (abs(rs["neutral"]) < 0.2)
    for rs in repair_scores
)
print(f"\n✓ Repair sign stability: {'PASS' if signs_ok else 'FAIL'}")

# ------------------------------------------------------------------
# 8. Visualization
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("Generating figure...")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: PCA
ax = axes[0, 0]
colors = plt.cm.tab10(np.linspace(0, 1, N_STAGES + 1))
for s in range(N_STAGES):
    mask = adata_normal.obs["cell_type"] == f"Stage_{s}"
    ax.scatter(adata_normal.obsm["X_pca"][mask, 0], adata_normal.obsm["X_pca"][mask, 1],
               label=f"Stage {s}", alpha=0.4, s=8, c=[colors[s]])
ax.scatter(adata_disease.obsm["X_pca"][:, 0], adata_disease.obsm["X_pca"][:, 1],
           label="Disease", c="red", marker="x", s=20, alpha=0.5)
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("A. Simulated Developmental Landscape")
ax.legend(loc="upper left", fontsize=7)

# Panel B: Cost map
ax = axes[0, 1]
pairs = [f"S{i}\u2192S{i+1}" for i in range(N_STAGES - 1)]
bars = ax.bar(pairs, distances, color="steelblue", edgecolor="black")
ax.set_ylabel("W2 Distance (Sinkhorn)")
ax.set_title("B. OT Cost Map")
for bar, val in zip(bars, distances):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
            f"{val:.3f}", ha="center", va="bottom", fontsize=9)

# Panel C: TDI
ax = axes[1, 0]
stages = list(tdi_dists.keys())
vals = [tdi_dists[s] for s in stages]
colors_bar = ["#e74c3c" if s == "Stage_3" else "#3498db" for s in stages]
bars = ax.bar(stages, vals, color=colors_bar, edgecolor="black")
ax.set_ylabel("EMD W2 Distance")
ax.set_title("C. TDI: Disease-to-Normal Distances")
ax.set_xticklabels(stages, rotation=45, ha="right")
for bar, val in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f"{val:.2f}", ha="center", va="bottom", fontsize=9)

# Panel D: Repair by category
ax = axes[1, 1]
cats = ["Diff. drivers\n(genes 0-29)", "Disease pert.\n(genes 30-59)", "Neutral\n(genes 60+)"]
means = [diff_mean, perturb_mean, neutral_mean]
colors_cat = ["#2ecc71", "#e74c3c", "#95a5a6"]
bars = ax.bar(cats, means, color=colors_cat, edgecolor="black")
ax.axhline(y=0, color="black", linewidth=0.5)
ax.set_ylabel("Mean Repair Delta")
ax.set_title("D. Repair Directionality")
for bar, val in zip(bars, means):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02 * np.sign(val),
            f"{val:+.3f}", ha="center", va="bottom" if val > 0 else "top", fontsize=10)

plt.tight_layout()
out_path = "/data1/yja/zhongzhuan/5.external/scTDRP/scripts/simulation_validation_v2.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"\nFigure saved: {out_path}")

# ------------------------------------------------------------------
# 9. Summary
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("VALIDATION SUMMARY")
print("=" * 60)
print(f"  [1] OT Cost Map Monotonicity: {'PASS' if monotonic else 'FAIL'}")
print(f"  [2] TDI Accuracy:              {'PASS' if correct else 'FAIL'}")
print(f"  [3] Repair Directionality:     {'PASS' if direction_ok else 'FAIL'}")
print(f"  [4] Epsilon Stability:         {'PASS' if signs_ok else 'FAIL'}")
print("=" * 60)
