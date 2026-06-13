#!/usr/bin/env python3
"""
scTDRP Concept Validation Simulation v3
========================================
Low-dimensional controlled simulation in a known 2D subspace:
  - X-axis: differentiation axis (5 stages at x=0,2,4,6,8)
  - Y-axis: disease perturbation axis (disease cells at y=+3)
No PCA ambiguity; directly validates OT logic.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ot

np.random.seed(42)

N_STAGES = 5
N_CELLS = 150
NOISE = 0.3
OFFSET = 2.5

# ------------------------------------------------------------------
# 1. Generate data in 2D: (differentiation, perturbation)
# ------------------------------------------------------------------
stage_names = [f"S{i}" for i in range(N_STAGES)]

normal_data = []
labels = []
for s in range(N_STAGES):
    x = np.random.normal(s * 2.0, NOISE, size=N_CELLS)
    y = np.random.normal(0, NOISE, size=N_CELLS)
    normal_data.append(np.column_stack([x, y]))
    labels.extend([stage_names[s]] * N_CELLS)

normal_data = np.vstack(normal_data)  # (750, 2)

# Disease: arrested at S3 (x=6) + perturbation in y-direction
disease_x = np.random.normal(6.0, NOISE, size=N_CELLS)
disease_y = np.random.normal(OFFSET, NOISE, size=N_CELLS)
disease_data = np.column_stack([disease_x, disease_y])
disease_label = ["Disease"] * N_CELLS

# ------------------------------------------------------------------
# 2. Stage distributions (uniform weights)
# ------------------------------------------------------------------
stage_dists = {}
for s in range(N_STAGES):
    mask = np.array(labels) == stage_names[s]
    X = normal_data[mask]
    w = np.ones(X.shape[0]) / X.shape[0]
    stage_dists[stage_names[s]] = (X, w)

disease_w = np.ones(N_CELLS) / N_CELLS

# ------------------------------------------------------------------
# 3. Test 1: OT Cost Map (adjacent stages)
# ------------------------------------------------------------------
print("=" * 60)
print("TEST 1: OT Cost Map Monotonicity")
print("=" * 60)

adj_costs = []
for i in range(N_STAGES - 1):
    X1, w1 = stage_dists[stage_names[i]]
    X2, w2 = stage_dists[stage_names[i+1]]
    M = ot.dist(X1, X2, metric="sqeuclidean")
    plan = ot.emd(w1, w2, M)
    cost = float(np.sum(plan * M))
    adj_costs.append(cost)
    print(f"  {stage_names[i]} -> {stage_names[i+1]}: W2 = {cost:.4f}")

# Cumulative S0 -> Si
cum_costs = []
for i in range(1, N_STAGES):
    X0, w0 = stage_dists["S0"]
    Xi, wi = stage_dists[stage_names[i]]
    M = ot.dist(X0, Xi, metric="sqeuclidean")
    plan = ot.emd(w0, wi, M)
    cost = float(np.sum(plan * M))
    cum_costs.append(cost)

mono = all(cum_costs[i] <= cum_costs[i+1] for i in range(len(cum_costs)-1))
print(f"\nCumulative S0->Si: {[f'{c:.2f}' for c in cum_costs]}")
print(f"✓ Monotonicity: {'PASS' if mono else 'FAIL'}")

# ------------------------------------------------------------------
# 4. Test 2: TDI (disease vs each stage)
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 2: TDI Accuracy")
print("=" * 60)

tdi = {}
for s in range(N_STAGES):
    Xs, ws = stage_dists[stage_names[s]]
    M = ot.dist(disease_data, Xs, metric="sqeuclidean")
    plan = ot.emd(disease_w, ws, M)
    cost = float(np.sum(plan * M))
    tdi[stage_names[s]] = cost
    marker = " <-- GROUND TRUTH (S3)" if s == 3 else ""
    print(f"  Disease -> S{s}: W2 = {cost:.4f}{marker}")

pred = min(tdi, key=tdi.get)
tdi_ok = (pred == "S3")
print(f"\n✓ Predicted: {pred} | Ground truth: S3 | {'PASS' if tdi_ok else 'FAIL'}")

# ------------------------------------------------------------------
# 5. Test 3: Cell-level TDI (sample 20 disease cells)
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 3: Cell-level TDI (20 sampled cells)")
print("=" * 60)

cell_tdi = []
for idx in np.random.choice(N_CELLS, 20, replace=False):
    Xc = disease_data[idx:idx+1]
    wc = np.ones(1)
    dists = {}
    for s in range(N_STAGES):
        Xs, ws = stage_dists[stage_names[s]]
        M = ot.dist(Xc, Xs, metric="sqeuclidean")
        plan = ot.emd(wc, ws, M)
        dists[stage_names[s]] = float(np.sum(plan * M))
    cell_tdi.append(min(dists, key=dists.get))

from collections import Counter
cnt = Counter(cell_tdi)
print(f"  Attribution counts: {dict(cnt)}")
cell_ok = cnt["S3"] >= 15  # at least 15/20 correct
print(f"✓ Cell-level accuracy: {cnt['S3']}/20 -> {'PASS' if cell_ok else 'FAIL'}")

# ------------------------------------------------------------------
# 6. Test 4: Repair pathway (disease -> terminal S4)
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 4: Repair Pathway Directionality")
print("=" * 60)

# We'll treat x-coordinate as "differentiation driver" and y as "disease perturbation"
# Terminal = S4
X_term, w_term = stage_dists["S4"]
M = ot.dist(disease_data, X_term, metric="sqeuclidean")
plan = ot.emd(disease_w, w_term, M)

# Gene-level repair delta: here each dimension is a "gene"
# delta_x = expected change in x to reach terminal
# delta_y = expected change in y to reach terminal
delta = np.zeros(2)
for i in range(disease_data.shape[0]):
    for j in range(X_term.shape[0]):
        delta += plan[i, j] * (X_term[j] - disease_data[i])

print(f"\nMean repair delta (transport-weighted):")
print(f"  X (differentiation driver): {delta[0]:+.4f}  (expected: positive, push toward S4)")
print(f"  Y (disease perturbation):   {delta[1]:+.4f}  (expected: negative, remove offset)")

repair_ok = (delta[0] > 0.5) and (delta[1] < -1.0)
print(f"✓ Directionality: {'PASS' if repair_ok else 'FAIL'}")

# Per-cell repair deltas
repair_x = []
repair_y = []
for i in range(disease_data.shape[0]):
    dx = 0.0
    dy = 0.0
    for j in range(X_term.shape[0]):
        dx += plan[i, j] * (X_term[j, 0] - disease_data[i, 0])
        dy += plan[i, j] * (X_term[j, 1] - disease_data[i, 1])
    repair_x.append(dx)
    repair_y.append(dy)

print(f"\nPer-cell repair stats:")
print(f"  dX: mean={np.mean(repair_x):+.3f}, std={np.std(repair_x):.3f}")
print(f"  dY: mean={np.mean(repair_y):+.3f}, std={np.std(repair_y):.3f}")

# ------------------------------------------------------------------
# 7. Test 5: Epsilon sensitivity (Sinkhorn on S3->S4)
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("TEST 5: Epsilon Sensitivity (Sinkhorn)")
print("=" * 60)

X3, w3 = stage_dists["S3"]
X4, w4 = stage_dists["S4"]
M = ot.dist(X3, X4, metric="sqeuclidean")
M /= M.max()

eps_vals = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]
sink_costs = []
for eps in eps_vals:
    try:
        p = ot.sinkhorn(w3, w4, M, eps)
        c = float(np.sum(p * M))
        sink_costs.append(c)
    except Exception as e:
        sink_costs.append(None)
        print(f"  eps={eps}: FAILED ({e})")
        continue
    print(f"  eps={eps:5.3f}: W2 = {c:.4f}")

# Check relative stability (ratio to EMD ground truth)
emd_plan = ot.emd(w3, w4, M)
emd_cost = float(np.sum(emd_plan * M))
print(f"  EMD (ground truth):  W2 = {emd_cost:.4f}")

stable = all(s is not None and 0.5 <= s/emd_cost <= 2.0 for s in sink_costs)
print(f"✓ Stability within 2x of EMD: {'PASS' if stable else 'FAIL'}")

# ------------------------------------------------------------------
# 8. Visualization
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("Generating figure...")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(11, 10))

# Panel A: landscape
ax = axes[0, 0]
colors = plt.cm.tab10(np.linspace(0, 1, N_STAGES))
for s in range(N_STAGES):
    mask = np.array(labels) == stage_names[s]
    ax.scatter(normal_data[mask, 0], normal_data[mask, 1],
               label=f"Stage {s}", alpha=0.4, s=10, c=[colors[s]])
ax.scatter(disease_data[:, 0], disease_data[:, 1],
           label="Disease", c="red", marker="x", s=15, alpha=0.4)
ax.set_xlabel("Differentiation axis (X)")
ax.set_ylabel("Perturbation axis (Y)")
ax.set_title("A. Simulated Landscape")
ax.legend(loc="upper left", fontsize=8)
ax.axhline(y=0, color="gray", linewidth=0.5)

# Panel B: cost map
ax = axes[0, 1]
pairs = [f"S{i}\u2192S{i+1}" for i in range(N_STAGES - 1)]
bars = ax.bar(pairs, adj_costs, color="steelblue", edgecolor="black")
ax.set_ylabel("W2 Distance")
ax.set_title("B. OT Cost Map")
for bar, val in zip(bars, adj_costs):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
            f"{val:.2f}", ha="center", va="bottom", fontsize=9)

# Panel C: TDI
ax = axes[1, 0]
vals = [tdi[stage_names[s]] for s in range(N_STAGES)]
colors_bar = ["#e74c3c" if s == 3 else "#3498db" for s in range(N_STAGES)]
bars = ax.bar(stage_names, vals, color=colors_bar, edgecolor="black")
ax.set_ylabel("W2 Distance (Disease -> Stage)")
ax.set_title("C. TDI")
for bar, val in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
            f"{val:.2f}", ha="center", va="bottom", fontsize=9)

# Panel D: epsilon sensitivity
ax = axes[1, 1]
valid = [(e, c) for e, c in zip(eps_vals, sink_costs) if c is not None]
if valid:
    ex, cy = zip(*valid)
    ax.plot(ex, cy, "o-", label="Sinkhorn", color="steelblue")
    ax.axhline(y=emd_cost, color="red", linestyle="--", label="EMD (ground truth)")
    ax.set_xlabel("Regularization epsilon")
    ax.set_ylabel("W2 Distance")
    ax.set_title("D. Sinkhorn Sensitivity (S3->S4)")
    ax.set_xscale("log")
    ax.legend()

plt.tight_layout()
out = "/data1/yja/zhongzhuan/5.external/scTDRP/scripts/simulation_validation_v3.png"
plt.savefig(out, dpi=300, bbox_inches="tight")
print(f"Figure saved: {out}")

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
print("\n" + "=" * 60)
print("VALIDATION SUMMARY")
print("=" * 60)
print(f"  [1] Cost Map Monotonicity:    {'PASS' if mono else 'FAIL'}")
print(f"  [2] TDI Accuracy (group):     {'PASS' if tdi_ok else 'FAIL'}")
print(f"  [3] TDI Accuracy (cell):      {'PASS' if cell_ok else 'FAIL'}")
print(f"  [4] Repair Directionality:    {'PASS' if repair_ok else 'FAIL'}")
print(f"  [5] Epsilon Stability:        {'PASS' if stable else 'FAIL'}")
print("=" * 60)
