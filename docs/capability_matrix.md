# Capability Matrix: scTDRP vs Competing Methods

## Problem Setting

| Aspect | Description |
|--------|-------------|
| **Input** | Normal reference developmental atlas (multi-stage) + Disease snapshot (single time-point) |
| **Goal** | Identify where disease cells deviate from normal trajectory and infer repair strategies |
| **Key Challenge** | No time-series data for disease; need cross-condition mapping |

## Method Comparison

| Capability | DGE+GSEA | PCA EucDist | GeneSet Score | Pseudotime Dev | Waddington-OT | MOSCOT | scTDRP-NoOT | **scTDRP** |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **Requires disease time-series** | ✗ | ✗ | ✗ | ✗ | ✓ | ✓ | ✗ | ✗ |
| **Arrest stage identification** | ✗ | ~ | ✗ | ~ | ~ | ~ | ✓ | ✓ |
| **Quantitative deviation score** | ✗ | ~ | ~ | ~ | ✗ | ✗ | ✓ | ✓ |
| **Repair direction inference** | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✓ |
| **Module-level therapeutic strategy** | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | ✓ |
| **Single-cell resolution** | ✗ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **Distribution-aware distance** | ✗ | ✗ | ✗ | ✗ | ✓ | ✓ | ✗ | ✓ |

**Legend:** ✓ = fully supported, ~ = partially supported, ✗ = not supported

## Detailed Method Descriptions

### 1. DGE + GSEA (Differential Expression + Gene Set Enrichment)
- **Approach**: Compare disease vs normal cells using DESeq2/edgeR, then perform GSEA on differentially expressed genes.
- **Limitations**: Treats disease and normal as static populations; ignores developmental continuity; cannot infer therapeutic direction (up/down regulation).

### 2. PCA Euclidean Distance (PCA_EucDist)
- **Approach**: Compute Euclidean distance from each disease cell to the nearest normal stage centroid in PCA space.
- **Limitations**: Ignores intra-stage heterogeneity; no repair inference; simplistic distance metric.

### 3. Gene-set Differentiation Score (GeneSet_Score)
- **Approach**: Compute AddModuleScore for known differentiation gene sets (e.g., P7 vs P8) and take their difference.
- **Limitations**: Requires prior knowledge of relevant gene sets; cannot identify arrest stage; no repair inference.

### 4. Pseudotime Deviation (Pseudotime_Dev)
- **Approach**: Project disease cells onto normal pseudotime axis (e.g., PC1) and measure deviation from terminal state.
- **Limitations**: Assumes linear trajectory; sensitive to PCA orientation; no repair inference.

### 5. Waddington-OT (WOT)
- **Approach**: Optimal transport between consecutive time points to reconstruct developmental trajectories.
- **Limitations**: Requires time-series data from the same process; designed for fate prediction, not disease-normal mapping; no repair inference.
- **Adaptation for comparison**: We treated normal stages as pseudo-time points and disease cells as an additional time point, computing transport from the terminal normal stage to disease.

### 6. MOSCOT (TemporalProblem)
- **Approach**: Fused Gromov-Wasserstein optimal transport for multi-condition/time alignment.
- **Limitations**: Requires time-series or multi-condition design; optimized for alignment quality, not deviation quantification; no repair inference.
- **Adaptation for comparison**: We used TemporalProblem on metacells with normal stages as time points and disease as the final time point.

### 7. scTDRP-NoOT (Ablation)
- **Approach**: Replace Wasserstein distance with Euclidean distance to stage centroids; keep all other components (metacells, TDI framework).
- **Purpose**: Isolates the contribution of OT vs simple geometric distance.
- **Limitations**: No repair pathway inference; less robust to distribution shape.

### 8. scTDRP (Full)
- **Approach**: Optimal transport-based framework with four components:
  1. **OT Cost Map**: Wasserstein distances between consecutive normal stages
  2. **TDI**: Minimal Wasserstein distance from each disease cell to all normal stages
  3. **Repair Pathway**: Transport plan from disease to terminal stage → gene-level repair deltas
  4. **Module Strategy**: Aggregation of gene-level deltas into functional module up/down strategies
- **Unique strengths**:
  - Does **not** require disease time-series
  - Outputs **actionable** module-level therapeutic hypotheses
  - Quantifies **both** where and how far disease cells are arrested

## Key Takeaways

1. **scTDRP is the only method that does not require disease time-series data** while still providing quantitative developmental mapping.
2. **scTDRP is the only method that infers repair directions** (which genes/modules to up/down regulate).
3. **OT-based distance (scTDRP) vs Euclidean distance (scTDRP-NoOT)**: In real data with complex heterogeneity, OT captures distribution shape and is more robust to outliers. In simple simulated data with clear linear trajectories, Euclidean distance can perform equally well for stage identification.
4. **WOT and MOSCOT are powerful for their original problems** (trajectory inference, multi-condition alignment) but are not designed for disease-normal developmental mapping and therapeutic inference.
