# Figure Legends

## Supplementary Figure 1 | Concept validation of scTDRP on simulated data

**(A)** Simulated developmental landscape with five normal differentiation stages (Stage 0–4, colored dots) distributed along the X-axis (differentiation driver) and disease cells (red crosses) arrested at Stage 3 with an orthogonal perturbation along the Y-axis (disease perturbation). Each stage contains 150 cells with Gaussian noise (σ = 0.3). Disease cells are displaced by 2.5 units in the perturbation direction.

**(B)** Optimal transport cost map showing Wasserstein-2 distances between adjacent normal stages. Distances are consistent (~3.9–4.2) across all stage transitions, reflecting uniform step size in the simulated differentiation axis.

**(C)** Transport-based Developmental Index (TDI) computed as the Wasserstein-2 distance from the disease cell population to each normal stage. The minimal distance is correctly identified at Stage 3 (red bar, 6.27), the ground-truth arrest stage. The symmetric distances to Stage 2 and Stage 4 (~10.13 each) reflect the geometric position of disease cells equidistant from the two neighboring stages in the differentiation axis.

**(D)** Sinkhorn regularization sensitivity. Wasserstein-2 distances between Stage 3 and Stage 4 computed by Sinkhorn algorithm with varying entropic regularization parameters ε (blue curve), compared to the exact EMD solution (red dashed line). Sinkhorn distances remain within 10% of the ground truth across three orders of magnitude of ε (0.001–1.0), demonstrating robustness to parameter choice.

## Figure 4 | scTDRP Validation in B-cell Acute Lymphoblastic Leukemia (B-ALL)

**(A)** TDI comparison of all normal B cells (n = 16,878, green) versus B-ALL malignant B cells (n = 3,569, red). Boxes show median and interquartile range; whiskers show 1.5× IQR. The overall difference is not significant (Mann-Whitney U, p = 0.30) because B-ALL cells map to Large Pre-B, the normal stage with intrinsically lowest TDI.

**(B)** Proportion of B-ALL malignant cells assigned to each normal B-cell developmental stage. The vast majority (82.5%) arrest at the Large Pre-B stage, consistent with the known biology of B-ALL precursor accumulation at the pre-B receptor checkpoint.

**(C)** TDI distribution of normal B cells stratified by their true annotated stage (Pro-B VDJ, Large Pre-B, Immature B, Mature B). Large Pre-B has the intrinsically lowest TDI (mean = 0.173), explaining why B-ALL cells—despite being transcriptionally deviant—do not appear globally elevated when compared against the full normal pool.

**(D)** Density distribution of single-cell TDI comparing B-ALL malignant cells (n = 3,569, red) to normal Large Pre-B cells (n = 3,488, green), the stage-matched reference. B-ALL cells show a markedly right-shifted distribution (Mann-Whitney U p < 1 × 10⁻³⁰⁰; Cohen's d = 1.91), demonstrating significant disease-specific deviation within the matched normal stage.
