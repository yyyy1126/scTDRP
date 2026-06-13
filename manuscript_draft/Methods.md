# Methods

## Data preprocessing and quality control

Single-cell RNA-seq datasets were processed using Scanpy [1]. Cells with fewer than 200 detected genes and genes expressed in fewer than 3 cells were removed. Mitochondrial gene percentage and total counts per cell were inspected to exclude low-quality cells. Raw counts were normalized to 10,000 reads per cell and log-transformed (log1p). Highly variable genes (HVGs) were selected using the Seurat v3 method [2] (top 2,000 HVGs unless otherwise specified). Principal component analysis (PCA) was performed on the HVG-subset expression matrix, and the top 50 principal components were retained for downstream optimal transport calculations.

## Normal developmental stage definition

Normal hematopoietic differentiation stages were defined based on established cell-type annotations from the literature and dataset metadata. For the erythroid lineage, cells were classified into seven ordered stages: megakaryocyte-erythroid progenitor (MEP), burst-forming unit-erythroid (BFU-E), colony-forming unit-erythroid (CFU-E), proerythroblast (Pro-E), basophilic erythroblast (Baso-E), polychromatic erythroblast (Poly-E), and orthochromatic erythroblast (Ortho-E). For the megakaryocytic lineage, stages included MEP, megakaryocyte precursor (Mk-Precursor), and megakaryocyte (Mk). Stage order was determined by known hematopoietic differentiation hierarchies [3,4].


## B-ALL dataset and preprocessing

**Data source.** B-cell acute lymphoblastic leukemia (B-ALL) scRNA-seq data were obtained from Witkowski et al. (2020; GEO accession GSE130116). The dataset comprised 24 bone marrow samples from 10 pediatric B-ALL patients (matched diagnosis, remission, and relapse) and 4 healthy bone marrow donors. For the present analysis, we selected 7 diagnostic B-ALL samples and 4 healthy donor samples (Table S2).

**Sample composition.** B-ALL patient bone marrow samples were separated into CD19+ B-cell and CD19−CD45+ non-B-cell fractions by flow-based sorting and then mixed at a ratio of 1:5 (CD19+:CD19−CD45+) prior to scRNA-seq. Total CD45+ bone marrow cells were sequenced from healthy donors.

**Raw data structure.** GSE130116 was generated with the 10x Genomics Chromium v2 platform using a shared barcode whitelist of 737,280 cell barcodes. Each sample was deposited as an independent raw (unfiltered) feature-barcode matrix (33,538 genes × 737,280 barcodes). Because the barcode pool was shared across all samples, we first mapped each barcode to its originating sample by identifying the single sample in which that barcode column contained non-zero counts. Barcodes with non-zero counts in multiple samples or in none were marked as ambiguous/unassigned and excluded.

**Quality control.** After sample-specific barcode extraction, we applied the following QC thresholds: total UMI counts ≥ 200, detected genes ≥ 100, and mitochondrial gene percentage < 20%. These thresholds were chosen to balance retention of genuine cells against exclusion of empty droplets, which are abundant in the deposited raw matrices. After QC, 7,139 B-ALL diagnosis cells and 57 healthy donor cells were retained. Healthy donor cells were too few for independent analysis; the normal B-cell reference trajectory was instead constructed from 16,878 B cells (Pro-B VDJ, Large Pre-B, Immature B, Mature B) curated from the blood map dataset.

**B-cell selection.** Because B-ALL diagnostic samples contained a mixture of CD19+ malignant B cells and CD19− non-B immune microenvironment cells, we performed B-cell selection using marker-gene scoring. For each cell, B-cell score (CD19, MS4A1, CD79A, CD79B, PAX5, VPREB1, RAG1, RAG2), T-cell score (CD3D, CD3E, CD4, CD8A, TRAC), and myeloid score (CD14, LYZ, S100A8, S100A9, CD33) were computed using Scanpy's `score_genes` function after log-normalization. Cells with B-cell score above the median were retained as the malignant B-cell compartment (n = 3,569), representing the leukemic blast population for downstream scTDRP analysis.

## GBM dataset and preprocessing

**Data source.** Glioblastoma (GBM) scRNA-seq data were obtained from Darmanis et al. (2017; GEO accession GSE84465). The dataset profiled 3,589 cells from four primary IDH-wild-type GBM patients and included both neoplastic cells and non-neoplastic tumor microenvironment cells.

**Normal reference.** We constructed a normal neural lineage reference from 579 non-neoplastic cells annotated as oligodendrocyte precursor cells (OPC, n = 406), astrocytes (n = 88), or oligodendrocytes (n = 85). The ordered reference trajectory was defined as OPC → Astrocyte → Oligodendrocyte, with Oligodendrocyte as the terminal stage.

**Disease cells.** Neoplastic cells annotated in the original study (n = 1,091) were used as the disease sample. No additional marker-based selection was applied.

**Note on metacell aggregation.** In the GBM validation, both the normal reference and disease cells were analyzed at single-cell resolution without metacell aggregation. Because the reference contained only three stages and adjacent stages (OPC and Oligodendrocyte) are transcriptionally close, metacell averaging blurred stage boundaries and shifted stage assignments toward the terminal stage. This observation suggests that metacell aggregation should be applied cautiously when the reference trajectory contains few stages or closely spaced cell states.

## Metacell construction

To enable computationally efficient optimal transport analysis while preserving biological heterogeneity, we aggregated cells within each normal stage into metacells using Leiden clustering [5]. Briefly, for each developmental stage, we performed Leiden clustering on the PCA-reduced expression matrix with stage-specific resolution parameters (Table S1). The expression profiles of cells within each cluster were averaged to generate a metacell, and the cluster size (number of constituent cells) was recorded as the metacell weight. For disease cells, Leiden clustering was performed across all malignant cells at resolution 1.0, and cluster averages were computed similarly. Metacell weights were incorporated into optimal transport calculations to ensure that large cell populations exert proportionally greater influence on the transport plan.

This metacell strategy serves two purposes. First, it reduces the computational complexity of Sinkhorn iterations from O(n × m) to O(n' × m'), where n' and m' are the numbers of metacells, enabling analysis of datasets with thousands of cells. Second, it provides a natural regularization by averaging out technical noise within biologically homogeneous subpopulations while preserving inter-subpopulation variance.

## Optimal transport cost map construction

We modeled each normal developmental stage as a probability distribution over gene expression space. The expression distribution of stage $s$ was represented by its metacell expression matrix $X_s \in \mathbb{R}^{n_s \times d}$ and weight vector $w_s \in \mathbb{R}^{n_s}$, where $n_s$ is the number of metacells in stage $s$ and $d$ is the number of features (PCA dimensions or selected genes).

For any two consecutive stages $s_i$ and $s_{i+1}$, we computed the entropically regularized optimal transport plan using the Sinkhorn algorithm [6,7] as implemented in the Python Optimal Transport (POT) library [8]. Given cost matrix $C \in \mathbb{R}^{n_i \times n_{i+1}}$ where $C_{jk} = \|x_j - y_k\|_2^2$ is the squared Euclidean distance between metacell $j$ in stage $s_i$ and metacell $k$ in stage $s_{i+1}$, the optimal transport plan $\Gamma^*$ is obtained by solving:

$$\Gamma^* = \arg\min_{\Gamma \in \Pi(w_i, w_{i+1})} \langle \Gamma, C \rangle + \varepsilon \cdot \text{KL}(\Gamma \| w_i \otimes w_{i+1})$$

where $\Pi(w_i, w_{i+1})$ denotes the set of couplings with marginal distributions $w_i$ and $w_{i+1}$, KL is the Kullback-Leibler divergence, and $\varepsilon > 0$ is the entropic regularization parameter (default $\varepsilon = 0.05 \cdot \text{median}(C)$).

The Wasserstein-2 distance between stages $s_i$ and $s_{i+1}$ is then:

$$W_2(s_i, s_{i+1}) = \sqrt{\langle \Gamma^*, C \rangle}$$

By iterating over all adjacent stage pairs, we constructed the **OT cost map** $M \in \mathbb{R}^{S \times S}$, where $S$ is the total number of stages and $M_{ij} = W_2(s_i, s_j)$ for consecutive stages (non-consecutive entries were left undefined or computed by path summation).

## Transport-based Developmental Index (TDI)

To determine the developmental position of each disease metacell relative to the normal trajectory, we computed the Transport-based Developmental Index (TDI) as the minimal Wasserstein distance from the disease metacell to any normal stage:

$$\text{TDI}(c) = \min_{s \in \{1, \ldots, S\}} W_2(c, s)$$

where $W_2(c, s)$ is the Wasserstein-2 distance between disease metacell $c$ and normal stage $s$, computed via entropically regularized optimal transport with the same regularization parameter as above. The **attributed stage** $s^*(c)$ was defined as the normal stage achieving this minimum:

$$s^*(c) = \arg\min_{s \in \{1, \ldots, S\}} W_2(c, s)$$

TDI values were computed in the PCA-reduced space (50 dimensions) to balance biological signal preservation with computational efficiency. For visualization and interpretation, TDI values were normalized to a [0, 1] scale, where 0 indicates identity with a normal stage and 1 indicates maximal developmental deviation.

## Repair pathway inference

The core therapeutic insight of scTDRP comes from analyzing the transport plan between disease cells and the terminal normal stage. For each disease metacell $c$ with expression vector $x_c \in \mathbb{R}^d$ and the terminal stage $t$ with metacell matrix $X_t \in \mathbb{R}^{n_t \times d}$, we computed the optimal transport plan $\Gamma^{c \to t}$ and extracted the gene-level repair quantity.

For each gene $g$ (in the original expression space, not PCA), the repair quantity $\Delta_g$ was computed as the expected expression difference between the terminal stage and the disease cell, weighted by the transport plan:

$$\Delta_g = \sum_{j=1}^{n_c} \sum_{k=1}^{n_t} \Gamma^{c \to t}_{jk} \cdot (X_{t,kg} - X_{c,jg})$$

where $X_{c,jg}$ is the expression of gene $g$ in the $j$-th constituent metacell of disease cluster $c$, and $X_{t,kg}$ is the expression of gene $g$ in the $k$-th metacell of the terminal stage. Positive $\Delta_g$ indicates that gene $g$ should be upregulated to push the disease cell toward the terminal stage, while negative $\Delta_g$ indicates downregulation.

This calculation was performed on the module-gene subspace (genes belonging to prior functional modules) to ensure biological interpretability. Genes were ranked by $\Delta_g$ magnitude to identify the strongest repair targets.

## Module-level repair strategy

To translate gene-level repair quantities into actionable therapeutic strategies, we aggregated $\Delta_g$ values by prior functional modules extracted from consensus non-negative matrix factorization (cNMF) [9,10] of the normal differentiation data. cNMF was run on the normal cell expression matrix with parameters $k = 8$ and local density threshold $\lambda = 2.0$. For each module $m$, the module repair score was computed as:

$$\text{Score}(m) = \frac{1}{|G_m|} \sum_{g \in G_m} \Delta_g$$

where $G_m$ is the set of genes belonging to module $m$. Modules with positive scores were flagged for upregulation, and modules with negative scores were flagged for downregulation. Statistical significance was assessed by permutation testing: gene labels were randomly shuffled 1,000 times, and the empirical p-value was computed as the fraction of permutations yielding a module score at least as extreme as the observed value.

## Statistical analysis and validation

All optimal transport computations were performed using POT v0.9.6 [8]. Statistical significance of stage attribution was assessed using permutation tests (1,000 permutations). For survival-relevant biomarkers identified in the erythroid lineage, Kaplan-Meier survival curves were compared using the log-rank test. All analyses were performed in Python 3.10 with Scanpy v1.9, AnnData v0.8, NumPy v1.23, and SciPy v1.10.

## Comparison with WaddingtonOT and Moscot

For the GBM cross-system validation, we additionally ran two independent optimal-transport frameworks on the identical PCA embedding and cell labels to benchmark scTDRP's numerical consistency.

**WaddingtonOT.** We used WOT v1.0 with `day_field="day"`, `local_pca=0` (so that the provided 50-dimensional PCA embedding was used directly), `epsilon=0.05`, `lambda1=1`, and `lambda2=50`. Disease cells were assigned day 0 and terminal-stage oligodendrocytes day 1; the resulting transport map was row-normalized to infer the expected terminal PCA coordinate for each disease cell, from which the repair direction was computed as the difference to the disease cell's original PCA coordinate. For stage assignment, all normal cells (OPC, Astrocyte, Oligodendrocyte) were assigned day 1 and the transport mass to each stage was summed per disease cell; the stage with the largest mass was taken as the WOT best-match stage.

**Moscot.** We used moscot v0.5.0 `TemporalProblem` with `time_key="day"` and `joint_attr="X_pca"`, solved with `epsilon=0.01`. Cell transitions from day 0 (disease) to day 1 (normal) were aggregated by stage to obtain stage probabilities. Repair direction was obtained by pulling the terminal-stage PCA coordinates back onto disease cells through the solved OT coupling.

## Software availability

scTDRP is implemented as a Python package with modular components for distance computation, repair pathway inference, and visualization. The source code, analysis notebooks, and processed data are available at https://github.com/yyyy1126/scTDRP and archived at Zenodo (DOI: https://doi.org/10.5281/zenodo.15640823).

## References

[1] Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19(1):15.
[2] Stuart T, et al. Comprehensive integration of single-cell data. Cell. 2019;177(7):1888-1902.
[3] Pellin D, et al. A comprehensive single cell transcriptional landscape of human hematopoietic progenitors. Nat Commun. 2019;10:2395.
[4] Jassinskaja M, et al. Deep phenotyping of human bone marrow explains normal blood formation and leukemia. bioRxiv. 2023.
[5] Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing well-connected communities. Sci Rep. 2019;9(1):5233.
[6] Cuturi M. Sinkhorn distances: Lightspeed computation of optimal transport. NeurIPS. 2013;26.
[7] Peyré G, Cuturi M. Computational optimal transport. Found Trends Mach Learn. 2019;11(5-6):355-607.
[8] Flamary R, et al. POT: Python Optimal Transport. J Mach Learn Res. 2021;22(78):1-8.
[9] Kotliar D, et al. Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq. eLife. 2019;8:e43803.
[10] Lareau C, et al. The dynamics of expression programs in normal hematopoiesis and acute myeloid leukemia. bioRxiv. 2023.
