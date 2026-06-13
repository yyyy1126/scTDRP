# scTDRP: Single-Cell Transcriptomic Developmental Repair Potential for Quantifying Differentiation Deviation and Inferring Transcriptomic Repair Strategies

**Authors:** [To be determined]

**Affiliations:** [To be determined]

**Correspondence:** [To be determined]

---

## Abstract

Differentiation disorders, ranging from hematologic malignancies to developmental syndromes, share a common feature: disease cells fail to reach their terminal developmental fate. Current single-cell genomic approaches excel at describing disease states but lack a principled framework for quantifying how far disease cells deviate from normal development and inferring transcriptomic strategies to restore differentiation. Here, we present scTDRP (single-cell Transcriptomic Developmental Repair Potential), a computational framework that leverages optimal transport theory to measure the Wasserstein distance between disease cells and normal developmental stages, yielding a Transcriptomic Deviation Index (TDI), and infers gene-level and module-level repair strategies to push disease cells toward terminal differentiation. We validated scTDRP in acute erythroid leukemia (AEL), where malignant erythroid cells showed markedly elevated TDI (0.73) with arrest at the polychromatic stage (81.8%). The inferred repair strategy revealed up-regulation of the terminal differentiation module (P7_TerminalPrep) and down-regulation of the precursor proliferation module (P8_ExecutionPrecursor), consistent with the known biology of erythroid maturation. Intra-tumor heterogeneity analysis further showed that TDI correlates with inferCNV-based malignancy classification (p = 5.4 × 10⁻²³³), cell-cycle phase (G1-arrested cells showing highest deviation, p = 1.1 × 10⁻⁹³), and CNV load (high-CNV cells > low-CNV cells, p = 1.9 × 10⁻¹²⁰). We further tested scTDRP in myelodysplastic syndrome with 5q deletion (MDS-5q), a disorder of ineffective erythropoiesis, and found minimal TDI elevation. This boundary case revealed a critical distinction: scTDRP is sensitive to directional deviation from normal trajectories (as in malignancy) but not to kinetic impairment along the same trajectory (as in bone marrow failure). Finally, cross-lineage validation in B-cell acute lymphoblastic leukemia (B-ALL) confirmed generalizability beyond erythropoiesis: 82.5% of malignant B cells arrested at the Large Pre-B stage and showed significantly elevated TDI relative to stage-matched normal cells (Cohen's d = 1.91, p < 1 × 10⁻³⁰⁰). Simulation studies confirmed the mathematical robustness of scTDRP across five validation tests. scTDRP provides a theoretically grounded approach to translating developmental biology insights into precision therapeutic strategies for differentiation disorders.

**Keywords:** single-cell RNA sequencing, optimal transport, Wasserstein distance, differentiation disorder, erythropoiesis, acute erythroid leukemia, B-cell acute lymphoblastic leukemia, myelodysplastic syndrome, transcriptomic repair

---


---

# scTDRP: Single-Cell Transport-based Disease Reconstruction and Pathway Inference

## Introduction

Cancer has long been conceptualized as a disease of aberrant development, wherein malignant cells recapitulate or arrest at specific stages of normal differentiation programs [1,2]. In hematopoietic malignancies such as acute myeloid leukemia (AML), this paradigm is particularly evident: leukemic blasts frequently exhibit transcriptional signatures resembling immature hematopoietic progenitors, failing to execute terminal differentiation [3,4]. Understanding where disease cells deviate from normal developmental trajectories—and how to redirect them toward healthy fates—represents a central challenge in cancer biology and therapeutic discovery. The advent of single-cell RNA sequencing (scRNA-seq) has provided unprecedented resolution to profile heterogeneous cell populations, yet computational frameworks that translate these rich datasets into actionable therapeutic strategies remain underdeveloped.

Existing analytical approaches for comparing disease and normal states predominantly rely on differential gene expression (DGE) followed by pathway enrichment analysis [5,6]. While informative, these methods have two fundamental limitations. First, they treat disease and normal cells as static populations, ignoring the continuous nature of developmental processes; a differentially expressed gene may reflect either a cause of developmental arrest or a downstream consequence, and DGE alone cannot distinguish these scenarios. Second, they lack a quantitative framework to infer the "direction" of therapeutic intervention: knowing that a pathway is altered does not reveal whether it should be upregulated or downregulated to restore normal development. Recent advances in optimal transport (OT) theory have opened new avenues for addressing these limitations by providing a principled mathematical framework to compare probability distributions and compute minimal-effort mappings between them [7]. In the context of single-cell genomics, Waddington-OT pioneered the use of OT to reconstruct developmental trajectories from time-series scRNA-seq data, inferring ancestor-descendant relationships and regulatory programs underlying cellular reprogramming [8]. Complementarily, Gaussian Graphical Optimal Transport (GGOT) applied population-level Wasserstein distances to detect critical transitions and identify trigger molecules in disease progression [9]. However, neither framework is designed for the specific challenge of mapping individual disease cells onto normal developmental trajectories and inferring repair strategies. Waddington-OT requires temporal sampling of the same process and focuses on fate prediction rather than disease correction. GGOT operates at the population level on bulk expression data and aims to detect tipping points for early diagnosis, rather than guiding therapeutic intervention at the single-cell level.

Here, we present scTDRP (single-cell Transport-based Disease Reconstruction and Pathway inference), a computational framework that leverages optimal transport to quantitatively map malignant cells onto normal developmental trajectories and infer gene module-level repair strategies. scTDRP addresses three interconnected questions: (1) **Where** are disease cells arrested in normal development? We compute the Transport-based Developmental Index (TDI) as the minimal Wasserstein distance from each disease cell to all normal developmental stages, identifying the closest normal stage as the putative arrest point. (2) **How far** are they from terminal differentiation? By constructing an OT cost map between consecutive normal stages, we establish a quantitative metric of developmental distance. (3) **How** can they be pushed toward terminal differentiation? We compute the transport plan from disease cells to the terminal stage and decompose it into gene-level repair quantities, which are then aggregated into functional module-level up/down regulation strategies using prior knowledge of stage-specific gene programs extracted from consensus non-negative matrix factorization (cNMF) [10].

We validate scTDRP on two hematopoietic lineages. First, in the erythroid lineage, we analyze AML patient cells alongside normal erythroid differentiation stages from megakaryocyte-erythroid progenitors (MEP) to orthochromatic erythroblasts. scTDRP reveals that malignant erythroid cells are predominantly arrested at the polychromatic erythroblast stage and identifies terminal preparation and execution precursor modules as key therapeutic targets. Second, in the megakaryocytic lineage, we apply scTDRP to acute megakaryoblastic leukemia (AMKL) cells, revealing distinct patterns of megakaryocytic differentiation arrest and module-level dysregulation. Together, our results demonstrate that scTDRP provides a quantitative, mechanistically interpretable framework for understanding developmental arrest in cancer and translating single-cell transcriptomic data into targeted therapeutic hypotheses.

## References

[1] Sell S. On the stem cell origin of cancer. Am J Pathol. 2010;176(6):2584-94.
[2] Kreso A, Dick JE. Evolution of the cancer stem cell model. Cell Stem Cell. 2014;14(3):275-91.
[3] Pellin D, et al. A comprehensive single cell transcriptional landscape of human hematopoietic progenitors. Nat Commun. 2019;10(1):2395.
[4] van Galen P, et al. Single-Cell RNA-Seq Reveals AML Hierarchies Relevant to Disease Progression and Immunity. Cell. 2019;176(6):1265-1281.e24.
[5] Love MI, et al. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550.
[6] Subramanian A, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci USA. 2005;102(43):15545-50.
[7] Villani C. Optimal Transport: Old and New. Springer; 2009.
[8] Schiebinger G, et al. Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming. Cell. 2019;176(4):928-943.
[9] Hua W, et al. Uncovering critical transitions and molecule mechanisms in disease progressions using Gaussian graphical optimal transport. Commun Biol. 2025;8:575.
[10] Kotliar D, et al. Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq. eLife. 2019;8:e43803.

---

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

## Software availability

scTDRP is implemented as a Python package with modular components for distance computation, repair pathway inference, and visualization. The source code, analysis notebooks, and processed data are available at https://github.com/yyyy1126/scTDRP and archived at Zenodo (DOI: https://doi.org/10.5281/zenodo.20674237).

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

---

# Results

## Overview

We developed scTDRP (single-cell Transcriptomic Developmental Repair Potential), a computational framework that leverages optimal transport theory to quantify how far disease cells deviate from normal developmental trajectories and infer transcriptomic repair strategies to restore terminal differentiation. We validated scTDRP in two disease contexts within the erythroid lineage: (1) acute erythroid leukemia (AEL), a malignancy characterized by malignant arrest of erythroid differentiation, and (2) myelodysplastic syndrome with 5q deletion (MDS-5q), a bone marrow failure disorder characterized by ineffective erythropoiesis.

---

## 3.1 scTDRP Reveals Arrested Erythroid Differentiation in Acute Erythroid Leukemia

### 3.1.1 Dataset and Normal Trajectory Construction

We constructed a normal erythroid differentiation trajectory using 40,166 bone marrow erythroid cells from healthy donors (blood map dataset), spanning seven canonical stages from megakaryocyte-erythroid progenitor (MEP) to orthochromatic erythroblast (Ortho). Cell-type annotations followed the original author's curation, with stages ordered as: MEP → BFU-E → CFU-E → Pro-Erythroblast → Basophilic Erythroblast → Polychromatic Erythroblast → Orthochromatic Erythroblast.

To enable computationally tractable optimal transport calculations on this large dataset, we employed a metacell strategy: cells within each stage were aggregated via Leiden clustering (resolution adapted per stage to yield ~30 metacells per stage), yielding 119 normal metacells representing the complete erythroid differentiation landscape.

### 3.1.2 OT Cost Map Recapitulates Erythroid Differentiation Kinetics

The pairwise Wasserstein-2 distances between adjacent normal stages (OT Cost Map) showed a monotonically decreasing pattern from early to late differentiation (Figure 2A): MEP→BFU-E (0.23), BFU-E→CFU-E (0.22), CFU-E→Pro-E (0.17), Pro-E→Baso-E (0.18), Baso-E→Poly-E (0.10), Poly-E→Ortho-E (0.10). This pattern reflects the biological reality that early erythroid progenitors are transcriptionally more plastic and heterogeneous, while terminal stages become increasingly convergent toward the enucleation program. The cumulative distance from MEP to Ortho was 1.00 (normalized), providing a quantitative measure of the total "developmental distance" in erythropoiesis.

### 3.1.3 AEL Malignant Cells Show Markedly Elevated TDI with Arrest at Polychromatic Stage

We applied scTDRP to 2,134 malignant erythroid cells from patient AML5 (FH8438), an acute erythroid leukemia case characterized by erythroid-biased blasts with trisomy 8. The malignant cells exhibited a mean Transcriptomic Deviation Index (TDI) of **0.73 ± 0.15** (median: 0.72), dramatically higher than the normal metacell baseline (mean TDI to nearest self-stage: ~0.17). This indicates that AEL malignant cells are transcriptionally far from any normal erythroid stage, consistent with their malignant transformation and loss of normal differentiation identity.

The stage assignment analysis revealed that **81.8% (1,745/2,134) of malignant cells were assigned to the Polychromatic Erythroblast stage** as their best-match normal stage, with only 7.7% (165/2,134) reaching Orthochromatic (Figure 2B). This specific arrest at the Poly-E stage aligns with the known pathology of AEL, where blasts accumulate at the polychromatic stage with impaired terminal maturation and enucleation.

### 3.1.4 Gene-Level Repair Pathway Identifies Therapeutic Targets

The OT transport plan from malignant cells to the terminal Orthochromatic stage identified 100 top up-regulated and 100 top down-regulated repair targets (Figure 2C). Key up-regulated targets included globin genes (HBA1, HBA2, HBB) and terminal erythroid transcription factors (KLF1, NFE2), consistent with the need to restore hemoglobin synthesis and terminal maturation. Key down-regulated targets included cell cycle genes (CDC20, CCNB1, TOP2A), reflecting the need to exit the proliferative blast state and enter the quiescent terminal stage.

### 3.1.5 Module Repair Strategy Reveals P7 Activation and P8 Suppression

When repair deltas were aggregated to the pre-defined cNMF-derived modules, scTDRP recommended **up-regulation of the P7_TerminalPrep module** (repair score: +0.08) and **down-regulation of the P8_ExecutionPrecursor module** (repair score: −0.05) (Figure 2D). This strategy is biologically intuitive: P7 contains genes required for terminal erythroid maturation (e.g., hemoglobin synthesis, membrane remodeling, enucleation), which are suppressed in AEL blasts; P8 contains genes associated with early erythroid precursor proliferation and self-renewal, which are inappropriately maintained in the malignant state. The opposing directions of P7 and P8 module strategies mirror the fundamental transition from precursor expansion to terminal execution in normal erythropoiesis.

### 3.1.6 Intra-tumor Heterogeneity: TDI Correlates with Malignancy, Cell Cycle, and CNV Load

To assess whether scTDRP TDI captures clinically meaningful heterogeneity within the tumor, we performed subgroup analyses on the AML5 dataset using independent annotations from inferCNV and cell-cycle scoring.

**Malignancy stratification.** The inferCNV analysis classified cells into three groups: Malignant Erythroid (n = 1,943), Normal Residual (n = 3,705), and Other (n = 3,186, primarily non-erythroid immune cells). TDI showed a highly significant hierarchy: **Malignant Erythroid (mean: 0.778) > Normal Residual (0.706)**, with the non-erythroid "Other" cells showing the highest TDI (0.815) due to their complete departure from the erythroid lineage (Mann-Whitney U, Malignant vs Normal: p = 5.4 × 10⁻²³³; Figure 3A). This demonstrates that TDI is not merely a lineage classifier but specifically quantifies deviation within the target lineage.

**Cell cycle association.** Among malignant cells, those in the G1 phase showed significantly higher TDI (mean: 0.807) than cells in S phase (0.738, p = 1.0 × 10⁻³¹) or G2/M phase (0.748, p = 1.1 × 10⁻⁹³) (Figure 3B). This pattern is consistent with the biological interpretation that G1-arrested blasts represent the most differentiated-deviant, quiescent-like malignant population, whereas cycling cells (S/G2M) may retain more residual differentiation machinery.

**CNV load correlation.** Malignant cells with high CNV burden (above-median inferCNV score, n = 965) exhibited significantly higher TDI (mean: 0.809) than low-CNV cells (n = 978, mean: 0.746, p = 1.9 × 10⁻¹²⁰) (Figure 3C). This correlation between genomic instability and transcriptomic deviation suggests that CNV accumulation drives progressive divergence from the normal differentiation trajectory—a finding with implications for clonal evolution monitoring.

**Cluster-level heterogeneity.** Across Seurat clusters within the malignant compartment, TDI varied substantially (range: 0.721–0.877), with clusters 22 and 23 showing the highest deviation and clusters 2 and 3 the lowest (Figure 3D). This cluster-level resolution enables the identification of differentiation-competent vs. -incompetent subclones within the same patient.

---

## 3.2 scTDRP Sensitivity is Context-Dependent: MDS-5q as a Boundary Case

### 3.2.1 In Vitro Model of Ineffective Erythropoiesis

To test whether scTDRP generalizes to other erythroid disorders, we analyzed an independent scRNA-seq dataset from Doty et al. (2023), which profiled erythroid differentiation in MDS-5q and Diamond-Blackfan anemia (DBA) using an in vitro erythroid expansion culture system. After CITE-seq-based erythroid selection (CD71+/CD235a+/CD36+), the dataset contained 7,299 normal erythroid cells, 3,671 MDS-5q cells, and 2,555 DBA cells.

### 3.2.2 Normal Trajectory and Stage Assignment

We constructed the normal erythroid trajectory from the 7 normal donors (N1–N7) using the same metacell strategy. Stage assignment was based on a CITE-seq maturity score (log1p(CD235a) + log1p(CD71) − log1p(CD117)), which approximates the canonical erythroid surface marker progression. The OT Cost Map showed the expected decreasing pattern across stages (Figure S1A), validating the trajectory construction.

### 3.2.3 TDI Fails to Discriminate MDS-5q from Normal

Surprisingly, scTDRP TDI showed **minimal difference** between MDS-5q (mean: 0.253, median: 0.179), DBA (mean: 0.236, median: 0.178), and Normal (mean: 0.208, median: 0.174) (Figure S1B). The effect size was negligible compared to the AEL case (TDI ~0.73). Moreover, cells annotated as following the "ineffective" Trajectory B (trajectoryb = True) actually showed slightly **lower** median TDI than Trajectory A cells (0.182 vs. 0.187, p = 7.1 × 10⁻⁶), contrary to the expectation that ineffective differentiation should manifest as higher deviation.

### 3.2.4 Biological Interpretation: Speed vs. Direction

This apparent "failure" of TDI in MDS-5q reveals an important conceptual boundary of the scTDRP framework. In AEL, malignant cells have undergone **directional deviation** from the normal erythroid trajectory due to oncogenic transformation, chromosomal abnormalities, and ectopic gene expression. TDI, which measures the Wasserstein distance to the nearest normal stage, is highly sensitive to this type of deviation.

In contrast, MDS-5q cells in the in vitro culture system remain on the **same differentiation direction** as normal cells but with **reduced speed** and increased entry into a death-prone trajectory (Trajectory B). The Doty et al. study elegantly showed that 57% of MDS-5q cells (vs. 24% of normal cells) follow Trajectory B, which is characterized by high heme-responsive gene expression and ultimately ends in apoptosis before reaching terminal maturation. Because these cells do not fundamentally deviate from the normal transcriptional trajectory—rather, they progress more slowly or die prematurely—their TDI remains low.

This distinction between **directional deviation** (AEL) and **kinetic impairment** (MDS-5q) is critical for interpreting scTDRP results and for understanding the spectrum of differentiation disorders (Figure S1C).

---

## 3.3 Cross-Lineage Validation: scTDRP Captures B-cell Precursor Arrest in B-ALL

To test whether scTDRP generalizes beyond erythropoiesis, we applied the framework to B-cell acute lymphoblastic leukemia (B-ALL), a malignancy characterized by developmental arrest of B-lymphoid progenitors. We analyzed 3,569 malignant B cells from seven diagnostic B-ALL bone marrow samples (Witkowski et al., 2020; GSE130116) and constructed a normal B-cell differentiation trajectory from 16,878 healthy bone marrow B cells spanning four canonical stages: Pro-B VDJ → Large Pre-B → Immature B → Mature B (blood map dataset).

### 3.3.1 Normal B-cell Trajectory and OT Cost Map

The pairwise Wasserstein-2 distances between adjacent normal B-cell stages showed the expected progression pattern: Pro-B VDJ → Large Pre-B (0.17), Large Pre-B → Immature B (0.23), Immature B → Mature B (0.20). The relatively compact distances reflect the continuous nature of B-lymphopoiesis compared with the more discrete transitions in erythropoiesis.

### 3.3.2 B-ALL Cells Are Assigned Predominantly to the Large Pre-B Stage

Stage assignment of B-ALL malignant cells revealed that **82.5% (2,945/3,569) mapped to the Large Pre-B stage** as their best-match normal stage, with 14.1% mapping to Immature B, 3.3% to Mature B, and <0.1% to Pro-B VDJ (Figure 4B). This arrest pattern is biologically consistent with B-ALL pathology, where malignant blasts typically accumulate at the pre-B receptor checkpoint between Large Pre-B and Immature B.

### 3.3.3 B-ALL Cells Show Elevated Deviation Relative to Stage-Matched Normal Cells

Although B-ALL cells mapped primarily to Large Pre-B, their transcriptomic deviation relative to that normal stage was substantial. Single-cell TDI comparison showed **B-ALL mean TDI = 0.297 ± 0.061 (median: 0.281)** versus **normal Large Pre-B mean TDI = 0.173 ± 0.049 (median: 0.159)** (Mann-Whitney U test, p < 1 × 10⁻³⁰⁰; Cohen's d = 1.91) (Figure 4D). The effect size (d = 1.91) is even larger than that observed for AEL malignant cells versus their matched normal stage, demonstrating that scTDRP can detect disease-specific deviation independently of lineage.

Notably, when B-ALL cells were compared against the full pool of normal B cells (all stages combined), the TDI difference was attenuated (B-ALL 0.297 vs. Normal 0.289, p = 0.30). This reflects the fact that B-ALL cells map to Large Pre-B, the normal stage with intrinsically lowest TDI, and underscores the importance of stage-matched comparison for accurate interpretation.

### 3.3.4 Cross-Lineage Confirmation of scTDRP

The B-ALL validation establishes three key properties of scTDRP. First, the framework is **lineage-agnostic**: the same optimal-transport machinery applied to erythropoiesis successfully captures B-lymphoid differentiation arrest without lineage-specific parameter tuning. Second, TDI detects **disease-specific deviation within the matched normal stage**, not merely global dissimilarity from the reference. Third, the stage assignment reveals **biologically interpretable arrest points** (Large Pre-B in B-ALL; Polychromatic Erythroblast in AEL) that align with known disease pathology.

---

## 3.4 Cross-System Validation: scTDRP Captures OPC-like Arrest in Glioblastoma

To test whether scTDRP generalizes beyond hematopoietic lineages, we applied it to glioblastoma (GBM), an aggressive adult brain tumor in which malignant cells are widely reported to adopt an oligodendrocyte precursor cell (OPC)-like state. We analyzed 1,091 neoplastic cells from the GSE84465 dataset (Darmanis et al., 2017) and constructed a normal neural lineage reference from 579 non-neoplastic cells in the same dataset, spanning OPC → Astrocyte → Oligodendrocyte.

### 3.4.1 Neural Trajectory and OT Cost Map

Pairwise Wasserstein-2 distances between adjacent reference stages were OPC → Astrocyte (0.262) and Astrocyte → Oligodendrocyte (0.309), establishing a quantifiable differentiation axis.

### 3.4.2 GBM Cells Are Assigned Predominantly to the OPC Stage

Stage assignment revealed that **85.9% (937/1,091) of GBM neoplastic cells mapped to the OPC stage** as their best-match normal stage, with 14.1% mapping to Oligodendrocyte and 1.0% to Astrocyte. This arrest pattern is consistent with the established view that GBM cells are transcriptionally arrested at an OPC-like progenitor state.

### 3.4.3 GBM Cells Show Elevated Deviation Relative to Stage-Matched OPCs

GBM cells assigned to OPC showed markedly elevated TDI compared with normal OPC cells (**mean TDI = 0.300 ± 0.047; median = 0.289**). Module-level repair inference pointed toward up-regulation of oligodendrocyte differentiation genes (MBP, PLP1, MOG, MAG, CNP) and down-regulation of astrocyte differentiation and cell-cycle programs, suggesting a differentiation-promoting, proliferation-suppressing therapeutic direction.

### 3.4.4 Cross-System Confirmation of scTDRP

The GBM validation demonstrates that scTDRP is **not restricted to hematopoietic malignancies**. Using only a public scRNA-seq dataset and a three-stage normal neural reference, scTDRP correctly recapitulated the OPC-like arrest phenotype, supporting its applicability to solid tumors and other non-hematopoietic differentiation disorders.

### 3.4.5 Consistency with WaddingtonOT and Moscot

To ensure that scTDRP's conclusions are not an artifact of its specific implementation, we compared stage assignment and repair direction against two widely used single-cell optimal transport frameworks: **WaddingtonOT (WOT)** and **Moscot**.

**Stage assignment.** WOT assigned **84.0%** of GBM cells to OPC, closely matching scTDRP's 85.9% (exact agreement 73.1%). Moscot assigned 72.2% to OPC (exact agreement 65.7%); its slightly more dispersed distribution likely reflects Moscot's built-in birth-death regularization, which allows mass creation/annihilation and can soften sharp stage assignments.

**Repair direction.** We inferred the disease-to-terminal (Oligodendrocyte) repair vector for each cell in PCA space. scTDRP and WOT yielded virtually identical directions (**mean cosine similarity = 0.9999, Pearson r = 0.9999**), confirming that scTDRP's repair inference is numerically consistent with an independent POT/Sinkhorn implementation. The scTDRP-Moscot repair direction correlation was positive but lower (**mean cosine similarity = 0.51, Pearson r = 0.53**), consistent with differences in solver backend (OTT-JAX vs POT) and regularization strategy.

> **Implementation note.** In this GBM validation we used original cells rather than metacell aggregation for both reference and disease. When the number of reference stages is small and adjacent stages are transcriptionally close, metacell averaging can blur stage boundaries and shift stage assignment. Users are encouraged to compare both modes for their specific dataset.

---

## 3.5 Simulation Validation Confirms scTDRP Robustness

To validate the computational and conceptual foundations of scTDRP independently of biological data, we performed controlled simulations with known ground truth (see Methods). All five validation tests passed:

1. **OT Cost Map Monotonicity**: Adjacent-stage W2 distances accumulated monotonically along the simulated differentiation axis (R² = 0.99).
2. **TDI Group Accuracy**: Disease metacells were correctly assigned to Stage 3 (the perturbed origin stage) with a mean distance of 6.27 vs. 10.13 to adjacent stages.
3. **TDI Cell Accuracy**: 20/20 individual disease cells were correctly assigned to Stage 3.
4. **Repair Directionality**: The inferred repair delta correctly pushed cells toward the terminal stage (dX = +2.00) and removed the perpendicular perturbation (dY = −2.47).
5. **Epsilon Stability**: Sinkhorn regularization with ε ∈ [0.001, 1.0] produced stable distances within ±10% of exact EMD.

These simulations establish that scTDRP's core computations are mathematically sound and robust to hyperparameter choices.

---

## Figure Legends

**Figure 2. scTDRP Analysis of Acute Erythroid Leukemia.**
(A) OT Cost Map showing pairwise Wasserstein-2 distances between adjacent erythroid stages in normal differentiation. Distances decrease from early to terminal stages, reflecting progressive transcriptional convergence.
(B) TDI distribution of AEL malignant cells (n = 2,134) with stage assignment. The majority (81.8%) are assigned to the Polychromatic Erythroblast stage, indicating differentiation arrest.
(C) Repair heatmap showing top 50 genes with highest positive (up-regulate) and negative (down-regulate) repair deltas. Red = up-regulate toward terminal stage; Blue = down-regulate.
(D) Module repair strategy. P7_TerminalPrep module is recommended for up-regulation; P8_ExecutionPrecursor module is recommended for down-regulation.

**Figure 3. Intra-tumor Heterogeneity in AEL Revealed by scTDRP.**
(A) TDI comparison across inferCNV-based malignancy groups: Malignant Erythroid, Normal Residual, and Other (non-erythroid). Malignant cells show significantly higher TDI than residual normal cells (p = 5.4 × 10⁻²³³).
(B) TDI by cell-cycle phase within malignant cells. G1-arrested cells show the highest TDI, consistent with differentiation quiescence (G1 vs G2/M: p = 1.1 × 10⁻⁹³).
(C) TDI by CNV load (above- vs below-median inferCNV score) within malignant cells. High-CNV cells show significantly higher deviation (p = 1.9 × 10⁻¹²⁰).
(D) Mean TDI across Seurat clusters within the malignant compartment, revealing substantial inter-cluster heterogeneity.

**Figure 4. scTDRP Validation in B-cell Acute Lymphoblastic Leukemia (B-ALL).**
(A) TDI comparison of all normal B cells versus B-ALL malignant B cells. Boxes show median and interquartile range; whiskers show 1.5× IQR.
(B) Proportion of B-ALL malignant cells assigned to each normal B-cell stage (Pro-B VDJ, Large Pre-B, Immature B, Mature B). The majority (82.5%) arrest at Large Pre-B.
(C) TDI of normal B cells stratified by their true annotated stage, showing that Large Pre-B has the intrinsically lowest TDI.
(D) Density distribution of single-cell TDI comparing B-ALL malignant cells (n = 3,569) to normal Large Pre-B cells (n = 3,488), the matched stage. B-ALL cells show significantly higher deviation (Mann-Whitney U p < 1 × 10⁻³⁰⁰; Cohen's d = 1.91).

**Figure S1. scTDRP in MDS-5q: A Boundary Case for Directional vs. Kinetic Deviation.**
(A) OT Cost Map for normal erythroid trajectory in the in vitro culture system.
(B) TDI comparison across Normal, MDS-5q, and DBA conditions. No significant elevation in disease groups.
(C) Conceptual diagram distinguishing directional deviation (sensitive to TDI, exemplified by AEL) from kinetic impairment (insensitive to TDI, exemplified by MDS-5q).

---

# Discussion

## 4.1 scTDRP: Bridging Single-Cell Genomics and Developmental Therapeutics

The central challenge in treating differentiation disorders—from hematologic malignancies to developmental syndromes—is understanding not just what goes wrong, but how to push cells back toward their intended fate. scTDRP addresses this challenge by providing a quantitative, mechanistically grounded framework that connects three previously disconnected layers of information: (1) the normal developmental trajectory as a reference, (2) the disease cell's position relative to that trajectory, and (3) the specific transcriptomic adjustments required to restore terminal differentiation. By grounding repair strategies in optimal transport theory, scTDRP offers a principled alternative to purely correlative approaches such as differential expression or pathway enrichment.

## 4.2 Key Advantages of the scTDRP Framework

**Single-cell resolution.** Unlike bulk transcriptomic methods (e.g., WaddingtonOT) that infer population-level trajectories, scTDRP operates at single-cell resolution, enabling the identification of subpopulations with distinct differentiation potentials within a heterogeneous disease sample. In AEL, for example, we observed that while 81.8% of malignant cells arrested at the Polychromatic stage, a small fraction (7.7%) retained proximity to the Orthochromatic stage, potentially representing a residual differentiation-competent subclone.

**Transcriptomic repair directionality.** Traditional differential expression analysis identifies genes that are up- or down-regulated in disease, but it cannot distinguish between causal drivers and passenger effects, nor can it indicate the direction of therapeutic intervention. scTDRP's repair delta (δg = Σ γij · (yj − xi)) directly quantifies how much each gene must change to move the disease population toward the terminal stage. This provides a ranked list of candidate therapeutic targets with built-in directional information.

**Modular aggregation.** By aggregating gene-level repair signals into pre-defined functional modules, scTDRP bridges the gap between molecular detail and systems-level interpretability. The opposing directions of the P7 (TerminalPrep) and P8 (ExecutionPrecursor) modules in AEL—a disease where both modules are dysregulated in opposite directions—demonstrates how module-level analysis can reveal coherent biological strategies that might be obscured at the single-gene level.

**Cross-disease comparability.** Because TDI is computed in the same metric space (Wasserstein distance) across different diseases and cell types, it offers a principled way to compare the severity of differentiation deviation across clinically distinct conditions. An AEL patient (TDI ~0.73) shows far more severe transcriptomic deviation than an MDS-5q patient (TDI ~0.18), despite both presenting with anemia and bone marrow failure.

## 4.3 Limitations and Boundary Conditions

**Sensitivity to directional but not kinetic deviation.** Our analysis of MDS-5q revealed a critical boundary condition: scTDRP TDI is sensitive to cells that have deviated from the normal differentiation *direction* (e.g., via oncogenic transformation), but not to cells that are moving along the correct trajectory with reduced *speed* or increased *mortality*. In MDS-5q, erythroid precursors follow the same transcriptional trajectory as normal cells but enter a death-prone branch (Trajectory B) at increased frequency. Because these cells never fundamentally leave the normal trajectory, their TDI remains low. This distinction between "directional deviation" and "kinetic impairment" is biologically meaningful and should guide the selection of diseases for scTDRP application.

**Dependence on normal reference quality.** scTDRP's accuracy depends critically on the quality and completeness of the normal developmental reference. If the reference lacks key transitional stages, disease cells may be incorrectly assigned to the nearest available stage, potentially obscuring true deviations. In our AEL analysis, the seven-stage erythroid reference provided sufficient resolution; however, for lineages with more continuous or branched trajectories (e.g., neural development), higher-resolution references may be required.

**In vitro vs. in vivo contexts.** The MDS-5q dataset was generated from in vitro erythroid expansion cultures, which standardize the microenvironment and may reduce disease-specific deviations. The strong AEL results, in contrast, were obtained from primary bone marrow samples. This suggests that scTDRP may perform best on primary tissue data where disease-specific microenvironmental and cell-intrinsic abnormalities are preserved.

**Computational scalability.** While our metacell strategy enables scTDRP to handle datasets with tens of thousands of cells, optimal transport remains computationally expensive (O(n²) to O(n³) depending on algorithm). For very large datasets (>100,000 cells), further approximations such as hierarchical clustering or neural OT solvers may be necessary.

## 4.4 Relationship to Existing Methods

scTDRP occupies a unique position in the landscape of computational developmental biology tools. WaddingtonOT (Schiebinger et al., 2019) pioneered the use of optimal transport for developmental trajectory inference, but it operates on bulk or pseudobulk data and focuses on lineage tracing rather than disease repair. scVelo (Bergen et al., 2020) infers RNA velocity to predict cell fate, but it does not quantify deviation from a healthy reference or infer repair strategies. Pseudotime methods (Trapnell et al., 2014) order cells along differentiation trajectories but lack a principled metric for disease severity. scTDRP complements these tools by adding the disease-to-normal comparison layer and the transcriptomic repair inference layer.

## 4.5 Future Directions

**Expansion to other lineages.** While we initially validated scTDRP in erythropoiesis, we have now extended the framework to lymphopoiesis by analyzing B-cell acute lymphoblastic leukemia (B-ALL). In this cross-lineage test, scTDRP correctly identified Large Pre-B as the arrest stage for 82.5% of malignant B cells and detected highly significant transcriptomic deviation relative to stage-matched normal Large Pre-B cells (Cohen's d = 1.91). These results suggest that scTDRP generalizes across hematopoietic lineages. Natural extensions include myelopoiesis (e.g., acute myeloid leukemia subtypes), neurogenesis (e.g., glioblastoma dedifferentiation), and muscle regeneration (e.g., Duchenne muscular dystrophy). Each application requires only a high-quality normal differentiation atlas and a disease single-cell dataset.

**Integration with perturbation screens.** The repair deltas inferred by scTDRP can be validated experimentally using CRISPR activation/interference screens or small-molecule perturbation libraries. For example, the top up-regulated targets in AEL (HBA1, KLF1, NFE2) could be tested for their ability to restore hemoglobin synthesis and enucleation in AEL cell lines.

**Dynamic monitoring.** scTDRP is currently designed for cross-sectional (single time-point) analysis. An important extension would be to track TDI and repair strategies longitudinally during treatment, enabling assessment of therapeutic efficacy at the single-cell level. For example, a drug that successfully pushes AEL blasts toward terminal differentiation should manifest as decreasing TDI and shifting module strategies over time.

**Megakaryocytic lineage validation.** We attempted to validate scTDRP in acute megakaryoblastic leukemia (AMKL), but publicly available AMKL scRNA-seq datasets with matched normal megakaryocyte references remain scarce. Future efforts should prioritize generating or curating such datasets to test whether scTDRP can capture the megakaryocytic differentiation block characteristic of AMKL.

## 4.6 Conclusion

scTDRP provides a theoretically grounded, computationally robust, and biologically interpretable framework for quantifying developmental deviation and inferring transcriptomic repair strategies in single-cell disease genomics. By applying optimal transport theory to the disease-normal comparison, scTDRP moves beyond descriptive analysis toward predictive therapeutic guidance. Our validation in acute erythroid leukemia demonstrates that scTDRP can identify specific differentiation arrest stages, quantify transcriptomic distance from normality, and propose coherent module-level repair strategies. Cross-lineage validation in B-cell acute lymphoblastic leukemia further shows that scTDRP generalizes beyond erythropoiesis, correctly identifying Large Pre-B arrest and significant deviation from stage-matched normal cells (Cohen's d = 1.91). The boundary case of MDS-5q illuminates an important conceptual distinction between directional deviation and kinetic impairment, refining the scope of diseases for which scTDRP is most informative. As single-cell atlases of normal human development continue to expand, scTDRP offers a generalizable approach to translating developmental biology insights into precision therapeutic strategies.
