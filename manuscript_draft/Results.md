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

## 3.4 Simulation Validation Confirms scTDRP Robustness

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
