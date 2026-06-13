#!/usr/bin/env python3
"""
Direct comparison: scTDRP vs Moscot on GBM cross-system validation.

We use Moscot's TemporalProblem to compute:
  1. Disease-to-normal-stage transition probabilities (analogous to scTDRP stage assignment).
  2. Disease-to-terminal-stage coupling and inferred repair direction.

Results are compared with scTDRP's TDI / best-match stage and repair delta.
"""

import os
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
np.random.seed(42)

DATA_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/data/gbm"
OUT_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/results_OT_comparison"
os.makedirs(OUT_DIR, exist_ok=True)

STAGE_ORDER = ["OPC", "Astrocyte", "Oligodendrocyte"]
TERMINAL_STAGE = "Oligodendrocyte"

import sys
sys.path.insert(0, "/data1/yja/zhongzhuan/5.external/scTDRP/src")
from scTDRP import TDRPAnalyzer
from scTDRP.utils import compute_expression_distribution


def load_and_preprocess():
    """Reproduce GBM preprocessing from run_gbm_validation.py."""
    df = pd.read_csv(
        os.path.join(DATA_DIR, "GSE84465_GBM_All_data.csv.gz"),
        index_col=0, sep=r"\s+"
    )
    meta = pd.read_csv(os.path.join(DATA_DIR, "GSE84465_metadata.csv"), index_col=0)
    meta = meta.rename(columns={
        "Unnamed: 2": "diagnosis",
        "_1": "plate_id", "_2": "well", "_3": "tissue",
        "_4": "patient_id", "_5": "tsne_cluster",
        "_6": "cell_type", "_7": "neoplastic", "_8": "selection",
    })
    meta["sample_id"] = meta["plate_id"].astype(str) + "." + meta["well"].astype(str)
    common = list(set(df.columns) & set(meta["sample_id"]))
    df = df[common]
    meta = meta.set_index("sample_id").loc[common]

    adata = sc.AnnData(X=df.T.values)
    adata.obs_names = df.columns
    adata.var_names = df.index
    adata.var_names_make_unique()
    for col in meta.columns:
        adata.obs[col] = meta[col].values

    normal_types = {"OPC", "Astocyte", "Oligodendrocyte"}
    keep = (adata.obs["neoplastic"] == "Neoplastic") | adata.obs["cell_type"].isin(normal_types)
    adata = adata[keep].copy()
    sc.pp.filter_genes(adata, min_cells=3)

    adata.obs["source"] = adata.obs["neoplastic"].map(
        {"Neoplastic": "disease", "Regular": "normal"}
    ).fillna("normal")
    stage_rename = {"Astocyte": "Astrocyte"}
    adata.obs["stage"] = adata.obs["cell_type"].map(
        lambda x: stage_rename.get(x, x) if pd.notna(x) else x
    ).where(adata.obs["source"] == "normal", "Disease")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

    return adata


def run_scTDRP(adata):
    """Run scTDRP on the same data and return cell-level results."""
    normal = adata[adata.obs["source"] == "normal"].copy()
    disease = adata[adata.obs["source"] == "disease"].copy()

    analyzer = TDRPAnalyzer(
        normal_adata=normal,
        stage_key="stage",
        terminal_stage=TERMINAL_STAGE,
        stage_order=STAGE_ORDER,
        use_rep="X_pca",
        n_top_genes=2000,
    )
    analyzer.prepare_data(flavor="seurat_v3")
    analyzer.build_ot_cost_map(metric="sqeuclidean", reg=0.01)
    tdi_df = analyzer.compute_tdi(disease, metric="sqeuclidean", reg=0.01)

    disease.obs["scTDRP_TDI"] = tdi_df["TDI"].values
    disease.obs["scTDRP_BestStage"] = tdi_df["Best_Match_Stage"].values

    # Repair direction in PCA space
    terminal = normal[normal.obs["stage"] == TERMINAL_STAGE].copy()
    repair = analyzer.infer_repair_pathway(
        disease, metric="sqeuclidean", reg=0.01, top_n=100
    )
    disease.obs["scTDRP_RepairDist"] = repair["wasserstein_distance"]
    return disease, analyzer, repair


def run_moscot_stage_assignment(adata):
    """Use Moscot TemporalProblem to get disease-to-stage transition probabilities."""
    from moscot.problems.time import TemporalProblem

    # Prepare AnnData with time = 0 (disease) and 1 (normal)
    ad = adata.copy()
    ad.obs["day"] = ad.obs["source"].map({"disease": 0, "normal": 1}).astype(int)
    ad.obs["stage_day"] = pd.Categorical(ad.obs["stage"].astype(str), categories=STAGE_ORDER + ["Disease"])

    tp = TemporalProblem(ad)
    tp = tp.prepare(time_key="day", joint_attr="X_pca")
    tp = tp.solve(epsilon=0.01)

    # Cell transition from disease (day=0) to normal stage (day=1)
    ct = tp.cell_transition(
        source=0,
        target=1,
        source_groups=None,
        target_groups="stage_day",
        forward=True,
        aggregation_mode="cell",
    )

    # ct is a DataFrame: rows = disease cells, columns = stages (including 'Disease' if any)
    stage_cols = [c for c in ct.columns if c in STAGE_ORDER]
    ct = ct[stage_cols]
    moscot_best = ct.idxmax(axis=1)
    moscot_entropy = -(ct * np.log(ct + 1e-12)).sum(axis=1)

    return ct, moscot_best, moscot_entropy


def run_moscot_repair(adata):
    """Use Moscot TemporalProblem to infer disease-to-terminal repair direction."""
    from moscot.problems.time import TemporalProblem

    ad = adata.copy()
    ad.obs["day"] = np.where(ad.obs["source"] == "disease", 0, 1).astype(int)

    # Keep only disease and terminal stage cells
    mask = (ad.obs["source"] == "disease") | (
        (ad.obs["source"] == "normal") & (ad.obs["stage"] == TERMINAL_STAGE)
    )
    ad2 = ad[mask].copy()

    tp = TemporalProblem(ad2)
    tp = tp.prepare(time_key="day", joint_attr="X_pca")
    tp = tp.solve(epsilon=0.01)

    # Pull terminal PCA back onto disease cells via the OT coupling
    # Result: each disease cell gets a weighted average of terminal-cell PC coordinates
    terminal_mask = ad2.obs["source"] == "normal"
    disease_mask = ad2.obs["source"] == "disease"
    terminal_pca = ad2.obsm["X_pca"][terminal_mask]

    prob = tp.problems[(0, 1)]
    moscot_target_pca = prob.solution.pull(terminal_pca)

    disease_pca = ad2.obsm["X_pca"][disease_mask]
    moscot_delta_pca = moscot_target_pca - disease_pca
    return moscot_delta_pca


def run_wot_stage_assignment(adata):
    """Use WaddingtonOT to compute disease-to-normal-stage mass distribution."""
    import wot

    ad = adata.copy()
    # Build WOT-compatible AnnData with PCA in .X
    X_pca = ad.obsm["X_pca"]
    ad_wot = sc.AnnData(X=X_pca)
    ad_wot.obs_names = ad.obs_names
    ad_wot.obs["day"] = np.where(ad.obs["source"] == "disease", 0.0, 1.0)
    ad_wot.obs["stage_wot"] = pd.Categorical(
        ad.obs["stage"].astype(str), categories=STAGE_ORDER + ["Disease"]
    )

    ot_model = wot.ot.OTModel(ad_wot, day_field="day", local_pca=0, epsilon=0.05, lambda1=1, lambda2=50)
    tmap_ad = ot_model.compute_transport_map(0, 1)
    plan = np.asarray(tmap_ad.X)
    plan_norm = plan / (plan.sum(axis=1, keepdims=True) + 1e-12)

    # Aggregate mass by target stage
    target_stages = ad_wot[tmap_ad.var_names].obs["stage_wot"].values
    stage_mass = pd.DataFrame(0.0, index=tmap_ad.obs_names, columns=STAGE_ORDER)
    for stage in STAGE_ORDER:
        cols = np.where(target_stages == stage)[0]
        if len(cols) > 0:
            stage_mass[stage] = plan_norm[:, cols].sum(axis=1)

    wot_best = stage_mass.idxmax(axis=1)
    return stage_mass, wot_best


def run_wot_repair(adata):
    """Use WaddingtonOT to infer disease-to-terminal repair direction."""
    import wot

    ad = adata.copy()
    mask = (ad.obs["source"] == "disease") | (
        (ad.obs["source"] == "normal") & (ad.obs["stage"] == TERMINAL_STAGE)
    )
    ad2 = ad[mask].copy()

    X_pca = ad2.obsm["X_pca"]
    ad_wot = sc.AnnData(X=X_pca)
    ad_wot.obs_names = ad2.obs_names
    ad_wot.obs["day"] = np.where(ad2.obs["source"] == "disease", 0.0, 1.0)

    ot_model = wot.ot.OTModel(ad_wot, day_field="day", local_pca=0, epsilon=0.05, lambda1=1, lambda2=50)
    tmap_ad = ot_model.compute_transport_map(0, 1)
    plan = np.asarray(tmap_ad.X)
    plan_norm = plan / (plan.sum(axis=1, keepdims=True) + 1e-12)

    terminal_mask = ad2.obs["source"] == "normal"
    disease_mask = ad2.obs["source"] == "disease"
    terminal_pca = ad2.obsm["X_pca"][terminal_mask]

    wot_target_pca = plan_norm @ terminal_pca
    disease_pca = ad2.obsm["X_pca"][disease_mask]
    wot_delta_pca = wot_target_pca - disease_pca
    return wot_delta_pca


def evaluate_stage_consistency(disease_obs, moscot_best, moscot_probs):
    """Compare scTDRP and Moscot stage assignments."""
    print("\n" + "=" * 70)
    print("Stage Assignment Comparison")
    print("=" * 70)

    print("\nscTDRP best-stage distribution:")
    print(disease_obs["scTDRP_BestStage"].value_counts(normalize=True))

    print("\nMoscot best-stage distribution:")
    print(moscot_best.value_counts(normalize=True))

    # Agreement
    agreement = (disease_obs["scTDRP_BestStage"].values == moscot_best.values).mean()
    print(f"\nExact agreement: {agreement*100:.1f}%")

    # Confusion matrix
    confusion = pd.crosstab(
        disease_obs["scTDRP_BestStage"], moscot_best,
        rownames=["scTDRP"], colnames=["Moscot"],
        normalize="index"
    )
    print("\nRow-normalized confusion (scTDRP -> Moscot):")
    print(confusion.round(3))

    # Mean Moscot probability for scTDRP-assigned stage
    probs_assigned = []
    for i, stage in enumerate(disease_obs["scTDRP_BestStage"]):
        if stage in moscot_probs.columns:
            probs_assigned.append(moscot_probs.iloc[i][stage])
        else:
            probs_assigned.append(np.nan)
    print(f"\nMean Moscot probability for scTDRP-assigned stage: {np.nanmean(probs_assigned):.3f}")

    return confusion


def plot_results(disease_obs, moscot_probs, confusion):
    """Generate comparison figures."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Stage distribution comparison
    stage_comp = pd.DataFrame({
        "scTDRP": disease_obs["scTDRP_BestStage"].value_counts(normalize=True),
        "Moscot": moscot_probs.idxmax(axis=1).value_counts(normalize=True),
    }).reindex(STAGE_ORDER).fillna(0)
    stage_comp.plot(kind="bar", ax=axes[0], color=["steelblue", "coral"])
    axes[0].set_ylabel("Fraction of disease cells")
    axes[0].set_title("Stage assignment: scTDRP vs Moscot")
    axes[0].tick_params(axis="x", rotation=45)
    axes[0].legend()

    # Confusion heatmap
    sns.heatmap(confusion.reindex(index=STAGE_ORDER, columns=STAGE_ORDER, fill_value=0),
                annot=True, fmt=".2f", cmap="Blues", ax=axes[1], vmin=0, vmax=1)
    axes[1].set_title("Stage assignment confusion")
    axes[1].set_xlabel("Moscot")
    axes[1].set_ylabel("scTDRP")

    plt.tight_layout()
    out = os.path.join(OUT_DIR, "scTDRP_vs_moscot_stage.pdf")
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"\nSaved: {out}")


def main():
    print("=" * 70)
    print("scTDRP vs WaddingtonOT / Moscot comparison on GBM")
    print("=" * 70)

    adata = load_and_preprocess()
    disease, analyzer, repair = run_scTDRP(adata)

    print("\nRunning Moscot stage assignment...")
    moscot_probs, moscot_best, moscot_entropy = run_moscot_stage_assignment(adata)

    print("\nRunning Moscot repair direction...")
    moscot_repair = run_moscot_repair(adata)

    print("\nRunning WaddingtonOT stage assignment...")
    wot_probs, wot_best = run_wot_stage_assignment(adata)

    print("\nRunning WaddingtonOT repair direction...")
    wot_repair = run_wot_repair(adata)

    disease.obs["Moscot_BestStage"] = moscot_best.values
    disease.obs["Moscot_Entropy"] = moscot_entropy.values
    disease.obs["WOT_BestStage"] = wot_best.values

    confusion = evaluate_stage_consistency(disease.obs, moscot_best, moscot_probs)
    plot_results(disease.obs, moscot_probs, confusion)

    # WOT stage consistency
    print("\n" + "=" * 70)
    print("WOT Stage Assignment")
    print("=" * 70)
    print("\nWOT best-stage distribution:")
    print(wot_best.value_counts(normalize=True))
    wot_agree_scTDRP = (disease.obs["scTDRP_BestStage"].values == wot_best.values).mean()
    wot_agree_moscot = (moscot_best.values == wot_best.values).mean()
    print(f"\nExact agreement scTDRP-WOT: {wot_agree_scTDRP*100:.1f}%")
    print(f"Exact agreement Moscot-WOT: {wot_agree_moscot*100:.1f}%")

    # Repair direction comparison in PCA space
    print("\n" + "=" * 70)
    print("Repair Direction Comparison (PCA space)")
    print("=" * 70)

    plan = repair["transport_plan"]
    plan_norm = plan / (plan.sum(axis=1, keepdims=True) + 1e-12)
    terminal = adata[
        (adata.obs["source"] == "normal") & (adata.obs["stage"] == TERMINAL_STAGE)
    ]
    scTDRP_target_pca = plan_norm @ terminal.obsm["X_pca"]
    scTDRP_delta_pca = scTDRP_target_pca - disease.obsm["X_pca"]

    # Moscot and WOT deltas already computed in PCA space
    moscot_delta_pca = moscot_repair
    wot_delta_pca = wot_repair

    def repair_cosine(a, b):
        cos_sims = []
        for i in range(a.shape[0]):
            ai = a[i]
            bi = b[i]
            norm = np.linalg.norm(ai) * np.linalg.norm(bi)
            cos = np.dot(ai, bi) / (norm + 1e-12)
            cos_sims.append(cos)
        return np.array(cos_sims)

    cos_moscot = repair_cosine(scTDRP_delta_pca, moscot_delta_pca)
    cos_wot = repair_cosine(scTDRP_delta_pca, wot_delta_pca)

    print("\n  scTDRP vs Moscot:")
    print(f"    Mean cosine similarity: {cos_moscot.mean():.4f}")
    print(f"    Median cosine similarity: {np.median(cos_moscot):.4f}")
    print(f"    Fraction positive (>0): {(cos_moscot > 0).mean()*100:.1f}%")
    pearson_m = np.corrcoef(scTDRP_delta_pca.flatten(), moscot_delta_pca.flatten())[0, 1]
    print(f"    Pearson correlation (all PC deltas): {pearson_m:.4f}")

    print("\n  scTDRP vs WaddingtonOT:")
    print(f"    Mean cosine similarity: {cos_wot.mean():.4f}")
    print(f"    Median cosine similarity: {np.median(cos_wot):.4f}")
    print(f"    Fraction positive (>0): {(cos_wot > 0).mean()*100:.1f}%")
    pearson_w = np.corrcoef(scTDRP_delta_pca.flatten(), wot_delta_pca.flatten())[0, 1]
    print(f"    Pearson correlation (all PC deltas): {pearson_w:.4f}")

    # Save outputs
    disease.obs.to_csv(os.path.join(OUT_DIR, "gbm_scTDRP_moscot_obs.csv"))
    moscot_probs.to_csv(os.path.join(OUT_DIR, "gbm_moscot_stage_probs.csv"))
    confusion.to_csv(os.path.join(OUT_DIR, "gbm_stage_confusion.csv"))
    wot_probs.to_csv(os.path.join(OUT_DIR, "gbm_wot_stage_probs.csv"))
    np.save(os.path.join(OUT_DIR, "gbm_scTDRP_delta_pca.npy"), scTDRP_delta_pca)
    np.save(os.path.join(OUT_DIR, "gbm_moscot_delta_pca.npy"), moscot_delta_pca)
    np.save(os.path.join(OUT_DIR, "gbm_wot_delta_pca.npy"), wot_delta_pca)

    # Plot repair direction cosine similarity
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(cos_moscot, bins=40, color="steelblue", edgecolor="black")
    axes[0].axvline(cos_moscot.mean(), color="red", linestyle="--", label=f"mean={cos_moscot.mean():.3f}")
    axes[0].set_xlabel("Cosine similarity")
    axes[0].set_ylabel("Number of cells")
    axes[0].set_title("scTDRP vs Moscot repair direction")
    axes[0].legend()

    axes[1].hist(cos_wot, bins=40, color="coral", edgecolor="black")
    axes[1].axvline(cos_wot.mean(), color="red", linestyle="--", label=f"mean={cos_wot.mean():.3f}")
    axes[1].set_xlabel("Cosine similarity")
    axes[1].set_ylabel("Number of cells")
    axes[1].set_title("scTDRP vs WaddingtonOT repair direction")
    axes[1].legend()

    plt.tight_layout()
    out_repair = os.path.join(OUT_DIR, "scTDRP_vs_OT_repair_cosine.pdf")
    plt.savefig(out_repair, dpi=300)
    plt.close()
    print(f"\nSaved: {out_repair}")

    print("\n" + "=" * 70)
    print("Comparison complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
