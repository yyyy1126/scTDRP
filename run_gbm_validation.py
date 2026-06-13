#!/usr/bin/env python3
"""
GBM cross-system validation for scTDRP.

Uses GSE84465 (Darmanis et al., 2017):
  - Regular (non-neoplastic) neural cells as the normal reference:
      OPC -> Astrocyte -> Oligodendrocyte
  - Neoplastic cells as the disease sample.

Expected biological findings:
  - GBM cells are transcriptionally arrested at an OPC-like progenitor state.
  - scTDRP should assign most neoplastic cells to the OPC stage.
  - Repair direction should point toward differentiated astrocyte/oligodendrocyte programs.
"""

import os
import json
import warnings
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
np.random.seed(42)

DATA_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/data/gbm"
OUT_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/results_gbm"
FIG_DIR = os.path.join(OUT_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

STAGE_ORDER = ["OPC", "Astrocyte", "Oligodendrocyte"]
TERMINAL_STAGE = "Oligodendrocyte"
EXPECTED_ARREST_STAGE = "OPC"

# Optional imports
import sys
try:
    from scTDRP import TDRPAnalyzer, build_stage_metacells, build_metacells
except Exception:
    sys.path.insert(0, "/data1/yja/zhongzhuan/5.external/scTDRP/src")
    from scTDRP import TDRPAnalyzer, build_stage_metacells, build_metacells


def load_gbm_data():
    """Load GSE84465 expression matrix and metadata, build AnnData."""
    print("=" * 70)
    print("[1] Loading GSE84465 expression matrix")
    print("=" * 70)

    expr_path = os.path.join(DATA_DIR, "GSE84465_GBM_All_data.csv.gz")
    meta_path = os.path.join(DATA_DIR, "GSE84465_metadata.csv")

    # The matrix is genes x cells, whitespace-delimited
    df = pd.read_csv(expr_path, index_col=0, sep=r"\s+")
    print(f"  Expression matrix: {df.shape}")

    meta = pd.read_csv(meta_path, index_col=0)
    # Column mapping inferred from series matrix order:
    # Unnamed: 2=diagnosis, _1=plate, _2=well, _3=tissue, _4=patient, _5=tsne cluster,
    # _6=cell type, _7=neoplastic, _8=selection
    meta = meta.rename(columns={
        "Unnamed: 2": "diagnosis",
        "_1": "plate_id",
        "_2": "well",
        "_3": "tissue",
        "_4": "patient_id",
        "_5": "tsne_cluster",
        "_6": "cell_type",
        "_7": "neoplastic",
        "_8": "selection",
    })

    # Build sample IDs from plate_id + well to match matrix columns
    meta["sample_id"] = meta["plate_id"].astype(str) + "." + meta["well"].astype(str)

    # Subset to common samples
    common_samples = list(set(df.columns) & set(meta["sample_id"]))
    print(f"  Common samples between matrix and metadata: {len(common_samples)}")
    df = df[common_samples]
    meta = meta.set_index("sample_id").loc[common_samples]

    if df.shape[1] != meta.shape[0]:
        raise ValueError(f"Matrix columns ({df.shape[1]}) != metadata rows ({meta.shape[0]})")

    # Build AnnData (cells x genes)
    adata = sc.AnnData(X=df.T.values)
    adata.obs_names = df.columns
    adata.var_names = df.index
    adata.var_names_make_unique()

    # Attach metadata
    meta = meta.loc[adata.obs_names]
    for col in meta.columns:
        adata.obs[col] = meta[col].values

    print(f"  AnnData: {adata.shape}")
    print("\n  Cell type distribution:")
    print(adata.obs["cell_type"].value_counts())
    print("\n  Neoplastic distribution:")
    print(adata.obs["neoplastic"].value_counts())

    return adata


def preprocess_gbm(adata):
    """Filter cells/genes, normalize, log1p, HVG, scale, PCA."""
    print("\n" + "=" * 70)
    print("[2] Preprocessing")
    print("=" * 70)

    # Keep only neoplastic cells and relevant normal neural cell types
    # Note: Darmanis metadata spells Astrocyte as "Astocyte"
    normal_types = {"OPC", "Astocyte", "Oligodendrocyte"}
    keep_mask = (
        (adata.obs["neoplastic"] == "Neoplastic") |
        adata.obs["cell_type"].isin(normal_types)
    )
    adata = adata[keep_mask].copy()
    print(f"  After cell-type filter: {adata.shape}")

    # Gene QC
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  After gene filter: {adata.shape}")

    # Label source and stage
    adata.obs["source"] = adata.obs["neoplastic"].map(
        {"Neoplastic": "disease", "Regular": "normal"}
    ).fillna("normal")
    stage_rename = {"Astocyte": "Astrocyte"}
    adata.obs["stage"] = adata.obs["cell_type"].map(
        lambda x: stage_rename.get(x, x) if pd.notna(x) else x
    ).where(adata.obs["source"] == "normal", np.nan)

    # Normalize and log1p (data are TPM-like, but standard log-normalization is fine)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG on all retained cells
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")

    # Scale and PCA
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

    print(f"  HVGs: {adata.var['highly_variable'].sum()}")
    return adata


def build_reference_and_disease(adata):
    """Split into normal reference and disease cells."""
    adata_normal = adata[adata.obs["source"] == "normal"].copy()
    adata_disease = adata[adata.obs["source"] == "disease"].copy()

    # Restrict normal to selected stages
    adata_normal = adata_normal[adata_normal.obs["stage"].isin(STAGE_ORDER)].copy()
    print(f"\n  Normal reference: {adata_normal.shape}")
    print(adata_normal.obs["stage"].value_counts())
    print(f"  Disease cells: {adata_disease.shape}")
    return adata_normal, adata_disease


def run_scTDRP_gbm(adata_normal, adata_disease):
    """Run scTDRP and return results."""
    print("\n" + "=" * 70)
    print("[3] Running scTDRP")
    print("=" * 70)

    # Use original normal cells as reference (no metacell aggregation)
    # Metacell averaging can blur stage boundaries when only a few stages exist.
    adata_meta_normal = adata_normal
    print(f"  Normal reference cells used directly: {adata_meta_normal.shape}")
    print(adata_meta_normal.obs["stage"].value_counts())

    # Use disease cells directly (n=1091 is small enough for OT)
    adata_meta_disease = adata_disease
    print(f"  Disease cells used directly: {adata_meta_disease.shape}")

    # Initialize analyzer
    analyzer = TDRPAnalyzer(
        normal_adata=adata_meta_normal,
        stage_key="stage",
        terminal_stage=TERMINAL_STAGE,
        stage_order=STAGE_ORDER,
        use_rep="X_pca",
        n_top_genes=2000,
    )
    analyzer.prepare_data(flavor="seurat_v3")
    analyzer.build_ot_cost_map(metric="sqeuclidean", reg=0.01)

    # Compute TDI
    tdi_df = analyzer.compute_tdi(adata_meta_disease, metric="sqeuclidean", reg=0.01)

    # Direct cell-level TDI (no metacell mapping needed)
    adata_disease.obs["TDI"] = tdi_df["TDI"].values
    adata_disease.obs["Best_Match_Stage"] = tdi_df["Best_Match_Stage"].values

    # Infer repair pathway in gene space using HVGs
    hvg_genes = adata_disease.var_names[adata_disease.var["highly_variable"]].tolist()
    analyzer.gene_list = hvg_genes
    analyzer.use_rep = "X"
    repair = analyzer.infer_repair_pathway(
        adata_meta_disease, metric="sqeuclidean", reg=0.01, top_n=100
    )

    # Module strategy: neural differentiation vs proliferation
    module_dict = {
        "OPC_Progenitor": ["PDGFRA", "CSPG4", "SOX10", "OLIG1", "OLIG2"],
        "Astrocyte_Differentiation": ["GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A3"],
        "Oligodendrocyte_Differentiation": ["MBP", "PLP1", "MOG", "MAG", "CNP"],
        "CellCycle": ["MKI67", "TOP2A", "CDK1", "CCNB1", "AURKA"],
    }
    module_strategy = analyzer.module_repair_strategy(module_dict, threshold=0.0, method="mean")

    return analyzer, tdi_df, repair, module_strategy, adata_meta_disease


def evaluate_and_plot(adata_disease, analyzer, module_strategy):
    """Compute summary statistics and generate figures."""
    print("\n" + "=" * 70)
    print("[4] Evaluation and visualization")
    print("=" * 70)

    stage_counts = adata_disease.obs["Best_Match_Stage"].value_counts(normalize=True)
    print("\n  Best-match stage distribution (neoplastic cells):")
    for stage, prop in stage_counts.items():
        print(f"    {stage}: {prop*100:.1f}%")

    frac_opc = stage_counts.get(EXPECTED_ARREST_STAGE, 0.0)
    print(f"\n  Fraction assigned to expected arrest stage ({EXPECTED_ARREST_STAGE}): {frac_opc*100:.1f}%")

    tdi_mean = adata_disease.obs["TDI"].mean()
    tdi_median = adata_disease.obs["TDI"].median()
    print(f"  Mean TDI: {tdi_mean:.4f}")
    print(f"  Median TDI: {tdi_median:.4f}")

    # Save cell-level results
    out_csv = os.path.join(OUT_DIR, "gbm_cell_tdi.csv")
    adata_disease.obs[["cell_type", "neoplastic", "patient_id", "TDI", "Best_Match_Stage"]].to_csv(out_csv)
    print(f"\n  Saved: {out_csv}")

    # Save module strategy
    if module_strategy is not None:
        module_csv = os.path.join(OUT_DIR, "gbm_module_strategy.csv")
        module_strategy.to_csv(module_csv, index=False)
        print(f"  Saved: {module_csv}")
        print("\n  Module-level repair strategy:")
        print(module_strategy.to_string(index=False))

    # Save h5ad (fill NaN stage to avoid h5py string conversion error)
    h5ad_path = os.path.join(OUT_DIR, "gbm_scTDRP_results.h5ad")
    adata_disease.obs["stage"] = adata_disease.obs["stage"].fillna("Disease")
    adata_disease.write_h5ad(h5ad_path)
    print(f"  Saved: {h5ad_path}")

    # Figures
    fig_stage = os.path.join(FIG_DIR, "gbm_stage_distribution.pdf")
    fig, ax = plt.subplots(figsize=(6, 4))
    stage_counts.plot(kind="bar", ax=ax, color="steelblue")
    ax.set_ylabel("Fraction of neoplastic cells")
    ax.set_title("scTDRP stage assignment (GBM)")
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    plt.savefig(fig_stage, dpi=300)
    plt.close()
    print(f"  Saved: {fig_stage}")

    fig_tdi = os.path.join(FIG_DIR, "gbm_tdi_distribution.pdf")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(adata_disease.obs["TDI"].dropna(), bins=30, color="steelblue", edgecolor="black")
    ax.axvline(tdi_median, color="red", linestyle="--", label=f"median={tdi_median:.3f}")
    ax.set_xlabel("TDI")
    ax.set_ylabel("Number of cells")
    ax.set_title("TDI distribution (GBM neoplastic cells)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(fig_tdi, dpi=300)
    plt.close()
    print(f"  Saved: {fig_tdi}")

    fig_cost = os.path.join(FIG_DIR, "gbm_ot_cost_map.pdf")
    analyzer.plot_ot_cost_map(save_path=fig_cost)
    print(f"  Saved: {fig_cost}")

    if module_strategy is not None:
        fig_module = os.path.join(FIG_DIR, "gbm_module_strategy.pdf")
        analyzer.plot_module_strategy(save_path=fig_module)
        print(f"  Saved: {fig_module}")


def main():
    adata = load_gbm_data()
    adata = preprocess_gbm(adata)
    adata_normal, adata_disease = build_reference_and_disease(adata)
    analyzer, tdi_df, repair, module_strategy, adata_meta_disease = run_scTDRP_gbm(
        adata_normal, adata_disease
    )
    evaluate_and_plot(adata_disease, analyzer, module_strategy)
    print("\n" + "=" * 70)
    print("GBM validation complete")
    print("=" * 70)


if __name__ == "__main__":
    main()
