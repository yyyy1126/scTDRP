#!/usr/bin/env python3
"""
Generate Figure 1: scTDRP method overview and erythroid leukemia results.

Panels:
  A. Workflow schematic
  B. OT cost map between normal erythroid stages
  C. TDI distribution (malignant vs stage-matched normal)
  D. Repair heatmap (top genes)
  E. Module repair strategy
"""

import os
import json
import warnings
from collections import Counter

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0
np.random.seed(42)

AEL_PATH = "/data1/yja/zhongzhuan/4.results/4.infercnv/infercnv_aml5/aml5_annotated_with_cnv_cycle.h5ad"
NORMAL_PATH = "/data1/yja/zhongzhuan/1.data/processed/erythroid_lineage_from_MEP.h5ad"
MODULES_PATH = "/data1/yja/zhongzhuan/5.external/scTDRP/modules.json"
OUT_DIR = "/data1/yja/zhongzhuan/5.external/scTDRP/figures_figure1"
os.makedirs(OUT_DIR, exist_ok=True)

STAGE_ORDER = [
    'MEP', 'BFU-E', 'CFU-E', 'Pro-Erythroblast',
    'Basophilic Erythroblast', 'Polychromatic Erythroblast',
    'Orthochromatic Erythroblast'
]
TERMINAL_STAGE = 'Orthochromatic Erythroblast'
EXPECTED_ARREST_STAGE = 'Polychromatic Erythroblast'

import sys
try:
    from scTDRP import TDRPAnalyzer, build_stage_metacells
except Exception:
    sys.path.insert(0, "/data1/yja/zhongzhuan/5.external/scTDRP/src")
    from scTDRP import TDRPAnalyzer, build_stage_metacells


def load_and_preprocess():
    """Load AEL and normal erythroid data with unified gene space and PCA."""
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    adata_ael = sc.read_h5ad(AEL_PATH)
    adata_normal = sc.read_h5ad(NORMAL_PATH)

    # Normal gene mapping: Ensembl -> HGNC symbol
    adata_normal.var['gene_symbol'] = adata_normal.var['gene_symbols'].astype(str)
    mean_expr = np.array(adata_normal.X.mean(axis=0)).flatten()
    if hasattr(mean_expr, 'toarray'):
        mean_expr = mean_expr.toarray().flatten()
    adata_normal.var['_mean_expr'] = mean_expr

    keep_idx = []
    for sym, group in adata_normal.var.groupby('gene_symbol'):
        if len(group) > 1:
            keep_idx.append(group['_mean_expr'].idxmax())
        else:
            keep_idx.append(group.index[0])
    adata_normal = adata_normal[:, keep_idx].copy()
    adata_normal.var_names = adata_normal.var['gene_symbol'].values
    adata_normal.var_names_make_unique()

    common_genes = list(set(adata_ael.var_names) & set(adata_normal.var_names))
    adata_ael = adata_ael[:, common_genes].copy()
    adata_normal = adata_normal[:, common_genes].copy()

    if 'logcounts' in adata_ael.layers:
        adata_ael.X = adata_ael.layers['logcounts']

    mask_malignant = adata_ael.obs['malignancy'] == 'Malignant Erythroid'
    adata_malignant = adata_ael[mask_malignant].copy()
    adata_residual = adata_ael[~mask_malignant].copy()

    adata_joint = sc.concat(
        [adata_normal, adata_malignant], label='_source',
        keys=['normal', 'disease'], index_unique='-', join='outer'
    )
    sc.pp.normalize_total(adata_joint, target_sum=1e4)
    sc.pp.log1p(adata_joint)
    sc.pp.highly_variable_genes(adata_joint, n_top_genes=2000, flavor='seurat_v3')

    hvg_mask = adata_joint.var['highly_variable'].values
    adata_normal.var['highly_variable'] = hvg_mask
    adata_malignant.var['highly_variable'] = hvg_mask
    adata_residual.var['highly_variable'] = hvg_mask

    jmask_normal = adata_joint.obs['_source'] == 'normal'
    jmask_disease = adata_joint.obs['_source'] == 'disease'

    X_fit = np.asarray(adata_joint[:, hvg_mask].X.toarray() if hasattr(adata_joint.X, 'toarray') else adata_joint[:, hvg_mask].X)
    scaler = StandardScaler()
    X_fit_scaled = scaler.fit_transform(X_fit)
    pca = PCA(n_components=50, svd_solver='arpack')
    pca.fit(X_fit_scaled)

    def project_to_pca(adata_subset):
        X = np.asarray(adata_subset[:, hvg_mask].X.toarray() if hasattr(adata_subset.X, 'toarray') else adata_subset[:, hvg_mask].X)
        X_scaled = scaler.transform(X)
        return pca.transform(X_scaled)

    adata_normal.obsm['X_pca'] = project_to_pca(adata_joint[jmask_normal])
    adata_malignant.obsm['X_pca'] = project_to_pca(adata_joint[jmask_disease])
    adata_residual.obsm['X_pca'] = project_to_pca(adata_residual)

    if 'AuthorCellType' in adata_normal.obs.columns:
        adata_normal.obs['stage'] = adata_normal.obs['AuthorCellType'].values
    else:
        adata_normal.obs['stage'] = 'Unknown'

    return adata_normal, adata_malignant, adata_residual


def run_scTDRP(adata_normal, adata_malignant):
    """Run scTDRP and return results."""
    adata_meta_normal = build_stage_metacells(
        adata_normal,
        stage_key='stage',
        stage_order=STAGE_ORDER,
        use_rep='X_pca',
        resolution_scale=1.0,
        target_metacells=20,
    )

    analyzer = TDRPAnalyzer(
        normal_adata=adata_meta_normal,
        stage_key='stage',
        terminal_stage=TERMINAL_STAGE,
        stage_order=STAGE_ORDER,
        use_rep='X_pca',
        n_top_genes=2000,
    )
    analyzer.prepare_data(flavor='seurat_v3')
    analyzer.build_ot_cost_map(metric='sqeuclidean', reg=0.01)
    tdi_df = analyzer.compute_tdi(adata_malignant, metric='sqeuclidean', reg=0.01)

    adata_malignant.obs['TDI'] = tdi_df['TDI'].values
    adata_malignant.obs['Best_Match_Stage'] = tdi_df['Best_Match_Stage'].values

    # Module repair
    with open(MODULES_PATH) as f:
        modules = json.load(f)
    module_genes = set()
    for glist in modules.values():
        module_genes.update(glist)
    module_genes = [g for g in module_genes if g in adata_malignant.var_names]

    analyzer.gene_list = module_genes
    analyzer.use_rep = 'X'
    repair = analyzer.infer_repair_pathway(adata_malignant, metric='sqeuclidean', reg=0.01, top_n=50)

    module_strategy = analyzer.module_repair_strategy(modules, threshold=0.0, method='mean')

    return analyzer, tdi_df, repair, module_strategy


def draw_workflow_panel(ax):
    """Draw panel A: workflow schematic."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    def box(x, y, w, h, text, color):
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.3",
                              facecolor=color, edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=11, fontweight='bold')

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", lw=2.5, color='#333333'))

    # Vertical workflow on the left side
    box(1.0, 7.5, 2.5, 1.3, "Normal\nreference", "#E8F4F8")
    box(1.0, 5.0, 2.5, 1.3, "Disease\ncells", "#FFE8E8")
    arrow(2.25, 7.5, 2.25, 6.3)
    arrow(2.25, 6.2, 2.25, 5.0)

    box(4.5, 6.0, 2.5, 1.3, "OT cost\nmap", "#FFF4E6")
    arrow(3.5, 6.3, 4.5, 6.6)

    box(7.5, 6.0, 2.5, 1.3, "TDI +\nstage assign", "#E8F8E8")
    arrow(7.0, 6.6, 7.5, 6.6)

    box(7.5, 3.0, 2.5, 1.3, "Repair\npathway", "#F0E8F8")
    arrow(8.75, 6.0, 8.75, 4.3)

    box(4.5, 3.0, 2.5, 1.3, "Module\nstrategy", "#F8F0E8")
    arrow(7.5, 3.6, 7.0, 3.6)

    # Side labels
    ax.text(5.0, 0.8, "Quantify deviation  →  Identify arrest  →  Infer repair",
            ha='center', fontsize=11, style='italic', color='#444444')

    ax.set_title("A. scTDRP workflow", fontsize=13, fontweight='bold', loc='left')


def draw_cost_map_panel(ax, analyzer):
    """Draw panel B: OT cost map."""
    cost_map = analyzer.cost_map
    n = len(STAGE_ORDER)
    mat = np.zeros((n, n))
    for i, s1 in enumerate(STAGE_ORDER):
        for j, s2 in enumerate(STAGE_ORDER):
            if (s1, s2) in cost_map:
                mat[i, j] = cost_map[(s1, s2)]
            elif i == j:
                mat[i, j] = 0
            else:
                mat[i, j] = np.nan

    # Build annotation matrix: only show adjacent-stage distances
    annot = np.full_like(mat, "", dtype=object)
    for (s1, s2), dist in cost_map.items():
        i = STAGE_ORDER.index(s1)
        j = STAGE_ORDER.index(s2)
        annot[i, j] = f"{dist:.3f}"

    sns.heatmap(mat, mask=np.isnan(mat), annot=annot, fmt="", cmap="YlOrRd",
                xticklabels=[s.replace(' Erythroblast', '') for s in STAGE_ORDER],
                yticklabels=[s.replace(' Erythroblast', '') for s in STAGE_ORDER],
                ax=ax, cbar_kws={'label': 'W2 distance'}, vmin=0)
    ax.set_title("B. OT cost map (normal erythropoiesis)", fontsize=12, fontweight='bold', loc='left')


def draw_tdi_panel(ax, adata_malignant, adata_normal, analyzer):
    """Draw panel C: TDI distribution."""
    # TDI for malignant cells
    mal_tdi = adata_malignant.obs['TDI'].dropna().values

    # Stage-matched normal TDI: for each normal stage, compute TDI to other stages
    stage_tdi = {}
    for stage in STAGE_ORDER:
        mask = adata_normal.obs['stage'] == stage
        if mask.sum() == 0:
            continue
        stage_cells = adata_normal[mask]
        tdi = analyzer.compute_tdi(stage_cells, metric='sqeuclidean', reg=0.01)
        stage_tdi[stage] = tdi['TDI'].values

    # Violin plot
    data_for_plot = []
    labels = []
    data_for_plot.append(mal_tdi)
    labels.append('AEL malignant')
    for stage in STAGE_ORDER:
        if stage in stage_tdi:
            data_for_plot.append(stage_tdi[stage])
            labels.append(stage.replace(' Erythroblast', ''))

    parts = ax.violinplot(data_for_plot, positions=range(len(data_for_plot)),
                          showmeans=True, showmedians=False, widths=0.7)
    colors = ['#D62728'] + ['#1F77B4'] * len(STAGE_ORDER)
    for pc, c in zip(parts['bodies'], colors):
        pc.set_facecolor(c)
        pc.set_alpha(0.6)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha='right')
    ax.set_ylabel('TDI')
    ax.set_title('C. TDI: malignant vs stage-matched normal', fontsize=12, fontweight='bold', loc='left')


def draw_repair_heatmap(ax, repair):
    """Draw panel D: repair heatmap (top up/down genes)."""
    deltas = repair['repair_deltas']
    sorted_genes = sorted(deltas.items(), key=lambda x: x[1])
    top_down = sorted_genes[:15]
    top_up = sorted_genes[-15:][::-1]

    genes = [g for g, _ in top_down + top_up]
    values = [d for _, d in top_down + top_up]

    vmax = max(abs(np.min(values)), abs(np.max(values)))
    cmap = sns.diverging_palette(250, 15, s=75, l=40, n=9, center="light", as_cmap=True)
    sns.heatmap(np.array(values).reshape(-1, 1), cmap=cmap, center=0,
                vmin=-vmax, vmax=vmax,
                yticklabels=genes, xticklabels=['Repair Δ'], ax=ax,
                cbar_kws={'label': 'Δ expression'})
    ax.set_title('D. Gene-level repair targets', fontsize=12, fontweight='bold', loc='left')


def draw_module_panel(ax, module_strategy):
    """Draw panel E: module strategy."""
    df = module_strategy.copy()
    colors = ['#2CA02C' if s == '上调 (Up-regulate)' else '#D62728' for s in df['Strategy']]
    ax.barh(df['Module'], df['Repair_Score'], color=colors, edgecolor='black')
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Repair score')
    ax.set_title('E. Module-level repair strategy', fontsize=12, fontweight='bold', loc='left')
    ax.invert_yaxis()


def main():
    print("Loading and preprocessing AEL data...")
    adata_normal, adata_malignant, adata_residual = load_and_preprocess()

    print("Running scTDRP...")
    analyzer, tdi_df, repair, module_strategy = run_scTDRP(adata_normal, adata_malignant)

    # Reset analyzer to PCA mode for TDI panel
    analyzer.use_rep = 'X_pca'
    analyzer.gene_list = None

    print("Generating Figure 1...")
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.35, wspace=0.45)

    ax_a = fig.add_subplot(gs[:, 0])  # workflow spans both rows
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[1, 2])

    draw_workflow_panel(ax_a)
    draw_cost_map_panel(ax_b, analyzer)
    draw_tdi_panel(ax_c, adata_malignant, adata_normal, analyzer)
    draw_repair_heatmap(ax_d, repair)
    draw_module_panel(ax_e, module_strategy)

    plt.suptitle('Figure 1. scTDRP: single-cell Transcriptomic Developmental Repair Potential',
                 fontsize=14, fontweight='bold', y=0.98)

    out_pdf = os.path.join(OUT_DIR, 'Figure1_scTDRP_overview.pdf')
    out_png = os.path.join(OUT_DIR, 'Figure1_scTDRP_overview.png')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
