"""可视化工具"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict


def plot_ot_cost_map(
    cost_map: Dict,
    figsize: tuple = (10, 6),
    color: str = "steelblue",
    save_path: str = None,
):
    """
    绘制OT代价地图
    
    Parameters
    ----------
    cost_map : dict
        {(s1, s2): distance}
    """
    stages = []
    distances = []
    for (s1, s2), dist in cost_map.items():
        stages.append(f"{s1}\n→\n{s2}")
        distances.append(dist)
    
    fig, ax = plt.subplots(figsize=figsize)
    bars = ax.bar(range(len(stages)), distances, color=color, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages, fontsize=9)
    ax.set_ylabel("Wasserstein-2 Distance", fontsize=12)
    ax.set_title("OT Cost Map: Normal Differentiation Trajectory", fontsize=14, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # 标注数值
    for bar, dist in zip(bars, distances):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{dist:.3f}", ha="center", va="bottom", fontsize=9)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_tdi_distribution(
    tdi_df: pd.DataFrame,
    figsize: tuple = (12, 5),
    bins: int = 50,
    color: str = "coral",
    save_path: str = None,
):
    """
    绘制TDI分布及最优归属阶段
    
    Parameters
    ----------
    tdi_df : pd.DataFrame
        TDRPAnalyzer.compute_tdi() 的输出
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # 左图：TDI直方图
    ax1 = axes[0]
    ax1.hist(tdi_df["TDI"], bins=bins, color=color, edgecolor="black", alpha=0.7)
    ax1.axvline(tdi_df["TDI"].mean(), color="darkred", linestyle="--", linewidth=2,
                label=f"Mean: {tdi_df['TDI'].mean():.3f}")
    ax1.axvline(tdi_df["TDI"].median(), color="navy", linestyle="--", linewidth=2,
                label=f"Median: {tdi_df['TDI'].median():.3f}")
    ax1.set_xlabel("Trajectory Deviation Index (TDI)", fontsize=12)
    ax1.set_ylabel("Number of Cells", fontsize=12)
    ax1.set_title("TDI Distribution", fontsize=13, fontweight="bold")
    ax1.legend()
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # 右图：最优归属阶段计数
    ax2 = axes[1]
    stage_counts = tdi_df["Best_Match_Stage"].value_counts()
    colors = sns.color_palette("husl", len(stage_counts))
    ax2.barh(stage_counts.index, stage_counts.values, color=colors, edgecolor="black")
    ax2.set_xlabel("Number of Cells", fontsize=12)
    ax2.set_title("Best Match Stage Distribution", fontsize=13, fontweight="bold")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_repair_heatmap(
    repair_results: Dict,
    top_n: int = 50,
    figsize: tuple = (8, 10),
    cmap: str = "RdBu_r",
    center: float = 0,
    save_path: str = None,
):
    """
    绘制Top修复靶点热图
    
    Parameters
    ----------
    repair_results : dict
        TDRPAnalyzer.infer_repair_pathway() 的输出
    """
    deltas = repair_results.get("repair_deltas", {})
    if not deltas:
        print("No gene-level repair deltas available for plotting.")
        return None
    
    df = pd.DataFrame(list(deltas.items()), columns=["Gene", "Delta"])
    df = df.sort_values("Delta", ascending=False)
    
    # 取Top上调和Top下调
    top_up = df.head(top_n // 2).copy()
    top_down = df.tail(top_n // 2).copy()
    plot_df = pd.concat([top_up, top_down])
    plot_df = plot_df.sort_values("Delta", ascending=True)
    
    fig, ax = plt.subplots(figsize=figsize)
    colors = plot_df["Delta"].values
    norm = plt.Normalize(vmin=colors.min(), vmax=colors.max())
    
    bars = ax.barh(range(len(plot_df)), plot_df["Delta"].values,
                   color=plt.cm.RdBu_r(norm(colors)), edgecolor="black", linewidth=0.3)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["Gene"].values, fontsize=8)
    ax.set_xlabel("Repair Delta (Disease → Terminal)", fontsize=12)
    ax.set_title(f"Top {top_n} Repair Targets", fontsize=13, fontweight="bold")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_module_strategy(
    module_strategy_df: pd.DataFrame,
    figsize: tuple = (10, 6),
    save_path: str = None,
):
    """
    绘制模块修复策略图
    
    Parameters
    ----------
    module_strategy_df : pd.DataFrame
        TDRPAnalyzer.module_repair_strategy() 的输出
    """
    df = module_strategy_df.sort_values("Repair_Score", ascending=True)
    
    fig, ax = plt.subplots(figsize=figsize)
    colors = ["#e74c3c" if s > 0 else "#3498db" if s < 0 else "#95a5a6"
              for s in df["Repair_Score"]]
    
    bars = ax.barh(df["Module"], df["Repair_Score"], color=colors, edgecolor="black")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Module Repair Score", fontsize=12)
    ax.set_title("Module-Level Repair Strategy", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    # 图例
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#e74c3c", label="Up-regulate"),
        Patch(facecolor="#3498db", label="Down-regulate"),
        Patch(facecolor="#95a5a6", label="Unchanged"),
    ]
    ax.legend(handles=legend_elements, loc="lower right")
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig


def plot_stage_distance_heatmap(
    tdi_df: pd.DataFrame,
    stage_cols: list = None,
    figsize: tuple = (10, 8),
    cmap: str = "YlOrRd",
    save_path: str = None,
):
    """
    绘制细胞到各阶段的距离热图（样本）
    
    Parameters
    ----------
    tdi_df : pd.DataFrame
        包含 Dist_* 列的TDI结果
    """
    if stage_cols is None:
        stage_cols = [c for c in tdi_df.columns if c.startswith("Dist_")]
    
    dist_matrix = tdi_df[stage_cols].values
    # 随机采样200个细胞避免过大
    if dist_matrix.shape[0] > 200:
        idx = np.random.choice(dist_matrix.shape[0], 200, replace=False)
        dist_matrix = dist_matrix[idx, :]
    
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(dist_matrix, cmap=cmap, ax=ax, cbar_kws={"label": "Wasserstein Distance"})
    ax.set_xticklabels([c.replace("Dist_", "") for c in stage_cols], rotation=45, ha="right")
    ax.set_ylabel("Cells (sampled)", fontsize=12)
    ax.set_title("Cell-to-Stage Distance Heatmap", fontsize=13, fontweight="bold")
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    return fig
