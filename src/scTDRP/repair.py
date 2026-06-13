"""修复路径推断"""

import numpy as np
import ot
from anndata import AnnData
from typing import Dict

from .utils import compute_expression_distribution
from .distance import compute_wasserstein_distance


def compute_gene_repair_delta(
    X_disease: np.ndarray,
    X_terminal: np.ndarray,
    transport_plan: np.ndarray,
    gene_names: list = None,
) -> Dict[str, float]:
    """
    基于传输计划计算每个基因的修复量
    
    Parameters
    ----------
    X_disease : np.ndarray (n_disease_cells, n_genes)
        疾病细胞表达矩阵
    X_terminal : np.ndarray (n_terminal_cells, n_genes)
        正常终末阶段细胞表达矩阵
    transport_plan : np.ndarray (n_disease_cells, n_terminal_cells)
        最优传输计划
    gene_names : list, optional
        基因名称列表
    
    Returns
    -------
    dict
        {gene_name: delta} 正值表示需上调，负值表示需下调
    """
    n_genes = X_disease.shape[1]
    if gene_names is None:
        gene_names = [f"gene_{i}" for i in range(n_genes)]
    
    repair_deltas = {}
    for g in range(n_genes):
        # 计算基因g的预期表达改变
        delta = 0.0
        for i in range(X_disease.shape[0]):
            for j in range(X_terminal.shape[0]):
                delta += transport_plan[i, j] * (X_terminal[j, g] - X_disease[i, g])
        gene_name = gene_names[g] if g < len(gene_names) else f"gene_{g}"
        repair_deltas[gene_name] = delta
    
    return repair_deltas


def infer_repair_pathway(
    adata_disease: AnnData,
    adata_terminal: AnnData,
    gene_list: list = None,
    use_rep: str = "X_pca",
    metric: str = "sqeuclidean",
    reg: float = 0.01,
    top_n: int = 100,
) -> Dict:
    """
    推断从疾病状态到正常终末状态的修复路径
    
    Parameters
    ----------
    adata_disease : AnnData
        疾病细胞
    adata_terminal : AnnData
        正常终末阶段细胞
    gene_list : list, optional
        分析用的基因列表
    use_rep : str
        表示层
    metric : str
        OT距离度量
    reg : float
        Sinkhorn正则化参数
    top_n : int
        返回Top N修复靶点
    
    Returns
    -------
    dict
        {
            "wasserstein_distance": float,
            "transport_plan": np.ndarray,
            "repair_deltas": dict,
            "top_up_targets": list,
            "top_down_targets": list,
        }
    """
    # 提取表达矩阵
    X_disease = compute_expression_distribution(adata_disease, gene_list=gene_list, use_rep=use_rep)
    X_terminal = compute_expression_distribution(adata_terminal, gene_list=gene_list, use_rep=use_rep)
    
    # 计算OT距离和传输计划
    dist, transport_plan = compute_wasserstein_distance(
        X_disease, X_terminal, metric=metric, reg=reg, method="sinkhorn"
    )
    
    # 获取基因名
    if gene_list is not None:
        gene_names = gene_list
    elif use_rep == "X" or use_rep is None:
        gene_names = adata_disease.var_names.tolist()
    else:
        gene_names = None  # PCA空间无法对应到具体基因
    
    # 计算基因修复量
    if gene_names is not None and use_rep in ("X", None):
        repair_deltas = compute_gene_repair_delta(X_disease, X_terminal, transport_plan, gene_names)
        
        # 排序
        sorted_genes = sorted(repair_deltas.items(), key=lambda x: x[1], reverse=True)
        top_up = [g for g, d in sorted_genes[:top_n] if d > 0]
        top_down = [g for g, d in sorted_genes[-top_n:] if d < 0]
    else:
        repair_deltas = {}
        top_up = []
        top_down = []
    
    return {
        "wasserstein_distance": dist,
        "transport_plan": transport_plan,
        "repair_deltas": repair_deltas,
        "top_up_targets": top_up,
        "top_down_targets": top_down,
    }


def compute_cell_level_tdi(
    X_cell: np.ndarray,
    stage_distributions: Dict[str, tuple],
    metric: str = "sqeuclidean",
    reg: float = 0.01,
) -> Dict[str, float]:
    """
    计算单个细胞到各正常阶段的OT距离（TDI的组成部分）
    
    Parameters
    ----------
    X_cell : np.ndarray (1, n_features) or (n_features,)
        单个细胞的表达向量
    stage_distributions : dict
        {stage: (X, weights)}
    
    Returns
    -------
    dict
        {stage_name: distance}
    """
    if X_cell.ndim == 1:
        X_cell = X_cell.reshape(1, -1)
    
    distances = {}
    for stage, (X_stage, w_stage) in stage_distributions.items():
        w_cell = np.ones(1)
        dist, _ = compute_wasserstein_distance(
            X_cell, X_stage, w_cell, w_stage, metric=metric, reg=reg
        )
        distances[stage] = dist
    
    return distances
