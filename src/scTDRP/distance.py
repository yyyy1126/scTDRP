"""Wasserstein 距离计算与阶段分布构建"""

import numpy as np
import ot
from anndata import AnnData
from typing import Dict, Tuple

from .utils import compute_expression_distribution, get_stage_cells


def build_stage_distributions(
    adata: AnnData,
    stage_key: str,
    stage_order: list = None,
    gene_list: list = None,
    use_rep: str = "X_pca",
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    构建各分化阶段的表达分布
    
    Parameters
    ----------
    adata : AnnData
        正常参考数据
    stage_key : str
        obs列名，标识分化阶段
    stage_order : list, optional
        阶段顺序，如 ["MEP", "BFU-E", "CFU-E", ...]
    gene_list : list, optional
        使用的基因列表
    use_rep : str
        使用的表示层
    
    Returns
    -------
    dict
        {stage_name: (expressions, weights)} 
        expressions: (n_cells, n_features)
        weights: (n_cells,) 均匀权重或根据质量调整
    """
    if stage_key not in adata.obs.columns:
        raise KeyError(f"Stage key '{stage_key}' not found")
    
    stages = adata.obs[stage_key].unique()
    if stage_order is not None:
        stages = [s for s in stage_order if s in stages]
    
    stage_distributions = {}
    for stage in stages:
        stage_adata = get_stage_cells(adata, stage_key, stage)
        X = compute_expression_distribution(stage_adata, gene_list=gene_list, use_rep=use_rep)
        # 均匀权重，可扩展为根据细胞质量调整
        weights = np.ones(X.shape[0]) / X.shape[0]
        stage_distributions[stage] = (X, weights)
    
    return stage_distributions


def compute_wasserstein_distance(
    X_source: np.ndarray,
    X_target: np.ndarray,
    w_source: np.ndarray = None,
    w_target: np.ndarray = None,
    metric: str = "sqeuclidean",
    reg: float = 0.01,
    method: str = "sinkhorn",
) -> Tuple[float, np.ndarray]:
    """
    计算两个分布之间的 Wasserstein-2 距离及最优传输计划
    
    Parameters
    ----------
    X_source : np.ndarray (n_samples_1, n_features)
    X_target : np.ndarray (n_samples_2, n_features)
    w_source : np.ndarray (n_samples_1,), optional
    w_target : np.ndarray (n_samples_2,), optional
    metric : str
        距离度量，如 "sqeuclidean", "euclidean", "cosine"
    reg : float
        Sinkhorn 正则化参数
    method : str
        "sinkhorn" 或 "emd"
    
    Returns
    -------
    distance : float
        Wasserstein 距离
    transport_plan : np.ndarray (n_samples_1, n_samples_2)
        最优传输计划
    """
    n1, n2 = X_source.shape[0], X_target.shape[0]
    
    if w_source is None:
        w_source = np.ones(n1) / n1
    if w_target is None:
        w_target = np.ones(n2) / n2
    
    w_source = w_source.astype(np.float64)
    w_target = w_target.astype(np.float64)
    
    # 计算代价矩阵
    M = ot.dist(X_source, X_target, metric=metric)
    M = M.astype(np.float64)
    
    # 归一化代价矩阵（Sinkhorn稳定性）
    M /= M.max() if M.max() > 0 else 1.0
    
    if method == "sinkhorn":
        transport_plan = ot.sinkhorn(w_source, w_target, M, reg)
        distance = np.sum(transport_plan * M)
    elif method == "emd":
        transport_plan = ot.emd(w_source, w_target, M)
        distance = np.sum(transport_plan * M)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return distance, transport_plan


def compute_ot_cost_map(
    stage_distributions: Dict[str, Tuple[np.ndarray, np.ndarray]],
    metric: str = "sqeuclidean",
    reg: float = 0.01,
) -> Dict[Tuple[str, str], float]:
    """
    计算相邻阶段之间的 OT 代价地图
    
    Returns
    -------
    dict
        {(stage_i, stage_j): wasserstein_distance}
    """
    stages = list(stage_distributions.keys())
    cost_map = {}
    
    for i in range(len(stages) - 1):
        s1, s2 = stages[i], stages[i + 1]
        X1, w1 = stage_distributions[s1]
        X2, w2 = stage_distributions[s2]
        dist, _ = compute_wasserstein_distance(X1, X2, w1, w2, metric=metric, reg=reg)
        cost_map[(s1, s2)] = dist
    
    return cost_map
