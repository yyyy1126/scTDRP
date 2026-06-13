"""工具函数"""

import numpy as np
import scanpy as sc
from anndata import AnnData


def validate_anndata(adata: AnnData, required_obs_keys: list = None) -> None:
    """验证 AnnData 对象的基本结构"""
    if not isinstance(adata, AnnData):
        raise TypeError(f"Expected AnnData, got {type(adata)}")
    if adata.X is None:
        raise ValueError("AnnData.X is None")
    if required_obs_keys:
        for key in required_obs_keys:
            if key not in adata.obs.columns:
                raise KeyError(f"Required obs key '{key}' not found in adata.obs")


def extract_highly_variable_genes(
    adata: AnnData,
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    layer: str = None,
) -> list:
    """提取高变基因列表"""
    adata_copy = adata.copy()
    if layer is not None:
        sc.pp.highly_variable_genes(
            adata_copy, n_top_genes=n_top_genes, flavor=flavor, layer=layer
        )
    else:
        sc.pp.highly_variable_genes(
            adata_copy, n_top_genes=n_top_genes, flavor=flavor
        )
    return adata_copy.var_names[adata_copy.var["highly_variable"]].tolist()


def compute_expression_distribution(
    adata: AnnData,
    gene_list: list = None,
    use_rep: str = "X_pca",
) -> np.ndarray:
    """
    计算细胞群体的表达分布矩阵
    
    Parameters
    ----------
    adata : AnnData
    gene_list : list, optional
        指定基因列表，默认使用全部基因
    use_rep : str
        使用的表示层，如 'X_pca', 'X' 等
    
    Returns
    -------
    np.ndarray
        表达矩阵 (n_cells, n_features)
    """
    if use_rep == "X" or use_rep is None:
        if gene_list is not None:
            X = adata[:, gene_list].X
        else:
            X = adata.X
    else:
        if use_rep not in adata.obsm:
            raise KeyError(f"{use_rep} not found in adata.obsm")
        X = adata.obsm[use_rep]
    
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=np.float64)


def get_stage_cells(adata: AnnData, stage_key: str, stage_value: str) -> AnnData:
    """获取特定阶段的细胞子集"""
    if stage_key not in adata.obs.columns:
        raise KeyError(f"Stage key '{stage_key}' not found in adata.obs")
    mask = adata.obs[stage_key] == stage_value
    if mask.sum() == 0:
        raise ValueError(f"No cells found for stage '{stage_value}'")
    return adata[mask].copy()


def normalize_distribution_weights(weights: np.ndarray) -> np.ndarray:
    """归一化分布权重"""
    weights = np.asarray(weights, dtype=np.float64)
    if weights.sum() == 0:
        raise ValueError("Weights sum to zero")
    return weights / weights.sum()
