"""模块级修复策略"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple


def aggregate_repair_to_module(
    repair_deltas: Dict[str, float],
    module_dict: Dict[str, List[str]],
    method: str = "mean",
) -> Dict[str, float]:
    """
    将基因级修复量聚合到模块级
    
    Parameters
    ----------
    repair_deltas : dict
        {gene_name: delta}
    module_dict : dict
        {module_name: [gene1, gene2, ...]}
    method : str
        "mean" 或 "sum"
    
    Returns
    -------
    dict
        {module_name: module_repair_score}
    """
    module_scores = {}
    
    for module_name, gene_list in module_dict.items():
        deltas = []
        for gene in gene_list:
            if gene in repair_deltas:
                deltas.append(repair_deltas[gene])
        
        if len(deltas) == 0:
            module_scores[module_name] = 0.0
            continue
        
        if method == "mean":
            module_scores[module_name] = np.mean(deltas)
        elif method == "sum":
            module_scores[module_name] = np.sum(deltas)
        elif method == "median":
            module_scores[module_name] = np.median(deltas)
        else:
            raise ValueError(f"Unknown method: {method}")
    
    return module_scores


def module_repair_strategy(
    repair_deltas: Dict[str, float],
    module_dict: Dict[str, List[str]],
    threshold: float = 0.0,
    method: str = "mean",
) -> pd.DataFrame:
    """
    生成模块级修复策略表
    
    Parameters
    ----------
    repair_deltas : dict
    module_dict : dict
    threshold : float
        判定显著性的阈值
    method : str
    
    Returns
    -------
    pd.DataFrame
        模块修复策略表，包含模块名、修复分值、策略（上调/下调/不变）
    """
    module_scores = aggregate_repair_to_module(repair_deltas, module_dict, method=method)
    
    records = []
    for module_name, score in module_scores.items():
        if score > threshold:
            strategy = "上调 (Up-regulate)"
        elif score < -threshold:
            strategy = "下调 (Down-regulate)"
        else:
            strategy = "不变 (Unchanged)"
        
        # 统计模块中覆盖的基因数
        n_genes = len(module_dict.get(module_name, []))
        n_covered = sum(1 for g in module_dict.get(module_name, []) if g in repair_deltas)
        
        records.append({
            "Module": module_name,
            "Repair_Score": score,
            "Strategy": strategy,
            "N_Genes_Total": n_genes,
            "N_Genes_Covered": n_covered,
            "Coverage_Rate": n_covered / n_genes if n_genes > 0 else 0,
        })
    
    df = pd.DataFrame(records)
    df = df.sort_values("Repair_Score", ascending=False)
    return df


def compare_module_repair_across_conditions(
    condition_repair_deltas: Dict[str, Dict[str, float]],
    module_dict: Dict[str, List[str]],
    method: str = "mean",
) -> pd.DataFrame:
    """
    比较多个条件下的模块修复策略
    
    Parameters
    ----------
    condition_repair_deltas : dict
        {condition_name: {gene: delta}}
    module_dict : dict
    method : str
    
    Returns
    -------
    pd.DataFrame
        长格式表格，用于绘制分组柱状图
    """
    records = []
    for condition, deltas in condition_repair_deltas.items():
        module_scores = aggregate_repair_to_module(deltas, module_dict, method=method)
        for module, score in module_scores.items():
            records.append({
                "Condition": condition,
                "Module": module,
                "Repair_Score": score,
            })
    
    return pd.DataFrame(records)
