"""数据读写工具"""

import os
import json
import numpy as np
import pandas as pd
from anndata import AnnData
import scanpy as sc


def load_h5ad(filepath: str) -> AnnData:
    """加载 h5ad 文件"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    return sc.read_h5ad(filepath)


def save_results(results: dict, output_dir: str, prefix: str = "scTDRP") -> None:
    """保存分析结果到指定目录"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存 TDI
    if "tdi" in results:
        tdi_path = os.path.join(output_dir, f"{prefix}_tdi.csv")
        results["tdi"].to_csv(tdi_path, index=False)
        print(f"TDI saved to {tdi_path}")
    
    # 保存修复靶点
    if "repair" in results and "repair_deltas" in results["repair"]:
        deltas = results["repair"]["repair_deltas"]
        if deltas:
            df = pd.DataFrame(list(deltas.items()), columns=["Gene", "Repair_Delta"])
            df = df.sort_values("Repair_Delta", ascending=False)
            repair_path = os.path.join(output_dir, f"{prefix}_repair_deltas.csv")
            df.to_csv(repair_path, index=False)
            print(f"Repair deltas saved to {repair_path}")
    
    # 保存模块策略
    if "module_strategy" in results:
        module_path = os.path.join(output_dir, f"{prefix}_module_strategy.csv")
        results["module_strategy"].to_csv(module_path, index=False)
        print(f"Module strategy saved to {module_path}")
    
    # 保存代价地图
    if "cost_map" in results:
        cost_path = os.path.join(output_dir, f"{prefix}_cost_map.json")
        # 将tuple key转为str key
        cost_map_serializable = {f"{k[0]}->{k[1]}": float(v) for k, v in results["cost_map"].items()}
        with open(cost_path, "w") as f:
            json.dump(cost_map_serializable, f, indent=2)
        print(f"Cost map saved to {cost_path}")


def load_module_dict(filepath: str) -> dict:
    """
    加载模块定义文件
    
    支持格式:
    - JSON: {"module_name": ["gene1", "gene2", ...]}
    - GMT: 标准GMT格式
    - CSV: 两列，Module 和 Gene
    """
    ext = os.path.splitext(filepath)[1].lower()
    
    if ext == ".json":
        with open(filepath, "r") as f:
            return json.load(f)
    
    elif ext == ".gmt":
        module_dict = {}
        with open(filepath, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    module_name = parts[0]
                    genes = parts[2:]
                    module_dict[module_name] = genes
        return module_dict
    
    elif ext == ".csv":
        df = pd.read_csv(filepath)
        if "Module" not in df.columns or "Gene" not in df.columns:
            raise ValueError("CSV must contain 'Module' and 'Gene' columns")
        module_dict = {}
        for module, group in df.groupby("Module"):
            module_dict[module] = group["Gene"].tolist()
        return module_dict
    
    else:
        raise ValueError(f"Unsupported file format: {ext}")


def save_transport_plan(transport_plan: np.ndarray, output_path: str) -> None:
    """保存传输计划为npy文件"""
    np.save(output_path, transport_plan)
    print(f"Transport plan saved to {output_path}")
