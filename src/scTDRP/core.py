"""
核心分析类 TDRPAnalyzer
整合OT代价地图、TDI计算、修复路径推断和模块级策略
"""

import numpy as np
import pandas as pd
from anndata import AnnData
from typing import Dict, List, Optional, Tuple
import warnings

from .utils import (
    validate_anndata,
    compute_expression_distribution,
    get_stage_cells,
    extract_highly_variable_genes,
)
from .distance import (
    build_stage_distributions,
    compute_wasserstein_distance,
    compute_ot_cost_map,
)
from .repair import (
    infer_repair_pathway,
    compute_cell_level_tdi,
)
from .module import module_repair_strategy
from .io import save_results
from . import plotting
from . import metacells


class TDRPAnalyzer:
    """
    scTDRP 主分析类
    
    Parameters
    ----------
    normal_adata : AnnData
        正常参考数据，需包含 stage_key 标注的分化阶段
    stage_key : str
        obs列名，标识分化阶段
    terminal_stage : str
        终末分化阶段名称
    stage_order : list, optional
        阶段顺序，如 ["MEP", "BFU-E", "CFU-E", ...]
    use_rep : str
        使用的表示层，如 "X_pca", "X"
    n_top_genes : int
        高变基因数量
    """
    
    def __init__(
        self,
        normal_adata: AnnData,
        stage_key: str,
        terminal_stage: str,
        stage_order: List[str] = None,
        use_rep: str = "X_pca",
        n_top_genes: int = 2000,
    ):
        validate_anndata(normal_adata, required_obs_keys=[stage_key])
        
        self.normal_adata = normal_adata.copy()
        self.stage_key = stage_key
        self.terminal_stage = terminal_stage
        self.use_rep = use_rep
        self.n_top_genes = n_top_genes
        
        # 确定阶段顺序
        available_stages = normal_adata.obs[stage_key].unique().tolist()
        if stage_order is not None:
            self.stage_order = [s for s in stage_order if s in available_stages]
        else:
            self.stage_order = available_stages
        
        if terminal_stage not in self.stage_order:
            raise ValueError(f"Terminal stage '{terminal_stage}' not found in stage_order")
        
        # 内部状态
        self.gene_list = None
        self.stage_distributions = None
        self.cost_map = None
        self.tdi_results = None
        self.repair_results = None
        self.module_strategy = None
        
    def prepare_data(self, flavor: str = "seurat_v3") -> None:
        """数据预处理：计算高变基因和PCA表示"""
        if self.use_rep == "X_pca":
            if "X_pca" not in self.normal_adata.obsm:
                # 提取高变基因
                if "highly_variable" not in self.normal_adata.var.columns:
                    self.gene_list = extract_highly_variable_genes(
                        self.normal_adata, n_top_genes=self.n_top_genes, flavor=flavor
                    )
                else:
                    self.gene_list = self.normal_adata.var_names[
                        self.normal_adata.var["highly_variable"]
                    ].tolist()
                
                # 计算PCA
                import scanpy as sc
                sc.pp.pca(self.normal_adata, n_comps=50, use_highly_variable=True)
        else:
            if self.use_rep not in self.normal_adata.obsm and self.use_rep != "X":
                raise KeyError(f"Representation '{self.use_rep}' not found in adata.obsm")
    
    def build_ot_cost_map(
        self,
        metric: str = "sqeuclidean",
        reg: float = 0.01,
    ) -> Dict[Tuple[str, str], float]:
        """
        构建正常分化轨迹的OT代价地图
        
        Returns
        -------
        dict
            {(stage_i, stage_j): wasserstein_distance}
        """
        if self.stage_distributions is None:
            self.stage_distributions = build_stage_distributions(
                self.normal_adata,
                stage_key=self.stage_key,
                stage_order=self.stage_order,
                gene_list=self.gene_list,
                use_rep=self.use_rep,
            )
        
        self.cost_map = compute_ot_cost_map(
            self.stage_distributions,
            metric=metric,
            reg=reg,
        )
        
        print("=" * 50)
        print("OT Cost Map (Normal Differentiation Trajectory)")
        print("=" * 50)
        for (s1, s2), dist in self.cost_map.items():
            print(f"{s1:15s} -> {s2:15s}: {dist:.4f}")
        print("=" * 50)
        
        return self.cost_map
    
    def compute_tdi(
        self,
        disease_adata: AnnData,
        metric: str = "sqeuclidean",
        reg: float = 0.01,
    ) -> pd.DataFrame:
        """
        计算疾病细胞的轨迹偏离指数（TDI）
        
        Returns
        -------
        pd.DataFrame
            每个细胞的TDI及最优归属阶段
        """
        if self.stage_distributions is None:
            raise RuntimeError("Please run build_ot_cost_map() first")
        
        X_disease = compute_expression_distribution(
            disease_adata, gene_list=self.gene_list, use_rep=self.use_rep
        )
        
        records = []
        for i in range(X_disease.shape[0]):
            X_cell = X_disease[i:i+1, :]
            distances = compute_cell_level_tdi(
                X_cell, self.stage_distributions, metric=metric, reg=reg
            )
            
            tdi = min(distances.values())
            best_stage = min(distances, key=distances.get)
            
            record = {
                "cell_idx": i,
                "TDI": tdi,
                "Best_Match_Stage": best_stage,
            }
            record.update({f"Dist_{stage}": d for stage, d in distances.items()})
            records.append(record)
        
        self.tdi_results = pd.DataFrame(records)
        
        # 统计摘要
        print(f"\nTDI Summary (n={len(self.tdi_results)}):")
        print(f"  Mean TDI: {self.tdi_results['TDI'].mean():.4f}")
        print(f"  Median TDI: {self.tdi_results['TDI'].median():.4f}")
        print(f"  Stage distribution:")
        stage_counts = self.tdi_results["Best_Match_Stage"].value_counts()
        for stage, count in stage_counts.items():
            print(f"    {stage}: {count} ({count/len(self.tdi_results)*100:.1f}%)")
        
        return self.tdi_results
    
    def infer_repair_pathway(
        self,
        disease_adata: AnnData,
        metric: str = "sqeuclidean",
        reg: float = 0.01,
        top_n: int = 100,
    ) -> Dict:
        """
        推断从疾病状态到正常终末状态的修复路径
        
        Returns
        -------
        dict
            包含距离、传输计划、修复靶点等
        """
        terminal_adata = get_stage_cells(
            self.normal_adata, self.stage_key, self.terminal_stage
        )
        
        self.repair_results = infer_repair_pathway(
            adata_disease=disease_adata,
            adata_terminal=terminal_adata,
            gene_list=self.gene_list,
            use_rep=self.use_rep,
            metric=metric,
            reg=reg,
            top_n=top_n,
        )
        
        print(f"\nRepair Pathway (Disease -> {self.terminal_stage}):")
        print(f"  Wasserstein Distance: {self.repair_results['wasserstein_distance']:.4f}")
        print(f"  Top Up-regulated Targets: {len(self.repair_results['top_up_targets'])}")
        print(f"  Top Down-regulated Targets: {len(self.repair_results['top_down_targets'])}")
        
        return self.repair_results
    
    def module_repair_strategy(
        self,
        module_dict: Dict[str, List[str]],
        threshold: float = 0.0,
        method: str = "mean",
    ) -> pd.DataFrame:
        """
        计算模块级修复策略
        
        Parameters
        ----------
        module_dict : dict
            {module_name: [gene1, gene2, ...]}
        threshold : float
            判定显著性的阈值
        method : str
            聚合方法: "mean", "sum", "median"
        
        Returns
        -------
        pd.DataFrame
            模块修复策略表
        """
        if self.repair_results is None or not self.repair_results.get("repair_deltas"):
            raise RuntimeError("Please run infer_repair_pathway() first with gene-level representation")
        
        self.module_strategy = module_repair_strategy(
            repair_deltas=self.repair_results["repair_deltas"],
            module_dict=module_dict,
            threshold=threshold,
            method=method,
        )
        
        print("\nModule Repair Strategy:")
        print(self.module_strategy.to_string(index=False))
        
        return self.module_strategy
    
    def run_full_pipeline(
        self,
        disease_adata: AnnData,
        module_dict: Dict[str, List[str]] = None,
        metric: str = "sqeuclidean",
        reg: float = 0.01,
        top_n: int = 100,
        output_dir: str = None,
    ) -> Dict:
        """
        运行完整分析流程
        
        Returns
        -------
        dict
            所有结果汇总
        """
        self.prepare_data()
        self.build_ot_cost_map(metric=metric, reg=reg)
        self.compute_tdi(disease_adata, metric=metric, reg=reg)
        self.infer_repair_pathway(disease_adata, metric=metric, reg=reg, top_n=top_n)
        
        if module_dict is not None:
            self.module_repair_strategy(module_dict)
        
        results = {
            "cost_map": self.cost_map,
            "tdi": self.tdi_results,
            "repair": self.repair_results,
            "module_strategy": self.module_strategy,
        }
        
        if output_dir is not None:
            save_results(results, output_dir)
        
        return results
    
    # ==================== 可视化接口 ====================
    
    def plot_ot_cost_map(self, **kwargs):
        """绘制OT代价地图"""
        if self.cost_map is None:
            raise RuntimeError("Please run build_ot_cost_map() first")
        return plotting.plot_ot_cost_map(self.cost_map, **kwargs)
    
    def plot_tdi_distribution(self, **kwargs):
        """绘制TDI分布"""
        if self.tdi_results is None:
            raise RuntimeError("Please run compute_tdi() first")
        return plotting.plot_tdi_distribution(self.tdi_results, **kwargs)
    
    def plot_repair_heatmap(self, top_n: int = 50, **kwargs):
        """绘制修复靶点热图"""
        if self.repair_results is None:
            raise RuntimeError("Please run infer_repair_pathway() first")
        return plotting.plot_repair_heatmap(self.repair_results, top_n=top_n, **kwargs)
    
    def plot_module_strategy(self, **kwargs):
        """绘制模块修复策略图"""
        if self.module_strategy is None:
            raise RuntimeError("Please run module_repair_strategy() first")
        return plotting.plot_module_strategy(self.module_strategy, **kwargs)

    # ==================== 元细胞接口 ====================

    def build_metacells(
        self,
        adata: AnnData,
        resolution: float = 1.0,
        n_neighbors: int = 15,
        label: Optional[str] = None,
    ) -> AnnData:
        """
        将细胞聚合为元细胞（Metacells）
        
        Parameters
        ----------
        adata : AnnData
            输入细胞，需包含 self.use_rep 对应的表示层
        resolution : float
            Leiden 聚类分辨率
        n_neighbors : int
            k-NN 邻居数
        label : str, optional
            元细胞 ID 前缀
        
        Returns
        -------
        AnnData
            元细胞对象
        """
        return metacells.build_metacells(
            adata=adata,
            use_rep=self.use_rep,
            resolution=resolution,
            n_neighbors=n_neighbors,
            label=label,
        )

    def build_stage_metacells(
        self,
        resolution: float = 1.0,
        resolution_scale: float = 1.0,
        n_neighbors: int = 15,
        target_metacells: int = 30,
    ) -> AnnData:
        """
        按分化阶段分别构建元细胞
        
        Returns
        -------
        AnnData
            跨阶段的元细胞对象，.obs['stage'] 为阶段标签
        """
        return metacells.build_stage_metacells(
            adata=self.normal_adata,
            stage_key=self.stage_key,
            stage_order=self.stage_order,
            use_rep=self.use_rep,
            resolution=resolution,
            resolution_scale=resolution_scale,
            n_neighbors=n_neighbors,
            target_metacells=target_metacells,
        )
