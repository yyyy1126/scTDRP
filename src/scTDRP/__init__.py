"""
scTDRP: single-cell Terminal Differentiation Repair Pathway

基于最优传输的单细胞终末分化轨迹偏离量化与修复路径推断方法。
"""

from .core import TDRPAnalyzer
from .distance import compute_wasserstein_distance, build_stage_distributions
from .repair import infer_repair_pathway, compute_gene_repair_delta
from .module import module_repair_strategy, aggregate_repair_to_module
from .utils import validate_anndata, extract_highly_variable_genes
from .metacells import build_metacells, build_stage_metacells, map_metacell_scores_to_cells

__version__ = "0.1.0"
__all__ = [
    "TDRPAnalyzer",
    "compute_wasserstein_distance",
    "build_stage_distributions",
    "infer_repair_pathway",
    "compute_gene_repair_delta",
    "module_repair_strategy",
    "aggregate_repair_to_module",
    "validate_anndata",
    "extract_highly_variable_genes",
    "build_metacells",
    "build_stage_metacells",
    "map_metacell_scores_to_cells",
]
