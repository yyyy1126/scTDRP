"""
scTDRP 红系脱核示例

使用正常骨髓红系数据作为参考，AML5恶性红系细胞作为疾病样本，
演示完整的 scTDRP 分析流程。
"""

import scanpy as sc
import pandas as pd
from scTDRP import TDRPAnalyzer

# ==================== 1. 数据加载 ====================
# 正常红系参考数据（假设为 BoneMarrowMap 中提取的红系细胞）
# 需包含 obs['cell_type'] 列，值为: MEP, BFU-E, CFU-E, Pro-E, Baso-E, Poly-E, Ortho-E
adata_normal = sc.read_h5ad("path/to/normal_erythroid.h5ad")

# AML5 恶性红系细胞
adata_aml5 = sc.read_h5ad("path/to/aml5_malignant_erythroid.h5ad")

# ==================== 2. 初始化分析器 ====================
analyzer = TDRPAnalyzer(
    normal_adata=adata_normal,
    stage_key="cell_type",
    terminal_stage="Ortho-E",
    stage_order=["MEP", "BFU-E", "CFU-E", "Pro-E", "Baso-E", "Poly-E", "Ortho-E"],
    use_rep="X_pca",
    n_top_genes=2000,
)

# ==================== 3. 运行完整流程 ====================
# 先验模块定义（对应原文的P7准备模块和P8执行模块）
module_dict = {
    "P7_TerminalPrep": [
        "ALAS2", "BPGM", "HBA1", "HBA2", "HBB", "BNIP3L", "NCOA4",
        "SLC4A1", "EPB42", "STOM", "TSPO2", "FECH", "BSG"
    ],
    "P8_ExecutionPrecursor": [
        "CDK1", "CCNA2", "CCNB1", "CCNB2", "PLK1", "KIF11",
        "KIF2C", "AURKB", "BUB1", "BUB3", "CDC20", "DLGAP5"
    ],
}

results = analyzer.run_full_pipeline(
    disease_adata=adata_aml5,
    module_dict=module_dict,
    metric="sqeuclidean",
    reg=0.01,
    top_n=100,
    output_dir="./scTDRP_results",
)

# ==================== 4. 可视化 ====================
analyzer.plot_ot_cost_map(save_path="./figures/01_ot_cost_map.pdf")
analyzer.plot_tdi_distribution(save_path="./figures/02_tdi_distribution.pdf")
analyzer.plot_repair_heatmap(top_n=50, save_path="./figures/03_repair_heatmap.pdf")
analyzer.plot_module_strategy(save_path="./figures/04_module_strategy.pdf")

# ==================== 5. 结果解读 ====================
print("\n" + "="*60)
print("结果解读")
print("="*60)

# TDI统计
tdi_mean = results["tdi"]["TDI"].mean()
tdi_median = results["tdi"]["TDI"].median()
print(f"\n1. AML5 恶性细胞的平均 TDI: {tdi_mean:.4f}")
print(f"   中位数 TDI: {tdi_median:.4f}")
print(f"   TDI 越高，说明偏离正常终末状态越远。")

# 最优归属阶段
stage_dist = results["tdi"]["Best_Match_Stage"].value_counts(normalize=True)
print(f"\n2. AML5 恶性细胞的最优归属阶段分布:")
for stage, prop in stage_dist.items():
    print(f"   {stage}: {prop*100:.1f}%")

# 修复靶点
repair = results["repair"]
print(f"\n3. 修复路径 (Disease -> Ortho-E):")
print(f"   Wasserstein 距离: {repair['wasserstein_distance']:.4f}")
print(f"   需上调的 Top 5 靶点: {repair['top_up_targets'][:5]}")
print(f"   需下调的 Top 5 靶点: {repair['top_down_targets'][:5]}")

# 模块策略
if results["module_strategy"] is not None:
    print(f"\n4. 模块级修复策略:")
    for _, row in results["module_strategy"].iterrows():
        print(f"   {row['Module']}: {row['Strategy']} (Score={row['Repair_Score']:.4f})")

print("="*60)
