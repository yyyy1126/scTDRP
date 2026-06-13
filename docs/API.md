# scTDRP API 文档

## TDRPAnalyzer

主分析类，整合完整分析流程。

### 初始化

```python
TDRPAnalyzer(
    normal_adata,      # AnnData: 正常参考数据
    stage_key,         # str: 阶段注释列名
    terminal_stage,    # str: 终末阶段名称
    stage_order=None,  # list: 阶段顺序
    use_rep="X_pca",   # str: 表示层
    n_top_genes=2000,  # int: 高变基因数
)
```

### 核心方法

#### prepare_data
```python
prepare_data(flavor="seurat_v3")
```
数据预处理，计算高变基因和PCA。

#### build_ot_cost_map
```python
build_ot_cost_map(metric="sqeuclidean", reg=0.01)
```
构建正常分化轨迹的OT代价地图。

**返回**: `dict` — {(stage_i, stage_j): distance}

#### compute_tdi
```python
compute_tdi(disease_adata, metric="sqeuclidean", reg=0.01)
```
计算疾病细胞的轨迹偏离指数（TDI）。

**返回**: `pd.DataFrame` — 每个细胞的TDI及最优归属阶段

#### infer_repair_pathway
```python
infer_repair_pathway(disease_adata, metric="sqeuclidean", reg=0.01, top_n=100)
```
推断从疾病状态到正常终末状态的修复路径。

**返回**: `dict` — 包含距离、传输计划、修复靶点

#### module_repair_strategy
```python
module_repair_strategy(module_dict, threshold=0.0, method="mean")
```
计算模块级修复策略。

**参数**:
- `module_dict`: {module_name: [gene1, gene2, ...]}
- `threshold`: 判定显著性的阈值
- `method`: "mean", "sum", "median"

**返回**: `pd.DataFrame`

#### run_full_pipeline
```python
run_full_pipeline(
    disease_adata,
    module_dict=None,
    metric="sqeuclidean",
    reg=0.01,
    top_n=100,
    output_dir=None,
)
```
运行完整分析流程。

---

## 底层函数

### distance.compute_wasserstein_distance
```python
compute_wasserstein_distance(
    X_source, X_target,
    w_source=None, w_target=None,
    metric="sqeuclidean", reg=0.01, method="sinkhorn"
)
```
计算Wasserstein距离及最优传输计划。

### repair.compute_gene_repair_delta
```python
compute_gene_repair_delta(X_disease, X_terminal, transport_plan, gene_names=None)
```
基于传输计划计算每个基因的修复量。

### module.aggregate_repair_to_module
```python
aggregate_repair_to_module(repair_deltas, module_dict, method="mean")
```
将基因级修复量聚合到模块级。
