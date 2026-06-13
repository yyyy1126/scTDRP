# scTDRP 执行任务清单

按顺序执行，每一步都有命令，复制粘贴即可。

---

## 任务1：安装环境（10分钟）

```bash
cd /data1/yja/zhongzhuan/5.external/scTDRP

# 新建conda环境
conda create -n scTDRP python=3.10 -y
conda activate scTDRP

# 安装依赖
pip install numpy scipy pandas scanpy anndata POT scikit-learn matplotlib seaborn
```

**验证安装成功**：
```bash
python -c "import scanpy, anndata, ot; print('OK')"
```

---

## 任务2：确认数据列名（5分钟）

```bash
conda activate scTDRP
cd /data1/yja/zhongzhuan/5.external/scTDRP
```

创建并运行探查脚本：

```bash
cat > check_columns.py << 'EOF'
import scanpy as sc

# === 正常红系数据 ===
print("="*60)
print("NORMAL: erythroid_lineage_from_MEP.h5ad")
print("="*60)
ad = sc.read_h5ad("../../1.data/processed/erythroid_lineage_from_MEP.h5ad")
print("Obs columns:", ad.obs.columns.tolist())
print("Shape:", ad.shape)

# 尝试找阶段列
for col in ad.obs.columns:
    uniq = ad.obs[col].unique()
    if len(uniq) < 20 and len(uniq) > 2:
        print(f"\nColumn '{col}' has {len(uniq)} unique values:")
        print(uniq.tolist())

# === AML5 数据 ===
print("\n" + "="*60)
print("DISEASE: aml5_seurat.h5ad")
print("="*60)
ad2 = sc.read_h5ad("../../1.data/raw/scRNAseq/Kuusanmaki2024_AML5/aml5_seurat.h5ad")
print("Obs columns:", ad2.obs.columns.tolist())
print("Shape:", ad2.shape)

for col in ad2.obs.columns:
    uniq = ad2.obs[col].unique()
    if len(uniq) < 20:
        print(f"\nColumn '{col}' has {len(uniq)} unique values:")
        print(uniq.tolist())
EOF

python check_columns.py > check_columns.log 2>&1
cat check_columns.log
```

**把 `check_columns.log` 的输出保存下来**，下一步需要据此修改脚本配置。

---

## 任务3：修改脚本配置（5分钟）

打开 `run_erythroid_scTDRP.py`，根据任务2的输出修改以下配置（在文件第28-50行左右）：

```python
# 正常数据中的阶段列名（根据check_columns.log修改）
STAGE_KEY = "cell_type"

# 阶段顺序（根据实际值修改）
STAGE_ORDER = ["MEP", "BFU-E", "CFU-E", "Pro-E", "Baso-E", "Poly-E", "Ortho-E"]

# 终末阶段名称
TERMINAL_STAGE = "Ortho-E"

# AML5中恶性细胞的筛选列和值（根据check_columns.log修改）
DISEASE_CONDITION_KEY = "condition"
DISEASE_CONDITION_VALUE = "malignant"
```

**常见需要修改的情况**：
- 如果正常数据的阶段列叫 `leiden` 或 `cluster` → 改 `STAGE_KEY`
- 如果阶段名称是 `ProE` 而不是 `Pro-E` → 改 `STAGE_ORDER`
- 如果 AML5 没有 condition 列，直接用全部细胞 → 注释掉筛选逻辑

---

## 任务4：运行分析（30分钟-2小时，取决于机器）

```bash
conda activate scTDRP
cd /data1/yja/zhongzhuan/5.external/scTDRP

# 后台运行（不挂终端）
nohup python run_erythroid_scTDRP.py > run.log 2>&1 &

# 查看进度
tail -f run.log
```

按 `Ctrl+C` 退出 tail，任务继续在后台跑。

---

## 任务5：检查结果（10分钟）

运行完成后：

```bash
# 看结果文件
ls -lh results/
ls -lh figures/

# 看关键数字
cat run.log | grep -A 20 "Results Summary"
```

**验证清单**：
- [ ] `results/scTDRP_tdi.csv` 存在且非空
- [ ] `results/scTDRP_repair_deltas.csv` 存在且非空
- [ ] `results/scTDRP_module_strategy.csv` 中 P7 显示"上调"，P8 显示"下调"
- [ ] `figures/` 下有 4 张 PDF
- [ ] `run.log` 最后显示 "Done!"

---

## 任务6：微调参数（可选，30分钟）

如果结果看起来不合理，尝试调整：

| 问题 | 解决方案 |
|:---|:---|
| TDI 分布太集中 | 减小 `OT_REG`（如 0.005） |
| TDI 分布太分散 | 增大 `OT_REG`（如 0.05） |
| 修复靶点全是管家基因 | 改用 `USE_REP = "X_pca"` 而非 `"X"` |
| 模块策略方向反了 | 检查 P7/P8 基因是否对应正确 |

修改后重新运行任务4。

---

## 任务7：准备论文图表（1小时）

4 张 figures/ 下的 PDF 可以直接放入论文：

| 图编号 | 文件 | 对应论文内容 |
|:---|:---|:---|
| Fig.1 | `01_ot_cost_map.pdf` | 正常分化OT代价地图 |
| Fig.2 | `02_tdi_distribution.pdf` | AML5恶性细胞TDI分布 |
| Fig.3 | `03_repair_heatmap.pdf` | Top修复靶点 |
| Fig.4 | `04_module_strategy.pdf` | P7/P8模块修复策略 |

如果需要调整样式，修改 `src/scTDRP/plotting.py` 中的颜色、字体等参数。

---

## 全部完成后

- [ ] 环境安装成功
- [ ] 数据列名已确认
- [ ] 脚本配置已修改
- [ ] 运行成功，4张图已生成
- [ ] 结果合理（P7需上调、P8需下调、TDI>0）

全部打勾后，这个项目就算跑通了。下一步是写论文方法学和结果部分。
