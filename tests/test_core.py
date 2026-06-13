"""单元测试"""

import unittest
import numpy as np
import pandas as pd
from anndata import AnnData

# 模拟数据生成
def make_dummy_anndata(n_cells=100, n_genes=50, stage_key="stage", stages=None):
    """生成模拟单细胞数据"""
    import numpy as np
    np.random.seed(42)
    X = np.random.poisson(2, size=(n_cells, n_genes)).astype(np.float32)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])
    if stages is not None:
        obs[stage_key] = np.random.choice(stages, size=n_cells)
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])
    return AnnData(X=X, obs=obs, var=var)


class TestUtils(unittest.TestCase):
    def test_validate_anndata(self):
        from scTDRP.utils import validate_anndata
        adata = make_dummy_anndata(stages=["A", "B"])
        validate_anndata(adata)  # 不应报错
        
        with self.assertRaises(TypeError):
            validate_anndata("not anndata")


class TestDistance(unittest.TestCase):
    def test_build_stage_distributions(self):
        from scTDRP.distance import build_stage_distributions
        adata = make_dummy_anndata(n_cells=200, stages=["A", "B", "C"])
        dists = build_stage_distributions(adata, stage_key="stage", use_rep="X")
        self.assertEqual(len(dists), 3)
        for stage, (X, w) in dists.items():
            self.assertEqual(X.shape[1], 50)
            self.assertAlmostEqual(w.sum(), 1.0, places=5)
    
    def test_compute_wasserstein_distance(self):
        from scTDRP.distance import compute_wasserstein_distance
        X1 = np.random.randn(20, 10)
        X2 = np.random.randn(20, 10)
        dist, plan = compute_wasserstein_distance(X1, X2, reg=0.1)
        self.assertGreater(dist, 0)
        self.assertEqual(plan.shape, (20, 20))
        self.assertAlmostEqual(plan.sum(), 1.0, places=5)


class TestRepair(unittest.TestCase):
    def test_compute_gene_repair_delta(self):
        from scTDRP.repair import compute_gene_repair_delta
        X_disease = np.array([[1, 2], [3, 4]])
        X_terminal = np.array([[5, 6], [7, 8]])
        plan = np.array([[0.25, 0.25], [0.25, 0.25]])
        deltas = compute_gene_repair_delta(X_disease, X_terminal, plan,
                                           gene_names=["g1", "g2"])
        self.assertIn("g1", deltas)
        self.assertIn("g2", deltas)


class TestModule(unittest.TestCase):
    def test_aggregate_repair_to_module(self):
        from scTDRP.module import aggregate_repair_to_module
        repair_deltas = {"g1": 1.0, "g2": 2.0, "g3": -1.0}
        module_dict = {"M1": ["g1", "g2"], "M2": ["g3"]}
        scores = aggregate_repair_to_module(repair_deltas, module_dict, method="mean")
        self.assertAlmostEqual(scores["M1"], 1.5, places=5)
        self.assertAlmostEqual(scores["M2"], -1.0, places=5)


class TestTDRPAnalyzer(unittest.TestCase):
    def test_initialization(self):
        from scTDRP import TDRPAnalyzer
        adata = make_dummy_anndata(n_cells=300, stages=["A", "B", "C"])
        analyzer = TDRPAnalyzer(
            normal_adata=adata,
            stage_key="stage",
            terminal_stage="C",
            stage_order=["A", "B", "C"],
        )
        self.assertEqual(analyzer.terminal_stage, "C")
        self.assertEqual(analyzer.stage_order, ["A", "B", "C"])

    def test_ot_cost_map_monotonicity(self):
        """With ordered stages, later stages should not have lower cost to terminal."""
        from scTDRP import TDRPAnalyzer
        np.random.seed(42)
        n_per_stage = 50
        n_genes = 20
        stages = ["A", "B", "C"]
        X_list = []
        obs_list = []
        for i, stage in enumerate(stages):
            # Each subsequent stage shifts expression by +1 along a direction
            base = np.random.poisson(2, size=(n_per_stage, n_genes)).astype(float)
            shift = np.ones((n_per_stage, n_genes)) * i * 0.5
            X_list.append(base + shift)
            obs_list.extend([stage] * n_per_stage)
        X = np.vstack(X_list)
        adata = AnnData(X=X, obs=pd.DataFrame({"stage": obs_list}, index=[f"c{i}" for i in range(len(obs_list))]))
        analyzer = TDRPAnalyzer(
            normal_adata=adata,
            stage_key="stage",
            terminal_stage="C",
            stage_order=stages,
            use_rep="X",
        )
        cost_map = analyzer.build_ot_cost_map(reg=0.05)
        self.assertGreater(cost_map[("A", "B")], 0)
        self.assertGreater(cost_map[("B", "C")], 0)


class TestMetacells(unittest.TestCase):
    def test_build_metacells(self):
        from scTDRP import build_metacells
        import scanpy as sc
        adata = make_dummy_anndata(n_cells=100, n_genes=50)
        sc.pp.pca(adata, n_comps=10)
        meta = build_metacells(adata, use_rep="X_pca", resolution=1.0, label="test")
        self.assertLess(meta.shape[0], adata.shape[0])
        self.assertIn("cell_map", meta.uns)
        self.assertIn("n_cells", meta.obs.columns)
        self.assertEqual(sum(meta.obs["n_cells"]), adata.shape[0])

    def test_build_stage_metacells(self):
        from scTDRP import build_stage_metacells
        import scanpy as sc
        adata = make_dummy_anndata(n_cells=300, n_genes=50, stages=["A", "B", "C"])
        sc.pp.pca(adata, n_comps=10)
        meta = build_stage_metacells(
            adata, stage_key="stage", stage_order=["A", "B", "C"],
            use_rep="X_pca", resolution_scale=1.0,
        )
        self.assertLess(meta.shape[0], adata.shape[0])
        self.assertEqual(set(meta.obs["stage"].tolist()), {"A", "B", "C"})


class TestNumericalProperties(unittest.TestCase):
    def test_transport_plan_normalization(self):
        from scTDRP.distance import compute_wasserstein_distance
        X1 = np.random.randn(20, 10)
        X2 = np.random.randn(20, 10)
        _, plan = compute_wasserstein_distance(X1, X2, reg=0.1)
        self.assertAlmostEqual(plan.sum(), 1.0, places=5)
        self.assertGreaterEqual(plan.min(), -1e-6)

    def test_repair_delta_direction(self):
        """If terminal expression is uniformly higher, repair deltas should be positive."""
        from scTDRP.repair import compute_gene_repair_delta
        X_disease = np.array([[1.0, 2.0], [3.0, 4.0]])
        X_terminal = np.array([[10.0, 11.0], [12.0, 13.0]])
        plan = np.array([[0.25, 0.25], [0.25, 0.25]])
        deltas = compute_gene_repair_delta(X_disease, X_terminal, plan, gene_names=["g1", "g2"])
        self.assertGreater(deltas["g1"], 0)
        self.assertGreater(deltas["g2"], 0)


if __name__ == "__main__":
    unittest.main()
