"""
Metacell construction utilities for scTDRP.

Metacells aggregate biologically similar cells into a smaller number of
representative profiles, reducing the computational cost of optimal transport
while preserving population structure.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from .utils import validate_anndata


def _resolve_representation(adata: AnnData, use_rep: str) -> str:
    """Ensure the requested representation exists, falling back to X_pca."""
    if use_rep in adata.obsm:
        return use_rep
    if "X_pca" in adata.obsm:
        return "X_pca"
    raise KeyError(
        f"Representation '{use_rep}' not found in adata.obsm and 'X_pca' is also unavailable. "
        "Please compute PCA or provide another embedding."
    )


def build_metacells(
    adata: AnnData,
    use_rep: str = "X_pca",
    resolution: float = 1.0,
    n_neighbors: int = 15,
    label: Optional[str] = None,
    copy: bool = True,
) -> AnnData:
    """
    Aggregate cells into metacells via Leiden clustering.

    Parameters
    ----------
    adata : AnnData
        Input cells. Must contain `use_rep` in `.obsm`.
    use_rep : str
        Embedding used for neighborhood graph and clustering. Defaults to 'X_pca'.
    resolution : float
        Leiden resolution parameter. Higher values yield more metacells.
    n_neighbors : int
        Number of neighbors for the k-NN graph.
    label : str, optional
        Prefix for metacell IDs.
    copy : bool
        If True, work on a copy of `adata`.

    Returns
    -------
    AnnData
        Metacell object where each observation is one metacell. `.obs` contains
        'n_cells' (number of aggregated cells) and 'metacell_id'. `.uns['cell_map']`
        maps each original cell to its metacell ID.
    """
    validate_anndata(adata)
    rep_key = _resolve_representation(adata, use_rep)

    ad = adata.copy() if copy else adata
    sc.pp.neighbors(ad, use_rep=rep_key, n_neighbors=n_neighbors)
    sc.tl.leiden(ad, resolution=resolution)

    clusters = ad.obs["leiden"].unique()
    n_meta = len(clusters)
    prefix = label if label is not None else "MC"

    X_list = []
    obs_list = []
    meta_rep = [] if rep_key in ad.obsm else None

    for cl in clusters:
        mask = ad.obs["leiden"] == cl
        cells = ad[mask]
        X_cl = cells.X.mean(axis=0)
        if hasattr(X_cl, "toarray"):
            X_cl = X_cl.toarray().flatten()
        X_list.append(np.asarray(X_cl).flatten())
        obs_list.append({
            "metacell_id": f"{prefix}_{cl}",
            "n_cells": int(mask.sum()),
            "leiden": cl,
        })
        if meta_rep is not None:
            meta_rep.append(ad.obsm[rep_key][mask].mean(axis=0))

    meta_X = np.vstack(X_list)
    meta_obs = pd.DataFrame(obs_list)
    meta_var = ad.var.copy()
    meta_ad = AnnData(X=meta_X, obs=meta_obs, var=meta_var)
    meta_ad.obs_names = meta_obs["metacell_id"].values.astype(str)

    if meta_rep is not None:
        meta_ad.obsm[rep_key] = np.vstack(meta_rep)

    # Store cell-to-metacell mapping for downstream deconvolution
    meta_ad.uns["cell_map"] = {
        cell: f"{prefix}_{cl}" for cell, cl in ad.obs["leiden"].items()
    }

    return meta_ad


def build_stage_metacells(
    adata: AnnData,
    stage_key: str,
    stage_order: Optional[List[str]] = None,
    use_rep: str = "X_pca",
    resolution: float = 1.0,
    resolution_scale: float = 1.0,
    n_neighbors: int = 15,
    target_metacells: int = 30,
    min_cells: int = 10,
    copy: bool = True,
) -> AnnData:
    """
    Build metacells separately for each developmental stage.

    This preserves stage identity while reducing within-stage cell counts,
    which is useful when constructing the normal reference distribution for
    optimal transport.

    Parameters
    ----------
    adata : AnnData
        Input cells with `stage_key` annotation.
    stage_key : str
        Column in `.obs` identifying developmental stages.
    stage_order : list, optional
        Ordered list of stage names to process. If None, uses all stages in
        `adata.obs[stage_key]`.
    use_rep : str
        Embedding used for clustering.
    resolution : float
        Base Leiden resolution. Ignored if `resolution_scale` is used.
    resolution_scale : float
        If > 0, automatically tune resolution per stage to yield approximately
        `target_metacells` metacells given the stage size.
    n_neighbors : int
        Number of neighbors for the k-NN graph.
    target_metacells : int
        Desired number of metacells per stage when `resolution_scale` is used.
    min_cells : int
        Minimum number of cells required to build metacells for a stage.
    copy : bool
        If True, work on a copy of `adata`.

    Returns
    -------
    AnnData
        Concatenated metacells across all stages. `.obs['stage']` contains the
        original stage label.
    """
    validate_anndata(adata, required_obs_keys=[stage_key])
    ad = adata.copy() if copy else adata

    if stage_order is None:
        stage_order = ad.obs[stage_key].unique().tolist()

    meta_list = []
    for stage in stage_order:
        if stage not in ad.obs[stage_key].values:
            continue
        stage_ad = ad[ad.obs[stage_key] == stage].copy()
        if stage_ad.shape[0] < min_cells:
            continue

        if resolution_scale > 0:
            # Heuristic: scale resolution inversely with stage size
            res = max(0.5, min(3.0, target_metacells * resolution_scale / max(1, stage_ad.shape[0] / 50)))
        else:
            res = resolution

        meta_stage = build_metacells(
            stage_ad,
            use_rep=use_rep,
            resolution=res,
            n_neighbors=n_neighbors,
            label=str(stage),
            copy=True,
        )
        meta_stage.obs[stage_key] = stage
        meta_list.append(meta_stage)

    if len(meta_list) == 0:
        raise ValueError("No metacells could be built. Check stage_key and min_cells.")

    # Determine categories for ordered stage label
    available_stages = [s for s in stage_order if any(meta.obs[stage_key].iloc[0] == s for meta in meta_list)]

    import anndata as ad_module
    meta_combined = ad_module.concat(meta_list, label="_stage_batch", index_unique="-")
    meta_combined.obs[stage_key] = pd.Categorical(
        [s for meta in meta_list for s in [meta.obs[stage_key].iloc[0]] * meta.shape[0]],
        categories=available_stages,
        ordered=True,
    )
    return meta_combined


def map_metacell_scores_to_cells(
    metacell_scores: Dict[str, np.ndarray],
    metacell_adata: AnnData,
    cell_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Map per-metacell scores back to individual cells using the cell_map stored
    in `.uns['cell_map']`.

    Parameters
    ----------
    metacell_scores : dict
        {score_name: array of length n_metacells}
    metacell_adata : AnnData
        Metacell object containing 'cell_map' in `.uns`.
    cell_names : list, optional
        If provided, return scores only for these cells in this order.

    Returns
    -------
    pd.DataFrame
        Cell-level scores with one column per score name.
    """
    if "cell_map" not in metacell_adata.uns:
        raise KeyError("Metacell object does not contain 'cell_map' in .uns")

    cell_map = metacell_adata.uns["cell_map"]
    if cell_names is None:
        cell_names = list(cell_map.keys())

    meta_to_score = {}
    for name, scores in metacell_scores.items():
        mapping = dict(zip(metacell_adata.obs_names, np.asarray(scores).flatten()))
        meta_to_score[name] = mapping

    records = []
    for cell in cell_names:
        meta_id = cell_map.get(cell)
        record = {"cell": cell}
        for name in metacell_scores:
            record[name] = meta_to_score[name].get(meta_id, np.nan)
        records.append(record)

    return pd.DataFrame(records).set_index("cell")
