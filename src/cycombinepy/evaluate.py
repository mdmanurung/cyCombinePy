"""Batch-effect evaluation: EMD, MAD, and a scib-metrics wrapper.

Ports of ``compute_emd`` / ``evaluate_emd`` / ``compute_mad`` / ``evaluate_mad``
from ``R/evaluate_performance.R``, plus a wrapper over ``scib_metrics`` for
scanpy-native benchmark metrics (kBET, iLISI/cLISI, graph connectivity, ...).
"""

from __future__ import annotations

from itertools import combinations
from typing import Iterable

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import wasserstein_distance

from cycombinepy._utils import (
    check_obs_values_not_missing,
    marker_matrix,
    resolve_markers,
)


def compute_emd(
    adata: AnnData,
    cell_key: str = "cycombine_som",
    batch_key: str = "batch",
    markers: Iterable[str] | None = None,
    layer: str | None = None,
) -> pd.DataFrame:
    """Per (cluster, marker, batch-pair) 1-D Earth Mover's distance.

    Returns a tidy DataFrame with columns ``cluster, marker, batch1, batch2,
    emd``. Uses :func:`scipy.stats.wasserstein_distance`, which is equivalent to
    the 1-D EMD that cyCombine computes via ``emdist::emd2d`` on single-column
    matrices.
    """
    label_values = check_obs_values_not_missing(
        adata,
        cell_key,
        context="compute_emd()",
    )
    batch_values = check_obs_values_not_missing(
        adata,
        batch_key,
        context="compute_emd()",
    )
    markers = resolve_markers(adata, markers)
    X = marker_matrix(
        adata,
        markers,
        layer=layer,
        require_finite=True,
        context="compute_emd()",
    )
    labels = label_values.astype(str).to_numpy()
    batches = batch_values.astype(str).to_numpy()

    # Two passes: first count, then fill preallocated arrays. This avoids the
    # list-of-dicts intermediate and the per-row DataFrame construction cost.
    n_markers = len(markers)
    uniq_labels = np.unique(labels)
    plan: list[tuple] = []  # (lab, b1, b2, A_idx, B_idx)
    for lab in uniq_labels:
        mask_l = labels == lab
        present = sorted(np.unique(batches[mask_l]).tolist())
        for b1, b2 in combinations(present, 2):
            A_idx = np.flatnonzero(mask_l & (batches == b1))
            B_idx = np.flatnonzero(mask_l & (batches == b2))
            if A_idx.size == 0 or B_idx.size == 0:
                continue
            plan.append((lab, b1, b2, A_idx, B_idx))

    n_rows = len(plan) * n_markers
    cluster_col = np.empty(n_rows, dtype=object)
    marker_col = np.empty(n_rows, dtype=object)
    batch1_col = np.empty(n_rows, dtype=object)
    batch2_col = np.empty(n_rows, dtype=object)
    emd_col = np.empty(n_rows, dtype=np.float64)

    k = 0
    for lab, b1, b2, A_idx, B_idx in plan:
        A = X[A_idx]
        B = X[B_idx]
        for j, marker in enumerate(markers):
            cluster_col[k] = lab
            marker_col[k] = marker
            batch1_col[k] = b1
            batch2_col[k] = b2
            emd_col[k] = wasserstein_distance(A[:, j], B[:, j])
            k += 1

    return pd.DataFrame(
        {
            "cluster": cluster_col,
            "marker": marker_col,
            "batch1": batch1_col,
            "batch2": batch2_col,
            "emd": emd_col,
        }
    )


def evaluate_emd(
    uncorrected: pd.DataFrame,
    corrected: pd.DataFrame,
) -> pd.DataFrame:
    """Join uncorrected vs corrected EMD and compute percent reduction."""
    keys = ["cluster", "marker", "batch1", "batch2"]
    merged = uncorrected.merge(
        corrected, on=keys, suffixes=("_uncorrected", "_corrected")
    )
    merged["reduction"] = merged["emd_uncorrected"] - merged["emd_corrected"]
    denom = merged["emd_uncorrected"].replace(0, np.nan)
    merged["reduction_pct"] = 100.0 * merged["reduction"] / denom
    return merged


def compute_mad(
    adata: AnnData,
    cell_key: str = "cycombine_som",
    batch_key: str = "batch",
    markers: Iterable[str] | None = None,
    layer: str | None = None,
) -> pd.DataFrame:
    """Per (cluster, marker, batch) Median Absolute Deviation.

    Returns a tidy DataFrame with columns ``cluster, marker, batch, mad``.
    Mirrors ``compute_mad`` in ``R/evaluate_performance.R``: MAD is the median
    of ``|x - median(x)|`` within each (cluster, batch) block.
    """
    label_values = check_obs_values_not_missing(
        adata,
        cell_key,
        context="compute_mad()",
    )
    batch_values = check_obs_values_not_missing(
        adata,
        batch_key,
        context="compute_mad()",
    )
    markers = resolve_markers(adata, markers)
    X = marker_matrix(
        adata,
        markers,
        layer=layer,
        require_finite=True,
        context="compute_mad()",
    )
    labels = label_values.astype(str).to_numpy()
    batches = batch_values.astype(str).to_numpy()

    n_markers = len(markers)
    uniq_labels = np.unique(labels)

    plan: list[tuple] = []  # (lab, b, block_idx)
    for lab in uniq_labels:
        mask_l = labels == lab
        for b in np.unique(batches[mask_l]):
            block_idx = np.flatnonzero(mask_l & (batches == b))
            if block_idx.size == 0:
                continue
            plan.append((lab, b, block_idx))

    n_rows = len(plan) * n_markers
    cluster_col = np.empty(n_rows, dtype=object)
    marker_col = np.empty(n_rows, dtype=object)
    batch_col = np.empty(n_rows, dtype=object)
    mad_col = np.empty(n_rows, dtype=np.float64)

    k = 0
    for lab, b, block_idx in plan:
        block = X[block_idx]
        med = np.median(block, axis=0)
        mad = np.median(np.abs(block - med), axis=0)
        for j, marker in enumerate(markers):
            cluster_col[k] = lab
            marker_col[k] = marker
            batch_col[k] = b
            mad_col[k] = mad[j]
            k += 1

    return pd.DataFrame(
        {
            "cluster": cluster_col,
            "marker": marker_col,
            "batch": batch_col,
            "mad": mad_col,
        }
    )


def evaluate_mad(
    uncorrected: pd.DataFrame,
    corrected: pd.DataFrame,
) -> pd.DataFrame:
    """Join uncorrected vs corrected MAD and compute percent reduction."""
    keys = ["cluster", "marker", "batch"]
    merged = uncorrected.merge(
        corrected, on=keys, suffixes=("_uncorrected", "_corrected")
    )
    merged["reduction"] = merged["mad_uncorrected"] - merged["mad_corrected"]
    denom = merged["mad_uncorrected"].replace(0, np.nan)
    merged["reduction_pct"] = 100.0 * merged["reduction"] / denom
    return merged


def scib_metrics(
    adata: AnnData,
    batch_key: str,
    label_key: str | None = None,
    embedding_key: str = "X_pca",
    layer: str | None = None,
) -> dict:
    """Run selected scib-metrics batch-correction metrics on an AnnData.

    Computes a PCA on ``adata.X`` (or ``adata.layers[layer]`` if supplied) and
    evaluates batch-mixing metrics from :mod:`scib_metrics`. To compare
    uncorrected and corrected data, call it once for each representation and
    compare the returned dictionaries.

    Returns a dict of scalar scores. Metrics that require a biological label are
    skipped if ``label_key`` is ``None``.
    """
    batch_values = check_obs_values_not_missing(
        adata,
        batch_key,
        context="scib_metrics()",
    )
    if label_key is not None:
        label_values = check_obs_values_not_missing(
            adata,
            label_key,
            context="scib_metrics()",
        )
    else:
        label_values = None
    markers = list(adata.var_names)
    X = marker_matrix(
        adata,
        markers,
        layer=layer,
        require_finite=True,
        context="scib_metrics()",
    )

    try:
        from scib_metrics import (  # type: ignore
            graph_connectivity,
            ilisi_knn,
            silhouette_batch,
        )
        from scib_metrics.nearest_neighbors import pynndescent  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "scib_metrics is required for cycombinepy.scib_metrics; install with "
            "`pip install scib-metrics`."
        ) from exc

    import scanpy as sc

    a = adata.copy()
    a.X = X

    sc.pp.pca(a, n_comps=min(20, a.n_vars - 1))
    X_emb = a.obsm[embedding_key]
    batches = batch_values.to_numpy()

    knn = pynndescent(X_emb, n_neighbors=30, random_state=0)
    scores: dict = {}
    scores["graph_connectivity"] = float(graph_connectivity(knn, batches))
    scores["ilisi"] = float(ilisi_knn(knn, batches))
    if label_key is not None:
        labels = label_values.to_numpy()
        scores["silhouette_batch"] = float(
            silhouette_batch(X_emb, labels, batches)
        )
    return scores
