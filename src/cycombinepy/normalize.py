"""Batch-wise normalization of expression data prior to clustering.

Port of ``normalize`` / ``quantile_norm`` / ``clr_norm*`` from
``R/02_batch_correct.R``.
"""

from __future__ import annotations

from typing import Iterable, Literal

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.interpolate import PchipInterpolator
from scipy.stats import rankdata

from cycombinepy._utils import (
    check_obs_key,
    marker_matrix,
    resolve_markers,
    set_marker_matrix,
)

NormMethod = Literal["scale", "rank", "CLR", "CLR_seu", "CLR_med", "qnorm", "none"]
TiesMethod = Literal["average", "min", "max", "dense", "ordinal"]


def _clr(x: np.ndarray) -> np.ndarray:
    """Seurat-flavor CLR: geom_mean = expm1(sum(log1p(x[x>0])) / n)."""
    pos = x[x > 0]
    if pos.size == 0:
        return x
    geom_mean = np.expm1(np.sum(np.log1p(pos)) / x.size)
    if geom_mean == 0:
        return x
    return np.log1p(x / geom_mean)


def _clr_mean(x: np.ndarray) -> np.ndarray:
    """CLR with geom_mean = expm1(mean(log1p(x[x>=0])))."""
    nn = x[x >= 0]
    if nn.size == 0:
        return x
    geom_mean = np.expm1(np.mean(np.log1p(nn)))
    if geom_mean == 0:
        return x
    return np.log((x + 1.0) / geom_mean)


def _clr_med(x: np.ndarray) -> np.ndarray:
    """CLR using the median as the center."""
    m = np.nanmedian(x)
    if m == 0:
        m = 1.0
    return np.log((x + 1.0) / m)


def _rank_norm(x: np.ndarray, ties_method: TiesMethod = "average") -> np.ndarray:
    """Percentile rank within a 1-D vector."""
    return rankdata(x, method=ties_method) / x.size


def _scale(x: np.ndarray) -> np.ndarray:
    """Z-score with an ``ddof=1`` std to match R's ``scale``."""
    mu = np.mean(x)
    sd = np.std(x, ddof=1)
    if sd == 0 or not np.isfinite(sd):
        return x - mu
    return (x - mu) / sd


def _apply_column_wise(block: np.ndarray, fn) -> np.ndarray:
    out = np.empty_like(block, dtype=float)
    for j in range(block.shape[1]):
        col = block[:, j]
        if np.sum(col) == 0:
            out[:, j] = col
        else:
            out[:, j] = fn(col)
    return out


def _quantile_norm(
    X: np.ndarray,
    batches: np.ndarray,
    n_quantiles: int = 5,
) -> np.ndarray:
    """Monotone-spline quantile normalization per marker across batches.

    Mirrors ``quantile_norm`` in ``R/02_batch_correct.R:126``. For each marker we
    compute a reference set of quantiles across all cells, then map each batch's
    quantiles onto the reference using a PCHIP (monotone cubic) spline.
    """
    out = X.astype(float, copy=True)
    q_levels = np.linspace(0.0, 1.0, n_quantiles)
    for j in range(X.shape[1]):
        col = X[:, j]
        ref_q = np.quantile(col, q_levels)
        for b in np.unique(batches):
            mask = batches == b
            batch_q = np.quantile(col[mask], q_levels)
            # PCHIP requires strictly increasing x. Fall back to identity if degenerate.
            if np.any(np.diff(batch_q) <= 0):
                continue
            spline = PchipInterpolator(batch_q, ref_q, extrapolate=True)
            out[mask, j] = spline(col[mask])
    return out


def normalize(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    method: NormMethod = "scale",
    batch_key: str = "batch",
    ties_method: TiesMethod = "average",
    layer: str | None = None,
    copy: bool = False,
) -> AnnData | None:
    """Batch-wise normalize marker columns of ``adata``.

    Port of ``normalize`` in ``R/02_batch_correct.R:27-111``. Each batch is
    processed independently, so that downstream clustering is less influenced by
    between-batch shifts.

    Parameters
    ----------
    adata
        AnnData containing expression in ``adata.X`` (or a layer).
    markers
        Var names to normalize. If ``None``, :func:`cycombinepy.get_markers` is used.
    method
        One of ``"scale"``, ``"rank"``, ``"CLR"``, ``"CLR_seu"``, ``"CLR_med"``,
        ``"qnorm"``, ``"none"``.
    batch_key
        Column in ``adata.obs`` identifying the batch.
    ties_method
        Tie-breaking rule for ``method="rank"``.
    layer
        If given, read / write that layer instead of ``adata.X``.
    copy
        If True, return a copy; otherwise modify in place and return None.
    """
    if method == "none":
        return adata.copy() if copy else None

    check_obs_key(adata, batch_key)
    markers = resolve_markers(adata, markers)
    if copy:
        adata = adata.copy()

    X = marker_matrix(adata, markers, layer=layer)
    batches = np.asarray(adata.obs[batch_key].values)

    if method == "qnorm":
        new_X = _quantile_norm(X, batches)
    else:
        if method == "scale":
            fn = _scale
        elif method == "rank":
            fn = lambda col: _rank_norm(col, ties_method)
        elif method == "CLR":
            fn = _clr_mean
        elif method == "CLR_seu":
            fn = _clr
        elif method == "CLR_med":
            fn = _clr_med
        else:
            raise ValueError(f"Unknown normalization method: {method!r}")

        new_X = np.empty_like(X, dtype=float)
        # Preserve row order: apply per-batch, scatter results back.
        for b in pd.unique(batches):
            mask = batches == b
            new_X[mask] = _apply_column_wise(X[mask], fn)

    set_marker_matrix(adata, markers, new_X, layer=layer)
    return adata if copy else None
