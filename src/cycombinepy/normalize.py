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
    check_obs_values_not_missing,
    marker_matrix,
    resolve_markers,
    set_marker_matrix,
)

NormMethod = Literal["scale", "rank", "CLR", "CLR_seu", "CLR_med", "qnorm", "none"]
TiesMethod = Literal["average", "min", "max", "dense", "ordinal"]


def _zero_col_mask(block: np.ndarray) -> np.ndarray:
    """Per-column boolean mask flagging columns whose entries all sum to zero.

    Mirrors the ``np.sum(col) == 0`` short-circuit used by the scalar
    per-column helpers.
    """
    return block.sum(axis=0) == 0


def _restore_zero_cols(
    out: np.ndarray, block: np.ndarray, zero_cols: np.ndarray
) -> np.ndarray:
    if zero_cols.any():
        out[:, zero_cols] = block[:, zero_cols]
    return out


def _scale_block(block: np.ndarray) -> np.ndarray:
    """Column-wise z-score (ddof=1) over a (n_cells, n_markers) block."""
    mu = block.mean(axis=0)
    sd = block.std(axis=0, ddof=1)
    # For bad (zero or non-finite) sd columns, divide by 1.0 → (x - mu) is returned.
    bad = (sd == 0) | ~np.isfinite(sd)
    safe_sd = np.where(bad, 1.0, sd)
    out = (block - mu) / safe_sd
    return _restore_zero_cols(out, block, _zero_col_mask(block))


def _rank_block(block: np.ndarray, ties_method: TiesMethod) -> np.ndarray:
    """Column-wise percentile rank over a block."""
    # rankdata supports an axis kwarg from SciPy 1.10+.
    ranks = rankdata(block, method=ties_method, axis=0)
    out = ranks / block.shape[0]
    return _restore_zero_cols(out, block, _zero_col_mask(block))


def _clr_block(block: np.ndarray) -> np.ndarray:
    """Column-wise Seurat-flavor CLR: geom_mean = expm1(sum(log1p(x[x>0])) / n).

    Matches :func:`_apply_column_wise` + scalar ``_clr``: if a column has no
    positive entries, or the computed geometric mean is 0, the column is
    returned unchanged.
    """
    pos_mask = block > 0
    # Replace non-positive entries with 0 before log1p so we don't emit
    # spurious RuntimeWarnings from log1p(negative) (log1p(0) == 0, so the
    # masked-out contribution to the sum is exactly 0).
    pos_vals = np.where(pos_mask, block, 0.0)
    sum_log1p_pos = np.log1p(pos_vals).sum(axis=0)
    # The scalar version divides by ``x.size`` (full column length), not pos count.
    n = block.shape[0]
    geom_mean = np.expm1(sum_log1p_pos / n)
    safe_gm = np.where(geom_mean == 0, 1.0, geom_mean)
    # Only compute log1p on block/safe_gm where the column is valid; use a
    # safe divisor to suppress warnings from invalid columns whose output we
    # discard anyway.
    with np.errstate(invalid="ignore"):
        transformed = np.log1p(block / safe_gm)
    valid = pos_mask.any(axis=0) & (geom_mean != 0)
    out = np.where(valid, transformed, block)
    return _restore_zero_cols(out, block, _zero_col_mask(block))


def _clr_mean_block(block: np.ndarray) -> np.ndarray:
    """Column-wise CLR with geom_mean = expm1(mean(log1p(x[x>=0])))."""
    nn_mask = block >= 0
    count = nn_mask.sum(axis=0)
    # Mask out entries < 0 before log1p to match the scalar version and avoid
    # log1p(negative) warnings.
    nn_vals = np.where(nn_mask, block, 0.0)
    sum_log1p = np.log1p(nn_vals).sum(axis=0)
    safe_count = np.where(count > 0, count, 1)
    mean_log1p = sum_log1p / safe_count
    geom_mean = np.expm1(mean_log1p)
    safe_gm = np.where(geom_mean == 0, 1.0, geom_mean)
    with np.errstate(invalid="ignore"):
        transformed = np.log((block + 1.0) / safe_gm)
    valid = (count > 0) & (geom_mean != 0)
    out = np.where(valid, transformed, block)
    return _restore_zero_cols(out, block, _zero_col_mask(block))


def _clr_med_block(block: np.ndarray) -> np.ndarray:
    """Column-wise CLR using the median as the center."""
    m = np.nanmedian(block, axis=0)
    m = np.where(m == 0, 1.0, m)
    out = np.log((block + 1.0) / m)
    return _restore_zero_cols(out, block, _zero_col_mask(block))


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
    out = X.astype(float, copy=True) if X.dtype != np.float64 else X.copy()
    q_levels = np.linspace(0.0, 1.0, n_quantiles)
    # Reference quantiles across all cells, computed once for every marker.
    ref_q_all = np.quantile(X, q_levels, axis=0)  # (n_q, n_markers)
    for b in np.unique(batches):
        mask = batches == b
        Xb = X[mask]
        # Per-batch quantiles for every marker in one call.
        batch_q_all = np.quantile(Xb, q_levels, axis=0)  # (n_q, n_markers)
        for j in range(X.shape[1]):
            bq = batch_q_all[:, j]
            # PCHIP requires strictly increasing x. Fall back to identity if degenerate.
            if np.any(np.diff(bq) <= 0):
                continue
            spline = PchipInterpolator(bq, ref_q_all[:, j], extrapolate=True)
            out[mask, j] = spline(Xb[:, j])
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

    batch_values = check_obs_values_not_missing(
        adata,
        batch_key,
        context="normalize()",
    )
    markers = resolve_markers(adata, markers)
    if copy:
        adata = adata.copy()

    X = marker_matrix(
        adata,
        markers,
        layer=layer,
        require_finite=True,
        context="normalize()",
    )
    batches = np.asarray(batch_values.values)

    if method == "qnorm":
        new_X = _quantile_norm(X, batches)
    else:
        if method == "scale":
            block_fn = _scale_block
        elif method == "rank":
            def block_fn(b):
                return _rank_block(b, ties_method)
        elif method == "CLR":
            block_fn = _clr_mean_block
        elif method == "CLR_seu":
            block_fn = _clr_block
        elif method == "CLR_med":
            block_fn = _clr_med_block
        else:
            raise ValueError(f"Unknown normalization method: {method!r}")

        new_X = np.empty_like(X, dtype=X.dtype)
        # Preserve row order: apply per-batch, scatter results back.
        for b in pd.unique(batches):
            mask = batches == b
            new_X[mask] = block_fn(X[mask])

    set_marker_matrix(adata, markers, new_X, layer=layer)
    return adata if copy else None
