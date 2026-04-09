"""Internal utilities shared across pycombine modules."""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from anndata import AnnData

# Mirror of cyCombine's cyCombine::non_markers default set
# (decoded from data/non_markers.rda).
DEFAULT_NON_MARKERS: tuple[str, ...] = (
    "LiveDead",
    "SSC",
    "FSC",
    "batch",
    "sample",
    "condition",
    "anchor",
    "som",
    "label",
    "id",
    "covar",
    "celltype",
    "model_prediction",
    "predicted_celltype",
    "cell",
    "cell_id",
)


def get_markers(
    adata: AnnData,
    non_markers: Iterable[str] | None = None,
) -> list[str]:
    """Return var_names that are not in the non-markers blacklist (case-insensitive).

    Mirrors ``cyCombine::get_markers`` in ``R/utils_helper.R`` using the default
    non-marker set from ``data/non_markers.rda``.
    """
    if non_markers is None:
        non_markers = DEFAULT_NON_MARKERS
    blacklist = {s.lower() for s in non_markers}
    return [v for v in adata.var_names if v.lower() not in blacklist]


def resolve_markers(
    adata: AnnData,
    markers: Iterable[str] | None,
) -> list[str]:
    """Normalize a user-supplied ``markers`` argument to a list of var_names."""
    if markers is None:
        return get_markers(adata)
    markers = list(markers)
    missing = [m for m in markers if m not in adata.var_names]
    if missing:
        raise KeyError(f"Markers not found in adata.var_names: {missing}")
    return markers


def check_obs_key(adata: AnnData, key: str, location: str = "adata.obs") -> None:
    """Raise if ``key`` is not a column of ``adata.obs``."""
    if key not in adata.obs.columns:
        raise KeyError(f'Column "{key}" was not found in {location}')


def as_dense(x) -> np.ndarray:
    """Return a dense ``np.ndarray`` from dense or sparse input."""
    if hasattr(x, "toarray"):
        return np.asarray(x.toarray())
    return np.asarray(x)


def marker_matrix(adata: AnnData, markers: list[str], layer: str | None = None) -> np.ndarray:
    """Extract a (n_cells, n_markers) dense float array for the given markers."""
    idx = [adata.var_names.get_loc(m) for m in markers]
    if layer is None:
        X = adata.X
    else:
        X = adata.layers[layer]
    X = as_dense(X)
    return np.asarray(X[:, idx], dtype=float)


def set_marker_matrix(
    adata: AnnData,
    markers: list[str],
    values: np.ndarray,
    layer: str | None = None,
) -> None:
    """Write ``values`` (n_cells, n_markers) back into ``adata.X`` or a layer."""
    idx = [adata.var_names.get_loc(m) for m in markers]
    if layer is None:
        X = as_dense(adata.X).astype(float, copy=True)
        X[:, idx] = values
        adata.X = X
    else:
        if layer not in adata.layers:
            adata.layers[layer] = as_dense(adata.X).astype(float, copy=True)
        L = as_dense(adata.layers[layer]).astype(float, copy=True)
        L[:, idx] = values
        adata.layers[layer] = L


def check_confound(batch, mod: np.ndarray | None = None) -> bool:
    """Return True if ``batch`` is confounded with ``mod``.

    Port of ``check_confound`` in ``R/utils_helper.R`` (adapted from ``sva::ComBat``).
    Tests for rank deficiency of ``[batch_dummies | mod]`` after dropping intercept-
    like columns.
    """
    batch = pd.Series(batch).astype("category")
    # one-hot (no intercept) batch model
    batchmod = pd.get_dummies(batch, drop_first=False).to_numpy(dtype=float)

    if mod is None:
        design = batchmod
    else:
        mod = np.asarray(mod, dtype=float)
        if mod.ndim == 1:
            mod = mod.reshape(-1, 1)
        design = np.hstack([batchmod, mod])

    # Drop all-ones columns (intercept-like), matching R's
    # `check <- apply(design, 2, function(x) all(x == 1))`.
    keep = ~np.all(design == 1, axis=0)
    design = design[:, keep]
    if design.size == 0:
        return False

    rank = np.linalg.matrix_rank(design)
    # Rank-deficient design after dropping intercept columns ⇒ confounded.
    # The R version branches on ncol vs n_batch for messaging, but the final
    # result is TRUE in every sub-branch.
    return bool(rank < design.shape[1])
