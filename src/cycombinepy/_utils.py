"""Internal utilities shared across cycombinepy modules."""

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
MAX_MISSING_OBS_INDICES = 5


def check_var_names_unique(adata: AnnData) -> None:
    """Raise if ``adata.var_names`` contains duplicate marker names."""
    if adata.var_names.is_unique:
        return
    duplicates = adata.var_names[adata.var_names.duplicated()].unique().tolist()
    raise ValueError(f"adata.var_names must be unique; duplicate names: {duplicates}")


def check_layer_key(
    adata: AnnData,
    layer: str | None,
    location: str = "adata.layers",
) -> None:
    """Raise if ``layer`` is not present in ``adata.layers``."""
    if layer is None:
        return
    if layer not in adata.layers:
        raise KeyError(f"Layer {layer!r} was not found in {location}")


def check_obs_values_not_missing(
    adata: AnnData,
    key: str,
    *,
    context: str,
) -> pd.Series:
    """Return ``adata.obs[key]`` after rejecting actual pandas missing values."""
    check_obs_key(adata, key)
    values = adata.obs[key]
    missing = values.isna()
    if missing.any():
        n_missing = int(missing.sum())
        missing_rows = values.index[missing][:MAX_MISSING_OBS_INDICES].tolist()
        raise ValueError(
            f"{context} requires non-missing values in adata.obs[{key!r}]; "
            f"found {n_missing} missing value(s); "
            f"first missing rows: {missing_rows}"
        )
    return values


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
    check_var_names_unique(adata)
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


def _marker_indices(adata: AnnData, markers: list[str]) -> np.ndarray:
    """Resolve marker names to an integer column index array."""
    check_var_names_unique(adata)
    get_loc = adata.var_names.get_loc
    return np.fromiter((get_loc(m) for m in markers), dtype=np.intp, count=len(markers))


def marker_matrix(
    adata: AnnData,
    markers: list[str],
    layer: str | None = None,
    dtype=float,
    require_finite: bool = False,
    context: str = "marker_matrix()",
) -> np.ndarray:
    """Extract a (n_cells, n_markers) dense array for the given markers.

    Parameters
    ----------
    dtype
        Target dtype. Defaults to ``float`` (float64) for backwards
        compatibility. Pass ``None`` to preserve the source dtype, which avoids
        an unnecessary upcast when the source is already float32/float64.
    require_finite
        If true, reject NaN or infinite values after dtype conversion.
    context
        Operation name included in finite-value validation errors.
    """
    idx = _marker_indices(adata, markers)
    check_layer_key(adata, layer)
    X = adata.X if layer is None else adata.layers[layer]
    # Fancy indexing always returns a fresh copy, so no extra .copy() is needed.
    if hasattr(X, "toarray"):
        # Only materialize the requested columns from sparse to dense.
        sub = np.asarray(X[:, idx].toarray())
    else:
        sub = np.asarray(X)[:, idx]
    if dtype is not None and sub.dtype != np.dtype(dtype):
        sub = sub.astype(dtype, copy=False)
    if require_finite:
        try:
            finite = np.isfinite(sub)
        except TypeError as exc:
            raise ValueError(f"{context} requires finite marker values") from exc
        if not np.all(finite):
            raise ValueError(f"{context} requires finite marker values")
    return sub


def _can_write_in_place(arr) -> bool:
    """True if ``arr`` is a writeable dense ndarray we can assign into."""
    return (
        isinstance(arr, np.ndarray)
        and arr.flags.writeable
        and np.issubdtype(arr.dtype, np.floating)
    )


def set_marker_matrix(
    adata: AnnData,
    markers: list[str],
    values: np.ndarray,
    layer: str | None = None,
) -> None:
    """Write ``values`` (n_cells, n_markers) back into ``adata.X`` or a layer.

    Writes in place when the target is already a writeable dense float ndarray
    (and the AnnData is not a view), avoiding full-matrix copies on the hot path.
    Falls back to a single dense-float materialization for sparse/non-float
    targets.
    """
    idx = _marker_indices(adata, markers)

    if layer is None:
        target = adata.X
        if _can_write_in_place(target) and not adata.is_view:
            # In-place column write: no allocation.
            if values.dtype != target.dtype:
                values = values.astype(target.dtype, copy=False)
            target[:, idx] = values
            return
        # Fallback: single materialization to a writable float ndarray.
        X = as_dense(target).astype(float, copy=True)
        X[:, idx] = values
        adata.X = X
        return

    existing = adata.layers.get(layer) if layer in adata.layers else None
    if existing is not None and _can_write_in_place(existing) and not adata.is_view:
        if values.dtype != existing.dtype:
            values = values.astype(existing.dtype, copy=False)
        existing[:, idx] = values
        return

    if existing is None:
        # Initialise the layer from X with a single dense-float materialization.
        L = as_dense(adata.X).astype(float, copy=True)
    else:
        L = as_dense(existing).astype(float, copy=True)
    L[:, idx] = values
    adata.layers[layer] = L


def check_confound(batch, mod: np.ndarray | None = None) -> bool:
    """Return True if ``batch`` is confounded with ``mod``.

    Port of ``check_confound`` in ``R/utils_helper.R`` (adapted from ``sva::ComBat``).
    Tests for rank deficiency of ``[batch_dummies | mod]`` after dropping intercept-
    like columns.
    """
    # Categorical codes avoid the per-call pandas get_dummies overhead.
    cat = pd.Categorical(batch)
    codes = cat.codes
    n_levels = len(cat.categories)
    if n_levels == 0:
        batchmod = np.zeros((len(codes), 0), dtype=np.float64)
    else:
        batchmod = np.eye(n_levels, dtype=np.float64)[codes]

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
