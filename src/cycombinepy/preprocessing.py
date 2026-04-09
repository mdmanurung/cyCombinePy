"""Preprocessing utilities: asinh transform with optional derandomization."""

from __future__ import annotations

from typing import Iterable

import numpy as np
from anndata import AnnData

from cycombinepy._utils import marker_matrix, resolve_markers, set_marker_matrix


def transform_asinh(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    cofactor: float = 5.0,
    derand: bool = True,
    reverse: bool = False,
    layer: str | None = None,
    copy: bool = False,
    seed: int | None = None,
) -> AnnData | None:
    """Asinh-transform marker columns of ``adata`` with an optional derandomization step.

    Port of ``transform_asinh`` in ``R/01_prepare_data.R:375``. Derandomization mirrors
    ``randomize_matrix`` in ``R/utils_helper.R``: take ``ceil(x)`` then subtract a
    uniform draw on ``[0, 0.9999]`` before the asinh.

    Parameters
    ----------
    adata
        AnnData with raw expression in ``adata.X`` (or ``adata.layers[layer]``).
    markers
        Var names to transform. If ``None``, uses :func:`cycombinepy.get_markers`.
    cofactor
        Asinh cofactor. Common values: 5 (CyTOF), 150 (flow), 6000 (spectral).
    derand
        If True, apply the derandomization before asinh.
    reverse
        If True, apply ``sinh(x) * cofactor`` instead (inverse transform).
    layer
        If given, read / write that layer instead of ``adata.X``.
    copy
        If True, return a copy; otherwise modify in place and return None.
    seed
        Seed for derandomization RNG.
    """
    markers = resolve_markers(adata, markers)
    if copy:
        adata = adata.copy()

    X = marker_matrix(adata, markers, layer=layer)

    if reverse:
        X = np.sinh(X) * cofactor
    else:
        if derand:
            rng = np.random.default_rng(seed)
            X = np.ceil(X) - rng.uniform(0.0, 0.9999, size=X.shape)
            X[X < 0] = 0.0
        X = np.arcsinh(X / cofactor)

    set_marker_matrix(adata, markers, X, layer=layer)
    return adata if copy else None
