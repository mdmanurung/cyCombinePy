"""High-level ``batch_correct`` orchestrator.

Port of ``batch_correct`` in ``R/02_batch_correct.R:66-210``. Runs the full
cyCombine pipeline: batch-wise normalize → SOM clustering → per-cluster ComBat
correction. Supports iterative correction with multiple SOM grid sizes by
passing ``xdim``/``ydim`` as sequences.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np
from anndata import AnnData

from pycombine._utils import marker_matrix, resolve_markers, set_marker_matrix
from pycombine.cluster import create_som
from pycombine.correct import CORRECTED_LAYER, correct_data
from pycombine.normalize import NormMethod, TiesMethod, normalize


def _as_list(v) -> list:
    if isinstance(v, (list, tuple, np.ndarray)):
        return list(v)
    return [v]


def batch_correct(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    batch_key: str = "batch",
    label_key: str = "cycombine_som",
    xdim: int | Sequence[int] = 8,
    ydim: int | Sequence[int] = 8,
    rlen: int = 10,
    seed: int = 473,
    n_clusters: int | None = None,
    norm_method: NormMethod = "scale",
    ties_method: TiesMethod = "average",
    covar: str | None = None,
    anchor: str | None = None,
    ref_batch=None,
    parametric: bool = True,
    out_layer: str = CORRECTED_LAYER,
    copy: bool = False,
) -> AnnData | None:
    """Full cyCombine pipeline: normalize → SOM → per-cluster ComBat.

    Parameters
    ----------
    adata
        Input AnnData. ``adata.X`` is assumed to already be on an appropriate
        scale (e.g. post-asinh for cytometry).
    markers
        Var names to normalize/cluster/correct. Defaults to
        :func:`pycombine.get_markers`.
    batch_key
        Column in ``adata.obs`` holding batch assignments.
    label_key
        Column in ``adata.obs`` to write cluster labels to.
    xdim, ydim
        SOM grid dimensions. Sequences trigger iterative correction: for each
        ``(x, y)`` pair, re-normalize, re-cluster, and re-correct.
    rlen
        SOM training passes (forwarded to FlowSOM if supported).
    seed
        FlowSOM random seed.
    n_clusters
        If set, metacluster the SOM nodes into this many clusters.
    norm_method
        Normalization method used for clustering. See :func:`pycombine.normalize`.
    ties_method
        Tie-breaking rule for ``norm_method="rank"``.
    covar, anchor, ref_batch, parametric
        Forwarded to :func:`pycombine.correct_data`.
    out_layer
        Layer name to store the corrected matrix in.
    copy
        If True, return a corrected copy; otherwise mutate in place.
    """
    if copy:
        adata = adata.copy()

    markers = resolve_markers(adata, markers)
    xdims = _as_list(xdim)
    ydims = _as_list(ydim)
    if len(xdims) != len(ydims):
        raise ValueError("xdim and ydim must have the same length")

    # Snapshot the starting matrix — clustering operates on a normalized copy but
    # correction is applied to the current "working" expression.
    working = marker_matrix(adata, markers).copy()

    for x, y in zip(xdims, ydims):
        # 1. Build a normalized view for clustering without touching `adata.X`.
        norm_view = adata.copy()
        set_marker_matrix(norm_view, markers, working)
        normalize(
            norm_view,
            markers=markers,
            method=norm_method,
            batch_key=batch_key,
            ties_method=ties_method,
        )

        # 2. Cluster on the normalized view, store labels on `adata.obs`.
        create_som(
            norm_view,
            markers=markers,
            xdim=x,
            ydim=y,
            n_clusters=n_clusters,
            seed=seed,
            rlen=rlen,
            label_key=label_key,
        )
        adata.obs[label_key] = norm_view.obs[label_key].values

        # 3. Correct the *working* (not normalized) values per cluster.
        work_view = adata.copy()
        set_marker_matrix(work_view, markers, working)
        correct_data(
            work_view,
            label_key=label_key,
            markers=markers,
            batch_key=batch_key,
            covar=covar,
            anchor=anchor,
            parametric=parametric,
            ref_batch=ref_batch,
            out_layer=out_layer,
        )
        working = marker_matrix(work_view, markers, layer=out_layer)

    set_marker_matrix(adata, markers, working, layer=out_layer)
    return adata if copy else None
