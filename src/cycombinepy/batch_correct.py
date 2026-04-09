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

from cycombinepy._utils import marker_matrix, resolve_markers, set_marker_matrix
from cycombinepy.cluster import create_som
from cycombinepy.correct import CORRECTED_LAYER, correct_data
from cycombinepy.normalize import NormMethod, TiesMethod, normalize


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
        :func:`cycombinepy.get_markers`.
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
        Normalization method used for clustering. See :func:`cycombinepy.normalize`.
    ties_method
        Tie-breaking rule for ``norm_method="rank"``.
    covar, anchor, ref_batch, parametric
        Forwarded to :func:`cycombinepy.correct_data`.
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

    # Working copy of the marker matrix that accumulates corrections between
    # iterations. Clustering sees a normalized view; correction sees the current
    # unnormalized working state.
    working = marker_matrix(adata, markers).copy()
    scratch = adata.copy()

    for x, y in zip(xdims, ydims):
        # Normalize + cluster on a fresh normalized view.
        set_marker_matrix(scratch, markers, working)
        normalize(
            scratch,
            markers=markers,
            method=norm_method,
            batch_key=batch_key,
            ties_method=ties_method,
        )
        create_som(
            scratch,
            markers=markers,
            xdim=x,
            ydim=y,
            n_clusters=n_clusters,
            seed=seed,
            rlen=rlen,
            label_key=label_key,
        )
        adata.obs[label_key] = scratch.obs[label_key].values

        # Correct the (unnormalized) working state per cluster.
        set_marker_matrix(scratch, markers, working)
        correct_data(
            scratch,
            label_key=label_key,
            markers=markers,
            batch_key=batch_key,
            covar=covar,
            anchor=anchor,
            parametric=parametric,
            ref_batch=ref_batch,
            out_layer=out_layer,
        )
        working = marker_matrix(scratch, markers, layer=out_layer)

    set_marker_matrix(adata, markers, working, layer=out_layer)
    return adata if copy else None
