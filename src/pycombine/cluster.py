"""SOM clustering wrapper around FlowSOM_Python.

Port of ``create_som`` in ``R/02_batch_correct.R:210``. Unlike the R version
(which supports kohonen / FlowSOM / FuseSOM / kmeans backends), the Python port
exclusively uses saeyslab's ``FlowSOM`` package.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
from anndata import AnnData

from pycombine._utils import marker_matrix, resolve_markers


def create_som(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    xdim: int = 8,
    ydim: int = 8,
    n_clusters: int | None = None,
    seed: int = 473,
    rlen: int = 10,
    layer: str | None = None,
    label_key: str = "cycombine_som",
    copy: bool = False,
) -> AnnData | None:
    """Train a FlowSOM on the marker columns of ``adata`` and store cluster labels.

    Cells are mapped to either the raw SOM node (default) or a metacluster if
    ``n_clusters`` is supplied. The resulting 1-indexed integer labels are written
    to ``adata.obs[label_key]`` as a ``category``.

    Parameters
    ----------
    adata
        AnnData to cluster.
    markers
        Var names used for clustering. Defaults to :func:`pycombine.get_markers`.
    xdim, ydim
        SOM grid dimensions. Default 8x8 matches cyCombine.
    n_clusters
        If set, consensus-metacluster the SOM nodes into this many clusters.
    seed
        Random seed passed to FlowSOM.
    rlen
        Number of training passes. (Kept for parity with R; forwarded if the
        installed FlowSOM version accepts it.)
    layer
        Read features from this layer rather than ``adata.X``.
    label_key
        Name of the ``adata.obs`` column to write labels into.
    copy
        If True, return a modified copy rather than mutating in place.
    """
    try:
        import flowsom as fs
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "create_som requires FlowSOM_Python. Install it with "
            "`pip install flowsom`."
        ) from exc

    markers = resolve_markers(adata, markers)
    if copy:
        adata = adata.copy()

    X = marker_matrix(adata, markers, layer=layer)

    # FlowSOM_Python accepts an AnnData; build a minimal one with just the markers
    # we care about so that the SOM is trained on the correct columns.
    import anndata as ad

    feat = ad.AnnData(X=X.astype(np.float32), var={"marker": markers})
    feat.var_names = list(markers)

    n_clus = n_clusters if n_clusters is not None else xdim * ydim
    try:
        fsom = fs.FlowSOM(
            feat,
            cols_to_use=list(markers),
            xdim=xdim,
            ydim=ydim,
            n_clusters=n_clus,
            rlen=rlen,
            seed=seed,
        )
    except TypeError:
        # Older versions may not accept ``rlen``.
        fsom = fs.FlowSOM(
            feat,
            cols_to_use=list(markers),
            xdim=xdim,
            ydim=ydim,
            n_clusters=n_clus,
            seed=seed,
        )

    if n_clusters is None:
        labels = np.asarray(fsom.get_cell_data().obs["clustering"]).astype(int)
    else:
        labels = np.asarray(fsom.get_cell_data().obs["metaclustering"]).astype(int)

    # Store as 1-indexed categorical labels for parity with R.
    labels = labels + 1 if labels.min() == 0 else labels
    adata.obs[label_key] = labels.astype(str)
    adata.obs[label_key] = adata.obs[label_key].astype("category")

    return adata if copy else None
