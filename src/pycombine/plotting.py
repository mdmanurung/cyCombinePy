"""Plotting helpers.

Pragmatic subset of ``R/utils_plotting.R``: density comparison, dimensionality
reduction, and an EMD heatmap summary.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from anndata import AnnData

from pycombine._utils import marker_matrix, resolve_markers


def plot_density(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    layer: str | None = "cycombine_corrected",
    batch_key: str = "batch",
    filename: str | os.PathLike | None = None,
):
    """Per-marker density plot colored by batch.

    If ``layer`` is set and present, the corrected distribution is overlaid
    alongside the uncorrected ``adata.X`` for easy before/after comparison.
    Mirrors ``plot_density`` in ``R/utils_plotting.R``.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    markers = resolve_markers(adata, markers)
    batches = adata.obs[batch_key].astype(str).values

    rows = []
    X = marker_matrix(adata, markers)
    for j, marker in enumerate(markers):
        for v, b in zip(X[:, j], batches):
            rows.append({"marker": marker, "value": v, batch_key: b, "kind": "uncorrected"})
    if layer is not None and layer in adata.layers:
        X2 = marker_matrix(adata, markers, layer=layer)
        for j, marker in enumerate(markers):
            for v, b in zip(X2[:, j], batches):
                rows.append({"marker": marker, "value": v, batch_key: b, "kind": "corrected"})

    df = pd.DataFrame(rows)
    g = sns.FacetGrid(
        df,
        col="marker",
        row="kind" if "corrected" in df["kind"].unique() else None,
        hue=batch_key,
        sharey=False,
        col_wrap=None if "corrected" in df["kind"].unique() else 4,
    )
    g.map(sns.kdeplot, "value", fill=True, alpha=0.3, common_norm=False)
    g.add_legend()

    if filename is not None:
        g.figure.savefig(filename, dpi=120, bbox_inches="tight")
    return g.figure


def plot_dimred(
    adata: AnnData,
    kind: str = "umap",
    color: str | list[str] = "batch",
    layer: str | None = None,
    seed: int = 0,
):
    """Thin wrapper around ``scanpy.pl.umap`` / ``scanpy.pl.pca``.

    If ``layer`` is provided, the PCA is computed on that layer so the plot
    reflects the corrected values.
    """
    import matplotlib.pyplot as plt
    import scanpy as sc

    a = adata.copy()
    if layer is not None:
        a.X = a.layers[layer]

    sc.pp.pca(a, n_comps=min(20, a.n_vars - 1))
    if kind == "pca":
        sc.pl.pca(a, color=color, show=False)
        return plt.gcf()
    sc.pp.neighbors(a, n_neighbors=15, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    sc.pl.umap(a, color=color, show=False)
    return plt.gcf()


def plot_emd_heatmap(emd_df: pd.DataFrame, filename: str | os.PathLike | None = None):
    """Heatmap of mean EMD per (cluster, marker) from :func:`compute_emd` output."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    pivot = emd_df.pivot_table(
        index="cluster", columns="marker", values="emd", aggfunc="mean"
    )
    fig, ax = plt.subplots(figsize=(max(4, 0.4 * pivot.shape[1]), max(3, 0.3 * pivot.shape[0])))
    sns.heatmap(pivot, annot=False, cmap="viridis", ax=ax)
    ax.set_title("Mean EMD per (cluster, marker)")
    if filename is not None:
        fig.savefig(filename, dpi=120, bbox_inches="tight")
    return fig
