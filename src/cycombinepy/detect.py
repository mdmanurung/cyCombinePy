"""Batch-effect detection utilities.

Port of ``detect_batch_effect_express`` / ``detect_batch_effect`` from
``R/detect_batch_effect.R``. The Python versions return a dict of matplotlib
figures (or save them to ``out_dir``) rather than printing the R plots.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from anndata import AnnData

from cycombinepy._utils import (
    check_obs_values_not_missing,
    marker_matrix,
    resolve_markers,
)
from cycombinepy.evaluate import compute_emd, compute_mad


def _ensure_single_cluster_label(adata: AnnData, key: str = "_cycombine_all") -> str:
    """Add a trivial cluster column so compute_emd/mad treat the data as one group."""
    if key not in adata.obs.columns:
        adata.obs[key] = "all"
    return key


def detect_batch_effect_express(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    batch_key: str = "batch",
    sample_key: str | None = "sample",
    downsample: int | None = None,
    out_dir: str | os.PathLike | None = None,
    seed: int = 472,
) -> dict:
    """Quick 3-panel batch-effect summary (EMD heatmap, density, MDS).

    Returns a ``dict`` of matplotlib figures keyed by ``"emd"``, ``"density"``
    and ``"mds"``. If ``out_dir`` is given, the figures are saved as PNGs and
    the dict is still returned for further inspection.

    Matches ``detect_batch_effect_express`` in ``R/detect_batch_effect.R``.
    """
    check_obs_values_not_missing(
        adata,
        batch_key,
        context="detect_batch_effect_express()",
    )
    if sample_key is not None and sample_key in adata.obs.columns:
        check_obs_values_not_missing(
            adata,
            sample_key,
            context="detect_batch_effect_express()",
        )
    markers = resolve_markers(adata, markers)

    if downsample is not None and adata.n_obs > downsample:
        rng = np.random.default_rng(seed)
        idx = rng.choice(adata.n_obs, downsample, replace=False)
        adata = adata[idx].copy()

    X = marker_matrix(
        adata,
        markers,
        require_finite=True,
        context="detect_batch_effect_express()",
    )

    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "Detection plots require matplotlib + seaborn; install with "
            "`pip install matplotlib seaborn`."
        ) from exc

    cluster_key = _ensure_single_cluster_label(adata)
    emd_df = compute_emd(
        adata, cell_key=cluster_key, batch_key=batch_key, markers=markers
    )

    # 1. EMD heatmap (mean per marker across batch pairs)
    pivot = emd_df.groupby("marker")["emd"].mean().to_frame("mean_emd")
    fig_emd, ax_emd = plt.subplots(figsize=(4, max(2, 0.3 * len(pivot))))
    sns.heatmap(pivot, annot=True, cmap="viridis", ax=ax_emd)
    ax_emd.set_title("Mean EMD per marker (between batches)")

    # 2. Density plots per marker, colored by batch
    long = pd.DataFrame(X, columns=markers)
    long[batch_key] = adata.obs[batch_key].values
    long = long.melt(id_vars=[batch_key], var_name="marker", value_name="value")
    g = sns.FacetGrid(long, col="marker", hue=batch_key, col_wrap=4, sharey=False)
    g.map(sns.kdeplot, "value", fill=True, alpha=0.3, common_norm=False)
    g.add_legend()
    fig_density = g.figure

    # 3. MDS of per-sample median expression
    fig_mds, ax_mds = plt.subplots(figsize=(4, 4))
    if sample_key is not None and sample_key in adata.obs.columns:
        df = pd.DataFrame(X, columns=markers)
        df[sample_key] = adata.obs[sample_key].values
        df[batch_key] = adata.obs[batch_key].values
        medians = df.groupby(sample_key)[list(markers)].median()
        sample_batch = df.groupby(sample_key)[batch_key].first()

        from sklearn.manifold import MDS

        mds = MDS(
            n_components=2,
            random_state=seed,
            normalized_stress="auto",
            init="random",
        )
        coords = mds.fit_transform(medians.values)
        for b in sample_batch.unique():
            mask = (sample_batch == b).values
            ax_mds.scatter(coords[mask, 0], coords[mask, 1], label=str(b))
        ax_mds.legend(title=batch_key)
        ax_mds.set_xlabel("MDS1")
        ax_mds.set_ylabel("MDS2")
        ax_mds.set_title("MDS of per-sample medians")
    else:
        ax_mds.text(0.5, 0.5, "No sample_key provided", ha="center", va="center")

    figs = {"emd": fig_emd, "density": fig_density, "mds": fig_mds}
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, fig in figs.items():
            fig.savefig(out_dir / f"detect_{name}.png", dpi=120, bbox_inches="tight")
    return figs


def detect_batch_effect(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    batch_key: str = "batch",
    sample_key: str | None = "sample",
    downsample: int | None = None,
    out_dir: str | os.PathLike | None = None,
    seed: int = 472,
) -> dict:
    """Extended batch-effect diagnostic: express output plus UMAP and MAD.

    Matches ``detect_batch_effect`` in ``R/detect_batch_effect.R``.
    """
    check_obs_values_not_missing(
        adata,
        batch_key,
        context="detect_batch_effect()",
    )
    if sample_key is not None and sample_key in adata.obs.columns:
        check_obs_values_not_missing(
            adata,
            sample_key,
            context="detect_batch_effect()",
        )
    markers = resolve_markers(adata, markers)
    marker_matrix(
        adata,
        markers,
        require_finite=True,
        context="detect_batch_effect()",
    )

    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover
        raise ImportError("scanpy + matplotlib required for detect_batch_effect.") from exc

    figs = detect_batch_effect_express(
        adata,
        markers=markers,
        batch_key=batch_key,
        sample_key=sample_key,
        downsample=downsample,
        out_dir=None,
        seed=seed,
    )

    # UMAP colored by batch
    a = adata.copy()
    sc.pp.pca(a, n_comps=min(20, a.n_vars - 1))
    sc.pp.neighbors(a, n_neighbors=15)
    sc.tl.umap(a, random_state=seed)
    fig_umap, ax = plt.subplots(figsize=(5, 5))
    for b in a.obs[batch_key].unique():
        mask = (a.obs[batch_key] == b).to_numpy()
        ax.scatter(
            a.obsm["X_umap"][mask, 0],
            a.obsm["X_umap"][mask, 1],
            s=3,
            label=str(b),
            alpha=0.6,
        )
    ax.legend(title=batch_key)
    ax.set_title("UMAP (uncorrected) colored by batch")
    figs["umap"] = fig_umap

    # Per-marker MAD summary
    cluster_key = _ensure_single_cluster_label(adata)
    mad_df = compute_mad(adata, cell_key=cluster_key, batch_key=batch_key, markers=markers)
    import seaborn as sns

    fig_mad, ax_mad = plt.subplots(figsize=(6, max(3, 0.3 * mad_df["marker"].nunique())))
    sns.barplot(data=mad_df, x="mad", y="marker", hue="batch", ax=ax_mad)
    ax_mad.set_title("MAD per marker per batch")
    figs["mad"] = fig_mad

    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, fig in figs.items():
            fig.savefig(out_dir / f"detect_{name}.png", dpi=120, bbox_inches="tight")

    return figs
