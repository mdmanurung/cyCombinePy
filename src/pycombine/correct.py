"""Per-cluster ComBat correction.

Port of ``correct_data`` in ``R/02_batch_correct.R:356-544``. The AnnData is split
by its SOM cluster label, each sub-group is corrected with
:func:`pycombine.combat.run_combat`, and results are stitched back in the original
row order. Values are capped to the per-cluster min/max of the input (matching R
lines 524-531).
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from anndata import AnnData

from pycombine._utils import (
    check_confound,
    check_obs_key,
    marker_matrix,
    resolve_markers,
    set_marker_matrix,
)
from pycombine.combat import run_combat

CORRECTED_LAYER = "cycombine_corrected"


def _build_model_matrix(
    df_sub: pd.DataFrame,
    covar: str | None,
    anchor: str | None,
) -> np.ndarray | None:
    """Build a design matrix (sans intercept) from covar and/or anchor columns.

    Uses :mod:`formulaic` if available so the output matches R's
    ``stats::model.matrix`` (treatment contrasts, drop first level). Falls back
    to pandas ``get_dummies`` otherwise.
    """
    terms = []
    if covar is not None:
        terms.append(covar)
    if anchor is not None:
        terms.append(anchor)
    if not terms:
        return None

    try:
        from formulaic import model_matrix

        formula = " + ".join(terms)
        mm = model_matrix(formula, df_sub.assign(**{t: df_sub[t].astype("category") for t in terms}))
        arr = np.asarray(mm, dtype=float)
        # Drop the intercept column so we hand inmoose a pure covariate block.
        if arr.shape[1] and np.all(arr[:, 0] == 1):
            arr = arr[:, 1:]
        return arr if arr.size else None
    except ImportError:
        pieces = []
        for t in terms:
            pieces.append(
                pd.get_dummies(df_sub[t].astype("category"), drop_first=True).to_numpy(dtype=float)
            )
        arr = np.hstack(pieces) if pieces else None
        return arr if arr is not None and arr.size else None


def _resolve_num_factors(
    series: pd.Series,
    batch: pd.Series,
    design: np.ndarray | None,
) -> int:
    """Return effective number of factor levels, mirroring R lines 455-506.

    - 1 if the covariate is confounded with batch
    - 1 if the cluster is heavily skewed to a single level
    - else the number of distinct levels.
    """
    if check_confound(batch, design):
        return 1
    counts = series.value_counts()
    total = counts.sum()
    n = counts.size
    if total < counts.max() + n * 5:
        return 1
    return n


def correct_data(
    adata: AnnData,
    label_key: str = "cycombine_som",
    markers: Iterable[str] | None = None,
    batch_key: str = "batch",
    covar: str | None = None,
    anchor: str | None = None,
    parametric: bool = True,
    ref_batch=None,
    layer: str | None = None,
    out_layer: str = CORRECTED_LAYER,
    copy: bool = False,
) -> AnnData | None:
    """Per-cluster ComBat batch correction.

    Parameters
    ----------
    adata
        AnnData with a cluster label in ``adata.obs[label_key]`` and a batch in
        ``adata.obs[batch_key]``.
    label_key
        Column in ``adata.obs`` with the SOM cluster id (from :func:`create_som`).
    markers
        Var names to correct. If ``None``, uses :func:`pycombine.get_markers`.
    batch_key
        Column in ``adata.obs`` giving the batch assignment.
    covar, anchor
        Optional ``adata.obs`` columns used as ComBat covariates. Skew- and
        confound-detection follow the R logic at lines 455-506.
    parametric
        Parametric vs. non-parametric ComBat prior.
    ref_batch
        Optional reference batch that is kept unchanged.
    layer
        If given, read the uncorrected matrix from this layer rather than ``X``.
    out_layer
        Name of the layer to store the corrected matrix in.
    copy
        If True, return a corrected copy; otherwise mutate in place.
    """
    check_obs_key(adata, batch_key)
    check_obs_key(adata, label_key)
    if covar is not None:
        check_obs_key(adata, covar)
    if anchor is not None:
        check_obs_key(adata, anchor)

    markers = resolve_markers(adata, markers)
    if copy:
        adata = adata.copy()

    X = marker_matrix(adata, markers, layer=layer)  # (n_cells, n_markers)
    n_cells = X.shape[0]

    labels = adata.obs[label_key].astype(str).to_numpy()
    batches = adata.obs[batch_key].astype(str).to_numpy()

    corrected = X.copy()

    for lab in pd.unique(labels):
        idx = np.where(labels == lab)[0]
        if idx.size == 0:
            continue
        sub_X = X[idx]  # (n_sub, n_markers)
        sub_batch = pd.Series(batches[idx])

        uniq_batches = sub_batch.unique()
        if uniq_batches.size <= 1:
            # Only one batch in this cluster — nothing to correct. (R lines 448-452)
            continue

        sub_df = adata.obs.iloc[idx]

        # Covar / anchor handling: determine effective level count
        num_covar = 1
        if covar is not None:
            cov_design = _build_model_matrix(sub_df, covar, None)
            num_covar = _resolve_num_factors(sub_df[covar], sub_batch, cov_design)

        num_anchor = 1
        if anchor is not None:
            anc_design = _build_model_matrix(sub_df, None, anchor)
            num_anchor = _resolve_num_factors(sub_df[anchor], sub_batch, anc_design)

        # If both are non-trivial, check that their combination is not confounded
        # with batch; if it is, drop anchor (R prioritises covar, lines 489-495).
        if num_covar > 1 and num_anchor > 1:
            joint = _build_model_matrix(sub_df, covar, anchor)
            if check_confound(sub_batch, joint):
                num_anchor = 1

        eff_covar = covar if num_covar > 1 else None
        eff_anchor = anchor if num_anchor > 1 else None
        mod = _build_model_matrix(sub_df, eff_covar, eff_anchor)

        # inmoose expects (n_features, n_samples)
        x_t = sub_X.T
        try:
            corrected_sub = run_combat(
                x_t,
                batch=sub_batch.values,
                mod=mod,
                parametric=parametric,
                ref_batch=ref_batch,
            ).T
        except Exception as exc:  # pragma: no cover
            # If ComBat fails inside a cluster (e.g. singular cov), leave untouched.
            # This matches the spirit of R's skip-on-confound handling.
            import warnings

            warnings.warn(
                f"ComBat failed for cluster {lab!r} ({exc}); leaving uncorrected.",
                RuntimeWarning,
            )
            continue

        # Cap to per-marker min/max within this cluster (R lines 524-531).
        lo = sub_X.min(axis=0)
        hi = sub_X.max(axis=0)
        corrected_sub = np.clip(corrected_sub, lo, hi)

        corrected[idx] = corrected_sub

    set_marker_matrix(adata, markers, corrected, layer=out_layer)
    return adata if copy else None
