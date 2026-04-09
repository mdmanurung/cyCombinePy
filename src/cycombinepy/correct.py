"""Per-cluster ComBat correction.

Port of ``correct_data`` in ``R/02_batch_correct.R:356-544``. The AnnData is split
by its SOM cluster label, each sub-group is corrected with
:func:`cycombinepy.combat.run_combat`, and results are stitched back in the original
row order. Values are capped to the per-cluster min/max of the input (matching R
lines 524-531).
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from anndata import AnnData

from cycombinepy._utils import (
    check_confound,
    check_obs_key,
    marker_matrix,
    resolve_markers,
    set_marker_matrix,
)
from cycombinepy.combat import run_combat

CORRECTED_LAYER = "cycombine_corrected"


def _build_model_matrix(
    df_sub: pd.DataFrame,
    covar: str | None,
    anchor: str | None,
) -> np.ndarray | None:
    """Build a design matrix (sans intercept) from covar and/or anchor columns.

    Uses :mod:`formulaic` to match R's ``stats::model.matrix`` (treatment
    contrasts, drop first level).
    """
    from formulaic import model_matrix

    terms = [t for t in (covar, anchor) if t is not None]
    if not terms:
        return None

    sub = df_sub[terms].astype("category")
    mm = np.asarray(model_matrix(" + ".join(terms), sub), dtype=float)
    # Drop the intercept column so we hand inmoose a pure covariate block.
    if mm.shape[1] and np.all(mm[:, 0] == 1):
        mm = mm[:, 1:]
    return mm if mm.size else None


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
        Var names to correct. If ``None``, uses :func:`cycombinepy.get_markers`.
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

    # Convert label/batch to categorical codes once; group rows by label via
    # a single stable argsort + np.split to avoid the per-cluster O(N)
    # boolean scans that the previous implementation did.
    label_cat = pd.Categorical(adata.obs[label_key].astype(str).to_numpy())
    batch_cat = pd.Categorical(adata.obs[batch_key].astype(str).to_numpy())
    label_codes = label_cat.codes
    batch_codes = batch_cat.codes
    batch_categories = np.asarray(batch_cat.categories)

    order = np.argsort(label_codes, kind="stable")
    # Boundaries between sorted label groups.
    sorted_codes = label_codes[order]
    boundaries = np.flatnonzero(np.diff(sorted_codes)) + 1
    cluster_index_groups = np.split(order, boundaries)

    # Pre-slice obs columns needed for covar/anchor as integer code arrays;
    # we only materialize a small per-cluster DataFrame when _build_model_matrix
    # is actually called.
    obs = adata.obs

    corrected = X.copy()

    for idx in cluster_index_groups:
        if idx.size == 0:
            continue
        sub_X = X[idx]  # (n_sub, n_markers)
        sub_batch_codes = batch_codes[idx]

        # Detect the set of distinct batches present in this cluster without
        # falling back to pandas.
        present_codes = np.unique(sub_batch_codes)
        if present_codes.size <= 1:
            # Only one batch in this cluster — nothing to correct. (R lines 448-452)
            continue

        sub_batch_values = batch_categories[sub_batch_codes]
        sub_batch = pd.Series(sub_batch_values)

        # Covar / anchor handling: determine effective level count.
        num_covar = 1
        num_anchor = 1
        sub_df = None  # lazy — only built when needed

        if covar is not None or anchor is not None:
            sub_df = obs.iloc[idx]

        if covar is not None:
            cov_design = _build_model_matrix(sub_df, covar, None)
            num_covar = _resolve_num_factors(sub_df[covar], sub_batch, cov_design)

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
        if eff_covar is None and eff_anchor is None:
            mod = None
        else:
            if sub_df is None:
                sub_df = obs.iloc[idx]
            mod = _build_model_matrix(sub_df, eff_covar, eff_anchor)

        # inmoose expects (n_features, n_samples) and is sensitive to float32
        # underflow in its EB priors → upcast to float64 at the ComBat boundary.
        x_t = np.ascontiguousarray(sub_X.T, dtype=np.float64)
        lab = label_cat.categories[label_codes[idx[0]]]
        try:
            corrected_sub = run_combat(
                x_t,
                batch=sub_batch_values,
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
        # Fuse the clip into a single maximum/minimum pair to avoid the extra
        # allocation from ``np.clip``.
        lo = sub_X.min(axis=0)
        hi = sub_X.max(axis=0)
        if corrected_sub.dtype != X.dtype:
            corrected_sub = corrected_sub.astype(X.dtype, copy=False)
        np.maximum(corrected_sub, lo, out=corrected_sub)
        np.minimum(corrected_sub, hi, out=corrected_sub)

        corrected[idx] = corrected_sub

    set_marker_matrix(adata, markers, corrected, layer=out_layer)
    return adata if copy else None
