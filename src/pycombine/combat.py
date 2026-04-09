"""Thin wrapper around ``inmoose.pycombat.pycombat_norm``.

Port of ``.combat`` in ``R/02_batch_correct.R:548``. We keep the R convention:
the input is a ``(n_features, n_samples)`` matrix so both the calling code in
:mod:`pycombine.correct` and ``inmoose`` see the same shape.
"""

from __future__ import annotations

import numpy as np


def run_combat(
    x: np.ndarray,
    batch,
    mod: np.ndarray | None = None,
    parametric: bool = True,
    ref_batch=None,
) -> np.ndarray:
    """Run ComBat on a ``(n_features, n_samples)`` matrix.

    Parameters
    ----------
    x
        Expression matrix with features in rows and cells (samples) in columns,
        matching both R's ``sva::ComBat`` input and ``inmoose.pycombat_norm``.
    batch
        Array-like of length ``n_samples`` giving the batch for each cell.
    mod
        Optional covariate model matrix ``(n_samples, n_covar)``. Built with
        :func:`formulaic.model_matrix` by the caller.
    parametric
        If True, use the parametric empirical-Bayes prior.
    ref_batch
        Optional reference batch that is kept unchanged.
    """
    from inmoose.pycombat import pycombat_norm

    batch_arr = np.asarray([str(b) for b in batch])
    result = pycombat_norm(
        counts=x,
        batch=batch_arr,
        covar_mod=mod,
        par_prior=parametric,
        ref_batch=ref_batch,
    )
    return np.asarray(result)
