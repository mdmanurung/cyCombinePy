"""Shared pytest fixtures for cycombinepy."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest


def make_normalization_exact_adata() -> ad.AnnData:
    X = np.array(
        [
            [1.0, 0.0],
            [1.0, 2.0],
            [3.0, 4.0],
            [10.0, 1.0],
            [10.0, 1.0],
            [20.0, 5.0],
        ]
    )
    obs = pd.DataFrame({"batch": ["a", "a", "a", "b", "b", "b"]})
    obs.index = obs.index.astype(str)
    out = ad.AnnData(X=X, obs=obs)
    out.var_names = ["CD0", "CD1"]
    return out


@pytest.fixture
def normalization_exact_adata() -> ad.AnnData:
    return make_normalization_exact_adata()


@pytest.fixture
def synthetic_adata() -> ad.AnnData:
    """A tiny 2-batch cytometry-like AnnData with a planted batch shift.

    - 600 cells, 6 markers + 1 scatter column.
    - 3 hidden "cell types" (simulated by different mean vectors).
    - Batch 2 has a per-marker constant shift on top of batch 1.
    """
    rng = np.random.default_rng(0)
    n_per_type = 100
    n_types = 3
    n_markers = 6
    means = rng.normal(0, 1, size=(n_types, n_markers))

    blocks_b1 = []
    blocks_b2 = []
    celltypes_b1 = []
    celltypes_b2 = []
    for t in range(n_types):
        blocks_b1.append(rng.normal(means[t], 0.3, size=(n_per_type, n_markers)))
        blocks_b2.append(rng.normal(means[t] + 1.5, 0.3, size=(n_per_type, n_markers)))
        celltypes_b1 += [f"type{t}"] * n_per_type
        celltypes_b2 += [f"type{t}"] * n_per_type

    X = np.vstack(blocks_b1 + blocks_b2)
    # Add a non-marker "FSC" column to exercise marker filtering.
    X = np.hstack([X, rng.uniform(100, 200, size=(X.shape[0], 1))])

    var_names = [f"CD{i}" for i in range(n_markers)] + ["FSC"]
    obs = pd.DataFrame(
        {
            "batch": ["b1"] * (n_per_type * n_types) + ["b2"] * (n_per_type * n_types),
            "sample": [f"s{i}" for i in range(X.shape[0])],
            "celltype": celltypes_b1 + celltypes_b2,
        }
    )
    obs.index = obs.index.astype(str)

    adata = ad.AnnData(X=X.astype(float), obs=obs)
    adata.var_names = var_names
    return adata
