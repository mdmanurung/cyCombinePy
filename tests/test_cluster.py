import numpy as np

from cycombinepy import create_som


def test_create_som_writes_labels(synthetic_adata):
    create_som(synthetic_adata, xdim=3, ydim=3, seed=0)
    assert "cycombine_som" in synthetic_adata.obs.columns
    labels = synthetic_adata.obs["cycombine_som"]
    # Categorical, at least 2 distinct labels, no NaNs
    assert labels.dtype.name == "category"
    assert labels.nunique() >= 2
    assert not labels.isna().any()


def test_create_som_respects_markers(synthetic_adata):
    # Passing an explicit subset should still produce labels for every cell.
    create_som(
        synthetic_adata,
        markers=["CD0", "CD1", "CD2"],
        xdim=3,
        ydim=3,
        seed=0,
    )
    assert len(synthetic_adata.obs["cycombine_som"]) == synthetic_adata.n_obs
