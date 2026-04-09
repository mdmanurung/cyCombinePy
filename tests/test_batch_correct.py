import numpy as np

from cycombinepy import batch_correct
from cycombinepy.correct import CORRECTED_LAYER


def _gap(X, batch):
    return float(np.linalg.norm(X[batch == "b1"].mean(0) - X[batch == "b2"].mean(0)))


def test_batch_correct_full_pipeline_reduces_gap(synthetic_adata):
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]
    idx = [synthetic_adata.var_names.get_loc(m) for m in markers]
    batch = synthetic_adata.obs["batch"].to_numpy()

    before = _gap(synthetic_adata.X[:, idx], batch)

    batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0)

    assert CORRECTED_LAYER in synthetic_adata.layers
    after = _gap(synthetic_adata.layers[CORRECTED_LAYER][:, idx], batch)
    assert after < before * 0.5, f"before={before} after={after}"


def test_batch_correct_iterative_grids(synthetic_adata):
    batch_correct(synthetic_adata, xdim=[3, 4], ydim=[3, 4], seed=0)
    assert CORRECTED_LAYER in synthetic_adata.layers


def test_batch_correct_copy_returns_new(synthetic_adata):
    out = batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0, copy=True)
    assert out is not None
    assert CORRECTED_LAYER in out.layers
    assert CORRECTED_LAYER not in synthetic_adata.layers
