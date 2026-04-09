import numpy as np

from pycombine import normalize


def test_normalize_scale_per_batch_has_zero_mean(synthetic_adata):
    normalize(synthetic_adata, method="scale")
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]
    idx = [synthetic_adata.var_names.get_loc(m) for m in markers]
    for b in ("b1", "b2"):
        mask = (synthetic_adata.obs["batch"] == b).to_numpy()
        block = synthetic_adata.X[mask][:, idx]
        np.testing.assert_allclose(block.mean(axis=0), 0, atol=1e-10)


def test_normalize_rank_maps_to_percentiles(synthetic_adata):
    normalize(synthetic_adata, method="rank")
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]
    idx = [synthetic_adata.var_names.get_loc(m) for m in markers]
    for b in ("b1", "b2"):
        mask = (synthetic_adata.obs["batch"] == b).to_numpy()
        block = synthetic_adata.X[mask][:, idx]
        assert block.min() > 0.0
        assert block.max() <= 1.0 + 1e-12


def test_normalize_none_is_noop(synthetic_adata):
    before = synthetic_adata.X.copy()
    normalize(synthetic_adata, method="none")
    np.testing.assert_array_equal(synthetic_adata.X, before)


def test_normalize_qnorm_runs(synthetic_adata):
    normalize(synthetic_adata, method="qnorm")
    assert np.all(np.isfinite(synthetic_adata.X))
