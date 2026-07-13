import numpy as np

from cycombinepy import transform_asinh


def _positive_marker_copy(adata):
    out = adata.copy()
    markers = [v for v in out.var_names if v != "FSC"]
    idx = [out.var_names.get_loc(marker) for marker in markers]
    out.X[:, idx] = np.abs(out.X[:, idx]) + 1.0
    return out, markers, idx


def test_transform_asinh_derand_same_seed_exact(synthetic_adata):
    adata1, markers, idx = _positive_marker_copy(synthetic_adata)
    adata2 = adata1.copy()

    transform_asinh(adata1, markers=markers, cofactor=5, derand=True, seed=42)
    transform_asinh(adata2, markers=markers, cofactor=5, derand=True, seed=42)

    np.testing.assert_array_equal(adata1.X[:, idx], adata2.X[:, idx])


def test_transform_asinh_derand_different_seed_changes_values(synthetic_adata):
    adata1, markers, idx = _positive_marker_copy(synthetic_adata)
    adata2 = adata1.copy()

    transform_asinh(adata1, markers=markers, cofactor=5, derand=True, seed=42)
    transform_asinh(adata2, markers=markers, cofactor=5, derand=True, seed=43)

    assert np.all(np.isfinite(adata1.X[:, idx]))
    assert np.all(np.isfinite(adata2.X[:, idx]))
    assert not np.allclose(adata1.X[:, idx], adata2.X[:, idx])
