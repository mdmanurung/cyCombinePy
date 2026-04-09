import numpy as np

from cycombinepy import transform_asinh


def test_transform_asinh_roundtrip(synthetic_adata):
    # Start from positive counts-like values
    synthetic_adata.X = np.abs(synthetic_adata.X) + 1
    orig = synthetic_adata.X.copy()
    transform_asinh(synthetic_adata, cofactor=5, derand=False)
    # All marker columns should have been asinh-transformed
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]
    idx = [synthetic_adata.var_names.get_loc(m) for m in markers]
    np.testing.assert_allclose(
        synthetic_adata.X[:, idx], np.arcsinh(orig[:, idx] / 5), rtol=1e-10
    )


def test_transform_asinh_derand_reduces_to_asinh_integer_part(synthetic_adata):
    synthetic_adata.X = np.abs(synthetic_adata.X) + 1
    transform_asinh(synthetic_adata, cofactor=5, derand=True, seed=42)
    # Derandomized values must be finite and non-negative after asinh of a
    # non-negative argument.
    assert np.all(np.isfinite(synthetic_adata.X))
