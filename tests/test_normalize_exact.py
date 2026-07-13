import numpy as np

from cycombinepy import normalize


def test_normalize_rank_average_ties_exact_per_batch(normalization_exact_adata):
    normalize(normalization_exact_adata, method="rank", ties_method="average")

    cd0 = normalization_exact_adata[:, "CD0"].X.ravel()
    np.testing.assert_allclose(cd0[:3], [0.5, 0.5, 1.0])
    np.testing.assert_allclose(cd0[3:], [0.5, 0.5, 1.0])


def test_normalize_scale_constant_column_stays_finite(normalization_exact_adata):
    normalization_exact_adata.X[:3, 1] = 7.0

    normalize(normalization_exact_adata, method="scale")

    assert np.isfinite(normalization_exact_adata.X).all()


def test_normalize_copy_leaves_input_unchanged(normalization_exact_adata):
    before = normalization_exact_adata.X.copy()

    normalized = normalize(normalization_exact_adata, method="rank", copy=True)

    np.testing.assert_array_equal(normalization_exact_adata.X, before)
    assert normalized is not normalization_exact_adata
    np.testing.assert_allclose(normalized[:, "CD0"].X.ravel(), [0.5, 0.5, 1.0] * 2)
