import numpy as np

from cycombinepy import batch_correct
from cycombinepy.correct import CORRECTED_LAYER
from cycombinepy.evaluate import (
    compute_emd,
    compute_mad,
    evaluate_emd,
    evaluate_mad,
)


def test_compute_emd_returns_expected_shape(synthetic_adata):
    batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0)
    df = compute_emd(synthetic_adata, cell_key="cycombine_som")
    assert set(df.columns) == {"cluster", "marker", "batch1", "batch2", "emd"}
    assert (df["emd"] >= 0).all()


def test_evaluate_emd_reports_reduction(synthetic_adata):
    batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0)
    uncorr = compute_emd(synthetic_adata, cell_key="cycombine_som")
    corr = compute_emd(
        synthetic_adata, cell_key="cycombine_som", layer=CORRECTED_LAYER
    )
    merged = evaluate_emd(uncorr, corr)
    # On average, the corrected EMD should be smaller than uncorrected.
    assert merged["reduction"].mean() > 0


def test_compute_and_evaluate_mad(synthetic_adata):
    batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0)
    uncorr = compute_mad(synthetic_adata, cell_key="cycombine_som")
    corr = compute_mad(
        synthetic_adata, cell_key="cycombine_som", layer=CORRECTED_LAYER
    )
    merged = evaluate_mad(uncorr, corr)
    assert "reduction" in merged.columns
    assert np.isfinite(merged["mad_corrected"]).any()
