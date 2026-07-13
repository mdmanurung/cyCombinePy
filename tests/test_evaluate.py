import numpy as np

from cycombinepy.correct import CORRECTED_LAYER
from cycombinepy.evaluate import (
    compute_emd,
    compute_mad,
    evaluate_emd,
    evaluate_mad,
)


def _label_by_celltype(adata):
    adata.obs["cycombine_som"] = adata.obs["celltype"].values


def _add_zero_corrected_marker_layer(adata):
    corrected = adata.X.copy()
    markers = [v for v in adata.var_names if v != "FSC"]
    idx = [adata.var_names.get_loc(m) for m in markers]
    corrected[:, idx] = 0.0
    adata.layers[CORRECTED_LAYER] = corrected


def test_compute_emd_returns_expected_shape(synthetic_adata):
    _label_by_celltype(synthetic_adata)
    df = compute_emd(synthetic_adata, cell_key="cycombine_som")
    assert set(df.columns) == {"cluster", "marker", "batch1", "batch2", "emd"}
    assert (df["emd"] >= 0).all()


def test_evaluate_emd_reports_reduction(synthetic_adata):
    _label_by_celltype(synthetic_adata)
    _add_zero_corrected_marker_layer(synthetic_adata)
    uncorr = compute_emd(synthetic_adata, cell_key="cycombine_som")
    corr = compute_emd(
        synthetic_adata, cell_key="cycombine_som", layer=CORRECTED_LAYER
    )
    merged = evaluate_emd(uncorr, corr)
    # On average, the corrected EMD should be smaller than uncorrected.
    assert merged["reduction"].mean() > 0


def test_compute_and_evaluate_mad(synthetic_adata):
    _label_by_celltype(synthetic_adata)
    _add_zero_corrected_marker_layer(synthetic_adata)
    uncorr = compute_mad(synthetic_adata, cell_key="cycombine_som")
    corr = compute_mad(
        synthetic_adata, cell_key="cycombine_som", layer=CORRECTED_LAYER
    )
    merged = evaluate_mad(uncorr, corr)
    assert "reduction" in merged.columns
    assert np.isfinite(merged["mad_corrected"]).any()
