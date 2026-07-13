import numpy as np
import pandas as pd
import pytest

from cycombinepy import check_confound, get_markers
from cycombinepy._utils import (
    check_obs_values_not_missing,
    marker_matrix,
    resolve_markers,
)


def test_get_markers_excludes_non_markers(synthetic_adata):
    markers = get_markers(synthetic_adata)
    assert "FSC" not in markers
    assert all(m.startswith("CD") for m in markers)
    assert len(markers) == 6


def test_resolve_markers_rejects_duplicate_var_names(synthetic_adata):
    adata = synthetic_adata.copy()
    adata.var_names = ["CD1", "CD1", "CD2", "CD3", "CD4", "CD5", "FSC"]

    with pytest.raises(ValueError, match=r"adata\.var_names must be unique.*CD1"):
        resolve_markers(adata, ["CD1"])


def test_marker_matrix_raises_clear_error_for_missing_layer(synthetic_adata):
    with pytest.raises(
        KeyError,
        match=r"Layer 'missing' was not found in adata\.layers",
    ):
        marker_matrix(synthetic_adata, ["CD0"], layer="missing")


def test_marker_matrix_rejects_non_finite_values_when_required(synthetic_adata):
    adata = synthetic_adata.copy()
    adata.X[0, 0] = np.nan

    with pytest.raises(
        ValueError,
        match=r"normalize\(\).*requires finite marker values",
    ):
        marker_matrix(adata, ["CD0"], require_finite=True, context="normalize()")


def test_check_obs_values_not_missing_allows_literal_nan_string(synthetic_adata):
    adata = synthetic_adata.copy()
    adata.obs["label"] = "nan"

    values = check_obs_values_not_missing(adata, "label", context="unit()")

    assert values.eq("nan").all()


def test_check_obs_values_not_missing_rejects_actual_missing_values(synthetic_adata):
    adata = synthetic_adata.copy()
    adata.obs["label"] = "ok"
    adata.obs.loc[adata.obs.index[0], "label"] = pd.NA

    with pytest.raises(ValueError, match=r"unit\(\).*requires non-missing values"):
        check_obs_values_not_missing(adata, "label", context="unit()")


def test_check_obs_values_not_missing_caps_reported_indices(synthetic_adata):
    adata = synthetic_adata[:8].copy()
    adata.obs_names = [f"row{i}" for i in range(adata.n_obs)]
    adata.obs["label"] = "ok"
    adata.obs.loc[adata.obs_names[:7], "label"] = pd.NA

    with pytest.raises(ValueError) as excinfo:
        check_obs_values_not_missing(adata, "label", context="unit()")

    message = str(excinfo.value)
    assert "found 7 missing value(s)" in message
    assert "row0" in message
    assert "row4" in message
    assert "row5" not in message
    assert "row6" not in message


def test_check_confound_true_when_covar_matches_batch():
    batch = np.array(["b1"] * 10 + ["b2"] * 10)
    # Perfectly confounded binary covariate
    cov = np.array([0] * 10 + [1] * 10).reshape(-1, 1).astype(float)
    assert check_confound(batch, cov) is True


def test_check_confound_false_when_covar_independent():
    batch = np.array(["b1", "b2"] * 10)
    cov = np.array([0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1])
    cov = cov.reshape(-1, 1).astype(float)
    assert check_confound(batch, cov) is False


def test_check_confound_none_returns_false():
    batch = np.array(["b1"] * 5 + ["b2"] * 5)
    assert check_confound(batch, None) is False
