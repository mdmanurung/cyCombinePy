import numpy as np
import pandas as pd

from cycombinepy import check_confound, get_markers


def test_get_markers_excludes_non_markers(synthetic_adata):
    markers = get_markers(synthetic_adata)
    assert "FSC" not in markers
    assert all(m.startswith("CD") for m in markers)
    assert len(markers) == 6


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
