import builtins

import anndata as ad
import numpy as np
import pandas as pd
import pytest


def test_normalize_rejects_missing_batch_without_mutating(synthetic_adata):
    from cycombinepy import normalize

    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "batch"] = np.nan
    before = synthetic_adata.X.copy()

    with pytest.raises(ValueError, match=r"normalize\(\) requires non-missing values"):
        normalize(synthetic_adata, method="scale")

    np.testing.assert_array_equal(synthetic_adata.X, before)


def test_correct_data_rejects_missing_label_before_string_coercion(
    synthetic_adata, monkeypatch
):
    from cycombinepy import correct_data

    synthetic_adata.obs["cycombine_som"] = ["1"] * synthetic_adata.n_obs
    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "cycombine_som"] = np.nan
    monkeypatch.setattr("cycombinepy.correct.run_combat", lambda x, **kwargs: x)

    with pytest.raises(
        ValueError,
        match=r"correct_data\(\) requires non-missing values",
    ):
        correct_data(synthetic_adata, label_key="cycombine_som")


def test_correct_data_rejects_missing_covar_before_string_coercion(
    synthetic_adata, monkeypatch
):
    from cycombinepy import correct_data

    synthetic_adata.obs["cycombine_som"] = ["1"] * synthetic_adata.n_obs
    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "celltype"] = pd.NA
    monkeypatch.setattr("cycombinepy.correct.run_combat", lambda x, **kwargs: x)

    with pytest.raises(
        ValueError,
        match=r"correct_data\(\) requires non-missing values.*'celltype'",
    ):
        correct_data(synthetic_adata, label_key="cycombine_som", covar="celltype")


def test_correct_data_rejects_missing_anchor_before_string_coercion(
    synthetic_adata, monkeypatch
):
    from cycombinepy import correct_data

    synthetic_adata.obs["cycombine_som"] = ["1"] * synthetic_adata.n_obs
    synthetic_adata.obs["anchor"] = ["a"] * synthetic_adata.n_obs
    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "anchor"] = pd.NA
    monkeypatch.setattr("cycombinepy.correct.run_combat", lambda x, **kwargs: x)

    with pytest.raises(
        ValueError,
        match=r"correct_data\(\) requires non-missing values.*'anchor'",
    ):
        correct_data(synthetic_adata, label_key="cycombine_som", anchor="anchor")


def test_correct_data_rejects_missing_batch_before_string_coercion(
    synthetic_adata, monkeypatch
):
    from cycombinepy import correct_data

    synthetic_adata.obs["cycombine_som"] = ["1"] * synthetic_adata.n_obs
    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "batch"] = pd.NA
    monkeypatch.setattr("cycombinepy.correct.run_combat", lambda x, **kwargs: x)

    with pytest.raises(
        ValueError,
        match=r"correct_data\(\) requires non-missing values.*'batch'",
    ):
        correct_data(synthetic_adata, label_key="cycombine_som")


def test_compute_emd_rejects_non_finite_marker_matrix(synthetic_adata):
    from cycombinepy.evaluate import compute_emd

    synthetic_adata.obs["cycombine_som"] = ["1"] * synthetic_adata.n_obs
    synthetic_adata.X[0, 0] = np.inf

    with pytest.raises(
        ValueError,
        match=r"compute_emd\(\) requires finite marker values",
    ):
        compute_emd(synthetic_adata, cell_key="cycombine_som")


def test_compute_mad_rejects_non_finite_marker_matrix(synthetic_adata):
    from cycombinepy.evaluate import compute_mad

    synthetic_adata.obs["cycombine_som"] = ["1"] * synthetic_adata.n_obs
    synthetic_adata.X[0, 0] = np.nan

    with pytest.raises(
        ValueError,
        match=r"compute_mad\(\) requires finite marker values",
    ):
        compute_mad(synthetic_adata, cell_key="cycombine_som")


def test_compute_mad_rejects_missing_batch_values(synthetic_adata):
    from cycombinepy.evaluate import compute_mad

    synthetic_adata.obs["cycombine_som"] = ["1"] * synthetic_adata.n_obs
    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "batch"] = pd.NA

    with pytest.raises(
        ValueError,
        match=r"compute_mad\(\) requires non-missing values.*'batch'",
    ):
        compute_mad(synthetic_adata, cell_key="cycombine_som")


def test_scib_metrics_rejects_missing_batch_before_optional_import(
    synthetic_adata, monkeypatch
):
    from cycombinepy.evaluate import scib_metrics

    real_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name.startswith(("scanpy", "scib_metrics")):
            raise AssertionError("optional import reached before validation")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    with pytest.raises(KeyError, match='Column "missing_batch"'):
        scib_metrics(synthetic_adata, batch_key="missing_batch")


def test_scib_metrics_rejects_non_finite_layer_before_optional_import(
    synthetic_adata, monkeypatch
):
    from cycombinepy.evaluate import scib_metrics

    synthetic_adata.layers["bad"] = synthetic_adata.X.copy()
    synthetic_adata.layers["bad"][0, 0] = np.inf
    real_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name.startswith(("scanpy", "scib_metrics")):
            raise AssertionError("optional import reached before validation")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    with pytest.raises(
        ValueError,
        match=r"scib_metrics\(\) requires finite marker values",
    ):
        scib_metrics(synthetic_adata, batch_key="batch", layer="bad")


def test_detect_batch_effect_express_rejects_missing_batch(synthetic_adata):
    from cycombinepy.detect import detect_batch_effect_express

    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "batch"] = pd.NA

    with pytest.raises(
        ValueError,
        match=r"detect_batch_effect_express\(\) requires non-missing values.*'batch'",
    ):
        detect_batch_effect_express(synthetic_adata)


def test_detect_batch_effect_express_rejects_missing_existing_sample_key(
    synthetic_adata,
):
    from cycombinepy.detect import detect_batch_effect_express

    synthetic_adata.obs.loc[synthetic_adata.obs.index[0], "sample"] = pd.NA

    with pytest.raises(
        ValueError,
        match=r"detect_batch_effect_express\(\) requires non-missing values.*'sample'",
    ):
        detect_batch_effect_express(synthetic_adata)


def test_detect_batch_effect_rejects_missing_batch_before_optional_import(
    synthetic_adata,
    monkeypatch,
):
    from cycombinepy.detect import detect_batch_effect

    real_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        if name.startswith(("matplotlib", "scanpy")):
            raise AssertionError("optional import reached before validation")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    with pytest.raises(KeyError, match='Column "missing_batch"'):
        detect_batch_effect(synthetic_adata, batch_key="missing_batch")


def test_read_fcs_dir_raises_for_missing_requested_markers(tmp_path, monkeypatch):
    from cycombinepy.io import read_fcs_dir

    fcs = tmp_path / "sample.fcs"
    fcs.write_bytes(b"fake-fcs")

    def fake_read(_path):
        a = ad.AnnData(X=np.ones((3, 2)), obs=pd.DataFrame(index=["0", "1", "2"]))
        a.var_names = ["CD0", "CD1"]
        return a

    monkeypatch.setattr("cycombinepy.io._read_fcs_one", fake_read)

    with pytest.raises(KeyError, match="Markers not found"):
        read_fcs_dir(tmp_path, markers=["CD0", "missing"], transform=False)
