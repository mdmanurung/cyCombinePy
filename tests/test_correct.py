import json
import sys
import types

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

from cycombinepy import (
    CombatCorrectionError,
    ConfoundedDesignError,
    correct_data,
    create_som,
    normalize,
)
from cycombinepy.correct import CORRECTED_LAYER


def _batch_mean_gap(X, batch):
    b1 = X[batch == "b1"].mean(axis=0)
    b2 = X[batch == "b2"].mean(axis=0)
    return float(np.linalg.norm(b1 - b2))


def _cluster_on_normalized(adata):
    """Normalize a copy and copy the cluster labels back to ``adata``.

    Mirrors the intended cyCombine workflow: cluster on batch-normalized data
    (so clusters represent cell populations, not batches), then apply ComBat to
    the *unnormalized* data per cluster.
    """
    tmp = adata.copy()
    normalize(tmp, method="scale")
    create_som(tmp, xdim=3, ydim=3, seed=0)
    adata.obs["cycombine_som"] = tmp.obs["cycombine_som"].values


def _fixed_som_adata():
    X = np.array(
        [
            [1.0, 10.0],
            [2.0, 20.0],
            [3.0, 30.0],
            [4.0, 40.0],
            [5.0, 50.0],
            [6.0, 60.0],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(
        {
            "batch": ["b1", "b1", "b2", "b1", "b1", "b1"],
            "cycombine_som": ["multi", "multi", "multi", "single", "single", "single"],
        },
        index=[f"cell{i}" for i in range(X.shape[0])],
    )
    adata = AnnData(X=X, obs=obs)
    adata.var_names = ["CD3", "CD4"]
    return adata


def _fixed_som_adata_with_confounded_covar():
    adata = _fixed_som_adata()
    adata.obs["celltype"] = ["t", "t", "b", "t", "t", "t"]
    return adata


def _fixed_som_adata_with_skewed_covar():
    adata = _fixed_som_adata()
    adata.obs["celltype"] = ["t", "b", "t", "t", "t", "t"]
    return adata


def _cluster_rows(report):
    clusters = report["clusters"]
    return {
        str(label): {
            "label": str(label),
            "n_cells": int(n_cells),
            "batches": json.loads(str(batches)),
            "status": str(status),
            "terms": json.loads(str(terms)),
            "exception": json.loads(str(exception)),
        }
        for label, n_cells, batches, status, terms, exception in zip(
            clusters["label"],
            clusters["n_cells"],
            clusters["batches"],
            clusters["status"],
            clusters["terms"],
            clusters["exception"],
        )
    }


def test_correct_data_writes_success_report_to_uns(monkeypatch):
    adata = _fixed_som_adata()
    calls = []

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        calls.append(
            {
                "shape": x.shape,
                "batch": list(batch),
                "mod": mod,
                "parametric": parametric,
                "ref_batch": ref_batch,
            }
        )
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    result = correct_data(adata, label_key="cycombine_som", markers=["CD3", "CD4"])

    assert result is None
    assert calls == [
        {
            "shape": (2, 3),
            "batch": ["b1", "b1", "b2"],
            "mod": None,
            "parametric": True,
            "ref_batch": None,
        }
    ]
    report = adata.uns["cycombinepy_correction"]
    assert report["schema_version"] == "1.0"
    assert report["function"] == "correct_data"
    assert report["status"] == "completed"
    assert report["output_written"] is True
    assert report["timestamp_utc"].endswith("Z")
    assert report["package"]["name"] == "cycombinepy"
    assert isinstance(report["package"]["version"], str)
    assert report["parameters"] == {
        "batch_key": "batch",
        "label_key": "cycombine_som",
        "covar": None,
        "anchor": None,
        "parametric": True,
        "ref_batch": None,
        "out_layer": CORRECTED_LAYER,
        "error_policy": "raise",
        "confound_policy": "raise",
    }
    assert report["markers"] == ["CD3", "CD4"]

    clusters = _cluster_rows(report)
    assert clusters["multi"] == {
        "label": "multi",
        "n_cells": 3,
        "batches": ["b1", "b2"],
        "status": "corrected",
        "terms": [],
        "exception": None,
    }
    assert clusters["single"] == {
        "label": "single",
        "n_cells": 3,
        "batches": ["b1"],
        "status": "skipped_single_batch",
        "terms": [],
        "exception": None,
    }


def test_correct_data_return_report_matches_uns(monkeypatch):
    adata = _fixed_som_adata()

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    report = correct_data(
        adata,
        label_key="cycombine_som",
        markers=["CD3", "CD4"],
        return_report=True,
    )

    assert report == adata.uns["cycombinepy_correction"]
    assert report["status"] == "completed"
    assert report["output_written"] is True


def test_correct_data_copy_return_report_returns_copy_and_report(monkeypatch):
    adata = _fixed_som_adata()

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    result = correct_data(
        adata,
        label_key="cycombine_som",
        markers=["CD3", "CD4"],
        copy=True,
        return_report=True,
    )

    assert "cycombinepy_correction" not in adata.uns
    assert CORRECTED_LAYER not in adata.layers
    assert isinstance(result, tuple)
    out, report = result
    assert out is not adata
    assert report == out.uns["cycombinepy_correction"]
    assert CORRECTED_LAYER in out.layers


def test_correct_data_report_survives_h5ad_roundtrip(monkeypatch, tmp_path):
    adata = _fixed_som_adata()

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    correct_data(adata, label_key="cycombine_som", markers=["CD3", "CD4"])
    path = tmp_path / "corrected.h5ad"

    adata.write_h5ad(path)
    restored = read_h5ad(path)

    assert "cycombinepy_correction" in restored.uns
    report = restored.uns["cycombinepy_correction"]
    assert report["schema_version"] == "1.0"
    assert report["function"] == "correct_data"
    assert report["status"] == "completed"
    assert bool(report["output_written"]) is True
    assert list(report["markers"]) == ["CD3", "CD4"]
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "corrected"
    assert clusters["multi"]["batches"] == ["b1", "b2"]
    assert clusters["single"]["status"] == "skipped_single_batch"
    assert clusters["single"]["batches"] == ["b1"]


def test_correct_data_combat_failure_default_raises_and_writes_failed_report(monkeypatch):
    adata = _fixed_som_adata()

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        raise RuntimeError("singular design")

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    with pytest.raises(CombatCorrectionError, match="ComBat failed.*multi") as exc_info:
        correct_data(adata, label_key="cycombine_som", markers=["CD3", "CD4"])

    assert CORRECTED_LAYER not in adata.layers
    report = adata.uns["cycombinepy_correction"]
    assert exc_info.value.report is report
    assert report["status"] == "failed"
    assert report["output_written"] is False
    assert report["parameters"]["error_policy"] == "raise"
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "failed"
    assert clusters["multi"]["exception"] == {
        "type": "RuntimeError",
        "message": "singular design",
    }


def test_correct_data_combat_failure_report_policy_records_and_writes_layer(monkeypatch):
    adata = _fixed_som_adata()

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        raise RuntimeError("singular design")

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    correct_data(
        adata,
        label_key="cycombine_som",
        markers=["CD3", "CD4"],
        error_policy="report",
    )

    assert CORRECTED_LAYER in adata.layers
    np.testing.assert_array_equal(adata.layers[CORRECTED_LAYER], adata.X)
    report = adata.uns["cycombinepy_correction"]
    assert report["status"] == "completed_with_failures"
    assert report["output_written"] is True
    assert report["parameters"]["error_policy"] == "report"
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "failed"
    assert clusters["multi"]["exception"]["type"] == "RuntimeError"
    assert clusters["single"]["status"] == "skipped_single_batch"


def test_correct_data_combat_failure_warn_policy_emits_warning_and_writes_layer(monkeypatch):
    adata = _fixed_som_adata()

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        raise RuntimeError("singular design")

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    with pytest.warns(RuntimeWarning, match="ComBat failed.*multi"):
        correct_data(
            adata,
            label_key="cycombine_som",
            markers=["CD3", "CD4"],
            error_policy="warn",
        )

    assert CORRECTED_LAYER in adata.layers
    report = adata.uns["cycombinepy_correction"]
    assert report["status"] == "completed_with_failures"
    assert report["parameters"]["error_policy"] == "warn"
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "failed"
    assert clusters["multi"]["exception"]["message"] == "singular design"


def test_correct_data_confounded_covar_default_raises_before_combat(monkeypatch):
    adata = _fixed_som_adata_with_confounded_covar()
    calls = []

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        calls.append(mod)
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    with pytest.raises(
        ConfoundedDesignError,
        match="Confounded design.*multi.*celltype",
    ) as exc_info:
        correct_data(
            adata,
            label_key="cycombine_som",
            markers=["CD3", "CD4"],
            covar="celltype",
        )

    assert calls == []
    assert CORRECTED_LAYER not in adata.layers
    report = adata.uns["cycombinepy_correction"]
    assert exc_info.value.report is report
    assert report["status"] == "failed"
    assert report["output_written"] is False
    assert report["parameters"]["confound_policy"] == "raise"
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "failed"
    assert clusters["multi"]["exception"]["type"] == "ConfoundedDesignError"
    assert clusters["multi"]["exception"]["dropped_terms"] == [
        {
            "term": "celltype",
            "role": "covar",
            "reason": "confounded_with_batch",
        }
    ]


def test_correct_data_skewed_covar_default_drops_and_reports_adjustment(monkeypatch):
    adata = _fixed_som_adata_with_skewed_covar()
    calls = []

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        calls.append(mod)
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    correct_data(
        adata,
        label_key="cycombine_som",
        markers=["CD3", "CD4"],
        covar="celltype",
    )

    assert len(calls) == 1
    assert calls[0] is None
    assert CORRECTED_LAYER in adata.layers
    report = adata.uns["cycombinepy_correction"]
    assert report["status"] == "completed_with_adjustments"
    assert report["output_written"] is True
    assert report["parameters"]["confound_policy"] == "raise"
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "corrected"
    assert clusters["multi"]["terms"] == []
    assert clusters["multi"]["exception"]["type"] == "ConfoundedDesignAdjustment"
    assert clusters["multi"]["exception"]["dropped_terms"] == [
        {
            "term": "celltype",
            "role": "covar",
            "reason": "skewed_to_single_level",
        }
    ]


def test_correct_data_reraises_transitive_model_matrix_import_failure(monkeypatch):
    adata = _fixed_som_adata_with_skewed_covar()
    fake_formulaic = types.ModuleType("formulaic")

    def fake_model_matrix(formula, sub):
        raise ModuleNotFoundError(
            "No module named 'formulaic_optional'",
            name="formulaic_optional",
        )

    fake_formulaic.model_matrix = fake_model_matrix
    monkeypatch.setitem(sys.modules, "formulaic", fake_formulaic)

    with pytest.raises(ModuleNotFoundError) as exc_info:
        correct_data(
            adata,
            label_key="cycombine_som",
            markers=["CD3", "CD4"],
            covar="celltype",
        )

    assert exc_info.value.name == "formulaic_optional"
    assert CORRECTED_LAYER not in adata.layers


def test_correct_data_confounded_covar_skip_policy_marks_cluster_and_writes_layer(monkeypatch):
    adata = _fixed_som_adata_with_confounded_covar()
    calls = []

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        calls.append(mod)
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    correct_data(
        adata,
        label_key="cycombine_som",
        markers=["CD3", "CD4"],
        covar="celltype",
        confound_policy="skip",
    )

    assert calls == []
    assert CORRECTED_LAYER in adata.layers
    np.testing.assert_array_equal(adata.layers[CORRECTED_LAYER], adata.X)
    report = adata.uns["cycombinepy_correction"]
    assert report["status"] == "completed_with_adjustments"
    assert report["output_written"] is True
    assert report["parameters"]["confound_policy"] == "skip"
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "skipped_confounded_design"
    assert clusters["multi"]["exception"]["type"] == "ConfoundedDesign"
    assert clusters["multi"]["exception"]["dropped_terms"][0]["term"] == "celltype"
    assert clusters["single"]["status"] == "skipped_single_batch"


def test_correct_data_confounded_covar_drop_policy_corrects_with_audited_drop(monkeypatch):
    adata = _fixed_som_adata_with_confounded_covar()
    calls = []

    def fake_run_combat(x, batch, mod, parametric, ref_batch):
        calls.append(mod)
        return x

    monkeypatch.setattr("cycombinepy.correct.run_combat", fake_run_combat)

    correct_data(
        adata,
        label_key="cycombine_som",
        markers=["CD3", "CD4"],
        covar="celltype",
        confound_policy="drop",
    )

    assert len(calls) == 1
    assert calls[0] is None
    assert CORRECTED_LAYER in adata.layers
    report = adata.uns["cycombinepy_correction"]
    assert report["status"] == "completed_with_adjustments"
    assert report["parameters"]["confound_policy"] == "drop"
    clusters = _cluster_rows(report)
    assert clusters["multi"]["status"] == "corrected"
    assert clusters["multi"]["terms"] == []
    assert clusters["multi"]["exception"]["type"] == "ConfoundedDesignAdjustment"
    assert clusters["multi"]["exception"]["dropped_terms"] == [
        {
            "term": "celltype",
            "role": "covar",
            "reason": "confounded_with_batch",
        }
    ]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"error_policy": "ignore"}, "error_policy"),
        ({"confound_policy": "ignore"}, "confound_policy"),
    ],
)
def test_correct_data_invalid_policy_values_raise_value_error(kwargs, message):
    adata = _fixed_som_adata()

    with pytest.raises(ValueError, match=message):
        correct_data(
            adata,
            label_key="cycombine_som",
            markers=["CD3", "CD4"],
            **kwargs,
        )


@pytest.mark.requires_flowsom
@pytest.mark.requires_inmoose
def test_correct_data_reduces_batch_gap(synthetic_adata):
    pytest.importorskip("flowsom")
    pytest.importorskip("inmoose")
    _cluster_on_normalized(synthetic_adata)
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]
    idx = [synthetic_adata.var_names.get_loc(m) for m in markers]
    batch = synthetic_adata.obs["batch"].to_numpy()

    before = _batch_mean_gap(synthetic_adata.X[:, idx], batch)

    correct_data(synthetic_adata, label_key="cycombine_som", markers=markers)
    corrected = synthetic_adata.layers[CORRECTED_LAYER][:, idx]

    after = _batch_mean_gap(corrected, batch)
    assert after < before * 0.5, f"batch gap did not shrink (before={before}, after={after})"


@pytest.mark.requires_flowsom
@pytest.mark.requires_inmoose
def test_correct_data_with_covar(synthetic_adata):
    pytest.importorskip("flowsom")
    pytest.importorskip("inmoose")
    _cluster_on_normalized(synthetic_adata)
    correct_data(
        synthetic_adata,
        label_key="cycombine_som",
        covar="celltype",
    )
    assert CORRECTED_LAYER in synthetic_adata.layers
