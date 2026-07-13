import importlib

import numpy as np
import pandas as pd
import pytest
from anndata import read_h5ad

from cycombinepy import batch_correct
from cycombinepy.correct import CORRECTED_LAYER, CombatCorrectionError


batch_correct_module = importlib.import_module("cycombinepy.batch_correct")


def _set_som_labels(adata, markers, xdim, ydim, n_clusters, seed, rlen, label_key):
    adata.obs[label_key] = [f"som_{xdim}_{ydim}"] * adata.n_obs


def _fake_correct_report(status="completed", output_written=True):
    return {
        "schema_version": "1.0",
        "function": "correct_data",
        "status": status,
        "output_written": output_written,
        "parameters": {"out_layer": CORRECTED_LAYER},
        "markers": ["CD3", "CD4"],
        "clusters": {
            "label": ["som"],
            "n_cells": [6],
            "batches": ['["b1","b2"]'],
            "status": ["corrected"],
            "terms": ["[]"],
            "exception": ["null"],
        },
    }


def _iteration_reports(report):
    import json

    return [json.loads(item) for item in report["iterations"]]


def _gap(X, batch):
    return float(np.linalg.norm(X[batch == "b1"].mean(0) - X[batch == "b2"].mean(0)))


@pytest.mark.requires_flowsom
@pytest.mark.requires_inmoose
def test_batch_correct_full_pipeline_reduces_gap(synthetic_adata):
    pytest.importorskip("flowsom")
    pytest.importorskip("inmoose")
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]
    idx = [synthetic_adata.var_names.get_loc(m) for m in markers]
    batch = synthetic_adata.obs["batch"].to_numpy()

    before = _gap(synthetic_adata.X[:, idx], batch)

    batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0)

    assert CORRECTED_LAYER in synthetic_adata.layers
    after = _gap(synthetic_adata.layers[CORRECTED_LAYER][:, idx], batch)
    assert after < before * 0.5, f"before={before} after={after}"


@pytest.mark.requires_flowsom
@pytest.mark.requires_inmoose
def test_batch_correct_iterative_grids(synthetic_adata):
    pytest.importorskip("flowsom")
    pytest.importorskip("inmoose")
    batch_correct(synthetic_adata, xdim=[3, 4], ydim=[3, 4], seed=0)
    assert CORRECTED_LAYER in synthetic_adata.layers


@pytest.mark.requires_flowsom
@pytest.mark.requires_inmoose
def test_batch_correct_copy_returns_new(synthetic_adata):
    pytest.importorskip("flowsom")
    pytest.importorskip("inmoose")
    out = batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0, copy=True)
    assert out is not None
    assert CORRECTED_LAYER in out.layers
    assert CORRECTED_LAYER not in synthetic_adata.layers


def test_batch_correct_forwards_error_and_confound_policies(monkeypatch, synthetic_adata):
    calls = []

    def fake_correct_data(adata, **kwargs):
        calls.append(kwargs)
        adata.layers[kwargs["out_layer"]] = adata.X.copy()
        return _fake_correct_report()

    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)
    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    batch_correct(
        synthetic_adata,
        markers=["CD3", "CD4"],
        xdim=3,
        ydim=3,
        error_policy="report",
        confound_policy="drop",
        uns_key="batch_report",
    )

    assert len(calls) == 1
    assert calls[0]["error_policy"] == "report"
    assert calls[0]["confound_policy"] == "drop"
    assert calls[0]["return_report"] is True
    assert calls[0]["uns_key"] != "batch_report"
    assert "batch_report" in synthetic_adata.uns
    assert synthetic_adata.uns["batch_report"]["function"] == "batch_correct"


def test_batch_correct_return_report_shapes(monkeypatch, synthetic_adata):
    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)

    def fake_correct_data(adata, **kwargs):
        adata.layers[kwargs["out_layer"]] = adata.X.copy()
        return _fake_correct_report()

    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    report = batch_correct(
        synthetic_adata,
        markers=["CD3", "CD4"],
        xdim=3,
        ydim=3,
        return_report=True,
    )

    assert report == synthetic_adata.uns["cycombinepy_correction"]
    assert report["function"] == "batch_correct"
    assert len(report["iterations"]) == 1
    assert report["markers"] == ["CD3", "CD4"]
    assert report["output_written"] is True
    assert CORRECTED_LAYER in synthetic_adata.layers


def test_batch_correct_copy_return_report_returns_copy_and_report(monkeypatch, synthetic_adata):
    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)

    def fake_correct_data(adata, **kwargs):
        adata.layers[kwargs["out_layer"]] = adata.X.copy()
        return _fake_correct_report()

    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    out, report = batch_correct(
        synthetic_adata,
        markers=["CD3", "CD4"],
        xdim=3,
        ydim=3,
        copy=True,
        return_report=True,
    )

    assert out is not synthetic_adata
    assert CORRECTED_LAYER not in synthetic_adata.layers
    assert "cycombinepy_correction" not in synthetic_adata.uns
    assert CORRECTED_LAYER in out.layers
    assert report == out.uns["cycombinepy_correction"]


def test_batch_correct_iterative_report_aggregates_status(monkeypatch, synthetic_adata):
    reports = [
        _fake_correct_report(status="completed"),
        _fake_correct_report(status="completed_with_adjustments"),
    ]

    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)

    def fake_correct_data(adata, **kwargs):
        adata.layers[kwargs["out_layer"]] = adata.X.copy()
        return reports.pop(0)

    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    report = batch_correct(
        synthetic_adata,
        markers=["CD3", "CD4"],
        xdim=[3, 4],
        ydim=[3, 4],
        return_report=True,
    )

    assert report["status"] == "completed_with_adjustments"
    assert len(report["iterations"]) == 2
    assert [item["status"] for item in _iteration_reports(report)] == [
        "completed",
        "completed_with_adjustments",
    ]


def test_batch_correct_fatal_correction_writes_failed_aggregate_report(
    monkeypatch, synthetic_adata
):
    exception_report = _fake_correct_report(status="failed", output_written=False)

    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)

    def fake_correct_data(adata, **kwargs):
        raise CombatCorrectionError("ComBat failed for cluster 'som'", exception_report)

    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    with pytest.raises(CombatCorrectionError, match="ComBat failed"):
        batch_correct(
            synthetic_adata,
            markers=["CD3", "CD4"],
            xdim=3,
            ydim=3,
            return_report=True,
        )

    assert CORRECTED_LAYER not in synthetic_adata.layers
    report = synthetic_adata.uns["cycombinepy_correction"]
    assert report["function"] == "batch_correct"
    assert report["status"] == "failed"
    assert report["output_written"] is False
    assert _iteration_reports(report)[0]["status"] == "failed"


def test_batch_correct_report_survives_h5ad_roundtrip(
    monkeypatch, synthetic_adata, tmp_path
):
    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)

    def fake_correct_data(adata, **kwargs):
        adata.layers[kwargs["out_layer"]] = adata.X.copy()
        return _fake_correct_report()

    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    batch_correct(
        synthetic_adata,
        markers=["CD3", "CD4"],
        xdim=3,
        ydim=3,
        return_report=True,
    )
    path = tmp_path / "batch_corrected.h5ad"

    synthetic_adata.write_h5ad(path)
    restored = read_h5ad(path)

    report = restored.uns["cycombinepy_correction"]
    assert report["function"] == "batch_correct"
    assert bool(report["output_written"]) is True
    assert list(report["markers"]) == ["CD3", "CD4"]
    assert _iteration_reports(report)[0]["function"] == "correct_data"


@pytest.mark.parametrize(
    ("ref_batch", "expected"),
    [
        (pd.Timestamp("2020-01-01"), "2020-01-01T00:00:00"),
        (pd.Interval(left=1, right=3, closed="right"), "(1, 3]"),
    ],
)
def test_batch_correct_pandas_ref_batch_reports_are_serialized(
    monkeypatch, synthetic_adata, ref_batch, expected
):
    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)

    def fake_correct_data(adata, **kwargs):
        adata.layers[kwargs["out_layer"]] = adata.X.copy()
        report = _fake_correct_report()
        report["parameters"]["ref_batch"] = kwargs["ref_batch"]
        return report

    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    report = batch_correct(
        synthetic_adata,
        markers=["CD3", "CD4"],
        xdim=3,
        ydim=3,
        ref_batch=ref_batch,
        return_report=True,
    )

    assert report["parameters"]["ref_batch"] == expected
    assert _iteration_reports(report)[0]["parameters"]["ref_batch"] == expected


def test_batch_correct_fatal_pandas_report_reraises_original_error(
    monkeypatch, synthetic_adata
):
    exception_report = _fake_correct_report(status="failed", output_written=False)
    exception_report["parameters"]["ref_batch"] = pd.Timestamp("2020-01-01")
    exception_report["clusters"]["exception"] = [pd.NA]

    monkeypatch.setattr(batch_correct_module, "create_som", _set_som_labels)

    def fake_correct_data(adata, **kwargs):
        raise CombatCorrectionError("ComBat failed for cluster 'som'", exception_report)

    monkeypatch.setattr(batch_correct_module, "correct_data", fake_correct_data)

    with pytest.raises(CombatCorrectionError, match="ComBat failed"):
        batch_correct(
            synthetic_adata,
            markers=["CD3", "CD4"],
            xdim=3,
            ydim=3,
            return_report=True,
        )

    report = synthetic_adata.uns["cycombinepy_correction"]
    assert report["status"] == "failed"
    assert _iteration_reports(report)[0]["parameters"]["ref_batch"] == (
        "2020-01-01T00:00:00"
    )
    assert _iteration_reports(report)[0]["clusters"]["exception"] == [None]
