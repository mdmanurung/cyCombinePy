import json

import numpy as np
import pytest

from cycombinepy import correct_data
from cycombinepy.correct import CORRECTED_LAYER


def _batch_mean_gap(X, batch):
    b1 = X[batch == "b1"].mean(axis=0)
    b2 = X[batch == "b2"].mean(axis=0)
    return float(np.linalg.norm(b1 - b2))


def _cluster_batch_gaps(adata, markers, layer=None):
    idx = [adata.var_names.get_loc(marker) for marker in markers]
    X = adata.layers[layer] if layer is not None else adata.X
    batch = adata.obs["batch"].to_numpy()
    labels = adata.obs["cycombine_som"].astype(str).to_numpy()
    return {
        label: _batch_mean_gap(X[labels == label][:, idx], batch[labels == label])
        for label in sorted(set(labels))
    }


def _cluster_rows(report):
    clusters = report["clusters"]
    return {
        str(label): {
            "status": str(status),
            "batches": json.loads(str(batches)),
        }
        for label, batches, status in zip(
            clusters["label"],
            clusters["batches"],
            clusters["status"],
        )
    }


@pytest.mark.requires_inmoose
def test_correct_data_reduces_balanced_additive_batch_gap_per_cluster(synthetic_adata):
    pytest.importorskip("inmoose")
    synthetic_adata.obs["cycombine_som"] = synthetic_adata.obs["celltype"].astype(
        "category"
    )
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]

    before = _cluster_batch_gaps(synthetic_adata, markers)

    report = correct_data(
        synthetic_adata,
        label_key="cycombine_som",
        markers=markers,
        return_report=True,
    )

    assert CORRECTED_LAYER in synthetic_adata.layers
    assert report["status"] == "completed"
    clusters = _cluster_rows(report)
    assert set(clusters) == set(before)
    assert all(row["status"] == "corrected" for row in clusters.values())
    assert all(row["batches"] == ["b1", "b2"] for row in clusters.values())

    after = _cluster_batch_gaps(synthetic_adata, markers, layer=CORRECTED_LAYER)
    for label, before_gap in before.items():
        assert after[label] < before_gap * 0.5, (
            f"{label} batch gap did not shrink enough "
            f"(before={before_gap}, after={after[label]})"
        )
