import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cycombinepy.evaluate import compute_emd, compute_mad


def test_compute_emd_two_batches_one_marker_exact():
    adata = ad.AnnData(
        X=np.array([[0.0], [0.0], [1.0], [1.0]]),
        obs=pd.DataFrame(
            {
                "batch": ["a", "a", "b", "b"],
                "cycombine_som": ["c0", "c0", "c0", "c0"],
            },
            index=["0", "1", "2", "3"],
        ),
    )
    adata.var_names = ["CD0"]

    out = compute_emd(adata)

    assert len(out) == 1
    row = out.iloc[0]
    assert row[["cluster", "marker", "batch1", "batch2"]].to_dict() == {
        "cluster": "c0",
        "marker": "CD0",
        "batch1": "a",
        "batch2": "b",
    }
    assert row["emd"] == pytest.approx(1.0)


def test_compute_mad_two_batches_one_marker_exact():
    adata = ad.AnnData(
        X=np.array([[1.0], [3.0], [5.0], [2.0], [2.0], [8.0]]),
        obs=pd.DataFrame(
            {
                "batch": ["a", "a", "a", "b", "b", "b"],
                "cycombine_som": ["c0", "c0", "c0", "c0", "c0", "c0"],
            },
            index=["0", "1", "2", "3", "4", "5"],
        ),
    )
    adata.var_names = ["CD0"]

    out = compute_mad(adata)
    out = out.sort_values(["cluster", "marker", "batch"]).reset_index(drop=True)

    assert out[["cluster", "marker", "batch"]].to_dict("records") == [
        {"cluster": "c0", "marker": "CD0", "batch": "a"},
        {"cluster": "c0", "marker": "CD0", "batch": "b"},
    ]
    assert out["mad"].to_numpy() == pytest.approx([2.0, 0.0])
