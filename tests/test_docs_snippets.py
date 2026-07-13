import numpy as np
import pandas as pd
from anndata import AnnData

import cycombinepy as pc
from cycombinepy.correct import CORRECTED_LAYER


def test_usage_modular_workflow_with_fixed_labels_avoids_optional_dependencies():
    X = np.array(
        [
            [1.0, 10.0],
            [2.0, 11.0],
            [3.0, 12.0],
            [8.0, 18.0],
            [9.0, 19.0],
            [10.0, 20.0],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(
        {
            "batch": ["b1", "b1", "b1", "b2", "b2", "b2"],
            "cycombine_som": ["cluster_a"] * 3 + ["cluster_b"] * 3,
        },
        index=[f"cell_{i}" for i in range(X.shape[0])],
    )
    adata = AnnData(X=X, obs=obs)
    adata.var_names = ["CD3", "CD4"]

    adata.layers["cycombine_normalized"] = adata.X.copy()
    pc.normalize(
        adata,
        method="scale",
        batch_key="batch",
        layer="cycombine_normalized",
    )
    pc.correct_data(
        adata,
        label_key="cycombine_som",
        batch_key="batch",
        out_layer=CORRECTED_LAYER,
    )

    assert "cycombine_normalized" in adata.layers
    assert CORRECTED_LAYER in adata.layers
    np.testing.assert_allclose(adata.layers[CORRECTED_LAYER], adata.X)

    emd = pc.compute_emd(adata, cell_key="cycombine_som")
    mad = pc.compute_mad(adata, cell_key="cycombine_som", layer=CORRECTED_LAYER)

    assert list(emd.columns) == ["cluster", "marker", "batch1", "batch2", "emd"]
    assert emd.empty
    assert set(mad.columns) == {"cluster", "marker", "batch", "mad"}
    assert set(mad["cluster"]) == {"cluster_a", "cluster_b"}
    assert np.isfinite(mad["mad"]).all()
