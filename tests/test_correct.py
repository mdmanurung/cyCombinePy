import numpy as np

from pycombine import correct_data, create_som, normalize
from pycombine.correct import CORRECTED_LAYER


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


def test_correct_data_reduces_batch_gap(synthetic_adata):
    _cluster_on_normalized(synthetic_adata)
    markers = [v for v in synthetic_adata.var_names if v != "FSC"]
    idx = [synthetic_adata.var_names.get_loc(m) for m in markers]
    batch = synthetic_adata.obs["batch"].to_numpy()

    before = _batch_mean_gap(synthetic_adata.X[:, idx], batch)

    correct_data(synthetic_adata, label_key="cycombine_som", markers=markers)
    corrected = synthetic_adata.layers[CORRECTED_LAYER][:, idx]

    after = _batch_mean_gap(corrected, batch)
    assert after < before * 0.5, f"batch gap did not shrink (before={before}, after={after})"


def test_correct_data_with_covar(synthetic_adata):
    _cluster_on_normalized(synthetic_adata)
    correct_data(
        synthetic_adata,
        label_key="cycombine_som",
        covar="celltype",
    )
    assert CORRECTED_LAYER in synthetic_adata.layers
