import matplotlib

matplotlib.use("Agg")

from cycombinepy import compute_emd
from cycombinepy.plotting import plot_density, plot_emd_heatmap


def _add_plotting_inputs(adata):
    adata.obs["cycombine_som"] = adata.obs["celltype"].astype("category")
    adata.layers["cycombine_corrected"] = adata.X.copy()


def test_plot_density_returns_figure(synthetic_adata):
    _add_plotting_inputs(synthetic_adata)
    fig = plot_density(synthetic_adata)
    assert fig is not None


def test_plot_emd_heatmap_returns_figure(synthetic_adata):
    _add_plotting_inputs(synthetic_adata)
    df = compute_emd(synthetic_adata, cell_key="cycombine_som")
    fig = plot_emd_heatmap(df)
    assert fig is not None
