import matplotlib

matplotlib.use("Agg")

from pycombine import batch_correct, compute_emd
from pycombine.plotting import plot_density, plot_emd_heatmap


def test_plot_density_returns_figure(synthetic_adata):
    batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0)
    fig = plot_density(synthetic_adata)
    assert fig is not None


def test_plot_emd_heatmap_returns_figure(synthetic_adata):
    batch_correct(synthetic_adata, xdim=3, ydim=3, seed=0)
    df = compute_emd(synthetic_adata, cell_key="cycombine_som")
    fig = plot_emd_heatmap(df)
    assert fig is not None
