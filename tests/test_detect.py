import matplotlib

matplotlib.use("Agg")

from cycombinepy.detect import detect_batch_effect_express


def test_detect_batch_effect_express_returns_figs(synthetic_adata, tmp_path):
    figs = detect_batch_effect_express(
        synthetic_adata, out_dir=tmp_path, downsample=200, seed=1
    )
    assert set(figs.keys()) >= {"emd", "density", "mds"}
    # Figures were written to disk
    assert (tmp_path / "detect_emd.png").exists()
    assert (tmp_path / "detect_density.png").exists()
    assert (tmp_path / "detect_mds.png").exists()


def test_detect_batch_effect_express_allows_missing_sample_key_column(synthetic_adata):
    figs = detect_batch_effect_express(
        synthetic_adata,
        sample_key="missing_sample",
        downsample=100,
        seed=1,
    )
    assert set(figs.keys()) >= {"emd", "density", "mds"}
