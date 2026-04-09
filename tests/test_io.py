"""Tests for FCS I/O. Skipped when the optional pytometry/readfcs deps are missing."""

from pathlib import Path

import pytest


readfcs = pytest.importorskip("readfcs")

FCS_PATH = Path(__file__).resolve().parents[1] / "inst" / "extdata" / "test.fcs"


@pytest.mark.skipif(not FCS_PATH.exists(), reason="test.fcs not available")
def test_read_fcs_dir_basic(tmp_path):
    # Copy the reference FCS so `read_fcs_dir` has a directory to scan.
    import shutil

    shutil.copy(FCS_PATH, tmp_path / "sample1.fcs")

    from cycombinepy.io import read_fcs_dir

    adata = read_fcs_dir(tmp_path, transform=False)
    assert adata.n_obs > 0
    assert adata.n_vars > 0
    assert "filename" in adata.obs.columns
    assert (adata.obs["filename"] == "sample1.fcs").all()


@pytest.mark.skipif(not FCS_PATH.exists(), reason="test.fcs not available")
def test_read_fcs_dir_with_metadata(tmp_path):
    import shutil

    shutil.copy(FCS_PATH, tmp_path / "sample1.fcs")

    import pandas as pd

    meta = pd.DataFrame(
        {"filename": ["sample1.fcs"], "Batch": ["A"], "Patient": ["P1"]}
    )
    from cycombinepy.io import read_fcs_dir

    adata = read_fcs_dir(
        tmp_path, metadata=meta, batch_key="Batch", sample_key="Patient", transform=False
    )
    assert (adata.obs["batch"] == "A").all()
    assert (adata.obs["sample"] == "P1").all()
