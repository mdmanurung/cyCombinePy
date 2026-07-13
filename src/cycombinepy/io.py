"""FCS I/O utilities.

Port of ``compile_fcs`` / ``convert_flowset`` / ``prepare_data`` from
``R/01_prepare_data.R``. Uses ``pytometry`` (with ``readfcs`` internally) to parse
FCS files into AnnData.
"""

from __future__ import annotations

import glob
import os
from pathlib import Path
from typing import Iterable, Literal

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData

from cycombinepy._utils import resolve_markers
from cycombinepy.preprocessing import transform_asinh


def _read_fcs_one(path: str | os.PathLike) -> AnnData:
    """Read a single FCS file into AnnData via pytometry/readfcs.

    The pytometry reader is imported lazily so that the top-level ``import
    cycombinepy.io`` works even if pytometry isn't installed.
    """
    try:
        # Import the submodule directly. pytometry's top-level package imports
        # optional dependencies that are not needed for reading FCS.
        from pytometry.io._readfcs import read_fcs as _pt_read_fcs  # type: ignore
    except ImportError:
        try:
            import readfcs  # type: ignore

            return readfcs.read(str(path))
        except ImportError as exc:
            raise ImportError(
                "read_fcs_dir requires pytometry or readfcs. Install with "
                "`pip install pytometry` or `pip install readfcs`."
            ) from exc

    return _pt_read_fcs(str(path))


def read_fcs_dir(
    data_dir: str | os.PathLike,
    pattern: str = "*.fcs",
    metadata: pd.DataFrame | str | os.PathLike | None = None,
    filename_col: str = "filename",
    sample_key: str | None = None,
    batch_key: str | None = None,
    condition_key: str | None = None,
    anchor_key: str | None = None,
    markers: Iterable[str] | None = None,
    transform: bool = True,
    cofactor: float = 5.0,
    derand: bool = True,
    downsample: int | None = None,
    sampling_type: Literal["random", "per_batch", "per_sample"] = "random",
    seed: int | None = None,
) -> AnnData:
    """Read all FCS files in ``data_dir`` into a single AnnData.

    Mirrors ``compile_fcs`` + ``convert_flowset`` + ``prepare_data`` from
    ``R/01_prepare_data.R``. Metadata (a DataFrame or a CSV/Excel path) is joined
    on the basename of each FCS file via ``filename_col``. Its columns are
    renamed to ``batch`` / ``sample`` / ``condition`` / ``anchor`` if the
    corresponding ``*_key`` argument points at them.

    Parameters
    ----------
    data_dir
        Directory containing FCS files.
    pattern
        Glob pattern for selecting files.
    metadata
        DataFrame or path to a CSV/TSV/XLSX table. Must contain ``filename_col``
        matching the FCS basenames.
    filename_col
        Column in ``metadata`` holding the FCS filenames.
    sample_key / batch_key / condition_key / anchor_key
        Columns of ``metadata`` to use for ``sample`` / ``batch`` / ``condition``
        / ``anchor`` respectively. Resulting ``adata.obs`` will use those
        canonical names.
    markers
        Restrict to these var_names after loading (optional).
    transform
        If True, apply :func:`cycombinepy.transform_asinh` with ``cofactor``.
    cofactor, derand
        Forwarded to ``transform_asinh``.
    downsample
        If given, downsample each unit (defined by ``sampling_type``) to this
        many cells.
    sampling_type
        How to downsample: uniformly at random, or per batch / per sample.
    seed
        RNG seed.
    """
    data_dir = Path(data_dir)
    files = sorted(glob.glob(str(data_dir / pattern)))
    if not files:
        raise FileNotFoundError(f"No FCS files matching {pattern} in {data_dir}")

    adatas = []
    for f in files:
        a = _read_fcs_one(f)
        a.obs["filename"] = os.path.basename(f)
        adatas.append(a)

    adata = ad.concat(adatas, join="outer", index_unique="-")

    # Join metadata if provided
    if metadata is not None:
        if isinstance(metadata, (str, os.PathLike)):
            p = Path(metadata)
            if p.suffix.lower() in (".xls", ".xlsx"):
                meta_df = pd.read_excel(p)
            elif p.suffix.lower() == ".tsv":
                meta_df = pd.read_csv(p, sep="\t")
            else:
                meta_df = pd.read_csv(p)
        else:
            meta_df = metadata.copy()

        if filename_col not in meta_df.columns:
            raise KeyError(
                f"metadata is missing the filename column {filename_col!r}"
            )

        meta_df = meta_df.set_index(filename_col)
        for col in meta_df.columns:
            adata.obs[col] = meta_df.loc[adata.obs["filename"].values, col].values

    # Canonicalize key columns
    rename_map = {
        sample_key: "sample",
        batch_key: "batch",
        condition_key: "condition",
        anchor_key: "anchor",
    }
    for src, dst in rename_map.items():
        if src is not None and src != dst and src in adata.obs.columns:
            adata.obs[dst] = adata.obs[src].values

    # Marker subset
    if markers is not None:
        keep = resolve_markers(adata, markers)
        adata = adata[:, keep].copy()

    # Downsample
    if downsample is not None:
        rng = np.random.default_rng(seed)
        if sampling_type == "random":
            if adata.n_obs > downsample:
                idx = rng.choice(adata.n_obs, downsample, replace=False)
                adata = adata[idx].copy()
        elif sampling_type in ("per_batch", "per_sample"):
            key = "batch" if sampling_type == "per_batch" else "sample"
            if key not in adata.obs.columns:
                raise KeyError(f"sampling_type={sampling_type!r} requires obs[{key!r}]")
            parts = []
            for value in adata.obs[key].unique():
                mask = (adata.obs[key] == value).to_numpy()
                subset = adata[mask]
                if subset.n_obs > downsample:
                    take = rng.choice(subset.n_obs, downsample, replace=False)
                    subset = subset[take].copy()
                parts.append(subset)
            adata = ad.concat(parts, join="outer")
        else:
            raise ValueError(f"Unknown sampling_type: {sampling_type!r}")

    # Asinh transform
    if transform:
        transform_asinh(adata, cofactor=cofactor, derand=derand, seed=seed)

    return adata
