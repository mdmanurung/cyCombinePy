"""``batch_correct`` orchestration.

Port of ``batch_correct`` in ``R/02_batch_correct.R:66-210``. Runs the full
cyCombine pipeline: batch-wise normalization, SOM clustering, and per-cluster ComBat
correction. Supports iterative correction with multiple SOM grid sizes by
passing ``xdim``/``ydim`` as sequences.
"""

from __future__ import annotations

import json
from typing import Iterable, Literal, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from anndata import AnnData

from cycombinepy._utils import marker_matrix, resolve_markers, set_marker_matrix
from cycombinepy.cluster import create_som
from cycombinepy.correct import (
    CONFOUND_POLICIES,
    CORRECTED_LAYER,
    CORRECTION_REPORT_SCHEMA_VERSION,
    ERROR_POLICIES,
    CombatCorrectionError,
    ConfoundedDesignError,
    correct_data,
)
from cycombinepy.normalize import NormMethod, TiesMethod, normalize


def _as_list(v) -> list:
    if isinstance(v, (list, tuple, np.ndarray)):
        return list(v)
    return [v]


def _is_missing_scalar(value) -> bool:
    if isinstance(value, (dict, list, tuple, np.ndarray, pd.Index, pd.Series)):
        return False
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    if isinstance(missing, (bool, np.bool_)):
        return bool(missing)
    return False


def _normalize_report_value(value):
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if _is_missing_scalar(value):
        return None
    if isinstance(value, float):
        return None if pd.isna(value) else value
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if pd.isna(value) else float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if isinstance(value, pd.Timedelta):
        return value.isoformat() if hasattr(value, "isoformat") else str(value)
    if isinstance(value, pd.Interval):
        return str(value)
    if isinstance(value, np.ndarray):
        return _normalize_report_value(value.tolist())
    if isinstance(value, (pd.Index, pd.Series)):
        return _normalize_report_value(value.tolist())
    if isinstance(value, dict):
        return {str(k): _normalize_report_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_normalize_report_value(item) for item in value]
    return str(value)


def _report_json(report: dict) -> str:
    return json.dumps(
        _normalize_report_value(report),
        separators=(",", ":"),
        sort_keys=True,
    )


def _fallback_report_json(report, exc: Exception) -> str:
    return json.dumps(
        {
            "schema_version": CORRECTION_REPORT_SCHEMA_VERSION,
            "function": "correct_data",
            "status": "failed",
            "serialization_error": {
                "type": type(exc).__name__,
                "message": str(exc),
            },
            "report": str(report),
        },
        separators=(",", ":"),
        sort_keys=True,
    )


def _aggregate_status(iteration_reports: list[dict]) -> str:
    statuses = [str(report.get("status", "")) for report in iteration_reports]
    if any(status == "completed_with_failures" for status in statuses):
        return "completed_with_failures"
    if any(status == "completed_with_adjustments" for status in statuses):
        return "completed_with_adjustments"
    return "completed"


def _build_scratch(
    adata: AnnData,
    markers: list[str],
    working: np.ndarray,
    needed_obs_keys: list[str],
) -> AnnData:
    """Build a minimal AnnData shim for iterative normalize/cluster/correct.

    Shares only the obs columns that downstream steps read, carries
    ``working`` as ``X`` (so no full-matrix copy), and uses ``markers`` as
    ``var_names`` so ``marker_matrix`` lookups still resolve. ``obs`` is a
    shallow per-column copy so that downstream label writes don't leak back.
    """
    obs = pd.DataFrame(
        {k: adata.obs[k].values for k in needed_obs_keys if k in adata.obs.columns},
        index=adata.obs.index,
    )
    shim = ad.AnnData(X=working, obs=obs)
    shim.var_names = list(markers)
    return shim


def batch_correct(
    adata: AnnData,
    markers: Iterable[str] | None = None,
    batch_key: str = "batch",
    label_key: str = "cycombine_som",
    xdim: int | Sequence[int] = 8,
    ydim: int | Sequence[int] = 8,
    rlen: int = 10,
    seed: int = 473,
    n_clusters: int | None = None,
    norm_method: NormMethod = "scale",
    ties_method: TiesMethod = "average",
    covar: str | None = None,
    anchor: str | None = None,
    ref_batch=None,
    parametric: bool = True,
    out_layer: str = CORRECTED_LAYER,
    copy: bool = False,
    error_policy: Literal["raise", "report", "warn"] = "raise",
    confound_policy: Literal["raise", "skip", "drop"] = "raise",
    return_report: bool = False,
    uns_key: str = "cycombinepy_correction",
) -> AnnData | tuple[AnnData, dict] | dict | None:
    """Run normalization, SOM clustering, and per-cluster ComBat.

    Parameters
    ----------
    adata
        Input AnnData. ``adata.X`` is assumed to already be on an appropriate
        scale (e.g. post-asinh for cytometry).
    markers
        Var names to normalize/cluster/correct. Defaults to
        :func:`cycombinepy.get_markers`.
    batch_key
        Column in ``adata.obs`` holding batch assignments.
    label_key
        Column in ``adata.obs`` to write cluster labels to.
    xdim, ydim
        SOM grid dimensions. Sequences trigger iterative correction: for each
        ``(x, y)`` pair, re-normalize, re-cluster, and re-correct.
    rlen
        SOM training passes (forwarded to FlowSOM if supported).
    seed
        FlowSOM random seed.
    n_clusters
        If set, metacluster the SOM nodes into this many clusters.
    norm_method
        Normalization method used for clustering. See :func:`cycombinepy.normalize`.
    ties_method
        Tie-breaking rule for ``norm_method="rank"``.
    covar, anchor, ref_batch, parametric
        Forwarded to :func:`cycombinepy.correct_data`.
    out_layer
        Layer name to store the corrected matrix in.
    copy
        If True, return a corrected copy; otherwise mutate in place.
    error_policy
        Forwarded to :func:`cycombinepy.correct_data`.
    confound_policy
        Forwarded to :func:`cycombinepy.correct_data`.
    return_report
        If True, return the aggregate batch orchestration report. With
        ``copy=True``, returns ``(adata, report)``.
    uns_key
        Key in ``adata.uns`` where the aggregate batch report is stored.
    """
    if error_policy not in ERROR_POLICIES:
        raise ValueError(
            "error_policy must be one of 'raise', 'report', or 'warn'; "
            f"got {error_policy!r}"
        )
    if confound_policy not in CONFOUND_POLICIES:
        raise ValueError(
            "confound_policy must be one of 'raise', 'skip', or 'drop'; "
            f"got {confound_policy!r}"
        )

    if copy:
        adata = adata.copy()

    markers = resolve_markers(adata, markers)
    xdims = _as_list(xdim)
    ydims = _as_list(ydim)
    if len(xdims) != len(ydims):
        raise ValueError("xdim and ydim must have the same length")

    # Working copy of the marker matrix that accumulates corrections between
    # iterations. Clustering sees a normalized view; correction sees the current
    # unnormalized working state. marker_matrix already returns a fresh array, so
    # no extra .copy() is needed here.
    working = marker_matrix(adata, markers)

    # Only carry the obs columns downstream helpers actually read. This avoids
    # the large `adata.copy()` that was previously run per call.
    needed_obs_keys = [k for k in (batch_key, covar, anchor) if k is not None]
    scratch = _build_scratch(adata, markers, working, needed_obs_keys)

    iteration_reports: list[dict] = []
    aggregate_report = {
        "schema_version": CORRECTION_REPORT_SCHEMA_VERSION,
        "function": "batch_correct",
        "status": "completed",
        "output_written": False,
        "parameters": {
            "batch_key": batch_key,
            "label_key": label_key,
            "xdim": [int(x) for x in xdims],
            "ydim": [int(y) for y in ydims],
            "rlen": int(rlen),
            "seed": int(seed),
            "n_clusters": None if n_clusters is None else int(n_clusters),
            "norm_method": norm_method,
            "ties_method": ties_method,
            "covar": covar,
            "anchor": anchor,
            "ref_batch": _normalize_report_value(ref_batch),
            "parametric": bool(parametric),
            "out_layer": out_layer,
            "error_policy": error_policy,
            "confound_policy": confound_policy,
        },
        "markers": [str(marker) for marker in markers],
        "iterations": [],
    }

    for iteration, (x, y) in enumerate(zip(xdims, ydims)):
        # Normalize + cluster on a fresh normalized view. Re-seat the current
        # working state as scratch.X; _can_write_in_place handles the in-place write.
        set_marker_matrix(scratch, markers, working)
        normalize(
            scratch,
            markers=markers,
            method=norm_method,
            batch_key=batch_key,
            ties_method=ties_method,
        )
        create_som(
            scratch,
            markers=markers,
            xdim=x,
            ydim=y,
            n_clusters=n_clusters,
            seed=seed,
            rlen=rlen,
            label_key=label_key,
        )
        adata.obs[label_key] = scratch.obs[label_key].values

        # Correct the (unnormalized) working state per cluster.
        set_marker_matrix(scratch, markers, working)
        iteration_uns_key = f"{uns_key}_iteration_{iteration}"
        try:
            iteration_report = correct_data(
                scratch,
                label_key=label_key,
                markers=markers,
                batch_key=batch_key,
                covar=covar,
                anchor=anchor,
                parametric=parametric,
                ref_batch=ref_batch,
                out_layer=out_layer,
                return_report=True,
                uns_key=iteration_uns_key,
                error_policy=error_policy,
                confound_policy=confound_policy,
            )
        except (CombatCorrectionError, ConfoundedDesignError) as exc:
            exception_report = getattr(exc, "report", None)
            if exception_report is not None:
                iteration_reports.append(exception_report)
                try:
                    aggregate_report["iterations"].append(_report_json(exception_report))
                except Exception as serialization_exc:
                    aggregate_report["iterations"].append(
                        _fallback_report_json(exception_report, serialization_exc)
                    )
            aggregate_report["status"] = "failed"
            aggregate_report["output_written"] = False
            adata.uns[uns_key] = aggregate_report
            raise
        if iteration_report is None:
            iteration_report = scratch.uns[iteration_uns_key]
        iteration_reports.append(iteration_report)
        aggregate_report["iterations"].append(_report_json(iteration_report))
        working = marker_matrix(scratch, markers, layer=out_layer)

    set_marker_matrix(adata, markers, working, layer=out_layer)
    aggregate_report["status"] = _aggregate_status(iteration_reports)
    aggregate_report["output_written"] = True
    adata.uns[uns_key] = aggregate_report

    if copy:
        if return_report:
            return adata, aggregate_report
        return adata
    if return_report:
        return aggregate_report
    return None
