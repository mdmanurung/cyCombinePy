"""Per-cluster ComBat correction.

Port of ``correct_data`` in ``R/02_batch_correct.R:356-544``. The AnnData is split
by its SOM cluster label, each sub-group is corrected with
:func:`cycombinepy.combat.run_combat`, and results are stitched back in the original
row order. Values are capped to the per-cluster min/max of the input (matching R
lines 524-531).
"""

from __future__ import annotations

from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import json
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from anndata import AnnData

from cycombinepy._utils import (
    check_confound,
    check_obs_values_not_missing,
    marker_matrix,
    resolve_markers,
    set_marker_matrix,
)
from cycombinepy.combat import run_combat

CORRECTED_LAYER = "cycombine_corrected"
CORRECTION_REPORT_UNS_KEY = "cycombinepy_correction"
CORRECTION_REPORT_SCHEMA_VERSION = "1.0"
ERROR_POLICIES = {"raise", "report", "warn"}
CONFOUND_POLICIES = {"raise", "skip", "drop"}
FATAL_CONFOUND_REASONS = {"confounded_with_batch", "joint_confounded_with_batch"}


class CombatCorrectionError(RuntimeError):
    """Raised when per-cluster ComBat fails under ``error_policy='raise'``."""

    def __init__(self, message: str, report: dict):
        super().__init__(message)
        self.report = report


class ConfoundedDesignError(ValueError):
    """Raised when covariates are confounded under ``confound_policy='raise'``."""

    def __init__(self, message: str, report: dict):
        super().__init__(message)
        self.report = report


def _utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _package_report() -> dict[str, str]:
    try:
        import cycombinepy

        package_version = getattr(cycombinepy, "__version__")
    except Exception:  # pragma: no cover - import-time fallback
        try:
            package_version = version("cycombinepy")
        except PackageNotFoundError:
            package_version = "unknown"
    return {"name": "cycombinepy", "version": str(package_version)}


def _json_report_value(value) -> str:
    return json.dumps(value, separators=(",", ":"), sort_keys=True)


def _build_model_matrix(
    df_sub: pd.DataFrame,
    covar: str | None,
    anchor: str | None,
) -> np.ndarray | None:
    """Build a design matrix (sans intercept) from covar and/or anchor columns.

    Uses :mod:`formulaic` to match R's ``stats::model.matrix`` (treatment
    contrasts, drop first level).
    """
    terms = [t for t in (covar, anchor) if t is not None]
    if not terms:
        return None

    sub = df_sub[terms].astype("category")
    try:
        from formulaic import model_matrix

        mm = np.asarray(model_matrix(" + ".join(terms), sub), dtype=float)
    except ModuleNotFoundError as exc:
        if exc.name != "formulaic":
            raise
        dummies = [
            pd.get_dummies(sub[term], drop_first=True, dtype=float) for term in terms
        ]
        dummies = [dummy for dummy in dummies if dummy.shape[1]]
        if not dummies:
            return None
        mm = pd.concat(dummies, axis=1).to_numpy(dtype=float)
    # Drop the intercept column so we hand inmoose a pure covariate block.
    if mm.shape[1] and np.all(mm[:, 0] == 1):
        mm = mm[:, 1:]
    return mm if mm.size else None


def _resolve_num_factors(
    series: pd.Series,
    batch: pd.Series,
    design: np.ndarray | None,
) -> int:
    """Return effective number of factor levels, mirroring R lines 455-506.

    - 1 if the covariate is confounded with batch
    - 1 if the cluster is heavily skewed to a single level
    - else the number of distinct levels.
    """
    if check_confound(batch, design):
        return 1
    counts = series.value_counts()
    total = counts.sum()
    n = counts.size
    if total < counts.max() + n * 5:
        return 1
    return n


def _resolve_factor_decision(
    term: str,
    role: str,
    series: pd.Series,
    batch: pd.Series,
    design: np.ndarray | None,
) -> tuple[int, dict[str, str] | None]:
    if check_confound(batch, design):
        return 1, {
            "term": term,
            "role": role,
            "reason": "confounded_with_batch",
        }
    counts = series.value_counts()
    total = counts.sum()
    n = counts.size
    if total < counts.max() + n * 5:
        return 1, {
            "term": term,
            "role": role,
            "reason": "skewed_to_single_level",
        }
    return n, None


def _confounded_design_message(
    cluster_label: str,
    dropped_terms: list[dict[str, str]],
) -> str:
    terms = ", ".join(
        f"{item['role']} {item['term']} ({item['reason']})"
        for item in dropped_terms
    )
    return f"Confounded design for cluster {cluster_label!r}: {terms}"


def _confounded_design_payload(
    cluster_label: str,
    dropped_terms: list[dict[str, str]],
    payload_type: str,
) -> dict[str, object]:
    return {
        "type": payload_type,
        "message": _confounded_design_message(cluster_label, dropped_terms),
        "dropped_terms": dropped_terms,
    }


def _has_fatal_confound_drop(dropped_terms: list[dict[str, str]]) -> bool:
    return any(item["reason"] in FATAL_CONFOUND_REASONS for item in dropped_terms)


def correct_data(
    adata: AnnData,
    label_key: str = "cycombine_som",
    markers: Iterable[str] | None = None,
    batch_key: str = "batch",
    covar: str | None = None,
    anchor: str | None = None,
    parametric: bool = True,
    ref_batch=None,
    layer: str | None = None,
    out_layer: str = CORRECTED_LAYER,
    copy: bool = False,
    return_report: bool = False,
    uns_key: str = CORRECTION_REPORT_UNS_KEY,
    error_policy: Literal["raise", "report", "warn"] = "raise",
    confound_policy: Literal["raise", "skip", "drop"] = "raise",
) -> AnnData | tuple[AnnData, dict] | dict | None:
    """Per-cluster ComBat batch correction.

    Parameters
    ----------
    adata
        AnnData with a cluster label in ``adata.obs[label_key]`` and a batch in
        ``adata.obs[batch_key]``.
    label_key
        Column in ``adata.obs`` with the SOM cluster id (from :func:`create_som`).
    markers
        Var names to correct. If ``None``, uses :func:`cycombinepy.get_markers`.
    batch_key
        Column in ``adata.obs`` giving the batch assignment.
    covar, anchor
        Optional ``adata.obs`` columns used as ComBat covariates. Skew- and
        confound-detection follow the R logic at lines 455-506.
    parametric
        Parametric vs. non-parametric ComBat prior.
    ref_batch
        Optional reference batch that is kept unchanged.
    layer
        If given, read the uncorrected matrix from this layer rather than ``X``.
    out_layer
        Name of the layer to store the corrected matrix in.
    copy
        If True, return a corrected copy; otherwise mutate in place.
    return_report
        If True, return the correction report dict. With ``copy=True``, returns
        ``(adata, report)``.
    uns_key
        Key in ``adata.uns`` where the correction report is stored.
    error_policy
        How to handle ComBat exceptions: raise, record in the report, or warn
        and record in the report.
    confound_policy
        How to handle requested covariates that are confounded or dropped by
        the effective design logic.
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

    check_obs_values_not_missing(adata, batch_key, context="correct_data()")
    check_obs_values_not_missing(adata, label_key, context="correct_data()")
    if covar is not None:
        check_obs_values_not_missing(adata, covar, context="correct_data()")
    if anchor is not None:
        check_obs_values_not_missing(adata, anchor, context="correct_data()")

    markers = resolve_markers(adata, markers)
    if copy:
        adata = adata.copy()

    X = marker_matrix(
        adata,
        markers,
        layer=layer,
        require_finite=True,
        context="correct_data()",
    )  # (n_cells, n_markers)

    report = {
        "schema_version": CORRECTION_REPORT_SCHEMA_VERSION,
        "function": "correct_data",
        "status": "completed",
        "output_written": False,
        "timestamp_utc": _utc_timestamp(),
        "package": _package_report(),
        "parameters": {
            "batch_key": batch_key,
            "label_key": label_key,
            "covar": covar,
            "anchor": anchor,
            "parametric": bool(parametric),
            "ref_batch": ref_batch,
            "out_layer": out_layer,
            "error_policy": error_policy,
            "confound_policy": confound_policy,
        },
        "markers": [str(marker) for marker in markers],
        "clusters": {
            "label": [],
            "n_cells": [],
            "batches": [],
            "status": [],
            "terms": [],
            "exception": [],
        },
    }
    any_cluster_failed = False
    any_adjustment = False

    # Convert label/batch to categorical codes once; group rows by label via
    # a single stable argsort + np.split to avoid the per-cluster O(N)
    # boolean scans that the previous implementation did.
    label_cat = pd.Categorical(adata.obs[label_key].astype(str).to_numpy())
    batch_cat = pd.Categorical(adata.obs[batch_key].astype(str).to_numpy())
    label_codes = label_cat.codes
    batch_codes = batch_cat.codes
    batch_categories = np.asarray(batch_cat.categories)

    order = np.argsort(label_codes, kind="stable")
    # Boundaries between sorted label groups.
    sorted_codes = label_codes[order]
    boundaries = np.flatnonzero(np.diff(sorted_codes)) + 1
    cluster_index_groups = np.split(order, boundaries)

    # Pre-slice obs columns needed for covar/anchor as integer code arrays;
    # we only materialize a small per-cluster DataFrame when _build_model_matrix
    # is actually called.
    obs = adata.obs

    corrected = X.copy()

    for idx in cluster_index_groups:
        if idx.size == 0:
            continue
        sub_X = X[idx]  # (n_sub, n_markers)
        sub_batch_codes = batch_codes[idx]

        # Detect the set of distinct batches present in this cluster without
        # falling back to pandas.
        present_codes = np.unique(sub_batch_codes)
        lab = str(label_cat.categories[label_codes[idx[0]]])
        cluster_pos = len(report["clusters"]["label"])
        report["clusters"]["label"].append(lab)
        report["clusters"]["n_cells"].append(int(idx.size))
        report["clusters"]["batches"].append(
            _json_report_value([str(batch_categories[code]) for code in present_codes])
        )
        report["clusters"]["status"].append("pending")
        report["clusters"]["terms"].append(_json_report_value([]))
        report["clusters"]["exception"].append(_json_report_value(None))
        if present_codes.size <= 1:
            # Only one batch in this cluster; nothing to correct. (R lines 448-452)
            report["clusters"]["status"][cluster_pos] = "skipped_single_batch"
            continue

        sub_batch_values = batch_categories[sub_batch_codes]
        sub_batch = pd.Series(sub_batch_values)

        # Covar / anchor handling: determine effective level count.
        num_covar = 1
        num_anchor = 1
        sub_df = None  # Built only when needed.
        dropped_terms: list[dict[str, str]] = []

        if covar is not None or anchor is not None:
            sub_df = obs.iloc[idx]

        if covar is not None:
            cov_design = _build_model_matrix(sub_df, covar, None)
            num_covar, cov_issue = _resolve_factor_decision(
                covar,
                "covar",
                sub_df[covar],
                sub_batch,
                cov_design,
            )
            if cov_issue is not None:
                dropped_terms.append(cov_issue)

        if anchor is not None:
            anc_design = _build_model_matrix(sub_df, None, anchor)
            num_anchor, anchor_issue = _resolve_factor_decision(
                anchor,
                "anchor",
                sub_df[anchor],
                sub_batch,
                anc_design,
            )
            if anchor_issue is not None:
                dropped_terms.append(anchor_issue)

        # If both are non-trivial, check that their combination is not confounded
        # with batch; if it is, drop anchor (R prioritises covar, lines 489-495).
        if num_covar > 1 and num_anchor > 1:
            joint = _build_model_matrix(sub_df, covar, anchor)
            if check_confound(sub_batch, joint):
                num_anchor = 1
                dropped_terms.append(
                    {
                        "term": anchor,
                        "role": "anchor",
                        "reason": "joint_confounded_with_batch",
                    }
                )

        if dropped_terms:
            if confound_policy == "raise" and _has_fatal_confound_drop(dropped_terms):
                report["clusters"]["status"][cluster_pos] = "failed"
                report["clusters"]["exception"][cluster_pos] = _json_report_value(
                    _confounded_design_payload(
                        lab,
                        dropped_terms,
                        "ConfoundedDesignError",
                    )
                )
                report["status"] = "failed"
                adata.uns[uns_key] = report
                raise ConfoundedDesignError(
                    _confounded_design_message(lab, dropped_terms),
                    report,
                )
            if confound_policy == "skip":
                report["clusters"]["status"][cluster_pos] = "skipped_confounded_design"
                report["clusters"]["exception"][cluster_pos] = _json_report_value(
                    _confounded_design_payload(
                        lab,
                        dropped_terms,
                        "ConfoundedDesign",
                    )
                )
                any_adjustment = True
                continue
            report["clusters"]["exception"][cluster_pos] = _json_report_value(
                _confounded_design_payload(
                    lab,
                    dropped_terms,
                    "ConfoundedDesignAdjustment",
                )
            )
            any_adjustment = True

        eff_covar = covar if num_covar > 1 else None
        eff_anchor = anchor if num_anchor > 1 else None
        if eff_covar is None and eff_anchor is None:
            mod = None
        else:
            if sub_df is None:
                sub_df = obs.iloc[idx]
            mod = _build_model_matrix(sub_df, eff_covar, eff_anchor)
        report["clusters"]["terms"][cluster_pos] = _json_report_value(
            [term for term in (eff_covar, eff_anchor) if term is not None]
        )

        # inmoose expects (n_features, n_samples) and is sensitive to float32
        # underflow in its EB priors, so upcast to float64 at the ComBat boundary.
        x_t = np.ascontiguousarray(sub_X.T, dtype=np.float64)
        try:
            corrected_sub = run_combat(
                x_t,
                batch=sub_batch_values,
                mod=mod,
                parametric=parametric,
                ref_batch=ref_batch,
            ).T
        except Exception as exc:  # pragma: no cover
            report["clusters"]["status"][cluster_pos] = "failed"
            report["clusters"]["exception"][cluster_pos] = _json_report_value(
                {
                    "type": type(exc).__name__,
                    "message": str(exc),
                }
            )
            any_cluster_failed = True
            if error_policy == "raise":
                report["status"] = "failed"
                adata.uns[uns_key] = report
                raise CombatCorrectionError(
                    f"ComBat failed for cluster {lab!r}: {exc}",
                    report,
                ) from exc
            if error_policy == "warn":
                import warnings

                warnings.warn(
                    f"ComBat failed for cluster {lab!r} ({exc}); leaving uncorrected.",
                    RuntimeWarning,
                )
            continue

        # Cap to per-marker min/max within this cluster (R lines 524-531).
        # Fuse the clip into a single maximum/minimum pair to avoid the extra
        # allocation from ``np.clip``.
        lo = sub_X.min(axis=0)
        hi = sub_X.max(axis=0)
        if corrected_sub.dtype != X.dtype:
            corrected_sub = corrected_sub.astype(X.dtype, copy=False)
        np.maximum(corrected_sub, lo, out=corrected_sub)
        np.minimum(corrected_sub, hi, out=corrected_sub)

        corrected[idx] = corrected_sub
        report["clusters"]["status"][cluster_pos] = "corrected"

    set_marker_matrix(adata, markers, corrected, layer=out_layer)
    report["output_written"] = True
    if any_cluster_failed:
        report["status"] = "completed_with_failures"
    elif any_adjustment:
        report["status"] = "completed_with_adjustments"
    adata.uns[uns_key] = report
    if copy:
        if return_report:
            return adata, report
        return adata
    if return_report:
        return report
    return None
