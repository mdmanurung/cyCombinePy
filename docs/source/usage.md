# Usage guide

This page describes the `cycombinepy` data model and correction workflow. For
complete runnable examples, see the
two tutorial notebooks:

- {doc}`notebooks/cycombine`: main batch-correction pipeline
- {doc}`notebooks/detect_batch_effects`: diagnostic plots

## Data conventions

`cycombinepy` operates on an
[`AnnData`](https://anndata.readthedocs.io/) object:

| Slot                                         | Meaning                                     |
| -------------------------------------------- | ------------------------------------------- |
| `adata.X`                                    | cells × markers expression (post-asinh)     |
| `adata.obs["batch"]`                         | batch assignment (**required**)             |
| `adata.obs["sample"]`                        | sample / patient id (optional)              |
| `adata.obs["condition"]`                     | biological condition (optional covariate)   |
| `adata.obs["anchor"]`                        | reference sample carried across batches     |
| `adata.obs["cycombine_som"]`                 | SOM cluster labels (set by `create_som`)    |
| `adata.layers["cycombine_corrected"]`        | corrected expression (set by `correct_data`)|
| `adata.uns["cycombinepy_correction"]`        | H5AD-safe correction report                 |

All functions read from `adata.X` by default; pass `layer=...` to read or
write a named layer instead. Columns other than `batch` are optional but
enable additional behavior (covariate-aware correction, anchor-based
alignment, diagnostic MDS plots, etc.).

Input validation is deliberately early. Public numerical entry points reject
missing required metadata, non-unique marker names, missing requested layers,
missing requested markers, and non-finite marker matrices before running
normalization, correction, evaluation, or plotting.

## The pipeline in one call

The most common entry point is `cycombinepy.batch_correct`, which runs the
full pipeline in one call:

```python
import cycombinepy as pc
from cycombinepy.correct import CORRECTED_LAYER

report = pc.batch_correct(
    adata,
    batch_key="batch",
    xdim=8, ydim=8,         # 64-node SOM
    rlen=10,
    norm_method="scale",    # batch-wise z-score before clustering
    covar="condition",      # optional: preserve a biological covariate
    seed=473,
    error_policy="raise",
    confound_policy="raise",
    return_report=True,
)

# The corrected matrix is now in adata.layers["cycombine_corrected"]
# (constant CORRECTED_LAYER). adata.X is unchanged, and the same report
# is stored in adata.uns["cycombinepy_correction"].
```

This call is equivalent to running three steps explicitly:

```python
from cycombinepy.correct import CORRECTED_LAYER

adata.layers["cycombine_normalized"] = adata.X.copy()
pc.normalize(adata, method="scale", batch_key="batch", layer="cycombine_normalized")
pc.create_som(
    adata,
    xdim=8,
    ydim=8,
    layer="cycombine_normalized",
    label_key="cycombine_som",
)
pc.correct_data(
    adata,
    label_key="cycombine_som",
    batch_key="batch",
    covar="condition",
    layer=None,
    out_layer=CORRECTED_LAYER,
    return_report=True,
)
```

Use the modular form to change components (for example,
`norm_method="rank"`), reuse an existing clustering, or inspect the
intermediate state.

## Loading data

### From a directory of FCS files

```python
from cycombinepy import io as pcio

adata = pcio.read_fcs_dir(
    "data/",
    metadata="metadata.csv",
    filename_col="Filename",
    batch_key="Batch",
    sample_key="Patient",
    condition_key="condition",
    cofactor=5,            # 5 CyTOF, 150 flow, 6000 spectral
    derand=True,
)
```

`read_fcs_dir` wraps
[`pytometry`](https://github.com/buettnerlab/pytometry) /
[`readfcs`](https://github.com/laminlabs/readfcs) and applies the
arcsinh transform in one step. Requires `pip install "cycombinepy[io]"`.

### From FCS files with no directory metadata

If you only have the FCS files and want to assign batches manually (this
is the path used in the {doc}`notebooks/cycombine` tutorial):

```python
import anndata as ad
import readfcs

b1 = readfcs.read("batch1.fcs"); b1.obs["batch"] = "batch1"
b2 = readfcs.read("batch2.fcs"); b2.obs["batch"] = "batch2"
adata = ad.concat([b1, b2], join="outer", index_unique="-")
pc.transform_asinh(adata, cofactor=5)
```

### From an existing AnnData

Set `adata.obs["batch"]`.

## Preprocessing

```python
pc.transform_asinh(adata, cofactor=5, derand=True)
```

- `cofactor`: marker-type-appropriate scale factor. Typical values:
  **5** (mass cytometry), **150** (conventional flow), **6000**
  (full-spectrum flow).
- `derand=True` applies the cyCombine-style derandomization trick
  (`ceil(x) - Uniform(0, 0.9999)`) to break ties at low intensity.

## Normalization

```python
pc.normalize(adata, method="scale", batch_key="batch")
```

`normalize(adata)` mutates `adata.X` unless `layer=` or `copy=True` is used.
For the modular correction workflow, normalize into a layer so `adata.X`
remains the unnormalized correction input.

Available methods:

| Method      | Notes                                                              |
| ----------- | ------------------------------------------------------------------ |
| `"scale"`   | Z-score per batch (default; matches cyCombine's `"scale"`)         |
| `"rank"`    | Percentile rank per batch (matches `"rank"` with ties handling)    |
| `"CLR"`     | Centered log-ratio, centered on the arithmetic mean of `log1p`     |
| `"CLR_seu"` | Seurat-flavor CLR                                                  |
| `"CLR_med"` | CLR centered on the median                                         |
| `"qnorm"`   | PCHIP monotone-spline quantile normalization                       |

Normalization is used for SOM clustering. The correction step itself runs on
the unnormalized expression.

## Clustering

```python
pc.create_som(adata, xdim=8, ydim=8, rlen=10, seed=473)
```

Trains a FlowSOM on the normalized marker view and writes 1-indexed
integer cluster ids to `adata.obs["cycombine_som"]`. Pass
`n_clusters=20` to metacluster the SOM nodes via consensus clustering.

## Per-cluster ComBat correction

```python
pc.correct_data(
    adata,
    label_key="cycombine_som",
    batch_key="batch",
    covar="condition",      # optional
    anchor="reference_id",  # optional
    ref_batch=None,
    parametric=True,
    error_policy="raise",
    confound_policy="raise",
    return_report=True,
)
```

Runs [ComBat](https://github.com/epigenelabs/inmoose) inside each SOM
cluster in isolation, capping the corrected values to the per-cluster
min/max of the input (matching R `cyCombine` lines 524-531). Clusters
with cells from only one batch are skipped and left unchanged. True
confounded covariate or anchor designs raise `ConfoundedDesignError` by
default; low-support or skewed terms are dropped and audited in the
correction report. Use `confound_policy="skip"` to leave affected
confounded clusters unchanged, or `confound_policy="drop"` to preserve
the legacy drop-and-audit behavior.

ComBat failures raise `CombatCorrectionError` by default and do not write a
mixed corrected/uncorrected output layer. Use `error_policy="report"` to leave
failed clusters unchanged and record the failure, or `error_policy="warn"` to
also emit a warning. With `return_report=True`, `correct_data` returns the same
report that it stores under `adata.uns["cycombinepy_correction"]`.

## Evaluating a correction

```python
from cycombinepy.correct import CORRECTED_LAYER

emd_before = pc.compute_emd(adata, cell_key="cycombine_som", layer=None)
emd_after  = pc.compute_emd(adata, cell_key="cycombine_som",
                            layer=CORRECTED_LAYER)
report = pc.evaluate_emd(emd_before, emd_after)
report.groupby("marker")["reduction_pct"].mean()
```

`compute_emd` uses
[`scipy.stats.wasserstein_distance`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wasserstein_distance.html)
to compute per-`(cluster, marker, batch-pair)` 1-D Earth Mover's Distance;
`evaluate_emd` merges uncorrected vs. corrected tables and reports
absolute and percent reduction. `compute_mad` / `evaluate_mad` provide
the same interface for Median Absolute Deviation.

## Detecting batch effects

Before correction, run the diagnostic plots to check whether the data show a
batch effect:

```python
figs = pc.detect_batch_effect_express(
    adata, batch_key="batch", sample_key="sample",
    downsample=10000, out_dir="before/",
)

# Or the extended report, which adds UMAP and per-marker MAD:
figs = pc.detect_batch_effect(
    adata, batch_key="batch", sample_key="sample",
    downsample=10000, out_dir="before/", seed=434,
)
```

Both functions return a dict of matplotlib figures and optionally save
them as PNGs under `out_dir`. See the {doc}`notebooks/detect_batch_effects`
notebook for the corresponding workflow.

## Plotting helpers

```python
from cycombinepy import plotting as pcpl
from cycombinepy.correct import CORRECTED_LAYER

pcpl.plot_dimred(adata, kind="umap", color="batch",
                 layer=CORRECTED_LAYER)
pcpl.plot_density(adata, batch_key="batch", layer=CORRECTED_LAYER)
pcpl.plot_emd_heatmap(emd_df)
```

These wrappers call `scanpy.pl.umap` and `seaborn.kdeplot` while handling the
before/after layer logic. For direct control of plot parameters, call
scanpy/seaborn directly.
