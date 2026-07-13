# cycombinepy

`cycombinepy` is an AnnData-based implementation inspired by
[cyCombine](https://github.com/biosurf/cyCombine) for batch correction of
single-cell cytometry data. It uses established Python packages for the main
numerical steps:

| Component              | Library                                                         |
| ---------------------- | --------------------------------------------------------------- |
| ComBat correction      | [`inmoose.pycombat`](https://github.com/epigenelabs/inmoose)    |
| SOM clustering         | [`FlowSOM`](https://github.com/saeyslab/FlowSOM_Python)         |
| FCS I/O                | [`pytometry`](https://github.com/buettnerlab/pytometry)         |
| Batch-effect metrics   | [`scib-metrics`](https://github.com/YosefLab/scib-metrics)      |

## Scope and provenance

cyCombinePy implements the core AnnData workflow inspired by the R cyCombine
package: normalization for clustering, FlowSOM clustering, and per-cluster
ComBat correction.

1. **Batch-wise normalize** each marker (`cycombinepy.normalize`)
2. **Self-organizing map** clustering of cells (`cycombinepy.create_som`)
3. **Per-cluster ComBat** correction with optional covariates
   (`cycombinepy.correct_data`)

Step 1 operates on a normalized view so that downstream clusters are less driven
by technical variation. Step 3 is applied to the
*unnormalized* data per cluster so rare populations are not over-corrected.

The API also validates requested marker names, missing observation
metadata, finite marker matrices, and requested layers before numerical
routines run. Correction functions write an H5AD-safe report to
`adata.uns["cycombinepy_correction"]`, and strict defaults fail closed for
ComBat errors or fully confounded covariate/anchor designs.

Out of scope for cyCombinePy are Seurat / SingleCellExperiment wrappers, panel
merging, `ComBat_seq`, alternative clustering backends, and `run_analysis`.

## Main entry points

- {func}`cycombinepy.batch_correct`: run normalization, FlowSOM clustering, and
  per-cluster ComBat in one call.
- {func}`cycombinepy.correct_data`: run audited per-cluster ComBat when SOM
  labels already exist.
- {func}`cycombinepy.detect_batch_effect_express` and
  {func}`cycombinepy.detect_batch_effect`: inspect marker-level and embedding
  batch effects before correction.
- {doc}`notebooks/cycombine`: source-only vignette for the correction
  workflow, including `return_report=True`, strict policies, and the
  layer-based modular API.
- {doc}`notebooks/detect_batch_effects`: source-only vignette for diagnostic
  plots and validation behavior.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
usage
citation
```

```{toctree}
:maxdepth: 2
:caption: Tutorials

notebooks/cycombine
notebooks/detect_batch_effects
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
```

## Indices

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
