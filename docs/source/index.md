# pycombine

`pycombine` is a Python port of
[cyCombine](https://github.com/biosurf/cyCombine) for batch correction of
single-cell cytometry data. It is AnnData-native and reuses mature Python
libraries for the numerical heavy lifting:

| Component              | Library                                                         |
| ---------------------- | --------------------------------------------------------------- |
| ComBat correction      | [`inmoose.pycombat`](https://github.com/epigenelabs/inmoose)    |
| SOM clustering         | [`FlowSOM`](https://github.com/saeyslab/FlowSOM_Python)         |
| FCS I/O                | [`pytometry`](https://github.com/buettnerlab/pytometry)         |
| Batch-effect metrics   | [`scib-metrics`](https://github.com/YosefLab/scib-metrics)      |

The pipeline ported over unchanged from the R package is:

1. **Batch-wise normalize** each marker (`pycombine.normalize`)
2. **Self-organizing map** clustering of cells (`pycombine.create_som`)
3. **Per-cluster ComBat** correction with optional covariates
   (`pycombine.correct_data`)

Step 1 operates on a normalized view so that downstream clusters represent
biology rather than technical variation. Step 3 is applied to the
*unnormalized* data per cluster so rare populations are not over-corrected.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
usage
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
