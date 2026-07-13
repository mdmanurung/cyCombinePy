# cyCombinePy API Workflows

## One-call correction

```python
import cycombinepy as pc

report = pc.batch_correct(
    adata,
    batch_key="batch",
    xdim=8,
    ydim=8,
    rlen=10,
    norm_method="scale",
    covar=None,
    error_policy="raise",
    confound_policy="raise",
    return_report=True,
)
```

The corrected matrix is written to
`adata.layers["cycombine_corrected"]`. The report is also stored in
`adata.uns["cycombinepy_correction"]`.

## Modular correction

Use the modular API when labels already exist or the user needs to inspect each
step.

```python
import cycombinepy as pc

adata.layers["cycombine_normalized"] = adata.X.copy()
pc.normalize(
    adata,
    method="scale",
    batch_key="batch",
    layer="cycombine_normalized",
)
pc.create_som(
    adata,
    xdim=8,
    ydim=8,
    layer="cycombine_normalized",
    label_key="cycombine_som",
)
report = pc.correct_data(
    adata,
    label_key="cycombine_som",
    batch_key="batch",
    covar=None,
    error_policy="raise",
    confound_policy="raise",
    return_report=True,
)
```

## Diagnostics and evaluation

Run diagnostics before correction:

```python
figs = pc.detect_batch_effect_express(
    adata,
    batch_key="batch",
    sample_key="sample",
    downsample=10000,
)
```

Evaluate EMD after correction:

```python
from cycombinepy.correct import CORRECTED_LAYER

emd_before = pc.compute_emd(adata, cell_key="cycombine_som")
emd_after = pc.compute_emd(
    adata,
    cell_key="cycombine_som",
    layer=CORRECTED_LAYER,
)
summary = pc.evaluate_emd(emd_before, emd_after)
```

## Common checks

- `adata.var_names` must be unique.
- Requested marker names must exist in `adata.var_names`.
- Requested layers must exist.
- Marker matrices must contain finite values.
- Fully confounded covariate or anchor designs raise
  `ConfoundedDesignError` under strict defaults.
- ComBat failures raise `CombatCorrectionError` under strict defaults.

