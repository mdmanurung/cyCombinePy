# pycombine

Python port of [cyCombine](https://github.com/biosurf/cyCombine) for batch
correction of single-cell cytometry data.

pycombine is AnnData-native and reuses existing Python libraries instead of
reimplementing primitives:

- **ComBat**: [`inmoose.pycombat`](https://github.com/epigenelabs/inmoose)
- **SOM clustering**: [`FlowSOM`](https://github.com/saeyslab/FlowSOM_Python)
- **FCS I/O**: [`pytometry`](https://github.com/buettnerlab/pytometry) /
  [`readfcs`](https://github.com/laminlabs/readfcs)
- **Batch-effect metrics**:
  [`scib-metrics`](https://github.com/YosefLab/scib-metrics)

## Pipeline

The cyCombine workflow ports over unchanged:

1. **Batch-wise normalize** expression per marker (`pycombine.normalize`)
2. **Self-organizing map** clustering of cells (`pycombine.create_som`)
3. **Per-cluster ComBat** correction with optional covariates and anchors
   (`pycombine.correct_data`)

Step 1 operates on a normalized view so clusters represent biology rather than
batch. Step 3 is applied to the unnormalized data per cluster so rare
populations aren't over-corrected.

## Quickstart

```python
import pycombine as pc

# 1. Load FCS files into AnnData
adata = pc.io.read_fcs_dir(
    "data/",
    metadata="metadata.csv",
    batch_key="Batch",
    sample_key="Patient",
    condition_key="condition",
    cofactor=5,           # asinh cofactor for CyTOF
)

# 2. Inspect batch effects before correction
figs = pc.detect_batch_effect_express(adata, out_dir="before/")

# 3. End-to-end batch correction
pc.batch_correct(
    adata,
    xdim=8, ydim=8,
    covar="condition",
)
# Corrected matrix is now in adata.layers["cycombine_corrected"]

# 4. Evaluate
from pycombine.correct import CORRECTED_LAYER
uncorr = pc.compute_emd(adata, cell_key="cycombine_som")
corr   = pc.compute_emd(adata, cell_key="cycombine_som", layer=CORRECTED_LAYER)
report = pc.evaluate_emd(uncorr, corr)
print(report.groupby("marker")["reduction_pct"].mean())
```

Or use the modular API:

```python
pc.transform_asinh(adata, cofactor=5)
pc.normalize(adata, method="scale")
pc.create_som(adata, xdim=8, ydim=8)
pc.correct_data(adata, label_key="cycombine_som", covar="condition")
```

## Public API

| Function | Purpose |
|---|---|
| `batch_correct` | Full pipeline orchestrator |
| `transform_asinh` | Asinh transform with derandomization |
| `normalize` | Batch-wise scale / rank / CLR / qnorm |
| `create_som` | FlowSOM clustering |
| `correct_data` | Per-cluster ComBat correction |
| `compute_emd`, `evaluate_emd` | Earth-Mover's-Distance batch evaluation |
| `compute_mad`, `evaluate_mad` | Median-Absolute-Deviation batch evaluation |
| `detect_batch_effect`, `detect_batch_effect_express` | Diagnostic plots |
| `get_markers`, `check_confound` | Utilities |

FCS I/O lives in `pycombine.io`, plotting in `pycombine.plotting`, and an
optional `scib_metrics` wrapper in `pycombine.evaluate`.

## Installation

```bash
pip install -e ".[all,dev]"
```

## Data structure conventions

- `adata.X`: cells × markers expression (post-asinh, pre-correction)
- `adata.obs["batch"]`: batch assignment (required)
- `adata.obs["sample"]`, `adata.obs["condition"]`, `adata.obs["anchor"]`:
  optional metadata
- `adata.obs["cycombine_som"]`: SOM cluster labels (written by `create_som`)
- `adata.layers["cycombine_corrected"]`: corrected expression (written by
  `correct_data` / `batch_correct`)

## Citation

If you use pycombine please cite the original cyCombine paper:

> Pedersen, C.B., Dam, S.H., Barnkob, M.B., *et al.* cyCombine allows for robust
> integration of single-cell cytometry datasets within and across technologies.
> *Nat Commun* **13**, 1698 (2022).
> <https://doi.org/10.1038/s41467-022-29383-5>
