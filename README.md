# cyCombinePy

AnnData-native implementation inspired by
[cyCombine](https://github.com/biosurf/cyCombine) for batch correction of
single-cell cytometry data.

cyCombinePy stores data in AnnData objects and calls established Python
packages for the main numerical steps:

- **ComBat**: [`inmoose.pycombat`](https://github.com/epigenelabs/inmoose)
- **SOM clustering**: [`FlowSOM`](https://github.com/saeyslab/FlowSOM_Python)
- **FCS I/O**: [`pytometry`](https://github.com/buettnerlab/pytometry) /
  [`readfcs`](https://github.com/laminlabs/readfcs)
- **Batch-effect metrics**:
  [`scib-metrics`](https://github.com/YosefLab/scib-metrics)

## Scope and provenance

cyCombinePy implements the core AnnData workflow inspired by the R cyCombine
package: normalization for clustering, FlowSOM clustering, and per-cluster
ComBat correction.

1. **Batch-wise normalize** expression per marker (`cycombinepy.normalize`)
2. **Self-organizing map** clustering of cells (`cycombinepy.create_som`)
3. **Per-cluster ComBat** correction with optional covariates and anchors
   (`cycombinepy.correct_data`)

Step 1 operates on a normalized view so clusters are less driven by batch. Step
3 is applied to the unnormalized data per cluster so rare populations are not
over-corrected.

Out of scope for cyCombinePy are:

- Seurat / SingleCellExperiment wrappers
- Panel merging
- `ComBat_seq`
- Alternative clustering backends
- `run_analysis`

## Quickstart

```python
import cycombinepy as pc
from cycombinepy.correct import CORRECTED_LAYER
from cycombinepy.io import read_fcs_dir

# 1. Load FCS files into AnnData
adata = read_fcs_dir(
    "data/",
    metadata="metadata.csv",
    batch_key="Batch",
    sample_key="Patient",
    condition_key="condition",
    cofactor=5,           # asinh cofactor for CyTOF
)

# 2. Inspect batch effects before correction
figs = pc.detect_batch_effect_express(adata, out_dir="before/")

# 3. Run batch correction
pc.batch_correct(
    adata,
    xdim=8, ydim=8,
    covar="condition",
)
# Corrected matrix is now in adata.layers["cycombine_corrected"]

# 4. Evaluate
uncorr = pc.compute_emd(adata, cell_key="cycombine_som")
corr   = pc.compute_emd(adata, cell_key="cycombine_som", layer=CORRECTED_LAYER)
report = pc.evaluate_emd(uncorr, corr)
print(report.groupby("marker")["reduction_pct"].mean())
```

Or use the modular API:

```python
from cycombinepy.correct import CORRECTED_LAYER

pc.transform_asinh(adata, cofactor=5)
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
pc.correct_data(
    adata,
    label_key="cycombine_som",
    batch_key="batch",
    covar="condition",
    layer=None,
    out_layer=CORRECTED_LAYER,
)
```

`normalize(adata)` mutates `adata.X` unless `layer=` or `copy=True` is used.
Keeping normalization in a layer leaves `adata.X` unnormalized for correction.

Minimal synthetic example without optional FCS, plotting, FlowSOM, or ComBat
dependencies:

<!-- cycombinepy-snippet:start quickstart-synthetic -->
```python
import anndata as ad
import numpy as np
import pandas as pd

import cycombinepy as pc
from cycombinepy.correct import CORRECTED_LAYER

X = np.array(
    [
        [1.0, 2.0, 3.0],
        [1.5, 2.5, 3.5],
        [6.0, 7.0, 8.0],
        [6.5, 7.5, 8.5],
    ],
    dtype=float,
)
obs = pd.DataFrame(
    {
        "batch": ["batch_a", "batch_a", "batch_b", "batch_b"],
        "cycombine_som": ["cluster_a", "cluster_a", "cluster_b", "cluster_b"],
    },
    index=[f"cell_{i}" for i in range(4)],
)
adata = ad.AnnData(X=X, obs=obs)
adata.var_names = ["CD3", "CD19", "CD45"]

# Keep adata.X as the unnormalized correction input.
adata.layers["cycombine_normalized"] = adata.X.copy()
pc.normalize(adata, method="scale", batch_key="batch", layer="cycombine_normalized")

# Fixed labels avoid the optional FlowSOM dependency in this minimal example.
pc.correct_data(
    adata,
    label_key="cycombine_som",
    batch_key="batch",
    layer=None,
    out_layer=CORRECTED_LAYER,
)

assert CORRECTED_LAYER in adata.layers
assert np.array_equal(adata.layers[CORRECTED_LAYER], adata.X)
```
<!-- cycombinepy-snippet:end quickstart-synthetic -->

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

FCS I/O lives in `cycombinepy.io`, plotting in `cycombinepy.plotting`, and an
optional `scib_metrics` wrapper in `cycombinepy.evaluate`.

## Installation

```bash
pip install cycombinepy
pip install "cycombinepy[io,plotting]"
pip install "cycombinepy[all]"
```

For local development from a checkout:

```bash
pip install -e ".[all,dev]"
```

## Claude Code skill

cyCombinePy ships an [Agent Skill](https://docs.claude.com/en/docs/claude-code/skills)
for [Claude Code](https://claude.com/claude-code). The skill teaches the agent
the package conventions: AnnData layout, arcsinh transformation, FlowSOM
clustering, per-cluster ComBat correction, strict correction reports, and
EMD/MAD evaluation.

The skill is bundled with the Python package, but Claude Code does not scan
installed Python packages. Install it once into your personal Claude Code skills
directory:

```bash
cycombinepy-install-skills --agent claude
```

This copies the skill to `~/.claude/skills/cycombinepy/`, where it is available
to Claude Code in every project. Re-run with `--force` after upgrading
cyCombinePy to refresh the installed copy:

```bash
cycombinepy-install-skills --agent claude --force
```

Once installed, ask Claude Code to work on cyCombinePy tasks such as "set up an
AnnData object for batch correction", "run FlowSOM and per-cluster ComBat", or
"evaluate correction with EMD and MAD". The skill is a small router
(`SKILL.md`) that points to reference workflows loaded only when needed.

If you would rather not copy files into your home directory, point Claude Code
at the bundled copy in place instead:

```bash
export CLAUDE_SKILLS_PATH="$(cycombinepy-install-skills --print-path)"
```

The same bundled skill can also be installed for Codex:

```bash
cycombinepy-install-skills --agent codex
```

Running `cycombinepy-install-skills` without `--agent` installs both targets:
`~/.claude/skills/cycombinepy/` and
`${CODEX_HOME:-~/.codex}/skills/cycombinepy/`.

## Data structure conventions

- `adata.X`: cells × markers expression (post-asinh, pre-correction)
- `adata.obs["batch"]`: batch assignment (required)
- `adata.obs["sample"]`, `adata.obs["condition"]`, `adata.obs["anchor"]`:
  optional metadata
- `adata.obs["cycombine_som"]`: SOM cluster labels (written by `create_som`)
- `adata.layers["cycombine_corrected"]`: corrected expression (written by
  `correct_data` / `batch_correct`)

## Citation

If you use cyCombinePy, please cite both this software and the original
cyCombine paper. Repository metadata is available in `CITATION.cff`.

> Pedersen, C.B., Dam, S.H., Barnkob, M.B., *et al.* cyCombine allows for robust
> integration of single-cell cytometry datasets within and across technologies.
> *Nat Commun* **13**, 1698 (2022).
> <https://doi.org/10.1038/s41467-022-29383-5>
